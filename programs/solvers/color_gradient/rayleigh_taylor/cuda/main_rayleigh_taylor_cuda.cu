#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#include "cuda/cuda_environment.h"
#include "cuda/cuda_error.h"
#include "cuda/cuda_memory.h"
#include "lbm/equation_of_state.h"
#include "lbm/isotropic_gradient.h"

// Lattice dimensions
const int Lx = 256;
const int Ly = 1024;
const int Q = 9;

// Simulation parameters
const double dx = 1.0;
const double dt = 1.0;
const double c_dx = 1.0e-5;
const double c_dt = c_dx / 347.0 / std::sqrt(3.0);

const int numSteps = 50000;
const int interval = 5000;
const double epsilon = 1.0e-10;

// Speed of sound
const double cs = dx / std::sqrt(3.0) / dt;
const double cs2 = cs * cs;
const double cs4 = cs2 * cs2;
const double cs6 = cs4 * cs2;

const double a_g = 9.81e2 / (c_dx / c_dt / c_dt);
const double rho1 = 4.0;
const double rho2 = 1.0;
const double c1 = 347.0 / (c_dx / c_dt);
const double c2 = 347.0 / (c_dx / c_dt);
const double c1_squared = c1 * c1;
const double c2_squared = c2 * c2;
const double radius = 10.0;
const double sigma = 0.0;
const double ch_width_init = 1.1 * dx;
const double ch_width_ope = 1.6 * dx;

const double nu = 1.0e-4 / (c_dx * c_dx / c_dt);
const double nu_b = 1.0e-4 / (c_dx * c_dx / c_dt);

const double p1_inf = rho1 * c1 * c1 - rho2 * c2 * c2 - sigma / radius;
const double p2_inf = 0.0;

// CUDA constant memory for lattice weights and discrete velocities
__constant__ double d_xi[Q][2];
__constant__ double d_w[Q];

// Stencil points for E4 in constant memory
__constant__ int d_e4_cx[8];
__constant__ int d_e4_cy[8];
__constant__ double d_e4_w[8];

__device__ inline double device_pressure(double rho_val, double phi_val, double c1_sq, double c2_sq, double p1_i, double p2_i) {
    double c_hat_squared = 0.5 * (c1_sq + c2_sq) + 0.5 * phi_val * (c1_sq - c2_sq);
    double c_bar_squared = 0.5 * (c1_sq - c2_sq) + 0.5 * phi_val * (c1_sq + c2_sq);
    double linear = p2_i - p1_i + rho_val * c_bar_squared;
    double mixing = 1.0 - phi_val * phi_val;
    if (mixing < 0.0) mixing = 0.0;
    double discriminant = linear * linear + rho_val * rho_val * mixing * c1_sq * c2_sq;
    return 0.5 * (rho_val * c_hat_squared - p1_i - p2_i + std::sqrt(discriminant));
}

__global__ void k_initialize(double* phi, double* rho, double* rho_mdt, double* p, double* p_mdt,
                             double* u, double* F, int lx, int ly, double c1_sq, double c2_sq,
                             double p1_i, double p2_i, double a_grav, double w_init, double w_ope) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= lx || j >= ly) return;

    int idx = i * ly + j;
    int y0 = ly / 2;
    double c0 = 0.1;
    double pi_val = 3.14159265358979323846;

    double phi_local = std::tanh((j - y0 - c0 * lx * std::cos(-2.0 * pi_val * i / lx)) / w_ope);
    phi[idx] = phi_local;

    u[idx * 2 + 0] = 0.0;
    u[idx * 2 + 1] = 0.0;

    double rho_local = rho1 * (0.5 + 0.5 * phi_local) + rho2 * (0.5 - 0.5 * phi_local);
    rho[idx] = rho_local;
    rho_mdt[idx] = rho_local;

    double p_local = rho_local * ((1.0 + phi_local) * 0.5 * c1_sq + (1.0 - phi_local) * 0.5 * c2_sq)
                   - (1.0 + phi_local) * 0.5 * p1_i - (1.0 - phi_local) * 0.5 * p2_i;
    p[idx] = p_local;
    p_mdt[idx] = p_local;

    F[idx * 2 + 0] = 0.0;
    F[idx * 2 + 1] = -rho_local * a_grav;
}

__global__ void k_cal_equilibrium(double* f_eq, const double* rho, const double* u, const double* p,
                                  int lx, int ly, double c_s2, double c_s4, double c_s6) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= lx || j >= ly) return;

    int idx = i * ly + j;
    double rho_local = rho[idx];
    double u_x = u[idx * 2 + 0];
    double u_y = u[idx * 2 + 1];
    double p_local = p[idx];

    for (int k = 0; k < Q; ++k) {
        double H0 = 1.0;
        double Hx = d_xi[k][0];
        double Hy = d_xi[k][1];

        double Hxx = Hx * Hx - c_s2;
        double Hyy = Hy * Hy - c_s2;
        double Hxy = Hx * Hy;

        double Hxxy = Hx * Hx * Hy - c_s2 * Hy;
        double Hyyx = Hy * Hy * Hx - c_s2 * Hx;
        double Hxxx = Hx * Hx * Hx - c_s2 * 3.0 * Hx;
        double Hyyy = Hy * Hy * Hy - c_s2 * 3.0 * Hy;

        double Hxxyy = Hx * Hx * Hy * Hy - c_s2 * (Hx * Hx + Hy * Hy) + c_s4;

        double E = d_w[k] * ((Hxx + Hyy) / (2.0 * c_s4) - Hxxyy / (4.0 * c_s6));
        double term1 = rho_local * d_w[k] * (H0 + u_x * Hx / c_s2 + u_y * Hy / c_s2
                     + 0.5 * (u_x * u_x * Hxx + 2.0 * u_x * u_y * Hxy + u_y * u_y * Hyy) / c_s4);
        f_eq[idx * Q + k] = term1 + (p_local - rho_local * c_s2) * (E + d_w[k] * (u_x * (Hyyx + Hxxx) + u_y * (Hyyy + Hxxy)) / (2.0 * c_s6));
    }
}

__global__ void k_init_distrib(double* f, double* g, const double* f_eq, const double* phi, int lx, int ly) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= lx || j >= ly) return;

    int idx = i * ly + j;
    double phi_val = phi[idx];
    for (int k = 0; k < Q; ++k) {
        double feq = f_eq[idx * Q + k];
        f[idx * Q + k] = feq;
        g[idx * Q + k] = feq * phi_val;
    }
}

__global__ void k_force(double* S, const double* u, const double* F, const double* p, const double* p_mdt,
                        const double* rho, const double* rho_mdt, int lx, int ly,
                        double c_s2, double c_s4, double c_s6, double d_t) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= lx || j >= ly) return;

    int idx = i * ly + j;
    double u_x = u[idx * 2 + 0];
    double u_y = u[idx * 2 + 1];
    double F_x = F[idx * 2 + 0];
    double F_y = F[idx * 2 + 1];

    double derive_x = 0.0;
    double derive_y = 0.0;
    for (int k = 0; k < Q; ++k) {
        int xi_x = static_cast<int>(d_xi[k][0]);
        int xi_y = static_cast<int>(d_xi[k][1]);
        int ind_x = (i + xi_x + lx) % lx;
        if (j == 0 && (k == 4 || k == 7 || k == 8)) {
            // wall boundary
        } else if (j == ly - 1 && (k == 2 || k == 5 || k == 6)) {
            // wall boundary
        } else {
            int ind_y = j + xi_y;
            int n_idx = ind_x * ly + ind_y;
            double factor = d_w[k] * (p[n_idx] - rho[n_idx] * c_s2);
            derive_x += factor * xi_x * u[n_idx * 2 + 0];
            derive_y += factor * xi_y * u[n_idx * 2 + 1];
        }
    }
    derive_x /= (d_t * c_s2);
    derive_y /= (d_t * c_s2);

    for (int k = 0; k < Q; ++k) {
        double xik0 = d_xi[k][0];
        double xik1 = d_xi[k][1];
        double Hxxk = xik0 * xik0 - c_s2;
        double Hyyk = xik1 * xik1 - c_s2;
        double Hxyk = xik0 * xik1;

        double term1 = (F_x * xik0 + F_y * xik1) / c_s2;
        double term23 = (u_x * F_x * Hxxk + u_y * F_y * Hyyk + (u_x * F_y + u_y * F_x) * Hxyk) / c_s4;
        double S_F = d_w[k] * (term1 + term23);

        double H_nu = (xik0 * xik0 - xik1 * xik1) * 0.5;
        double H_b = (xik0 * xik0 + xik1 * xik1) * 0.5 - c_s2;
        double S_Sp = d_w[k] * (derive_y * (3.0 * H_nu - H_b) + derive_x * (-3.0 * H_nu - H_b)) / (2.0 * c_s4);

        double H_xxyy = xik0 * xik0 * xik1 * xik1 - c_s2 * (xik0 * xik0 + xik1 * xik1) + c_s4;
        double E = d_w[k] * ((Hxxk + Hyyk) / (2.0 * c_s4) - H_xxyy / (4.0 * c_s6));
        double S_t = (p[idx] - p_mdt[idx] - (rho[idx] - rho_mdt[idx]) * c_s2) * E;

        S[idx * Q + k] = S_F + S_Sp + S_t;
    }
}

__global__ void k_collide(double* omega_1, const double* f, const double* f_eq, const double* S,
                          const double* rho, const double* p, int lx, int ly,
                          double nu_val, double nu_b_val, double c_s2, double c_s4, double d_t) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= lx || j >= ly) return;

    int idx = i * ly + j;
    double rho_local = rho[idx];
    double p_local = p[idx];

    double tau_nu = rho_local * nu_val / (p_local * d_t) + 0.5;
    double tau_b = rho_local * nu_b_val / (p_local * d_t) + 0.5;
    double sum_nu_neq = 0.0, sum_b_neq = 0.0, sum_xy_neq = 0.0;

    for (int k = 0; k < Q; ++k) {
        double xi_x = d_xi[k][0], xi_y = d_xi[k][1];
        double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
        double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - c_s2;
        double H_xy = xi_x * xi_y;

        double f_neq = f[idx * Q + k] - f_eq[idx * Q + k] + 0.5 * S[idx * Q + k];
        sum_nu_neq += f_neq * H_nu;
        sum_b_neq += f_neq * H_b;
        sum_xy_neq += f_neq * H_xy;
    }

    for (int k = 0; k < Q; ++k) {
        double xi_x = d_xi[k][0], xi_y = d_xi[k][1];
        double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
        double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - c_s2;
        double H_xy = xi_x * xi_y;

        double f_nu_neq = H_nu / c_s4 * sum_nu_neq;
        double f_b_neq = H_b / c_s4 * sum_b_neq;
        double f_xy_neq = H_xy / c_s4 * sum_xy_neq;
        omega_1[idx * Q + k] = d_w[k] * (1.0 - 1.0 / tau_nu) * (f_nu_neq + f_xy_neq)
                            + d_w[k] * (1.0 - 1.0 / tau_b) * f_b_neq;
    }
}

__global__ void k_collide_surface(double* omega_2, const double* phi, const double* rho, const double* p,
                                 int lx, int ly, double sigma_val, double nu_val, double nu_b_val,
                                 double c_s2, double c_s4, double d_t, double eps) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= lx || j >= ly) return;

    int idx = i * ly + j;

    double Cx = 0.0, Cy = 0.0;
    for (int n = 0; n < 8; ++n) {
        int jp = j + d_e4_cy[n];
        if (jp < 0 || jp >= ly) continue;
        int ip = (i + d_e4_cx[n] + lx) % lx;
        double val = phi[ip * ly + jp];
        Cx += d_e4_w[n] * d_e4_cx[n] * val;
        Cy += d_e4_w[n] * d_e4_cy[n] * val;
    }

    double norm_C = std::sqrt(Cx * Cx + Cy * Cy);
    double tau_nu = rho[idx] * nu_val / (p[idx] * d_t) + 0.5;
    double tau_b = rho[idx] * nu_b_val / (p[idx] * d_t) + 0.5;

    for (int k = 0; k < Q; ++k) {
        if (norm_C > eps) {
            double xi_x = d_xi[k][0], xi_y = d_xi[k][1];
            double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
            double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - c_s2;
            double H_xy = xi_x * xi_y;
            omega_2[idx * Q + k] = sigma_val * d_w[k] / (4.0 * norm_C * c_s4)
                                * ((2.0 * Cx * Cy * H_xy + (Cx * Cx - Cy * Cy) * H_nu) / tau_nu
                                - ((Cx * Cx + Cy * Cy) * H_b) / tau_b);
        } else {
            omega_2[idx * Q + k] = 0.0;
        }
    }
}

__global__ void k_recolor(double* omega_3, const double* phi, const double* p,
                          int lx, int ly, double w_ope, double c_s2, double d_t, double eps) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= lx || j >= ly) return;

    int idx = i * ly + j;

    double Cx = 0.0, Cy = 0.0;
    for (int n = 0; n < 8; ++n) {
        int jp = j + d_e4_cy[n];
        if (jp < 0 || jp >= ly) continue;
        int ip = (i + d_e4_cx[n] + lx) % lx;
        double val = phi[ip * ly + jp];
        Cx += d_e4_w[n] * d_e4_cx[n] * val;
        Cy += d_e4_w[n] * d_e4_cy[n] * val;
    }

    double grad_phi_x = Cx / d_t;
    double grad_phi_y = Cy / d_t;
    double norm_grad = std::sqrt(grad_phi_x * grad_phi_x + grad_phi_y * grad_phi_y);

    double phi_val = phi[idx];
    double p_val = p[idx];

    for (int k = 0; k < Q; ++k) {
        if (norm_grad > eps) {
            double xi_x = d_xi[k][0], xi_y = d_xi[k][1];
            omega_3[idx * Q + k] = d_w[k] * p_val * (1.0 - phi_val * phi_val) / (2.0 * w_ope)
                                * (xi_x * grad_phi_x + xi_y * grad_phi_y) / (c_s2 * norm_grad);
        } else {
            omega_3[idx * Q + k] = 0.0;
        }
    }
}

__global__ void k_stream(double* f_out, double* g_out, const double* f_eq, const double* omega_1,
                         const double* omega_2, const double* omega_3, const double* S,
                         const double* phi, const double* rho, const double* p,
                         double* rho_mdt, double* p_mdt, int lx, int ly) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= lx || j >= ly) return;

    int idx = i * ly + j;
    rho_mdt[idx] = rho[idx];
    p_mdt[idx] = p[idx];

    for (int k = 0; k < Q; ++k) {
        int ip = (i + static_cast<int>(d_xi[k][0]) + lx) % lx;
        int jp, kp;
        if (j == 0 && (k == 4 || k == 7 || k == 8)) {
            jp = j;
            kp = k - 2;
        } else if (j == ly - 1 && (k == 2 || k == 5 || k == 6)) {
            jp = j;
            kp = k + 2;
        } else {
            jp = j + static_cast<int>(d_xi[k][1]);
            kp = k;
        }
        int target_idx = (ip * ly + jp) * Q + kp;
        double f_val = f_eq[idx * Q + k] + omega_1[idx * Q + k] + omega_2[idx * Q + k] + 0.5 * S[idx * Q + k];
        f_out[target_idx] = f_val;
        g_out[target_idx] = f_val * phi[idx] + omega_3[idx * Q + k];
    }
}

__global__ void k_macroscopic(double* rho, double* u, double* F, const double* f,
                              int lx, int ly, double a_grav, double d_t) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= lx || j >= ly) return;

    int idx = i * ly + j;
    double sum_f = 0.0;
    double sum_xi_x = 0.0;
    double sum_xi_y = 0.0;

    for (int k = 0; k < Q; ++k) {
        double f_local = f[idx * Q + k];
        sum_f += f_local;
        sum_xi_x += f_local * d_xi[k][0];
        sum_xi_y += f_local * d_xi[k][1];
    }

    rho[idx] = sum_f;
    double F_x = 0.0;
    double F_y = -sum_f * a_grav;
    F[idx * 2 + 0] = F_x;
    F[idx * 2 + 1] = F_y;

    u[idx * 2 + 0] = (sum_xi_x + F_x * d_t * 0.5) / sum_f;
    u[idx * 2 + 1] = (sum_xi_y + F_y * d_t * 0.5) / sum_f;
}

__global__ void k_phase_field(double* phi, double* p, const double* f, const double* g,
                              const double* rho, int lx, int ly,
                              double c1_sq, double c2_sq, double p1_i, double p2_i) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= lx || j >= ly) return;

    int idx = i * ly + j;
    double sum_f = 0.0;
    double sum_g = 0.0;
    for (int k = 0; k < Q; ++k) {
        sum_f += f[idx * Q + k];
        sum_g += g[idx * Q + k];
    }

    double phi_val = sum_g / sum_f;
    phi[idx] = phi_val;

    double rho_local = rho[idx];
    p[idx] = device_pressure(rho_local, phi_val, c1_sq, c2_sq, p1_i, p2_i);
}

void outputDataCSV(int timestep, const double* h_rho, const double* h_u, const double* h_phi, const double* h_p) {
    std::ofstream fileDensity("density_" + std::to_string(timestep) + ".csv");
    std::ofstream fileVelocity("velocity_" + std::to_string(timestep) + ".csv");
    std::ofstream filePhase("phase_" + std::to_string(timestep) + ".csv");
    std::ofstream filePressure("pressure_" + std::to_string(timestep) + ".csv");

    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            int idx = i * Ly + j;
            fileDensity << h_rho[idx];
            fileVelocity << h_u[idx * 2 + 0] << "," << h_u[idx * 2 + 1];
            filePhase << h_phi[idx];
            filePressure << h_p[idx];

            if (i < Lx - 1) {
                fileDensity << ",";
                fileVelocity << ",";
                filePhase << ",";
                filePressure << ",";
            }
        }
        fileDensity << "\n";
        fileVelocity << "\n";
        filePhase << "\n";
        filePressure << "\n";
    }
}

int main(int argc, char** argv) {
    cglbm::lbm::GradientStencil gradient_stencil = cglbm::lbm::GradientStencil::E4;
    if (argc > 1 && !cglbm::lbm::stencil_from_name(argv[1], &gradient_stencil)) {
        std::cerr << "Unknown gradient stencil '" << argv[1] << "'; expected E4, E6 or E8." << std::endl;
        return 2;
    }

    std::cout << "colour gradient stencil = " << cglbm::lbm::stencil_name(gradient_stencil) << std::endl;
    std::cout << cglbm::cuda::describe() << std::endl;

    if (!cglbm::cuda::available()) {
        std::cerr << "Error: CUDA device not available." << std::endl;
        return 1;
    }

    cglbm::cuda::set_device(0);

    // Copy constants to device constant memory
    const double h_xi[Q][2] = {
        {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1},
        {1, 1}, {-1, 1}, {-1, -1}, {1, -1}
    };
    const double h_w[Q] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
    CGLBM_CUDA_CHECK(cudaMemcpyToSymbol(d_xi, h_xi, sizeof(h_xi)));
    CGLBM_CUDA_CHECK(cudaMemcpyToSymbol(d_w, h_w, sizeof(h_w)));

    const int h_e4_cx[8] = {1, 0, -1, 0, 1, -1, -1, 1};
    const int h_e4_cy[8] = {0, 1, 0, -1, 1, 1, -1, -1};
    const double h_e4_w[8] = {1./3., 1./3., 1./3., 1./3., 1./12., 1./12., 1./12., 1./12.};
    CGLBM_CUDA_CHECK(cudaMemcpyToSymbol(d_e4_cx, h_e4_cx, sizeof(h_e4_cx)));
    CGLBM_CUDA_CHECK(cudaMemcpyToSymbol(d_e4_cy, h_e4_cy, sizeof(h_e4_cy)));
    CGLBM_CUDA_CHECK(cudaMemcpyToSymbol(d_e4_w, h_e4_w, sizeof(h_e4_w)));

    // Allocate device arrays
    const std::size_t N = Lx * Ly;
    cglbm::cuda::DeviceBuffer<double> d_phi(N);
    cglbm::cuda::DeviceBuffer<double> d_rho(N);
    cglbm::cuda::DeviceBuffer<double> d_rho_mdt(N);
    cglbm::cuda::DeviceBuffer<double> d_p(N);
    cglbm::cuda::DeviceBuffer<double> d_p_mdt(N);
    cglbm::cuda::DeviceBuffer<double> d_u(N * 2);
    cglbm::cuda::DeviceBuffer<double> d_F(N * 2);

    cglbm::cuda::DeviceBuffer<double> d_f(N * Q);
    cglbm::cuda::DeviceBuffer<double> d_f_stream(N * Q);
    cglbm::cuda::DeviceBuffer<double> d_g(N * Q);
    cglbm::cuda::DeviceBuffer<double> d_g_stream(N * Q);
    cglbm::cuda::DeviceBuffer<double> d_f_eq(N * Q);
    cglbm::cuda::DeviceBuffer<double> d_omega_1(N * Q);
    cglbm::cuda::DeviceBuffer<double> d_omega_2(N * Q);
    cglbm::cuda::DeviceBuffer<double> d_omega_3(N * Q);
    cglbm::cuda::DeviceBuffer<double> d_S(N * Q);

    dim3 threads(16, 16);
    dim3 blocks((Lx + threads.x - 1) / threads.x, (Ly + threads.y - 1) / threads.y);

    k_initialize<<<blocks, threads>>>(d_phi.data(), d_rho.data(), d_rho_mdt.data(), d_p.data(),
                                      d_p_mdt.data(), d_u.data(), d_F.data(), Lx, Ly,
                                      c1_squared, c2_squared, p1_inf, p2_inf, a_g,
                                      ch_width_init, ch_width_ope);
    CGLBM_CUDA_CHECK(cudaGetLastError());

    k_cal_equilibrium<<<blocks, threads>>>(d_f_eq.data(), d_rho.data(), d_u.data(), d_p.data(),
                                           Lx, Ly, cs2, cs4, cs6);
    CGLBM_CUDA_CHECK(cudaGetLastError());

    k_init_distrib<<<blocks, threads>>>(d_f.data(), d_g.data(), d_f_eq.data(), d_phi.data(), Lx, Ly);
    CGLBM_CUDA_CHECK(cudaGetLastError());
    cglbm::cuda::synchronize();

    std::vector<double> h_rho(N);
    std::vector<double> h_u(N * 2);
    std::vector<double> h_phi(N);
    std::vector<double> h_p(N);

    d_rho.copy_to_host(h_rho.data());
    d_u.copy_to_host(h_u.data());
    d_phi.copy_to_host(h_phi.data());
    d_p.copy_to_host(h_p.data());
    outputDataCSV(0, h_rho.data(), h_u.data(), h_phi.data(), h_p.data());

    int steps_to_run = (argc > 2) ? std::atoi(argv[2]) : numSteps;
    int print_interval = (steps_to_run < interval) ? steps_to_run : interval;

    std::cout << "Starting CUDA simulation: " << Lx << "x" << Ly << " lattice, "
              << steps_to_run << " steps..." << std::endl;
    double t_start = cglbm::cuda::wall_time();

    double* cur_f = d_f.data();
    double* next_f = d_f_stream.data();
    double* cur_g = d_g.data();
    double* next_g = d_g_stream.data();

    for (int step = 1; step <= steps_to_run; ++step) {
        k_force<<<blocks, threads>>>(d_S.data(), d_u.data(), d_F.data(), d_p.data(), d_p_mdt.data(),
                                     d_rho.data(), d_rho_mdt.data(), Lx, Ly, cs2, cs4, cs6, dt);

        k_collide<<<blocks, threads>>>(d_omega_1.data(), cur_f, d_f_eq.data(), d_S.data(),
                                       d_rho.data(), d_p.data(), Lx, Ly, nu, nu_b, cs2, cs4, dt);

        k_collide_surface<<<blocks, threads>>>(d_omega_2.data(), d_phi.data(), d_rho.data(), d_p.data(),
                                               Lx, Ly, sigma, nu, nu_b, cs2, cs4, dt, epsilon);

        k_recolor<<<blocks, threads>>>(d_omega_3.data(), d_phi.data(), d_p.data(), Lx, Ly,
                                       ch_width_ope, cs2, dt, epsilon);

        k_stream<<<blocks, threads>>>(next_f, next_g, d_f_eq.data(), d_omega_1.data(),
                                      d_omega_2.data(), d_omega_3.data(), d_S.data(),
                                      d_phi.data(), d_rho.data(), d_p.data(),
                                      d_rho_mdt.data(), d_p_mdt.data(), Lx, Ly);

        // Swap stream buffers
        std::swap(cur_f, next_f);
        std::swap(cur_g, next_g);

        k_macroscopic<<<blocks, threads>>>(d_rho.data(), d_u.data(), d_F.data(), cur_f,
                                           Lx, Ly, a_g, dt);

        k_phase_field<<<blocks, threads>>>(d_phi.data(), d_p.data(), cur_f, cur_g,
                                           d_rho.data(), Lx, Ly, c1_squared, c2_squared,
                                           p1_inf, p2_inf);

        k_cal_equilibrium<<<blocks, threads>>>(d_f_eq.data(), d_rho.data(), d_u.data(), d_p.data(),
                                               Lx, Ly, cs2, cs4, cs6);

        if (step % print_interval == 0) {
            cglbm::cuda::synchronize();
            std::cout << "Step " << step << " / " << steps_to_run << std::endl;
            d_rho.copy_to_host(h_rho.data());
            d_u.copy_to_host(h_u.data());
            d_phi.copy_to_host(h_phi.data());
            d_p.copy_to_host(h_p.data());
            outputDataCSV(step, h_rho.data(), h_u.data(), h_phi.data(), h_p.data());
        }
    }

    cglbm::cuda::synchronize();
    double t_end = cglbm::cuda::wall_time();
    double total_time = t_end - t_start;
    double mlups = (static_cast<double>(steps_to_run) * N) / (total_time * 1.0e6);
    std::cout << "Finished simulation in " << total_time << " s (" << mlups << " MLUPs)." << std::endl;

    return 0;
}
