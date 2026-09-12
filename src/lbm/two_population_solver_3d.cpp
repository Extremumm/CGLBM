#include "lbm/two_population_solver_3d.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "lbm/field.h"
#include "lbm/output_writer.h"

namespace cglbm {
namespace lbm {

TwoPopulationSolver3D::TwoPopulationSolver3D(CaseConfig config)
    : config_(std::move(config)), lattice_(&lattice_3d(config_.lattice_3d)), nx_(config_.nx),
      ny_(config_.ny), nz_(config_.nz) {
    if (nx_ <= 0 || ny_ <= 0 || nz_ <= 0) {
        throw std::invalid_argument("lattice must have at least one node on each axis");
    }
    if (!config_.initial_phase_3d) {
        throw std::invalid_argument("case '" + config_.name + "' has no initial phase field");
    }
    const int reach = stencil_reach_3d(config_.stencil_3d);
    if (config_.boundary == Boundary::WallY && ny_ < 2 * reach + 1) {
        throw std::invalid_argument("a wall-bounded case needs room for the gradient stencil");
    }
    const Physics& physics = config_.physics;
    if (physics.rho1 <= 0.0 || physics.rho2 <= 0.0) {
        throw std::invalid_argument("both densities must be positive");
    }

    dt_ = config_.units.dt;
    cs2_ = config_.units.cs2();
    wall_y_ = config_.boundary == Boundary::WallY;
    parallel_ = config_.parallel;

    alpha2_ = physics.alpha2;
    alpha1_ = 1.0 - (1.0 - alpha2_) * physics.rho2 / physics.rho1;
    if (!(alpha1_ >= 0.0 && alpha1_ <= 1.0) || !(alpha2_ >= 0.0 && alpha2_ <= 1.0)) {
        throw std::invalid_argument(
            "alpha must lie in [0, 1] for both fluids; lower physics.alpha2 for this density "
            "ratio");
    }

    const Lattice3D& lattice = *lattice_;

    // The two moments the enhanced equilibrium's amplitude is fixed by. Summed
    // here rather than tabulated per lattice: `lambda` is exactly the quantity
    // a port carries over by hand and gets wrong, and no conservation check in
    // this file would notice. See `lattice3d.h`.
    double second_moment = 0.0;  // S = sum_q w_q e_x^2 e_y^2
    double third_moment = 0.0;   // T = sum_q w_q e_x^2 e_y^2 (3|e|^2 - 5)
    for (int q = 0; q < lattice.q; q++) {
        const double ex = lattice.xi[q][0];
        const double ey = lattice.xi[q][1];
        const double weighted = lattice.w[q] * ex * ex * ey * ey;
        second_moment += weighted;
        third_moment += weighted * (3.0 * speed_squared_3d(lattice, q) - 5.0);
    }
    if (!(std::fabs(third_moment) > 0.0)) {
        throw std::invalid_argument("this lattice cannot carry the enhanced equilibrium");
    }

    const double alpha[2] = {alpha1_, alpha2_};
    for (int fluid = 0; fluid < 2; fluid++) {
        cs_squared_[fluid] = sound_speed_squared_3d(lattice, alpha[fluid]);
        enhancement_[fluid] = (cs_squared_[fluid] - 3.0 * second_moment) / (3.0 * third_moment);
        rest_weights_[fluid].resize(static_cast<std::size_t>(lattice.q));
        for (int q = 0; q < lattice.q; q++) {
            rest_weights_[fluid][static_cast<std::size_t>(q)] =
                rest_weight_3d(lattice, q, alpha[fluid]);
        }
    }
    const double nu2 = physics.nu2 < 0.0 ? physics.nu : physics.nu2;
    mu_[0] = physics.rho1 * physics.nu;
    mu_[1] = physics.rho2 * nu2;

    rho_1_ = Field3D(nx_, ny_, nz_);
    rho_2_ = Field3D(nx_, ny_, nz_);
    rho_ = Field3D(nx_, ny_, nz_);
    p_ = Field3D(nx_, ny_, nz_);
    u_ = Field3D(nx_, ny_, nz_, 3);
    phi_n_ = Field3D(nx_, ny_, nz_);
    force_ = Field3D(nx_, ny_, nz_, 3);
    grad_x_ = Field3D(nx_, ny_, nz_);
    grad_y_ = Field3D(nx_, ny_, nz_);
    grad_z_ = Field3D(nx_, ny_, nz_);
    normal_x_ = Field3D(nx_, ny_, nz_);
    normal_y_ = Field3D(nx_, ny_, nz_);
    normal_z_ = Field3D(nx_, ny_, nz_);
    f1_ = Field3D(nx_, ny_, nz_, lattice.q);
    f2_ = Field3D(nx_, ny_, nz_, lattice.q);

    const double* drive = config_.physics.body_force;
    has_body_force_ = drive[0] != 0.0 || drive[1] != 0.0 || drive[2] != 0.0;

    if (config_.mhd.enabled) {
        mhd_ = std::make_unique<QuasiStaticMhd3D>(
            nx_, ny_, nz_, wall_y_, config_.mhd, parallel_);
    }
}

MacroscopicState3D TwoPopulationSolver3D::state() const {
    return MacroscopicState3D{&rho_, &u_, &phi_n_, &p_};
}

void TwoPopulationSolver3D::equilibrium(
    int fluid, double rho_k, double u_x, double u_y, double u_z, double* out) const {
    // The enhancement runs over 3|e|^2 - (d + 2), which is 5 in three
    // dimensions, and its amplitude was solved for at construction from the
    // lattice's own moments -- (3 (c_s^k)^2 - 1) on D3Q19, half of that on
    // D3Q27. It leaves rho and rho u untouched because
    // sum_q w_q e e (3|e|^2 - 5) vanishes on both.
    const Lattice3D& lattice = *lattice_;
    const double enhancement = enhancement_[fluid];
    const double* phi = rest_weights_[fluid].data();
    const double u_squared = u_x * u_x + u_y * u_y + u_z * u_z;
    for (int q = 0; q < lattice.q; q++) {
        const double e_squared = speed_squared_3d(lattice, q);
        const double eu =
            lattice.xi[q][0] * u_x + lattice.xi[q][1] * u_y + lattice.xi[q][2] * u_z;
        const double first = 3.0 * eu * (1.0 + enhancement * (3.0 * e_squared - 5.0));
        out[q] = rho_k * phi[static_cast<std::size_t>(q)] +
                 rho_k * lattice.w[q] * (first + 4.5 * eu * eu - 1.5 * u_squared);
    }
}

void TwoPopulationSolver3D::densities() {
    const bool parallel = parallel_;
    const double rho1_bulk = config_.physics.rho1;
    const double rho2_bulk = config_.physics.rho2;
    const int q_count = lattice_->q;
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                double sum1 = 0.0, sum2 = 0.0;
                for (int q = 0; q < q_count; q++) {
                    sum1 += f1_(i, j, k, q);
                    sum2 += f2_(i, j, k, q);
                }
                rho_1_(i, j, k) = sum1;
                rho_2_(i, j, k) = sum2;
                rho_(i, j, k) = sum1 + sum2;
                p_(i, j, k) = sum1 * cs_squared_[0] + sum2 * cs_squared_[1];
                const double scaled1 = sum1 / rho1_bulk;
                const double scaled2 = sum2 / rho2_bulk;
                const double total = scaled1 + scaled2;
                phi_n_(i, j, k) = total > 0.0 ? (scaled1 - scaled2) / total : 0.0;
            }
        }
    }
}

void TwoPopulationSolver3D::gradient_at(const double* field,
                                        int i,
                                        int j,
                                        int k,
                                        double* grad_x,
                                        double* grad_y,
                                        double* grad_z) const {
    if (wall_y_) {
        gradient_wall_y_3d(
            field, nx_, ny_, nz_, i, j, k, config_.stencil_3d, grad_x, grad_y, grad_z);
    } else {
        gradient_periodic_3d(
            field, nx_, ny_, nz_, i, j, k, config_.stencil_3d, grad_x, grad_y, grad_z);
    }
}

void TwoPopulationSolver3D::update_colour_gradient() {
    const bool parallel = parallel_;
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                double gx = 0.0, gy = 0.0, gz = 0.0;
                gradient_at(phi_n_.data(), i, j, k, &gx, &gy, &gz);
                grad_x_(i, j, k) = gx;
                grad_y_(i, j, k) = gy;
                grad_z_(i, j, k) = gz;
                const double norm = std::sqrt(gx * gx + gy * gy + gz * gz);
                if (norm > kInterfaceGradientFloor) {
                    normal_x_(i, j, k) = -gx / norm;
                    normal_y_(i, j, k) = -gy / norm;
                    normal_z_(i, j, k) = -gz / norm;
                } else {
                    normal_x_(i, j, k) = 0.0;
                    normal_y_(i, j, k) = 0.0;
                    normal_z_(i, j, k) = 0.0;
                }
            }
        }
    }
}

void TwoPopulationSolver3D::surface_force() {
    const bool parallel = parallel_;
    const double sigma = config_.physics.sigma;
    const double gravity = config_.physics.gravity;
    const bool has_gravity = gravity != 0.0;
    const bool has_body_force = has_body_force_;
    const double* drive = config_.physics.body_force;
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                // K = -(div n - n_a n_b d_a n_b), the surface divergence. The
                // second term vanishes for a unit normal and keeps the estimate
                // bounded where the discrete one is not quite a unit vector.
                double d[3][3];
                gradient_at(normal_x_.data(), i, j, k, &d[0][0], &d[0][1], &d[0][2]);
                gradient_at(normal_y_.data(), i, j, k, &d[1][0], &d[1][1], &d[1][2]);
                gradient_at(normal_z_.data(), i, j, k, &d[2][0], &d[2][1], &d[2][2]);
                const double n[3] = {normal_x_(i, j, k), normal_y_(i, j, k), normal_z_(i, j, k)};

                double divergence = 0.0, projection = 0.0;
                for (int a = 0; a < 3; a++) {
                    divergence += d[a][a];  // d_a n_a
                    for (int b = 0; b < 3; b++) {
                        projection += n[a] * n[b] * d[b][a];  // n_a n_b d_b n_a
                    }
                }
                const double curvature = -(divergence - projection);

                const double gx = grad_x_(i, j, k);
                const double gy = grad_y_(i, j, k);
                const double gz = grad_z_(i, j, k);
                const bool on_interface =
                    std::sqrt(gx * gx + gy * gy + gz * gz) > kInterfaceGradientFloor;
                const double magnitude = on_interface ? -0.5 * sigma * curvature : 0.0;
                force_(i, j, k, 0) = magnitude * gx;
                force_(i, j, k, 1) = magnitude * gy;
                force_(i, j, k, 2) = magnitude * gz;
                if (has_gravity) {
                    force_(i, j, k, 1) -= rho_(i, j, k) * gravity;
                }
                if (has_body_force) {
                    for (int a = 0; a < 3; a++) {
                        force_(i, j, k, a) += drive[a];
                    }
                }
            }
        }
    }

    if (mhd_) {
        // `u_` is the bare momentum here, not the reported velocity: see the
        // header. The potential solve is warm-started from the previous step,
        // so this costs a handful of iterations rather than a fresh solve.
        mhd_->solve(u_, phi_n_);
        mhd_->add_lorentz_force(force_);
    }
}

void TwoPopulationSolver3D::update_velocity(bool with_force) {
    const bool parallel = parallel_;
    const int q_count = lattice_->q;
    const Lattice3D& lattice = *lattice_;
    const double half_step = with_force ? 0.5 * dt_ : 0.0;
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                double sum[3] = {0.0, 0.0, 0.0};
                for (int q = 0; q < q_count; q++) {
                    const double total = f1_(i, j, k, q) + f2_(i, j, k, q);
                    for (int a = 0; a < 3; a++) {
                        sum[a] += total * lattice.xi[q][a];
                    }
                }
                const double density = rho_(i, j, k);
                for (int a = 0; a < 3; a++) {
                    u_(i, j, k, a) = (sum[a] + half_step * force_(i, j, k, a)) / density;
                }
            }
        }
    }
}

void TwoPopulationSolver3D::collide() {
    const bool parallel = parallel_;
    const int q_count = lattice_->q;
    const Lattice3D& lattice = *lattice_;
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                const double u_x = u_(i, j, k, 0);
                const double u_y = u_(i, j, k, 1);
                const double u_z = u_(i, j, k, 2);
                const double density = rho_(i, j, k);
                const double fraction = 0.5 * (1.0 + phi_n_(i, j, k));
                const double mu = fraction * mu_[0] + (1.0 - fraction) * mu_[1];
                const double tau = mu / (p_(i, j, k) * dt_) + 0.5;
                const double omega = 1.0 / tau;
                const double guo = (1.0 - 0.5 * omega) * dt_;
                const double share1 = rho_1_(i, j, k) / density;

                const double f[3] = {force_(i, j, k, 0), force_(i, j, k, 1), force_(i, j, k, 2)};

                double eq1[kQ3D27], eq2[kQ3D27];
                equilibrium(0, rho_1_(i, j, k), u_x, u_y, u_z, eq1);
                equilibrium(1, rho_2_(i, j, k), u_x, u_y, u_z, eq2);

                for (int q = 0; q < q_count; q++) {
                    const double* e = lattice.xi[q];
                    const double eu = e[0] * u_x + e[1] * u_y + e[2] * u_z;
                    double source = 0.0;
                    source += ((e[0] - u_x) + eu * e[0] / cs2_) * f[0];
                    source += ((e[1] - u_y) + eu * e[1] / cs2_) * f[1];
                    source += ((e[2] - u_z) + eu * e[2] / cs2_) * f[2];
                    source *= lattice.w[q] / cs2_;

                    f1_(i, j, k, q) += -omega * (f1_(i, j, k, q) - eq1[q]) + guo * share1 * source;
                    f2_(i, j, k, q) +=
                        -omega * (f2_(i, j, k, q) - eq2[q]) + guo * (1.0 - share1) * source;
                }
            }
        }
    }
}

double TwoPopulationSolver3D::rest_weight(int i, int j, int k, int q) const {
    const double phi_1 = rest_weights_[0][static_cast<std::size_t>(q)];
    const double phi_2 = rest_weights_[1][static_cast<std::size_t>(q)];
    return (rho_1_(i, j, k) * phi_1 + rho_2_(i, j, k) * phi_2) / rho_(i, j, k);
}

void TwoPopulationSolver3D::recolor() {
    const bool parallel = parallel_;
    const double beta = config_.physics.beta;
    const int q_count = lattice_->q;
    const Lattice3D& lattice = *lattice_;
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                const double density = rho_(i, j, k);
                const double share1 = rho_1_(i, j, k) / density;
                const double gx = grad_x_(i, j, k);
                const double gy = grad_y_(i, j, k);
                const double gz = grad_z_(i, j, k);
                const double norm = std::sqrt(gx * gx + gy * gy + gz * gz);
                const double strength =
                    norm > kGradientEpsilon
                        ? beta * rho_1_(i, j, k) * rho_2_(i, j, k) / (density * norm)
                        : 0.0;
                for (int q = 0; q < q_count; q++) {
                    const double total = f1_(i, j, k, q) + f2_(i, j, k, q);
                    double push = 0.0;
                    if (q != 0 && strength != 0.0) {
                        const double* e = lattice.xi[q];
                        const double speed = std::sqrt(speed_squared_3d(lattice, q));
                        push = rest_weight(i, j, k, q) * strength *
                               (e[0] * gx + e[1] * gy + e[2] * gz) / speed;
                    }
                    f1_(i, j, k, q) = share1 * total + push;
                    f2_(i, j, k, q) = total - f1_(i, j, k, q);
                }
            }
        }
    }
}

void TwoPopulationSolver3D::stream() {
    const Lattice3D& lattice = *lattice_;
    const int q_count = lattice.q;
    Field3D next1(nx_, ny_, nz_, q_count);
    Field3D next2(nx_, ny_, nz_, q_count);
    const bool parallel = parallel_;
    const bool wall = wall_y_;
#pragma omp parallel for collapse(3) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                for (int q = 0; q < q_count; q++) {
                    const int dx = static_cast<int>(lattice.xi[q][0]);
                    const int dy = static_cast<int>(lattice.xi[q][1]);
                    const int dz = static_cast<int>(lattice.xi[q][2]);
                    const int ip = (i + dx + nx_) % nx_;
                    const int kp = (k + dz + nz_) % nz_;
                    int jp = j + dy;
                    int qp = q;
                    if (wall && (jp < 0 || jp >= ny_)) {
                        // Half-way bounce-back: the population stays put and
                        // comes back along the reversed direction.
                        jp = j;
                        qp = lattice.opposite[q];
                    } else if (!wall) {
                        jp = (j + dy + ny_) % ny_;
                    }
                    next1(ip, jp, kp, qp) = f1_(i, j, k, q);
                    next2(ip, jp, kp, qp) = f2_(i, j, k, q);
                }
            }
        }
    }
    f1_ = std::move(next1);
    f2_ = std::move(next2);
}

void TwoPopulationSolver3D::initialize() {
    const Physics& physics = config_.physics;
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < nz_; k++) {
                double indicator = config_.initial_phase_3d(config_, i, j, k);
                if (indicator > 1.0) {
                    indicator = 1.0;
                } else if (indicator < -1.0) {
                    indicator = -1.0;
                }
                const double fraction = 0.5 * (1.0 + indicator);
                double eq1[kQ3D27], eq2[kQ3D27];
                equilibrium(0, fraction * physics.rho1, 0.0, 0.0, 0.0, eq1);
                equilibrium(1, (1.0 - fraction) * physics.rho2, 0.0, 0.0, 0.0, eq2);
                for (int q = 0; q < lattice_->q; q++) {
                    f1_(i, j, k, q) = eq1[q];
                    f2_(i, j, k, q) = eq2[q];
                }
            }
        }
    }
    densities();
    update_colour_gradient();
    if (mhd_) {
        update_velocity(false);
    }
    surface_force();
    update_velocity(true);
}

void TwoPopulationSolver3D::refresh() {
    densities();
    update_velocity(true);
}

void TwoPopulationSolver3D::step() {
    densities();
    update_colour_gradient();
    if (mhd_) {
        // The magnetic force is evaluated at the bare momentum, before Guo's
        // half step folds the force back in; see the header.
        update_velocity(false);
    }
    surface_force();
    update_velocity(true);
    collide();
    recolor();
    stream();
}

void TwoPopulationSolver3D::interface_axes(double* radii) const {
    const int centre[3] = {nx_ / 2, ny_ / 2, nz_ / 2};
    const int extent[3] = {nx_, ny_, nz_};
    for (int axis = 0; axis < 3; axis++) {
        radii[axis] = std::numeric_limits<double>::quiet_NaN();
        const int reach = extent[axis] - centre[axis] - 1;
        double previous = phi_n_(centre[0], centre[1], centre[2]);
        for (int m = 1; m <= reach; m++) {
            int node[3] = {centre[0], centre[1], centre[2]};
            node[axis] += m;
            const double value = phi_n_(node[0], node[1], node[2]);
            if (previous >= 0.0 && value < 0.0) {
                // The crossing sits between m - 1 and m, at the linear
                // interpolant of the two values bracketing it.
                radii[axis] = (m - 1) + previous / (previous - value);
                break;
            }
            previous = value;
        }
    }
}

void TwoPopulationSolver3D::write_midplane(CsvWriter& writer, int timestep) const {
    const int k = nz_ / 2;
    Field density(nx_, ny_), phase(nx_, ny_), pressure(nx_, ny_), velocity(nx_, ny_, 2);
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            density(i, j) = rho_(i, j, k);
            phase(i, j) = phi_n_(i, j, k);
            pressure(i, j) = p_(i, j, k);
            velocity(i, j, 0) = u_(i, j, k, 0);
            velocity(i, j, 1) = u_(i, j, k, 1);
        }
    }
    writer.write_grids(timestep, MacroscopicState{&density, &velocity, &phase, &pressure});
}

void write_vtk_3d(const std::string& prefix,
                  int timestep,
                  const Field3D& density,
                  const Field3D& velocity,
                  const Field3D& phase,
                  const Field3D& pressure) {
    const std::string name = prefix + "_" + std::to_string(timestep) + ".vtk";
    std::ofstream file(name);
    if (!file) {
        throw OutputError("cannot open '" + name + "' for writing");
    }
    const int nx = density.nx(), ny = density.ny(), nz = density.nz();
    file.precision(8);
    file << "# vtk DataFile Version 3.0\nCGLBM " << timestep
         << "\nASCII\nDATASET STRUCTURED_POINTS\n"
         << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n"
         << "ORIGIN 0 0 0\nSPACING 1 1 1\n"
         << "POINT_DATA " << static_cast<long>(nx) * ny * nz << "\n";
    const Field3D* scalars[3] = {&density, &phase, &pressure};
    const char* names[3] = {"density", "phase", "pressure"};
    for (int s = 0; s < 3; s++) {
        file << "SCALARS " << names[s] << " double 1\nLOOKUP_TABLE default\n";
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    file << (*scalars[s])(i, j, k) << "\n";
                }
            }
        }
    }
    file << "VECTORS velocity double\n";
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                file << velocity(i, j, k, 0) << " " << velocity(i, j, k, 1) << " "
                     << velocity(i, j, k, 2) << "\n";
            }
        }
    }
    file.close();
    if (!file) {
        throw OutputError("failed while writing '" + name + "'");
    }
}

void TwoPopulationSolver3D::run() {
    CsvWriter writer(config_.output_precision);
    const bool track = config_.track_interface;
    double radii[3];

    initialize();
    write_midplane(writer, 0);
    if (config_.write_vtk_field) {
        write_vtk_3d("field", 0, rho_, u_, phi_n_, p_);
    }
    if (track) {
        writer.open_droplet_track();
        interface_axes(radii);
        writer.write_axes(0, radii);
    }

    for (int timestep = 1; timestep <= config_.steps; timestep++) {
        step();
        if (track) {
            // `step()` leaves phi_n_ as it stood before the collision it
            // performed, so the fields are recovered from the streamed
            // populations before the track is taken: the line written against
            // `timestep` is then the state at `timestep`, not one step behind.
            densities();
            interface_axes(radii);
            writer.write_axes(timestep, radii);
        }
        if (timestep % config_.interval == 0) {
            std::cout << "Step " << timestep;
            if (mhd_) {
                // What the potential cost and whether it got there. A run that
                // stops converging says so here rather than in the flow field
                // three thousand steps later.
                std::cout << "  potential: " << mhd_->iterations() << " iterations, residual "
                          << mhd_->residual() << ", charge imbalance "
                          << mhd_->charge_imbalance();
                if (!mhd_->converged()) {
                    std::cout << "  [did not reach tolerance]";
                }
            }
            std::cout << std::endl;
            refresh();
            write_midplane(writer, timestep);
            if (config_.write_vtk_field) {
                write_vtk_3d("field", timestep, rho_, u_, phi_n_, p_);
            }
        }
    }
}

}  // namespace lbm
}  // namespace cglbm
