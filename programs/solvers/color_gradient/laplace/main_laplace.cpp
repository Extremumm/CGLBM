#include <iostream>
#include <cmath>
#include <cstdlib>

#include "lbm/equation_of_state.h"
#include "lbm/isotropic_gradient.h"
#include "lbm/mixture.h"
#include "lbm/surface_force.h"

// Usage: laplace [E4|E6|E8] [density_ratio] [viscosity_ratio]
//
//   density_ratio    rho1/rho2, default 20. The light fluid is fixed (rho2, nu2)
//                    and the droplet is made heavier.
//   viscosity_ratio  mu1/mu2, default equal to the density ratio, i.e. the same
//                    kinematic viscosity in both fluids as in the shipped case.
//                    At large density ratios give it explicitly: a single
//                    kinematic viscosity makes the droplet's relaxation time
//                    grow with the density ratio (tau ~ 5e4 at 1e4), which the
//                    scheme does not survive. `laplace E8 1e4 1` is the
//                    high-density-ratio benchmark of docs/numerics.md.

// In this version, all the intermediate variables are calculated for clarity. The code is not optimized for performance.
// Periodic boundary conditions p170 of the book (Graduate Texts in Physics) Timm Krüger, Halim Kusumaatmaja, Alexandr Kuzmin, Orest Shardt, Goncalo Silva, Erlend Magnus Viggen (auth.) - The Lattice Boltzmann Method_ Principles and Practice-Springer
// Initial conditions depend on the problem

//How to define Reynolds number in LBM? Convert the Reynolds number to physical units from the lattice units

// Define your lattice dimensions and parameters
const int Lx = 128; // Number of lattice nodes in the x-direction
const int Ly = 128; // Number of lattice nodes in the y-direction
const int Q = 9;  // Number of discrete velocities


// Define your simulation parameters without dimensions
const double dx = 1.; // Lattice spacing : 1e-3 m
const double dt = 1.; // Time step : 16638.0e-6 s

const double c_dx = 1.e-5 ; // m : conversion factor from lattice units to physical units
const double c_dt = c_dx/347./sqrt(3.); // s : conversion factor from lattice units to physical units

const int numSteps = 30000; // Number of simulation steps
const int interval = 1000; // Output interval

// Speed of sound and related constants
const double cs = dx / sqrt(3.0) / dt; // Speed of sound in the lattice
const double cs2 = cs * cs;
const double cs4 = cs2 * cs2;
const double cs6 = cs4 * cs2;

//parameters
double rho1 = 20.; // kg * m-3 Density for component 1, overridable from the command line
const double rho2 = 1.;  // kg * m-3 Density for component 2
const double c1 = 347./(c_dx/c_dt);//dx / sqrt(3.0) / dt; // Speed of sound for component 1 in lattice units
const double c2 = 347./(c_dx/c_dt);//dx / sqrt(3.0) / dt; // Speed of sound for component 2 in lattice units
const double radius = 10.; // Radius of the droplet
const double sigma = 1./(c_dx * c_dx * c_dx / c_dt / c_dt); //1.e-3 / (347*347*3*1.e-4);// 1e-3; //s1.e-3; // surface tension
//const double sigma = (rho1 - rho2) * cs2 * radius;
const double ch_width_init = 1.1 * dx; // Characteristic width of the interface
const double ch_width_ope = 1.6 * dx;

// viscosities to calculate Relaxation times in collide step
// Component 2, the surrounding fluid, keeps the viscosity the case always had;
// component 1 gets its own, from the viscosity ratio (see setComponents).
const double nu2   = 1.e-2/(c_dx * c_dx / c_dt) ; //1.0e-3*sqrt(3)/(347*1.0e-3); // kinematic viscosity in lattice units p284 of the book (Graduate Texts in Physics) Timm Krüger, Halim Kusumaatmaja, Alexandr Kuzmin, Orest Shardt, Goncalo Silva, Erlend Magnus Viggen (auth.) - The Lattice Boltzmann Method_ Principles and Practice-Springer
const double nu_b2 = 1.e-2/(c_dx * c_dx / c_dt) ; //1.0e-3*sqrt(3)/(347*1.0e-3); // bulk viscosity in lattice units
double nu1   = nu2;   // kinematic viscosity of component 1, lattice units
double nu_b1 = nu_b2; // bulk viscosity of component 1, lattice units

//The pressure at infinity is used at the equation of state to calculate the pressure and equilibrium distribution function
double p1_inf = rho1 * c1 * c1 - rho2 * c2 * c2 - sigma / radius; // Pressure at infinity for component 1
const double p2_inf = 0.; // Pressure at infinity for component 2

// Derive everything that depends on the density and viscosity ratios.
//
// p1_inf is chosen so that component 1 sits at density rho1 under the Laplace
// pressure p2 + sigma/R while component 2 sits at rho2 under p2 = rho2 c2^2.
// The mixture viscosity is rho (Y1 nu1 + Y2 nu2), i.e. the volume-weighted
// dynamic viscosity alpha1 mu1 + alpha2 mu2 (src/lbm/mixture.h); with
// viscosity_ratio == density_ratio, nu1 == nu2 and the original single-
// viscosity case is recovered exactly.
void setComponents(double density_ratio, double viscosity_ratio) {
    rho1 = density_ratio * rho2;
    p1_inf = rho1 * c1 * c1 - rho2 * c2 * c2 - sigma / radius;
    nu1 = viscosity_ratio * nu2 * rho2 / rho1;
    nu_b1 = viscosity_ratio * nu_b2 * rho2 / rho1;
}

cglbm::lbm::ComponentPair componentPair() {
    return {c1 * c1, c2 * c2, p1_inf, p2_inf};
}

double S[Lx][Ly][Q]; // Force term
double F[Lx][Ly][2]; // Volumic force : kg * m^-2 * s^-2. Here the surface tension, see calSurfaceForce
//for Laplace equation test, the gravity is on the perpendicular direction to the interface: 0
//for instability test, the gravity is on the parallel direction to the interface

double rho[Lx][Ly]; // Density
double rho_mdt[Lx][Ly]; // Density at the previous time step
double u[Lx][Ly][2]; // Velocity components
double momentum[Lx][Ly][2]; // First moment of f, before the half-force correction of the velocity
double p[Lx][Ly]; // Pressure
double p_mdt[Lx][Ly]; // Pressure at the previous time step
double phi[Lx][Ly]; // Phase field function: mass-fraction difference Y1 - Y2
double psi[Lx][Ly]; // Normalised colour field: volume-fraction difference alpha1 - alpha2
double stress_xx[Lx][Ly]; // Capillary stress sigma/2 (|grad psi| I - grad psi grad psi / |grad psi|), xx
double stress_xy[Lx][Ly]; // Capillary stress, xy
double stress_yy[Lx][Ly]; // Capillary stress, yy
double normal_x[Lx][Ly]; // Unit interface normal grad(psi)/|grad(psi)|, x component
double normal_y[Lx][Ly]; // Unit interface normal grad(psi)/|grad(psi)|, y component

double omega_1[Lx][Ly][Q]; // Collision term
double omega_3[Lx][Ly][Q]; // Recoloring term
double f[Lx][Ly][Q]; // Sum of distribution functions
double g[Lx][Ly][Q]; // Difference of distribution functions
double f_eq[Lx][Ly][Q]; // Sum of distribution functions at equilibrium

const double xi[Q][2] = {
    {0, 0},
    {1, 0}, {0, 1}, {-1, 0}, {0, -1},
    {1, 1}, {-1, 1}, {-1, -1}, {1, -1}
};

//weights
const double w1 = 4./9.;
const double w2 = 1./9.;
const double w3 = 1./36.;
const double w[Q] = {w1, w2, w2, w2, w2, w3, w3, w3, w3};



// Isotropy order of the colour gradient, overridable from the command line
// (e.g. `laplace E4`). See src/lbm/isotropic_gradient.h.
//
// Leclaire, Reggio & Trepanier, Computers & Fluids 48, 98 (2011) show that the
// nearest-neighbour gradient is what limits the colour-gradient model at large
// density contrast, and that a higher-order isotropic stencil lets Laplace's
// law hold to O(10^4).

cglbm::lbm::GradientStencil gradient_stencil = cglbm::lbm::GradientStencil::E8;

void calEquilibrium() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double rho_local = rho[i][j];  // Local density
            double u_x = u[i][j][0];
            double u_y = u[i][j][1];
            double p_local = p[i][j];
            for (int k = 0; k < Q; k++) {
                double H0 = 1.0;
                double Hx = xi[k][0];
                double Hy = xi[k][1];

                double Hxx = Hx * Hx - cs2 ; // xi[k][0]*xi[k][0] - cs2;
                double Hyy = Hy * Hy - cs2 ; // xi[k][1]*xi[k][1] - cs2;
                double Hxy = xi[k][0]*xi[k][1]; // Hyx = Hxy

                double Hxxy = Hx * Hx * Hy - cs2 * Hy;
                double Hyyx = Hy * Hy * Hx - cs2 * Hx;
                double Hxxx = pow(Hx, 3) - cs2 * 3. * Hx;
                double Hyyy = pow(Hy, 3) - cs2 * 3. * Hy;

                double Hxxyy = Hx * Hx * Hy * Hy - cs2 * (Hx * Hx + Hy * Hy) + cs4;

                double E = w[k]*((Hxx + Hyy) / (2.*cs4) - Hxxyy / (4.*cs6));
                double term1 = rho_local * w[k] * (H0 + u_x * Hx / cs2 + u_y * Hy / cs2 + 0.5 * (u_x * u_x * Hxx + u_x * u_y * Hxy * 2. + u_y * u_y * Hyy) / cs4);
                f_eq[i][j][k] = term1 + (p_local - rho_local * cs2) * (E +  w[k] * (u_x * (Hyyx + Hxxx) + u_y * (Hyyy + Hxxy))/(2.*cs6));
                // f_eq[i][j][k] = term1 + (p_local - rho_local * cs2) * (E +  w[k] * (u_x * (Hyyx) + u_y * (Hxxy))/(2.*cs6));
                // std::cout << "k = " << k << " f_eq = " << f_eq[i][j][k] << std::endl;
            }
        }
    }
}

// Density and momentum. The velocity needs the force, which needs the phase
// field, so it is completed later by calVelocity.
void calMacroscopic() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double sum_f = 0.0;
            double sum_xi_x = 0.0;
            double sum_xi_y = 0.0;
            for (int k = 0; k < Q; k++) {
                double f_local = f[i][j][k];
                sum_f += f_local;
                sum_xi_x += f_local * xi[k][0];
                sum_xi_y += f_local * xi[k][1];
            }
            rho[i][j] = sum_f;
            momentum[i][j][0] = sum_xi_x;
            momentum[i][j][1] = sum_xi_y;
        }
    }
}

// Velocity with the half-force correction of the Guo et al. forcing scheme.
void calVelocity() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            u[i][j][0] = (momentum[i][j][0] + F[i][j][0]*dt*0.5) / rho[i][j];//add force term Guo al.
            u[i][j][1] = (momentum[i][j][1] + F[i][j][1]*dt*0.5) / rho[i][j];
        }
    }
}

// Function to calculate the phase field function
void calPhaseField() {
    const cglbm::lbm::ComponentPair components = componentPair();
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double sum_f = 0.0;
            double sum_g = 0.0;
            for (int k = 0; k < Q; k++) {
                sum_f += f[i][j][k];
                sum_g += g[i][j][k];
            }
            phi[i][j] = sum_g / sum_f;
            if (std::fabs(phi[i][j]) > 1.0 + 1.0e-6) {
                std::cout << "Error : phi = " << phi[i][j] << std::endl;
            }
            double rho_local = rho[i][j];
            double phi_local = sum_g / sum_f;
            // Two-component equation of state: each component keeps its own sound
            // speed and pressure at infinity, which is what decouples the density
            // ratio from the sound-speed ratio.
            //   T. Lafarge, P. Boivin, N. Odier, B. Cuenot, Phys. Fluids 33, 082110
            //   (2021), doi:10.1063/5.0061638 -- see src/lbm/equation_of_state.h.
            //
            // A linear mixing rule used to overwrite this value. On the shipped
            // Laplace case that cost 17 % of the pressure jump at a density ratio of
            // 10 and 28 % at 20; it is kept as pressure_linear_mixing() for
            // comparison, not used here.
            const double p_local = cglbm::lbm::pressure(rho_local, phi_local, components);
            p[i][j] = p_local;
            // The colour field the interface is located by. phi is a mass
            // fraction, and at a density ratio r its zero contour sits
            // (W/2) ln(r) outside the density interface -- 7.4 nodes at 1e4.
            // psi = alpha1 - alpha2 is centred on the density interface; see
            // src/lbm/mixture.h.
            psi[i][j] = cglbm::lbm::normalised_phase(phi_local, p_local, components);
            }
    }
}

// Function to perform the collision step
void collide() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double rho_local = rho[i][j];  // Local density
            double p_local = p[i][j];  // Local pressure
            // Calculate relaxation times, tau = mu / p + 1/2 with the mixture
            // viscosity rho (Y1 nu1 + Y2 nu2) = alpha1 mu1 + alpha2 mu2
            double nu_local   = cglbm::lbm::mixture_kinematic_viscosity(phi[i][j], nu1, nu2);
            double nu_b_local = cglbm::lbm::mixture_kinematic_viscosity(phi[i][j], nu_b1, nu_b2);
            double tau_nu = rho_local * nu_local   / (p_local * dt) + 0.5;  //shear relaxation time
            double tau_b  = rho_local * nu_b_local / (p_local * dt) + 0.5; //bulk relaxation time
            double sum_nu_neq = 0.0, sum_b_neq = 0.0, sum_xy_neq = 0.0;

            // Calculate $f_{k, i}^{r, neq}$ with $k\inn\{\nu, b, x y\}$
            for (int k = 0; k < Q; k++) {
                double xi_x = xi[k][0], xi_y = xi[k][1];
                double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
                double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - cs2;
                double H_xy = xi_x * xi_y;

                double f_neq = f[i][j][k] - f_eq[i][j][k] + 0.5*S[i][j][k]; // Calculate non-equilibrium part for each direction
                sum_nu_neq += f_neq * H_nu;
                sum_b_neq  += f_neq * H_b;
                sum_xy_neq += f_neq * H_xy;
            }
            // Calculate $\Omega_i^{(1)}$
            for (int k = 0; k < Q; k++) {
                double xi_x = xi[k][0], xi_y = xi[k][1];
                double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
                double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - cs2;
                double H_xy = xi_x * xi_y;

                double f_nu_neq = H_nu / cs4 * sum_nu_neq; //optimal way to calculate f_nu_neq without using vector, only using scalar
                double f_b_neq  = H_b  / cs4 * sum_b_neq;
                double f_xy_neq = H_xy / cs4 * sum_xy_neq;
                omega_1[i][j][k] = w[k]*(1.0 - 1.0 / tau_nu) * (f_nu_neq + f_xy_neq) + w[k]*(1.0 - 1.0 / tau_b) * f_b_neq;
            }
        }
    }
}

double S_F[Lx][Ly][Q]; // Force term for the potential volume forces eg.gravity
double S_Sp[Lx][Ly][Q]; // corrective term corrects the error stemming from the third order moment, improperly resolved on the D2Q9 lattice:
double S_t[Lx][Ly][Q]; // temporal correction term
// Function to calculate the force term
void force() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double u_x = u[i][j][0];
            double u_y = u[i][j][1];
            double F_x = F[i][j][0];
            double F_y = F[i][j][1];
            for (int k = 0; k < Q; k++) {
                double xik0 = xi[k][0];
                double xik1 = xi[k][1];
                double Hxxk = xik0*xik0 - cs2;
                double Hyyk = xik1*xik1 - cs2;
                double Hxyk = xik0*xik1;

                double term1 = (F_x * xik0 + F_y * xik1) / cs2;
                double term23 = (u_x * F_x * Hxxk + u_y * F_y * Hyyk + (u_x * F_y + u_y * F_x) * Hxyk) / cs4;

                S_F[i][j][k] = w[k] * (term1 + term23);
            }
            double derive_x = 0., derive_y = 0.;
            for (int k = 0; k < Q; k++) {
                int xi_x = (int)xi[k][0];
                int xi_y = (int)xi[k][1];
                int ind_x = (i + xi_x + Lx) % Lx; // periodic boundary conditions
                int ind_y = (j + xi_y + Ly) % Ly; // periodic boundary conditions
                derive_x += w[k] * xi_x * (p[ind_x][ind_y] - rho[ind_x][ind_y] * cs2) * u[ind_x][ind_y][0];
                derive_y += w[k] * xi_y * (p[ind_x][ind_y] - rho[ind_x][ind_y] * cs2) * u[ind_x][ind_y][1];
            }
            derive_x /= (dt*cs2);
            derive_y /= (dt*cs2);
            for (int k = 0; k < Q; k++) {
                double H_nu = (xi[k][0] * xi[k][0] - xi[k][1] * xi[k][1]) / 2.;
                double H_b  = (xi[k][0] * xi[k][0] + xi[k][1] * xi[k][1]) / 2. -  cs2 ;
                S_Sp[i][j][k] = w[k] * (derive_y*(3.*H_nu-H_b) + derive_x*(-3.*H_nu-H_b)) / (2.*cs4);
            }
            for (int k = 0; k < Q; k++) {
                double H_xx = xi[k][0] * xi[k][0] - cs2;
                double H_yy = xi[k][1] * xi[k][1] - cs2;
                double H_xxyy = xi[k][0] * xi[k][0] * xi[k][1] * xi[k][1] - cs2 * (xi[k][0] * xi[k][0] + xi[k][1] * xi[k][1]) + cs4;
                double E = w[k] * ((H_xx+H_yy) / (2.*cs4) - H_xxyy / (4.*cs6)) ; //w[k] * ((H_xx+H_yy) / (2.*cs4) - H_xxyy / (4.*cs6));
                S_t[i][j][k] = (p[i][j] - p_mdt[i][j] - (rho[i][j] - rho_mdt[i][j]) * cs2) * E;
            }

            for (int k = 0; k < Q; k++) {
                S[i][j][k] = S_F[i][j][k] + S_Sp[i][j][k] + S_t[i][j][k];
            }
        }
    }
}

// Surface tension, as the divergence of the capillary stress
// T = sigma/2 (|grad psi| I - grad psi grad psi / |grad psi|), which enters
// through S_F like any body force (src/lbm/surface_force.h).
//
// This replaces the perturbation operator Omega^(2), which wrote the same
// stress into the non-equilibrium populations scaled by 1/tau. Streaming then
// mixes it between neighbours of different tau before it is relaxed, and
// across an interface with a viscosity contrast the jump came out wrong:
// 0.71 sigma/R on this case at density ratio 20, where tau changes twentyfold
// across the droplet boundary. The force does not depend on tau, and being a
// divergence it conserves momentum exactly.
//
//   B. Lafaurie, C. Nardone, R. Scardovelli, S. Zaleski, G. Zanetti,
//   J. Comput. Phys. 113, 134 (1994)
void calSurfaceForce() {
    // the stress first: the force differentiates it
    for (int i=0 ; i<Lx ; i++){
        for (int j=0 ; j<Ly ; j++){
            // Colour gradient, already carrying the 1/c_s^2 factor.
            double Cx = 0.0, Cy = 0.0;
            cglbm::lbm::gradient_periodic(&psi[0][0], Lx, Ly, i, j, gradient_stencil, &Cx, &Cy);
            cglbm::lbm::unit_normal(Cx, Cy, &normal_x[i][j], &normal_y[i][j]);
            cglbm::lbm::capillary_stress(sigma, Cx, Cy, &stress_xx[i][j], &stress_xy[i][j], &stress_yy[i][j]);
        }
    }
    for (int i=0 ; i<Lx ; i++){
        for (int j=0 ; j<Ly ; j++){
            cglbm::lbm::surface_force(&stress_xx[0][0], &stress_xy[0][0], &stress_yy[0][0], Lx, Ly, i, j,
                                      gradient_stencil, cglbm::lbm::Boundary::Periodic, &F[i][j][0], &F[i][j][1]);
        }
    }
}

void recolor(){
    for (int i=0 ; i<Lx ; i++){
        for (int j=0 ; j<Ly ; j++){
            // Direction of the colour gradient, from calSurfaceForce. It is taken
            // on psi, which has the same direction as grad(phi) but is centred
            // on the density interface. The amplitude still uses phi: the tanh
            // profile it maintains in the mass fraction is the same tanh profile
            // in the volume fraction, shifted (src/lbm/mixture.h).
            double n_x = normal_x[i][j];
            double n_y = normal_y[i][j];
            if (n_x != 0.0 || n_y != 0.0) { // no interface where the gradient vanishes
                for (int k = 0; k < Q; k++) {
                        double xi_x = xi[k][0],  xi_y = xi[k][1];
                        omega_3[i][j][k] = w[k] * p[i][j] * (1 - phi[i][j] * phi[i][j]) / (2. * ch_width_ope) * (xi_x * n_x + xi_y * n_y) / cs2;
                    }
            }
            else {
                for (int k = 0; k < Q; k++) {
                    omega_3[i][j][k] = 0.0;
                }
            }
        }
    }
}



// Function to perform streaming step
void stream() {
    for (int i=0 ; i<Lx ; i++){
        // std::cout << "i = " << i << std::endl;
        for (int j=0 ; j<Ly ; j++){
            // std::cout << "j = " << j << std::endl;
            for (int k=0 ; k<Q ; k++){
                int ip = (i + (int)xi[k][0] + Lx) % Lx;
                int jp = (j + (int)xi[k][1] + Ly) % Ly;
                f[ip][jp][k] = f_eq[i][j][k] + omega_1[i][j][k] + 0.5*S[i][j][k]; // surface tension is in S, via F
                g[ip][jp][k] = f[ip][jp][k] * phi[i][j] + omega_3[i][j][k];
            }
            // Save the macroscopic variables at the previous time step for the force step
            rho_mdt[i][j] = rho[i][j];
            p_mdt[i][j] = p[i][j];
        }
    }
}

#include <fstream>
//#include <iostream> // Already included in the head
#include <string>
// Function to output data in VTK format
void outputVTK(const std::string& filename, int timestep) {
    std::ofstream vtkfile;
    std::string fullFilename = filename + std::to_string(timestep) + ".vtk";
    vtkfile.open(fullFilename);

    // VTK file header and data structure
    vtkfile << "# vtk DataFile Version 4.0" << std::endl;
    vtkfile << "Lattice Boltzmann Method Data" << std::endl;
    vtkfile << "ASCII" << std::endl;
    vtkfile << "DATASET STRUCTURED_POINTS" << std::endl;
    vtkfile << "DIMENSIONS " << Lx << " " << Ly << " 1" << std::endl;
    vtkfile << "ORIGIN 0 0 0" << std::endl;
    vtkfile << "SPACING 1 1 1" << std::endl;

    // Output density
    vtkfile << "POINT_DATA " << Lx * Ly << std::endl;
    vtkfile << "SCALARS density float" << std::endl;
    vtkfile << "LOOKUP_TABLE default" << std::endl;
    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            vtkfile << rho[i][j] << std::endl;
        }
    }

    // Output velocity
    vtkfile << "VECTORS velocity float" << std::endl;
    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            vtkfile << u[i][j][0] << " " << u[i][j][1] << " 0.0" << std::endl;
        }
    }

    // Output phase field
    vtkfile << "SCALARS phase_field float" << std::endl;
    vtkfile << "LOOKUP_TABLE default" << std::endl;
    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            vtkfile << phi[i][j] << std::endl;
        }
    }

    //output pressure
    vtkfile << "SCALARS pressure float" << std::endl;
    vtkfile << "LOOKUP_TABLE default" << std::endl;
    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            vtkfile << p[i][j] << std::endl;
        }
    }

    vtkfile.close();
}

// #include <iomanip>  // For controlling float print precision
// void outputDataCSV(const std::string& baseName, int timestep) {
void outputDataCSV(int timestep) {
    std::ofstream fileDensity, fileVelocity, filePhase, filePressure;
    std::string densityFilename = "density_" + std::to_string(timestep) + ".csv";
    std::string velocityFilename = "velocity_" + std::to_string(timestep) + ".csv";
    std::string phaseFilename = "phase_" + std::to_string(timestep) + ".csv";
    std::string pressureFilename = "pressure_" + std::to_string(timestep) + ".csv";

    fileDensity.open(densityFilename);
    fileVelocity.open(velocityFilename);
    filePhase.open(phaseFilename);
    filePressure.open(pressureFilename);

    // Set precision for float values
    // fileDensity << std::fixed << std::setprecision(8);
    // fileVelocity << std::fixed << std::setprecision(8);
    // filePhase << std::fixed << std::setprecision(8);
    // filePressure << std::fixed << std::setprecision(8);

    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            fileDensity << rho[i][j];
            fileVelocity << u[i][j][0] << "," << u[i][j][1];
            filePhase << phi[i][j];
            filePressure << p[i][j];
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

    fileDensity.close();
    fileVelocity.close();
    filePhase.close();
    filePressure.close();
}

// Function to initialize the simulation for ellipsoidal droplet

void initialize() {
    // Initialize distribution functions, density, velocity, phase field function, pressure at infinity, etc.
    // Define the center of the bubble/droplet
    int x0 = Lx / 2;
    int y0 = Ly / 2;
    // Define the radius of the droplet/bubble and the interface width
    double r =  radius * dx ; // Lx / 8. * dx; // 16 lattice units

    // In order to get f_eq, we need to calculate the macroscopic variables rho(with rho1, rho2 at different nodes), u, p
    //
    // The tanh profile is given to the volume fraction, the pressure steps from
    // p2 outside to p2 + sigma/R inside across it, and each component takes the
    // density its own branch gives at that pressure (src/lbm/mixture.h). The
    // state is then in equilibrium with the equation of state node by node.
    // Setting rho linear in phi instead, as this case used to, leaves the
    // interface far from it -- at density ratio 1000 its pressure starts about
    // 300 times the ambient one and the run blows up within ten steps.
    const cglbm::lbm::ComponentPair components = componentPair();
    const double p_outside = rho2 * c2 * c2 - p2_inf; // ambient pressure, in component 2
    for (int i=0; i<Lx ; i++){
        for (int j=0; j<Ly ; j++){
            double distance = sqrt((i - x0) * (i - x0) + (j - y0) * (j - y0));
            double alpha1 = 0.5 * (1.0 - tanh((distance - r) / ch_width_init)); // volume fraction of the droplet
            double p_local = p_outside + sigma / radius * alpha1;
            const cglbm::lbm::MixtureState state = cglbm::lbm::mixture_from_volume_fraction(alpha1, p_local, components);
            double phi_local = state.phi;
            phi[i][j] = phi_local;  // Local phase field
            u[i][j][0] = 0.0; // static flow field
            u[i][j][1] = 0.0; // static flow field
            double rho_local = state.rho;
            rho[i][j] = rho_local;  // Loc-al density
            rho_mdt[i][j] = rho_local; // At zero time step, the value of previous step is the one of current step

            // Pressure from the equation of state; it returns p_local to rounding.
            p_local = cglbm::lbm::pressure(rho_local, phi_local, components);
            p[i][j] = p_local;
            p_mdt[i][j] = p_local;
            psi[i][j] = cglbm::lbm::normalised_phase(phi_local, p_local, components);
        }
    }
    calSurfaceForce();
    calEquilibrium();
    for (int i=0 ; i<Lx ; i++){
        for (int j=0 ; j<Ly ; j++){
            for (int k=0 ; k<Q ; k++){
                f[i][j][k] = f_eq[i][j][k];
                g[i][j][k] = f[i][j][k] * phi[i][j];
            }
        }
    }
}

const double damping = 0.001;
const int numBoundary = 1;
//Damp the acoustic wave on the boundary
void applyAbsorbingBoundary() {
    for (int i = 2; i < Lx; i++) {
        for (int j = 2; j < Ly; j++) {
            // if ( i<numBoundary&j<numBoundary || i>Lx-numBoundary&j>Ly-numBoundary || i<numBoundary&j>Ly-numBoundary || i>Lx-numBoundary&j<numBoundary) {
            if ( i<numBoundary || i>Lx-numBoundary|| j>Ly-numBoundary || j<numBoundary) {
                for (int k=0 ; k<Q ; k++){
                f[i][j][k] *= damping; // Apply damping
                g[i][j][k] *= damping; // Apply damping
                }
            }
        }
    }
}

// Function to run the simulation
void runSimulation() {
    // double r3 = pow((sqrt(a2) + sqrt(b2)) / 2.* radius, 3);
    std::cout << "sigma = " << sigma << std::endl;
    std::cout << "radius = " << radius << std::endl;
    std::cout << "p1_inf = " << p1_inf << std::endl;
    std::cout << "p2_inf = " << p2_inf << std::endl;
    std::cout << "rho1 = " << rho1 << std::endl;
    std::cout << "rho2 = " << rho2 << std::endl;
    std::cout << "nu1 = " << nu1 << std::endl;
    std::cout << "nu2 = " << nu2 << std::endl;
    // std::cout << "T_theo = " << 2*3.1415 * sqrt((rho1+rho2)*r3/6./sigma) << std::endl; //T_{\text {theo }}=2 \pi \sqrt{\frac{\left(\rho_1+\rho_2\right) r^3}{6 \sigma}}
    initialize();
    //outputVTK("lbm_output_", 0);
    outputDataCSV(0);
    for (int i = 1; i < numSteps+1; i++) {
        // std::cout << "Step " << i << std::endl;
        force();
        collide();
        recolor();
        stream();
        // if (i < 2){
        //     applyAbsorbingBoundary();
        // }

        calMacroscopic();
        calPhaseField();
        calSurfaceForce();
        calVelocity();
        calEquilibrium();

        //Output or visualization code here
        if (i % interval == 0) {  // Output every 100 steps
            //outputVTK("lbm_output_", i);
            std::cout << "Step " << i << std::endl;
            outputDataCSV(i);
        }
        // outputVTK("lbm_output_", i);
        // outputDataCSV(i);
    }
}

// A strictly positive number, or false.
bool parsePositive(const char* text, double* value) {
    char* end = nullptr;
    const double parsed = std::strtod(text, &end);
    if (end == text || *end != '\0' || !(parsed > 0.0) || !std::isfinite(parsed)) {
        return false;
    }
    *value = parsed;
    return true;
}

int main(int argc, char** argv) {
    // Optional first argument selects the colour-gradient stencil: E4, E6, E8.
    if (argc > 1 && !cglbm::lbm::stencil_from_name(argv[1], &gradient_stencil)) {
        std::cerr << "Unknown gradient stencil '" << argv[1]
                  << "'; expected E4, E6 or E8." << std::endl;
        return 2;
    }
    // Optional second and third: the density ratio rho1/rho2 and the dynamic
    // viscosity ratio mu1/mu2, which defaults to the density ratio.
    double density_ratio = rho1 / rho2;
    if (argc > 2 && !parsePositive(argv[2], &density_ratio)) {
        std::cerr << "Invalid density ratio '" << argv[2] << "'; expected a positive number." << std::endl;
        return 2;
    }
    double viscosity_ratio = density_ratio;
    if (argc > 3 && !parsePositive(argv[3], &viscosity_ratio)) {
        std::cerr << "Invalid viscosity ratio '" << argv[3] << "'; expected a positive number." << std::endl;
        return 2;
    }
    setComponents(density_ratio, viscosity_ratio);
    std::cout << "colour gradient stencil = " << cglbm::lbm::stencil_name(gradient_stencil)
              << std::endl;
    std::cout << "density ratio = " << density_ratio << std::endl;
    std::cout << "viscosity ratio = " << viscosity_ratio << std::endl;
    runSimulation();
    return 0;
}
