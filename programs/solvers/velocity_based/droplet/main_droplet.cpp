#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "lbm/isotropic_gradient.h"
#include "lbm/surface_force.h"
#include "lbm/velocity_based.h"

// A droplet at a large density ratio, with the velocity-based scheme of
// src/lbm/velocity_based.h: static (the Laplace law) or launched through a
// quiescent lighter fluid (an interface that moves).
//
// Usage: droplet [E4|E6|E8] [density_ratio] [velocity] [viscosity_ratio]
//
//   density_ratio    rho1/rho2, default 1e4.
//   velocity         initial speed of the droplet along x, lattice units per
//                    step; the surrounding fluid starts at rest. Default 0: a
//                    static droplet, the Laplace benchmark.
//   viscosity_ratio  mu1/mu2, default 1.
//
// The colour-gradient solvers stream the density and the momentum, which jump
// by the density ratio across the interface; at 1e4 their droplets diverge as
// soon as they move at 1e-3. Here the streamed moments are p/(rho cs^2) and u,
// both continuous, and the density follows from a bounded volume fraction.
// The momentum is updated link by link, with equal and opposite exchanges, so
// that it is conserved to rounding; see src/lbm/velocity_based.h and
// docs/numerics.md.
//
// Periodic in both directions.

namespace vb = cglbm::lbm::velocity_based;

const int Lx = 128; // Number of lattice nodes in the x-direction
const int Ly = 128; // Number of lattice nodes in the y-direction
const int Q = vb::kQ; // Number of discrete velocities

const double dx = 1.; // Lattice spacing
const double dt = 1.; // Time step

const double c_dx = 1.e-5 ; // m : conversion factor from lattice units to physical units
const double c_dt = c_dx/347./sqrt(3.); // s : conversion factor from lattice units to physical units

const int numSteps = 10000; // Number of simulation steps
const int interval = 1000; // Output interval

const double cs2 = vb::kSoundSpeedSquared;

//parameters
double rho1 = 1.e4; // kg * m-3 Density of the droplet, overridable from the command line
const double rho2 = 1.;  // kg * m-3 Density of the surrounding fluid
const double radius = 10.; // Radius of the droplet
const double sigma = 1./(c_dx * c_dx * c_dx / c_dt / c_dt); // surface tension, as in the laplace case
const double width = 1.6 * dx; // Interface width W: c = (1 + tanh(x/W))/2
double velocity = 0.; // Initial velocity of the droplet, lattice units

// viscosities: the surrounding fluid gets a tenth of the laplace case's, which
// puts its relaxation time at 1; the droplet's follows from the ratio
const double nu2 = 1.e-3/(c_dx * c_dx / c_dt); // kinematic viscosity of component 2, lattice units
const double mu2 = rho2 * nu2; // dynamic viscosity of component 2
double mu1 = mu2; // dynamic viscosity of component 1
// The trace of the non-equilibrium relaxes at its own rate. A heavy fluid with
// a modest dynamic viscosity has a shear tau within 1e-3 of 1/2, and without
// this its acoustic modes are not damped.
const double tau_bulk = 1.0;
// Weight of the populations' own non-equilibrium in the hybrid regularised
// collision; the rest is its finite-difference estimate. It damps the grid-scale
// modes of a heavy fluid whose shear tau is within 1e-4 of 1/2.
const double hybrid_weight = 0.98;

double g[Lx][Ly][Q]; // Hydrodynamic populations: moments P and u
double h[Lx][Ly][Q]; // Phase populations: moment c
double g_new[Lx][Ly][Q];
double h_new[Lx][Ly][Q];

double c[Lx][Ly]; // Volume fraction of component 1
double rho[Lx][Ly]; // Density rho2 + c (rho1 - rho2)
double psi[Lx][Ly]; // 2c - 1, the field the surface tension is built on
double P[Lx][Ly]; // Pressure number p / (rho cs^2); p is the pressure relative to the far field
double u[Lx][Ly][2]; // Velocity
double u_x[Lx][Ly]; // Velocity components, as scalar fields for the gradient stencils
double u_y[Lx][Ly];
double a[Lx][Ly][2]; // Acceleration: surface tension and pressure
double u_lattice[Lx][Ly][2]; // Velocity from the momentum exchange, before the half acceleration
// The previous step, which the momentum exchange reads
double rho_old[Lx][Ly];
double P_old[Lx][Ly];
double u_old[Lx][Ly][2];
double a_old[Lx][Ly][2];
double normal_x[Lx][Ly]; // Unit interface normal grad(psi)/|grad(psi)|
double normal_y[Lx][Ly];
double stress_xx[Lx][Ly]; // Capillary stress sigma/2 (|grad psi| I - grad psi grad psi / |grad psi|)
double stress_xy[Lx][Ly];
double stress_yy[Lx][Ly];

// Isotropy order of the gradients (surface tension, normals).
cglbm::lbm::GradientStencil gradient_stencil = cglbm::lbm::GradientStencil::E8;

double dynamicViscosity(double c_local) {
    const double bounded = (c_local < 0.0) ? 0.0 : ((c_local > 1.0) ? 1.0 : c_local);
    return bounded * mu1 + (1.0 - bounded) * mu2;
}

// c from the phase populations, then everything that follows from it and P.
void calMacroscopic() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double sum_h = 0.0;
            double sum_g = 0.0;
            for (int k = 0; k < Q; k++) {
                sum_h += h[i][j][k];
                sum_g += g[i][j][k];
            }
            c[i][j] = sum_h;
            const double bounded = (sum_h < 0.0) ? 0.0 : ((sum_h > 1.0) ? 1.0 : sum_h);
            rho[i][j] = rho2 + bounded * (rho1 - rho2);
            psi[i][j] = 2.0 * bounded - 1.0;
            P[i][j] = sum_g;
        }
    }
}

// Surface force div(T) and pressure force, both per unit mass.
void calAcceleration() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double gx = 0.0, gy = 0.0;
            cglbm::lbm::gradient_periodic(&psi[0][0], Lx, Ly, i, j, gradient_stencil, &gx, &gy);
            cglbm::lbm::unit_normal(gx, gy, &normal_x[i][j], &normal_y[i][j]);
            cglbm::lbm::capillary_stress(sigma, gx, gy, &stress_xx[i][j], &stress_xy[i][j], &stress_yy[i][j]);
        }
    }
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double fx = 0.0, fy = 0.0;
            cglbm::lbm::surface_force(&stress_xx[0][0], &stress_xy[0][0], &stress_yy[0][0], Lx, Ly, i, j,
                                      gradient_stencil, cglbm::lbm::Boundary::Periodic, &fx, &fy);
            double px = 0.0, py = 0.0;
            vb::pressure_force(&P[0][0], &rho[0][0], Lx, Ly, i, j, &px, &py);
            a[i][j][0] = fx / rho[i][j] + px;
            a[i][j][1] = fy / rho[i][j] + py;
        }
    }
}

// The momentum after streaming: what each node held after its collision, plus
// the exchanges over its eight links. The populations have just streamed, so
// g[x][k] is what the neighbour x - xi_k sent along k, and g[x - xi_k][opp k]
// what x sent back.
void calMomentum() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            // post-collision velocity u + a/2 of the previous step
            double jx = rho_old[i][j] * (u_old[i][j][0] + 0.5 * a_old[i][j][0]);
            double jy = rho_old[i][j] * (u_old[i][j][1] + 0.5 * a_old[i][j][1]);
            for (int k = 1; k < Q; k++) {
                const int id = (i - vb::kVelocity[k][0] + Lx) % Lx;
                const int jd = (j - vb::kVelocity[k][1] + Ly) % Ly;
                const int back = vb::kOpposite[k];
                vb::LinkEnd donor;
                donor.outgoing = g[i][j][k] - vb::kWeight[k] * P_old[id][jd];
                donor.phase = h[i][j][k];
                donor.rho = rho[id][jd];
                donor.mu = dynamicViscosity(c[id][jd]);
                donor.ux = u_old[id][jd][0];
                donor.uy = u_old[id][jd][1];
                vb::LinkEnd receiver;
                receiver.outgoing = g[id][jd][back] - vb::kWeight[k] * P_old[i][j];
                receiver.phase = h[id][jd][back];
                receiver.rho = rho[i][j];
                receiver.mu = dynamicViscosity(c[i][j]);
                receiver.ux = u_old[i][j][0];
                receiver.uy = u_old[i][j][1];
                double link_x = 0.0, link_y = 0.0;
                vb::link_momentum(k, donor, receiver, rho1, rho2, &link_x, &link_y);
                jx += link_x;
                jy += link_y;
            }
            u_lattice[i][j][0] = jx / rho[i][j];
            u_lattice[i][j][1] = jy / rho[i][j];
        }
    }
}

// The populations take the exchanged velocity; the macroscopic one adds the
// half acceleration of the forcing scheme.
void calVelocity() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            vb::set_velocity(g[i][j], u_lattice[i][j][0], u_lattice[i][j][1]);
            u[i][j][0] = u_lattice[i][j][0] + 0.5 * a[i][j][0] * dt;
            u[i][j][1] = u_lattice[i][j][1] + 0.5 * a[i][j][1] * dt;
            u_x[i][j] = u[i][j][0];
            u_y[i][j] = u[i][j][1];
        }
    }
}

// Keep the fields the next momentum exchange reads.
void saveStep() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            rho_old[i][j] = rho[i][j];
            P_old[i][j] = P[i][j];
            for (int d = 0; d < 2; d++) {
                u_old[i][j][d] = u[i][j][d];
                a_old[i][j][d] = a[i][j][d];
            }
        }
    }
}

// Collision of g, construction of h, and streaming of both.
void collideAndStream() {
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            double eq[Q], source[Q], post[Q], phase[Q];
            vb::hydrodynamic_equilibrium(P[i][j], u[i][j][0], u[i][j][1], eq);
            vb::forcing(u[i][j][0], u[i][j][1], a[i][j][0], a[i][j][1], source);
            const double tau = dynamicViscosity(c[i][j]) / (rho[i][j] * cs2 * dt) + 0.5; // shear relaxation time
            vb::VelocityGradient gradient;
            cglbm::lbm::gradient_periodic(&u_x[0][0], Lx, Ly, i, j, cglbm::lbm::GradientStencil::E4,
                                          &gradient.dux_dx, &gradient.dux_dy);
            cglbm::lbm::gradient_periodic(&u_y[0][0], Lx, Ly, i, j, cglbm::lbm::GradientStencil::E4,
                                          &gradient.duy_dx, &gradient.duy_dy);
            vb::collide_hybrid(g[i][j], eq, source, tau, tau_bulk, hybrid_weight, gradient, post);
            vb::phase_populations(c[i][j], u[i][j][0], u[i][j][1], normal_x[i][j], normal_y[i][j], width, phase);
            for (int k = 0; k < Q; k++) {
                const int ip = (i + vb::kVelocity[k][0] + Lx) % Lx;
                const int jp = (j + vb::kVelocity[k][1] + Ly) % Ly;
                g_new[ip][jp][k] = post[k];
                h_new[ip][jp][k] = phase[k];
            }
        }
    }
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            for (int k = 0; k < Q; k++) {
                g[i][j][k] = g_new[i][j][k];
                h[i][j][k] = h_new[i][j][k];
            }
        }
    }
}

void initialize() {
    const double x0 = Lx / 2;
    const double y0 = Ly / 2;
    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            const double distance = sqrt((i - x0) * (i - x0) + (j - y0) * (j - y0));
            const double c_local = 0.5 * (1.0 - tanh((distance - radius) / width));
            // the droplet moves, its surroundings are at rest; the Laplace jump
            // is in place from the start
            const double p_local = sigma / radius * c_local;
            const double rho_local = rho2 + c_local * (rho1 - rho2);
            u[i][j][0] = velocity * c_local;
            u[i][j][1] = 0.0;
            u_x[i][j] = u[i][j][0];
            u_y[i][j] = 0.0;
            double eq[Q], gamma[Q];
            vb::hydrodynamic_equilibrium(p_local / (rho_local * cs2), u[i][j][0], 0.0, eq);
            vb::velocity_equilibrium(u[i][j][0], 0.0, gamma);
            for (int k = 0; k < Q; k++) {
                g[i][j][k] = eq[k];
                h[i][j][k] = c_local * gamma[k];
            }
        }
    }
    calMacroscopic();
    calAcceleration();
}

// Four ASCII grids per output, as the colour-gradient solvers write them. The
// phase field is psi = 2c - 1, +1 in the droplet; the pressure is relative to
// the far field.
void outputDataCSV(int timestep) {
    std::ofstream fileDensity("density_" + std::to_string(timestep) + ".csv");
    std::ofstream fileVelocity("velocity_" + std::to_string(timestep) + ".csv");
    std::ofstream filePhase("phase_" + std::to_string(timestep) + ".csv");
    std::ofstream filePressure("pressure_" + std::to_string(timestep) + ".csv");
    fileDensity.precision(10);
    fileVelocity.precision(10);
    filePhase.precision(10);
    filePressure.precision(10);
    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            fileDensity << rho[i][j];
            fileVelocity << u[i][j][0] << "," << u[i][j][1];
            filePhase << psi[i][j];
            filePressure << P[i][j] * rho[i][j] * cs2;
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

void runSimulation() {
    std::cout << "sigma = " << sigma << std::endl;
    std::cout << "radius = " << radius << std::endl;
    std::cout << "rho1 = " << rho1 << std::endl;
    std::cout << "rho2 = " << rho2 << std::endl;
    std::cout << "mu1 = " << mu1 << std::endl;
    std::cout << "mu2 = " << mu2 << std::endl;
    std::cout << "velocity = " << velocity << std::endl;
    initialize();
    outputDataCSV(0);
    for (int n = 1; n < numSteps + 1; n++) {
        saveStep();
        collideAndStream();
        calMacroscopic();
        calMomentum();
        calAcceleration();
        calVelocity();
        if (n % interval == 0) {
            std::cout << "Step " << n << std::endl;
            outputDataCSV(n);
        }
    }
}

// A number, or false.
bool parseNumber(const char* text, double* value) {
    char* end = nullptr;
    const double parsed = std::strtod(text, &end);
    if (end == text || *end != '\0' || !std::isfinite(parsed)) {
        return false;
    }
    *value = parsed;
    return true;
}

int main(int argc, char** argv) {
    if (argc > 1 && !cglbm::lbm::stencil_from_name(argv[1], &gradient_stencil)) {
        std::cerr << "Unknown gradient stencil '" << argv[1] << "'; expected E4, E6 or E8." << std::endl;
        return 2;
    }
    double density_ratio = rho1 / rho2;
    if (argc > 2 && !(parseNumber(argv[2], &density_ratio) && density_ratio > 0.0)) {
        std::cerr << "Invalid density ratio '" << argv[2] << "'; expected a positive number." << std::endl;
        return 2;
    }
    if (argc > 3 && !parseNumber(argv[3], &velocity)) {
        std::cerr << "Invalid velocity '" << argv[3] << "'; expected a number." << std::endl;
        return 2;
    }
    double viscosity_ratio = 1.0;
    if (argc > 4 && !(parseNumber(argv[4], &viscosity_ratio) && viscosity_ratio > 0.0)) {
        std::cerr << "Invalid viscosity ratio '" << argv[4] << "'; expected a positive number." << std::endl;
        return 2;
    }
    rho1 = density_ratio * rho2;
    mu1 = viscosity_ratio * mu2;
    std::cout << "gradient stencil = " << cglbm::lbm::stencil_name(gradient_stencil) << std::endl;
    std::cout << "density ratio = " << density_ratio << std::endl;
    std::cout << "viscosity ratio = " << viscosity_ratio << std::endl;
    runSimulation();
    return 0;
}
