#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "lbm/isotropic_gradient.h"
#include "lbm/velocity_based_solver.h"

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

const int Lx = 128;  // Number of lattice nodes in the x-direction
const int Ly = 128;  // Number of lattice nodes in the y-direction

const double dx = 1.;  // Lattice spacing

const double c_dx = 1.e-5;  // m : conversion factor from lattice units to physical units
const double c_dt =
    c_dx / 347. / sqrt(3.);  // s : conversion factor from lattice units to physical units

const int numSteps = 10000;  // Number of simulation steps
const int interval = 1000;   // Output interval

// parameters
double rho1 = 1.e4;         // kg * m-3 Density of the droplet, overridable from the command line
const double rho2 = 1.;     // kg * m-3 Density of the surrounding fluid
const double radius = 10.;  // Radius of the droplet
const double sigma =
    1. / (c_dx * c_dx * c_dx / c_dt / c_dt);  // surface tension, as in the laplace case
const double width = 1.6 * dx;                // Interface width W: c = (1 + tanh(x/W))/2
double velocity = 0.;                         // Initial velocity of the droplet, lattice units

// viscosities: the surrounding fluid gets a tenth of the laplace case's, which
// puts its relaxation time at 1; the droplet's follows from the ratio
const double nu2 =
    1.e-3 / (c_dx * c_dx / c_dt);  // kinematic viscosity of component 2, lattice units
const double mu2 = rho2 * nu2;     // dynamic viscosity of component 2
double mu1 = mu2;                  // dynamic viscosity of component 1

// Isotropy order of the gradients (surface tension, normals).
cglbm::lbm::GradientStencil gradient_stencil = cglbm::lbm::GradientStencil::E8;

// Four ASCII grids per output, as the colour-gradient solvers write them. The
// phase field is psi = 2c - 1, +1 in the droplet; the pressure is relative to
// the far field.
void outputDataCSV(const vb::Solver& solver, int timestep) {
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
            fileDensity << solver.density(i, j);
            fileVelocity << solver.velocity_x(i, j) << "," << solver.velocity_y(i, j);
            filePhase << solver.phase(i, j);
            filePressure << solver.pressure(i, j);
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

    vb::SolverParameters parameters;
    parameters.nx = Lx;
    parameters.ny = Ly;
    parameters.rho1 = rho1;
    parameters.rho2 = rho2;
    parameters.mu1 = mu1;
    parameters.mu2 = mu2;
    parameters.surface_tension = sigma;
    parameters.width = width;
    parameters.stencil = gradient_stencil;
    vb::Solver solver(parameters);

    const double x0 = Lx / 2;
    const double y0 = Ly / 2;
    solver.initialize([&](int i, int j) {
        const double distance = sqrt((i - x0) * (i - x0) + (j - y0) * (j - y0));
        vb::NodeState state;
        state.c = 0.5 * (1.0 - tanh((distance - radius) / width));
        // the droplet moves, its surroundings are at rest; the Laplace jump
        // is in place from the start
        state.ux = velocity * state.c;
        state.p = sigma / radius * state.c;
        return state;
    });
    outputDataCSV(solver, 0);
    for (int n = 1; n < numSteps + 1; n++) {
        solver.step();
        if (n % interval == 0) {
            std::cout << "Step " << n << std::endl;
            outputDataCSV(solver, n);
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
        std::cerr << "Unknown gradient stencil '" << argv[1] << "'; expected E4, E6 or E8."
                  << std::endl;
        return 2;
    }
    double density_ratio = rho1 / rho2;
    if (argc > 2 && !(parseNumber(argv[2], &density_ratio) && density_ratio > 0.0)) {
        std::cerr << "Invalid density ratio '" << argv[2] << "'; expected a positive number."
                  << std::endl;
        return 2;
    }
    if (argc > 3 && !parseNumber(argv[3], &velocity)) {
        std::cerr << "Invalid velocity '" << argv[3] << "'; expected a number." << std::endl;
        return 2;
    }
    double viscosity_ratio = 1.0;
    if (argc > 4 && !(parseNumber(argv[4], &viscosity_ratio) && viscosity_ratio > 0.0)) {
        std::cerr << "Invalid viscosity ratio '" << argv[4] << "'; expected a positive number."
                  << std::endl;
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
