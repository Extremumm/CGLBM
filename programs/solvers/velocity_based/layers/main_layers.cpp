#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "lbm/isotropic_gradient.h"
#include "lbm/velocity_based_solver.h"

// Two layers at a large density ratio sheared across their interfaces, with
// the velocity-based scheme of src/lbm/velocity_based.h: a Kolmogorov flow.
//
// Usage: layers [E4|E6|E8] [density_ratio] [amplitude] [viscosity_ratio]
//
//   density_ratio    rho1/rho2, default 1e4.
//   amplitude        peak velocity U of the lighter layer, lattice units per
//                    step. Default 0.01.
//   viscosity_ratio  mu1/mu2, default 1.
//
// Component 1 fills 0 < y < Ly/2, component 2 the rest; both boundaries are
// periodic. A body force G = U mu2 k^2 sin(k y), k = 2 pi / Ly, along x drives
//
//     u = U (mu2 / mu) sin(k y),
//
// the steady flow of a layered fluid, whose shear stress (G / k) cos(k y) is
// largest on the two interfaces, where the velocity vanishes. Both layers
// start in it. The heavy layer, whose viscous time Ly^2 / nu1 is some 10^8
// steps at 1e4, is held by its inertia; the lighter one, some 10^4 steps,
// stays there only if the stress crosses the interface intact, and any slip
// shows as an offset of the whole light layer.
//
// The flow is uniform along x, which is 8 nodes wide.

namespace vb = cglbm::lbm::velocity_based;

const int Lx = 8; // Number of lattice nodes in the x-direction
const int Ly = 128; // Number of lattice nodes in the y-direction

const double dx = 1.; // Lattice spacing

const double c_dx = 1.e-5 ; // m : conversion factor from lattice units to physical units
const double c_dt = c_dx/347./sqrt(3.); // s : conversion factor from lattice units to physical units

const int numSteps = 10000; // Number of simulation steps
const int interval = 1000; // Output interval

//parameters, as in the droplet case
double rho1 = 1.e4; // kg * m-3 Density of the lower layer, overridable from the command line
const double rho2 = 1.;  // kg * m-3 Density of the upper layer
const double sigma = 1./(c_dx * c_dx * c_dx / c_dt / c_dt); // surface tension; flat interfaces feel none
const double width = 1.6 * dx; // Interface width W: c = (1 + tanh(x/W))/2
double amplitude = 0.01; // Peak velocity of the lighter layer, lattice units

const double nu2 = 1.e-3/(c_dx * c_dx / c_dt); // kinematic viscosity of component 2, lattice units
const double mu2 = rho2 * nu2; // dynamic viscosity of component 2
double mu1 = mu2; // dynamic viscosity of component 1

// Isotropy order of the gradients (surface tension, normals).
cglbm::lbm::GradientStencil gradient_stencil = cglbm::lbm::GradientStencil::E8;

// Four ASCII grids per output, as the droplet case writes them.
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

// Volume fraction of component 1 across the layers.
double volumeFraction(int j) {
    double distance = std::fabs(j - Ly / 4.0);
    distance = std::fmin(distance, Ly - distance); // periodic
    return 0.5 * (1.0 - tanh((distance - Ly / 4.0) / width));
}

void runSimulation() {
    std::cout << "rho1 = " << rho1 << std::endl;
    std::cout << "rho2 = " << rho2 << std::endl;
    std::cout << "mu1 = " << mu1 << std::endl;
    std::cout << "mu2 = " << mu2 << std::endl;
    std::cout << "amplitude = " << amplitude << std::endl;

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

    const double k = 2.0 * M_PI / Ly;
    solver.initialize([&](int, int j) {
        vb::NodeState state;
        state.c = volumeFraction(j);
        const double mu = state.c * mu1 + (1.0 - state.c) * mu2;
        state.ux = amplitude * mu2 / mu * sin(k * j);
        return state;
    });
    solver.set_body_force([&](int, int j, double* fx, double* fy) {
        *fx = amplitude * mu2 * k * k * sin(k * j);
        *fy = 0.0;
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
        std::cerr << "Unknown gradient stencil '" << argv[1] << "'; expected E4, E6 or E8." << std::endl;
        return 2;
    }
    double density_ratio = rho1 / rho2;
    if (argc > 2 && !(parseNumber(argv[2], &density_ratio) && density_ratio > 0.0)) {
        std::cerr << "Invalid density ratio '" << argv[2] << "'; expected a positive number." << std::endl;
        return 2;
    }
    if (argc > 3 && !parseNumber(argv[3], &amplitude)) {
        std::cerr << "Invalid amplitude '" << argv[3] << "'; expected a number." << std::endl;
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
