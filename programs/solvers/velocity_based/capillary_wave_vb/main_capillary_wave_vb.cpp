#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "lbm/isotropic_gradient.h"
#include "lbm/velocity_based_solver.h"

// A capillary wave on a heavy layer, with the velocity-based scheme of
// src/lbm/velocity_based.h: the same case as the colour-gradient
// `capillary_wave`, against the same exact viscous normal mode.
//
// Usage: capillary_wave_vb [E4|E6|E8] [density_ratio] [mu1] [steps]
//
//   density_ratio  rho1/rho2, default 1000.
//   mu1            dynamic viscosity of the heavy layer, lattice units,
//                  default 2. The light fluid's is 0.05.
//   steps          default 25000.
//
// A band of component 1 fills Ly/4 < y < 3Ly/4 of a doubly periodic box. The
// lower interface is displaced by `amplitude cos(2 pi x / Lx)` and released;
// `mode.csv` records every 50 steps the first Fourier coefficient of its
// height, read off the heavy fluid's volume in each column of the lower half.

namespace vb = cglbm::lbm::velocity_based;

namespace {

const int Lx = 64;
const int Ly = 128;
const double kPi = 3.14159265358979323846;

const double c_dx = 1.e-5;  // m : conversion factor from lattice units to physical units
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

const double rho2 = 1.;
const double mu2 = 0.05;
const double sigma = 1. / (c_dx * c_dx * c_dx / c_dt / c_dt);  // as in the laplace case
const double amplitude = 0.3;  // initial displacement of the lower interface, nodes
const int interval = 5000;     // grid output interval
const int track_interval = 50;

// Four ASCII grids per output, as the other programs write them.
void write_grids(const vb::Solver& solver, int timestep) {
    std::ofstream density("density_" + std::to_string(timestep) + ".csv");
    std::ofstream velocity("velocity_" + std::to_string(timestep) + ".csv");
    std::ofstream phase("phase_" + std::to_string(timestep) + ".csv");
    std::ofstream pressure("pressure_" + std::to_string(timestep) + ".csv");
    for (auto* file : {&density, &velocity, &phase, &pressure}) {
        file->precision(10);
    }
    for (int j = 0; j < Ly; ++j) {
        for (int i = 0; i < Lx; ++i) {
            const char* separator = i < Lx - 1 ? "," : "\n";
            density << solver.density(i, j) << separator;
            velocity << solver.velocity_x(i, j) << "," << solver.velocity_y(i, j) << separator;
            phase << solver.phase(i, j) << separator;
            pressure << solver.pressure(i, j) << separator;
        }
    }
}

// First cosine coefficient of the lower interface's height.
double mode_amplitude(const vb::Solver& solver) {
    double result = 0.0;
    for (int i = 0; i < Lx; ++i) {
        double volume = 0.0;
        for (int j = 0; j < Ly / 2; ++j) {
            volume += solver.volume_fraction(i, j);
        }
        result += 2.0 / Lx * (Ly / 2.0 - volume) * std::cos(2.0 * kPi * i / Lx);
    }
    return result;
}

// A number, or false.
bool parse_number(const char* text, double* value) {
    char* end = nullptr;
    const double parsed = std::strtod(text, &end);
    if (end == text || *end != '\0' || !std::isfinite(parsed)) {
        return false;
    }
    *value = parsed;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    cglbm::lbm::GradientStencil stencil = cglbm::lbm::GradientStencil::E8;
    double density_ratio = 1000.0;
    double mu1 = 2.0;
    double steps = 25000.0;
    if (argc > 1 && !cglbm::lbm::stencil_from_name(argv[1], &stencil)) {
        std::cerr << "Unknown gradient stencil '" << argv[1] << "'; expected E4, E6 or E8."
                  << std::endl;
        return 2;
    }
    if (argc > 2 && !(parse_number(argv[2], &density_ratio) && density_ratio > 1.0)) {
        std::cerr << "Invalid density ratio '" << argv[2] << "'; expected a number above 1."
                  << std::endl;
        return 2;
    }
    if (argc > 3 && !(parse_number(argv[3], &mu1) && mu1 > 0.0)) {
        std::cerr << "Invalid viscosity '" << argv[3] << "'; expected a positive number."
                  << std::endl;
        return 2;
    }
    if (argc > 4 && !(parse_number(argv[4], &steps) && steps >= 0.0)) {
        std::cerr << "Invalid step count '" << argv[4] << "'." << std::endl;
        return 2;
    }

    vb::SolverParameters parameters;
    parameters.nx = Lx;
    parameters.ny = Ly;
    parameters.rho1 = density_ratio * rho2;
    parameters.rho2 = rho2;
    parameters.mu1 = mu1;
    parameters.mu2 = mu2;
    parameters.surface_tension = sigma;
    parameters.stencil = stencil;

    std::cout << "gradient stencil = " << cglbm::lbm::stencil_name(stencil) << "\n"
              << "nx = " << Lx << "\n"
              << "ny = " << Ly << "\n"
              << "steps = " << static_cast<int>(steps) << "\n"
              << "rho1 = " << parameters.rho1 << "\n"
              << "rho2 = " << rho2 << "\n"
              << "mu1 = " << mu1 << "\n"
              << "mu2 = " << mu2 << "\n"
              << "sigma = " << sigma << "\n"
              << "width = " << parameters.width << "\n"
              << "amplitude = " << amplitude << std::endl;

    vb::Solver solver(parameters);
    solver.initialize([&](int i, int j) {
        const double lower = Ly / 4.0 + amplitude * std::cos(2.0 * kPi * i / Lx);
        const double upper = 3.0 * Ly / 4.0;
        const double phi = j < Ly / 2 ? std::tanh((j - lower) / parameters.width)
                                      : std::tanh((upper - j) / parameters.width);
        vb::NodeState state;
        state.c = 0.5 * (1.0 + phi);
        return state;
    });

    std::ofstream track("mode.csv");
    if (!track) {
        std::cerr << "capillary_wave_vb: cannot open mode.csv" << std::endl;
        return 1;
    }
    track.precision(12);
    track << "timestep,amplitude\n";
    write_grids(solver, 0);
    const int total = static_cast<int>(steps);
    for (int timestep = 0; timestep <= total; ++timestep) {
        if (timestep > 0) {
            solver.step();
        }
        if (timestep % track_interval == 0) {
            track << timestep << "," << mode_amplitude(solver) << "\n";
        }
        if (timestep > 0 && timestep % interval == 0) {
            std::cout << "Step " << timestep << std::endl;
            write_grids(solver, timestep);
        }
    }
    if (!track) {
        std::cerr << "capillary_wave_vb: cannot write mode.csv" << std::endl;
        return 1;
    }
    return 0;
}
