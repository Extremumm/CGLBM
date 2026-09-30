#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "lbm/isotropic_gradient.h"
#include "lbm/velocity_based_solver.h"

// Mode-2 oscillation of a two-dimensional droplet, with the velocity-based
// scheme of src/lbm/velocity_based.h: the same case as the colour-gradient
// `oscillation`, against the same exact viscous normal mode.
//
// Usage: oscillation_vb [E4|E6|E8] [density_ratio] [mu1] [steps]
//                       [--fourth-order-phase]
//
//   density_ratio  rho1/rho2, default 1000.
//   mu1            dynamic viscosity of the droplet, lattice units, default 2.
//                  The surrounding fluid's is 0.05.
//   steps          default 24000.
//   --fourth-order-phase
//                  builds the phase populations as
//                  SolverParameters::fourth_order_phase describes.
//
// The droplet is laid down as r = R' (1 + eps cos 2 theta), with R' shrunk so
// that it holds the area of a circle of radius 20, and released with the
// Laplace jump of that circle in place. `mode.csv` records every 50 steps the
// deformation
//
//     D = sum c (x^2 - y^2) / sum c (x^2 + y^2),
//
// the second moment of the volume fraction about the domain centre, which is
// 2 eps for a small deformation eps.

namespace vb = cglbm::lbm::velocity_based;

namespace {

const int Lx = 128;
const int Ly = 128;

const double c_dx = 1.e-5;  // m : conversion factor from lattice units to physical units
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

const double rho2 = 1.;
const double mu2 = 0.05;
const double radius = 20.;
const double sigma = 1. / (c_dx * c_dx * c_dx / c_dt / c_dt);  // as in the laplace case
const double deformation_0 = 0.03;  // initial mode-2 deformation, a fraction of the radius
const int interval = 6000;          // grid output interval
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

double deformation(const vb::Solver& solver) {
    double difference = 0.0;
    double sum = 0.0;
    for (int i = 0; i < Lx; ++i) {
        for (int j = 0; j < Ly; ++j) {
            const double x = i - Lx / 2;
            const double y = j - Ly / 2;
            const double c = solver.volume_fraction(i, j);
            difference += c * (x * x - y * y);
            sum += c * (x * x + y * y);
        }
    }
    return difference / sum;
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

// Takes --fourth-order-phase out of argv, leaving the positional arguments.
bool take_fourth_order_phase(int* argc, char** argv) {
    bool found = false;
    int kept = 1;
    for (int n = 1; n < *argc; ++n) {
        if (std::string(argv[n]) == "--fourth-order-phase") {
            found = true;
        } else {
            argv[kept++] = argv[n];
        }
    }
    *argc = kept;
    return found;
}

}  // namespace

int main(int argc, char** argv) {
    cglbm::lbm::GradientStencil stencil = cglbm::lbm::GradientStencil::E8;
    double density_ratio = 1000.0;
    double mu1 = 2.0;
    double steps = 24000.0;
    const bool fourth_order_phase = take_fourth_order_phase(&argc, argv);
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
    parameters.fourth_order_phase = fourth_order_phase;

    std::cout << "gradient stencil = " << cglbm::lbm::stencil_name(stencil) << "\n"
              << "nx = " << Lx << "\n"
              << "ny = " << Ly << "\n"
              << "steps = " << static_cast<int>(steps) << "\n"
              << "rho1 = " << parameters.rho1 << "\n"
              << "rho2 = " << rho2 << "\n"
              << "mu1 = " << mu1 << "\n"
              << "mu2 = " << mu2 << "\n"
              << "sigma = " << sigma << "\n"
              << "radius = " << radius << "\n"
              << "width = " << parameters.width << "\n"
              << "deformation = " << deformation_0 << "\n"
              << "fourth_order_phase = " << fourth_order_phase << std::endl;

    vb::Solver solver(parameters);
    const double scale = 1.0 / std::sqrt(1.0 + 0.5 * deformation_0 * deformation_0);
    solver.initialize([&](int i, int j) {
        const double x = i - Lx / 2;
        const double y = j - Ly / 2;
        const double distance = std::sqrt(x * x + y * y);
        const double cos_2theta = distance > 0.0 ? (x * x - y * y) / (distance * distance) : 0.0;
        const double r = radius * scale * (1.0 + deformation_0 * cos_2theta);
        vb::NodeState state;
        state.c = 0.5 * (1.0 - std::tanh((distance - r) / parameters.width));
        state.p = sigma / radius * state.c;
        return state;
    });

    std::ofstream track("mode.csv");
    if (!track) {
        std::cerr << "oscillation_vb: cannot open mode.csv" << std::endl;
        return 1;
    }
    track.precision(12);
    track << "timestep,deformation\n";
    write_grids(solver, 0);
    const int total = static_cast<int>(steps);
    for (int timestep = 0; timestep <= total; ++timestep) {
        if (timestep > 0) {
            solver.step();
        }
        if (timestep % track_interval == 0) {
            track << timestep << "," << deformation(solver) << "\n";
        }
        if (timestep > 0 && timestep % interval == 0) {
            std::cout << "Step " << timestep << std::endl;
            write_grids(solver, timestep);
        }
    }
    if (!track) {
        std::cerr << "oscillation_vb: cannot write mode.csv" << std::endl;
        return 1;
    }
    return 0;
}
