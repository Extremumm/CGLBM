#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "lbm/isotropic_gradient.h"
#include "lbm/velocity_based_output.h"
#include "lbm/velocity_based_solver.h"
#include "omp/omp_environment.h"

// A capillary wave on a heavy layer, with the velocity-based scheme of
// src/lbm/velocity_based.h: the same case as the colour-gradient
// `capillary_wave`, against the same exact viscous normal mode.
//
// Usage: capillary_wave_vb [E4|E6|E8|E10|E12] [density_ratio] [mu1] [steps] [amplitude]
//                          [--wavelength=N] [--fourth-order-phase]
//                          [--sixth-order-phase] [--viscosity=arithmetic|harmonic|laminate]
//                          [--threads=N]
//
//   density_ratio  rho1/rho2, default 1000.
//   mu1            dynamic viscosity of the heavy layer, lattice units,
//                  default 2. The light fluid's is 0.05.
//   steps          default 25000.
//   amplitude      initial displacement of the lower interface, nodes,
//                  default 0.3, where the wave is linear. At 1000 the
//                  colour-gradient solver's wave diverges from 4 upwards,
//                  where the interface moves at about 2e-3; this one does not
//                  (docs/numerics.md, "Oscillations against exact normal
//                  modes").
//   --wavelength=N the wavelength and width of the box, nodes, default 64; the
//                  box is twice as tall.
//   --fourth-order-phase
//                  builds the phase populations as
//                  SolverParameters::fourth_order_phase describes.
//   --sixth-order-phase
//                  and as SolverParameters::sixth_order_phase does.
//   --viscosity=M  how the two viscosities are mixed across the interface,
//                  SolverParameters::interface_viscosity; arithmetic by default.
//   --threads=N    runs the solver's loops on N OpenMP threads; the fields
//                  are the same to the last bit as on one.
//
// A band of component 1 fills Ly/4 < y < 3Ly/4 of a doubly periodic box. The
// lower interface is displaced by `amplitude cos(2 pi x / Lx)` and released;
// `mode.csv` records every 50 steps the first Fourier coefficient of its
// height, read off the heavy fluid's volume in each column of the lower half.

namespace vb = cglbm::lbm::velocity_based;

namespace {

// Set from --wavelength before anything is built.
int Lx = 64;
int Ly = 128;
const double kPi = 3.14159265358979323846;

const double c_dx = 1.e-5;  // m : conversion factor from lattice units to physical units
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

const double rho2 = 1.;
const double mu2 = 0.05;
const double sigma = 1. / (c_dx * c_dx * c_dx / c_dt / c_dt);  // as in the laplace case
const int interval = 5000;                                     // grid output interval
const int track_interval = 50;

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

// Takes the options out of argv, leaving the positional arguments.
bool take_options(
    int* argc, char** argv, int* wavelength, bool* fourth_order_phase, bool* sixth_order_phase) {
    const std::string wavelength_prefix = "--wavelength=";
    int kept = 1;
    for (int n = 1; n < *argc; ++n) {
        const std::string argument = argv[n];
        if (argument == "--fourth-order-phase") {
            *fourth_order_phase = true;
        } else if (argument == "--sixth-order-phase") {
            *sixth_order_phase = true;
        } else if (argument.rfind(wavelength_prefix, 0) == 0) {
            const std::string value = argument.substr(wavelength_prefix.size());
            double parsed = 0.0;
            if (!vb::parse_number(value.c_str(), &parsed) || parsed < 8.0 ||
                parsed != std::floor(parsed)) {
                std::cerr << "capillary_wave_vb: --wavelength expects a whole number of nodes, "
                             "at least 8, got '"
                          << value << "'." << std::endl;
                return false;
            }
            *wavelength = static_cast<int>(parsed);
        } else {
            argv[kept++] = argv[n];
        }
    }
    *argc = kept;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    cglbm::lbm::GradientStencil stencil = cglbm::lbm::GradientStencil::E8;
    double density_ratio = 1000.0;
    double mu1 = 2.0;
    double steps = 25000.0;
    double amplitude = 0.3;
    bool fourth_order_phase = false;
    bool sixth_order_phase = false;
    if (!take_options(&argc, argv, &Lx, &fourth_order_phase, &sixth_order_phase)) {
        return 2;
    }
    vb::InterfaceViscosity mixing = vb::InterfaceViscosity::Arithmetic;
    if (!vb::take_viscosity_option(&argc, argv, &mixing)) {
        std::cerr << "capillary_wave_vb: --viscosity expects arithmetic, harmonic or laminate."
                  << std::endl;
        return 2;
    }
    int threads = 0;
    if (!cglbm::omp::take_threads_option(&argc, argv, &threads)) {
        std::cerr << "capillary_wave_vb: --threads expects a whole number of at least 1."
                  << std::endl;
        return 2;
    }
    Ly = 2 * Lx;
    if (argc > 1 && !cglbm::lbm::stencil_from_name(argv[1], &stencil)) {
        std::cerr << "Unknown gradient stencil '" << argv[1]
                  << "'; expected E4, E6, E8, E10 or E12." << std::endl;
        return 2;
    }
    if (argc > 2 && !(vb::parse_number(argv[2], &density_ratio) && density_ratio > 1.0)) {
        std::cerr << "Invalid density ratio '" << argv[2] << "'; expected a number above 1."
                  << std::endl;
        return 2;
    }
    if (argc > 3 && !(vb::parse_number(argv[3], &mu1) && mu1 > 0.0)) {
        std::cerr << "Invalid viscosity '" << argv[3] << "'; expected a positive number."
                  << std::endl;
        return 2;
    }
    if (argc > 4 && !(vb::parse_number(argv[4], &steps) && steps >= 0.0)) {
        std::cerr << "Invalid step count '" << argv[4] << "'." << std::endl;
        return 2;
    }
    // at most a quarter of the layer, so the interface stays inside it
    if (argc > 5 && !(vb::parse_number(argv[5], &amplitude) && std::fabs(amplitude) < Ly / 8.0)) {
        std::cerr << "Invalid amplitude '" << argv[5] << "'; expected a number below " << Ly / 8
                  << " in magnitude." << std::endl;
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
    parameters.sixth_order_phase = sixth_order_phase;
    parameters.interface_viscosity = mixing;
    if (threads > 0) {
        cglbm::omp::set_thread_count(threads);
        parameters.parallel = true;
    }

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
              << "amplitude = " << amplitude << "\n"
              << "fourth_order_phase = " << fourth_order_phase << "\n"
              << "sixth_order_phase = " << sixth_order_phase << "\n"
              << "viscosity = " << vb::interface_viscosity_name(mixing) << std::endl;

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
    vb::write_fields(solver, 0);
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
            vb::write_fields(solver, timestep);
        }
    }
    if (!track) {
        std::cerr << "capillary_wave_vb: cannot write mode.csv" << std::endl;
        return 1;
    }
    return 0;
}
