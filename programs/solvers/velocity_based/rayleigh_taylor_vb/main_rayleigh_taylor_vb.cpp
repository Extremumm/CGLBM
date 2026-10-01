#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "lbm/isotropic_gradient.h"
#include "lbm/velocity_based_output.h"
#include "lbm/velocity_based_solver.h"
#include "omp/omp_environment.h"

// The Rayleigh-Taylor instability at a large density ratio, between resting
// walls, with the velocity-based scheme of src/lbm/velocity_based.h: a heavy
// fluid above a light one, its interface displaced by a small cosine and
// released, against the exact growth rate of the viscous normal mode.
//
// Usage: rayleigh_taylor_vb [E4|E6|E8|E10|E12] [density_ratio] [mu1] [steps] [amplitude]
//                           [--wavelength=N] [--gravity=G] [--fourth-order-phase]
//                           [--viscosity=arithmetic|harmonic|laminate] [--threads=N]
//
//   density_ratio  rho1/rho2, default 1000: an Atwood number of 0.998.
//   mu1            dynamic viscosity of the heavy fluid, lattice units,
//                  default 2. The light fluid's is 0.05, as in
//                  capillary_wave_vb.
//   steps          default 12000.
//   amplitude      initial displacement of the interface, nodes, default
//                  0.05, so that it grows twenty-fold before k a reaches 0.1.
//   --wavelength=N the wavelength and width of the box, nodes, default 64; the
//                  box is four times as tall.
//   --gravity=G    gravitational acceleration, lattice units, default 1e-5.
//   --fourth-order-phase
//                  builds the phase populations as
//                  SolverParameters::fourth_order_phase describes.
//   --viscosity=M  how the two viscosities are mixed across the interface,
//                  SolverParameters::interface_viscosity; arithmetic by default.
//   --threads=N    runs the solver's loops on N OpenMP threads; the fields
//                  are the same to the last bit as on one.
//
// Component 1, the heavy one, fills the upper half of the box and component 2
// the lower; the walls sit at y = -1/2 and y = Ly - 1/2, and the interface at
// y = (Ly - 1)/2 + amplitude cos(2 pi x / Lx). The surface tension is the
// other cases', which holds every wavelength below 33 nodes at the default
// gravity. Gravity acts on rho - rho2: the light fluid is weightless and the
// heavy fluid carries the whole hydrostatic head, which the initial pressure
// balances column by column. The reference changes only the pressure, but
// in this scheme it has to be the light fluid: the populations carry
// p / (rho cs^2), which does not travel with the interface, so a node the
// heavy fluid sweeps over turns its pressure into rho1/rho2 times as much,
// and the light fluid's pressure has to stay small.
//
// mode.csv records every 50 steps the first Fourier coefficient of the
// interface's height, from the heavy fluid's volume in each column, and the
// heights of the bubble's and the spike's tips, where the column through
// each crosses c = 1/2, relative to the undisturbed interface.

namespace vb = cglbm::lbm::velocity_based;

namespace {

int Lx = 64;
int Ly = 256;
const double kPi = 3.14159265358979323846;

const double c_dx = 1.e-5;                        // m, as in the other cases
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

const double rho2 = 1.;
const double mu2 = 0.05;
const double sigma = 1. / (c_dx * c_dx * c_dx / c_dt / c_dt);  // as in the laplace case
const int interval = 2000;                                     // grid output interval
const int track_interval = 50;

double centre() {
    return (Ly - 1) / 2.0;
}

// Height of the interface in column i: the wall at Ly - 1/2 less the heavy
// fluid's volume above it.
double column_height(const vb::Solver& solver, int i) {
    double volume = 0.0;
    for (int j = 0; j < Ly; ++j) {
        volume += solver.volume_fraction(i, j);
    }
    return Ly - 0.5 - volume;
}

double mode_amplitude(const vb::Solver& solver) {
    double result = 0.0;
    for (int i = 0; i < Lx; ++i) {
        result += 2.0 / Lx * (column_height(solver, i) - centre()) * std::cos(2.0 * kPi * i / Lx);
    }
    return result;
}

// Where column i, read from the bottom, first reaches c = 1/2.
double crossing(const vb::Solver& solver, int i) {
    for (int j = 1; j < Ly; ++j) {
        const double below = solver.volume_fraction(i, j - 1);
        const double above = solver.volume_fraction(i, j);
        if (below < 0.5 && above >= 0.5) {
            return j - 1 + (0.5 - below) / (above - below);
        }
    }
    return Ly - 0.5;
}

// Takes the options out of argv, leaving the positional arguments.
bool take_options(
    int* argc, char** argv, int* wavelength, double* gravity, bool* fourth_order_phase) {
    const std::string wavelength_prefix = "--wavelength=";
    const std::string gravity_prefix = "--gravity=";
    int kept = 1;
    for (int n = 1; n < *argc; ++n) {
        const std::string argument = argv[n];
        if (argument == "--fourth-order-phase") {
            *fourth_order_phase = true;
        } else if (argument.rfind(wavelength_prefix, 0) == 0) {
            const std::string value = argument.substr(wavelength_prefix.size());
            double parsed = 0.0;
            if (!vb::parse_number(value.c_str(), &parsed) || parsed < 8.0 ||
                parsed != std::floor(parsed)) {
                std::cerr << "rayleigh_taylor_vb: --wavelength expects a whole number of "
                             "nodes, at least 8, got '"
                          << value << "'." << std::endl;
                return false;
            }
            *wavelength = static_cast<int>(parsed);
        } else if (argument.rfind(gravity_prefix, 0) == 0) {
            const std::string value = argument.substr(gravity_prefix.size());
            if (!vb::parse_number(value.c_str(), gravity) || !(*gravity > 0.0)) {
                std::cerr << "rayleigh_taylor_vb: --gravity expects a positive number, got '"
                          << value << "'." << std::endl;
                return false;
            }
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
    double steps = 12000.0;
    double amplitude = 0.05;
    double gravity = 1.0e-5;
    bool fourth_order_phase = false;
    if (!take_options(&argc, argv, &Lx, &gravity, &fourth_order_phase)) {
        return 2;
    }
    vb::InterfaceViscosity mixing = vb::InterfaceViscosity::Arithmetic;
    if (!vb::take_viscosity_option(&argc, argv, &mixing)) {
        std::cerr << "rayleigh_taylor_vb: --viscosity expects arithmetic, harmonic or laminate."
                  << std::endl;
        return 2;
    }
    int threads = 0;
    if (!cglbm::omp::take_threads_option(&argc, argv, &threads)) {
        std::cerr << "rayleigh_taylor_vb: --threads expects a whole number of at least 1."
                  << std::endl;
        return 2;
    }
    Ly = 4 * Lx;
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
    if (argc > 5 && !(vb::parse_number(argv[5], &amplitude) && std::fabs(amplitude) < Lx / 4.0)) {
        std::cerr << "Invalid amplitude '" << argv[5] << "'; expected a number below " << Lx / 4
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
    parameters.interface_viscosity = mixing;
    parameters.boundary = cglbm::lbm::Boundary::WallY;
    parameters.gravity_y = -gravity;
    parameters.gravity_reference_density = rho2;
    if (threads > 0) {
        cglbm::omp::set_thread_count(threads);
        parameters.parallel = true;
    }

    // full precision, so that the exact solution can be rebuilt from the log
    std::cout.precision(17);
    std::cout << "gradient stencil = " << cglbm::lbm::stencil_name(stencil) << "\n"
              << "nx = " << Lx << "\n"
              << "ny = " << Ly << "\n"
              << "steps = " << static_cast<int>(steps) << "\n"
              << "rho1 = " << parameters.rho1 << "\n"
              << "rho2 = " << rho2 << "\n"
              << "mu1 = " << mu1 << "\n"
              << "mu2 = " << mu2 << "\n"
              << "sigma = " << sigma << "\n"
              << "gravity = " << gravity << "\n"
              << "width = " << parameters.width << "\n"
              << "amplitude = " << amplitude << "\n"
              << "fourth_order_phase = " << fourth_order_phase << "\n"
              << "viscosity = " << vb::interface_viscosity_name(mixing) << std::endl;

    auto fraction = [&](int i, int j) {
        const double interface = centre() + amplitude * std::cos(2.0 * kPi * i / Lx);
        return 0.5 * (1.0 + std::tanh((j - interface) / parameters.width));
    };
    // The hydrostatic pressure of each column: zero in the light fluid, and
    // falling through the heavy one by its excess weight, integrated by the
    // trapezoidal rule so that the lattice gradient balances it.
    std::vector<double> pressure(static_cast<std::size_t>(Lx) * Ly, 0.0);
    for (int i = 0; i < Lx; ++i) {
        double p = 0.0;
        double previous = 0.0;
        for (int j = 0; j < Ly; ++j) {
            const double excess = fraction(i, j) * (parameters.rho1 - rho2);
            if (j > 0) {
                p -= 0.5 * gravity * (previous + excess);
            }
            pressure[static_cast<std::size_t>(i) * Ly + j] = p;
            previous = excess;
        }
    }

    vb::Solver solver(parameters);
    solver.initialize([&](int i, int j) {
        vb::NodeState state;
        state.c = fraction(i, j);
        state.p = pressure[static_cast<std::size_t>(i) * Ly + j];
        return state;
    });

    std::ofstream track("mode.csv");
    if (!track) {
        std::cerr << "rayleigh_taylor_vb: cannot open mode.csv" << std::endl;
        return 1;
    }
    track.precision(12);
    track << "timestep,amplitude,bubble,spike\n";
    vb::write_fields(solver, 0);
    const int total = static_cast<int>(steps);
    for (int timestep = 0; timestep <= total; ++timestep) {
        if (timestep > 0) {
            solver.step();
        }
        if (timestep % track_interval == 0) {
            // a positive amplitude raises the light fluid at x = 0, the bubble,
            // and lowers the heavy fluid at Lx / 2, the spike
            track << timestep << "," << mode_amplitude(solver) << ","
                  << crossing(solver, 0) - centre() << "," << crossing(solver, Lx / 2) - centre()
                  << "\n";
        }
        if (timestep > 0 && timestep % interval == 0) {
            std::cout << "Step " << timestep << std::endl;
            vb::write_fields(solver, timestep);
        }
    }
    if (!track) {
        std::cerr << "rayleigh_taylor_vb: cannot write mode.csv" << std::endl;
        return 1;
    }
    return 0;
}
