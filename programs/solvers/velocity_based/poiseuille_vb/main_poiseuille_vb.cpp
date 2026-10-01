#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "lbm/isotropic_gradient.h"
#include "lbm/velocity_based_output.h"
#include "lbm/velocity_based_solver.h"
#include "omp/omp_environment.h"

// Layered Poiseuille flow: two fluids side by side in a channel between
// resting walls, driven along it by a uniform pressure gradient, with the
// velocity-based scheme of src/lbm/velocity_based.h. The standard wall-bounded
// test of the phase-field models at large density ratios (Zu & He 2013,
// Fakhari et al. 2017, Subhedar 2022), against its exact profile.
//
// Usage: poiseuille_vb [E4|E6|E8|E10|E12] [density_ratio] [viscosity_ratio] [steps]
//                      [--ny=N] [--start=rest|exact] [--threads=N]
//
//   density_ratio    rho1/rho2, default 1000.
//   viscosity_ratio  mu1/mu2, default 1000: the same kinematic viscosity in
//                    both fluids, so that both reach the steady flow in the
//                    same few 10^4 steps.
//   steps            default 30000.
//   --ny=N           the channel's width in nodes, default 64; even.
//   --start=rest     starts both fluids at rest (the default);
//   --start=exact    starts them in the exact profile, for a heavy fluid
//                    whose viscous time ny^2 / nu1 no run can wait out.
//   --threads=N      runs the solver's loops on N OpenMP threads; the fields
//                    are the same to the last bit as on one.
//
// Component 1, the heavy one, fills the lower half of the channel and
// component 2 the upper; the walls sit halfway between the first and last
// rows and their ghosts, at y = -1/2 and y = ny - 1/2, and the interface on
// the centre line between them. A uniform force G per unit volume along x
// drives the flow, and with h = ny / 2 and y' = y - (ny - 1) / 2 its steady
// profile is
//
//     u(y') = G h^2 / (2 mu) [-(y'/h)^2 - (y'/h) (mu2 - mu1) / (mu1 + mu2)
//                             + 2 mu / (mu1 + mu2)],
//
// with mu = mu1 below the interface and mu2 above it: a parabola in each
// fluid, the two meeting with the same velocity and the same shear stress.
// G is set so that the fastest point moves at 0.01. The density ratio does
// not enter the steady flow, only how the scheme gets there and holds it.
//
// profile.csv holds, at the end, each row's velocity along x averaged over
// the channel's length, the exact one, and the volume fraction of component
// 1; error.csv, every 500 steps, the relative L2 distance between the two.
// The flow is uniform along x, which is 4 nodes long.

namespace vb = cglbm::lbm::velocity_based;

namespace {

const int Lx = 4;
int Ly = 64;

const double rho2 = 1.;
const double mu2 = rho2 / 6.;  // relaxation time 1 in component 2
const double peak = 0.01;      // fastest velocity of the exact profile
const int interval = 10000;    // grid output interval
const int error_interval = 500;

// The exact profile at height y' from the centre line, for a unit drive.
double exact_per_drive(double y, double h, double mu1) {
    const double mu = y < 0.0 ? mu1 : mu2;
    const double s = y / h;
    return h * h / (2.0 * mu) * (-s * s - s * (mu2 - mu1) / (mu1 + mu2) + 2.0 * mu / (mu1 + mu2));
}

// Takes the options out of argv, leaving the positional arguments.
bool take_options(int* argc, char** argv, int* ny, bool* start_exact, double* width) {
    const std::string ny_prefix = "--ny=";
    const std::string width_prefix = "--width=";
    int kept = 1;
    for (int n = 1; n < *argc; ++n) {
        const std::string argument = argv[n];
        if (argument == "--start=rest") {
            *start_exact = false;
        } else if (argument == "--start=exact") {
            *start_exact = true;
        } else if (argument.rfind(ny_prefix, 0) == 0) {
            const std::string value = argument.substr(ny_prefix.size());
            double parsed = 0.0;
            if (!vb::parse_number(value.c_str(), &parsed) || parsed < 16.0 ||
                parsed != std::floor(parsed) || static_cast<int>(parsed) % 2 != 0) {
                std::cerr << "poiseuille_vb: --ny expects an even whole number of nodes, "
                             "at least 16, got '"
                          << value << "'." << std::endl;
                return false;
            }
            *ny = static_cast<int>(parsed);
        } else if (argument.rfind(width_prefix, 0) == 0) {
            const std::string value = argument.substr(width_prefix.size());
            if (!vb::parse_number(value.c_str(), width) || !(*width > 0.0)) {
                std::cerr << "poiseuille_vb: --width expects a positive number of nodes, got '"
                          << value << "'." << std::endl;
                return false;
            }
        } else if (argument.rfind("--start=", 0) == 0) {
            std::cerr << "poiseuille_vb: --start expects rest or exact, got '" << argument << "'."
                      << std::endl;
            return false;
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
    double viscosity_ratio = 1000.0;
    double steps = 30000.0;
    bool start_exact = false;
    vb::InterfaceViscosity mixing = vb::InterfaceViscosity::Arithmetic;
    double width = vb::SolverParameters().width;
    if (!take_options(&argc, argv, &Ly, &start_exact, &width)) {
        return 2;
    }
    if (!vb::take_viscosity_option(&argc, argv, &mixing)) {
        std::cerr << "poiseuille_vb: --viscosity expects arithmetic, harmonic or laminate."
                  << std::endl;
        return 2;
    }
    int threads = 0;
    if (!cglbm::omp::take_threads_option(&argc, argv, &threads)) {
        std::cerr << "poiseuille_vb: --threads expects a whole number of at least 1." << std::endl;
        return 2;
    }
    if (argc > 1 && !cglbm::lbm::stencil_from_name(argv[1], &stencil)) {
        std::cerr << "Unknown gradient stencil '" << argv[1]
                  << "'; expected E4, E6, E8, E10 or E12." << std::endl;
        return 2;
    }
    if (argc > 2 && !(vb::parse_number(argv[2], &density_ratio) && density_ratio > 0.0)) {
        std::cerr << "Invalid density ratio '" << argv[2] << "'; expected a positive number."
                  << std::endl;
        return 2;
    }
    if (argc > 3 && !(vb::parse_number(argv[3], &viscosity_ratio) && viscosity_ratio > 0.0)) {
        std::cerr << "Invalid viscosity ratio '" << argv[3] << "'; expected a positive number."
                  << std::endl;
        return 2;
    }
    if (argc > 4 && !(vb::parse_number(argv[4], &steps) && steps >= 0.0)) {
        std::cerr << "Invalid step count '" << argv[4] << "'." << std::endl;
        return 2;
    }

    const double mu1 = viscosity_ratio * mu2;
    const double h = Ly / 2.0;
    const double centre = (Ly - 1) / 2.0;
    // the drive that puts the exact profile's fastest point at `peak`
    double fastest = 0.0;
    for (int n = 0; n <= 20000; ++n) {
        fastest = std::fmax(fastest, exact_per_drive(-h + n * (2.0 * h / 20000), h, mu1));
    }
    const double drive = peak / fastest;
    auto exact = [&](int j) { return drive * exact_per_drive(j - centre, h, mu1); };

    vb::SolverParameters parameters;
    parameters.nx = Lx;
    parameters.ny = Ly;
    parameters.rho1 = density_ratio * rho2;
    parameters.rho2 = rho2;
    parameters.mu1 = mu1;
    parameters.mu2 = mu2;
    parameters.stencil = stencil;
    parameters.boundary = cglbm::lbm::Boundary::WallY;
    parameters.interface_viscosity = mixing;
    parameters.width = width;
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
              << "width = " << parameters.width << "\n"
              << "viscosity = " << vb::interface_viscosity_name(mixing) << "\n"
              << "drive = " << drive << "\n"
              << "start = " << (start_exact ? "exact" : "rest") << std::endl;

    vb::Solver solver(parameters);
    solver.initialize([&](int, int j) {
        vb::NodeState state;
        state.c = 0.5 * (1.0 - std::tanh((j - centre) / parameters.width));
        state.ux = start_exact ? exact(j) : 0.0;
        return state;
    });
    solver.set_body_force([&](int, int, double* fx, double* fy) {
        *fx = drive;
        *fy = 0.0;
    });

    auto mean_velocity = [&](int j) {
        double sum = 0.0;
        for (int i = 0; i < Lx; ++i) {
            sum += solver.velocity_x(i, j);
        }
        return sum / Lx;
    };
    auto relative_error = [&]() {
        double difference = 0.0;
        double norm = 0.0;
        for (int j = 0; j < Ly; ++j) {
            const double e = exact(j);
            const double u = mean_velocity(j);
            difference += (u - e) * (u - e);
            norm += e * e;
        }
        return std::sqrt(difference / norm);
    };

    std::ofstream error("error.csv");
    if (!error) {
        std::cerr << "poiseuille_vb: cannot open error.csv" << std::endl;
        return 1;
    }
    error.precision(12);
    error << "timestep,relative_l2_error\n";
    vb::write_fields(solver, 0);
    const int total = static_cast<int>(steps);
    for (int timestep = 0; timestep <= total; ++timestep) {
        if (timestep > 0) {
            solver.step();
        }
        if (timestep % error_interval == 0) {
            error << timestep << "," << relative_error() << "\n";
        }
        if (timestep > 0 && timestep % interval == 0) {
            std::cout << "Step " << timestep << std::endl;
            vb::write_fields(solver, timestep);
        }
    }

    std::ofstream profile("profile.csv");
    if (!profile) {
        std::cerr << "poiseuille_vb: cannot open profile.csv" << std::endl;
        return 1;
    }
    profile.precision(12);
    profile << "j,y,ux,exact,c\n";
    for (int j = 0; j < Ly; ++j) {
        double c = 0.0;
        for (int i = 0; i < Lx; ++i) {
            c += solver.volume_fraction(i, j) / Lx;
        }
        profile << j << "," << j - centre << "," << mean_velocity(j) << "," << exact(j) << "," << c
                << "\n";
    }
    if (!error || !profile) {
        std::cerr << "poiseuille_vb: cannot write the output" << std::endl;
        return 1;
    }
    return 0;
}
