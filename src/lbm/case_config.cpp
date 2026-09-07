#include "lbm/case_config.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "omp/omp_environment.h"

namespace cglbm {
namespace lbm {
namespace {

/// pi, spelled out rather than taken from `M_PI`, which is not in the standard
/// headers under `-std=c++17`. Same value to the last bit.
constexpr double kPi = 3.14159265358979323846;

/// Read `--name=value` out of `argument`, returning false when it is some other
/// option.
bool option_value(const std::string& argument, const std::string& name, std::string* value) {
    const std::string prefix = "--" + name + "=";
    if (argument.compare(0, prefix.size(), prefix) != 0) {
        return false;
    }
    *value = argument.substr(prefix.size());
    return true;
}

/// Parse a positive integer, reporting `name` when it is not one.
bool positive_integer(const std::string& text, const char* name, int* out) {
    try {
        const int value = std::stoi(text);
        if (value <= 0) {
            std::cerr << "--" << name << " must be positive, got '" << text << "'." << std::endl;
            return false;
        }
        *out = value;
        return true;
    } catch (const std::exception&) {
        std::cerr << "--" << name << " expects an integer, got '" << text << "'." << std::endl;
        return false;
    }
}

void print_usage(const std::string& program_name) {
    std::cout << "usage: " << program_name << " [E4|E6|E8] [options]\n"
              << "\n"
              << "Run the colour-gradient case built into this program. Output goes to the\n"
              << "current working directory under fixed names, so give every run a directory\n"
              << "of its own.\n"
              << "\n"
              << "  --stencil=E4|E6|E8  isotropy order of the colour gradient\n"
              << "  --initial-state=equilibrium|eos|linear\n"
              << "                      how rho and p are laid down at t = 0; `equilibrium`\n"
              << "                      starts the interface in mechanical equilibrium, which\n"
              << "                      is what a density ratio past ~100 needs\n"
              << "  --initial-profile=colour|normalised\n"
              << "                      field the initial tanh profile is prescribed in;\n"
              << "                      `normalised` is what makes the droplet come out at\n"
              << "                      the radius the case asked for\n"
              << "  --interface-field=colour|normalised\n"
              << "                      field the colour gradient is taken of; `normalised`\n"
              << "                      divides each component by its bulk density, so its\n"
              << "                      zero contour is the interface at any density ratio\n"
              << "  --surface-tension=perturbation|csf\n"
              << "                      capillary stress in Omega^(2), or a body force from\n"
              << "                      an explicit curvature (Ba et al.); `csf` with\n"
              << "                      `--interface-field=normalised` is what reaches a\n"
              << "                      density ratio of a few hundred\n"
              << "  --recolouring=width|latva-kokko\n"
              << "                      segregation strength from the pressure and an\n"
              << "                      interface width, or from the density and beta\n"
              << "  --beta=X            segregation strength of `latva-kokko`, in (0, 1]\n"
              << "  --rho1=X, --rho2=X  bulk densities of the two components\n"
              << "  --sigma=X           surface tension\n"
              << "  --radius=X          prescribed interface radius\n"
              << "  --width=X           interface width the recolouring holds. Wider resolves\n"
              << "                      a steeper density profile, at the cost of a curvature\n"
              << "                      that varies more across it. 1.6 is what every shipped\n"
              << "                      case uses\n"
              << "  --width-init=X      interface width of the initial profile\n"
              << "  --nu=X, --nu-b=X    kinematic shear and bulk viscosity of component 1\n"
              << "  --nu2=X, --nu-b2=X  the same for component 2; unset means equal to\n"
              << "                      component 1's. Setting nu2 = nu rho1/rho2 matches\n"
              << "                      the dynamic viscosities, which makes tau uniform\n"
              << "                      across the interface\n"
              << "  --nx=N, --ny=N, --nz=N\n"
              << "                      lattice size, overriding the case default\n"
              << "  --steps=N           number of time steps\n"
              << "  --interval=N        write the CSV grids every N steps\n"
              << "  --precision=N       significant digits in the CSV output (default 6)\n"
              << "  --threads=N         OpenMP threads; implies parallel execution\n"
              << "  --help              show this message\n";
}

}  // namespace

double density_at_pressure(
    double phi, double target_pressure, double rho1, double rho2, const ComponentPair& components) {
    const double linear = rho1 * (0.5 + 0.5 * phi) + rho2 * (0.5 - 0.5 * phi);

    // The root lies between the two bulk densities when the target lies between
    // the two bulk pressures; the bracket is widened generously so that a case
    // with a large surface-tension offset still encloses it.
    double lo = 0.25 * std::min(rho1, rho2);
    double hi = 4.0 * std::max(rho1, rho2);
    double f_lo = pressure(lo, phi, components) - target_pressure;
    double f_hi = pressure(hi, phi, components) - target_pressure;
    if (f_lo * f_hi > 0.0) {
        return linear;  // no sign change: leave the profile alone
    }

    // Bisection. The equation of state is a square root of a quadratic, so it
    // is smooth and cheap; 100 halvings take the bracket below any tolerance
    // that matters here, and this runs once per node at t = 0 only.
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double mid = 0.5 * (lo + hi);
        const double f_mid = pressure(mid, phi, components) - target_pressure;
        if (f_mid == 0.0) {
            return mid;
        }
        if (f_lo * f_mid < 0.0) {
            hi = mid;
            f_hi = f_mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }
    (void) f_hi;
    return 0.5 * (lo + hi);
}

double normalised_phase(double phi, double rho1, double rho2) {
    // The recolouring step can overshoot |phi| = 1 by a rounding error, and the
    // denominator below is only guaranteed positive on [-1, 1].
    if (phi > 1.0) {
        phi = 1.0;
    } else if (phi < -1.0) {
        phi = -1.0;
    }
    // Each component's density divided by its own bulk value. The common
    // factor rho/2 cancels between numerator and denominator.
    const double heavy = rho2 * (1.0 + phi);
    const double light = rho1 * (1.0 - phi);
    return (heavy - light) / (heavy + light);
}

double phase_from_normalised(double phi_n, double rho1, double rho2) {
    if (phi_n > 1.0) {
        phi_n = 1.0;
    } else if (phi_n < -1.0) {
        phi_n = -1.0;
    }
    const double sum = rho1 + rho2;
    const double difference = rho2 - rho1;
    return (sum * phi_n - difference) / (sum - difference * phi_n);
}

double matched_p1_inf(const Physics& physics) {
    return physics.rho1 * physics.c1 * physics.c1 - physics.rho2 * physics.c2 * physics.c2 -
           physics.sigma / physics.radius;
}

PhaseFieldInit droplet_interface() {
    return [](const CaseConfig& config, int i, int j) {
        const int x0 = config.nx / 2;
        const int y0 = config.ny / 2;
        const double r = config.physics.radius * config.units.dx;
        const double distance = std::sqrt(static_cast<double>((i - x0) * (i - x0)) +
                                          static_cast<double>((j - y0) * (j - y0)));
        return -std::tanh((distance - r) / config.physics.ch_width_init);
    };
}

PhaseFieldInit3D droplet_interface_3d() {
    return [](const CaseConfig& config, int i, int j, int k) {
        const double x = i - config.nx / 2;
        const double y = j - config.ny / 2;
        const double z = k - config.nz / 2;
        const double distance = std::sqrt(x * x + y * y + z * z);
        const double radius = config.physics.radius * config.units.dx;
        return -std::tanh((distance - radius) / config.physics.ch_width_init);
    };
}

PhaseFieldInit cosine_layer(double amplitude, bool inverted) {
    return [amplitude, inverted](const CaseConfig& config, int i, int j) {
        // phi(x, y, 0) = tanh( (y - y0 - A L cos(-2 pi x / L)) / W )
        const int y0 = config.ny / 2;
        const double displacement = amplitude * config.nx * std::cos(-2. * kPi * i / config.nx);
        const double profile = std::tanh(((j - y0) - displacement) / config.physics.ch_width_ope);
        return inverted ? -profile : profile;
    };
}

namespace {

bool nonnegative_double(const std::string& value, const char* name, double* out);

/// Parse a strictly positive, finite number.
bool positive_double(const std::string& value, const char* name, double* out) {
    double parsed = 0.0;
    if (!nonnegative_double(value, name, &parsed)) {
        return false;
    }
    if (parsed <= 0.0) {
        std::cerr << name << " must be greater than zero; got '" << value << "'." << std::endl;
        return false;
    }
    *out = parsed;
    return true;
}

/// Parse a viscosity: a finite, non-negative number.
bool nonnegative_double(const std::string& value, const char* name, double* out) {
    try {
        size_t consumed = 0;
        const double parsed = std::stod(value, &consumed);
        if (consumed != value.size() || !std::isfinite(parsed) || parsed < 0.0) {
            std::cerr << name << " must be a non-negative number; got '" << value << "'."
                      << std::endl;
            return false;
        }
        *out = parsed;
        return true;
    } catch (const std::exception&) {
        std::cerr << name << " must be a non-negative number; got '" << value << "'." << std::endl;
        return false;
    }
}

/// Read "colour"/"color" or "normalised"/"normalized" into `field`.
bool interface_field_from_name(const std::string& value, InterfaceField* field) {
    if (value == "colour" || value == "color") {
        *field = InterfaceField::Colour;
        return true;
    }
    if (value == "normalised" || value == "normalized") {
        *field = InterfaceField::BulkNormalised;
        return true;
    }
    return false;
}

}  // namespace

CommandLineResult
parse_command_line(CaseConfig& config, int argc, char** argv, const std::string& program_name) {
    int threads = 0;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        std::string value;

        if (argument == "--help" || argument == "-h") {
            print_usage(program_name);
            return CommandLineResult::Finished;
        }
        if (option_value(argument, "stencil", &value)) {
            if (!stencil_from_name(value.c_str(), &config.stencil)) {
                std::cerr << "Unknown gradient stencil '" << value << "'; expected E4, E6 or E8."
                          << std::endl;
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "interface-field", &value)) {
            if (!interface_field_from_name(value, &config.interface_field)) {
                std::cerr << "Unknown interface field '" << value
                          << "'; expected colour or normalised." << std::endl;
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "initial-profile", &value)) {
            if (!interface_field_from_name(value, &config.initial_profile_field)) {
                std::cerr << "Unknown initial profile field '" << value
                          << "'; expected colour or normalised." << std::endl;
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "surface-tension", &value)) {
            if (value == "perturbation") {
                config.surface_tension = SurfaceTension::Perturbation;
            } else if (value == "csf") {
                config.surface_tension = SurfaceTension::ContinuumSurfaceForce;
            } else {
                std::cerr << "Unknown surface tension form '" << value
                          << "'; expected perturbation or csf." << std::endl;
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "recolouring", &value) ||
            option_value(argument, "recoloring", &value)) {
            if (value == "width") {
                config.recolouring = Recolouring::InterfaceWidth;
            } else if (value == "latva-kokko") {
                config.recolouring = Recolouring::LatvaKokko;
            } else {
                std::cerr << "Unknown recolouring '" << value << "'; expected width or latva-kokko."
                          << std::endl;
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "rho1", &value)) {
            if (!positive_double(value, "rho1", &config.physics.rho1)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "rho2", &value)) {
            if (!positive_double(value, "rho2", &config.physics.rho2)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "sigma", &value)) {
            if (!nonnegative_double(value, "sigma", &config.physics.sigma)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "radius", &value)) {
            if (!positive_double(value, "radius", &config.physics.radius)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "width", &value)) {
            if (!nonnegative_double(value, "width", &config.physics.ch_width_ope)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "width-init", &value)) {
            if (!nonnegative_double(value, "width-init", &config.physics.ch_width_init)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "nu", &value)) {
            if (!nonnegative_double(value, "nu", &config.physics.nu)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "nu2", &value)) {
            if (!nonnegative_double(value, "nu2", &config.physics.nu2)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "nu-b", &value)) {
            if (!nonnegative_double(value, "nu-b", &config.physics.nu_b)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "nu-b2", &value)) {
            if (!nonnegative_double(value, "nu-b2", &config.physics.nu_b2)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "beta", &value)) {
            const double beta = std::atof(value.c_str());
            if (!(beta > 0.0) || beta > 1.0) {
                std::cerr << "beta must lie in (0, 1]; got '" << value << "'." << std::endl;
                return CommandLineResult::Error;
            }
            config.physics.beta = beta;
            continue;
        }
        if (option_value(argument, "initial-state", &value)) {
            if (value == "equilibrium") {
                config.initial_state = InitialState::MechanicalEquilibrium;
            } else if (value == "eos") {
                config.initial_state = InitialState::EquationOfStateP;
            } else if (value == "linear") {
                config.initial_state = InitialState::LinearDensity;
            } else {
                std::cerr << "Unknown initial state '" << value
                          << "'; expected equilibrium, eos or linear." << std::endl;
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "nx", &value)) {
            if (!positive_integer(value, "nx", &config.nx)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "nz", &value)) {
            if (!positive_integer(value, "nz", &config.nz)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "ny", &value)) {
            if (!positive_integer(value, "ny", &config.ny)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "steps", &value)) {
            if (!positive_integer(value, "steps", &config.steps)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "interval", &value)) {
            if (!positive_integer(value, "interval", &config.interval)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "precision", &value)) {
            if (!positive_integer(value, "precision", &config.output_precision)) {
                return CommandLineResult::Error;
            }
            continue;
        }
        if (option_value(argument, "threads", &value)) {
            if (!positive_integer(value, "threads", &threads)) {
                return CommandLineResult::Error;
            }
            config.parallel = true;
            continue;
        }
        // The historical form: a bare stencil name as the first argument.
        if (!argument.empty() && argument[0] != '-') {
            if (stencil_from_name(argument.c_str(), &config.stencil)) {
                continue;
            }
        }
        std::cerr << "Unknown argument '" << argument << "'. Try --help." << std::endl;
        return CommandLineResult::Error;
    }

    if (threads > 0) {
        omp::set_thread_count(threads);
    }
    if (config.matched_pressure_offset) {
        // The case tied p1_inf to the densities, sigma and the radius, any of
        // which may just have been overridden.
        config.physics.p1_inf = matched_p1_inf(config.physics);
    }
    return CommandLineResult::Run;
}

std::string describe(const CaseConfig& config) {
    std::ostringstream out;
    out.precision(17);
    const Physics& physics = config.physics;

    out << "case = " << config.name << "\n"
        << "nx = " << config.nx << "\n"
        << "ny = " << config.ny << "\n"
        << "nz = " << config.nz << "\n"
        << "steps = " << config.steps << "\n"
        << "interval = " << config.interval << "\n"
        << "boundary = " << (config.boundary == Boundary::WallY ? "wall_y" : "periodic_y") << "\n"
        << "stencil = " << stencil_name(config.stencil) << "\n"
        << "initial_state = "
        << (config.initial_state == InitialState::MechanicalEquilibrium
                ? "equilibrium"
                : (config.initial_state == InitialState::EquationOfStateP ? "eos" : "linear"))
        << "\n"
        << "initial_profile = "
        << (config.initial_profile_field == InterfaceField::BulkNormalised ? "normalised"
                                                                           : "colour")
        << "\n"
        << "interface_field = "
        << (config.interface_field == InterfaceField::BulkNormalised ? "normalised" : "colour")
        << "\n"
        << "surface_tension = "
        << (config.surface_tension == SurfaceTension::ContinuumSurfaceForce ? "csf"
                                                                            : "perturbation")
        << "\n"
        << "recolouring = "
        << (config.recolouring == Recolouring::LatvaKokko ? "latva-kokko" : "width") << "\n"
        << "dx = " << config.units.dx << "\n"
        << "dt = " << config.units.dt << "\n"
        << "rho1 = " << physics.rho1 << "\n"
        << "rho2 = " << physics.rho2 << "\n"
        << "c1 = " << physics.c1 << "\n"
        << "c2 = " << physics.c2 << "\n"
        << "nu = " << physics.nu << "\n"
        << "nu_b = " << physics.nu_b << "\n"
        << "nu2 = " << (physics.nu2 < 0.0 ? physics.nu : physics.nu2) << "\n"
        << "nu_b2 = " << (physics.nu_b2 < 0.0 ? physics.nu_b : physics.nu_b2) << "\n"
        << "sigma = " << physics.sigma << "\n"
        << "radius = " << physics.radius << "\n"
        << "gravity = " << physics.gravity << "\n"
        << "ch_width_init = " << physics.ch_width_init << "\n"
        << "ch_width_ope = " << physics.ch_width_ope << "\n"
        << "beta = " << physics.beta << "\n"
        << "alpha2 = " << physics.alpha2 << "\n"
        << "p1_inf = " << physics.p1_inf << "\n"
        << "p2_inf = " << physics.p2_inf << "\n"
        << "output_precision = " << config.output_precision << "\n"
        << "parallel = " << (config.parallel ? "true" : "false");
    return out.str();
}

}  // namespace lbm
}  // namespace cglbm
