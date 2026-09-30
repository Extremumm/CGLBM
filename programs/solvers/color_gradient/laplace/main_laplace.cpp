/// Laplace's law across a static droplet.
///
/// A droplet of radius R held by surface tension sigma carries a pressure jump
/// dp = sigma / R across its interface. The case starts from a droplet whose
/// pressure field satisfies that exactly -- `matched_p1_inf` puts the -sigma/R
/// offset into the pressure at infinity of component 1 -- and lets the scheme
/// relax. What it relaxes to is the measurement.
///
/// The domain is periodic on both axes: the droplet floats in an unbounded
/// fluid, so nothing but the scheme sets the jump.
///
/// See programs/solvers/color_gradient/laplace/tests/ for what is asserted of
/// the result, and docs/numerics.md for the gap between the relaxed jump and
/// the analytic one.
///
/// `--translate=U`, which this program reads before the shared options, starts
/// the whole box moving at U along x instead: the droplet should then be
/// carried along unchanged, and how far it falls behind the flow is the
/// Galilean-invariance test of tests/test_laplace_translating.py.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "lbm/case_config.h"
#include "lbm/solver.h"

namespace {

using cglbm::lbm::CaseConfig;

// Conversion factors from lattice to physical units. Following CONTRIBUTING,
// this is the only place a physical unit appears: everything handed to the
// solver below is already in lattice units.
const double c_dx = 1.e-5;                        // m
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

CaseConfig laplace_case() {
    CaseConfig config;
    config.name = "laplace";
    config.nx = 128;
    config.ny = 128;
    config.steps = 30000;
    config.interval = 1000;

    config.physics.rho1 = 20.;                 // kg m^-3, the droplet
    config.physics.rho2 = 1.;                  // kg m^-3, the surrounding fluid
    config.physics.c1 = 347. / (c_dx / c_dt);  // speed of sound, lattice units
    config.physics.c2 = 347. / (c_dx / c_dt);
    config.physics.radius = 10.;
    config.physics.sigma = 1. / (c_dx * c_dx * c_dx / c_dt / c_dt);
    // Kinematic viscosities, p. 284 of Kruger et al.
    config.physics.nu = 1.e-2 / (c_dx * c_dx / c_dt);
    config.physics.nu_b = 1.e-2 / (c_dx * c_dx / c_dt);
    config.physics.ch_width_init = 1.1 * config.units.dx;
    config.physics.ch_width_ope = 1.6 * config.units.dx;
    config.physics.p2_inf = 0.;
    config.physics.p1_inf = cglbm::lbm::matched_p1_inf(config.physics);
    // Keep it matched if --rho1, --sigma or --radius move on the command line.
    config.matched_pressure_offset = true;

    config.boundary = cglbm::lbm::Boundary::PeriodicY;
    config.initial_phase = cglbm::lbm::droplet_interface();
    // No walls here, so the gradient stencil is free to reach further.
    // Leclaire, Reggio & Trepanier, Computers & Fluids 48, 98 (2011) show the
    // nearest-neighbour gradient is what limits the model at large density
    // contrast, and that a higher-order isotropic one cuts spurious currents.
    config.stencil = cglbm::lbm::GradientStencil::E8;

    // The interface is the surface where the two components occupy equal
    // volume, and that is the zero of the bulk-normalised phase field, not of
    // the colour field -- at a density ratio of 1000 the colour field's zero
    // sits at phi = 0.998, deep inside the light fluid. Ba et al. Eq. (21).
    // Both the initial profile and the gradient the tension rides on are
    // therefore taken of phi_N.
    config.interface_field = cglbm::lbm::InterfaceField::BulkNormalised;
    // With the interface in the right place, the capillary stress of
    // Omega^(2) has to be injected where tau is largest, and it is divided by
    // tau: the tension collapses. Ba et al.'s continuum-surface-force operator
    // reaches the momentum equation through Guo's forcing instead, whose
    // factor stays bounded, and it is what holds Laplace's law to a couple of
    // per cent up to a density ratio of about 500 -- which is where this
    // configuration stops: at 1000 it diverges after 8.7e4 steps. Past it,
    // give the two fluids the same dynamic viscosity (--nu, --nu2) and add
    // --viscosity-mixing=dynamic, which keeps tau uniform across the
    // interface, and --surface-tension=stress, which conserves momentum:
    // tests/test_laplace_high_density_ratio.py runs that at 1e4. See
    // docs/numerics.md.
    config.surface_tension = cglbm::lbm::SurfaceTension::ContinuumSurfaceForce;
    config.warn_phase_out_of_range = true;

    return config;
}

/// Remove `--translate=U` from the arguments and store U; false if the value
/// is not a number.
bool take_translation(int* argc, char** argv, double* speed) {
    const std::string prefix = "--translate=";
    int kept = 1;
    for (int n = 1; n < *argc; ++n) {
        const std::string argument = argv[n];
        if (argument.rfind(prefix, 0) != 0) {
            argv[kept++] = argv[n];
            continue;
        }
        const std::string value = argument.substr(prefix.size());
        char* end = nullptr;
        *speed = std::strtod(value.c_str(), &end);
        if (value.empty() || *end != '\0' || !std::isfinite(*speed)) {
            std::cerr << "laplace: --translate expects a number, got '" << value << "'."
                      << std::endl;
            return false;
        }
    }
    *argc = kept;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    double speed = 0.0;
    if (!take_translation(&argc, argv, &speed)) {
        return 2;
    }
    CaseConfig config = laplace_case();
    if (speed != 0.0) {
        config.initial_velocity = [speed](const CaseConfig&, int, int, double* u_x, double* u_y) {
            *u_x = speed;
            *u_y = 0.0;
        };
    }
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "laplace")) {
    case cglbm::lbm::CommandLineResult::Finished:
        std::cout << "  --translate=U       start the whole box moving at U along x (0)"
                  << std::endl;
        return 0;
    case cglbm::lbm::CommandLineResult::Error:
        return 2;
    case cglbm::lbm::CommandLineResult::Run:
        break;
    }

    std::cout << cglbm::lbm::describe(config) << "\ntranslate = " << speed << std::endl;
    try {
        cglbm::lbm::Solver solver(config);
        solver.run();
    } catch (const std::exception& error) {
        std::cerr << "laplace: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
