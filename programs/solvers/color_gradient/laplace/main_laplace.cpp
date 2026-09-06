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

#include <cmath>
#include <iostream>

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
    // scheme stops: at 1000 it diverges after 8.7e4 steps. See
    // docs/numerics.md.
    config.surface_tension = cglbm::lbm::SurfaceTension::ContinuumSurfaceForce;
    config.warn_phase_out_of_range = true;

    return config;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = laplace_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "laplace")) {
    case cglbm::lbm::CommandLineResult::Finished:
        return 0;
    case cglbm::lbm::CommandLineResult::Error:
        return 2;
    case cglbm::lbm::CommandLineResult::Run:
        break;
    }

    std::cout << cglbm::lbm::describe(config) << std::endl;
    try {
        cglbm::lbm::Solver solver(config);
        solver.run();
    } catch (const std::exception& error) {
        std::cerr << "laplace: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
