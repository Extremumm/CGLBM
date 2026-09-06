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

    config.boundary = cglbm::lbm::Boundary::PeriodicY;
    config.initial_phase = cglbm::lbm::droplet_interface();
    // No walls here, so the gradient stencil is free to reach further.
    // Leclaire, Reggio & Trepanier, Computers & Fluids 48, 98 (2011) show the
    // nearest-neighbour gradient is what limits the model at large density
    // contrast, and that a higher-order isotropic one cuts spurious currents.
    config.stencil = cglbm::lbm::GradientStencil::E8;
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
