/// Capillary oscillation of a perturbed interface.
///
/// A flat interface is displaced by one cosine wavelength and released. Surface
/// tension pulls it back, and the interface oscillates at a period set by the
/// two densities, the wavelength and sigma. `interface.csv` records the
/// crossing along the mid-plane at every step, which is what the period is
/// measured from.
///
/// The domain is periodic along x and closed by resting walls along y, applied
/// as half-way bounce-back.

#include <cmath>
#include <iostream>

#include "lbm/case_config.h"
#include "lbm/solver.h"

namespace {

using cglbm::lbm::CaseConfig;

// Conversion factors from lattice to physical units; everything below is in
// lattice units.
const double c_dx = 1.e-5;                        // m
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

CaseConfig capillary_case() {
    CaseConfig config;
    config.name = "capillary";
    config.nx = 128;
    config.ny = 128;
    config.steps = 10000;
    config.interval = 100;

    config.physics.rho1 = 4.;
    config.physics.rho2 = 1.;
    config.physics.c1 = 347. / (c_dx / c_dt);
    config.physics.c2 = 347. / (c_dx / c_dt);
    config.physics.radius = 10.;
    config.physics.sigma = 0.02 / (c_dx * c_dx * c_dx / c_dt / c_dt);
    config.physics.nu = 1.e-4 / (c_dx * c_dx / c_dt);
    config.physics.nu_b = 1.e-4 / (c_dx * c_dx / c_dt);
    config.physics.ch_width_init = 1.1 * config.units.dx;
    config.physics.ch_width_ope = 1.6 * config.units.dx;
    config.physics.p2_inf = 0.;
    config.physics.p1_inf = cglbm::lbm::matched_p1_inf(config.physics);

    config.boundary = cglbm::lbm::Boundary::WallY;
    config.initial_phase = cglbm::lbm::cosine_layer(0.2, /*inverted=*/true);
    // E4 reaches a single node, so it stays exact next to a wall; a wider
    // stencil becomes one-sided over two nodes there and has not been
    // validated against this case.
    config.stencil = cglbm::lbm::GradientStencil::E4;
    config.track_interface = true;

    return config;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = capillary_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "capillary")) {
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
        std::cerr << "capillary: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
