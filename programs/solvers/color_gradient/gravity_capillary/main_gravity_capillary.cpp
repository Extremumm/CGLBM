/// A perturbed interface under both gravity and surface tension.
///
/// Same geometry as the `capillary` case, with the dense component below and a
/// uniform acceleration along -y. Gravity and surface tension both pull the
/// interface back to flat, so it rings as a capillary-gravity wave, tracked in
/// `interface.csv`. At this lattice spacing gravity is the weaker of the two
/// by far: at the 1.28 mm wavelength the Bond number (rho1 - rho2) g /
/// (sigma k^2) is 6e-5, so the wave is nearly the `capillary` one. It used to be
/// laid down with the dense component on top, the Rayleigh-Taylor arrangement,
/// which the tension held down only because gravity is this weak.
///
/// The domain is periodic along x and closed by resting walls along y.

#include <cmath>
#include <iostream>
#include <limits>

#include "lbm/case_config.h"
#include "lbm/solver.h"

namespace {

using cglbm::lbm::CaseConfig;

const double c_dx = 1.e-5;                        // m
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

CaseConfig gravity_capillary_case() {
    CaseConfig config;
    config.name = "gravity_capillary";
    config.nx = 128;
    config.ny = 128;
    config.steps = 50000;
    config.interval = 1000;

    config.physics.rho1 = 4.;
    config.physics.rho2 = 1.;
    config.physics.c1 = 347. / (c_dx / c_dt);
    config.physics.c2 = 347. / (c_dx / c_dt);
    // A flat interface carries no Laplace jump: an infinite radius keeps
    // matched_p1_inf from giving the dense layer a droplet's.
    config.physics.radius = std::numeric_limits<double>::infinity();
    config.physics.sigma = 0.02 / (c_dx * c_dx * c_dx / c_dt / c_dt);
    config.physics.nu = 1.e-4 / (c_dx * c_dx / c_dt);
    config.physics.nu_b = 1.e-4 / (c_dx * c_dx / c_dt);
    config.physics.gravity = 9.81 / (c_dx / c_dt / c_dt);  // m s^-2, lattice units
    config.physics.ch_width_init = 1.1 * config.units.dx;
    config.physics.ch_width_ope = 1.6 * config.units.dx;
    config.physics.p2_inf = 0.;
    config.physics.p1_inf = cglbm::lbm::matched_p1_inf(config.physics);
    // Keep it matched if --rho1 or --sigma move on the command line.
    config.matched_pressure_offset = true;

    config.boundary = cglbm::lbm::Boundary::WallY;
    // Dense component below: the stable arrangement under gravity.
    config.initial_phase = cglbm::lbm::cosine_layer(0.2, /*inverted=*/true);
    config.stencil = cglbm::lbm::GradientStencil::E4;
    config.track_interface = true;

    return config;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = gravity_capillary_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "gravity_capillary")) {
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
        std::cerr << "gravity_capillary: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
