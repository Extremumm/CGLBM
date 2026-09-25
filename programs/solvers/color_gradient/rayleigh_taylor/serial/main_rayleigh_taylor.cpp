/// Rayleigh-Taylor instability, serial reference resolution.
///
/// The dense component sits on top of the light one and the interface is
/// perturbed by one cosine wavelength. With no surface tension to hold it, the
/// perturbation grows: the dense fluid falls in a spike, the light one rises in
/// a bubble, and the shear along the flanks rolls them up.
///
/// This is the reference resolution, run serially. `rayleigh_taylor_omp` is the
/// same case at production resolution across OpenMP threads.
///
/// The domain is periodic along x and closed by resting walls along y.

#include <cmath>
#include <iostream>

#include "lbm/case_config.h"
#include "lbm/solver.h"

namespace {

using cglbm::lbm::CaseConfig;

const double c_dx = 1.e-5;                        // m
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

CaseConfig rayleigh_taylor_case() {
    CaseConfig config;
    config.name = "rayleigh_taylor";
    config.nx = 128;
    config.ny = 1028;
    config.steps = 5000000;
    config.interval = 10000;

    config.physics.rho1 = 4.;
    config.physics.rho2 = 1.;
    config.physics.c1 = 347. / (c_dx / c_dt);
    config.physics.c2 = 347. / (c_dx / c_dt);
    config.physics.radius = 10.;
    // No surface tension: nothing opposes the instability.
    config.physics.sigma = 0. / (c_dx * c_dx * c_dx / c_dt / c_dt);
    config.physics.nu = 1.e-4 / (c_dx * c_dx / c_dt);
    config.physics.nu_b = 1.e-4 / (c_dx * c_dx / c_dt);
    config.physics.gravity = 9.81e2 / (c_dx / c_dt / c_dt);  // m s^-2, lattice units
    config.physics.ch_width_init = 1.1 * config.units.dx;
    config.physics.ch_width_ope = 1.6 * config.units.dx;
    config.physics.p2_inf = 0.;
    config.physics.p1_inf = cglbm::lbm::matched_p1_inf(config.physics);

    config.boundary = cglbm::lbm::Boundary::WallY;
    // Dense component on top: the unstable arrangement.
    config.initial_phase = cglbm::lbm::cosine_layer(0.2, /*inverted=*/false);
    config.stencil = cglbm::lbm::GradientStencil::E4;

    return config;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = rayleigh_taylor_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "rayleigh_taylor")) {
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
        std::cerr << "rayleigh_taylor: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
