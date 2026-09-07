/// Rayleigh-Taylor instability in three dimensions, single mode.
///
/// The dense component sits on top of the light one, the interface carries one
/// cosine wavelength along x and one along z, and nothing holds it up: the
/// dense fluid falls in a spike, the light one rises in a bubble, and the shear
/// along the flanks rolls them up. It is the three-dimensional counterpart of
/// `rayleigh_taylor`, and what it adds to the three-dimensional solver's
/// coverage is the two ingredients `laplace_3d` and `oscillation_3d` leave
/// untouched: the walls along y and the body force.
///
/// The interface starts at `y0 + A nx cos(2 pi x / nx) cos(2 pi z / nz)`, the
/// single mode of a square cell. Where that product is +1 the interface is
/// high, so the light fluid has risen: a bubble, at the cell's four corners and
/// again at its centre. Where it is -1 the dense fluid has fallen: a spike, at
/// the midpoint of each of the cell's four edges. The mode's wavenumber is
/// `k = sqrt(kx^2 + kz^2) = 2 pi sqrt(2) / nx`, larger by sqrt(2) than the
/// two-dimensional case of the same width, which is where the third dimension
/// shows up.
///
/// **This case is a demonstration, not a validation, and the distinction is
/// worth stating.** The inviscid single-mode growth rate is `n = sqrt(A g k)`
/// with `A = (rho1 - rho2)/(rho1 + rho2)`, and it would be the obvious thing to
/// score against. At this resolution it cannot be: the viscous correction is
/// `-nu k^2`, which is 25 % of `n` here, the interface is 1.6 nodes wide
/// against a wavelength of 32, and the perturbation must start at an amplitude
/// comparable to that width to be visible at all. Lowering the viscosity enough
/// to make the comparison sharp puts tau below 0.52, and raising gravity enough
/// to do it makes the hydrostatic pressure a third of the bulk pressure. What
/// the tests beside this file assert is therefore what can be asserted: that
/// the walls hold every gram of both components, that the mode keeps the
/// symmetry it started with, and that the spike falls and the bubble rises.
/// `oscillation_3d` is where a number is checked against a formula.
///
/// Output is the z = nz/2 slice, which cuts through both. The displacement
/// along it is `-A nx cos(2 pi x / nx)`, so the slice holds a spike at x = 0 and
/// a bubble at x = nx/2, and the tests read the two off it. The droplet track of
/// `oscillation_3d` would say nothing here, so it is off.
///
/// The domain is periodic along x and z, and closed by resting walls along y.

#include <iostream>

#include "lbm/case_config.h"
#include "lbm/two_population_solver_3d.h"

namespace {

using cglbm::lbm::CaseConfig;

/// Dynamic viscosity shared by the two fluids, in lattice units.
///
/// Matched, so tau is uniform across the interface; 0.02 puts it at 0.55, low
/// enough to let the instability grow and high enough to stay well away from
/// the 0.5 the scheme is bounded by.
const double mu = 0.02;

/// Initial displacement of the interface, as a fraction of nx.
///
/// 1.5 nodes at nx = 32, which is about the width of the interface itself: a
/// smaller perturbation is not resolved, and a larger one starts outside the
/// linear regime the mode is chosen from.
const double amplitude = 0.047;

CaseConfig rayleigh_taylor_3d_case() {
    CaseConfig config;
    config.name = "rayleigh_taylor_3d";
    config.nx = 32;
    config.ny = 128;
    config.nz = 32;
    config.steps = 6000;
    config.interval = 500;

    config.physics.rho1 = 3.;  // the dense component, on top
    config.physics.rho2 = 1.;
    // Nothing opposes the instability.
    config.physics.sigma = 0.;
    config.physics.gravity = 1.e-4;

    config.physics.nu = mu / config.physics.rho1;
    config.physics.nu_b = config.physics.nu;
    config.physics.nu2 = mu / config.physics.rho2;
    config.physics.nu_b2 = config.physics.nu2;

    config.physics.alpha2 = 0.2;
    config.physics.beta = 0.7;
    config.physics.ch_width_ope = 1.6;

    config.boundary = cglbm::lbm::Boundary::WallY;
    // Dense component on top: the unstable arrangement, as in the
    // two-dimensional case.
    config.initial_phase_3d = cglbm::lbm::cosine_layer_3d(amplitude, /*inverted=*/false);
    // E4 reaches a single node, so it stays exact next to a wall; E6 reaches
    // two and becomes one-sided there, which is the same reason the
    // two-dimensional wall-bounded cases use the narrow stencil.
    config.stencil_3d = cglbm::lbm::GradientStencil3D::E4;
    config.parallel = true;

    return config;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = rayleigh_taylor_3d_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "rayleigh_taylor_3d")) {
    case cglbm::lbm::CommandLineResult::Finished:
        return 0;
    case cglbm::lbm::CommandLineResult::Error:
        return 2;
    case cglbm::lbm::CommandLineResult::Run:
        break;
    }

    std::cout << cglbm::lbm::describe(config) << std::endl;
    try {
        cglbm::lbm::TwoPopulationSolver3D solver(config);
        solver.run();
    } catch (const std::exception& error) {
        std::cerr << "rayleigh_taylor_3d: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
