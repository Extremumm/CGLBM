/// Rayleigh-Taylor instability in three dimensions, single mode.
///
/// The dense component sits on top of the light one, the interface carries one
/// cosine wavelength along x and one along z, and nothing holds it up: the
/// dense fluid falls in a spike and the light one rises in a bubble. It is the
/// three-dimensional counterpart of `rayleigh_taylor`, and what it adds to the
/// three-dimensional solver's coverage is the two ingredients `laplace_3d` and
/// `oscillation_3d` leave untouched: the walls along y and the body force.
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
/// score against. At this resolution it cannot be: the viscous correction
/// `-nu k^2` is 28 % of `n` at the mean of the two kinematic viscosities, and
/// the interface is 4.9 nodes wide between phi_N = +0.9 and -0.9, a sixth of
/// the wavelength it is being perturbed at. Lowering the viscosity enough to
/// make the comparison sharp puts tau below 0.52, and raising gravity enough to
/// do it makes the hydrostatic pressure a third of the bulk pressure. What the
/// tests beside this file assert is therefore what can be asserted: that the
/// walls hold every gram of both components, that the mode keeps the symmetry
/// it started with, and that the spike falls and the bubble rises.
/// `oscillation_3d` is where a number is checked against a formula.
///
/// Output is the z = nz/2 slice, which cuts through both. The displacement
/// along it is `-A nx cos(2 pi x / nx)`, so the slice holds a spike at x = 0
/// and a bubble at x = nx/2, and the tests read the two off it. The droplet
/// track of `oscillation_3d` would say nothing here, so it is off.
///
/// **What it does.** The interface starts 3.0 nodes peak to peak and reaches 55
/// by step 3000 -- the spike down to y = 32 and the bubble up to y = 86 in a
/// column of 128 -- growing with an e-folding time of 815 steps, fitted over
/// steps 500 to 2000. The inviscid rate would give 268; the difference is the
/// `-nu k^2` and the interface width the header above declines to score
/// against.
///
/// The spike front and the peak-to-peak separation both grow monotonically from
/// step 500 on. The bubble front alone does not, and the reason is measurement
/// rather than flow: the bubble's top is flat, so the highest column changes
/// from one output to the next and `nanmax` over the row wobbles by up to a
/// third of a node while the front climbs 21 of them. The tests beside this
/// file therefore score the separation, which is the quantity the instability
/// actually grows.
///
/// The domain is periodic along x and z, and closed by resting walls along y.

#include <iostream>

#include "lbm/case_config.h"
#include "lbm/two_population_solver_3d.h"

namespace {

using cglbm::lbm::CaseConfig;

/// Dynamic viscosity shared by the two fluids, in lattice units.
///
/// Matched, so tau is uniform across the interface; 0.04 puts it at 0.6.
///
/// As in `oscillation_3d`, this is a floor rather than a preference, and it is
/// set by the positivity bound of the enhanced equilibrium on D3Q19: an axial
/// population is `rho_k (phi_axial + u/2)` and goes negative once `|u|` passes
/// `(c_s^k)^2 / 3`, which is 4.4e-2 at this density ratio. Measured, mu = 0.02
/// diverges between steps 250 and 500, 0.025 by step 500 and 0.03 by step 750;
/// 0.04 runs the whole 3000 steps with the phase field inside [-1, 1]
/// throughout. docs/numerics.md, *The positivity bound on D3Q19*, has the rest.
///
/// It is not free. The viscous correction to the inviscid growth rate,
/// `-nu k^2` at the mean of the two kinematic viscosities, is 55 % of
/// `sqrt(A g k)` here against 28 % at the mu this case used to claim -- which
/// is why the header above is careful that this is a demonstration and not a
/// validation of the growth rate.
const double mu = 0.04;

/// Initial displacement of the interface, as a fraction of nx.
///
/// 1.5 nodes at nx = 32: a twentieth of the wavelength, so the mode starts
/// inside the linear regime it was chosen from, and a third of the width of the
/// diffuse interface, which is fine because the interface position is read from
/// where the profile crosses zero rather than from any single node. Measured,
/// the initial contour follows the prescribed cosine to 0.012 nodes.
const double amplitude = 0.047;

CaseConfig rayleigh_taylor_3d_case() {
    CaseConfig config;
    config.name = "rayleigh_taylor_3d";
    config.nx = 32;
    config.ny = 128;
    config.nz = 32;
    // The spike reaches the bottom wall at about 4000 steps and the two layers
    // then overturn, which is a real outcome and a different case; 3000 leaves
    // the spike 22 nodes clear of the wall with the mode still recognisable.
    config.steps = 3000;
    config.interval = 250;

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
