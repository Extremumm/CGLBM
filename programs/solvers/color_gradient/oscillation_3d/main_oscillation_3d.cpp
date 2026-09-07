/// Free oscillation of a deformed droplet in three dimensions.
///
/// `laplace_3d` measures what the three-dimensional scheme does when nothing
/// moves. This is the same scheme in motion, and it is checked against a result
/// that only exists in three dimensions.
///
/// A droplet is laid down as the spheroid `r(theta) = R [1 + eps P_2(cos
/// theta)]`, released from rest, and left to oscillate. Lamb gives the
/// frequency of that mode for two inviscid fluids:
///
///     omega_n^2 = n (n - 1) (n + 1) (n + 2) sigma
///                 / { R^3 [ (n + 1) rho_in + n rho_out ] }
///
/// which for the second harmonic, the lowest one a droplet has, is
///
///     omega_2^2 = 24 sigma / { R^3 (3 rho_in + 2 rho_out) }
///
/// The two-dimensional analogue of the same calculation gives
/// `omega^2 = 6 sigma / [R^3 (rho_in + rho_out)]`, so this is not a formula a
/// two-dimensional run could have been scored against: the mode shape, the
/// numerator and the way the two densities are weighted all change with the
/// dimension. It is the dynamic counterpart of `2 sigma / R` against
/// `sigma / R`.
///
/// The case is deliberately gentler than `laplace_3d`: a density ratio of 10
/// rather than 1000, because Lamb's result is inviscid and linear and a ratio
/// that large would need an interface far better resolved than R = 10 before
/// the comparison meant anything. Everything else is that case's -- sigma = 0.1,
/// R = 10, 48^3, matched dynamic viscosities -- and the surface tension it
/// produces has been measured there: the static jump at exactly these
/// parameters is 1.0235 times `2 sigma / R`, so what the frequency is scored
/// against is not in doubt.
///
/// The run appends the droplet's three semi-axes to `interface.csv` every step;
/// the mode-2 signal is `r_z - r_x`, a decaying sinusoid that
/// `pycglbm.oscillation` fits for a frequency and a decay rate.
///
/// **What it measures.** The frequency comes out below Lamb's, and by an amount
/// that falls as the droplet is better resolved:
///
///     R = 8  (48^3)   omega_0 / omega_Lamb = 0.834
///     R = 10 (48^3)                          0.872
///     R = 13 (64^3)                          0.894
///
/// which is `1 - 1.3/R` to within the scatter -- first order in the interface
/// width over the radius, the interface being about 1.6 nodes wide throughout.
/// The static tension is right to 2 % at all three, and the damping accounts
/// for under 1 % of it, so this is the diffuse interface and not the tension or
/// the viscosity. R = 10 is shipped because it is `laplace_3d`'s radius and
/// costs a quarter of what that case does; `--radius=13 --nx=64 --ny=64
/// --nz=64` is the better-resolved point above.
///
/// Reference
///  - H. Lamb, *Hydrodynamics*, 6th ed., Cambridge University Press (1932),
///    art. 275.

#include <iostream>

#include "lbm/case_config.h"
#include "lbm/two_population_solver_3d.h"

namespace {

using cglbm::lbm::CaseConfig;

/// Dynamic viscosity shared by the two fluids, in lattice units.
///
/// Matched, as in `laplace_3d`, so that tau is uniform across the interface --
/// at 0.55 here. It also sets how fast the oscillation dies: measured, the
/// amplitude falls to 0.46 of itself per period, which leaves four periods
/// worth fitting and a damping ratio of 0.12. Lower would be a cleaner
/// oscillation and a tau closer to the 0.5 the scheme is bounded by.
const double mu = 0.02;

/// Amplitude of the initial deformation, as a fraction of the radius.
///
/// Large enough that the semi-axes move over several nodes -- the interface is
/// about 1.6 nodes wide -- and small enough to stay in the linear regime Lamb's
/// frequency describes. The second-harmonic amplitude of a shape at eps = 0.1
/// is 5 % above the linear one.
const double deformation = 0.1;

CaseConfig oscillation_3d_case() {
    CaseConfig config;
    config.name = "oscillation_3d";
    config.nx = 48;
    config.ny = 48;
    config.nz = 48;
    config.steps = 4000;   // about 4.8 periods
    config.interval = 1000;

    config.physics.rho1 = 10.;  // the droplet
    config.physics.rho2 = 1.;   // the surrounding fluid
    config.physics.radius = 10.;
    config.physics.sigma = 0.1;

    config.physics.nu = mu / config.physics.rho1;
    config.physics.nu_b = config.physics.nu;
    config.physics.nu2 = mu / config.physics.rho2;
    config.physics.nu_b2 = config.physics.nu2;

    config.physics.alpha2 = 0.2;
    config.physics.beta = 0.7;
    config.physics.ch_width_init = 1.1;

    config.boundary = cglbm::lbm::Boundary::PeriodicY;
    config.initial_phase_3d = cglbm::lbm::oscillating_droplet_3d(deformation);
    config.stencil_3d = cglbm::lbm::GradientStencil3D::E6;
    config.track_interface = true;
    config.parallel = true;

    return config;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = oscillation_3d_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "oscillation_3d")) {
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
        std::cerr << "oscillation_3d: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
