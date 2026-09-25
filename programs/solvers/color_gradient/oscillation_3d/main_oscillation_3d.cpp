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
/// R = 10, 48^3, matched dynamic viscosities.
///
/// The run appends the droplet's three semi-axes to `interface.csv` every step;
/// the mode-2 signal is `r_z - r_x`, a decaying sinusoid that
/// `pycglbm.oscillation` fits for a frequency and a decay rate.
///
/// **What it measures.** The frequency comes out below Lamb's, and by an amount
/// that falls as the droplet is better resolved:
///
///     R = 8  (48^3)   omega_0 / omega_Lamb = 0.776
///     R = 10 (48^3)                          0.822
///     R = 13 (64^3)                          0.855
///
/// The gap closes as 1/R and closely: `gap * R` is 1.80, 1.78 and 1.88 at the
/// three radii. What they share is the interface, which the segregation
/// operator holds at the same 5.2 nodes between phi_N = +0.9 and -0.9 whatever
/// the radius is: 0.65 of R at R = 8, 0.52 at R = 10, 0.40 at R = 13. A deficit
/// first order in the interface width over the radius is what a diffuse
/// interface should give.
///
/// It is not the tension: the static jump at these parameters is 1.0235 times
/// `2 sigma / R` -- the same to four digits as before the equilibrium fix,
/// since a static jump does not see the third moment -- and quadrupling sigma
/// moves the ratio only to 0.863 while doubling the frequency. It is not the
/// density ratio: removing it makes the agreement worse, 0.798, at half the
/// damping.
///
/// It is *partly* the viscosity, and that is a change from what this case used
/// to claim. Doubling mu costs 10.3 points, not the 1.6 that was reported while
/// the scheme was running at 4.7 times the viscosity it had been asked for.
/// Since the damping ratio also falls with R -- 0.375, 0.346, 0.318 -- roughly
/// a fifth of the trend across the three radii is viscous rather than
/// geometric, and four fifths is resolution. Lowering mu to separate them is
/// not available; see the note on `mu` below. docs/numerics.md has the tables.
///
/// R = 10 is shipped because it is `laplace_3d`'s radius and costs a third of
/// what that case does; `--radius=13 --nx=64 --ny=64 --nz=64` is the
/// better-resolved point above, at about three times the cost.
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
/// at 0.65 here. It also sets how fast the oscillation dies: measured, the
/// amplitude falls to 0.10 of itself per period, for a damping ratio of 0.35.
///
/// **This is not a free choice, and 0.02 -- the obvious one -- does not run.**
/// The floor under it is the positivity bound of the enhanced equilibrium on
/// D3Q19, which is what makes this case harder than `laplace_3d` rather than
/// gentler. As `(c_s^k)^2` shrinks the equilibrium puts essentially all the
/// momentum on the six axial directions, leaving `rho_k (phi_axial + u/2)`
/// there, so an axial population goes negative once `|u|` passes
/// `(c_s^k)^2 / 3` -- 1.3e-2 at this density ratio, and half the room the same
/// construction leaves on D2Q9. `laplace_3d` never approaches its own bound
/// because its droplet does not move; here the droplet is the thing moving,
/// and the release transient is the worst moment of the run. At mu = 0.02 that
/// transient reaches |u| = 0.028, the populations go negative, and the run
/// diverges at step 86. At 0.045 it survives the release and runs away at step
/// 350. 0.05 holds but leaves a trace of the excursion; 0.06 holds it at
/// |u| = 0.007 and the populations at exactly zero from step 100 on.
///
/// The cost is that this is a heavily damped oscillation: the correction from
/// the observed frequency to the undamped one is 6.6 %, against 0.8 % at the
/// viscosity the case used to claim. Buying that back means getting under the
/// bound some other way than lowering mu: a smaller `deformation`, which
/// scales the release velocities directly, or a larger `(c_s^k)^2` from a lower
/// density ratio or an `alpha2` nearer zero. Both change what the case is, so
/// neither is done here.
///
/// docs/numerics.md, *The positivity bound on D3Q19*, has the measurements.
const double mu = 0.06;

/// Amplitude of the initial deformation, as a fraction of the radius.
///
/// Large enough that the semi-axes move over several nodes -- 1.5 of them
/// between the poles and the equator at t = 0, against an interface 5.2 nodes
/// wide -- and small enough to stay in the linear regime Lamb's frequency
/// describes. The second-harmonic amplitude of a shape at eps = 0.1 is 5 %
/// above the linear one.
const double deformation = 0.1;

CaseConfig oscillation_3d_case() {
    CaseConfig config;
    config.name = "oscillation_3d";
    config.nx = 48;
    config.ny = 48;
    config.nz = 48;
    config.steps = 4000;  // about 4.8 periods
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
