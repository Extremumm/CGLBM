/// Laplace's law across a static droplet in three dimensions.
///
/// The three-dimensional counterpart of `laplace_high_ratio`, and the same
/// model and parameters: a sphere of the heavy fluid at a density ratio of
/// 1000, held by surface tension, with the dynamic viscosities matched so that
/// tau is uniform across the interface. What differs is the law being tested.
/// For a sphere the curvature is 2/R rather than 1/R, so
///
///     dp = 2 sigma / R
///
/// and a curvature operator ported from two dimensions without thought would
/// give half of it. `programs/unit_testing/lbm/two_population_3d` checks the
/// discrete curvature against -2/R directly, so that this case is not the first
/// place such a mistake would show.
///
/// The lattice is D3Q19 and the colour gradient uses the sixth-order isotropic
/// stencil: in two dimensions the isotropy of that gradient was worth an order
/// of magnitude on the spurious currents, and there is no reason to expect less
/// here.
///
/// This case runs its lattice loops across OpenMP threads by default, unlike
/// the two-dimensional ones. Three dimensions costs about a hundred times more
/// work per step and the loops write only their own node, so the result does
/// not depend on the thread count -- which the tests beside this file check.
///
/// Output is the z = nz/2 slice, in the same four CSV files a two-dimensional
/// run writes, so the existing post-processing reads it unchanged.
///
/// The run length is a compromise. The jump overshoots and then creeps down --
/// 1.0318 at 1.5e4 steps, 1.0298 at 2e4, 1.0262 at 3e4 -- so 1.5e4 is within
/// 0.6 % of converged and costs twelve minutes on a two dozen cores, where 3e4
/// costs half an hour. Raise `--steps` if that 0.6 % matters.
///
/// It used to be within 0.25 %. The enhanced-equilibrium fix took 417 times the
/// viscosity out of the heavy fluid, which is what had been holding the
/// approach so tight; the limit the jump converges to barely moved, from 1.0240
/// to 1.0262, because a static jump does not see the third moment the fix
/// repairs. What it does see is the spurious currents, now 7.2e-5 at 1.5e4
/// steps against 1.9e-5 before, settling to 3.2e-5 by 3e4.

#include <iostream>

#include "lbm/case_config.h"
#include "lbm/two_population_solver_3d.h"

namespace {

using cglbm::lbm::CaseConfig;

/// Dynamic viscosity shared by the two fluids, in lattice units.
const double mu = 0.1667;

CaseConfig laplace_3d_case() {
    CaseConfig config;
    config.name = "laplace_3d";
    config.nx = 48;
    config.ny = 48;
    config.nz = 48;
    config.steps = 15000;
    config.interval = 2500;

    config.physics.rho1 = 1000.;  // the droplet
    config.physics.rho2 = 1.;     // the surrounding fluid
    config.physics.radius = 10.;
    config.physics.sigma = 0.1;

    // Equal dynamic viscosity, hence tau uniform across the interface. See
    // docs/report for what happens when it is not.
    config.physics.nu = mu / config.physics.rho1;
    config.physics.nu_b = config.physics.nu;
    config.physics.nu2 = mu / config.physics.rho2;
    config.physics.nu_b2 = config.physics.nu2;

    config.physics.alpha2 = 0.2;
    config.physics.beta = 0.7;
    config.physics.ch_width_init = 1.1;

    config.boundary = cglbm::lbm::Boundary::PeriodicY;
    config.initial_phase_3d = cglbm::lbm::droplet_interface_3d();
    config.stencil_3d = cglbm::lbm::GradientStencil3D::E6;
    config.parallel = true;

    return config;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = laplace_3d_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "laplace_3d")) {
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
        std::cerr << "laplace_3d: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
