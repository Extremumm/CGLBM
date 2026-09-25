/// Laplace's law across a static droplet, at a density ratio of 1000.
///
/// The same benchmark as `laplace`, run with `cglbm::lbm::TwoPopulationSolver`
/// instead of `cglbm::lbm::Solver`. The two are different colour-gradient
/// models and the difference is which of them can be here at all: `laplace`
/// runs at a density ratio of 20 and its scheme diverges at 1000, while this
/// one holds Laplace's law to a fraction of a per cent there.
///
/// The parameters follow Ba et al. (2016) Table I so the result can be read
/// against theirs directly: R = 25 in 100^2, rho_2 = 1, sigma = 0.1,
/// alpha_2 = 0.2, beta = 0.7, periodic on both axes.
///
/// One choice is worth spelling out because it is not the paper's and it
/// matters more than anything else here. Ba et al. give both fluids the same
/// *kinematic* viscosity, 0.1667. Since tau = mu / (p dt) + 1/2 and the two
/// bulk pressures are equal in this model, that puts tau at 348 in the heavy
/// fluid against 0.85 in the light one, and a BGK collision at tau = 348 leaves
/// a viscous stress large enough to absorb more than half the capillary force:
/// measured, the pressure jump comes out at 43 % of the force actually applied.
/// Ba et al. carry that with an MRT collision, which relaxes the ghost moments
/// at their own rate; this program instead matches the *dynamic* viscosities,
/// which makes tau uniform at 0.85 and needs no MRT. It is a different physical
/// case -- a viscosity ratio of 1000 rather than 1 -- and it is stated here
/// rather than buried, because it is the reason the numbers below are good.
///
/// See docs/numerics.md for the measurements and programs/.../tests/.

#include <iostream>

#include "lbm/case_config.h"
#include "lbm/two_population_solver.h"

namespace {

using cglbm::lbm::CaseConfig;

/// Dynamic viscosity shared by the two fluids, in lattice units.
const double mu = 0.1667;

CaseConfig laplace_high_ratio_case() {
    CaseConfig config;
    config.name = "laplace_high_ratio";
    config.nx = 100;
    config.ny = 100;
    config.steps = 40000;
    config.interval = 5000;

    config.physics.rho1 = 1000.;  // the droplet
    config.physics.rho2 = 1.;     // the surrounding fluid
    config.physics.radius = 25.;
    config.physics.sigma = 0.1;

    // Equal dynamic viscosity, hence tau uniform across the interface. The
    // kinematic viscosities differ by the density ratio, which is the price.
    config.physics.nu = mu / config.physics.rho1;
    config.physics.nu_b = config.physics.nu;
    config.physics.nu2 = mu / config.physics.rho2;
    config.physics.nu_b2 = config.physics.nu2;

    // Rest-particle weight of the light fluid; the heavy one follows from the
    // density ratio. Ba et al. use 0.2.
    config.physics.alpha2 = 0.2;
    // Segregation strength of the Latva-Kokko recolouring, Ba et al.'s value.
    config.physics.beta = 0.7;
    config.physics.ch_width_init = 1.1;

    config.boundary = cglbm::lbm::Boundary::PeriodicY;
    config.initial_phase = cglbm::lbm::droplet_interface();
    // The colour gradient is differentiated twice here -- once for the normal,
    // once for the curvature -- and its isotropy is what Leclaire et al. (2011)
    // identify as the thing that carries this model past a ratio of 10^3.
    config.stencil = cglbm::lbm::GradientStencil::E8;

    return config;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = laplace_high_ratio_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "laplace_high_ratio")) {
    case cglbm::lbm::CommandLineResult::Finished:
        return 0;
    case cglbm::lbm::CommandLineResult::Error:
        return 2;
    case cglbm::lbm::CommandLineResult::Run:
        break;
    }

    std::cout << cglbm::lbm::describe(config) << std::endl;
    try {
        cglbm::lbm::TwoPopulationSolver solver(config);
        solver.run();
    } catch (const std::exception& error) {
        std::cerr << "laplace_high_ratio: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
