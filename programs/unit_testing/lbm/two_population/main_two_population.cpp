// Checks the two-population colour-gradient solver of src/lbm/two_population_solver.h.
//
// That model carries one distribution per fluid and gets its density ratio
// from the rest-particle weight of the equilibrium rather than from an equation
// of state, which is what lets it run where `Solver` cannot. The invariants it
// stands on are different, and they are what this reports, as `key = value`
// lines the pytest beside this file parses:
//
//  1. the equilibrium's moments -- mass, momentum and the momentum flux
//     rho_k (c_s^k)^2 + rho_k u u -- which is what makes it a fluid at all;
//  2. the rest weights against the density ratio they are supposed to encode;
//  3. each fluid's mass conserved separately, not just the total;
//  4. both distributions non-negative, which is the property the recolouring
//     had to be adapted for and the one that bounds the phase field;
//  5. phi_N inside [-1, 1];
//  6. a uniform fluid at rest staying at rest.
//
//   main_two_population [stencil]

#include <cmath>
#include <iostream>
#include <string>

#include "lbm/case_config.h"
#include "lbm/two_population_solver.h"

namespace {

using cglbm::lbm::Boundary;
using cglbm::lbm::CaseConfig;
using cglbm::lbm::Field;
using cglbm::lbm::GradientStencil;
using cglbm::lbm::TwoPopulationSolver;

/// A small static droplet at the given density ratio.
CaseConfig droplet_case(double ratio, GradientStencil stencil) {
    CaseConfig config;
    config.name = "two_population_unit";
    config.nx = 48;
    config.ny = 48;
    config.steps = 50;
    config.interval = 1000000;

    config.physics.rho1 = ratio;
    config.physics.rho2 = 1.;
    config.physics.radius = 12.;
    config.physics.sigma = 0.1;
    // Equal dynamic viscosity, so tau is uniform across the interface.
    const double mu = 0.1667;
    config.physics.nu = mu / config.physics.rho1;
    config.physics.nu_b = config.physics.nu;
    config.physics.nu2 = mu / config.physics.rho2;
    config.physics.nu_b2 = config.physics.nu2;
    config.physics.alpha2 = 0.2;
    config.physics.beta = 0.7;
    config.physics.ch_width_init = 1.1;

    config.boundary = Boundary::PeriodicY;
    config.initial_phase = cglbm::lbm::droplet_interface();
    config.stencil = stencil;
    return config;
}

double component_mass(const TwoPopulationSolver& solver, int fluid) {
    const Field& density = solver.component_density(fluid);
    double total = 0.0;
    for (int i = 0; i < density.nx(); ++i) {
        for (int j = 0; j < density.ny(); ++j) {
            total += density(i, j);
        }
    }
    return total;
}

double min_population(const TwoPopulationSolver& solver) {
    double worst = 0.0;
    for (int fluid = 0; fluid < 2; ++fluid) {
        const Field& f = solver.population(fluid);
        for (int i = 0; i < f.nx(); ++i) {
            for (int j = 0; j < f.ny(); ++j) {
                for (int k = 0; k < cglbm::lbm::kQ; ++k) {
                    worst = std::min(worst, f(i, j, k));
                }
            }
        }
    }
    return worst;
}

double max_abs_phase(const TwoPopulationSolver& solver) {
    const Field& phase = solver.phase();
    double worst = 0.0;
    for (int i = 0; i < phase.nx(); ++i) {
        for (int j = 0; j < phase.ny(); ++j) {
            worst = std::max(worst, std::fabs(phase(i, j)));
        }
    }
    return worst;
}

double max_speed(const TwoPopulationSolver& solver) {
    const Field& u = solver.velocity();
    double worst = 0.0;
    for (int i = 0; i < u.nx(); ++i) {
        for (int j = 0; j < u.ny(); ++j) {
            worst = std::max(worst, std::hypot(u(i, j, 0), u(i, j, 1)));
        }
    }
    return worst;
}

/// Run a droplet and report what has to hold at every step.
void report_droplet(const std::string& prefix, double ratio, GradientStencil stencil) {
    const CaseConfig config = droplet_case(ratio, stencil);
    TwoPopulationSolver solver(config);
    solver.initialize();

    const double mass1 = component_mass(solver, 0);
    const double mass2 = component_mass(solver, 1);
    double worst_phase = max_abs_phase(solver);
    double worst_population = min_population(solver);
    for (int step = 0; step < config.steps; ++step) {
        solver.step();
        solver.refresh();
        worst_phase = std::max(worst_phase, max_abs_phase(solver));
        worst_population = std::min(worst_population, min_population(solver));
    }

    std::cout.precision(17);
    std::cout << prefix << "_mass1_drift = " << std::fabs(component_mass(solver, 0) - mass1) / mass1
              << "\n"
              << prefix << "_mass2_drift = " << std::fabs(component_mass(solver, 1) - mass2) / mass2
              << "\n"
              << prefix << "_min_population = " << worst_population << "\n"
              << prefix << "_max_abs_phase = " << worst_phase << "\n"
              << prefix << "_max_speed = " << max_speed(solver) << std::endl;
}

/// A uniform fluid with no interface and no gravity must stay exactly at rest.
void report_rest_state() {
    CaseConfig config = droplet_case(20.0, GradientStencil::E8);
    config.initial_phase = [](const CaseConfig&, int, int) { return 1.0; };
    TwoPopulationSolver solver(config);
    solver.initialize();
    for (int step = 0; step < config.steps; ++step) {
        solver.step();
    }
    solver.refresh();
    std::cout.precision(17);
    std::cout << "rest_max_speed = " << max_speed(solver) << std::endl;
}

/// The equilibrium's zeroth, first and second moments, and the rest weights.
///
/// The second moment is the one that carries the model: it must come out as
/// `rho_k (c_s^k)^2 delta + rho_k u u`, with each fluid's own sound speed, and
/// that is where the density ratio enters the pressure.
void report_equilibrium() {
    std::cout.precision(17);
    double worst_mass = 0.0, worst_momentum = 0.0, worst_stress = 0.0, worst_ratio = 0.0;
    for (double ratio : {1.0, 20.0, 1000.0, 100000.0}) {
        CaseConfig config = droplet_case(ratio, GradientStencil::E8);
        TwoPopulationSolver solver(config);
        // (1 - alpha_2) / (1 - alpha_1) must be the density ratio, Ba Eq. (7).
        const double encoded = (1.0 - solver.alpha2()) / (1.0 - solver.alpha1());
        worst_ratio = std::max(worst_ratio, std::fabs(encoded - ratio) / ratio);

        // A uniform patch of one fluid at a known velocity: initialise it, then
        // read the moments back out of the distribution.
        for (int fluid = 0; fluid < 2; ++fluid) {
            const double cs_squared =
                0.6 * (1.0 - (fluid == 0 ? solver.alpha1() : solver.alpha2()));
            const double rho_k = fluid == 0 ? ratio : 1.0;
            const double u_x = 0.031, u_y = -0.017;
            double eq[cglbm::lbm::kQ];
            solver.equilibrium_for_test(fluid, rho_k, u_x, u_y, eq);

            double mass = 0.0, mx = 0.0, my = 0.0, pxx = 0.0, pxy = 0.0;
            for (int k = 0; k < cglbm::lbm::kQ; ++k) {
                mass += eq[k];
                mx += eq[k] * cglbm::lbm::kXi[k][0];
                my += eq[k] * cglbm::lbm::kXi[k][1];
                pxx += eq[k] * cglbm::lbm::kXi[k][0] * cglbm::lbm::kXi[k][0];
                pxy += eq[k] * cglbm::lbm::kXi[k][0] * cglbm::lbm::kXi[k][1];
            }
            worst_mass = std::max(worst_mass, std::fabs(mass - rho_k) / rho_k);
            worst_momentum = std::max(
                worst_momentum,
                std::max(std::fabs(mx - rho_k * u_x), std::fabs(my - rho_k * u_y)) / rho_k);
            const double pxx_exact = rho_k * (cs_squared + u_x * u_x);
            const double pxy_exact = rho_k * u_x * u_y;
            worst_stress =
                std::max(worst_stress,
                         std::max(std::fabs(pxx - pxx_exact), std::fabs(pxy - pxy_exact)) / rho_k);
        }
    }
    std::cout << "equilibrium_mass_error = " << worst_mass << "\n"
              << "equilibrium_momentum_error = " << worst_momentum << "\n"
              << "equilibrium_stress_error = " << worst_stress << "\n"
              << "alpha_density_ratio_error = " << worst_ratio << std::endl;
}

}  // namespace

int main(int argc, char** argv) {
    GradientStencil stencil = GradientStencil::E8;
    if (argc > 1 && !cglbm::lbm::stencil_from_name(argv[1], &stencil)) {
        std::cerr << "Unknown gradient stencil '" << argv[1] << "'; expected E4, E6 or E8."
                  << std::endl;
        return 2;
    }

    std::cout << "stencil = " << cglbm::lbm::stencil_name(stencil) << std::endl;
    try {
        report_equilibrium();
        report_droplet("r20", 20.0, stencil);
        report_droplet("r1000", 1000.0, stencil);
        report_droplet("r100000", 100000.0, stencil);
        report_rest_state();
    } catch (const std::exception& error) {
        std::cerr << "two_population: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
