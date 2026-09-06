// Checks the shared colour-gradient kernel of src/lbm/solver.h.
//
// The solvers under programs/solvers each run one physical case for tens of
// thousands of steps, which makes them validation tests and slow ones. This
// program instead drives `cglbm::lbm::Solver` directly on a small lattice for a
// handful of steps and reports the invariants the scheme must hold at every
// step, as `key = value` lines the pytest beside this file parses:
//
//  1. mass conservation -- streaming moves populations, it does not create
//     them, and none of the three collision operators has a monopole;
//  2. the phase field stays in [-1, 1], the range the equation of state and
//     the recolouring operator are defined on;
//  3. a uniform state with no interface and no gravity stays exactly at rest,
//     which catches a spurious force in the source term;
//  4. the same, on the wall-bounded boundary, where bounce-back must not leak
//     mass through the walls.
//
//   main_lbm_solver [stencil]

#include <cmath>
#include <iostream>
#include <string>

#include "lbm/case_config.h"
#include "lbm/solver.h"

namespace {

using cglbm::lbm::Boundary;
using cglbm::lbm::CaseConfig;
using cglbm::lbm::GradientStencil;
using cglbm::lbm::Solver;

/// A small case with the shape of the shipped droplet runs, at a density ratio
/// of 20 so the equation of state is genuinely exercised.
CaseConfig small_case(Boundary boundary, GradientStencil stencil) {
    const double c_dx = 1.e-5;
    const double c_dt = c_dx / 347. / std::sqrt(3.);

    CaseConfig config;
    config.name = "unit";
    config.nx = 32;
    config.ny = 32;
    config.steps = 20;
    config.interval = 1000000;  // never reached: this program writes no files

    config.physics.rho1 = 20.;
    config.physics.rho2 = 1.;
    config.physics.c1 = 347. / (c_dx / c_dt);
    config.physics.c2 = 347. / (c_dx / c_dt);
    config.physics.radius = 6.;
    config.physics.sigma = 1. / (c_dx * c_dx * c_dx / c_dt / c_dt);
    config.physics.nu = 1.e-2 / (c_dx * c_dx / c_dt);
    config.physics.nu_b = 1.e-2 / (c_dx * c_dx / c_dt);
    config.physics.ch_width_init = 1.1;
    config.physics.ch_width_ope = 1.6;
    config.physics.p2_inf = 0.;
    config.physics.p1_inf = cglbm::lbm::matched_p1_inf(config.physics);

    config.boundary = boundary;
    config.initial_phase = cglbm::lbm::droplet_interface();
    config.stencil = stencil;
    return config;
}

/// Total density over the lattice.
double total_mass(const Solver& solver) {
    const cglbm::lbm::Field& rho = solver.density();
    double total = 0.0;
    for (int i = 0; i < rho.nx(); ++i) {
        for (int j = 0; j < rho.ny(); ++j) {
            total += rho(i, j);
        }
    }
    return total;
}

double max_abs_phase(const Solver& solver) {
    const cglbm::lbm::Field& phi = solver.phase();
    double worst = 0.0;
    for (int i = 0; i < phi.nx(); ++i) {
        for (int j = 0; j < phi.ny(); ++j) {
            worst = std::max(worst, std::fabs(phi(i, j)));
        }
    }
    return worst;
}

double max_speed(const Solver& solver) {
    const cglbm::lbm::Field& u = solver.velocity();
    double worst = 0.0;
    for (int i = 0; i < u.nx(); ++i) {
        for (int j = 0; j < u.ny(); ++j) {
            worst = std::max(worst, std::hypot(u(i, j, 0), u(i, j, 1)));
        }
    }
    return worst;
}

/// Run `config` for its step count, reporting the invariants under `prefix`.
void report_case(const std::string& prefix, const CaseConfig& config) {
    Solver solver(config);
    solver.initialize();

    const double mass_initial = total_mass(solver);
    double worst_phase = max_abs_phase(solver);
    for (int step = 0; step < config.steps; ++step) {
        solver.step();
        worst_phase = std::max(worst_phase, max_abs_phase(solver));
    }
    const double mass_final = total_mass(solver);

    std::cout.precision(17);
    std::cout << prefix << "_mass_initial = " << mass_initial << "\n"
              << prefix << "_mass_final = " << mass_final << "\n"
              << prefix << "_mass_drift = "
              << std::fabs(mass_final - mass_initial) / mass_initial << "\n"
              << prefix << "_max_abs_phase = " << worst_phase << "\n"
              << prefix << "_max_speed = " << max_speed(solver) << std::endl;
}

/// A lattice with no interface at all: phi = +1 everywhere, no gravity.
///
/// Nothing in the scheme should move it, so any velocity here is a spurious
/// force rather than physics.
void report_rest_state(const std::string& prefix, Boundary boundary) {
    CaseConfig config = small_case(boundary, GradientStencil::E4);
    config.steps = 10;
    config.initial_phase = [](const CaseConfig&, int, int) { return 1.0; };

    Solver solver(config);
    solver.initialize();
    const double mass_initial = total_mass(solver);
    for (int step = 0; step < config.steps; ++step) {
        solver.step();
    }

    std::cout.precision(17);
    std::cout << prefix << "_max_speed = " << max_speed(solver) << "\n"
              << prefix << "_mass_drift = "
              << std::fabs(total_mass(solver) - mass_initial) / mass_initial << std::endl;
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
        report_case("periodic", small_case(Boundary::PeriodicY, stencil));
        report_case("wall", small_case(Boundary::WallY, GradientStencil::E4));
        report_rest_state("rest_periodic", Boundary::PeriodicY);
        report_rest_state("rest_wall", Boundary::WallY);
    } catch (const std::exception& error) {
        std::cerr << "lbm_solver: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
