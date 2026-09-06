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
              << prefix << "_mass_drift = " << std::fabs(mass_final - mass_initial) / mass_initial
              << "\n"
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
              << prefix
              << "_mass_drift = " << std::fabs(total_mass(solver) - mass_initial) / mass_initial
              << std::endl;
}

/// Report the bulk-normalised phase field at the points its contract names.
///
/// phi_N must fix both bulks, be the identity at equal densities, and cross
/// zero where the two components are present in equal proportion of their own
/// bulk densities -- which is the point the raw colour field gets wrong, and
/// gets wrong by more the larger the density ratio.
void report_normalised_phase() {
    std::cout.precision(17);
    for (double ratio : {1.0, 20.0, 1000.0, 100000.0}) {
        const std::string tag = "phin_r" + std::to_string(static_cast<long>(ratio));
        std::cout << tag << "_at_plus_one = " << cglbm::lbm::normalised_phase(1.0, ratio, 1.0)
                  << "\n"
                  << tag << "_at_minus_one = " << cglbm::lbm::normalised_phase(-1.0, ratio, 1.0)
                  << "\n"
                  // The colour field's zero, which should be the interface only
                  // when the densities are equal.
                  << tag << "_at_zero = " << cglbm::lbm::normalised_phase(0.0, ratio, 1.0)
                  << "\n"
                  // Where phi_N crosses zero, i.e. the true interface.
                  << tag << "_at_interface = "
                  << cglbm::lbm::normalised_phase((ratio - 1.0) / (ratio + 1.0), ratio, 1.0)
                  << "\n"
                  // Overshoot must be clamped, not propagated.
                  << tag << "_clamped_above = " << cglbm::lbm::normalised_phase(1.5, ratio, 1.0)
                  << std::endl;
    }
    // Identity at equal densities, checked across the range rather than at a point.
    double worst_identity = 0.0;
    for (int n = -10; n <= 10; ++n) {
        const double phi = 0.1 * n;
        worst_identity =
            std::max(worst_identity, std::fabs(cglbm::lbm::normalised_phase(phi, 3.0, 3.0) - phi));
    }
    std::cout << "phin_identity_error = " << worst_identity << std::endl;
}

/// The configuration the shipped Laplace case runs: the interface located by
/// the bulk-normalised phase field and the tension applied as a body force
/// built from an explicit curvature.
CaseConfig ba_case(GradientStencil stencil) {
    CaseConfig config = small_case(Boundary::PeriodicY, stencil);
    config.interface_field = cglbm::lbm::InterfaceField::BulkNormalised;
    config.surface_tension = cglbm::lbm::SurfaceTension::ContinuumSurfaceForce;
    return config;
}

/// The Latva-Kokko segregation operator, at a density ratio it can hold.
///
/// It is unusable above a ratio of about 10 -- see docs/numerics.md -- so this
/// runs at 2, where it is the better of the two operators. The point here is
/// only that the code path conserves mass and keeps phi in range.
CaseConfig latva_kokko_case(GradientStencil stencil) {
    CaseConfig config = ba_case(stencil);
    config.recolouring = cglbm::lbm::Recolouring::LatvaKokko;
    config.physics.rho1 = 2.;
    config.physics.rho2 = 1.;
    config.physics.p1_inf = cglbm::lbm::matched_p1_inf(config.physics);
    return config;
}

/// A case whose two components have different kinematic viscosities.
///
/// `nu2 = nu * rho1 / rho2` matches the *dynamic* viscosities, which is what
/// makes tau uniform across the interface instead of spanning the density
/// ratio. The point of running it here is that it takes the interpolation
/// branch in `Solver::viscosity_at`, which no other case does.
CaseConfig viscosity_ratio_case(GradientStencil stencil) {
    CaseConfig config = ba_case(stencil);
    config.physics.nu2 = config.physics.nu * config.physics.rho1 / config.physics.rho2;
    config.physics.nu_b2 = config.physics.nu_b * config.physics.rho1 / config.physics.rho2;
    return config;
}

/// Check that this equilibrium already carries the "enhanced" third-order term.
///
/// Leclaire et al. (2013) added a term to the colour-gradient equilibrium to
/// fix the third-order velocity moment, which a two-component lattice gets
/// wrong whenever the components' sound speeds differ; Ba et al. (2016) restate
/// it as their Eq. (14) and build their high-density-ratio model on it. The
/// equilibrium here comes from a Hermite expansion instead and never mentions
/// alpha, so whether it contains the same correction is not obvious by
/// inspection. It does, identically. Writing Ba's parameters as
///
///     (c_s^k)^2 = 3/5 (1 - alpha_k),   p_k = rho_k (c_s^k)^2
///
/// their extra term over the standard equilibrium is
///
///     rho_k W_i (3 e_i.u) 1/2 (3 (c_s^k)^2 - 1) (3 |e_i|^2 - 4)
///       = (p_k - rho_k cs^2) W_i (e_i.u) 4.5 (3 |e_i|^2 - 4)
///
/// and the third-order Hermite term of `Solver::equilibrium` reduces to the
/// same thing, because `H_xxx + H_yyx = e_x (|e|^2 - 4 cs^2)` and
/// `1 / (2 cs^6) = 13.5`. This reports the largest disagreement over a sweep of
/// sound speed, density and velocity; it is a round-off number or the two have
/// diverged.
void report_enhanced_equilibrium() {
    const double cs2 = 1. / 3.;
    const double cs6 = cs2 * cs2 * cs2;
    double worst = 0.;
    for (double csk_squared : {0.05, 0.16, 1. / 3., 0.45, 0.58}) {
        for (double rho : {1., 7.3, 1000.}) {
            for (double u_x : {-0.11, 0., 0.07}) {
                for (double u_y : {0.03, -0.19}) {
                    const double p = rho * csk_squared;
                    const double alpha = 1. - (5. / 3.) * csk_squared;
                    for (int k = 0; k < cglbm::lbm::kQ; k++) {
                        const double x = cglbm::lbm::kXi[k][0];
                        const double y = cglbm::lbm::kXi[k][1];
                        const double Hxxy = x * x * y - cs2 * y;
                        const double Hyyx = y * y * x - cs2 * x;
                        const double Hxxx = std::pow(x, 3) - cs2 * 3. * x;
                        const double Hyyy = std::pow(y, 3) - cs2 * 3. * y;
                        const double ours = (p - rho * cs2) * cglbm::lbm::kW[k] *
                                            (u_x * (Hyyx + Hxxx) + u_y * (Hyyy + Hxxy)) /
                                            (2. * cs6);
                        const double theirs = rho * cglbm::lbm::kW[k] * 3. * (x * u_x + y * u_y) *
                                              0.5 * (3. * (0.6 * (1. - alpha)) - 1.) *
                                              (3. * (x * x + y * y) - 4.);
                        worst = std::max(worst, std::fabs(ours - theirs));
                    }
                }
            }
        }
    }
    std::cout.precision(17);
    std::cout << "enhanced_equilibrium_difference = " << worst << std::endl;
}

/// Run a droplet at a series of density ratios and report whether it survives.
///
/// The scheme used to diverge before step 200 at any ratio above about 100,
/// because the interface was started out of mechanical equilibrium and the
/// resulting pressure discontinuity drove a transient that reached Mach 1.4.
/// Started in equilibrium it survives the opening transient at 10^5. This is
/// the guard on *that*, and on nothing more: a few hundred steps is a smoke
/// test, not a statement about the density ratio the scheme can hold. Run long
/// enough, 10^3 and above diverge -- see docs/numerics.md. A regression in the
/// initialisation shows up here in a fraction of a second rather than in a
/// validation run.
void report_density_ratios(int steps) {
    std::cout.precision(17);
    for (double ratio : {1.e3, 1.e5}) {
        CaseConfig config = small_case(Boundary::PeriodicY, GradientStencil::E8);
        config.nx = 64;
        config.ny = 64;
        config.physics.radius = 8.;
        config.physics.rho1 = ratio;
        config.physics.rho2 = 1.;
        config.physics.p1_inf = cglbm::lbm::matched_p1_inf(config.physics);
        config.steps = steps;

        Solver solver(config);
        solver.initialize();
        bool finite = true;
        for (int step = 0; step < steps && finite; ++step) {
            solver.step();
            const cglbm::lbm::Field& rho = solver.density();
            for (int i = 0; i < rho.nx() && finite; ++i) {
                for (int j = 0; j < rho.ny(); ++j) {
                    if (!std::isfinite(rho(i, j)) || rho(i, j) <= 0.0) {
                        finite = false;
                        break;
                    }
                }
            }
        }
        const std::string tag =
            "ratio_1e" + std::to_string(static_cast<int>(std::lround(std::log10(ratio))));
        std::cout << tag << "_finite = " << (finite ? 1 : 0) << "\n";
        if (finite) {
            std::cout << tag << "_max_speed = " << max_speed(solver) << "\n"
                      << tag << "_max_abs_phase = " << max_abs_phase(solver) << "\n"
                      << tag << "_mass_drift = 0" << std::endl;
        } else {
            std::cout << tag << "_max_speed = inf\n"
                      << tag << "_max_abs_phase = inf\n"
                      << tag << "_mass_drift = inf" << std::endl;
        }
    }
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
        // The same invariants on the Ba et al. configuration the Laplace case
        // runs, and on the Latva-Kokko segregation operator. Neither shares a
        // code path with the two above: the tension is a body force rather than
        // Omega^(2), and the recolouring reads rho rather than p.
        report_case("csf", ba_case(stencil));
        report_case("latva_kokko", latva_kokko_case(stencil));
        report_case("viscosity_ratio", viscosity_ratio_case(stencil));
        report_rest_state("rest_periodic", Boundary::PeriodicY);
        report_rest_state("rest_wall", Boundary::WallY);
        report_normalised_phase();
        report_enhanced_equilibrium();
        report_density_ratios(300);
    } catch (const std::exception& error) {
        std::cerr << "lbm_solver: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
