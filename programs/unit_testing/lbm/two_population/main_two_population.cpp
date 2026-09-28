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
//  6. a uniform fluid at rest staying at rest;
//  7. the MRT collision: BGK again when every rate is the viscous one, each
//     moment relaxed at its own rate otherwise, and mass and momentum
//     untouched by the relaxation;
//  8. the third-moment source: no mass, no momentum, and exactly the trace
//     and normal-stress difference Ba et al. Eqs. (17)-(18) prescribe;
//  9. items 3 to 5 again with the MRT collision and the source term on;
// 10. the decay rate of a Taylor-Green vortex in one fluid, whose strain is
//     purely normal and so measures exactly the stress the third-moment
//     source corrects, with and without it.
//
//   main_two_population [stencil]

#include <cmath>
#include <iostream>
#include <string>

#include "lbm/case_config.h"
#include "lbm/mrt.h"
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
///
/// `mrt` switches on the MRT collision and the third-moment source, with Ba et
/// al.'s *kinematic* viscosity in both fluids, which is the case they are for:
/// tau reaches 348 in the heavy fluid at a ratio of 1000.
void report_droplet(const std::string& prefix,
                    double ratio,
                    GradientStencil stencil,
                    bool mrt = false) {
    CaseConfig config = droplet_case(ratio, stencil);
    if (mrt) {
        config.collision = cglbm::lbm::Collision::MRT;
        config.third_moment_correction = true;
        config.physics.nu = 0.1667;
        config.physics.nu_b = 0.1667;
        config.physics.nu2 = 0.1667;
        config.physics.nu_b2 = 0.1667;
    }
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

/// A deterministic spread of values in [-1, 1], so the tests need no seed.
double pseudo_random(int n) {
    return std::sin(12.9898 * n + 78.233 * n * n);
}

/// The MRT kernel against BGK, and against its own definition.
///
/// With every rate equal to one `omega`, `mrt_collide` must be BGK with Guo's
/// forcing, `f - omega (f - f_eq) + (1 - omega/2) S dt`. With different rates,
/// each moment of the result must be `m - s (m - m_eq) + (1 - s/2) m_S dt`,
/// and the mass must change only by what the source carries.
void report_mrt() {
    using cglbm::lbm::kMoment;
    using cglbm::lbm::kMomentNorm;
    using cglbm::lbm::kQ;
    double worst_bgk = 0.0, worst_moment = 0.0, worst_inverse = 0.0;
    for (int trial = 0; trial < 20; ++trial) {
        // As in the solver: the equilibrium carries the populations' mass, and
        // Guo's source none. BGK would relax a mass difference; MRT does not.
        double f[kQ], f_eq[kQ], source[kQ];
        double mass = 0.0, mass_eq = 0.0, mass_source = 0.0;
        for (int i = 0; i < kQ; ++i) {
            f[i] = 0.1 + 0.05 * pseudo_random(100 * trial + i);
            f_eq[i] = 0.1 + 0.05 * pseudo_random(100 * trial + i + 30);
            source[i] = 1.0e-3 * pseudo_random(100 * trial + i + 60);
            mass += f[i];
            mass_eq += f_eq[i];
            mass_source += source[i];
        }
        for (int i = 0; i < kQ; ++i) {
            f_eq[i] += (mass - mass_eq) / kQ;
            source[i] -= mass_source / kQ;
        }
        const double omega = 0.3 + 1.5 * (0.5 + 0.5 * pseudo_random(trial + 7));

        double rates[kQ], out[kQ];
        cglbm::lbm::mrt_rates(omega, omega, omega, omega, rates);
        cglbm::lbm::mrt_collide(f, f_eq, source, rates, 1.0, out);
        for (int i = 0; i < kQ; ++i) {
            const double bgk = f[i] - omega * (f[i] - f_eq[i]) + (1.0 - 0.5 * omega) * source[i];
            worst_bgk = std::max(worst_bgk, std::fabs(out[i] - bgk));
        }

        cglbm::lbm::mrt_rates(omega, 1.25, 1.14, 1.6, rates);
        cglbm::lbm::mrt_collide(f, f_eq, source, rates, 1.0, out);
        for (int a = 0; a < kQ; ++a) {
            double m = 0.0, m_eq = 0.0, m_source = 0.0, m_out = 0.0;
            for (int i = 0; i < kQ; ++i) {
                m += kMoment[a][i] * f[i];
                m_eq += kMoment[a][i] * f_eq[i];
                m_source += kMoment[a][i] * source[i];
                m_out += kMoment[a][i] * out[i];
            }
            const double expected = m - rates[a] * (m - m_eq) + (1.0 - 0.5 * rates[a]) * m_source;
            worst_moment = std::max(worst_moment, std::fabs(m_out - expected));
        }
    }
    // M M^T = diag(kMomentNorm): the basis is orthogonal, so the transpose
    // over the norms is the inverse.
    for (int a = 0; a < kQ; ++a) {
        for (int b = 0; b < kQ; ++b) {
            double dot = 0.0;
            for (int i = 0; i < kQ; ++i) {
                dot += kMoment[a][i] * kMoment[b][i];
            }
            worst_inverse =
                std::max(worst_inverse, std::fabs(dot - (a == b ? kMomentNorm[a] : 0.0)));
        }
    }
    std::cout.precision(17);
    std::cout << "mrt_bgk_difference = " << worst_bgk << "\n"
              << "mrt_moment_error = " << worst_moment << "\n"
              << "mrt_orthogonality_error = " << worst_inverse << std::endl;
}

/// The third-moment source: what it adds to each moment.
///
/// It must add nothing to the mass or the momentum, `(1 - s_e/2) div Q` to the
/// trace of the second moment, `(1 - s_nu/2) (d_x Q_x - d_y Q_y)` to its normal
/// difference, and nothing to the shear stress.
void report_third_moment_source() {
    using cglbm::lbm::kQ;
    using cglbm::lbm::kXi;
    const double dqx_dx = 3.7e-4, dqy_dy = -1.1e-4, s_e = 1.25, s_nu = 0.8;
    double source[kQ] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    cglbm::lbm::add_third_moment_source(dqx_dx, dqy_dy, s_e, s_nu, 1.0, source);
    double mass = 0.0, jx = 0.0, jy = 0.0, trace = 0.0, normal = 0.0, shear = 0.0;
    for (int i = 0; i < kQ; ++i) {
        const double ex = kXi[i][0], ey = kXi[i][1];
        mass += source[i];
        jx += source[i] * ex;
        jy += source[i] * ey;
        trace += source[i] * (ex * ex + ey * ey);
        normal += source[i] * (ex * ex - ey * ey);
        shear += source[i] * ex * ey;
    }
    const double trace_exact = (1.0 - 0.5 * s_e) * (dqx_dx + dqy_dy);
    const double normal_exact = (1.0 - 0.5 * s_nu) * (dqx_dx - dqy_dy);
    std::cout.precision(17);
    std::cout << "third_moment_mass = " << std::fabs(mass) << "\n"
              << "third_moment_momentum = " << std::max(std::fabs(jx), std::fabs(jy)) << "\n"
              << "third_moment_trace_error = "
              << std::fabs(trace - trace_exact) / std::fabs(trace_exact) << "\n"
              << "third_moment_normal_error = "
              << std::fabs(normal - normal_exact) / std::fabs(normal_exact) << "\n"
              << "third_moment_shear = " << std::fabs(shear) << std::endl;
}

/// The MRT collision in the solver conserves mass and gives the momentum
/// exactly the force, node by node.
///
/// Run with the source term and Ba et al.'s equal kinematic viscosities, so
/// that tau reaches 348 in the droplet and every rate differs from every other.
/// With Guo's half-force velocity the force is `F = 2 (rho u - sum f e)`, so
/// the post-collision momentum must be `2 rho u - sum f e`.
void report_mrt_conservation() {
    CaseConfig config = droplet_case(1000.0, GradientStencil::E8);
    config.collision = cglbm::lbm::Collision::MRT;
    config.third_moment_correction = true;
    config.physics.nu = config.physics.nu_b = 0.1667;
    config.physics.nu2 = config.physics.nu_b2 = 0.1667;
    TwoPopulationSolver solver(config);
    solver.initialize();
    for (int step = 0; step < 20; ++step) {
        solver.step();
    }
    solver.refresh();
    double worst_mass = 0.0, worst_momentum = 0.0;
    for (int i = 0; i < config.nx; ++i) {
        for (int j = 0; j < config.ny; ++j) {
            double post[cglbm::lbm::kQ];
            solver.collide_node_for_test(i, j, post);
            double mass = 0.0, mass_post = 0.0, jx = 0.0, jy = 0.0, jx_post = 0.0, jy_post = 0.0;
            for (int k = 0; k < cglbm::lbm::kQ; ++k) {
                const double f = solver.population(0)(i, j, k) + solver.population(1)(i, j, k);
                mass += f;
                mass_post += post[k];
                jx += f * cglbm::lbm::kXi[k][0];
                jy += f * cglbm::lbm::kXi[k][1];
                jx_post += post[k] * cglbm::lbm::kXi[k][0];
                jy_post += post[k] * cglbm::lbm::kXi[k][1];
            }
            const double rho = solver.density()(i, j);
            const double expected_x = 2.0 * rho * solver.velocity()(i, j, 0) - jx;
            const double expected_y = 2.0 * rho * solver.velocity()(i, j, 1) - jy;
            worst_mass = std::max(worst_mass, std::fabs(mass_post - mass) / rho);
            worst_momentum = std::max(
                worst_momentum,
                std::max(std::fabs(jx_post - expected_x), std::fabs(jy_post - expected_y)) / rho);
        }
    }
    std::cout.precision(17);
    std::cout << "mrt_solver_mass_error = " << worst_mass << "\n"
              << "mrt_solver_momentum_error = " << worst_momentum << std::endl;
}

/// Decay rate of an axis-aligned Taylor-Green vortex in a uniform fluid 1,
/// over the rate `2 nu k^2` the Navier-Stokes equations give.
///
/// Its strain is `d_x u_x = -d_y u_y` with no shear at all, so the decay is
/// set entirely by the normal viscous stress -- the one D2Q9 gets wrong for a
/// fluid off the lattice sound speed. `shear` instead runs a shear wave
/// `u_x(y)`, which only the off-diagonal stress sees and which the enhanced
/// equilibrium already gets right.
double taylor_green_rate(double ratio, bool corrected, bool shear) {
    constexpr double pi = 3.14159265358979323846;
    CaseConfig config = droplet_case(ratio, GradientStencil::E4);
    config.nx = 32;
    config.ny = 32;
    config.physics.sigma = 0.0;
    config.initial_phase = [](const CaseConfig&, int, int) { return 1.0; };
    const double amplitude = 1.0e-4 / std::sqrt(ratio);
    config.initial_velocity = [amplitude,
                               shear](const CaseConfig& c, int i, int j, double* u_x, double* u_y) {
        const double k = 2.0 * pi / c.nx;
        const double x = i + 0.5, y = j + 0.5;
        *u_x = shear ? amplitude * std::sin(k * y) : amplitude * std::sin(k * x) * std::cos(k * y);
        *u_y = shear ? 0.0 : -amplitude * std::cos(k * x) * std::sin(k * y);
    };
    config.third_moment_correction = corrected;
    // tau = 0.8 in fluid 1: mu = (tau - 1/2) p with p = rho1 (c_s^1)^2.
    const double cs_squared = 0.6 * 0.8 / ratio;
    const double mu = 0.3 * ratio * cs_squared;
    config.physics.nu = config.physics.nu_b = mu / ratio;
    config.physics.nu2 = config.physics.nu_b2 = mu / ratio;

    TwoPopulationSolver solver(config);
    solver.initialize();
    auto kinetic = [&solver]() {
        double energy = 0.0;
        const Field& u = solver.velocity();
        for (int i = 0; i < u.nx(); ++i) {
            for (int j = 0; j < u.ny(); ++j) {
                energy += u(i, j, 0) * u(i, j, 0) + u(i, j, 1) * u(i, j, 1);
            }
        }
        return energy;
    };
    solver.refresh();
    const double start = kinetic();
    const int steps = 2000;
    for (int step = 0; step < steps; ++step) {
        solver.step();
    }
    solver.refresh();
    const double k = 2.0 * pi / config.nx;
    const double nu = mu / ratio;
    const double rate = -std::log(kinetic() / start) / (2.0 * steps);
    return rate / ((shear ? 1.0 : 2.0) * nu * k * k);
}

void report_taylor_green() {
    std::cout.precision(17);
    for (double ratio : {1.0, 10.0}) {
        const std::string prefix = ratio == 1.0 ? "tg_r1" : "tg_r10";
        std::cout << prefix << "_uncorrected = " << taylor_green_rate(ratio, false, false) << "\n"
                  << prefix << "_corrected = " << taylor_green_rate(ratio, true, false) << "\n";
    }
    std::cout << "tg_r1000_shear = " << taylor_green_rate(1000.0, false, true) << std::endl;
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
        report_mrt();
        report_third_moment_source();
        report_mrt_conservation();
        report_droplet("mrt_r1000", 1000.0, stencil, true);
        report_droplet("mrt_r100000", 100000.0, stencil, true);
        // Independent of the stencil argument -- the vortex has no interface
        // -- and the slowest item here, so it runs once, with E8.
        if (stencil == GradientStencil::E8) {
            report_taylor_green();
        }
    } catch (const std::exception& error) {
        std::cerr << "two_population: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
