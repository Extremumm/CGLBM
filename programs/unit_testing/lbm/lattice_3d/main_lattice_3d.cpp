// Checks the two three-dimensional lattices against the conditions their
// weights were solved for, and the equilibrium built on each against its
// moments.
//
// `programs/unit_testing/lbm/two_population_3d` already does this for D3Q19.
// What is here is what only exists once there are two lattices, and it is
// arranged around the one failure this design is meant to make impossible.
//
// The two-population model carries the density ratio in the rest weight, so
// four numbers change with the lattice: the rest weights of the moving shells,
// the sound speed `(c_s^k)^2`, and the amplitude `lambda` of the enhanced
// equilibrium. Three of them are visible in a table. The fourth is not: it is
// derived, it appears in no conservation law, and carrying it over from another
// lattice leaves mass, momentum and the second moment all exactly right while
// the shear viscosity is wrong by a factor that grows with the density ratio.
// That is not hypothetical -- it is what the port from D2Q9 to D3Q19 did, and
// the correct D3Q27 value is *precisely* the incorrect D3Q19 one, so the same
// slip in the other direction would look right.
//
// So `lambda` is never written down. The solver sums
//
//     S = sum_q w_q e_x^2 e_y^2,   T = sum_q w_q e_x^2 e_y^2 (3|e_q|^2 - 5)
//
// over whichever lattice it was handed and sets `lambda = (c_s^2 - 3S)/(3T)`,
// and what this program measures is the moment that fixes it, `M3_xxy`, on both
// lattices at four density ratios. It also reports `lambda` against its closed
// form for each, so that a reader can see the factor of two rather than infer
// it.
//
// Three more things are checked, all of them lattice-specific:
//
//  1. the weight conditions themselves -- `sum phi = 1`, the second moment
//     giving `(c_s^k)^2`, fourth-order isotropy on both lattices and the sixth-
//     order relation that closes D3Q27's third shell;
//  2. the velocity at which the equilibrium first goes negative, which is what
//     bounds the heavy fluid at a high density ratio. The prediction is
//     `(c_s^k)^2/3` on D3Q19 and `(c_s^k)^2/2` on D3Q27, and it is measured by
//     minimising over the directions rather than asserted;
//  3. that the solver actually runs on D3Q27 at a density ratio of 1000 --
//     mass conserved per fluid and the phase field in range. The spurious
//     velocity of both lattices is reported beside it, but it is *not* asserted
//     against each other: six hundred steps is the transient, not the relaxed
//     state, and measured there D3Q27 is marginally the worse of the two
//     (8.2e-3 against 7.8e-3). The comparison the lattice was added for is the
//     converged one, which is `laplace_3d` over 1.5e4 steps, and there D3Q27
//     wins: 6.05e-5 against 7.17e-5.
//
//   main_lattice_3d

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "lbm/case_config.h"
#include "lbm/two_population_solver_3d.h"

namespace {

using cglbm::lbm::Boundary;
using cglbm::lbm::CaseConfig;
using cglbm::lbm::Field3D;
using cglbm::lbm::GradientStencil3D;
using cglbm::lbm::Lattice3D;
using cglbm::lbm::Lattice3DKind;
using cglbm::lbm::TwoPopulationSolver3D;

const Lattice3DKind kLattices[2] = {Lattice3DKind::D3Q19, Lattice3DKind::D3Q27};

/// Largest velocity set on this lattice, so a stack buffer is always enough.
constexpr int kMaxQ = cglbm::lbm::kQ3D27;

CaseConfig droplet_case(double ratio, Lattice3DKind kind) {
    CaseConfig config;
    config.name = "lattice_3d_unit";
    config.nx = 32;
    config.ny = 32;
    config.nz = 32;
    config.steps = 30;
    config.interval = 1000000;

    config.physics.rho1 = ratio;
    config.physics.rho2 = 1.;
    config.physics.radius = 8.;
    config.physics.sigma = 0.1;
    const double mu = 0.1667;
    config.physics.nu = mu / config.physics.rho1;
    config.physics.nu_b = config.physics.nu;
    config.physics.nu2 = mu / config.physics.rho2;
    config.physics.nu_b2 = config.physics.nu2;
    config.physics.alpha2 = 0.2;
    config.physics.beta = 0.7;
    config.physics.ch_width_init = 1.1;

    config.boundary = Boundary::PeriodicY;
    config.initial_phase_3d = cglbm::lbm::droplet_interface_3d();
    config.stencil_3d = GradientStencil3D::E6;
    config.lattice_3d = kind;
    return config;
}

/// The conditions the weights and the rest weights were solved for.
///
/// `sum w = 1` and `sum w e e = c_s^2 I` are what makes the lattice a lattice;
/// the isotropy relations are the ones that fix the shells against each other,
/// and they are stated as ratios between moments rather than against `c_s^4`
/// because every component of these velocities is 0 or +-1, so `e_x^4 = e_x^2`
/// identically and the absolute form would only ever be true at `c_s^2 = 1/3`.
void report_lattice_conditions() {
    std::cout.precision(17);
    for (Lattice3DKind kind : kLattices) {
        const Lattice3D& lattice = cglbm::lbm::lattice_3d(kind);
        const std::string tag = std::string("lattice_") + lattice.name;

        double sum = 0, m2 = 0, m4 = 0, m22 = 0, m42 = 0, m222 = 0, cross = 0;
        for (int q = 0; q < lattice.q; ++q) {
            const double w = lattice.w[q];
            const double x = lattice.xi[q][0], y = lattice.xi[q][1], z = lattice.xi[q][2];
            sum += w;
            m2 += w * x * x;
            cross += w * x * y;
            m4 += w * x * x * x * x;
            m22 += w * x * x * y * y;
            m42 += w * x * x * x * x * y * y;
            m222 += w * x * x * y * y * z * z;
        }
        std::cout << tag << "_velocities = " << lattice.q << "\n"
                  << tag << "_weight_sum_error = " << std::fabs(sum - 1.0) << "\n"
                  << tag << "_weight_cs2_error = " << std::fabs(m2 - 1.0 / 3.0) << "\n"
                  << tag << "_weight_off_diagonal = " << std::fabs(cross) << "\n"
                  << tag << "_weight_fourth_error = " << std::fabs(m4 - 3.0 * m22) << "\n";
        // Sixth order is not reachable on either lattice as a whole. The one
        // relation D3Q27 does satisfy is the one that closed its third shell,
        // and D3Q19 -- having no third shell -- cannot.
        std::cout << tag << "_weight_sixth_partial_error = " << std::fabs(m42 - 3.0 * m222)
                  << std::endl;

        // The rest weights, at an alpha well away from the standard one so that
        // the relations are tested rather than the tabulated numbers.
        for (double alpha : {0.2, 0.8}) {
            double phi_sum = 0, p2 = 0, p4 = 0, p22 = 0, p_cross = 0;
            for (int q = 0; q < lattice.q; ++q) {
                const double phi = cglbm::lbm::rest_weight_3d(lattice, q, alpha);
                const double x = lattice.xi[q][0], y = lattice.xi[q][1];
                phi_sum += phi;
                p2 += phi * x * x;
                p_cross += phi * x * y;
                p4 += phi * x * x * x * x;
                p22 += phi * x * x * y * y;
            }
            const double cs2 = cglbm::lbm::sound_speed_squared_3d(lattice, alpha);
            const std::string sub = tag + "_alpha" + (alpha < 0.5 ? "02" : "08");
            std::cout << sub << "_rest_sum_error = " << std::fabs(phi_sum - 1.0) << "\n"
                      << sub << "_rest_cs2_error = " << std::fabs(p2 - cs2) << "\n"
                      << sub << "_rest_off_diagonal = " << std::fabs(p_cross) << "\n"
                      << sub << "_rest_fourth_error = " << std::fabs(p4 - 3.0 * p22) << "\n"
                      << sub << "_cs2 = " << cs2 << std::endl;
        }
    }
}

/// The equilibrium's moments on both lattices, and the amplitude behind the
/// third one.
void report_equilibrium() {
    std::cout.precision(17);
    for (Lattice3DKind kind : kLattices) {
        const Lattice3D& lattice = cglbm::lbm::lattice_3d(kind);
        const std::string tag = std::string("equilibrium_") + lattice.name;
        double worst_mass = 0, worst_momentum = 0, worst_stress = 0, worst_third = 0;
        double worst_ratio = 0, worst_amplitude = 0, worst_pressure = 0;

        for (double ratio : {1.0, 20.0, 1000.0, 100000.0}) {
            TwoPopulationSolver3D solver(droplet_case(ratio, kind));
            const double encoded = (1.0 - solver.alpha2()) / (1.0 - solver.alpha1());
            worst_ratio = std::max(worst_ratio, std::fabs(encoded - ratio) / ratio);

            // The two bulk pressures are identically equal whichever lattice
            // carries them: the coefficient in (c_s^k)^2 cancels.
            const double p1 = ratio * solver.sound_speed_squared(0);
            const double p2 = 1.0 * solver.sound_speed_squared(1);
            worst_pressure = std::max(worst_pressure, std::fabs(p1 - p2) / p2);

            for (int fluid = 0; fluid < 2; ++fluid) {
                const double cs_squared = solver.sound_speed_squared(fluid);
                const double rho_k = fluid == 0 ? ratio : 1.0;
                const double u[3] = {0.031, -0.017, 0.023};
                double eq[kMaxQ];
                solver.equilibrium_for_test(fluid, rho_k, u, eq);

                double mass = 0, momentum[3] = {0, 0, 0}, pxx = 0, pxy = 0, pzz = 0;
                double qxxy = 0, qzzy = 0;
                for (int q = 0; q < lattice.q; ++q) {
                    const double x = lattice.xi[q][0], y = lattice.xi[q][1], z = lattice.xi[q][2];
                    mass += eq[q];
                    momentum[0] += eq[q] * x;
                    momentum[1] += eq[q] * y;
                    momentum[2] += eq[q] * z;
                    pxx += eq[q] * x * x;
                    pzz += eq[q] * z * z;
                    pxy += eq[q] * x * y;
                    qxxy += eq[q] * x * x * y;
                    qzzy += eq[q] * z * z * y;
                }
                worst_mass = std::max(worst_mass, std::fabs(mass - rho_k) / rho_k);
                for (int a = 0; a < 3; ++a) {
                    worst_momentum =
                        std::max(worst_momentum, std::fabs(momentum[a] - rho_k * u[a]) / rho_k);
                }
                worst_stress =
                    std::max({worst_stress,
                              std::fabs(pxx - rho_k * (cs_squared + u[0] * u[0])) / rho_k,
                              std::fabs(pzz - rho_k * (cs_squared + u[2] * u[2])) / rho_k,
                              std::fabs(pxy - rho_k * u[0] * u[1]) / rho_k});
                // The moment the amplitude exists to set, and the only one that
                // sees it. Both mixed forms must come out the same, since the
                // correction is isotropic in the plane.
                const double target = rho_k * cs_squared * u[1];
                const double scale = rho_k * std::fabs(u[1]);
                worst_third = std::max({worst_third,
                                        std::fabs(qxxy - target) / scale,
                                        std::fabs(qzzy - target) / scale});

                // The closed form, for the record rather than for the solver:
                // (3 c_s^2 - 1) on D3Q19 and half of it on D3Q27.
                const double expected = kind == Lattice3DKind::D3Q27
                                            ? 0.5 * (3.0 * cs_squared - 1.0)
                                            : (3.0 * cs_squared - 1.0);
                worst_amplitude = std::max(
                    worst_amplitude,
                    std::fabs(solver.enhancement(fluid) - expected) / std::fabs(expected));
            }
        }
        std::cout << tag << "_mass_error = " << worst_mass << "\n"
                  << tag << "_momentum_error = " << worst_momentum << "\n"
                  << tag << "_stress_error = " << worst_stress << "\n"
                  << tag << "_third_moment_error = " << worst_third << "\n"
                  << tag << "_amplitude_error = " << worst_amplitude << "\n"
                  << tag << "_density_ratio_error = " << worst_ratio << "\n"
                  << tag << "_pressure_match_error = " << worst_pressure << std::endl;
    }
}

/// Where the equilibrium first goes negative, measured rather than asserted.
///
/// At small `u` direction `q` carries `rho_k [phi_q + 3 w_q (e_q.u) A_q]` with
/// `A_q = 1 + lambda (3|e_q|^2 - 5)`, so the tightest direction bounds `|u|` by
/// `phi_q / (3 w_q |e_q| |A_q|)`. In the limit the heavy fluid at a high
/// density ratio sits in, that is `(c_s^k)^2/3` on D3Q19 and `(c_s^k)^2/2` on
/// D3Q27; the ratio between them is what a run at a ratio of 1000 spends.
void report_positivity_bound() {
    std::cout.precision(17);
    const double ratio = 1000.0;
    for (Lattice3DKind kind : kLattices) {
        const Lattice3D& lattice = cglbm::lbm::lattice_3d(kind);
        TwoPopulationSolver3D solver(droplet_case(ratio, kind));
        const std::string tag = std::string("positivity_") + lattice.name;

        // Fluid 0 is the heavy one, which is the one with the small sound speed
        // and so the small moving weights.
        const double cs_squared = solver.sound_speed_squared(0);
        const double alpha = solver.alpha1();
        const double amplitude = solver.enhancement(0);

        double bound = std::numeric_limits<double>::infinity();
        for (int q = 1; q < lattice.q; ++q) {
            const double e_squared = cglbm::lbm::speed_squared_3d(lattice, q);
            const double factor = 1.0 + amplitude * (3.0 * e_squared - 5.0);
            if (factor == 0.0) {
                continue;  // this direction carries no momentum, so no bound
            }
            const double phi = cglbm::lbm::rest_weight_3d(lattice, q, alpha);
            bound = std::min(bound,
                             phi / (3.0 * lattice.w[q] * std::sqrt(e_squared) *
                                    std::fabs(factor)));
        }
        const double predicted =
            kind == Lattice3DKind::D3Q27 ? cs_squared / 2.0 : cs_squared / 3.0;
        std::cout << tag << "_bound = " << bound << "\n"
                  << tag << "_bound_over_cs2 = " << bound / cs_squared << "\n"
                  << tag << "_predicted = " << predicted << "\n"
                  << tag << "_bound_error = " << std::fabs(bound - predicted) / predicted
                  << std::endl;
    }
}

double component_mass(const TwoPopulationSolver3D& solver, int fluid) {
    const Field3D& density = solver.component_density(fluid);
    double total = 0.0;
    for (int i = 0; i < density.nx(); ++i) {
        for (int j = 0; j < density.ny(); ++j) {
            for (int k = 0; k < density.nz(); ++k) {
                total += density(i, j, k);
            }
        }
    }
    return total;
}

double max_speed(const TwoPopulationSolver3D& solver) {
    const Field3D& u = solver.velocity();
    double worst = 0.0;
    for (int i = 0; i < u.nx(); ++i) {
        for (int j = 0; j < u.ny(); ++j) {
            for (int k = 0; k < u.nz(); ++k) {
                worst = std::max(worst,
                                 std::sqrt(u(i, j, k, 0) * u(i, j, k, 0) +
                                           u(i, j, k, 1) * u(i, j, k, 1) +
                                           u(i, j, k, 2) * u(i, j, k, 2)));
            }
        }
    }
    return worst;
}

double max_abs_phase(const TwoPopulationSolver3D& solver) {
    const Field3D& phase = solver.phase();
    double worst = 0.0;
    for (int i = 0; i < phase.nx(); ++i) {
        for (int j = 0; j < phase.ny(); ++j) {
            for (int k = 0; k < phase.nz(); ++k) {
                worst = std::max(worst, std::fabs(phase(i, j, k)));
            }
        }
    }
    return worst;
}

/// A droplet at a density ratio of 1000, relaxed on each lattice.
///
/// Both lattices are run side by side and both numbers reported. Six hundred
/// steps is enough for the interface to settle into the profile the recolouring
/// holds and for the run to show whether it is sound; it is nowhere near enough
/// for the spurious current to have converged, so that number is a record here
/// and a result only in `laplace_3d`.
void report_high_ratio_run() {
    std::cout.precision(17);
    for (Lattice3DKind kind : kLattices) {
        CaseConfig config = droplet_case(1000.0, kind);
        config.steps = 600;
        TwoPopulationSolver3D solver(config);
        const std::string tag = std::string("droplet_") + cglbm::lbm::lattice_3d(kind).name;

        solver.initialize();
        const double mass1 = component_mass(solver, 0);
        const double mass2 = component_mass(solver, 1);
        for (int step = 0; step < config.steps; ++step) {
            solver.step();
        }
        solver.refresh();

        std::cout << tag << "_mass1_drift = "
                  << std::fabs(component_mass(solver, 0) - mass1) / mass1 << "\n"
                  << tag << "_mass2_drift = "
                  << std::fabs(component_mass(solver, 1) - mass2) / mass2 << "\n"
                  << tag << "_max_abs_phase = " << max_abs_phase(solver) << "\n"
                  << tag << "_spurious_velocity = " << max_speed(solver) << std::endl;
    }
}

}  // namespace

int main() {
    report_lattice_conditions();
    report_equilibrium();
    report_positivity_bound();
    report_high_ratio_run();
    return 0;
}
