// Checks the three-dimensional colour-gradient solver of
// src/lbm/two_population_solver_3d.h and the stencils it stands on.
//
// The invariants are those of the two-dimensional model -- each fluid's mass
// conserved separately, both distributions non-negative, a uniform fluid at
// rest -- plus three that only exist here:
//
//  1. the D3Q19 gradient stencils are isotropic to the order they claim, which
//     is what the shell weights were solved for;
//  2. the rest weights phi_i^k satisfy sum = 1 and give (c_s^k)^2 = (1-alpha)/2,
//     which is the three-dimensional analogue of the D2Q9 relation and fixes
//     how the density ratio enters;
//  3. the curvature of a sphere comes out as -2/R, not -1/R. Laplace's law
//     reads 2 sigma / R in three dimensions, and getting this wrong would show
//     up only as a factor of two in a validation run.
//
//   main_two_population_3d [stencil]

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "lbm/case_config.h"
#include "lbm/two_population_solver_3d.h"

namespace {

using cglbm::lbm::Boundary;
using cglbm::lbm::CaseConfig;
using cglbm::lbm::Field3D;
using cglbm::lbm::GradientStencil3D;
using cglbm::lbm::TwoPopulationSolver3D;

CaseConfig droplet_case(double ratio, GradientStencil3D stencil) {
    CaseConfig config;
    config.name = "two_population_3d_unit";
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
    config.stencil_3d = stencil;
    return config;
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

double min_population(const TwoPopulationSolver3D& solver) {
    double worst = 0.0;
    for (int fluid = 0; fluid < 2; ++fluid) {
        const Field3D& f = solver.population(fluid);
        for (int i = 0; i < f.nx(); ++i) {
            for (int j = 0; j < f.ny(); ++j) {
                for (int k = 0; k < f.nz(); ++k) {
                    for (int q = 0; q < cglbm::lbm::kQ3D; ++q) {
                        worst = std::min(worst, f(i, j, k, q));
                    }
                }
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

/// The stencil weights were solved for these conditions; check they hold.
void report_stencil_isotropy() {
    std::cout.precision(17);
    for (GradientStencil3D stencil : {GradientStencil3D::E4, GradientStencil3D::E6}) {
        int count = 0;
        const cglbm::lbm::StencilPoint3D* points = cglbm::lbm::stencil_points_3d(stencil, &count);
        double m2 = 0, m4 = 0, m22 = 0, m6 = 0, m42 = 0, m222 = 0, cross = 0;
        for (int p = 0; p < count; ++p) {
            const double w = points[p].weight;
            const double x = points[p].cx, y = points[p].cy, z = points[p].cz;
            m2 += w * x * x;
            cross += w * x * y;
            m4 += w * x * x * x * x;
            m22 += w * x * x * y * y;
            m6 += w * std::pow(x, 6);
            m42 += w * std::pow(x, 4) * y * y;
            m222 += w * x * x * y * y * z * z;
        }
        const std::string tag = std::string("stencil_") + cglbm::lbm::stencil_name_3d(stencil);
        std::cout << tag << "_points = " << count << "\n"
                  << tag << "_second_moment = " << m2 << "\n"
                  << tag << "_off_diagonal = " << std::fabs(cross) << "\n"
                  << tag << "_fourth_error = " << std::fabs(m4 - 3.0 * m22) << "\n"
                  << tag << "_sixth_error = "
                  << std::max(std::fabs(m6 - 15.0 * m222), std::fabs(m42 - 3.0 * m222))
                  << std::endl;
    }
}

/// The equilibrium's moments, and the rest weights behind them.
void report_equilibrium() {
    std::cout.precision(17);
    double worst_mass = 0, worst_momentum = 0, worst_stress = 0, worst_ratio = 0, worst_cs = 0;
    for (double ratio : {1.0, 20.0, 1000.0, 100000.0}) {
        CaseConfig config = droplet_case(ratio, GradientStencil3D::E6);
        TwoPopulationSolver3D solver(config);
        const double encoded = (1.0 - solver.alpha2()) / (1.0 - solver.alpha1());
        worst_ratio = std::max(worst_ratio, std::fabs(encoded - ratio) / ratio);

        for (int fluid = 0; fluid < 2; ++fluid) {
            const double alpha = fluid == 0 ? solver.alpha1() : solver.alpha2();
            const double cs_squared = 0.5 * (1.0 - alpha);  // the D3Q19 relation
            const double rho_k = fluid == 0 ? ratio : 1.0;
            const double u[3] = {0.031, -0.017, 0.023};
            double eq[cglbm::lbm::kQ3D];
            solver.equilibrium_for_test(fluid, rho_k, u, eq);

            double mass = 0, momentum[3] = {0, 0, 0}, pxx = 0, pxy = 0, pzz = 0;
            for (int q = 0; q < cglbm::lbm::kQ3D; ++q) {
                mass += eq[q];
                for (int a = 0; a < 3; ++a) {
                    momentum[a] += eq[q] * cglbm::lbm::kXi3D[q][a];
                }
                pxx += eq[q] * cglbm::lbm::kXi3D[q][0] * cglbm::lbm::kXi3D[q][0];
                pzz += eq[q] * cglbm::lbm::kXi3D[q][2] * cglbm::lbm::kXi3D[q][2];
                pxy += eq[q] * cglbm::lbm::kXi3D[q][0] * cglbm::lbm::kXi3D[q][1];
            }
            worst_mass = std::max(worst_mass, std::fabs(mass - rho_k) / rho_k);
            for (int a = 0; a < 3; ++a) {
                worst_momentum =
                    std::max(worst_momentum, std::fabs(momentum[a] - rho_k * u[a]) / rho_k);
            }
            worst_stress = std::max({worst_stress,
                                     std::fabs(pxx - rho_k * (cs_squared + u[0] * u[0])) / rho_k,
                                     std::fabs(pzz - rho_k * (cs_squared + u[2] * u[2])) / rho_k,
                                     std::fabs(pxy - rho_k * u[0] * u[1]) / rho_k});
            // sum_q phi_q^k must be 1, which the mass moment at u = 0 shows.
            double rest[cglbm::lbm::kQ3D];
            const double zero[3] = {0, 0, 0};
            solver.equilibrium_for_test(fluid, 1.0, zero, rest);
            double total = 0, second = 0;
            for (int q = 0; q < cglbm::lbm::kQ3D; ++q) {
                total += rest[q];
                second += rest[q] * cglbm::lbm::kXi3D[q][0] * cglbm::lbm::kXi3D[q][0];
            }
            worst_cs = std::max({worst_cs, std::fabs(total - 1.0), std::fabs(second - cs_squared)});
        }
    }
    std::cout << "equilibrium_mass_error = " << worst_mass << "\n"
              << "equilibrium_momentum_error = " << worst_momentum << "\n"
              << "equilibrium_stress_error = " << worst_stress << "\n"
              << "rest_weight_error = " << worst_cs << "\n"
              << "alpha_density_ratio_error = " << worst_ratio << std::endl;
}

/// The curvature of a sphere must be -2/R, not -1/R.
///
/// Built from an analytic phase field rather than a relaxed droplet, so that
/// what is measured is the discrete curvature operator alone.
void report_sphere_curvature() {
    const int n = 64;
    const double radius = 20.0, width = 2.0;
    std::vector<double> phase(static_cast<std::size_t>(n) * n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                const double dx = i - n / 2, dy = j - n / 2, dz = k - n / 2;
                const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
                phase[(static_cast<std::size_t>(i) * n + j) * n + k] =
                    -std::tanh((r - radius) / width);
            }
        }
    }
    std::vector<double> nx(phase.size()), ny(phase.size()), nz(phase.size());
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                double gx, gy, gz;
                cglbm::lbm::gradient_periodic_3d(
                    phase.data(), n, n, n, i, j, k, GradientStencil3D::E6, &gx, &gy, &gz);
                const double norm = std::sqrt(gx * gx + gy * gy + gz * gz);
                const std::size_t index = (static_cast<std::size_t>(i) * n + j) * n + k;
                if (norm > 1.0e-6) {
                    nx[index] = -gx / norm;
                    ny[index] = -gy / norm;
                    nz[index] = -gz / norm;
                }
            }
        }
    }
    // Evaluate on the +x ray, at the node nearest the interface.
    const int i = n / 2 + static_cast<int>(radius);
    double d[3][3];
    cglbm::lbm::gradient_periodic_3d(
        nx.data(), n, n, n, i, n / 2, n / 2, GradientStencil3D::E6, &d[0][0], &d[0][1], &d[0][2]);
    cglbm::lbm::gradient_periodic_3d(
        ny.data(), n, n, n, i, n / 2, n / 2, GradientStencil3D::E6, &d[1][0], &d[1][1], &d[1][2]);
    cglbm::lbm::gradient_periodic_3d(
        nz.data(), n, n, n, i, n / 2, n / 2, GradientStencil3D::E6, &d[2][0], &d[2][1], &d[2][2]);
    const std::size_t index = (static_cast<std::size_t>(i) * n + n / 2) * n + n / 2;
    const double normal[3] = {nx[index], ny[index], nz[index]};
    double divergence = 0.0, projection = 0.0;
    for (int a = 0; a < 3; ++a) {
        divergence += d[a][a];
        for (int b = 0; b < 3; ++b) {
            projection += normal[a] * normal[b] * d[b][a];
        }
    }
    const double curvature = -(divergence - projection);
    std::cout.precision(17);
    std::cout << "sphere_curvature = " << curvature << "\n"
              << "sphere_curvature_exact = " << -2.0 / radius << std::endl;
}

void report_droplet(const std::string& prefix, double ratio, GradientStencil3D stencil) {
    const CaseConfig config = droplet_case(ratio, stencil);
    TwoPopulationSolver3D solver(config);
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

void report_rest_state() {
    CaseConfig config = droplet_case(20.0, GradientStencil3D::E6);
    config.initial_phase_3d = [](const CaseConfig&, int, int, int) { return 1.0; };
    TwoPopulationSolver3D solver(config);
    solver.initialize();
    for (int step = 0; step < config.steps; ++step) {
        solver.step();
    }
    solver.refresh();
    std::cout.precision(17);
    std::cout << "rest_max_speed = " << max_speed(solver) << std::endl;
}

}  // namespace

int main(int argc, char** argv) {
    GradientStencil3D stencil = GradientStencil3D::E6;
    if (argc > 1 && !cglbm::lbm::stencil_from_name_3d(argv[1], &stencil)) {
        std::cerr << "Unknown gradient stencil '" << argv[1] << "'; expected E4 or E6."
                  << std::endl;
        return 2;
    }

    std::cout << "stencil = " << cglbm::lbm::stencil_name_3d(stencil) << std::endl;
    try {
        report_stencil_isotropy();
        report_equilibrium();
        report_sphere_curvature();
        report_droplet("r20", 20.0, stencil);
        report_droplet("r1000", 1000.0, stencil);
        report_rest_state();
    } catch (const std::exception& error) {
        std::cerr << "two_population_3d: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
