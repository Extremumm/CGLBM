// Checks src/lbm/central_moments.
//
// Reported as `key = value` lines for the pytest beside this file:
//
//  1. central_moments and populations_from_central_moments are each other's
//     inverse;
//  2. generalized_equilibrium, built from its central moments, is the Hermite
//     form g_eq,4 + (p - rho cs^2)(E + Phi) derived from Saito et al. (2023),
//     at the ideal-gas pressure and at a heavy fluid's, rho = 1000, p = 1/3;
//  3. the equilibrium of `Solver` differs from it only in k_21, k_12 and k_22,
//     by the terms central_moments.h lists;
//  4. collide_central_moment conserves mass and momentum, relaxes the shear
//     and the trace at their rates, sets the higher moments to equilibrium,
//     and at unit rates returns the generalized equilibrium;
//  5. a uniform fluid at rho = 20, p = 1/3 and tau = 100, noise 1e-6, stays
//     quiet with the trace at rate 1, and blows up with it at 1 / tau.
//
//   main_lbm_central_moments

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "lbm/central_moments.h"
#include "lbm/d2q9.h"

namespace {

namespace lbm = cglbm::lbm;

const double kCs2 = 1.0 / 3.0;

/// The Hermite form of the generalized equilibrium; see central_moments.h.
void hermite_equilibrium(double rho, double ux, double uy, double p, bool generalized, double* f) {
    for (int i = 0; i < lbm::kQ; ++i) {
        const double x = lbm::kXi[i][0];
        const double y = lbm::kXi[i][1];
        const double w = lbm::kW[i];
        const double hxx = x * x - kCs2;
        const double hyy = y * y - kCs2;
        const double hxy = x * y;
        const double hxxy = hxx * y;
        const double hxyy = x * hyy;
        const double hxxyy = hxx * hyy;
        const double cs4 = kCs2 * kCs2;
        const double cs6 = cs4 * kCs2;
        const double cs8 = cs6 * kCs2;
        double g = 1.0 + (ux * x + uy * y) / kCs2 +
                   (ux * ux * hxx + uy * uy * hyy + 2.0 * ux * uy * hxy) / (2.0 * cs4);
        double phi = (ux * hxyy + uy * hxxy) / (2.0 * cs6);
        if (generalized) {
            g += (ux * ux * uy * hxxy + ux * uy * uy * hxyy) / (2.0 * cs6) +
                 ux * ux * uy * uy * hxxyy / (4.0 * cs8);
            phi += (ux * ux + uy * uy) * hxxyy / (4.0 * cs8);
        }
        const double e = (hxx + hyy) / (2.0 * cs4) - hxxyy / (4.0 * cs6);
        f[i] = rho * w * g + (p - rho * kCs2) * w * (e + phi);
    }
}

double max_difference(const double* a, const double* b) {
    double largest = 0.0;
    for (int i = 0; i < lbm::kQ; ++i) {
        largest = std::max(largest, std::fabs(a[i] - b[i]));
    }
    return largest;
}

void report_round_trip() {
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> population(0.01, 1.0);
    std::uniform_real_distribution<double> speed(-0.3, 0.3);
    double error = 0.0;
    for (int trial = 0; trial < 100; ++trial) {
        double f[lbm::kQ], back[lbm::kQ], k[3][3];
        for (double& v : f) {
            v = population(rng);
        }
        const double ux = speed(rng);
        const double uy = speed(rng);
        lbm::central_moments(f, ux, uy, k);
        lbm::populations_from_central_moments(k, ux, uy, back);
        error = std::max(error, max_difference(f, back));
    }
    std::cout << "round_trip_error = " << error << "\n";
}

void report_equilibrium() {
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> speed(-0.2, 0.2);
    double hermite_error = 0.0;
    double moment_error = 0.0;
    double solver_low_order_error = 0.0;
    double solver_k22_offset_error = 0.0;
    for (int trial = 0; trial < 100; ++trial) {
        const bool heavy = trial % 2 == 1;
        const double rho = heavy ? 1000.0 : 1.3;
        const double p = heavy ? 1.0 / 3.0 : 1.3 * kCs2;
        const double ux = speed(rng);
        const double uy = speed(rng);
        double built[lbm::kQ], hermite[lbm::kQ], solver[lbm::kQ];
        lbm::generalized_equilibrium(rho, ux, uy, p, built);
        hermite_equilibrium(rho, ux, uy, p, true, hermite);
        hermite_equilibrium(rho, ux, uy, p, false, solver);
        hermite_error = std::max(hermite_error, max_difference(built, hermite) / rho);

        double k[3][3], ks[3][3];
        lbm::central_moments(built, ux, uy, k);
        const double expected[3][3] = {{rho, 0, p}, {0, 0, 0}, {p, 0, p * kCs2}};
        lbm::central_moments(solver, ux, uy, ks);
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                moment_error = std::max(moment_error, std::fabs(k[a][b] - expected[a][b]) / rho);
                if (a + b <= 2) {
                    solver_low_order_error =
                        std::max(solver_low_order_error, std::fabs(ks[a][b] - k[a][b]) / rho);
                }
            }
        }
        // what the Solver's equilibrium keeps: the rho u^2 u terms of g_eq,2 in
        // k_21 and k_12, and in k_22 -(p - rho cs^2)|u|^2 + 3 rho u_x^2 u_y^2
        const double k22 =
            p * kCs2 - (p - rho * kCs2) * (ux * ux + uy * uy) + 3.0 * rho * ux * ux * uy * uy;
        solver_k22_offset_error = std::max({solver_k22_offset_error,
                                            std::fabs(ks[2][2] - k22) / rho,
                                            std::fabs(ks[2][1] + rho * ux * ux * uy) / rho,
                                            std::fabs(ks[1][2] + rho * ux * uy * uy) / rho});
    }
    std::cout << "hermite_equilibrium_error = " << hermite_error << "\n";
    std::cout << "equilibrium_moment_error = " << moment_error << "\n";
    std::cout << "solver_equilibrium_low_order_error = " << solver_low_order_error << "\n";
    std::cout << "solver_equilibrium_k22_error = " << solver_k22_offset_error << "\n";
}

void report_collision() {
    std::mt19937 rng(13);
    std::uniform_real_distribution<double> population(0.02, 0.2);
    std::uniform_real_distribution<double> rate(0.1, 1.9);
    double conservation = 0.0;
    double relaxation = 0.0;
    double unit_rate = 0.0;
    for (int trial = 0; trial < 100; ++trial) {
        double f[lbm::kQ];
        for (double& v : f) {
            v = population(rng);
        }
        double rho = 0.0, jx = 0.0, jy = 0.0;
        for (int i = 0; i < lbm::kQ; ++i) {
            rho += f[i];
            jx += f[i] * lbm::kXi[i][0];
            jy += f[i] * lbm::kXi[i][1];
        }
        const double ux = jx / rho;
        const double uy = jy / rho;
        const double p = 0.9 * rho * kCs2;
        const double omega_shear = rate(rng);
        const double omega_bulk = rate(rng);
        double post[lbm::kQ];
        lbm::collide_central_moment(f, ux, uy, p, omega_shear, omega_bulk, post);
        double rho_post = 0.0, jx_post = 0.0, jy_post = 0.0;
        for (int i = 0; i < lbm::kQ; ++i) {
            rho_post += post[i];
            jx_post += post[i] * lbm::kXi[i][0];
            jy_post += post[i] * lbm::kXi[i][1];
        }
        conservation = std::max({conservation,
                                 std::fabs(rho_post - rho),
                                 std::fabs(jx_post - jx),
                                 std::fabs(jy_post - jy)});
        double k[3][3], kp[3][3];
        lbm::central_moments(f, ux, uy, k);
        lbm::central_moments(post, ux, uy, kp);
        relaxation =
            std::max({relaxation,
                      std::fabs((kp[2][0] - kp[0][2]) - (1.0 - omega_shear) * (k[2][0] - k[0][2])),
                      std::fabs(kp[1][1] - (1.0 - omega_shear) * k[1][1]),
                      std::fabs((kp[2][0] + kp[0][2]) -
                                ((1.0 - omega_bulk) * (k[2][0] + k[0][2]) + omega_bulk * 2.0 * p)),
                      std::fabs(kp[2][1]),
                      std::fabs(kp[1][2]),
                      std::fabs(kp[2][2] - p * kCs2)});
        double unit[lbm::kQ], equilibrium[lbm::kQ];
        lbm::collide_central_moment(f, ux, uy, p, 1.0, 1.0, unit);
        lbm::generalized_equilibrium(rho, ux, uy, p, equilibrium);
        unit_rate = std::max(unit_rate, max_difference(unit, equilibrium));
    }
    std::cout << "collision_conservation_error = " << conservation << "\n";
    std::cout << "collision_relaxation_error = " << relaxation << "\n";
    std::cout << "collision_unit_rate_error = " << unit_rate << "\n";
}

/// Largest speed after `steps` steps of a uniform periodic fluid, 16 x 16,
/// started from noise of 1e-6, under collide_central_moment and streaming.
double uniform_fluid_speed(double omega_bulk, int steps) {
    const int n = 16;
    const double rho0 = 20.0;
    const double p_inf = rho0 * kCs2 - 1.0 / 3.0;
    const double tau = 100.0;
    std::mt19937 rng(17);
    std::normal_distribution<double> noise(0.0, 1.0e-6);
    std::vector<double> f(n * n * lbm::kQ), next(n * n * lbm::kQ);
    for (int m = 0; m < n * n; ++m) {
        const double rho = rho0 * (1.0 + noise(rng));
        lbm::generalized_equilibrium(
            rho, noise(rng), noise(rng), rho * kCs2 - p_inf, &f[m * lbm::kQ]);
    }
    double speed = 0.0;
    for (int step = 0; step < steps; ++step) {
        speed = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                const double* node = &f[(i * n + j) * lbm::kQ];
                double rho = 0.0, jx = 0.0, jy = 0.0;
                for (int k = 0; k < lbm::kQ; ++k) {
                    rho += node[k];
                    jx += node[k] * lbm::kXi[k][0];
                    jy += node[k] * lbm::kXi[k][1];
                }
                const double ux = jx / rho;
                const double uy = jy / rho;
                speed = std::max(speed, std::hypot(ux, uy));
                double post[lbm::kQ];
                lbm::collide_central_moment(
                    node, ux, uy, rho * kCs2 - p_inf, 1.0 / tau, omega_bulk, post);
                for (int k = 0; k < lbm::kQ; ++k) {
                    const int ip = (i + static_cast<int>(lbm::kXi[k][0]) + n) % n;
                    const int jp = (j + static_cast<int>(lbm::kXi[k][1]) + n) % n;
                    next[(ip * n + jp) * lbm::kQ + k] = post[k];
                }
            }
        }
        f.swap(next);
        if (!std::isfinite(speed) || speed > 1.0) {
            return 1.0;
        }
    }
    return speed;
}

void report_uniform_stability() {
    std::cout << "uniform_speed_unit_bulk = " << uniform_fluid_speed(1.0, 400) << "\n";
    std::cout << "uniform_speed_slow_bulk = " << uniform_fluid_speed(0.01, 400) << "\n";
}

}  // namespace

int main() {
    std::cout.precision(12);
    report_round_trip();
    report_equilibrium();
    report_collision();
    report_uniform_stability();
    std::cout << std::flush;
    return 0;
}
