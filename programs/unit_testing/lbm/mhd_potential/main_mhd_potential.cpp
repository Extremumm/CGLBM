// Checks the inductionless MHD potential solve, both discretisations of it,
// and the current it produces.
//
// The coupling stands entirely on one elliptic solve, and it is worth being
// precise about what can be demanded of it. Five things are measured here, and
// they are deliberately of different kinds.
//
//  1. **The linear solve is exact.** A potential is chosen, the finite-volume
//     operator is applied to it to manufacture a right-hand side, and the
//     solver is asked to recover the potential. No truncation error enters:
//     the answer is the tolerance the conjugate gradient reached, and nothing
//     else. Done with a uniform conductivity and with one that jumps by a
//     factor of 10^4 across an interface, because the second is what the
//     preconditioner has to survive.
//
//  2. **Charge is conserved.** `sum_faces J_f` at every node, for a real
//     velocity field. This is the property the whole face-based arrangement
//     exists for, and it is not a truncation statement either: with the
//     finite-volume potential the face current is built from the same
//     difference the operator was, so the imbalance *is* the linear residual.
//
//  3. **The insulating wall holds.** With walls along y the current may not
//     cross them, and the flux form makes that exact rather than approximate.
//
//  4. **The discretisation is second order.** Against an analytic potential,
//     under refinement. This is the one measurement where truncation is the
//     answer rather than the noise, and it is the only one that says the
//     operator solves the right equation rather than solving some equation
//     well.
//
//  5. **The lattice-Boltzmann march reaches the same answer.** The D3Q7
//     pseudo-time model is a different discretisation, so it cannot agree to
//     round-off; what is measured is that it converges at all, how many sweeps
//     it takes from cold, and how far its steady state sits from the
//     finite-volume one. Both are reported rather than asserted against each
//     other, because neither is the reference for the other -- the analytic
//     solution in (4) is, and both are run against it.
//
//   main_mhd_potential

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "lbm/field3d.h"
#include "lbm/quasi_static_mhd_3d.h"

namespace {

using cglbm::lbm::Field3D;
using cglbm::lbm::MhdPhysics;
using cglbm::lbm::PotentialSolver;
using cglbm::lbm::QuasiStaticMhd3D;

constexpr double kPi = 3.14159265358979323846;

/// A conductivity that jumps across a diffuse interface at `x = nx/2`.
///
/// The profile is the tanh the colour-gradient solver maintains, so the
/// conductivity the module sees here is shaped like the one it sees in a run.
Field3D interface_phase(int nx, int ny, int nz, double width) {
    Field3D phase(nx, ny, nz);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                phase(i, j, k) = std::tanh((i - 0.5 * nx) / width);
            }
        }
    }
    return phase;
}

Field3D uniform_phase(int nx, int ny, int nz, double value) {
    Field3D phase(nx, ny, nz);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                phase(i, j, k) = value;
            }
        }
    }
    return phase;
}

/// A smooth, triply periodic, zero-mean potential to manufacture from.
Field3D manufactured_potential(int nx, int ny, int nz) {
    Field3D field(nx, ny, nz);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                field(i, j, k) = std::sin(2.0 * kPi * i / nx) * std::cos(2.0 * kPi * j / ny) +
                                 0.4 * std::cos(2.0 * kPi * k / nz);
            }
        }
    }
    double mean = 0.0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                mean += field(i, j, k);
            }
        }
    }
    mean /= static_cast<double>(nx) * ny * nz;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                field(i, j, k) -= mean;
            }
        }
    }
    return field;
}

double max_difference(const Field3D& a, const Field3D& b) {
    double worst = 0.0;
    for (int i = 0; i < a.nx(); ++i) {
        for (int j = 0; j < a.ny(); ++j) {
            for (int k = 0; k < a.nz(); ++k) {
                worst = std::max(worst, std::fabs(a(i, j, k) - b(i, j, k)));
            }
        }
    }
    return worst;
}

double max_abs(const Field3D& a) {
    double worst = 0.0;
    for (int i = 0; i < a.nx(); ++i) {
        for (int j = 0; j < a.ny(); ++j) {
            for (int k = 0; k < a.nz(); ++k) {
                worst = std::max(worst, std::fabs(a(i, j, k)));
            }
        }
    }
    return worst;
}

MhdPhysics base_physics() {
    MhdPhysics physics;
    physics.enabled = true;
    physics.b[0] = 0.0;
    physics.b[1] = 0.0;
    physics.b[2] = 0.05;
    physics.conductivity1 = 1.0;
    physics.conductivity2 = 1.0;
    physics.tolerance = 1e-12;
    physics.max_iterations = 4000;
    return physics;
}

/// (1) The linear solve, against a potential it was handed the residual of.
void report_manufactured_solve() {
    std::cout.precision(17);
    const int n = 24;
    struct Case {
        const char* tag;
        double sigma1;
        double sigma2;
        bool interface;
    };
    const Case cases[2] = {{"uniform", 1.0, 1.0, false}, {"jump", 1.0e4, 1.0, true}};

    for (const Case& item : cases) {
        MhdPhysics physics = base_physics();
        physics.conductivity1 = item.sigma1;
        physics.conductivity2 = item.sigma2;
        QuasiStaticMhd3D module(n, n, n, false, physics, false);

        const Field3D phase =
            item.interface ? interface_phase(n, n, n, 1.6) : uniform_phase(n, n, n, 1.0);
        // The conductivity is set by a solve; a zero field leaves the
        // coefficients in place and costs nothing.
        Field3D velocity(n, n, n, 3);
        module.solve(velocity, phase);

        const Field3D exact = manufactured_potential(n, n, n);
        Field3D rhs(n, n, n);
        module.apply_operator(exact, rhs);
        module.solve_potential(rhs);

        const std::string tag = std::string("manufactured_") + item.tag;
        std::cout << tag << "_error = " << max_difference(module.potential(), exact) / max_abs(exact)
                  << "\n"
                  << tag << "_iterations = " << module.iterations() << "\n"
                  << tag << "_residual = " << module.residual() << "\n"
                  << tag << "_converged = " << (module.converged() ? 1 : 0) << std::endl;
    }
}

/// A swirling, divergence-free velocity with structure on all three axes.
Field3D swirl(int nx, int ny, int nz, double amplitude) {
    Field3D velocity(nx, ny, nz, 3);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                const double x = 2.0 * kPi * i / nx;
                const double y = 2.0 * kPi * j / ny;
                const double z = 2.0 * kPi * k / nz;
                velocity(i, j, k, 0) = amplitude * std::sin(x) * std::cos(y) * std::cos(z);
                velocity(i, j, k, 1) = -amplitude * std::cos(x) * std::sin(y) * std::cos(z);
                velocity(i, j, k, 2) = amplitude * std::cos(x) * std::cos(y) * std::sin(z);
            }
        }
    }
    return velocity;
}

/// (2) and (3): charge conservation, periodic and against an insulating wall.
void report_charge_conservation() {
    std::cout.precision(17);
    const int n = 24;
    for (bool wall : {false, true}) {
        for (PotentialSolver solver :
             {PotentialSolver::FiniteVolume, PotentialSolver::LatticeBoltzmann}) {
            MhdPhysics physics = base_physics();
            physics.conductivity1 = 1.0e4;
            physics.conductivity2 = 1.0;
            physics.solver = solver;
            if (solver == PotentialSolver::LatticeBoltzmann) {
                physics.tolerance = 1e-9;
                physics.max_iterations = 40000;
            }
            QuasiStaticMhd3D module(n, n, n, wall, physics, false);

            const Field3D phase = interface_phase(n, n, n, 1.6);
            const Field3D velocity = swirl(n, n, n, 0.02);
            module.solve(velocity, phase);

            // The current the force is actually built from, for scale: an
            // imbalance is only meaningful beside the current carrying it.
            const Field3D& current = module.current();
            double largest = 0.0;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    for (int k = 0; k < n; ++k) {
                        largest = std::max(largest,
                                           std::sqrt(current(i, j, k, 0) * current(i, j, k, 0) +
                                                     current(i, j, k, 1) * current(i, j, k, 1) +
                                                     current(i, j, k, 2) * current(i, j, k, 2)));
                    }
                }
            }

            const std::string tag = std::string("charge_") + (wall ? "wall_" : "periodic_") +
                                    cglbm::lbm::potential_solver_name(solver);
            std::cout << tag << "_imbalance = " << module.charge_imbalance() << "\n"
                      << tag << "_current = " << largest << "\n"
                      << tag << "_relative = "
                      << (largest > 0.0 ? module.charge_imbalance() / largest : 0.0) << "\n"
                      << tag << "_iterations = " << module.iterations() << "\n"
                      << tag << "_converged = " << (module.converged() ? 1 : 0) << std::endl;
        }
    }
}

/// (4) and (5): both discretisations against an analytic potential.
///
/// With a uniform conductivity, `B = B z^` and
///
///     u = (U sin(k y), U sin(k x), 0),
///
/// the drive is `div (sigma u x B) = sigma B U k [cos(k x) - cos(k y)]`, whose
/// potential is
///
///     phi = -(U B / k) [cos(k x) - cos(k y)],
///
/// exactly. The lattice spacing is one, so refining the lattice refines `k`
/// with it and the error should fall as `n^-2`.
void report_convergence() {
    std::cout.precision(17);
    for (PotentialSolver solver :
         {PotentialSolver::FiniteVolume, PotentialSolver::LatticeBoltzmann}) {
        const char* name = cglbm::lbm::potential_solver_name(solver);
        double previous = 0.0;
        for (int n : {16, 32, 64}) {
            const int nz = 4;
            MhdPhysics physics = base_physics();
            physics.solver = solver;
            const double field = physics.b[2];
            if (solver == PotentialSolver::LatticeBoltzmann) {
                physics.tolerance = 1e-10;
                physics.max_iterations = 200000;
            }
            QuasiStaticMhd3D module(n, n, nz, false, physics, false);

            const double amplitude = 0.02;
            const double k = 2.0 * kPi / n;
            Field3D velocity(n, n, nz, 3);
            Field3D exact(n, n, nz);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    for (int kk = 0; kk < nz; ++kk) {
                        velocity(i, j, kk, 0) = amplitude * std::sin(k * j);
                        velocity(i, j, kk, 1) = amplitude * std::sin(k * i);
                        velocity(i, j, kk, 2) = 0.0;
                        exact(i, j, kk) = -(amplitude * field / k) *
                                          (std::cos(k * i) - std::cos(k * j));
                    }
                }
            }
            module.solve(velocity, uniform_phase(n, n, nz, 1.0));

            const double error = max_difference(module.potential(), exact) / max_abs(exact);
            const std::string tag = std::string("analytic_") + name + "_n" + std::to_string(n);
            std::cout << tag << "_error = " << error << "\n"
                      << tag << "_iterations = " << module.iterations() << "\n"
                      << tag << "_converged = " << (module.converged() ? 1 : 0) << std::endl;
            if (previous > 0.0) {
                std::cout << tag << "_order = " << std::log2(previous / error) << std::endl;
            }
            previous = error;
        }
    }
}

/// (5), directly: how far apart the two discretisations land.
void report_solver_agreement() {
    std::cout.precision(17);
    const int n = 24;
    Field3D potentials[2];
    int sweeps[2] = {0, 0};
    int converged[2] = {0, 0};
    int index = 0;
    for (PotentialSolver solver :
         {PotentialSolver::FiniteVolume, PotentialSolver::LatticeBoltzmann}) {
        MhdPhysics physics = base_physics();
        physics.solver = solver;
        if (solver == PotentialSolver::LatticeBoltzmann) {
            physics.tolerance = 1e-10;
            physics.max_iterations = 200000;
        }
        QuasiStaticMhd3D module(n, n, n, false, physics, false);
        module.solve(swirl(n, n, n, 0.02), uniform_phase(n, n, n, 1.0));
        potentials[index] = module.potential();
        sweeps[index] = module.iterations();
        converged[index] = module.converged() ? 1 : 0;
        ++index;
    }
    std::cout << "agreement_difference = "
              << max_difference(potentials[0], potentials[1]) / max_abs(potentials[0]) << "\n"
              << "agreement_fv_iterations = " << sweeps[0] << "\n"
              << "agreement_lbm_sweeps = " << sweeps[1] << "\n"
              << "agreement_fv_converged = " << converged[0] << "\n"
              << "agreement_lbm_converged = " << converged[1] << std::endl;
}

}  // namespace

int main() {
    report_manufactured_solve();
    report_charge_conservation();
    report_convergence();
    report_solver_agreement();
    return 0;
}
