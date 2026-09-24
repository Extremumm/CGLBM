// Checks src/lbm/velocity_based.
//
// Reported as `key = value` lines for the pytest beside this file:
//
//  1. the equilibria, the forcing and the collision have the moments the
//     scheme is built on, and the collision conserves P and u;
//  2. the carrier Gamma_i(u) of the phase field is non-negative at the speeds
//     the solver runs at;
//  3. the pressure correction turns a uniform pressure across a density jump of
//     1e4 into no force at all, exactly;
//  4. the memoryless phase transport keeps its tanh profile, conserves c and
//     stays within [0, 1].
//
//   main_lbm_velocity_based

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "lbm/isotropic_gradient.h"
#include "lbm/surface_force.h"
#include "lbm/velocity_based.h"

namespace {

namespace vb = cglbm::lbm::velocity_based;

const double kCs2 = vb::kSoundSpeedSquared;

struct Moments {
    double m0;
    double mx;
    double my;
    double mxx;
    double myy;
    double mxy;
};

Moments moments(const double* f) {
    Moments m = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int k = 0; k < vb::kQ; ++k) {
        const double ex = vb::kVelocity[k][0];
        const double ey = vb::kVelocity[k][1];
        m.m0 += f[k];
        m.mx += f[k] * ex;
        m.my += f[k] * ey;
        m.mxx += f[k] * ex * ex;
        m.myy += f[k] * ey * ey;
        m.mxy += f[k] * ex * ey;
    }
    return m;
}

void report_moments() {
    const double P = 0.037;
    const double ux = 0.021;
    const double uy = -0.013;
    double eq[vb::kQ];
    vb::hydrodynamic_equilibrium(P, ux, uy, eq);
    const Moments e = moments(eq);
    double err = std::fabs(e.m0 - P);
    err = std::max(err, std::fabs(e.mx - ux));
    err = std::max(err, std::fabs(e.my - uy));
    err = std::max(err, std::fabs(e.mxx - (P * kCs2 + ux * ux)));
    err = std::max(err, std::fabs(e.myy - (P * kCs2 + uy * uy)));
    err = std::max(err, std::fabs(e.mxy - ux * uy));
    std::cout << "equilibrium_moment_error = " << err << "\n";

    // forcing: no mass, first moment a, second moment u a + a u
    const double ax = 3.0e-4;
    const double ay = -1.7e-4;
    double s[vb::kQ];
    vb::forcing(ux, uy, ax, ay, s);
    const Moments f = moments(s);
    double ferr = std::fabs(f.m0);
    ferr = std::max(ferr, std::fabs(f.mx - ax));
    ferr = std::max(ferr, std::fabs(f.my - ay));
    ferr = std::max(ferr, std::fabs(f.mxx - 2.0 * ux * ax));
    ferr = std::max(ferr, std::fabs(f.myy - 2.0 * uy * ay));
    ferr = std::max(ferr, std::fabs(f.mxy - (ux * ay + uy * ax)));
    std::cout << "forcing_moment_error = " << ferr << "\n";

    // collision: an arbitrary state keeps its P and u (plus half the force)
    double state[vb::kQ];
    for (int k = 0; k < vb::kQ; ++k) {
        state[k] = eq[k] + 1.0e-3 * std::sin(1.3 * k + 0.4);
    }
    const Moments before = moments(state);
    double post[vb::kQ];
    vb::hydrodynamic_equilibrium(before.m0, before.mx, before.my, eq);
    vb::collide(state, eq, s, 0.73, 1.0, post);
    const Moments after = moments(post);
    double cerr = std::fabs(after.m0 - before.m0);
    cerr = std::max(cerr, std::fabs(after.mx - (before.mx + 0.5 * ax)));
    cerr = std::max(cerr, std::fabs(after.my - (before.my + 0.5 * ay)));
    std::cout << "collision_conservation_error = " << cerr << "\n";

    // with both relaxation times at 1 the non-equilibrium is removed entirely
    vb::collide(state, eq, s, 1.0, 1.0, post);
    double relax = 0.0;
    for (int k = 0; k < vb::kQ; ++k) {
        relax = std::max(relax, std::fabs(post[k] - (eq[k] + 0.5 * s[k])));
    }
    std::cout << "unit_tau_error = " << relax << "\n";

    // phase populations: moments c and c u + A n
    const double c = 0.3;
    const double nx = 0.6;
    const double ny = -0.8;
    const double width = 1.6;
    double phase[vb::kQ];
    vb::phase_populations(c, ux, uy, nx, ny, width, phase);
    const Moments ph = moments(phase);
    const double A = 0.5 * kCs2 * 2.0 * c * (1.0 - c) / width;
    double perr = std::fabs(ph.m0 - c);
    perr = std::max(perr, std::fabs(ph.mx - (c * ux + A * nx)));
    perr = std::max(perr, std::fabs(ph.my - (c * uy + A * ny)));
    std::cout << "phase_moment_error = " << perr << "\n";
}

/// Smallest carrier population over all directions for |u| up to `speed`.
void report_carrier() {
    for (double speed : {0.01, 0.05, 0.1}) {
        double lowest = 1.0;
        for (int n = 0; n < 72; ++n) {
            const double angle = 2.0 * M_PI * n / 72.0;
            double gamma[vb::kQ];
            vb::velocity_equilibrium(speed * std::cos(angle), speed * std::sin(angle), gamma);
            for (int k = 0; k < vb::kQ; ++k) {
                lowest = std::min(lowest, gamma[k]);
            }
        }
        std::cout << "carrier_min_" << static_cast<int>(std::lround(speed * 100)) << " = " << lowest
                  << "\n";
    }
}

/// A uniform pressure across a density jump of 1e4 must exert no force.
void report_pressure_correction() {
    const int n = 32;
    const double p = 0.04;
    std::vector<double> rho(static_cast<std::size_t>(n) * n);
    std::vector<double> P(rho.size());
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double r = std::hypot(i - 16.0, j - 16.0);
            const double c = 0.5 * (1.0 - std::tanh((r - 7.0) / 1.6));
            const std::size_t k = static_cast<std::size_t>(i) * n + j;
            rho[k] = 1.0 + c * (1.0e4 - 1.0);
            P[k] = p / (rho[k] * kCs2);
        }
    }
    // lattice pressure force per unit mass: -cs^2 grad_lat(P); corrected total
    // must vanish, since p is uniform
    double worst = 0.0;
    double scale = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double lx = 0.0;
            double ly = 0.0;
            for (int k = 1; k < vb::kQ; ++k) {
                const int ip = (i + vb::kVelocity[k][0] + n) % n;
                const int jp = (j + vb::kVelocity[k][1] + n) % n;
                lx +=
                    vb::kWeight[k] * vb::kVelocity[k][0] * P[static_cast<std::size_t>(ip) * n + jp];
                ly +=
                    vb::kWeight[k] * vb::kVelocity[k][1] * P[static_cast<std::size_t>(ip) * n + jp];
            }
            double ax = 0.0;
            double ay = 0.0;
            vb::pressure_correction(P.data(), rho.data(), n, n, i, j, &ax, &ay);
            // -cs^2 grad_lat P = -(lx, ly), since grad_lat carries 1/cs^2
            worst = std::max(worst, std::hypot(ax - lx, ay - ly));
            scale = std::max(scale, std::hypot(lx, ly));
        }
    }
    std::cout << "uniform_pressure_residual = " << worst / scale << "\n";

    // and a uniform density needs no correction
    std::fill(rho.begin(), rho.end(), 3.0);
    double uniform = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double ax = 0.0;
            double ay = 0.0;
            vb::pressure_correction(P.data(), rho.data(), n, n, i, j, &ax, &ay);
            uniform = std::max(uniform, std::hypot(ax, ay));
        }
    }
    std::cout << "uniform_density_correction = " << uniform << "\n";
}

/// The phase transport alone, at rest, on a slab: profile, mass, bounds.
void report_phase_transport() {
    const int nx = 64;
    const int ny = 4;
    const double width = 1.6;
    const std::size_t size = static_cast<std::size_t>(nx) * ny;
    std::vector<double> c(size);
    std::vector<double> psi(size);
    std::vector<double> next(size);
    std::vector<double> exact(size);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const double value = 0.5 * (1.0 - std::tanh((std::fabs(i - 32.0) - 10.0) / width));
            exact[static_cast<std::size_t>(i) * ny + j] = value;
            c[static_cast<std::size_t>(i) * ny + j] = value;
        }
    }
    double mass0 = 0.0;
    for (double v : c) {
        mass0 += v;
    }
    double lowest = 1.0;
    double highest = 0.0;
    for (int step = 0; step < 2000; ++step) {
        for (std::size_t k = 0; k < size; ++k) {
            psi[k] = 2.0 * c[k] - 1.0;
        }
        std::fill(next.begin(), next.end(), 0.0);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double gx = 0.0;
                double gy = 0.0;
                cglbm::lbm::gradient_periodic(
                    psi.data(), nx, ny, i, j, cglbm::lbm::GradientStencil::E8, &gx, &gy);
                double nxv = 0.0;
                double nyv = 0.0;
                cglbm::lbm::unit_normal(gx, gy, &nxv, &nyv);
                double pop[vb::kQ];
                vb::phase_populations(
                    c[static_cast<std::size_t>(i) * ny + j], 0.0, 0.0, nxv, nyv, width, pop);
                for (int k = 0; k < vb::kQ; ++k) {
                    const int ip = (i + vb::kVelocity[k][0] + nx) % nx;
                    const int jp = (j + vb::kVelocity[k][1] + ny) % ny;
                    next[static_cast<std::size_t>(ip) * ny + jp] += pop[k];
                }
            }
        }
        c.swap(next);
        for (double v : c) {
            lowest = std::min(lowest, v);
            highest = std::max(highest, v);
        }
    }
    double mass = 0.0;
    double deviation = 0.0;
    for (std::size_t k = 0; k < size; ++k) {
        mass += c[k];
        deviation = std::max(deviation, std::fabs(c[k] - exact[k]));
    }
    std::cout << "phase_mass_error = " << std::fabs(mass - mass0) / mass0 << "\n";
    std::cout << "phase_profile_deviation = " << deviation << "\n";
    std::cout << "phase_min = " << lowest << "\n";
    std::cout << "phase_max = " << highest << "\n";
}

}  // namespace

int main() {
    std::cout.precision(12);
    report_moments();
    report_carrier();
    report_pressure_correction();
    report_phase_transport();
    std::cout << std::flush;
    return 0;
}
