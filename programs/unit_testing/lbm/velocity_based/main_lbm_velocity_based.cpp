// Checks src/lbm/velocity_based.
//
// Reported as `key = value` lines for the pytest beside this file:
//
//  1. the equilibria, the forcing and the collision have the moments the
//     scheme is built on, and the collision conserves P and u;
//  2. the carrier Gamma_i(u) of the phase field is non-negative at the speeds
//     the solver runs at;
//  3. the pressure force turns a uniform pressure across a density jump of 1e4
//     into no force at all, and sums to zero as a force density;
//  4. the link momentum exchange is equal and opposite at the two ends of a
//     link, with the same dissipation coefficient at both, sums to zero over a
//     lattice, carries a uniform velocity with the mass alone, and is the
//     lattice's own where there is a single component; the dissipation force
//     sums to zero and leaves a uniform velocity alone;
//  5. the hybrid collision reduces to the regularised one, conserves P and u,
//     and leaves a resolved non-equilibrium alone; set_velocity changes the
//     first moment only;
//  6. the memoryless phase transport keeps its tanh profile, conserves c and
//     stays within [0, 1]; at |u| = 0.2 the limited sharpening keeps every
//     phase population between 0 and its carrier.
//
//   main_lbm_velocity_based

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
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

/// A uniform pressure across a density jump of 1e4 must exert no force, and
/// the force density must sum to zero whatever the fields.
void report_pressure_force() {
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
    // scale: the force a single light node would feel from the jump in P
    double worst = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double ax = 0.0;
            double ay = 0.0;
            vb::pressure_force(P.data(), rho.data(), n, n, i, j, &ax, &ay);
            worst = std::max(worst, std::hypot(ax, ay));
        }
    }
    std::cout << "uniform_pressure_force = " << worst / (p / kCs2) << "\n";

    std::mt19937 random(7);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    double total_x = 0.0;
    double total_y = 0.0;
    double scale = 0.0;
    for (std::size_t k = 0; k < rho.size(); ++k) {
        rho[k] = 1.0 + 1.0e4 * unit(random);
        P[k] = unit(random) / rho[k];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double ax = 0.0;
            double ay = 0.0;
            vb::pressure_force(P.data(), rho.data(), n, n, i, j, &ax, &ay);
            const double r = rho[static_cast<std::size_t>(i) * n + j];
            total_x += r * ax;
            total_y += r * ay;
            scale = std::max(scale, r * std::hypot(ax, ay));
        }
    }
    std::cout << "pressure_force_total = " << std::hypot(total_x, total_y) / scale << "\n";
}

vb::LinkEnd random_end(std::mt19937& random, double rho1, double rho2) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    vb::LinkEnd end;
    const double c = unit(random);
    end.outgoing = 0.05 * (unit(random) - 0.5);
    end.phase = c * 0.1 * unit(random);
    end.rho = rho2 + c * (rho1 - rho2);
    end.mu = 0.1 + unit(random);
    end.ux = 0.02 * (unit(random) - 0.5);
    end.uy = 0.02 * (unit(random) - 0.5);
    return end;
}

void report_link_momentum() {
    const double rho1 = 1.0e4;
    const double rho2 = 1.0;
    std::mt19937 random(11);

    // equal and opposite: the same link seen from its other end
    double antisymmetry = 0.0;
    double dissipation_asymmetry = 0.0;
    double smallest_dissipation = 1.0;
    for (int trial = 0; trial < 200; ++trial) {
        for (int k = 1; k < vb::kQ; ++k) {
            const vb::LinkEnd a = random_end(random, rho1, rho2);
            const vb::LinkEnd b = random_end(random, rho1, rho2);
            double jx = 0.0, jy = 0.0, kx = 0.0, ky = 0.0, da = 0.0, db = 0.0;
            vb::link_momentum(k, a, b, rho1, rho2, &jx, &jy, &da);
            vb::link_momentum(vb::kOpposite[k], b, a, rho1, rho2, &kx, &ky, &db);
            const double scale = std::max(std::hypot(jx, jy), 1e-300);
            antisymmetry = std::max(antisymmetry, std::hypot(jx + kx, jy + ky) / scale);
            dissipation_asymmetry =
                std::max(dissipation_asymmetry, std::fabs(da - db) / std::max(da, 1e-300));
            smallest_dissipation = std::min(smallest_dissipation, std::min(da, db));
        }
    }
    std::cout << "link_antisymmetry = " << antisymmetry << "\n";
    std::cout << "link_dissipation_asymmetry = " << dissipation_asymmetry << "\n";
    std::cout << "link_dissipation_min = " << smallest_dissipation << "\n";

    // over a periodic lattice, indexed as the droplet solver indexes it
    const int n = 12;
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<vb::LinkEnd> node(static_cast<std::size_t>(n) * n);
    std::vector<double> sent(node.size() * vb::kQ);
    std::vector<double> phase(node.size() * vb::kQ);
    for (std::size_t m = 0; m < node.size(); ++m) {
        node[m] = random_end(random, rho1, rho2);
        for (int k = 0; k < vb::kQ; ++k) {
            sent[m * vb::kQ + k] = 0.05 * (unit(random) - 0.5);
            phase[m * vb::kQ + k] = 0.1 * unit(random);
        }
    }
    double total_x = 0.0, total_y = 0.0, scale = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const std::size_t x = static_cast<std::size_t>(i) * n + j;
            for (int k = 1; k < vb::kQ; ++k) {
                const int id = (i - vb::kVelocity[k][0] + n) % n;
                const int jd = (j - vb::kVelocity[k][1] + n) % n;
                const std::size_t y = static_cast<std::size_t>(id) * n + jd;
                vb::LinkEnd donor = node[y];
                donor.outgoing = sent[y * vb::kQ + k];
                donor.phase = phase[y * vb::kQ + k];
                vb::LinkEnd receiver = node[x];
                receiver.outgoing = sent[x * vb::kQ + vb::kOpposite[k]];
                receiver.phase = phase[x * vb::kQ + vb::kOpposite[k]];
                double jx = 0.0, jy = 0.0, d = 0.0;
                vb::link_momentum(k, donor, receiver, rho1, rho2, &jx, &jy, &d);
                total_x += jx;
                total_y += jy;
                scale = std::max(scale, std::hypot(jx, jy));
            }
        }
    }
    std::cout << "lattice_momentum_total = " << std::hypot(total_x, total_y) / scale << "\n";

    // a uniform velocity is carried by the mass flux alone, whatever the
    // densities and phase populations on the two sides
    double galilean = 0.0;
    const double ux = 0.013;
    const double uy = -0.007;
    double gamma[vb::kQ];
    vb::velocity_equilibrium(ux, uy, gamma);
    for (int trial = 0; trial < 200; ++trial) {
        for (int k = 1; k < vb::kQ; ++k) {
            vb::LinkEnd donor = random_end(random, rho1, rho2);
            vb::LinkEnd receiver = random_end(random, rho1, rho2);
            donor.ux = receiver.ux = ux;
            donor.uy = receiver.uy = uy;
            donor.outgoing = gamma[k] - vb::kWeight[k];
            receiver.outgoing = gamma[vb::kOpposite[k]] - vb::kWeight[k];
            double jx = 0.0, jy = 0.0, d = 0.0;
            vb::link_momentum(k, donor, receiver, rho1, rho2, &jx, &jy, &d);
            const double volume = gamma[k] - gamma[vb::kOpposite[k]];
            const double mass = (rho1 - rho2) * (donor.phase - receiver.phase) + rho2 * volume;
            // relative to the terms that cancel: the mass flux and the lighter
            // side's share of it
            const double rho_link = std::min(donor.rho, receiver.rho);
            const double cancelling =
                (std::fabs(mass) + rho_link * std::fabs(volume)) * std::hypot(ux, uy);
            galilean = std::max(galilean, std::hypot(jx - mass * ux, jy - mass * uy) / cancelling);
        }
    }
    std::cout << "link_uniform_velocity_error = " << galilean << "\n";

    // one component on both sides (c = 1): no excess mass flux, no added
    // viscosity, and the lattice exchange with its advection as the linear
    // part of the mass flux times the mean velocity
    double single = 0.0;
    double single_dissipation = 0.0;
    for (int trial = 0; trial < 200; ++trial) {
        for (int k = 1; k < vb::kQ; ++k) {
            vb::LinkEnd donor = random_end(random, rho1, rho2);
            vb::LinkEnd receiver = random_end(random, rho1, rho2);
            donor.rho = receiver.rho = rho1;
            donor.mu = receiver.mu = 0.3;
            double gd[vb::kQ], gr[vb::kQ];
            vb::velocity_equilibrium(donor.ux, donor.uy, gd);
            vb::velocity_equilibrium(receiver.ux, receiver.uy, gr);
            donor.phase = gd[k];
            receiver.phase = gr[vb::kOpposite[k]];
            double jx = 0.0, jy = 0.0, d = 0.0;
            vb::link_momentum(k, donor, receiver, rho1, rho2, &jx, &jy, &d);
            single_dissipation = std::max(single_dissipation, std::fabs(d));
            auto advective = [&](double vx, double vy) {
                const double eu = vb::kVelocity[k][0] * vx + vb::kVelocity[k][1] * vy;
                return vb::kWeight[k] *
                       (0.5 * eu * eu / (kCs2 * kCs2) - 0.5 * (vx * vx + vy * vy) / kCs2);
            };
            const double lattice =
                rho1 * (donor.outgoing + receiver.outgoing - advective(donor.ux, donor.uy) -
                        advective(receiver.ux, receiver.uy));
            const double linear = vb::kWeight[k] *
                                  (vb::kVelocity[k][0] * (donor.ux + receiver.ux) +
                                   vb::kVelocity[k][1] * (donor.uy + receiver.uy)) /
                                  kCs2;
            const double ex =
                rho1 * linear * 0.5 * (donor.ux + receiver.ux) + lattice * vb::kVelocity[k][0];
            const double ey =
                rho1 * linear * 0.5 * (donor.uy + receiver.uy) + lattice * vb::kVelocity[k][1];
            single = std::max(single,
                              std::hypot(jx - ex, jy - ey) / std::max(std::hypot(ex, ey), 1e-300));
        }
    }
    std::cout << "link_single_component_error = " << single << "\n";
    std::cout << "link_single_component_dissipation = " << single_dissipation << "\n";

    // the dissipation force, with the same coefficient at the two ends of
    // every link, sums to zero, and vanishes for a uniform velocity
    std::vector<double> coefficient(node.size() * vb::kQ, 0.0);
    std::vector<double> vel_x(node.size());
    std::vector<double> vel_y(node.size());
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const std::size_t x = static_cast<std::size_t>(i) * n + j;
            vel_x[x] = 0.02 * (unit(random) - 0.5);
            vel_y[x] = 0.02 * (unit(random) - 0.5);
            for (int k = 1; k < vb::kQ; ++k) {
                if (k == 1 || k == 2 || k == 5 || k == 6) {
                    // one value per link, stored at both of its ends
                    const int id = (i - vb::kVelocity[k][0] + n) % n;
                    const int jd = (j - vb::kVelocity[k][1] + n) % n;
                    const std::size_t y = static_cast<std::size_t>(id) * n + jd;
                    const double value = 1.0e3 * unit(random);
                    coefficient[x * vb::kQ + k] = value;
                    coefficient[y * vb::kQ + vb::kOpposite[k]] = value;
                }
            }
        }
    }
    double force_x = 0.0, force_y = 0.0, force_scale = 0.0, uniform = 0.0;
    std::vector<double> still_x(node.size(), 0.013);
    std::vector<double> still_y(node.size(), -0.007);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double fx = 0.0, fy = 0.0;
            vb::dissipation_force(
                coefficient.data(), vel_x.data(), vel_y.data(), n, n, i, j, &fx, &fy);
            force_x += fx;
            force_y += fy;
            force_scale = std::max(force_scale, std::hypot(fx, fy));
            vb::dissipation_force(
                coefficient.data(), still_x.data(), still_y.data(), n, n, i, j, &fx, &fy);
            uniform = std::max(uniform, std::hypot(fx, fy));
        }
    }
    std::cout << "dissipation_force_total = " << std::hypot(force_x, force_y) / force_scale << "\n";
    std::cout << "dissipation_force_uniform = " << uniform / force_scale << "\n";
}

/// The hybrid collision is collide at sigma = 1, conserves P and u, and does
/// not depend on sigma when the populations' non-equilibrium is the one the
/// velocity gradient predicts.
void report_collide_hybrid() {
    const double P = 0.021;
    const double ux = 0.011;
    const double uy = -0.017;
    const double tau = 0.5003;
    const double tau_bulk = 1.0;
    double eq[vb::kQ];
    double source[vb::kQ];
    vb::hydrodynamic_equilibrium(P, ux, uy, eq);
    vb::forcing(ux, uy, 2.0e-5, -1.0e-5, source);
    vb::VelocityGradient gradient;
    gradient.dux_dx = 3.0e-4;
    gradient.dux_dy = -2.0e-4;
    gradient.duy_dx = 5.0e-4;
    gradient.duy_dy = -1.0e-4;

    double state[vb::kQ];
    for (int k = 0; k < vb::kQ; ++k) {
        state[k] = eq[k] + 1.0e-4 * std::sin(2.1 * k + 0.7);
    }
    double plain[vb::kQ];
    double hybrid[vb::kQ];
    vb::collide(state, eq, source, tau, tau_bulk, plain);
    vb::collide_hybrid(state, eq, source, tau, tau_bulk, 1.0, gradient, hybrid);
    double same = 0.0;
    for (int k = 0; k < vb::kQ; ++k) {
        same = std::max(same, std::fabs(plain[k] - hybrid[k]));
    }
    std::cout << "hybrid_unit_weight_error = " << same << "\n";

    // the equilibrium of the state's own moments, as in a collision
    const Moments before = moments(state);
    double own[vb::kQ];
    vb::hydrodynamic_equilibrium(before.m0, before.mx, before.my, own);
    vb::collide_hybrid(state, own, source, tau, tau_bulk, 0.7, gradient, hybrid);
    const Moments after = moments(hybrid);
    double cerr = std::fabs(after.m0 - before.m0);
    cerr = std::max(cerr, std::fabs(after.mx - (before.mx + 0.5 * 2.0e-5)));
    cerr = std::max(cerr, std::fabs(after.my - (before.my - 0.5 * 1.0e-5)));
    std::cout << "hybrid_conservation_error = " << cerr << "\n";

    // populations whose non-equilibrium is the Chapman-Enskog one
    const double divergence = gradient.dux_dx + gradient.duy_dy;
    const double fxx =
        -tau * kCs2 * (2.0 * gradient.dux_dx - divergence) - tau_bulk * kCs2 * divergence;
    const double fyy =
        -tau * kCs2 * (2.0 * gradient.duy_dy - divergence) - tau_bulk * kCs2 * divergence;
    const double fxy = -tau * kCs2 * (gradient.dux_dy + gradient.duy_dx);
    for (int k = 0; k < vb::kQ; ++k) {
        const double ex = vb::kVelocity[k][0];
        const double ey = vb::kVelocity[k][1];
        state[k] = eq[k] - 0.5 * source[k] +
                   vb::kWeight[k] / (2.0 * kCs2 * kCs2) *
                       ((ex * ex - kCs2) * fxx + (ey * ey - kCs2) * fyy + 2.0 * ex * ey * fxy);
    }
    double at_one[vb::kQ];
    vb::collide_hybrid(state, eq, source, tau, tau_bulk, 1.0, gradient, at_one);
    vb::collide_hybrid(state, eq, source, tau, tau_bulk, 0.3, gradient, hybrid);
    double independent = 0.0;
    for (int k = 0; k < vb::kQ; ++k) {
        independent = std::max(independent, std::fabs(at_one[k] - hybrid[k]));
    }
    std::cout << "hybrid_resolved_error = " << independent << "\n";
}

/// set_velocity replaces the first moment and nothing else.
void report_set_velocity() {
    double pop[vb::kQ];
    for (int k = 0; k < vb::kQ; ++k) {
        pop[k] = 0.1 + 0.01 * std::sin(1.7 * k + 0.3);
    }
    const Moments before = moments(pop);
    vb::set_velocity(pop, 0.017, -0.004);
    const Moments after = moments(pop);
    double err = std::fabs(after.m0 - before.m0);
    err = std::max(err, std::fabs(after.mx - 0.017));
    err = std::max(err, std::fabs(after.my + 0.004));
    err = std::max(err, std::fabs(after.mxx - before.mxx));
    err = std::max(err, std::fabs(after.myy - before.myy));
    err = std::max(err, std::fabs(after.mxy - before.mxy));
    std::cout << "set_velocity_error = " << err << "\n";
}

/// The phase transport alone, at rest, on a slab: profile, mass, bounds.
/// At |u| = 0.2 the carrier against the flow is smaller than the unlimited
/// sharpening; the limited populations stay within [0, Gamma_k].
void report_phase_limiter() {
    const double speed = 0.2;
    const double width = 1.6;
    double below = 0.0;
    double above = 0.0;
    double mass = 0.0;
    double unlimited_below = 0.0;
    for (double c : {0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999}) {
        for (int a = 0; a < 16; ++a) {
            for (int b = 0; b < 16; ++b) {
                const double ux = speed * std::cos(2.0 * M_PI * a / 16.0);
                const double uy = speed * std::sin(2.0 * M_PI * a / 16.0);
                const double nx = std::cos(2.0 * M_PI * (b + 0.5) / 16.0);
                const double ny = std::sin(2.0 * M_PI * (b + 0.5) / 16.0);
                double h[vb::kQ];
                double gamma[vb::kQ];
                vb::phase_populations(c, ux, uy, nx, ny, width, h);
                vb::velocity_equilibrium(ux, uy, gamma);
                double sum = 0.0;
                for (int k = 0; k < vb::kQ; ++k) {
                    below = std::min(below, h[k]);
                    above = std::min(above, gamma[k] - h[k]);
                    sum += h[k];
                    // what the unlimited sharpening would give
                    const double en = vb::kVelocity[k][0] * nx + vb::kVelocity[k][1] * ny;
                    const double sharpening = 0.5 * kCs2 * 2.0 * c * (1.0 - c) / width;
                    unlimited_below = std::min(
                        unlimited_below, c * gamma[k] + vb::kWeight[k] * sharpening * en / kCs2);
                }
                mass = std::max(mass, std::fabs(sum - c));
            }
        }
    }
    std::cout << "phase_limited_min = " << below << "\n";
    std::cout << "phase_limited_complement_min = " << above << "\n";
    std::cout << "phase_limited_mass_error = " << mass << "\n";
    std::cout << "phase_unlimited_min = " << unlimited_below << "\n";
}

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
    report_pressure_force();
    report_link_momentum();
    report_collide_hybrid();
    report_set_velocity();
    report_phase_limiter();
    report_phase_transport();
    std::cout << std::flush;
    return 0;
}
