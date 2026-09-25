#include "lbm/velocity_based.h"

namespace cglbm {
namespace lbm {
namespace velocity_based {

namespace {

constexpr double kCs2 = kSoundSpeedSquared;
constexpr double kCs4 = kSoundSpeedSquared * kSoundSpeedSquared;

}  // namespace

void velocity_equilibrium(double ux, double uy, double* gamma) {
    const double usq = ux * ux + uy * uy;
    for (int k = 0; k < kQ; ++k) {
        const double eu = kVelocity[k][0] * ux + kVelocity[k][1] * uy;
        gamma[k] = kWeight[k] * (1.0 + eu / kCs2 + 0.5 * eu * eu / kCs4 - 0.5 * usq / kCs2);
    }
}

void hydrodynamic_equilibrium(double pressure_number, double ux, double uy, double* equilibrium) {
    velocity_equilibrium(ux, uy, equilibrium);
    for (int k = 0; k < kQ; ++k) {
        equilibrium[k] += kWeight[k] * (pressure_number - 1.0);
    }
}

void phase_populations(double c,
                       double ux,
                       double uy,
                       double normal_x,
                       double normal_y,
                       double width,
                       double* populations) {
    double gamma[kQ];
    velocity_equilibrium(ux, uy, gamma);
    // the sharpening vanishes outside [0, 1], so a rounding overshoot of c is
    // not amplified
    const double bounded = (c < 0.0) ? 0.0 : ((c > 1.0) ? 1.0 : c);
    const double mobility = 0.5 * kCs2;
    const double sharpening = mobility * 2.0 * bounded * (1.0 - bounded) / width;
    double flux[kQ];
    double theta = 1.0;
    for (int k = 0; k < kQ; ++k) {
        const double en = kVelocity[k][0] * normal_x + kVelocity[k][1] * normal_y;
        flux[k] = kWeight[k] * sharpening * en / kCs2;
        // keep 0 <= c Gamma_k + theta flux_k <= Gamma_k
        if (flux[k] < 0.0 && bounded * gamma[k] < -theta * flux[k]) {
            theta = bounded * gamma[k] / -flux[k];
        } else if (flux[k] > 0.0 && (1.0 - bounded) * gamma[k] < theta * flux[k]) {
            theta = (1.0 - bounded) * gamma[k] / flux[k];
        }
    }
    if (theta < 0.0) {
        theta = 0.0;
    }
    for (int k = 0; k < kQ; ++k) {
        populations[k] = c * gamma[k] + theta * flux[k];
    }
}

void forcing(double ux, double uy, double ax, double ay, double* source) {
    for (int k = 0; k < kQ; ++k) {
        const double ex = kVelocity[k][0];
        const double ey = kVelocity[k][1];
        const double first = (ex * ax + ey * ay) / kCs2;
        const double second = ((ex * ex - kCs2) * 2.0 * ux * ax + (ey * ey - kCs2) * 2.0 * uy * ay +
                               2.0 * ex * ey * (ux * ay + uy * ax)) /
                              (2.0 * kCs4);
        source[k] = kWeight[k] * (first + second);
    }
}

namespace {

/// Relax the second-order non-equilibrium (pxx, pyy, pxy) and rebuild the
/// post-collision populations around the equilibrium.
void relax(double pxx,
           double pyy,
           double pxy,
           const double* equilibrium,
           const double* source,
           double tau_shear,
           double tau_bulk,
           double* post_collision) {
    const double trace = (1.0 - 1.0 / tau_bulk) * 0.5 * (pxx + pyy);
    const double deviator = (1.0 - 1.0 / tau_shear) * 0.5 * (pxx - pyy);
    const double qxx = trace + deviator;
    const double qyy = trace - deviator;
    const double qxy = (1.0 - 1.0 / tau_shear) * pxy;
    for (int k = 0; k < kQ; ++k) {
        const double ex = kVelocity[k][0];
        const double ey = kVelocity[k][1];
        const double neq = kWeight[k] / (2.0 * kCs4) *
                           ((ex * ex - kCs2) * qxx + (ey * ey - kCs2) * qyy + 2.0 * ex * ey * qxy);
        post_collision[k] = equilibrium[k] + neq + 0.5 * source[k];
    }
}

/// Second-order moments of the non-equilibrium g - g^eq + S/2.
void non_equilibrium(const double* populations,
                     const double* equilibrium,
                     const double* source,
                     double* pxx,
                     double* pyy,
                     double* pxy) {
    *pxx = 0.0;
    *pyy = 0.0;
    *pxy = 0.0;
    for (int k = 0; k < kQ; ++k) {
        const double ex = kVelocity[k][0];
        const double ey = kVelocity[k][1];
        const double neq = populations[k] - equilibrium[k] + 0.5 * source[k];
        *pxx += neq * ex * ex;
        *pyy += neq * ey * ey;
        *pxy += neq * ex * ey;
    }
}

/// The u u part of Gamma_k(u): w_k ((xi_k . u)^2 / (2 cs^4) - u^2 / (2 cs^2)).
double advective_part(int k, double ux, double uy) {
    const double eu = kVelocity[k][0] * ux + kVelocity[k][1] * uy;
    return kWeight[k] * (0.5 * eu * eu / kCs4 - 0.5 * (ux * ux + uy * uy) / kCs2);
}

/// Gamma_k(u) for one direction.
double carrier(int k, double ux, double uy) {
    const double eu = kVelocity[k][0] * ux + kVelocity[k][1] * uy;
    return kWeight[k] * (1.0 + eu / kCs2 + 0.5 * eu * eu / kCs4 - 0.5 * (ux * ux + uy * uy) / kCs2);
}

}  // namespace

void collide(const double* populations,
             const double* equilibrium,
             const double* source,
             double tau_shear,
             double tau_bulk,
             double* post_collision) {
    double pxx, pyy, pxy;
    non_equilibrium(populations, equilibrium, source, &pxx, &pyy, &pxy);
    relax(pxx, pyy, pxy, equilibrium, source, tau_shear, tau_bulk, post_collision);
}

void collide_hybrid(const double* populations,
                    const double* equilibrium,
                    const double* source,
                    double tau_shear,
                    double tau_bulk,
                    double sigma,
                    const VelocityGradient& gradient,
                    double* post_collision) {
    double pxx, pyy, pxy;
    non_equilibrium(populations, equilibrium, source, &pxx, &pyy, &pxy);
    const double divergence = gradient.dux_dx + gradient.duy_dy;
    const double fxx =
        -tau_shear * kCs2 * (2.0 * gradient.dux_dx - divergence) - tau_bulk * kCs2 * divergence;
    const double fyy =
        -tau_shear * kCs2 * (2.0 * gradient.duy_dy - divergence) - tau_bulk * kCs2 * divergence;
    const double fxy = -tau_shear * kCs2 * (gradient.dux_dy + gradient.duy_dx);
    relax(sigma * pxx + (1.0 - sigma) * fxx,
          sigma * pyy + (1.0 - sigma) * fyy,
          sigma * pxy + (1.0 - sigma) * fxy,
          equilibrium,
          source,
          tau_shear,
          tau_bulk,
          post_collision);
}

void pressure_force(const double* pressure_number,
                    const double* rho,
                    int nx,
                    int ny,
                    int i,
                    int j,
                    double* ax,
                    double* ay) {
    double sx = 0.0;
    double sy = 0.0;
    for (int k = 1; k < kQ; ++k) {
        const int im = ((i - kVelocity[k][0]) % nx + nx) % nx;
        const int jm = ((j - kVelocity[k][1]) % ny + ny) % ny;
        const int m = im * ny + jm;
        // w_i p / cs^2 of the upstream node, p = rho cs^2 P
        const double weight = kWeight[k] * rho[m] * pressure_number[m];
        sx += weight * kVelocity[k][0];
        sy += weight * kVelocity[k][1];
    }
    const double rho_here = rho[i * ny + j];
    *ax = sx / rho_here;
    *ay = sy / rho_here;
}

void link_momentum(int k,
                   const LinkEnd& donor,
                   const LinkEnd& receiver,
                   double rho1,
                   double rho2,
                   double* jx,
                   double* jy,
                   double* dissipation) {
    const int opposite = kOpposite[k];
    const double ex = kVelocity[k][0];
    const double ey = kVelocity[k][1];
    // the lighter end sets how much momentum the lattice's exchange may carry
    const double rho_link = donor.rho < receiver.rho ? donor.rho : receiver.rho;

    // lattice exchange per unit mass, without its advective part
    const double exchange = donor.outgoing + receiver.outgoing;
    const double advective =
        advective_part(k, donor.ux, donor.uy) + advective_part(k, receiver.ux, receiver.uy);
    const double lattice = rho_link * (exchange - advective);

    // advection: the mass the phase populations carry, less the quadratic part
    // of the carriers' volume at the lighter density, times the mean velocity
    const double volume =
        carrier(k, donor.ux, donor.uy) - carrier(opposite, receiver.ux, receiver.uy);
    const double linear =
        kWeight[k] * (ex * (donor.ux + receiver.ux) + ey * (donor.uy + receiver.uy)) / kCs2;
    const double mass = (rho1 - rho2) * (donor.phase - receiver.phase) + rho2 * volume;
    const double excess = mass - rho_link * volume;
    const double carried = rho_link * linear + excess;
    const double mid_x = 0.5 * (donor.ux + receiver.ux);
    const double mid_y = 0.5 * (donor.uy + receiver.uy);

    // viscous stress: raise the link viscosity to the harmonic mean of mu
    const double target = 2.0 * donor.mu * receiver.mu / (donor.mu + receiver.mu);
    const double lattice_viscosity =
        rho_link * 0.5 * (donor.mu / donor.rho + receiver.mu / receiver.rho);
    const double beta = target > lattice_viscosity ? target - lattice_viscosity : 0.0;

    *jx = lattice * ex + carried * mid_x;
    *jy = lattice * ey + carried * mid_y;
    *dissipation = 0.5 * (excess < 0.0 ? -excess : excess) + beta * 2.0 * kWeight[k] / kCs2;
}

void dissipation_force(const double* coefficients,
                       const double* ux,
                       const double* uy,
                       int nx,
                       int ny,
                       int i,
                       int j,
                       double* fx,
                       double* fy) {
    const int here = i * ny + j;
    double sx = 0.0;
    double sy = 0.0;
    for (int k = 1; k < kQ; ++k) {
        const int id = ((i - kVelocity[k][0]) % nx + nx) % nx;
        const int jd = ((j - kVelocity[k][1]) % ny + ny) % ny;
        const int there = id * ny + jd;
        const double coefficient = coefficients[here * kQ + k];
        sx += coefficient * (ux[there] - ux[here]);
        sy += coefficient * (uy[there] - uy[here]);
    }
    *fx = sx;
    *fy = sy;
}

void set_velocity(double* populations, double ux, double uy) {
    double mx = 0.0;
    double my = 0.0;
    for (int k = 0; k < kQ; ++k) {
        mx += populations[k] * kVelocity[k][0];
        my += populations[k] * kVelocity[k][1];
    }
    const double dx = ux - mx;
    const double dy = uy - my;
    for (int k = 1; k < kQ; ++k) {
        populations[k] += kWeight[k] * (kVelocity[k][0] * dx + kVelocity[k][1] * dy) / kCs2;
    }
}

}  // namespace velocity_based
}  // namespace lbm
}  // namespace cglbm
