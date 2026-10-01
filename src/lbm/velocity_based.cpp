#include "lbm/velocity_based.h"

#include "lbm/isotropic_gradient.h"
#include "lbm/surface_force.h"

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

namespace {

/// D2Q9 weight of direction k at lattice temperature `temperature`.
double phase_weight(int k, double temperature) {
    if (temperature == kCs2) {
        return kWeight[k];
    }
    if (k == 0) {
        return 1.0 - 5.0 * temperature / 3.0;
    }
    return k <= 4 ? temperature / 3.0 : temperature / 12.0;
}

/// What the second-order Hermite carrier at `temperature` misses of
/// temperature I + u u on D2Q9: the diagonal components on the axis pairs
/// against the rest population, the off-diagonal one on the diagonals, with
/// no mass and no momentum. Zero at cs^2.
double moment_correction(int k, double ux, double uy, double temperature) {
    const double ax =
        ux * ux * (1.5 - 0.5 / temperature) + uy * uy * (0.5 - 1.0 / (6.0 * temperature));
    const double ay =
        uy * uy * (1.5 - 0.5 / temperature) + ux * ux * (0.5 - 1.0 / (6.0 * temperature));
    const double b = ux * uy * (1.0 - 1.0 / (3.0 * temperature));
    const int ex = kVelocity[k][0];
    const int ey = kVelocity[k][1];
    if (ex == 0 && ey == 0) {
        return -ax - ay;
    }
    if (ey == 0) {
        return 0.5 * ax;
    }
    if (ex == 0) {
        return 0.5 * ay;
    }
    return 0.25 * b * ex * ey;
}

/// Gamma_k(u) at `temperature`, for one direction.
double carrier(int k, double ux, double uy, double temperature) {
    const double eu = kVelocity[k][0] * ux + kVelocity[k][1] * uy;
    if (temperature == kCs2) {
        return kWeight[k] *
               (1.0 + eu / kCs2 + 0.5 * eu * eu / kCs4 - 0.5 * (ux * ux + uy * uy) / kCs2);
    }
    return phase_weight(k, temperature) *
               (1.0 + eu / temperature + 0.5 * eu * eu / (temperature * temperature) -
                0.5 * (ux * ux + uy * uy) / temperature) +
           moment_correction(k, ux, uy, temperature);
}

}  // namespace

void phase_carrier(double ux, double uy, double temperature, double* gamma) {
    if (temperature == kCs2) {
        velocity_equilibrium(ux, uy, gamma);
        return;
    }
    for (int k = 0; k < kQ; ++k) {
        gamma[k] = carrier(k, ux, uy, temperature);
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
                       double* populations,
                       double temperature,
                       double flux_x,
                       double flux_y) {
    double gamma[kQ];
    phase_carrier(ux, uy, temperature, gamma);
    // the sharpening vanishes outside [0, 1], so a rounding overshoot of c is
    // not amplified
    const double bounded = (c < 0.0) ? 0.0 : ((c > 1.0) ? 1.0 : c);
    const double mobility = 0.5 * temperature;
    const double sharpening = mobility * 2.0 * bounded * (1.0 - bounded) / width;
    double flux[kQ];
    double theta = 1.0;
    for (int k = 0; k < kQ; ++k) {
        const double carried = kVelocity[k][0] * (sharpening * normal_x + flux_x) +
                               kVelocity[k][1] * (sharpening * normal_y + flux_y);
        flux[k] = phase_weight(k, temperature) * carried / temperature;
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

double
field_at(const double* field, int nx, int ny, int i, int j, Boundary boundary, double parity) {
    const int ip = (i % nx + nx) % nx;
    if (boundary != Boundary::WallY) {
        return field[ip * ny + (j % ny + ny) % ny];
    }
    if (j < 0) {
        return parity * field[ip * ny + (-1 - j)];
    }
    if (j >= ny) {
        return parity * field[ip * ny + (2 * ny - 1 - j)];
    }
    return field[ip * ny + j];
}

void field_gradient(const double* field,
                    int nx,
                    int ny,
                    int i,
                    int j,
                    GradientStencil stencil,
                    Boundary boundary,
                    double parity,
                    double* grad_x,
                    double* grad_y) {
    if (boundary == Boundary::WallY) {
        gradient_mirror_y(field, nx, ny, i, j, stencil, parity, grad_x, grad_y);
    } else {
        gradient_periodic(field, nx, ny, i, j, stencil, grad_x, grad_y);
    }
}

double lattice_laplacian(const double* field, int nx, int ny, int i, int j, Boundary boundary) {
    const double here = field[i * ny + j];
    double sum = 0.0;
    for (int k = 1; k < kQ; ++k) {
        sum += kWeight[k] *
               (field_at(field, nx, ny, i + kVelocity[k][0], j + kVelocity[k][1], boundary) - here);
    }
    return 2.0 * sum / kCs2;
}

void interface_normal(const double* psi,
                      const double* laplacian_psi,
                      int nx,
                      int ny,
                      int i,
                      int j,
                      double* normal_x,
                      double* normal_y,
                      Boundary boundary) {
    double gx = 0.0, gy = 0.0, hx = 0.0, hy = 0.0;
    field_gradient(psi, nx, ny, i, j, GradientStencil::E4, boundary, 1.0, &gx, &gy);
    field_gradient(laplacian_psi, nx, ny, i, j, GradientStencil::E4, boundary, 1.0, &hx, &hy);
    unit_normal(gx - hx / 6.0, gy - hy / 6.0, normal_x, normal_y);
}

void phase_correction_flux(const double* laplacian_c,
                           int nx,
                           int ny,
                           int i,
                           int j,
                           double temperature,
                           double* flux_x,
                           double* flux_y,
                           Boundary boundary) {
    double hx = 0.0, hy = 0.0;
    field_gradient(laplacian_c, nx, ny, i, j, GradientStencil::E4, boundary, 1.0, &hx, &hy);
    *flux_x = -temperature / 24.0 * hx;
    *flux_y = -temperature / 24.0 * hy;
}

void sixth_order_normal(const double* psi,
                        const double* laplacian_psi,
                        const double* laplacian2_psi,
                        int nx,
                        int ny,
                        int i,
                        int j,
                        double* normal_x,
                        double* normal_y,
                        Boundary boundary) {
    double gx = 0.0, gy = 0.0, hx = 0.0, hy = 0.0, kx = 0.0, ky = 0.0;
    field_gradient(psi, nx, ny, i, j, GradientStencil::E4, boundary, 1.0, &gx, &gy);
    field_gradient(laplacian_psi, nx, ny, i, j, GradientStencil::E4, boundary, 1.0, &hx, &hy);
    field_gradient(laplacian2_psi, nx, ny, i, j, GradientStencil::E4, boundary, 1.0, &kx, &ky);
    auto at = [&](int ii, int jj) { return field_at(psi, nx, ny, ii, jj, boundary); };
    // (f(3) - 4 f(2) + 5 f(1) - 5 f(-1) + 4 f(-2) - f(-3)) / 2 = f^(5) + O(h^2)
    auto fifth = [&](int di, int dj) {
        return 0.5 * (at(i + 3 * di, j + 3 * dj) - 4.0 * at(i + 2 * di, j + 2 * dj) +
                      5.0 * at(i + di, j + dj) - 5.0 * at(i - di, j - dj) +
                      4.0 * at(i - 2 * di, j - 2 * dj) - at(i - 3 * di, j - 3 * dj));
    };
    unit_normal(gx - hx / 6.0 + kx / 36.0 + fifth(1, 0) / 180.0,
                gy - hy / 6.0 + ky / 36.0 + fifth(0, 1) / 180.0,
                normal_x,
                normal_y);
}

void sixth_order_flux(const double* laplacian_c,
                      const double* laplacian2_c,
                      int nx,
                      int ny,
                      int i,
                      int j,
                      double temperature,
                      double* flux_x,
                      double* flux_y,
                      Boundary boundary) {
    double hx = 0.0, hy = 0.0, kx = 0.0, ky = 0.0;
    field_gradient(laplacian_c, nx, ny, i, j, GradientStencil::E4, boundary, 1.0, &hx, &hy);
    field_gradient(laplacian2_c, nx, ny, i, j, GradientStencil::E4, boundary, 1.0, &kx, &ky);
    // D_x d_yy L c and D_y d_xx L c, from L c directly
    auto l = [&](int di, int dj) {
        return field_at(laplacian_c, nx, ny, i + di, j + dj, boundary);
    };
    const double mixed_x =
        0.5 * ((l(1, 1) - 2.0 * l(1, 0) + l(1, -1)) - (l(-1, 1) - 2.0 * l(-1, 0) + l(-1, -1)));
    const double mixed_y =
        0.5 * ((l(1, 1) - 2.0 * l(0, 1) + l(-1, 1)) - (l(1, -1) - 2.0 * l(0, -1) + l(-1, -1)));
    *flux_x = temperature * (-hx / 24.0 + 7.0 * kx / 480.0 - mixed_x / 360.0);
    *flux_y = temperature * (-hy / 24.0 + 7.0 * ky / 480.0 - mixed_y / 360.0);
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

/// relax, with the deviator relaxed at `tau_along` where it stretches along
/// the interface of unit normal n and at `tau_across` where it shears across
/// it. In the frame of n and its tangent t the deviator has two components,
/// D_nn = a cos 2t + b sin 2t and D_nt = b cos 2t - a sin 2t, with
/// a = (pxx - pyy) / 2, b = pxy and t the angle of n; each relaxes at its own
/// rate and the two are rotated back. Without a normal, both take tau_along.
void relax_laminate(double pxx,
                    double pyy,
                    double pxy,
                    const double* equilibrium,
                    const double* source,
                    double tau_across,
                    double tau_along,
                    double tau_bulk,
                    double normal_x,
                    double normal_y,
                    double* post_collision) {
    const double cos2 = normal_x * normal_x - normal_y * normal_y;
    const double sin2 = 2.0 * normal_x * normal_y;
    const double a = 0.5 * (pxx - pyy);
    const double b = pxy;
    const bool interface = cos2 != 0.0 || sin2 != 0.0;
    const double along = (1.0 - 1.0 / tau_along) * (interface ? a * cos2 + b * sin2 : a);
    const double across =
        (1.0 - 1.0 / (interface ? tau_across : tau_along)) * (interface ? b * cos2 - a * sin2 : b);
    const double deviator = interface ? along * cos2 - across * sin2 : along;
    const double qxy = interface ? along * sin2 + across * cos2 : across;
    const double trace = (1.0 - 1.0 / tau_bulk) * 0.5 * (pxx + pyy);
    const double qxx = trace + deviator;
    const double qyy = trace - deviator;
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

void collide_filtered(const double* populations,
                      const double* equilibrium,
                      const double* source,
                      double tau_shear,
                      double tau_bulk,
                      double sigma,
                      double* previous,
                      double* post_collision) {
    double pxx, pyy, pxy;
    non_equilibrium(populations, equilibrium, source, &pxx, &pyy, &pxy);
    const double mean_xx = 0.5 * (pxx + previous[0]);
    const double mean_yy = 0.5 * (pyy + previous[1]);
    const double mean_xy = 0.5 * (pxy + previous[2]);
    previous[0] = pxx;
    previous[1] = pyy;
    previous[2] = pxy;
    relax(sigma * pxx + (1.0 - sigma) * mean_xx,
          sigma * pyy + (1.0 - sigma) * mean_yy,
          sigma * pxy + (1.0 - sigma) * mean_xy,
          equilibrium,
          source,
          tau_shear,
          tau_bulk,
          post_collision);
}

void collide_laminate(const double* populations,
                      const double* equilibrium,
                      const double* source,
                      double tau_across,
                      double tau_along,
                      double tau_bulk,
                      double sigma,
                      double normal_x,
                      double normal_y,
                      double* previous,
                      double* post_collision) {
    double pxx, pyy, pxy;
    non_equilibrium(populations, equilibrium, source, &pxx, &pyy, &pxy);
    const double mean_xx = 0.5 * (pxx + previous[0]);
    const double mean_yy = 0.5 * (pyy + previous[1]);
    const double mean_xy = 0.5 * (pxy + previous[2]);
    previous[0] = pxx;
    previous[1] = pyy;
    previous[2] = pxy;
    relax_laminate(sigma * pxx + (1.0 - sigma) * mean_xx,
                   sigma * pyy + (1.0 - sigma) * mean_yy,
                   sigma * pxy + (1.0 - sigma) * mean_xy,
                   equilibrium,
                   source,
                   tau_across,
                   tau_along,
                   tau_bulk,
                   normal_x,
                   normal_y,
                   post_collision);
}

void pressure_force(const double* pressure_number,
                    const double* rho,
                    int nx,
                    int ny,
                    int i,
                    int j,
                    double* ax,
                    double* ay,
                    Boundary boundary) {
    double sx = 0.0;
    double sy = 0.0;
    for (int k = 1; k < kQ; ++k) {
        const int im = ((i - kVelocity[k][0]) % nx + nx) % nx;
        const int jm = j - kVelocity[k][1];
        // p / cs^2 = rho P of the upstream node
        double upstream = 0.0;
        if (boundary == Boundary::WallY && (jm < 0 || jm >= ny)) {
            // extrapolated through the wall from the two rows beside it
            const int first = jm < 0 ? 0 : ny - 1;
            const int second = jm < 0 ? 1 : ny - 2;
            const int a = im * ny + first;
            const int b = im * ny + second;
            upstream = 2.0 * rho[a] * pressure_number[a] - rho[b] * pressure_number[b];
        } else {
            const int m = im * ny + (jm % ny + ny) % ny;
            upstream = rho[m] * pressure_number[m];
        }
        const double weight = kWeight[k] * upstream;
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
                   double* dissipation,
                   double temperature) {
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
    const double volume = carrier(k, donor.ux, donor.uy, temperature) -
                          carrier(opposite, receiver.ux, receiver.uy, temperature);
    const double linear = phase_weight(k, temperature) *
                          (ex * (donor.ux + receiver.ux) + ey * (donor.uy + receiver.uy)) /
                          temperature;
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

LinkEnd wall_image(const LinkEnd& receiver) {
    LinkEnd image = receiver;
    image.ux = -receiver.ux;
    image.uy = -receiver.uy;
    return image;
}

void dissipation_force(const double* coefficients,
                       const double* ux,
                       const double* uy,
                       int nx,
                       int ny,
                       int i,
                       int j,
                       double* fx,
                       double* fy,
                       Boundary boundary) {
    const int here = i * ny + j;
    double sx = 0.0;
    double sy = 0.0;
    for (int k = 1; k < kQ; ++k) {
        int jd = j - kVelocity[k][1];
        if (boundary == Boundary::WallY && (jd < 0 || jd >= ny)) {
            // the node's own image, with no dissipation on the link
            continue;
        }
        jd = (jd % ny + ny) % ny;
        const int id = ((i - kVelocity[k][0]) % nx + nx) % nx;
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
