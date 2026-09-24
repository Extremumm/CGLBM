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
    velocity_equilibrium(ux, uy, populations);
    // the sharpening vanishes outside [0, 1], so a rounding overshoot of c is
    // not amplified
    const double bounded = (c < 0.0) ? 0.0 : ((c > 1.0) ? 1.0 : c);
    const double mobility = 0.5 * kCs2;
    const double sharpening = mobility * 2.0 * bounded * (1.0 - bounded) / width;
    for (int k = 0; k < kQ; ++k) {
        const double en = kVelocity[k][0] * normal_x + kVelocity[k][1] * normal_y;
        populations[k] = c * populations[k] + kWeight[k] * sharpening * en / kCs2;
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

void collide(const double* populations,
             const double* equilibrium,
             const double* source,
             double tau_shear,
             double tau_bulk,
             double* post_collision) {
    double pxx = 0.0;
    double pyy = 0.0;
    double pxy = 0.0;
    for (int k = 0; k < kQ; ++k) {
        const double ex = kVelocity[k][0];
        const double ey = kVelocity[k][1];
        const double neq = populations[k] - equilibrium[k] + 0.5 * source[k];
        pxx += neq * ex * ex;
        pyy += neq * ey * ey;
        pxy += neq * ex * ey;
    }
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

void pressure_correction(const double* pressure_number,
                         const double* rho,
                         int nx,
                         int ny,
                         int i,
                         int j,
                         double* ax,
                         double* ay) {
    const double rho_here = rho[i * ny + j];
    double sx = 0.0;
    double sy = 0.0;
    for (int k = 1; k < kQ; ++k) {
        const int ip = ((i + kVelocity[k][0]) % nx + nx) % nx;
        const int jp = ((j + kVelocity[k][1]) % ny + ny) % ny;
        const int m = ip * ny + jp;
        const double weight = kWeight[k] * pressure_number[m] * (1.0 - rho[m] / rho_here);
        sx += weight * kVelocity[k][0];
        sy += weight * kVelocity[k][1];
    }
    *ax = sx;
    *ay = sy;
}

void viscous_correction(double nu,
                        double dux_dx,
                        double dux_dy,
                        double duy_dx,
                        double duy_dy,
                        double dlnrho_dx,
                        double dlnrho_dy,
                        double* ax,
                        double* ay) {
    const double sxx = 2.0 * dux_dx;
    const double syy = 2.0 * duy_dy;
    const double sxy = dux_dy + duy_dx;
    *ax = nu * (sxx * dlnrho_dx + sxy * dlnrho_dy);
    *ay = nu * (sxy * dlnrho_dx + syy * dlnrho_dy);
}

}  // namespace velocity_based
}  // namespace lbm
}  // namespace cglbm
