#include "lbm/central_moments.h"

namespace cglbm {
namespace lbm {

namespace {

constexpr double kCs2 = 1.0 / 3.0;

/// One axis of the inverse: the populations at xi = -1, 0, 1 whose moments of
/// order 0, 1, 2 are M0, M1, M2 are (M2 - M1) / 2, M0 - M2 and (M2 + M1) / 2.
/// Row xi + 1, column the order.
constexpr double kAxisInverse[3][3] = {{0.0, -0.5, 0.5}, {1.0, 0.0, -1.0}, {0.0, 0.5, 0.5}};

constexpr double kBinomial[3][3] = {{1, 0, 0}, {1, 1, 0}, {1, 2, 1}};

}  // namespace

void central_moments(const double* f, double ux, double uy, double k[3][3]) {
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            k[a][b] = 0.0;
        }
    }
    for (int i = 0; i < kQ; ++i) {
        const double cx = kXi[i][0] - ux;
        const double cy = kXi[i][1] - uy;
        const double px[3] = {1.0, cx, cx * cx};
        const double py[3] = {1.0, cy, cy * cy};
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                k[a][b] += f[i] * px[a] * py[b];
            }
        }
    }
}

void populations_from_central_moments(const double k[3][3], double ux, double uy, double* f) {
    // raw moments m_ab = sum_i f_i xi_x^a xi_y^b, from xi = (xi - u) + u
    const double powx[3] = {1.0, ux, ux * ux};
    const double powy[3] = {1.0, uy, uy * uy};
    double m[3][3];
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            double sum = 0.0;
            for (int r = 0; r <= a; ++r) {
                for (int s = 0; s <= b; ++s) {
                    sum += kBinomial[a][r] * kBinomial[b][s] * powx[a - r] * powy[b - s] * k[r][s];
                }
            }
            m[a][b] = sum;
        }
    }
    for (int i = 0; i < kQ; ++i) {
        const double* lx = kAxisInverse[static_cast<int>(kXi[i][0]) + 1];
        const double* ly = kAxisInverse[static_cast<int>(kXi[i][1]) + 1];
        double sum = 0.0;
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                sum += lx[a] * ly[b] * m[a][b];
            }
        }
        f[i] = sum;
    }
}

void generalized_equilibrium(double rho, double ux, double uy, double p, double* f) {
    const double k[3][3] = {{rho, 0.0, p}, {0.0, 0.0, 0.0}, {p, 0.0, p * kCs2}};
    populations_from_central_moments(k, ux, uy, f);
}

void collide_central_moment(const double* f,
                            double ux,
                            double uy,
                            double p,
                            double omega_shear,
                            double omega_bulk,
                            double* post) {
    double k[3][3];
    central_moments(f, ux, uy, k);
    const double difference = (1.0 - omega_shear) * (k[2][0] - k[0][2]);
    const double trace = (1.0 - omega_bulk) * (k[2][0] + k[0][2]) + omega_bulk * 2.0 * p;
    k[2][0] = 0.5 * (trace + difference);
    k[0][2] = 0.5 * (trace - difference);
    k[1][1] *= 1.0 - omega_shear;
    k[2][1] = 0.0;
    k[1][2] = 0.0;
    k[2][2] = p * kCs2;
    populations_from_central_moments(k, ux, uy, post);
}

}  // namespace lbm
}  // namespace cglbm
