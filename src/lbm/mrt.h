#ifndef CGLBM_LBM_MRT_H
#define CGLBM_LBM_MRT_H

#include "lbm/d2q9.h"

/// Multiple-relaxation-time collision on D2Q9, in the moment basis of
/// Lallemand & Luo.
///
/// The basis is, one row per moment, in the velocity order of `kXi`:
///
///     rho, e, epsilon, j_x, q_x, j_y, q_y, p_xx, p_xy
///
/// Row `e` is `3|e_i|^2 - 4`, three times the trace of the stress less four
/// times the mass; `p_xx` is `e_x^2 - e_y^2` and `p_xy` is `e_x e_y`, the
/// deviatoric stress; `epsilon` and `q` are the fourth- and third-order moments
/// the Navier-Stokes equations do not contain. The rows are orthogonal, so the
/// inverse is the transpose divided by their squared norms.
///
/// References
///  - P. Lallemand, L.-S. Luo, "Theory of the lattice Boltzmann method:
///    Dispersion, dissipation, isotropy, Galilean invariance, and stability",
///    Phys. Rev. E 61, 6546 (2000), Eq. (7).
///  - Y. Ba, H. Liu, Q. Li, Q. Kang, J. Sun, Phys. Rev. E 94, 023310 (2016),
///    Eqs. (11)-(18) and Appendix A: the same basis for a colour-gradient
///    model, and the source term for the diagonal third moment.

namespace cglbm {
namespace lbm {

/// Lallemand & Luo's moment matrix `M`, row `a`, column `i`.
inline constexpr double kMoment[kQ][kQ] = {{1, 1, 1, 1, 1, 1, 1, 1, 1},       // rho
                                           {-4, -1, -1, -1, -1, 2, 2, 2, 2},  // e
                                           {4, -2, -2, -2, -2, 1, 1, 1, 1},   // epsilon
                                           {0, 1, 0, -1, 0, 1, -1, -1, 1},    // j_x
                                           {0, -2, 0, 2, 0, 1, -1, -1, 1},    // q_x
                                           {0, 0, 1, 0, -1, 1, 1, -1, -1},    // j_y
                                           {0, 0, -2, 0, 2, 1, 1, -1, -1},    // q_y
                                           {0, 1, -1, 1, -1, 0, 0, 0, 0},     // p_xx
                                           {0, 0, 0, 0, 0, 1, -1, 1, -1}};    // p_xy

/// Squared norm of each row of `kMoment`: `M^-1 = M^T diag(1 / kMomentNorm)`.
inline constexpr double kMomentNorm[kQ] = {9., 36., 36., 6., 12., 6., 12., 4., 4.};

/// Rows of `kMoment`, by name.
enum MomentRow {
    kMomentRho = 0,
    kMomentEnergy = 1,
    kMomentEnergySquare = 2,
    kMomentJx = 3,
    kMomentQx = 4,
    kMomentJy = 5,
    kMomentQy = 6,
    kMomentPxx = 7,
    kMomentPxy = 8
};

/// Relaxation rates in the order of `kMoment`'s rows.
///
/// The shear stress relaxes at `s_nu`, the rate the viscosity sets. The two
/// momenta take it too, and their rate is immaterial: with Guo's half-force
/// velocity `j - j^eq = -F dt / 2`, and any rate `s` gives
/// `s F / 2 + (1 - s / 2) F = F`. The mass is conserved and its rate is zero.
void mrt_rates(double s_nu, double s_e, double s_eps, double s_q, double* rates);

/// One collision in moment space, into `out[kQ]`:
///
///     out = f - M^-1 S M (f - f_eq) + M^-1 (I - S / 2) M source dt
///
/// with `S = diag(rates)` and `source` Guo's forcing term in velocity space.
/// Every rate equal to `1 / tau` is BGK with Guo's forcing, exactly.
void mrt_collide(const double* f,
                 const double* f_eq,
                 const double* source,
                 const double* rates,
                 double dt,
                 double* out);

/// The source term that restores the diagonal third moment, in velocity space.
///
/// D2Q9 gives `sum_i f_i^eq e_x^3 = rho_k u_x` where a fluid with its own sound
/// speed needs `3 p_k u_x`; the difference, summed over the fluids, is
/// `Q = sum_k (1 - 3 (c_s^k)^2) rho_k u`. Its divergence goes back into the
/// trace of the stress (row `e`) and into the normal-stress difference (row
/// `p_xx`), each with the `1 - s / 2` of a source term -- Ba et al. Eqs. (17),
/// (18):
///
///     C_e    = 3 (1 - s_e / 2) (d_x Q_x + d_y Q_y) dt
///     C_p_xx =   (1 - s_nu / 2) (d_x Q_x - d_y Q_y) dt
///
/// The two brackets are passed separately, `divergence` and
/// `normal_difference`, because they need not come from the same stencil:
/// see `SourceStencil`. It carries no mass and no momentum. Added to `out[kQ]`.
void add_third_moment_source(
    double divergence, double normal_difference, double s_e, double s_nu, double dt, double* out);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_MRT_H
