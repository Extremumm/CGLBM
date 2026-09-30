#ifndef CGLBM_LBM_CENTRAL_MOMENTS_H
#define CGLBM_LBM_CENTRAL_MOMENTS_H

#include "lbm/d2q9.h"

/// Central moments on D2Q9, and the colour-gradient collision of Saito et al.
/// (2023) built on them.
///
/// The colour-gradient equilibrium of `Solver`, that of Ba et al. (2016) and
/// Wen et al. (2019),
///
///     f_eq = g_eq,2 + (p - rho cs^2) (E + Phi_3),
///
/// has the right moments up to the second order, and in the third the
/// correction Phi_3 = w (u_x H_xyy + u_y H_xxy) / (2 cs^6) removes the part of
/// the non-ideal pressure that moves with the fluid. Its higher central moments
/// still depend on the velocity. Saito et al. write the equilibrium as
///
///     f_eq = g_eq,N + (p - rho cs^2) (E + Phi),
///
/// with g_eq,N the Maxwellian projected onto every Hermite polynomial the
/// lattice carries, and choose Phi, order by order, so that no central moment
/// depends on the velocity. On D2Q9 (the tensor Hermite polynomials H_ab,
/// a, b <= 2) that is
///
///     g_eq,4 = rho w [1 + (u_x H_10 + u_y H_01) / cs^2
///                     + (u_x^2 H_20 + u_y^2 H_02 + 2 u_x u_y H_11) / (2 cs^4)
///                     + (u_x^2 u_y H_21 + u_x u_y^2 H_12) / (2 cs^6)
///                     + u_x^2 u_y^2 H_22 / (4 cs^8)],
///     E      = w [(H_20 + H_02) / (2 cs^4) - H_22 / (4 cs^6)],
///     Phi    = w [(u_x H_12 + u_y H_21) / (2 cs^6) + |u|^2 H_22 / (4 cs^8)],
///
/// and its central moments k_ab = sum_i f_i (xi_x - u_x)^a (xi_y - u_y)^b are
///
///     k_00 = rho,  k_20 = k_02 = p,  k_22 = p cs^2,  all others 0,
///
/// whatever the velocity (checked symbolically in the derivation and by the unit
/// test). With Phi = 0, k_12 = -(p - rho cs^2) u_x, k_21 = -(p - rho cs^2) u_y
/// and k_22 = p cs^2 + (p - rho cs^2) |u|^2. With Phi_3 alone the third-order
/// moments vanish but k_22 = p cs^2 - (p - rho cs^2) |u|^2, and the equilibrium
/// of `Solver`, on g_eq,2, adds k_21 = -rho u_x^2 u_y, k_12 = -rho u_x u_y^2 and
/// 3 rho u_x^2 u_y^2 to k_22. In a heavy fluid p - rho cs^2 is about -rho cs^2:
/// at a density ratio of 1000 and |u| = 0.01, k_22 is 30 % above p cs^2.
///
/// The collision is done in central moments. The shear moments k_20 - k_02 and
/// k_11 relax at the rate the viscosity sets, and the trace and the third- and
/// fourth-order moments go to their equilibria, a rate of 1, as Saito et al.
/// set every rate but the shear's. At second order that is the regularised
/// collision `Solver` otherwise uses with a bulk relaxation time of 1; what
/// differs is the third and fourth order, which here keep the products of the
/// velocity with the relaxed stress that a raw-moment regularisation drops.
///
/// The trace has to go with k_22. The raw-moment regularisation leaves k_22 a
/// third of the trace's non-equilibrium; set to equilibrium while the trace
/// relaxes slowly, it takes that part from the rest population, and a uniform
/// fluid is linearly unstable: noise of 1e-6 overflows within a few hundred
/// steps at a trace relaxation time of 5 or 100, ideal gas or not. Tied to the
/// trace instead, k_22 is unstable at tau = 0.51; with both at rate 1 every
/// relaxation time from 0.51 to 100 is stable (research copy).
///
/// Reference
///  - S. Saito, N. Takada, S. Baba, S. Someya, H. Ito, "Generalized equilibria
///    for color-gradient lattice Boltzmann model based on higher-order Hermite
///    polynomials: A simplified implementation with central moments", Phys.
///    Rev. E 108, 065305 (2023), arXiv:2309.07801. Eqs. (25), (28), (29), (56)
///    and (57)-(63), for D3Q27; the D2Q9 form above is its two-dimensional
///    counterpart.
///  - Z. X. Wen, Q. Li, Y. Yu, K. H. Luo, "Improved three-dimensional
///    color-gradient lattice Boltzmann model for immiscible two-phase flows",
///    Phys. Rev. E 100, 023301 (2019), and Y. Ba et al., Phys. Rev. E 94,
///    023310 (2016): the third-order correction Phi_3, from the third-order
///    Hermite expansion of Q. Li et al., Phys. Rev. E 85, 016710 (2012).

namespace cglbm {
namespace lbm {

/// k[a][b] = sum_i f_i (xi_ix - u_x)^a (xi_iy - u_y)^b for a, b in {0, 1, 2}, in
/// the velocity order of `kXi`.
void central_moments(const double* f, double ux, double uy, double k[3][3]);

/// The D2Q9 populations whose central moments about (ux, uy) are `k`: the
/// inverse of central_moments. The nine monomials xi_x^a xi_y^b with a, b <= 2
/// are a basis on D2Q9, so the populations are fixed by them.
void populations_from_central_moments(const double k[3][3], double ux, double uy, double* f);

/// g_eq,4 + (p - rho cs^2) (E + Phi): the generalized equilibrium, in
/// population space.
void generalized_equilibrium(double rho, double ux, double uy, double p, double* f);

/// Post-collision populations of Saito et al.'s central-moment collision.
///
/// `f` is what the collision relaxes (in `Solver`, the populations with half
/// the source term added) and (ux, uy) the velocity its first central moments
/// vanish about. The shear moments relax at `omega_shear` towards 0, the trace
/// at `omega_bulk` towards 2 p, k_21 and k_12 are set to 0 and k_22 to
/// p cs^2; k_00, k_10 and k_01 are kept.
void collide_central_moment(const double* f,
                            double ux,
                            double uy,
                            double p,
                            double omega_shear,
                            double omega_bulk,
                            double* post);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_CENTRAL_MOMENTS_H
