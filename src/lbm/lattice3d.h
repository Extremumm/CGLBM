#ifndef CGLBM_LBM_LATTICE3D_H
#define CGLBM_LBM_LATTICE3D_H

#include "lbm/d3q19.h"
#include "lbm/d3q27.h"

/// One description of a three-dimensional lattice, so the solver has one body.
///
/// `TwoPopulationSolver3D` runs on D3Q19 or D3Q27. The two differ in eight
/// velocities, four numbers, and nothing else -- but three of those four are
/// derived quantities that a port carries over by hand, and getting one wrong
/// is invisible to every conservation check the solver makes. That happened
/// once already on the way from D2Q9 to D3Q19: the enhanced equilibrium's
/// amplitude was carried over unchanged, which left the heavy fluid with 417
/// times its intended shear viscosity while mass, momentum and the second
/// moment all stayed exactly right. See `two_population_solver_3d.h`.
///
/// So the rule here is that a lattice states only what is arbitrary about it --
/// its velocities, its weights, and the two coefficients that follow from
/// solving the rest-weight conditions -- and everything else is *computed* from
/// those by `TwoPopulationSolver3D`, identically for both. There is no place
/// left to write a number that is right for one lattice and wrong for the
/// other.
///
/// # The rest weights
///
/// The two-population model gives each fluid its own rest weight `alpha_k` and
/// carries the density ratio in it. Fixing the remaining weights `phi_q^k`
/// takes two conditions,
///
///     sum_q phi_q = 1                          (mass)
///     sum_q phi_q e_x^4 = 3 sum_q phi_q e_x^2 e_y^2   (fourth-order isotropy)
///
/// and, on D3Q27, one more, because there are three moving shells rather than
/// two. Taking the sixth-order relation `sum phi e_x^4 e_y^2 = 3 sum phi
/// e_x^2 e_y^2 e_z^2` -- the one part of sixth-order isotropy the lattice can
/// satisfy -- closes it, and reproduces the standard weights at
/// `alpha = 8/27`. The solutions are
///
///     D3Q19:  phi = alpha, (1-a)/12, (1-a)/24
///             (c_s^k)^2 = (1 - alpha) / 2
///     D3Q27:  phi = alpha, 2(1-a)/19, (1-a)/38, (1-a)/152
///             (c_s^k)^2 = 9 (1 - alpha) / 19
///
/// with `(1-a)` short for `1 - alpha_k`. Both leave the relation that matters,
///
///     rho_1 / rho_2 = (1 - alpha_2) / (1 - alpha_1),
///
/// unchanged -- the coefficient in `(c_s^k)^2` cancels between the two fluids --
/// and with it the identity of the two bulk pressures `rho_k (c_s^k)^2`.
///
/// # What the solver derives, and does not store
///
/// *The enhanced equilibrium's amplitude.* The correction
/// `3 (e.u) lambda (3|e|^2 - 5)` exists to drag the equilibrium's third moment
/// from the lattice's `c_s^2 = 1/3` to the fluid's own `(c_s^k)^2`, which is
/// what sets the shear viscosity. Taking `M3_xxy` of the equilibrium gives
///
///     M3_xxy / (rho_k u_y) = 3 S + 3 lambda T,
///     S = sum_q w_q e_x^2 e_y^2,   T = sum_q w_q e_x^2 e_y^2 (3|e_q|^2 - 5),
///
/// so `lambda = ((c_s^k)^2 - 3 S) / (3 T)`. `S` is `1/9` on both lattices; `T`
/// is `1/9` on D3Q19 and `2/9` on D3Q27, because the corner shell carries
/// `3|e|^2 - 5 = 4` where the edge shell carries 1. The amplitude is therefore
/// `3 (c_s^k)^2 - 1` on D3Q19 and *half* that on D3Q27 -- which is, exactly,
/// the wrong D3Q19 value that the port from two dimensions left behind. Two
/// lattices whose correct amplitudes differ by the same factor as an old bug is
/// as good a reason as exists not to write either of them down: the solver sums
/// `S` and `T` over the lattice it was handed.
///
/// *The velocity at which the equilibrium goes negative.* At small `u` the
/// equilibrium of direction `q` is `rho_k [phi_q + 3 w_q (e_q.u) A_q]` with
/// `A_q = 1 + lambda (3|e_q|^2 - 5)`, so it stays positive while
///
///     |u| < min_q  phi_q / (3 w_q |e_q| |A_q|).
///
/// In the high-ratio limit, where `(c_s^k)^2` is small, that is `(c_s^k)^2 / 3`
/// on D3Q19, from the axial shell, and on D3Q27 the axial, edge and corner
/// shells give `(c_s^k)^2 / 2`, `sqrt(2) (c_s^k)^2` and `(c_s^k)^2 / sqrt(3)`,
/// the first of which binds. Half again as much room, and it is the bound the
/// heavy fluid runs into first at a density ratio of 1000.
/// `programs/unit_testing/lbm/lattice_3d` measures all of it.

namespace cglbm {
namespace lbm {

/// Which three-dimensional lattice a case runs on.
enum class Lattice3DKind {
    D3Q19,  ///< 19 velocities, fourth-order isotropic
    D3Q27   ///< 27 velocities, the corners included
};

/// The velocities, weights and rest-weight coefficients of one lattice.
///
/// Everything here is a property of the velocity set. Nothing derived from a
/// fluid belongs in it.
struct Lattice3D {
    /// Name, as accepted by :func:`lattice_3d_from_name`.
    const char* name;

    /// Number of discrete velocities.
    int q;

    /// Discrete velocities, `xi[q][axis]`.
    const double (*xi)[3];

    /// Lattice weights of the standard equilibrium.
    const double* w;

    /// Index of the reversed velocity, for bounce-back.
    const int* opposite;

    /// `(c_s^k)^2 = cs2_coefficient * (1 - alpha_k)`.
    double cs2_coefficient;

    /// `phi_q^k = shell_coefficient[|e_q|^2] * (1 - alpha_k)` for `q != 0`.
    ///
    /// Indexed by the squared speed, which is 0, 1, 2 or 3. Entry 0 is unused
    /// -- the rest weight is `alpha_k` itself -- and entry 3 is zero on D3Q19,
    /// which has no corner shell.
    double shell_coefficient[4];
};

/// Squared speed of direction `q`: 0, 1, 2 or 3.
inline double speed_squared_3d(const Lattice3D& lattice, int q) {
    const double* e = lattice.xi[q];
    return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
}

/// The rest weight `phi_q^k` of a fluid with rest weight `alpha`.
inline double rest_weight_3d(const Lattice3D& lattice, int q, double alpha) {
    if (q == 0) {
        return alpha;
    }
    const int shell = static_cast<int>(speed_squared_3d(lattice, q) + 0.5);
    return lattice.shell_coefficient[shell] * (1.0 - alpha);
}

/// The sound speed `(c_s^k)^2` of a fluid with rest weight `alpha`.
inline double sound_speed_squared_3d(const Lattice3D& lattice, double alpha) {
    return lattice.cs2_coefficient * (1.0 - alpha);
}

inline constexpr Lattice3D kLatticeD3Q19 = {
    "D3Q19",
    kQ3D,
    kXi3D,
    kW3D,
    kOpposite3D,
    1. / 2.,
    {0.0, 1. / 12., 1. / 24., 0.0},
};

inline constexpr Lattice3D kLatticeD3Q27 = {
    "D3Q27",
    kQ3D27,
    kXi3D27,
    kW3D27,
    kOpposite3D27,
    9. / 19.,
    {0.0, 2. / 19., 1. / 38., 1. / 152.},
};

/// The lattice `kind` names.
inline const Lattice3D& lattice_3d(Lattice3DKind kind) {
    return kind == Lattice3DKind::D3Q27 ? kLatticeD3Q27 : kLatticeD3Q19;
}

/// Parse "D3Q19" or "D3Q27", case-insensitively.
bool lattice_3d_from_name(const char* name, Lattice3DKind* kind);

/// Name of `kind`, as accepted by :func:`lattice_3d_from_name`.
const char* lattice_3d_name(Lattice3DKind kind);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_LATTICE3D_H
