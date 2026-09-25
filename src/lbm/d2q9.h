#ifndef CGLBM_LBM_D2Q9_H
#define CGLBM_LBM_D2Q9_H

#include <cmath>

/// The D2Q9 velocity set and the lattice constants derived from it.
///
/// These values were repeated verbatim at the top of every solver. They are
/// collected here so that a scheme change touches one definition rather than
/// five, and so that a program cannot silently disagree with the library about
/// what `cs2` means.
///
/// Reference
///  - T. Kruger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E. M. Viggen,
///    "The Lattice Boltzmann Method: Principles and Practice", Springer (2017),
///    chapter 3. The D2Q9 weights and the speed of sound c_s = dx / (sqrt(3) dt).

namespace cglbm {
namespace lbm {

/// Number of discrete velocities.
inline constexpr int kQ = 9;

/// Discrete velocities, in the order the solvers assume.
///
/// The bounce-back reflection of the wall-bounded cases relies on this order:
/// k and k + 2 are opposite for k = 2 and k = 5, 6, which is what lets the
/// streaming step reflect with `kp = k +/- 2`.
///
///     6   2   5
///      \  |  /
///     3 - 0 - 1
///      /  |  \.
///     7   4   8
inline constexpr double kXi[kQ][2] = {
    {0, 0},  //
    {1, 0},
    {0, 1},
    {-1, 0},
    {0, -1},  //
    {1, 1},
    {-1, 1},
    {-1, -1},
    {1, -1}  //
};

/// Lattice weights: 4/9 at rest, 1/9 on the axes, 1/36 on the diagonals.
inline constexpr double kW[kQ] = {
    4. / 9.,  //
    1. / 9.,
    1. / 9.,
    1. / 9.,
    1. / 9.,  //
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.  //
};

/// Below this, the norm of a colour gradient is treated as zero.
///
/// The surface-tension and recolouring operators both divide by that norm,
/// which vanishes away from the interface -- everywhere except a few nodes.
inline constexpr double kGradientEpsilon = 1.0e-10;

/// Below this, a node is treated as carrying no interface at all.
///
/// `kGradientEpsilon` only keeps a division safe. The curvature is a harder
/// case: it is built from the *unit* normal, so it does not shrink with the
/// gradient it came from, and in the bulk -- where the gradient is nothing but
/// round-off -- it is a large random number. The capillary force
/// `sigma K |grad phi_N| / 2` then fails to vanish there.
///
/// This bites in proportion to the density ratio, because phi_N amplifies noise
/// in phi by `d phi_N / d phi = rho1 / rho2` at `phi = +1`. Measured in the
/// heavy bulk of the Laplace droplet, `|grad phi_N|` starts at 7e-12 and grows;
/// at a density ratio of 10^4 it crosses 1e-10 by step 9000 and the run is
/// destroyed by step 19000. The physical gradient at the interface is 0.66, so
/// a floor of 1e-6 sits five orders below anything real and three above the
/// noise at a ratio of 10^5.
///
/// Ba et al. state the same rule in words: the interfacial tension "is applied
/// only at the lattice sites where two fluids coexist" (after their Eq. 24).
inline constexpr double kInterfaceGradientFloor = 1.0e-6;

/// Lattice spacing, time step, and the sound speed they imply.
///
/// `dx` and `dt` are 1 in every case shipped here; they are kept explicit
/// because the equations below carry them, and a lattice unit is only a unit
/// while that stays visible.
struct LatticeUnits {
    double dx = 1.0;
    double dt = 1.0;

    /// Speed of sound, c_s = dx / (sqrt(3) dt).
    double cs() const {
        return dx / std::sqrt(3.0) / dt;
    }
    double cs2() const {
        const double c = cs();
        return c * c;
    }
    double cs4() const {
        const double c2 = cs2();
        return c2 * c2;
    }
    double cs6() const {
        return cs4() * cs2();
    }
};

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_D2Q9_H
