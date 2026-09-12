#ifndef CGLBM_LBM_D3Q27_H
#define CGLBM_LBM_D3Q27_H

/// The D3Q27 lattice: D3Q19 plus the eight corner neighbours.
///
/// The full three-dimensional stencil -- the rest particle, the six axial
/// neighbours, the twelve edge neighbours and the eight corners. It costs 42 %
/// more memory and streaming than D3Q19, and buys two things the colour-gradient
/// model at a high density ratio actually needs.
///
///  - *Isotropy.* D3Q19 leaves the fourth-order moment of the lattice weights
///    isotropic but not the sixth, and the colour gradient is differentiated
///    twice per step. Saito et al. measured the spurious velocity around a
///    static droplet at 1.2e-2 on D3Q19 against 5.8e-3 on D3Q27, everything
///    else held equal.
///
///  - *Headroom in the equilibrium.* The two-population model carries the
///    density ratio in the rest weight, which leaves the heavy fluid with a
///    very small `(c_s^k)^2` and correspondingly small moving weights. The
///    equilibrium goes negative once the velocity passes a bound set by the
///    tightest of them; see `lattice3d.h` for the arithmetic, which gives
///    `(c_s^k)^2 / 3` on D3Q19 and `(c_s^k)^2 / 2` here, half again as much.
///
/// The first nineteen velocities are those of `d3q19.h`, in the same order, so
/// that the two lattices agree wherever they overlap and the D3Q19 tables are a
/// prefix of these. Nothing depends on that, but a reader comparing the two
/// should not have to.
///
/// Reference
///  - S. Saito, Y. Abe, K. Koyama, "Lattice Boltzmann modeling and simulation
///    of liquid jet breakup", Phys. Rev. E 96, 013317 (2017),
///    doi:10.1103/PhysRevE.96.013317, arXiv:1705.03141. Eq. (1) is this
///    velocity set and Eq. (12) these weights; Sec. III A reports the spurious
///    velocity against D3Q19.
///  - Y. Qian, D. d'Humieres, P. Lallemand, "Lattice BGK models for
///    Navier-Stokes equation", Europhysics Letters 17(6), 479 (1992). The DnQm
///    family and its weights.

namespace cglbm {
namespace lbm {

/// Number of discrete velocities.
inline constexpr int kQ3D27 = 27;

/// Discrete velocities: rest, then axial, then edge, then corner.
inline constexpr double kXi3D27[kQ3D27][3] = {
    {0, 0, 0},  //
    {1, 0, 0},
    {-1, 0, 0},
    {0, 1, 0},
    {0, -1, 0},
    {0, 0, 1},   //
    {0, 0, -1},  //
    {1, 1, 0},
    {-1, -1, 0},
    {1, -1, 0},
    {-1, 1, 0},  //
    {1, 0, 1},
    {-1, 0, -1},
    {1, 0, -1},
    {-1, 0, 1},  //
    {0, 1, 1},
    {0, -1, -1},
    {0, 1, -1},
    {0, -1, 1},  //
    {1, 1, 1},
    {-1, -1, -1},
    {1, 1, -1},
    {-1, -1, 1},
    {1, -1, 1},
    {-1, 1, -1},
    {-1, 1, 1},
    {1, -1, -1}  //
};

/// Index of the reversed velocity, `kXi3D27[kOpposite3D27[i]] == -kXi3D27[i]`.
inline constexpr int kOpposite3D27[kQ3D27] = {
    0,  2,  1,  4,  3,  6,  5,  8,  7,  10, 9,  12, 11, 14,
    13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};

/// Lattice weights: 8/27 at rest, 2/27 axial, 1/54 edge, 1/216 corner.
inline constexpr double kW3D27[kQ3D27] = {
    8. / 27.,  //
    2. / 27.,
    2. / 27.,
    2. / 27.,
    2. / 27.,
    2. / 27.,
    2. / 27.,  //
    1. / 54.,
    1. / 54.,
    1. / 54.,
    1. / 54.,
    1. / 54.,
    1. / 54.,  //
    1. / 54.,
    1. / 54.,
    1. / 54.,
    1. / 54.,
    1. / 54.,
    1. / 54.,  //
    1. / 216.,
    1. / 216.,
    1. / 216.,
    1. / 216.,  //
    1. / 216.,
    1. / 216.,
    1. / 216.,
    1. / 216.  //
};

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_D3Q27_H
