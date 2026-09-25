#ifndef CGLBM_LBM_D3Q19_H
#define CGLBM_LBM_D3Q19_H

#include <cmath>

/// The D3Q19 lattice: the three-dimensional counterpart of `d2q9.h`.
///
/// Nineteen velocities -- the rest particle, the six axial neighbours and the
/// twelve edge neighbours -- with the corner neighbours of D3Q27 left out. It
/// is the standard choice for colour-gradient models in three dimensions:
/// isotropic enough for the hydrodynamics, and a third cheaper than D3Q27 in
/// both memory and streaming.
///
/// The velocities are ordered so that `kOpposite3D[i]` is the reversed
/// direction of `i`, which is what the bounce-back boundary needs. In two
/// dimensions the solvers exploited a numeric pattern for that (`k +- 2`); here
/// the table is written out, because the pattern would not survive anyone
/// reordering the velocities.
///
/// Reference
///  - Y. Qian, D. d'Humieres, P. Lallemand, "Lattice BGK models for
///    Navier-Stokes equation", Europhysics Letters 17(6), 479 (1992). The
///    DnQm family and its weights.

namespace cglbm {
namespace lbm {

/// Number of discrete velocities.
inline constexpr int kQ3D = 19;

/// Discrete velocities, rest first, then the six axial, then the twelve edges.
inline constexpr double kXi3D[kQ3D][3] = {
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
    {0, -1, 1}  //
};

/// Index of the reversed velocity, `kXi3D[kOpposite3D[i]] == -kXi3D[i]`.
inline constexpr int kOpposite3D[kQ3D] = {
    0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};

/// Lattice weights: 1/3 at rest, 1/18 on the axes, 1/36 on the edges.
inline constexpr double kW3D[kQ3D] = {
    1. / 3.,  //
    1. / 18.,
    1. / 18.,
    1. / 18.,
    1. / 18.,
    1. / 18.,
    1. / 18.,  //
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,  //
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.  //
};

/// Squared speed of direction `i`: 0, 1 or 2.
inline constexpr double speed_squared_3d(int i) {
    return kXi3D[i][0] * kXi3D[i][0] + kXi3D[i][1] * kXi3D[i][1] + kXi3D[i][2] * kXi3D[i][2];
}

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_D3Q19_H
