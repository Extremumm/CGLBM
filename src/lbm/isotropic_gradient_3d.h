#ifndef CGLBM_LBM_ISOTROPIC_GRADIENT_3D_H
#define CGLBM_LBM_ISOTROPIC_GRADIENT_3D_H

/// Isotropic gradient stencils in three dimensions.
///
/// The counterpart of `isotropic_gradient.h`, and it matters for the same
/// reason: the colour gradient is differentiated twice per step, once for the
/// interface normal and once for its divergence, and anisotropy there becomes
/// an anisotropic surface tension and spurious currents at the interface. In
/// two dimensions, going from the fourth-order stencil to the eighth-order one
/// was worth an order of magnitude on those currents.
///
/// A gradient of the form
///
///     grad psi(x) = sum_l W(|c_l|^2) c_l psi(x + c_l)
///
/// is isotropic to order 2n when the lattice tensors `sum_l W c_a c_b ...` of
/// rank up to 2n match their continuum counterparts. Solving those conditions
/// for the shell weights gives, exactly:
///
///     E4:  W(1) = 1/6,   W(2) = 1/12
///          18 neighbours, `sum W c_x^2 = 1`, fourth-order isotropic.
///     E6:  W(1) = 2/15,  W(2) = 1/15,  W(3) = 1/60,  W(4) = 1/120
///          32 neighbours, additionally sixth-order isotropic.
///
/// `E4` is the D3Q19 neighbourhood, so it costs nothing beyond the streaming
/// pattern the solver already touches. `E6` adds the eight corner neighbours
/// and the six second-axial ones, reaching two nodes along each axis; a field
/// that is not periodic then needs two valid nodes outside the region where the
/// gradient is evaluated.
///
/// Reference
///  - M. Sbragaglia, R. Benzi, L. Biferale, S. Succi, K. Sugiyama, F. Toschi,
///    "Generalized lattice Boltzmann method with multirange pseudopotential",
///    Phys. Rev. E 75, 026702 (2007). The construction of the shell weights.

namespace cglbm {
namespace lbm {

/// Isotropy order of the three-dimensional gradient stencil.
enum class GradientStencil3D {
    E4,  ///< 18 neighbours, isotropic to fourth order
    E6   ///< 32 neighbours, isotropic to sixth order
};

/// One term of a stencil: a lattice offset and its weight.
struct StencilPoint3D {
    int cx;
    int cy;
    int cz;
    double weight;
};

/// The points of `stencil`, with `count` set to how many there are.
///
/// The returned array has static storage; it must not be freed.
const StencilPoint3D* stencil_points_3d(GradientStencil3D stencil, int* count);

/// Largest offset the stencil reaches along any axis: 1 for E4, 2 for E6.
int stencil_reach_3d(GradientStencil3D stencil);

/// Name of the stencil, as accepted by :func:`stencil_from_name_3d`.
const char* stencil_name_3d(GradientStencil3D stencil);

/// Parse "E4" or "E6", case-insensitively.
bool stencil_from_name_3d(const char* name, GradientStencil3D* stencil);

/// Gradient of a triply periodic field at node (i, j, k).
///
/// `field` is `nx * ny * nz` values indexed `field[(i * ny + j) * nz + k]`,
/// which is the memory order of a `depth == 1` `Field3D`.
void gradient_periodic_3d(const double* field,
                          int nx,
                          int ny,
                          int nz,
                          int i,
                          int j,
                          int k,
                          GradientStencil3D stencil,
                          double* grad_x,
                          double* grad_y,
                          double* grad_z);

/// Gradient of a field periodic in x and z and bounded by walls along y.
///
/// A neighbour that would leave the domain through a y wall contributes
/// nothing, which is the same rule the two-dimensional wall gradient applies.
void gradient_wall_y_3d(const double* field,
                        int nx,
                        int ny,
                        int nz,
                        int i,
                        int j,
                        int k,
                        GradientStencil3D stencil,
                        double* grad_x,
                        double* grad_y,
                        double* grad_z);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_ISOTROPIC_GRADIENT_3D_H
