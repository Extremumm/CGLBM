#ifndef CGLBM_LBM_ISOTROPIC_GRADIENT_H
#define CGLBM_LBM_ISOTROPIC_GRADIENT_H

/// Discrete gradient stencils of increasing isotropy order.
///
/// The colour-gradient model differentiates the phase field twice per step: in
/// the surface-tension operator, where the result is divided by its own norm,
/// and in the recolouring operator, which follows its direction. Anisotropy in
/// that gradient therefore turns directly into an anisotropic surface tension
/// and into spurious currents at the interface.
///
/// The gradient of a lattice field is approximated as
///
///     grad phi(x) = sum_l W(|c_l|^2) c_l phi(x + c_l)
///
/// where the shell weights `W` are chosen so that the lattice tensors
/// `sum_l W c_a c_b`, `sum_l W c_a c_b c_g c_d`, ... match their isotropic
/// continuum counterparts up to a given order. `E4` reproduces the eight
/// nearest neighbours of D2Q9 and is what the solvers used originally; `E6` and
/// `E8` reach further and are isotropic to higher order.
///
/// References
///  - M. Sbragaglia, R. Benzi, L. Biferale, S. Succi, K. Sugiyama, F. Toschi,
///    "Generalized lattice Boltzmann method with multirange pseudopotential",
///    Phys. Rev. E 75, 026702 (2007). Construction of the isotropic shell
///    weights; gives W(1) = 1/3, W(2) = 1/12 for the standard D2Q9 case and
///    W(1) = 4/15, W(2) = 1/10, W(4) = 1/120 for the sixth-order set.
///  - S. Leclaire, M. Reggio, J.-Y. Trepanier, "Isotropic color gradient for
///    simulating very high-density ratios with a two-phase flow lattice
///    Boltzmann model", Computers & Fluids 48(1), 98-112 (2011),
///    doi:10.1016/j.compfluid.2011.04.001. Shows that replacing the
///    nearest-neighbour colour gradient by a higher-order isotropic one lets
///    Laplace's law hold up to density ratios of O(10^4) and cuts the spurious
///    currents by about an order of magnitude.
///  - Eighth-order weights W(1) = 4/21, W(2) = 4/45, W(4) = 1/60,
///    W(5) = 2/315, W(8) = 1/5040, as tabulated in arXiv:2505.23647,
///    "Higher-order tuning of interface physics in multiphase lattice
///    Boltzmann", appendix B.

namespace cglbm {
namespace lbm {

/// Isotropy order of the gradient stencil.
enum class GradientStencil {
    E4,  ///< 8 neighbours, isotropic to 4th order. The original D2Q9 stencil.
    E6,  ///< 12 neighbours, isotropic to 6th order.
    E8   ///< 24 neighbours, isotropic to 8th order.
};

/// One term of a stencil: a lattice offset and its weight.
struct StencilPoint {
    int cx;
    int cy;
    double weight;
};

/// The points of `stencil`, with `count` set to how many there are.
///
/// The returned array has static storage; it must not be freed.
const StencilPoint* stencil_points(GradientStencil stencil, int* count);

/// Largest offset the stencil reaches along either axis: 1 for E4, 2 otherwise.
///
/// A field whose boundary is not periodic needs at least this many valid nodes
/// outside the region where the gradient is evaluated.
int stencil_reach(GradientStencil stencil);

/// Name of the stencil, as accepted by :func:`stencil_from_name`.
const char* stencil_name(GradientStencil stencil);

/// Parse "E4", "E6" or "E8", case-insensitively.
///
/// Returns false and leaves `stencil` untouched when `name` matches none.
bool stencil_from_name(const char* name, GradientStencil* stencil);

/// Gradient of a doubly periodic field at node (i, j).
///
/// `field` is `nx * ny` values indexed `field[i * ny + j]`, i.e. the memory
/// order of the solvers' `[Lx][Ly]` arrays. Both axes wrap.
void gradient_periodic(const double* field,
                       int nx,
                       int ny,
                       int i,
                       int j,
                       GradientStencil stencil,
                       double* grad_x,
                       double* grad_y);

/// Gradient of a field periodic along x and bounded by resting walls along y.
///
/// A neighbour that would leave the domain through a y wall contributes
/// nothing, which is the rule the wall-bounded solvers apply inline: within
/// `stencil_reach` nodes of a wall the gradient becomes one-sided. For E4 this
/// reproduces the original eight-neighbour treatment exactly, term for term.
void gradient_wall_y(const double* field,
                     int nx,
                     int ny,
                     int i,
                     int j,
                     GradientStencil stencil,
                     double* grad_x,
                     double* grad_y);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_ISOTROPIC_GRADIENT_H
