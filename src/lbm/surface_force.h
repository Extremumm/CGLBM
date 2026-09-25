#ifndef CGLBM_LBM_SURFACE_FORCE_H
#define CGLBM_LBM_SURFACE_FORCE_H

#include "lbm/isotropic_gradient.h"

/// Surface tension as the divergence of a capillary stress.
///
///     T = sigma/2 * (|g| I - g g / |g|),   g = grad(psi),   F = div(T)
///
/// with psi the normalised colour field of mixture.h, +1 in component 1. The
/// stress is a tension sigma along the interface, spread over it with the
/// weight |grad(psi)|/2, which integrates to one since psi changes by 2. In the
/// continuum div(T) = sigma/2 * kappa * g with kappa = -div(g/|g|): the force
/// points into component 1 across a droplet of it, and a static droplet of
/// radius R carries the Laplace jump sigma / R. It enters the scheme through
/// the ordinary forcing term, like gravity.
///
/// Why a force rather than the perturbation operator Omega^(2) the solvers
/// used before: Omega^(2) writes this same stress into the non-equilibrium
/// populations scaled by 1/tau, and streaming carries it from node to node
/// before it is relaxed. The two cancel exactly only where tau is uniform.
/// Across an interface with a viscosity contrast they do not: on the Laplace
/// case at density ratio 1, with the droplet 100 times more viscous than its
/// surroundings, the jump comes out at 0.86 sigma/R, and at density ratio 20
/// with a single kinematic viscosity -- a factor 20 in tau across the
/// interface -- at 0.71 sigma/R once the stress sits on the density interface.
/// The force does not involve tau: 1.03 and 1.02 on the same two cases.
///
/// Why the divergence form rather than the curvature form sigma/2 kappa g of
/// Brackbill et al. and Lishchuk et al.: the two agree in the continuum, but
/// only the divergence of a stress conserves momentum on the lattice. With a
/// centred stencil the discrete divergence of any field sums to zero over a
/// periodic lattice, so the interface exerts no net force on the fluid,
/// whatever shape it takes. The curvature form does not have that property,
/// and on a droplet translating at 0.01 lattice units per step at density
/// ratio 100 it drained 36 % of the total momentum in 3000 steps.
///
/// References
///  - B. Lafaurie, C. Nardone, R. Scardovelli, S. Zaleski, G. Zanetti,
///    "Modelling merging and fragmentation in multiphase flows with SURFER",
///    J. Comput. Phys. 113, 134-147 (1994). The capillary stress tensor.
///  - J. U. Brackbill, D. B. Kothe, C. Zemach, "A continuum method for modeling
///    surface tension", J. Comput. Phys. 100, 335-354 (1992). The same force
///    in curvature form.
///  - C. Kublik, R. Tsai, "Integration over curves and surfaces defined by the
///    closest point mapping", Res. Math. Sci. 3, 3 (2016), arXiv:1504.05478,
///    proposition 2; C. Kublik, N. M. Tanushev, R. Tsai, J. Comput. Phys. 247,
///    279-311 (2013). The Jacobian 1 + d kappa between a curve and the curve at
///    distance d from it, which layer_weight uses.

namespace cglbm {
namespace lbm {

/// Below this norm the colour gradient is taken as zero: there is no interface.
constexpr double kInterfaceGradientThreshold = 1.0e-10;

/// Unit normal grad(psi) / |grad(psi)|, or (0, 0) away from any interface.
void unit_normal(double grad_x, double grad_y, double* normal_x, double* normal_y);

/// Capillary stress sigma/2 (|g| I - g g / |g|) at a node, from its colour
/// gradient g; zero away from any interface.
void capillary_stress(double sigma,
                      double grad_x,
                      double grad_y,
                      double* stress_xx,
                      double* stress_xy,
                      double* stress_yy);

/// Weight 1/J that carries the tension of each layer of a diffuse interface
/// onto the interface psi = 0.
///
/// The capillary stress spreads the tension over the interface width, and a
/// layer at distance d from psi = 0 is bent to its own radius: across a
/// droplet of radius R the pressure jump comes out as sigma times the mean of
/// 1/r over the layers rather than sigma / R, 3 % high at R = 10 with W = 1.6.
/// In 2D the length of psi = 0 is J = 1 + d div(n) times the length of the
/// layer at distance d, with n the unit normal into component 1 (Kublik & Tsai
/// 2016). A stress divided by J pulls in every layer as the interface itself
/// does: the jump across a circle is sigma / R whatever its width, a flat
/// interface is unchanged, and the force is still the divergence of a stress,
/// which conserves momentum.
///
/// d = -W atanh(psi) is the distance from psi = 0 on the equilibrium profile
/// psi = -tanh(d / W), positive in component 2; J is kept within [1/4, 4], a
/// bound that only a layer many widths away from a small droplet reaches.
double layer_weight(double psi, double divergence_of_normal, double width);

/// Surface force div(T) at node (i, j).
///
/// The three components of T are stored for every node, `nx * ny` values each
/// indexed `[i * ny + j]` as in isotropic_gradient.h, and are differentiated
/// with the same stencil as the colour field.
void surface_force(const double* stress_xx,
                   const double* stress_xy,
                   const double* stress_yy,
                   int nx,
                   int ny,
                   int i,
                   int j,
                   GradientStencil stencil,
                   Boundary boundary,
                   double* force_x,
                   double* force_y);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_SURFACE_FORCE_H
