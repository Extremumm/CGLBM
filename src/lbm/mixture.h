#ifndef CGLBM_LBM_MIXTURE_H
#define CGLBM_LBM_MIXTURE_H

#include "lbm/equation_of_state.h"

/// Mixture quantities of the two-component equation of state: volume fractions,
/// the normalised colour field, a consistent initial state and the mixture
/// viscosity.
///
/// The phase field of the colour-gradient method is a *mass* fraction. With
/// rho_1 and rho_2 the partial densities carried by f and g,
///
///     phi = (rho_1 - rho_2) / (rho_1 + rho_2),   Y_1 = (1 + phi) / 2,
///
/// and the equation of state of equation_of_state.h is the pressure-equilibrium
/// mixture of the two ideal-gas branches,
///
///     1 / rho = Y_1 / rho_1(p) + Y_2 / rho_2(p),   rho_k(p) = (p + p_k_inf) / c_k^2,
///
/// so that component k fills the volume fraction alpha_k = rho Y_k / rho_k(p),
/// and alpha_1 + alpha_2 = 1 wherever the pressure is the one the equation of
/// state returns.
///
/// Mass and volume fractions differ by a logistic shift,
///
///     logit(Y_1) = logit(alpha_1) + ln(rho_1(p) / rho_2(p)),
///
/// so a tanh profile of width W in one is a tanh profile of the same width in
/// the other, displaced by (W/2) ln(rho_1/rho_2). The density, linear in
/// alpha_1 at fixed pressure, follows the volume fraction. With W = 1.6 the
/// phi = 0 contour therefore sits 2.4 lattice units outside the density
/// interface at a density ratio of 20, and 7.4 units outside at 10^4.
/// Anything meant to act *on the interface* -- the surface tension, the
/// interface normal -- has to be built on the volume fraction rather than on
/// phi: that is the normalised colour field of :func:`normalised_phase`.
///
/// References
///  - S. Leclaire, M. Reggio, J.-Y. Trepanier, Computers & Fluids 48(1), 98-112
///    (2011), and Y. Ba, H. Liu, Q. Li, Q. Kang, J. Sun, Phys. Rev. E 94,
///    023310 (2016): the colour gradient taken on the normalised density
///    (rho_1/rho_1^0 - rho_2/rho_2^0) / (rho_1/rho_1^0 + rho_2/rho_2^0), which
///    is what lets their models reach large density ratios. With an equation
///    of state the reference densities become rho_k(p), and the normalised
///    density becomes the volume-fraction difference alpha_1 - alpha_2.
///  - T. Lafarge, P. Boivin, N. Odier, B. Cuenot, Phys. Fluids 33, 082110
///    (2021): the equation of state these fractions are defined by.

namespace cglbm {
namespace lbm {

/// Density of pure component 1 at pressure `p`: the inverse of its branch.
double component1_density(double p, const ComponentPair& components);

/// Density of pure component 2 at pressure `p`: the inverse of its branch.
double component2_density(double p, const ComponentPair& components);

/// Volume fraction alpha_1 of component 1 at a node, rho Y_1 / rho_1(p).
///
/// Not clamped: away from the pressure the equation of state returns for
/// (rho, phi), alpha_1 + alpha_2 differs from one, and that is worth seeing.
double volume_fraction(double rho, double phi, double p, const ComponentPair& components);

/// Normalised colour field psi = (alpha_1 - alpha_2) / (alpha_1 + alpha_2).
///
/// +1 in component 1, -1 in component 2, and centred on the density interface
/// rather than on phi = 0. The density cancels from the ratio, so psi depends
/// on phi and p only and stays within [-1, 1] even where (rho, phi, p) are not
/// exactly consistent; phi is clamped to [-1, 1] first, as the recolouring can
/// overshoot it by a rounding error.
double normalised_phase(double phi, double p, const ComponentPair& components);

/// A node of the mixture, as the solvers store it.
struct MixtureState {
    double rho;  ///< mixture density
    double phi;  ///< phase field, the mass-fraction difference Y_1 - Y_2
};

/// The mixture holding volume fraction `alpha1` of component 1 at pressure `p`.
///
/// Both components sit at the same pressure, each at its own density
/// rho_k(p), so that :func:`pressure` returns `p` for the result: an initial
/// state that does not start the run with a pressure transient. A density
/// linear in phi (`InitialState::LinearDensity`) puts the interface nodes far
/// from that equilibrium -- at a density ratio of 1000 their pressure comes out
/// about 300 times the ambient one. `Solver` itself starts from
/// `InitialState::MechanicalEquilibrium`, which solves the same problem by
/// bisection on the density; this function is what the unit tests check the
/// equation of state against.
MixtureState mixture_from_volume_fraction(double alpha1, double p, const ComponentPair& components);

/// Kinematic viscosity of the mixture, nu = Y_1 nu_1 + Y_2 nu_2.
///
/// Weighting the kinematic viscosities by mass fraction is the same as
/// weighting the dynamic viscosities by volume fraction, since
/// rho Y_k = alpha_k rho_k:
///
///     rho nu = alpha_1 mu_1 + alpha_2 mu_2.
///
/// With nu_1 == nu_2 it returns that common value, so a case written with a
/// single kinematic viscosity is reproduced exactly. `phi` is clamped to
/// [-1, 1], which keeps the result between nu_1 and nu_2.
double mixture_kinematic_viscosity(double phi, double nu1, double nu2);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_MIXTURE_H
