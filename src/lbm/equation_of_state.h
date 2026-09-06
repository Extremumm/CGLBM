#ifndef CGLBM_LBM_EQUATION_OF_STATE_H
#define CGLBM_LBM_EQUATION_OF_STATE_H

/// The two-component equation of state of the improved colour-gradient method.
///
/// Classical colour-gradient models tie the density ratio to the ratio of the
/// components' sound speeds, because the density contrast is carried by the
/// rest-particle weight of the equilibrium. That coupling is what limits them
/// to modest density ratios. Lafarge et al. remove it by giving each component
/// its own ideal-gas branch, with its own sound speed `c_k` and its own
/// pressure at infinity `p_k_inf`, and mixing the two branches through the
/// phase field. Density ratio and sound-speed ratio then become independent.
///
/// The mixture pressure is the positive root of the quadratic that matches both
/// branches at phi = +1 and phi = -1:
///
///     p = 1/2 ( rho c_hat^2 - p1_inf - p2_inf
///               + sqrt( (p2_inf - p1_inf + rho c_bar^2)^2
///                       + rho^2 (1 - phi^2) c1^2 c2^2 ) )
///
/// with the phase-weighted combinations
///
///     c_hat^2 = (c1^2 + c2^2)/2 + phi (c1^2 - c2^2)/2
///     c_bar^2 = (c1^2 - c2^2)/2 + phi (c1^2 + c2^2)/2
///
/// At phi = +1 this collapses to `rho c1^2 - p1_inf` and at phi = -1 to
/// `rho c2^2 - p2_inf`, so each bulk recovers its own ideal gas. Between them
/// the square root interpolates, and that interpolation is what carries the
/// density contrast across the interface.
///
/// Reference
///  - T. Lafarge, P. Boivin, N. Odier, B. Cuenot, "Improved color-gradient
///    method for lattice Boltzmann modeling of two-phase flows", Physics of
///    Fluids 33(8), 082110 (2021), doi:10.1063/5.0061638. The generalisation of
///    the colour-gradient method to an arbitrary equation of state, removing
///    the non-physical link between the density and sound-speed ratios.

namespace cglbm {
namespace lbm {

/// Sound speeds and reference pressures of the two components.
struct ComponentPair {
    double c1_squared;  ///< squared sound speed of component 1 (phi = +1)
    double c2_squared;  ///< squared sound speed of component 2 (phi = -1)
    double p1_inf;      ///< pressure at infinity of component 1
    double p2_inf;      ///< pressure at infinity of component 2
};

/// Mixture pressure at a node, from its density and phase field.
///
/// `phi` is expected in [-1, 1]; values slightly outside, as the recolouring
/// step can produce, are handled without a domain error in the square root.
double pressure(double rho, double phi, const ComponentPair& components);

/// The single-component pressure the model reduces to when the two branches
/// are collapsed onto one sound speed.
///
/// This is *not* the model's equation of state -- it is the linear mixing rule
/// that a debugging substitution once left in `calPhaseField`, kept here so the
/// old behaviour can still be reproduced deliberately and compared against.
/// It agrees with :func:`pressure` in the bulk whenever c1 == c2, and differs
/// across the interface, which is exactly where the pressure jump is set. On
/// the shipped Laplace case it costs 17 % of the jump at a density ratio of 10
/// and 28 % at 20.
double pressure_linear_mixing(double rho, double phi, double cs_squared,
                              const ComponentPair& components);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_EQUATION_OF_STATE_H
