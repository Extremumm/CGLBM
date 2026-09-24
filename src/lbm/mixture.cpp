#include "lbm/mixture.h"

#include <algorithm>

namespace cglbm {
namespace lbm {

namespace {

/// Mass fraction of component 1, with phi clamped to its physical range.
double mass_fraction_1(double phi) {
    return 0.5 * (1.0 + std::min(1.0, std::max(-1.0, phi)));
}

/// Specific volume 1/rho_k(p) of a component, kept finite and positive.
///
/// p + p_k_inf <= 0 would mean the component is stretched past the point where
/// its branch gives a density at all; it is then treated as occupying an
/// arbitrarily large volume, which is the limit the branch approaches.
double specific_volume(double p, double p_inf, double c_squared) {
    const double floor = 1.0e-300;
    return c_squared / std::max(p + p_inf, floor);
}

}  // namespace

double component1_density(double p, const ComponentPair& components) {
    return (p + components.p1_inf) / components.c1_squared;
}

double component2_density(double p, const ComponentPair& components) {
    return (p + components.p2_inf) / components.c2_squared;
}

double volume_fraction(double rho, double phi, double p, const ComponentPair& components) {
    return rho * 0.5 * (1.0 + phi) / component1_density(p, components);
}

double normalised_phase(double phi, double p, const ComponentPair& components) {
    const double y1 = mass_fraction_1(phi);
    const double y2 = 1.0 - y1;
    // alpha_k = rho Y_k v_k(p); the common factor rho cancels from the ratio.
    const double a1 = y1 * specific_volume(p, components.p1_inf, components.c1_squared);
    const double a2 = y2 * specific_volume(p, components.p2_inf, components.c2_squared);
    return (a1 - a2) / (a1 + a2);
}

MixtureState
mixture_from_volume_fraction(double alpha1, double p, const ComponentPair& components) {
    const double alpha2 = 1.0 - alpha1;
    const double partial1 = alpha1 * component1_density(p, components);
    const double partial2 = alpha2 * component2_density(p, components);
    const double rho = partial1 + partial2;
    return {rho, (partial1 - partial2) / rho};
}

double mixture_kinematic_viscosity(double phi, double nu1, double nu2) {
    const double y1 = mass_fraction_1(phi);
    return y1 * nu1 + (1.0 - y1) * nu2;
}

}  // namespace lbm
}  // namespace cglbm
