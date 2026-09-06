#include "lbm/equation_of_state.h"

#include <cmath>

namespace cglbm {
namespace lbm {

double pressure(double rho, double phi, const ComponentPair& components) {
    const double c1_squared = components.c1_squared;
    const double c2_squared = components.c2_squared;
    const double p1_inf = components.p1_inf;
    const double p2_inf = components.p2_inf;

    // Phase-weighted sound speeds: c_hat picks out the arithmetic mean and
    // c_bar the difference, so that phi = +/-1 selects one branch exactly.
    const double c_hat_squared =
        0.5 * (c1_squared + c2_squared) + 0.5 * phi * (c1_squared - c2_squared);
    const double c_bar_squared =
        0.5 * (c1_squared - c2_squared) + 0.5 * phi * (c1_squared + c2_squared);

    const double linear = p2_inf - p1_inf + rho * c_bar_squared;

    // 1 - phi^2 vanishes in either bulk. The recolouring step can overshoot
    // |phi| = 1 by a rounding error, which would make this negative and the
    // square root a NaN, so it is clamped rather than trusted.
    double mixing = 1.0 - phi * phi;
    if (mixing < 0.0) {
        mixing = 0.0;
    }

    const double discriminant = linear * linear + rho * rho * mixing * c1_squared * c2_squared;
    return 0.5 * (rho * c_hat_squared - p1_inf - p2_inf + std::sqrt(discriminant));
}

double pressure_linear_mixing(double rho, double phi, double cs_squared,
                              const ComponentPair& components) {
    return rho * cs_squared - 0.5 * (1.0 + phi) * components.p1_inf -
           0.5 * (1.0 - phi) * components.p2_inf;
}

}  // namespace lbm
}  // namespace cglbm
