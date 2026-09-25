#include "lbm/surface_force.h"

#include <cmath>

namespace cglbm {
namespace lbm {

void unit_normal(double grad_x, double grad_y, double* normal_x, double* normal_y) {
    const double norm = std::sqrt(grad_x * grad_x + grad_y * grad_y);
    if (norm > kInterfaceGradientThreshold) {
        *normal_x = grad_x / norm;
        *normal_y = grad_y / norm;
    } else {
        *normal_x = 0.0;
        *normal_y = 0.0;
    }
}

void capillary_stress(double sigma,
                      double grad_x,
                      double grad_y,
                      double* stress_xx,
                      double* stress_xy,
                      double* stress_yy) {
    const double norm = std::sqrt(grad_x * grad_x + grad_y * grad_y);
    if (norm > kInterfaceGradientThreshold) {
        const double half_sigma = 0.5 * sigma;
        *stress_xx = half_sigma * (norm - grad_x * grad_x / norm);
        *stress_xy = -half_sigma * grad_x * grad_y / norm;
        *stress_yy = half_sigma * (norm - grad_y * grad_y / norm);
    } else {
        *stress_xx = 0.0;
        *stress_xy = 0.0;
        *stress_yy = 0.0;
    }
}

double layer_weight(double psi, double divergence_of_normal, double width) {
    // atanh is finite for |psi| < 1; past this the layer carries no stress anyway
    constexpr double kLargest = 1.0 - 1.0e-12;
    const double bounded = (psi < -kLargest) ? -kLargest : ((psi > kLargest) ? kLargest : psi);
    const double distance = -width * std::atanh(bounded);
    const double jacobian = 1.0 + distance * divergence_of_normal;
    const double kept = (jacobian < 0.25) ? 0.25 : ((jacobian > 4.0) ? 4.0 : jacobian);
    return 1.0 / kept;
}

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
                   double* force_y) {
    // F_a = d_b T_ab. Each call returns both derivatives of one component;
    // d(T_xx)/dy and d(T_yy)/dx are not needed and are discarded.
    double dxx_dx = 0.0;
    double dxx_dy = 0.0;
    double dxy_dx = 0.0;
    double dxy_dy = 0.0;
    double dyy_dx = 0.0;
    double dyy_dy = 0.0;
    gradient(stress_xx, nx, ny, i, j, stencil, boundary, &dxx_dx, &dxx_dy);
    gradient(stress_xy, nx, ny, i, j, stencil, boundary, &dxy_dx, &dxy_dy);
    gradient(stress_yy, nx, ny, i, j, stencil, boundary, &dyy_dx, &dyy_dy);
    *force_x = dxx_dx + dxy_dy;
    *force_y = dxy_dx + dyy_dy;
}

}  // namespace lbm
}  // namespace cglbm
