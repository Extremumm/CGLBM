// Checks the isotropic gradient stencils of src/lbm.
//
// Three properties are reported, all as `key = value` lines the pytest beside
// this file parses:
//
//  1. the lattice tensor sum_l W c_a c_b must equal delta_ab, which makes the
//     stencil exact for a linear field;
//  2. the isotropy defect of the rank-4, rank-6 and rank-8 lattice tensors,
//     which is what separates E4, E6 and E8;
//  3. the angular error of the gradient of a radial interface profile -- the
//     quantity that actually matters, since the surface-tension operator
//     divides the colour gradient by its own norm.
//
//   main_lbm_gradient [stencil]

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "lbm/isotropic_gradient.h"

namespace {

using cglbm::lbm::GradientStencil;
using cglbm::lbm::StencilPoint;

/// sum_l W c_x^px c_y^py over the stencil.
double moment(GradientStencil stencil, int px, int py) {
    int count = 0;
    const StencilPoint* points = cglbm::lbm::stencil_points(stencil, &count);
    double total = 0.0;
    for (int n = 0; n < count; ++n) {
        total += points[n].weight * std::pow(static_cast<double>(points[n].cx), px) *
                 std::pow(static_cast<double>(points[n].cy), py);
    }
    return total;
}

/// Relative departure from isotropy of the rank-`order` lattice tensor.
///
/// In two dimensions an isotropic tensor of rank 2n satisfies
/// sum W c_x^2n = (2n-1)!! / (2n-3)!! * ... ; the ratios below are the standard
/// pairwise conditions, which are simpler to check and equivalent:
///     rank 4: <x^4> = 3 <x^2 y^2>
///     rank 6: <x^6> = 5 <x^4 y^2>
///     rank 8: <x^8> = 7 <x^6 y^2>
double isotropy_defect(GradientStencil stencil, int order) {
    double lhs = 0.0;
    double rhs = 0.0;
    if (order == 4) {
        lhs = moment(stencil, 4, 0);
        rhs = 3.0 * moment(stencil, 2, 2);
    } else if (order == 6) {
        lhs = moment(stencil, 6, 0);
        rhs = 5.0 * moment(stencil, 4, 2);
    } else {
        lhs = moment(stencil, 8, 0);
        rhs = 7.0 * moment(stencil, 6, 2);
    }
    const double scale = std::fabs(lhs) + std::fabs(rhs);
    return (scale > 0.0) ? std::fabs(lhs - rhs) / scale : 0.0;
}

/// Largest angular error, in degrees, of the gradient of a radial interface.
///
/// The field is the tanh profile a colour-gradient droplet actually settles
/// into. For a radial field the gradient must point exactly along the radius;
/// whatever angle it picks up instead is lattice anisotropy, and it feeds
/// straight into the surface tension through Omega^(2).
double radial_angle_error(GradientStencil stencil, int n, double radius, double width) {
    std::vector<double> field(static_cast<std::size_t>(n) * n, 0.0);
    const double centre = 0.5 * n;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double dx = i - centre;
            const double dy = j - centre;
            const double r = std::sqrt(dx * dx + dy * dy);
            field[static_cast<std::size_t>(i) * n + j] = -std::tanh((r - radius) / width);
        }
    }

    double worst = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double dx = i - centre;
            const double dy = j - centre;
            const double r = std::sqrt(dx * dx + dy * dy);
            // only sample the interface band, where the gradient is meaningful
            if (r < radius - 2.0 * width || r > radius + 2.0 * width || r < 1.0) {
                continue;
            }
            double gx = 0.0;
            double gy = 0.0;
            cglbm::lbm::gradient_periodic(field.data(), n, n, i, j, stencil, &gx, &gy);
            const double norm = std::sqrt(gx * gx + gy * gy);
            if (norm < 1.0e-12) {
                continue;
            }
            // angle between -grad(phi) and the outward radius
            const double cosine = -(gx * dx + gy * dy) / (norm * r);
            const double clamped = (cosine > 1.0) ? 1.0 : ((cosine < -1.0) ? -1.0 : cosine);
            const double degrees = std::acos(clamped) * 180.0 / M_PI;
            if (degrees > worst) {
                worst = degrees;
            }
        }
    }
    return worst;
}

/// Largest error of the gradient of the linear field phi = a x + b y.
double linear_field_error(GradientStencil stencil, int n, double a, double b) {
    std::vector<double> field(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            field[static_cast<std::size_t>(i) * n + j] = a * i + b * j;
        }
    }
    double worst = 0.0;
    const int reach = cglbm::lbm::stencil_reach(stencil);
    // away from the edges: the field is linear, not periodic, so the wrap at
    // the boundary is not part of what is being tested here
    for (int i = reach; i < n - reach; ++i) {
        for (int j = reach; j < n - reach; ++j) {
            double gx = 0.0;
            double gy = 0.0;
            cglbm::lbm::gradient_periodic(field.data(), n, n, i, j, stencil, &gx, &gy);
            worst = std::max(worst, std::max(std::fabs(gx - a), std::fabs(gy - b)));
        }
    }
    return worst;
}

}  // namespace

int main(int argc, char** argv) {
    GradientStencil stencil = GradientStencil::E4;
    if (argc > 1 && !cglbm::lbm::stencil_from_name(argv[1], &stencil)) {
        std::cerr << "Unknown stencil '" << argv[1] << "'" << std::endl;
        return 2;
    }

    int count = 0;
    cglbm::lbm::stencil_points(stencil, &count);

    std::cout.precision(12);
    std::cout << "stencil = " << cglbm::lbm::stencil_name(stencil) << "\n";
    std::cout << "points = " << count << "\n";
    std::cout << "reach = " << cglbm::lbm::stencil_reach(stencil) << "\n";

    // normalisation: sum W c_a c_b = delta_ab
    std::cout << "moment_xx = " << moment(stencil, 2, 0) << "\n";
    std::cout << "moment_yy = " << moment(stencil, 0, 2) << "\n";
    std::cout << "moment_xy = " << moment(stencil, 1, 1) << "\n";
    // and the odd moments must vanish, or the stencil would not be centred
    std::cout << "moment_x = " << moment(stencil, 1, 0) << "\n";
    std::cout << "moment_y = " << moment(stencil, 0, 1) << "\n";

    std::cout << "isotropy_defect_4 = " << isotropy_defect(stencil, 4) << "\n";
    std::cout << "isotropy_defect_6 = " << isotropy_defect(stencil, 6) << "\n";
    std::cout << "isotropy_defect_8 = " << isotropy_defect(stencil, 8) << "\n";

    std::cout << "linear_error = " << linear_field_error(stencil, 32, 0.37, -0.11) << "\n";
    std::cout << "radial_angle_error_deg = " << radial_angle_error(stencil, 128, 20.0, 1.6)
              << std::endl;

    return 0;
}
