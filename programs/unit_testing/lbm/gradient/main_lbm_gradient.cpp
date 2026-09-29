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
//     divides the colour gradient by its own norm;
//  4. that matched_face_value differences to kMatchedDerivative where the
//     density is uniform, and keeps each face within its two nodes where it
//     is not.
//
//   main_lbm_gradient [stencil]

#include <algorithm>
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

/// matched_face_value on a periodic line of n nodes, at the face between
/// `a` and `a + 1`.
double face_on_line(const std::vector<double>& psi, const std::vector<double>& rho, int a) {
    const int n = static_cast<int>(psi.size());
    double line_psi[6];
    double line_rho[6];
    for (int m = 0; m < 6; ++m) {
        const int node = ((a - 2 + m) % n + n) % n;
        line_psi[m] = psi[node];
        line_rho[m] = rho[node];
    }
    return cglbm::lbm::matched_face_value(line_psi, line_rho);
}

/// The face form of the matched derivative, on three lines of 64 nodes
/// carrying psi = (p - rho c_s^2) u for a smooth u:
///
///  - one density throughout: the face difference is kMatchedDerivative;
///  - the density varying by half a per cent: still the plain interpolation;
///  - an interface of width 1.6 between densities 1000 and 1: each face whose
///    six nodes span more than a per cent in density lies between its two
///    nodes, where the plain interpolation overshoots, and the face
///    differences still sum to zero over the line.
void report_matched_face() {
    const int n = 64;
    const double two_pi = 2.0 * std::acos(-1.0);
    std::vector<double> u(n), psi(n), rho(n);
    for (int x = 0; x < n; ++x) {
        u[x] = 1e-3 * (std::sin(two_pi * x / n) + 0.3 * std::cos(3.0 * two_pi * x / n));
    }
    auto fill = [&](auto&& density) {
        for (int x = 0; x < n; ++x) {
            rho[x] = density(x);
            psi[x] = (1.0 / 3.0 - rho[x] / 3.0) * u[x];
        }
    };

    fill([](int) { return 1000.0; });
    double uniform_error = 0.0;
    for (int x = 0; x < n; ++x) {
        double derivative = 0.0;
        for (int m = 1; m <= cglbm::lbm::kMatchedReach; ++m) {
            derivative +=
                cglbm::lbm::kMatchedDerivative[m - 1] * (psi[(x + m) % n] - psi[(x - m + n) % n]);
        }
        const double faces = face_on_line(psi, rho, x) - face_on_line(psi, rho, x - 1);
        uniform_error = std::max(uniform_error, std::fabs(faces - derivative));
    }

    fill([&](int x) { return 1000.0 * (1.0 + 0.0025 * std::sin(two_pi * x / n)); });
    double gentle_error = 0.0;
    for (int x = 0; x < n; ++x) {
        const double face =
            cglbm::lbm::kMatchedFace[0] * (psi[x] + psi[(x + 1) % n]) +
            cglbm::lbm::kMatchedFace[1] * (psi[(x - 1 + n) % n] + psi[(x + 2) % n]) +
            cglbm::lbm::kMatchedFace[2] * (psi[(x - 2 + n) % n] + psi[(x + 3) % n]);
        gentle_error = std::max(gentle_error, std::fabs(face_on_line(psi, rho, x) - face));
    }

    // a heavy band between x = 16 and x = 48
    fill([](int x) {
        const double alpha = 0.5 * (std::tanh((x - 16.0) / 1.6) - std::tanh((x - 48.0) / 1.6));
        return 1.0 + 999.0 * alpha;
    });
    double scale = 0.0;
    for (int x = 0; x < n; ++x) {
        scale = std::max(scale, std::fabs(psi[x]));
    }
    double overshoot = 0.0;
    double outside = 0.0;
    double sum = 0.0;
    for (int x = 0; x < n; ++x) {
        const double lo = std::min(psi[x], psi[(x + 1) % n]);
        const double hi = std::max(psi[x], psi[(x + 1) % n]);
        const double plain =
            cglbm::lbm::kMatchedFace[0] * (psi[x] + psi[(x + 1) % n]) +
            cglbm::lbm::kMatchedFace[1] * (psi[(x - 1 + n) % n] + psi[(x + 2) % n]) +
            cglbm::lbm::kMatchedFace[2] * (psi[(x - 2 + n) % n] + psi[(x + 3) % n]);
        const double face = face_on_line(psi, rho, x);
        sum += face - face_on_line(psi, rho, x - 1);
        double lightest = rho[x];
        double heaviest = rho[x];
        for (int m = -2; m <= 3; ++m) {
            lightest = std::min(lightest, rho[(x + m + n) % n]);
            heaviest = std::max(heaviest, rho[(x + m + n) % n]);
        }
        if (heaviest <= cglbm::lbm::kMatchedFaceDensityVariation * lightest) {
            continue;  // the heavy band's interior, which keeps the plain value
        }
        overshoot = std::max(overshoot, std::max(lo - plain, plain - hi) / scale);
        outside = std::max(outside, std::max(lo - face, face - hi) / scale);
    }

    std::cout << "matched_face_uniform_error = " << uniform_error / scale << "\n";
    std::cout << "matched_face_gentle_error = " << gentle_error << "\n";
    std::cout << "matched_face_plain_overshoot = " << overshoot << "\n";
    std::cout << "matched_face_outside = " << outside << "\n";
    std::cout << "matched_face_sum = " << sum / scale << "\n";
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
    report_matched_face();

    return 0;
}
