// Checks src/lbm/mixture and src/lbm/surface_force.
//
// Everything is reported as `key = value` lines, parsed by the pytest beside
// this file:
//
//  1. the initial state built from a volume fraction is in equilibrium with
//     the equation of state, and gives back its volume fraction;
//  2. the phi = 0 contour of a tanh volume-fraction profile sits
//     (W/2) ln(rho_1/rho_2) outside the density interface;
//  3. the mixture viscosity reduces to a single viscosity and mixes the
//     dynamic viscosities by volume;
//  4. the surface force of a circular interface integrates to the Laplace
//     jump sigma/R, sums to zero over any closed interface, and vanishes on a
//     flat one.
//
//   main_lbm_mixture [density_ratio] [stencil]

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "lbm/equation_of_state.h"
#include "lbm/isotropic_gradient.h"
#include "lbm/mixture.h"
#include "lbm/surface_force.h"

namespace {

using cglbm::lbm::ComponentPair;
using cglbm::lbm::GradientStencil;
using cglbm::lbm::MixtureState;

const double kSoundSpeedSquared = 1.0 / 3.0;
const double kLightDensity = 1.0;
const double kLaplaceJump = 0.0276834;  // sigma / R of the laplace case
const double kWidth = 1.6;              // interface width of the solvers

/// The components of the Laplace case at a given density ratio.
ComponentPair laplace_components(double density_ratio) {
    const double heavy = density_ratio * kLightDensity;
    return {kSoundSpeedSquared,
            kSoundSpeedSquared,
            heavy * kSoundSpeedSquared - kLightDensity * kSoundSpeedSquared - kLaplaceJump,
            0.0};
}

double ambient_pressure() {
    return kLightDensity * kSoundSpeedSquared;
}

/// Largest error of the round trip alpha -> (rho, phi) -> (p, alpha, psi).
void report_round_trip(const ComponentPair& components) {
    const double alphas[] = {0.0, 1.0e-6, 1.0e-3, 0.1, 0.5, 0.9, 0.999, 1.0 - 1.0e-6, 1.0};
    const double pressures[] = {ambient_pressure(), ambient_pressure() + kLaplaceJump};

    double pressure_error = 0.0;
    double alpha_error = 0.0;
    double psi_error = 0.0;
    for (double p : pressures) {
        for (double alpha : alphas) {
            const MixtureState state =
                cglbm::lbm::mixture_from_volume_fraction(alpha, p, components);
            const double p_back = cglbm::lbm::pressure(state.rho, state.phi, components);
            pressure_error = std::max(pressure_error, std::fabs(p_back - p) / p);
            const double alpha_back =
                cglbm::lbm::volume_fraction(state.rho, state.phi, p, components);
            alpha_error = std::max(alpha_error, std::fabs(alpha_back - alpha));
            const double psi = cglbm::lbm::normalised_phase(state.phi, p, components);
            psi_error = std::max(psi_error, std::fabs(psi - (2.0 * alpha - 1.0)));
        }
    }
    std::cout << "round_trip_pressure_error = " << pressure_error << "\n";
    std::cout << "round_trip_alpha_error = " << alpha_error << "\n";
    std::cout << "round_trip_psi_error = " << psi_error << "\n";

    // the two bulks are the pure components
    const double p = ambient_pressure();
    const MixtureState pure1 = cglbm::lbm::mixture_from_volume_fraction(1.0, p, components);
    const MixtureState pure2 = cglbm::lbm::mixture_from_volume_fraction(0.0, p, components);
    std::cout << "pure1_phi = " << pure1.phi << "\n";
    std::cout << "pure2_phi = " << pure2.phi << "\n";
    std::cout << "pure1_rho_error = "
              << std::fabs(pure1.rho - cglbm::lbm::component1_density(p, components)) << "\n";
    std::cout << "pure2_rho_error = "
              << std::fabs(pure2.rho - cglbm::lbm::component2_density(p, components)) << "\n";

    // the recolouring can push phi past +-1 by a rounding error; psi must not follow
    double overshoot = 0.0;
    for (double phi : {1.0 + 1.0e-9, -1.0 - 1.0e-9, 1.0 + 1.0e-3, -1.0 - 1.0e-3}) {
        const double psi = cglbm::lbm::normalised_phase(phi, p, components);
        overshoot = std::max(overshoot, std::fabs(psi) - 1.0);
    }
    std::cout << "psi_overshoot = " << overshoot << "\n";
}

/// Where phi crosses zero across a tanh volume-fraction profile of width W.
///
/// The profile is alpha_1(x) = (1 - tanh(x / W)) / 2, so the density
/// interface is at x = 0 and component 1 is at negative x. phi(x) is built
/// with the library, and its zero is found by bisection.
void report_phase_offset(const ComponentPair& components) {
    const double p = ambient_pressure();
    auto phi_at = [&](double x) {
        const double alpha = 0.5 * (1.0 - std::tanh(x / kWidth));
        return cglbm::lbm::mixture_from_volume_fraction(alpha, p, components).phi;
    };
    double lo = -50.0;
    double hi = 50.0;
    for (int n = 0; n < 200; ++n) {
        const double mid = 0.5 * (lo + hi);
        if (phi_at(mid) > 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    const double ratio = cglbm::lbm::component1_density(p, components) /
                         cglbm::lbm::component2_density(p, components);
    std::cout << "phase_zero_offset = " << 0.5 * (lo + hi) << "\n";
    std::cout << "phase_zero_offset_expected = " << 0.5 * kWidth * std::log(ratio) << "\n";
}

void report_viscosity(const ComponentPair& components) {
    const double nu = 1.66383;
    double single_error = 0.0;
    for (double phi : {-1.2, -1.0, -0.5, 0.0, 0.3, 1.0, 1.2}) {
        single_error = std::max(
            single_error, std::fabs(cglbm::lbm::mixture_kinematic_viscosity(phi, nu, nu) - nu));
    }
    std::cout << "viscosity_single_error = " << single_error << "\n";

    // rho nu_mix = alpha_1 mu_1 + alpha_2 mu_2, on a state halfway through
    const double p = ambient_pressure();
    const double rho1 = cglbm::lbm::component1_density(p, components);
    const double rho2 = cglbm::lbm::component2_density(p, components);
    const double nu1 = 1.0e-3;
    const double nu2 = 0.7;
    double mixing_error = 0.0;
    for (double alpha : {0.0, 0.01, 0.5, 0.99, 1.0}) {
        const MixtureState state = cglbm::lbm::mixture_from_volume_fraction(alpha, p, components);
        const double mu = state.rho * cglbm::lbm::mixture_kinematic_viscosity(state.phi, nu1, nu2);
        const double expected = alpha * rho1 * nu1 + (1.0 - alpha) * rho2 * nu2;
        mixing_error = std::max(mixing_error, std::fabs(mu - expected) / expected);
    }
    std::cout << "viscosity_mixing_error = " << mixing_error << "\n";
}

/// Surface force of a psi field: the stress at every node, then its divergence.
/// With `weighted`, each layer's stress carries layer_weight.
void force_field(const std::vector<double>& psi,
                 int n,
                 double sigma,
                 GradientStencil stencil,
                 cglbm::lbm::Boundary boundary,
                 std::vector<double>* force_x,
                 std::vector<double>* force_y,
                 bool weighted = false) {
    const std::size_t size = static_cast<std::size_t>(n) * n;
    std::vector<double> gx(size);
    std::vector<double> gy(size);
    std::vector<double> normal_x(size);
    std::vector<double> normal_y(size);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const std::size_t k = static_cast<std::size_t>(i) * n + j;
            cglbm::lbm::gradient(psi.data(), n, n, i, j, stencil, boundary, &gx[k], &gy[k]);
            cglbm::lbm::unit_normal(gx[k], gy[k], &normal_x[k], &normal_y[k]);
        }
    }
    std::vector<double> sxx(size);
    std::vector<double> sxy(size);
    std::vector<double> syy(size);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const std::size_t k = static_cast<std::size_t>(i) * n + j;
            double weight = 1.0;
            if (weighted) {
                double dnx_dx = 0.0, dnx_dy = 0.0, dny_dx = 0.0, dny_dy = 0.0;
                cglbm::lbm::gradient(
                    normal_x.data(), n, n, i, j, stencil, boundary, &dnx_dx, &dnx_dy);
                cglbm::lbm::gradient(
                    normal_y.data(), n, n, i, j, stencil, boundary, &dny_dx, &dny_dy);
                weight = cglbm::lbm::layer_weight(psi[k], dnx_dx + dny_dy, kWidth);
            }
            cglbm::lbm::capillary_stress(sigma * weight, gx[k], gy[k], &sxx[k], &sxy[k], &syy[k]);
        }
    }
    force_x->assign(size, 0.0);
    force_y->assign(size, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const std::size_t k = static_cast<std::size_t>(i) * n + j;
            cglbm::lbm::surface_force(sxx.data(),
                                      sxy.data(),
                                      syy.data(),
                                      n,
                                      n,
                                      i,
                                      j,
                                      stencil,
                                      boundary,
                                      &(*force_x)[k],
                                      &(*force_y)[k]);
        }
    }
}

/// The layer weight is r/R on every layer of a circle, and 1 on a flat interface.
void report_layer_weight() {
    const double radius = 10.0;
    double circle = 0.0;
    double flat = 0.0;
    for (double r = radius - 3.0 * kWidth; r <= radius + 3.0 * kWidth; r += 0.1) {
        const double psi = -std::tanh((r - radius) / kWidth);
        // the normal points into component 1, towards the centre: div n = -1/r
        const double weight = cglbm::lbm::layer_weight(psi, -1.0 / r, kWidth);
        circle = std::max(circle, std::fabs(weight - r / radius));
        flat = std::max(flat, std::fabs(cglbm::lbm::layer_weight(psi, 0.0, kWidth) - 1.0));
    }
    std::cout << "layer_weight_circle_error = " << circle << "\n";
    std::cout << "layer_weight_flat_error = " << flat << "\n";
}

/// Pressure jump the force of a psi circle supports: p_out - p_in is the sum of
/// F_x along the +x ray, since grad p = F at rest.
double ray_jump(const std::vector<double>& fx, int n) {
    double ray_sum = 0.0;
    for (int i = n / 2; i < n; ++i) {
        ray_sum += fx[static_cast<std::size_t>(i) * n + n / 2];
    }
    return -ray_sum;
}

/// Integrated force of a circle, net force of an ellipse, force on a flat interface.
void report_surface_force(GradientStencil stencil) {
    const int n = 128;
    const double radius = 20.0;
    const double sigma = 0.3;
    const double centre = 0.5 * n;
    const std::size_t size = static_cast<std::size_t>(n) * n;
    std::vector<double> psi(size);
    std::vector<double> fx;
    std::vector<double> fy;

    // circle of radius R: the force summed across the interface is sigma/R
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double r = std::hypot(i - centre, j - centre);
            psi[static_cast<std::size_t>(i) * n + j] = -std::tanh((r - radius) / kWidth);
        }
    }
    force_field(psi, n, sigma, stencil, cglbm::lbm::Boundary::Periodic, &fx, &fy);
    // jump: p_out - p_in = sum of F_x along the +x ray, since grad p = F at rest
    double ray_sum = 0.0;
    int inward = 0;
    int outward = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const std::size_t k = static_cast<std::size_t>(i) * n + j;
            const double dx = i - centre;
            const double dy = j - centre;
            if (std::fabs(std::hypot(dx, dy) - radius) < 2.0 * kWidth) {
                // the force must point towards the centre of the droplet
                if (fx[k] * dx + fy[k] * dy < 0.0) {
                    ++inward;
                } else {
                    ++outward;
                }
            }
            if (j == n / 2 && i >= n / 2) {
                ray_sum += fx[k];
            }
        }
    }
    std::cout << "force_jump_over_laplace = " << -ray_sum / (sigma / radius) << "\n";
    std::cout << "force_inward_nodes = " << inward << "\n";
    std::cout << "force_outward_nodes = " << outward << "\n";

    // the same circles with each layer weighted: sigma / R at R = 20 and 10
    for (double r : {20.0, 10.0}) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                psi[static_cast<std::size_t>(i) * n + j] =
                    -std::tanh((std::hypot(i - centre, j - centre) - r) / kWidth);
            }
        }
        force_field(psi, n, sigma, stencil, cglbm::lbm::Boundary::Periodic, &fx, &fy, false);
        std::cout << "force_jump_over_laplace_R" << r << " = " << ray_jump(fx, n) / (sigma / r)
                  << "\n";
        force_field(psi, n, sigma, stencil, cglbm::lbm::Boundary::Periodic, &fx, &fy, true);
        std::cout << "weighted_force_jump_over_laplace_R" << r << " = "
                  << ray_jump(fx, n) / (sigma / r) << "\n";
    }

    // an off-centre ellipse has no symmetry to cancel the force: a divergence
    // still sums to zero over the periodic lattice
    double total_x = 0.0;
    double total_y = 0.0;
    double scale = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double x = (i - 0.37 * n) / 25.0;
            const double y = (j - 0.61 * n) / 13.0;
            const double r = 13.0 * (std::hypot(x, y) - 1.0);
            psi[static_cast<std::size_t>(i) * n + j] = -std::tanh(r / kWidth);
        }
    }
    force_field(psi, n, sigma, stencil, cglbm::lbm::Boundary::Periodic, &fx, &fy);
    for (std::size_t k = 0; k < size; ++k) {
        total_x += fx[k];
        total_y += fy[k];
        scale += std::hypot(fx[k], fy[k]);
    }
    std::cout << "ellipse_net_force = " << std::hypot(total_x, total_y) / scale << "\n";
    force_field(psi, n, sigma, stencil, cglbm::lbm::Boundary::Periodic, &fx, &fy, true);
    total_x = total_y = scale = 0.0;
    for (std::size_t k = 0; k < size; ++k) {
        total_x += fx[k];
        total_y += fy[k];
        scale += std::hypot(fx[k], fy[k]);
    }
    std::cout << "weighted_ellipse_net_force = " << std::hypot(total_x, total_y) / scale << "\n";

    // a flat interface bounded by walls exerts no force at all
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            psi[static_cast<std::size_t>(i) * n + j] = -std::tanh((j - centre) / kWidth);
        }
    }
    force_field(psi, n, sigma, stencil, cglbm::lbm::Boundary::WallY, &fx, &fy);
    double flat = 0.0;
    for (std::size_t k = 0; k < size; ++k) {
        flat = std::max(flat, std::hypot(fx[k], fy[k]));
    }
    std::cout << "flat_force = " << flat << "\n";
    force_field(psi, n, sigma, stencil, cglbm::lbm::Boundary::WallY, &fx, &fy, true);
    flat = 0.0;
    for (std::size_t k = 0; k < size; ++k) {
        flat = std::max(flat, std::hypot(fx[k], fy[k]));
    }
    std::cout << "weighted_flat_force = " << flat << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    double density_ratio = 20.0;
    if (argc > 1) {
        density_ratio = std::atof(argv[1]);
        if (!(density_ratio > 0.0)) {
            std::cerr << "Invalid density ratio '" << argv[1] << "'" << std::endl;
            return 2;
        }
    }
    GradientStencil stencil = GradientStencil::E8;
    if (argc > 2 && !cglbm::lbm::stencil_from_name(argv[2], &stencil)) {
        std::cerr << "Unknown stencil '" << argv[2] << "'" << std::endl;
        return 2;
    }

    const ComponentPair components = laplace_components(density_ratio);

    std::cout.precision(12);
    std::cout << "density_ratio = " << density_ratio << "\n";
    std::cout << "stencil = " << cglbm::lbm::stencil_name(stencil) << "\n";
    report_round_trip(components);
    report_phase_offset(components);
    report_viscosity(components);
    report_layer_weight();
    report_surface_force(stencil);
    std::cout << std::flush;
    return 0;
}
