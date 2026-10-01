#include "lbm/velocity_based_solver.h"

#include <algorithm>
#include <stdexcept>

#include "lbm/surface_force.h"
#include "lbm/velocity_based.h"

namespace cglbm {
namespace lbm {
namespace velocity_based {

namespace {

constexpr double kCs2 = kSoundSpeedSquared;

double bounded_fraction(double c) {
    return (c < 0.0) ? 0.0 : ((c > 1.0) ? 1.0 : c);
}

/// Parities of the fields the stencils read, across a wall normal to y.
constexpr double kEven = 1.0;
constexpr double kOdd = -1.0;

}  // namespace

Solver::Solver(const SolverParameters& parameters)
    : parameters_(parameters), walls_(parameters.boundary == Boundary::WallY) {
    if (!(parameters_.phase_temperature > 0.0 && parameters_.phase_temperature < 0.6)) {
        throw std::invalid_argument("phase_temperature must lie in (0, 0.6)");
    }
    // the stencils reach four nodes, and a mirror reflects them once
    if (walls_ && parameters_.ny < 8) {
        throw std::invalid_argument("a lattice between walls must be at least 8 nodes tall");
    }
    const std::size_t nodes = static_cast<std::size_t>(parameters_.nx) * parameters_.ny;
    for (std::vector<double>* populations : {&g_, &h_, &g_new_, &h_new_, &dissipation_}) {
        populations->assign(nodes * kQ, 0.0);
    }
    // The populations start at equilibrium, and so does their non-equilibrium.
    previous_non_equilibrium_.assign(nodes * 3, 0.0);
    for (std::vector<double>* field : {&c_,
                                       &rho_,
                                       &psi_,
                                       &pressure_number_,
                                       &ux_,
                                       &uy_,
                                       &ax_,
                                       &ay_,
                                       &lattice_ux_,
                                       &lattice_uy_,
                                       &body_x_,
                                       &body_y_,
                                       &rho_old_,
                                       &pressure_number_old_,
                                       &ux_old_,
                                       &uy_old_,
                                       &ax_old_,
                                       &ay_old_,
                                       &grad_x_,
                                       &grad_y_,
                                       &normal_x_,
                                       &normal_y_,
                                       &laplacian_psi_,
                                       &laplacian_c_,
                                       &laplacian2_psi_,
                                       &laplacian2_c_,
                                       &phase_normal_x_,
                                       &phase_normal_y_,
                                       &phase_flux_x_,
                                       &phase_flux_y_,
                                       &stress_xx_,
                                       &stress_xy_,
                                       &stress_yy_}) {
        field->assign(nodes, 0.0);
    }
}

double Solver::arithmetic_viscosity(double c) const {
    const double bounded = bounded_fraction(c);
    return bounded * parameters_.mu1 + (1.0 - bounded) * parameters_.mu2;
}

double Solver::harmonic_viscosity(double c) const {
    const double bounded = bounded_fraction(c);
    return parameters_.mu1 * parameters_.mu2 /
           (bounded * parameters_.mu2 + (1.0 - bounded) * parameters_.mu1);
}

double Solver::dynamic_viscosity(double c) const {
    if (parameters_.interface_viscosity == InterfaceViscosity::Arithmetic) {
        return arithmetic_viscosity(c);
    }
    return harmonic_viscosity(c);
}

double Solver::pressure(int i, int j) const {
    const int m = node(i, j);
    return pressure_number_[m] * rho_[m] * kCs2;
}

void Solver::initialize(const std::function<NodeState(int i, int j)>& state) {
    for (int i = 0; i < parameters_.nx; ++i) {
        for (int j = 0; j < parameters_.ny; ++j) {
            const int m = node(i, j);
            const NodeState s = state(i, j);
            const double rho = parameters_.rho2 + s.c * (parameters_.rho1 - parameters_.rho2);
            double eq[kQ];
            double gamma[kQ];
            hydrodynamic_equilibrium(s.p / (rho * kCs2), s.ux, s.uy, eq);
            phase_carrier(s.ux, s.uy, parameters_.phase_temperature, gamma);
            for (int k = 0; k < kQ; ++k) {
                g_[m * kQ + k] = eq[k];
                h_[m * kQ + k] = s.c * gamma[k];
                dissipation_[m * kQ + k] = 0.0;
            }
            ux_[m] = s.ux;
            uy_[m] = s.uy;
            lattice_ux_[m] = s.ux;
            lattice_uy_[m] = s.uy;
        }
    }
    macroscopic();
    acceleration();
}

void Solver::set_body_force(
    const std::function<void(int i, int j, double* fx, double* fy)>& force) {
    for (int i = 0; i < parameters_.nx; ++i) {
        for (int j = 0; j < parameters_.ny; ++j) {
            const int m = node(i, j);
            force(i, j, &body_x_[m], &body_y_[m]);
        }
    }
}

void Solver::step() {
    save_step();
    collide_and_stream();
    macroscopic();
    momentum();
    acceleration();
    velocity();
}

// c from the phase populations, then everything that follows from it and P.
void Solver::macroscopic() {
    const int nodes = parameters_.nx * parameters_.ny;
    [[maybe_unused]] const bool parallel = parameters_.parallel;
#pragma omp parallel for if (parallel)
    for (int m = 0; m < nodes; ++m) {
        double sum_h = 0.0;
        double sum_g = 0.0;
        for (int k = 0; k < kQ; ++k) {
            sum_h += h_[m * kQ + k];
            sum_g += g_[m * kQ + k];
        }
        c_[m] = sum_h;
        const double bounded = bounded_fraction(sum_h);
        rho_[m] = parameters_.rho2 + bounded * (parameters_.rho1 - parameters_.rho2);
        psi_[m] = 2.0 * bounded - 1.0;
        pressure_number_[m] = sum_g;
    }
}

// Surface tension, pressure, link dissipation and body force, per unit mass.
void Solver::acceleration() {
    const int nx = parameters_.nx;
    const int ny = parameters_.ny;
    const GradientStencil stencil = parameters_.stencil;
    const Boundary boundary = parameters_.boundary;
    [[maybe_unused]] const bool parallel = parameters_.parallel;
#pragma omp parallel for if (parallel)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int m = node(i, j);
            field_gradient(
                psi_.data(), nx, ny, i, j, stencil, boundary, kEven, &grad_x_[m], &grad_y_[m]);
            unit_normal(grad_x_[m], grad_y_[m], &normal_x_[m], &normal_y_[m]);
        }
    }
    const bool sixth = parameters_.sixth_order_phase;
    if (parameters_.fourth_order_phase || sixth) {
#pragma omp parallel for if (parallel)
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                const int m = node(i, j);
                laplacian_psi_[m] = lattice_laplacian(psi_.data(), nx, ny, i, j, boundary);
                laplacian_c_[m] = lattice_laplacian(c_.data(), nx, ny, i, j, boundary);
            }
        }
        if (sixth) {
#pragma omp parallel for if (parallel)
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    const int m = node(i, j);
                    laplacian2_psi_[m] =
                        lattice_laplacian(laplacian_psi_.data(), nx, ny, i, j, boundary);
                    laplacian2_c_[m] =
                        lattice_laplacian(laplacian_c_.data(), nx, ny, i, j, boundary);
                }
            }
        }
        const double temperature = parameters_.phase_temperature;
#pragma omp parallel for if (parallel)
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                const int m = node(i, j);
                if (1.0 - psi_[m] * psi_[m] < kFourthOrderNormalBand) {
                    phase_normal_x_[m] = normal_x_[m];
                    phase_normal_y_[m] = normal_y_[m];
                } else if (sixth) {
                    sixth_order_normal(psi_.data(),
                                       laplacian_psi_.data(),
                                       laplacian2_psi_.data(),
                                       nx,
                                       ny,
                                       i,
                                       j,
                                       &phase_normal_x_[m],
                                       &phase_normal_y_[m],
                                       boundary);
                } else {
                    interface_normal(psi_.data(),
                                     laplacian_psi_.data(),
                                     nx,
                                     ny,
                                     i,
                                     j,
                                     &phase_normal_x_[m],
                                     &phase_normal_y_[m],
                                     boundary);
                }
                if (sixth) {
                    sixth_order_flux(laplacian_c_.data(),
                                     laplacian2_c_.data(),
                                     nx,
                                     ny,
                                     i,
                                     j,
                                     temperature,
                                     &phase_flux_x_[m],
                                     &phase_flux_y_[m],
                                     boundary);
                } else {
                    phase_correction_flux(laplacian_c_.data(),
                                          nx,
                                          ny,
                                          i,
                                          j,
                                          temperature,
                                          &phase_flux_x_[m],
                                          &phase_flux_y_[m],
                                          boundary);
                }
            }
        }
    } else {
        phase_normal_x_ = normal_x_;
        phase_normal_y_ = normal_y_;
    }
#pragma omp parallel for if (parallel)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int m = node(i, j);
            double dnx_dx = 0.0, dnx_dy = 0.0, dny_dx = 0.0, dny_dy = 0.0;
            field_gradient(
                normal_x_.data(), nx, ny, i, j, stencil, boundary, kEven, &dnx_dx, &dnx_dy);
            field_gradient(
                normal_y_.data(), nx, ny, i, j, stencil, boundary, kOdd, &dny_dx, &dny_dy);
            const double weight = layer_weight(psi_[m], dnx_dx + dny_dy, parameters_.width);
            capillary_stress(parameters_.surface_tension * weight,
                             grad_x_[m],
                             grad_y_[m],
                             &stress_xx_[m],
                             &stress_xy_[m],
                             &stress_yy_[m]);
        }
    }
    const bool gravity = parameters_.gravity_x != 0.0 || parameters_.gravity_y != 0.0;
#pragma omp parallel for if (parallel)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int m = node(i, j);
            // the divergence of the capillary stress, whose shear component is
            // odd across a wall and whose normal ones are even
            double dxx_dx = 0.0, dxx_dy = 0.0, dxy_dx = 0.0, dxy_dy = 0.0;
            double dyy_dx = 0.0, dyy_dy = 0.0;
            field_gradient(
                stress_xx_.data(), nx, ny, i, j, stencil, boundary, kEven, &dxx_dx, &dxx_dy);
            field_gradient(
                stress_xy_.data(), nx, ny, i, j, stencil, boundary, kOdd, &dxy_dx, &dxy_dy);
            field_gradient(
                stress_yy_.data(), nx, ny, i, j, stencil, boundary, kEven, &dyy_dx, &dyy_dy);
            const double fx = dxx_dx + dxy_dy;
            const double fy = dxy_dx + dyy_dy;
            double dx = 0.0, dy = 0.0;
            dissipation_force(dissipation_.data(),
                              lattice_ux_.data(),
                              lattice_uy_.data(),
                              nx,
                              ny,
                              i,
                              j,
                              &dx,
                              &dy,
                              boundary);
            double px = 0.0, py = 0.0;
            pressure_force(pressure_number_.data(), rho_.data(), nx, ny, i, j, &px, &py, boundary);
            double bx = body_x_[m];
            double by = body_y_[m];
            if (gravity) {
                const double buoyant = rho_[m] - parameters_.gravity_reference_density;
                bx += buoyant * parameters_.gravity_x;
                by += buoyant * parameters_.gravity_y;
            }
            ax_[m] = (fx + dx + bx) / rho_[m] + px;
            ay_[m] = (fy + dy + by) / rho_[m] + py;
        }
    }
}

// The momentum after streaming: what each node held after its collision, plus
// the exchanges over its eight links. The populations have just streamed, so
// g at x along k is what the neighbour x - xi_k sent, and g at x - xi_k along
// the opposite direction what x sent back. Through a wall both are what x sent
// into it, bounced back to x along k, and the neighbour is x's own image.
void Solver::momentum() {
    const int nx = parameters_.nx;
    const int ny = parameters_.ny;
    [[maybe_unused]] const bool parallel = parameters_.parallel;
#pragma omp parallel for if (parallel)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int m = node(i, j);
            // post-collision velocity u + a/2 of the previous step
            double jx = rho_old_[m] * (ux_old_[m] + 0.5 * ax_old_[m]);
            double jy = rho_old_[m] * (uy_old_[m] + 0.5 * ay_old_[m]);
            dissipation_[m * kQ] = 0.0;
            for (int k = 1; k < kQ; ++k) {
                const bool wall = through_wall(j, k);
                const int d =
                    wall ? m
                         : node((i - kVelocity[k][0] + nx) % nx, (j - kVelocity[k][1] + ny) % ny);
                const int back = kOpposite[k];
                LinkEnd receiver;
                receiver.outgoing = (wall ? g_[m * kQ + k] : g_[d * kQ + back]) -
                                    kWeight[k] * pressure_number_old_[m];
                receiver.phase = wall ? h_[m * kQ + k] : h_[d * kQ + back];
                receiver.rho = rho_[m];
                receiver.mu = dynamic_viscosity(c_[m]);
                receiver.ux = ux_old_[m];
                receiver.uy = uy_old_[m];
                LinkEnd donor;
                if (wall) {
                    donor = wall_image(receiver);
                } else {
                    donor.outgoing = g_[m * kQ + k] - kWeight[k] * pressure_number_old_[d];
                    donor.phase = h_[m * kQ + k];
                    donor.rho = rho_[d];
                    donor.mu = dynamic_viscosity(c_[d]);
                    donor.ux = ux_old_[d];
                    donor.uy = uy_old_[d];
                }
                double link_x = 0.0, link_y = 0.0;
                link_momentum(k,
                              donor,
                              receiver,
                              parameters_.rho1,
                              parameters_.rho2,
                              &link_x,
                              &link_y,
                              &dissipation_[m * kQ + k],
                              parameters_.phase_temperature);
                jx += link_x;
                jy += link_y;
            }
            lattice_ux_[m] = jx / rho_[m];
            lattice_uy_[m] = jy / rho_[m];
        }
    }
}

// The populations take the exchanged velocity; the macroscopic one adds the
// half acceleration of the forcing scheme.
void Solver::velocity() {
    const int nodes = parameters_.nx * parameters_.ny;
    [[maybe_unused]] const bool parallel = parameters_.parallel;
#pragma omp parallel for if (parallel)
    for (int m = 0; m < nodes; ++m) {
        set_velocity(&g_[m * kQ], lattice_ux_[m], lattice_uy_[m]);
        ux_[m] = lattice_ux_[m] + 0.5 * ax_[m];
        uy_[m] = lattice_uy_[m] + 0.5 * ay_[m];
    }
}

void Solver::save_step() {
    rho_old_ = rho_;
    pressure_number_old_ = pressure_number_;
    ux_old_ = ux_;
    uy_old_ = uy_;
    ax_old_ = ax_;
    ay_old_ = ay_;
}

// Collision of g, construction of h, and streaming of both. A population that
// would cross a wall comes back to its own node along the opposite direction,
// into a slot no other node streams to.
void Solver::collide_and_stream() {
    const int nx = parameters_.nx;
    const int ny = parameters_.ny;
    const bool laminate = parameters_.interface_viscosity == InterfaceViscosity::Laminate;
    [[maybe_unused]] const bool parallel = parameters_.parallel;
#pragma omp parallel for if (parallel)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int m = node(i, j);
            double eq[kQ], source[kQ], post[kQ], phase[kQ];
            hydrodynamic_equilibrium(pressure_number_[m], ux_[m], uy_[m], eq);
            forcing(ux_[m], uy_[m], ax_[m], ay_[m], source);
            const double tau = dynamic_viscosity(c_[m]) / (rho_[m] * kCs2) + 0.5;
            if (laminate) {
                const double tau_along = arithmetic_viscosity(c_[m]) / (rho_[m] * kCs2) + 0.5;
                collide_laminate(&g_[m * kQ],
                                 eq,
                                 source,
                                 tau,
                                 tau_along,
                                 parameters_.tau_bulk,
                                 parameters_.filter_weight,
                                 normal_x_[m],
                                 normal_y_[m],
                                 &previous_non_equilibrium_[m * 3],
                                 post);
            } else {
                collide_filtered(&g_[m * kQ],
                                 eq,
                                 source,
                                 tau,
                                 parameters_.tau_bulk,
                                 parameters_.filter_weight,
                                 &previous_non_equilibrium_[m * 3],
                                 post);
            }
            phase_populations(c_[m],
                              ux_[m],
                              uy_[m],
                              phase_normal_x_[m],
                              phase_normal_y_[m],
                              parameters_.width,
                              phase,
                              parameters_.phase_temperature,
                              phase_flux_x_[m],
                              phase_flux_y_[m]);
            for (int k = 0; k < kQ; ++k) {
                if (through_wall(j, kOpposite[k])) {
                    g_new_[m * kQ + kOpposite[k]] = post[k];
                    h_new_[m * kQ + kOpposite[k]] = phase[k];
                    continue;
                }
                const int target =
                    node((i + kVelocity[k][0] + nx) % nx, (j + kVelocity[k][1] + ny) % ny);
                g_new_[target * kQ + k] = post[k];
                h_new_[target * kQ + k] = phase[k];
            }
        }
    }
    std::swap(g_, g_new_);
    std::swap(h_, h_new_);
}

}  // namespace velocity_based
}  // namespace lbm
}  // namespace cglbm
