#include "lbm/velocity_based_solver.h"

#include <algorithm>

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

}  // namespace

Solver::Solver(const SolverParameters& parameters) : parameters_(parameters) {
    const std::size_t nodes = static_cast<std::size_t>(parameters_.nx) * parameters_.ny;
    for (std::vector<double>* populations : {&g_, &h_, &g_new_, &h_new_, &dissipation_}) {
        populations->assign(nodes * kQ, 0.0);
    }
    for (std::vector<double>* field : {&c_,          &rho_,
                                       &psi_,        &pressure_number_,
                                       &ux_,         &uy_,
                                       &ax_,         &ay_,
                                       &lattice_ux_, &lattice_uy_,
                                       &body_x_,     &body_y_,
                                       &rho_old_,    &pressure_number_old_,
                                       &ux_old_,     &uy_old_,
                                       &ax_old_,     &ay_old_,
                                       &grad_x_,     &grad_y_,
                                       &normal_x_,   &normal_y_,
                                       &stress_xx_,  &stress_xy_,
                                       &stress_yy_}) {
        field->assign(nodes, 0.0);
    }
}

double Solver::dynamic_viscosity(double c) const {
    const double bounded = bounded_fraction(c);
    return bounded * parameters_.mu1 + (1.0 - bounded) * parameters_.mu2;
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
            velocity_equilibrium(s.ux, s.uy, gamma);
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
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int m = node(i, j);
            gradient_periodic(psi_.data(), nx, ny, i, j, stencil, &grad_x_[m], &grad_y_[m]);
            unit_normal(grad_x_[m], grad_y_[m], &normal_x_[m], &normal_y_[m]);
        }
    }
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int m = node(i, j);
            double dnx_dx = 0.0, dnx_dy = 0.0, dny_dx = 0.0, dny_dy = 0.0;
            gradient_periodic(normal_x_.data(), nx, ny, i, j, stencil, &dnx_dx, &dnx_dy);
            gradient_periodic(normal_y_.data(), nx, ny, i, j, stencil, &dny_dx, &dny_dy);
            const double weight = layer_weight(psi_[m], dnx_dx + dny_dy, parameters_.width);
            capillary_stress(parameters_.surface_tension * weight,
                             grad_x_[m],
                             grad_y_[m],
                             &stress_xx_[m],
                             &stress_xy_[m],
                             &stress_yy_[m]);
        }
    }
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int m = node(i, j);
            double fx = 0.0, fy = 0.0;
            surface_force(stress_xx_.data(),
                          stress_xy_.data(),
                          stress_yy_.data(),
                          nx,
                          ny,
                          i,
                          j,
                          stencil,
                          Boundary::PeriodicY,
                          &fx,
                          &fy);
            double dx = 0.0, dy = 0.0;
            dissipation_force(dissipation_.data(),
                              lattice_ux_.data(),
                              lattice_uy_.data(),
                              nx,
                              ny,
                              i,
                              j,
                              &dx,
                              &dy);
            double px = 0.0, py = 0.0;
            pressure_force(pressure_number_.data(), rho_.data(), nx, ny, i, j, &px, &py);
            ax_[m] = (fx + dx + body_x_[m]) / rho_[m] + px;
            ay_[m] = (fy + dy + body_y_[m]) / rho_[m] + py;
        }
    }
}

// The momentum after streaming: what each node held after its collision, plus
// the exchanges over its eight links. The populations have just streamed, so
// g at x along k is what the neighbour x - xi_k sent, and g at x - xi_k along
// the opposite direction what x sent back.
void Solver::momentum() {
    const int nx = parameters_.nx;
    const int ny = parameters_.ny;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int m = node(i, j);
            // post-collision velocity u + a/2 of the previous step
            double jx = rho_old_[m] * (ux_old_[m] + 0.5 * ax_old_[m]);
            double jy = rho_old_[m] * (uy_old_[m] + 0.5 * ay_old_[m]);
            dissipation_[m * kQ] = 0.0;
            for (int k = 1; k < kQ; ++k) {
                const int d =
                    node((i - kVelocity[k][0] + nx) % nx, (j - kVelocity[k][1] + ny) % ny);
                const int back = kOpposite[k];
                LinkEnd donor;
                donor.outgoing = g_[m * kQ + k] - kWeight[k] * pressure_number_old_[d];
                donor.phase = h_[m * kQ + k];
                donor.rho = rho_[d];
                donor.mu = dynamic_viscosity(c_[d]);
                donor.ux = ux_old_[d];
                donor.uy = uy_old_[d];
                LinkEnd receiver;
                receiver.outgoing = g_[d * kQ + back] - kWeight[k] * pressure_number_old_[m];
                receiver.phase = h_[d * kQ + back];
                receiver.rho = rho_[m];
                receiver.mu = dynamic_viscosity(c_[m]);
                receiver.ux = ux_old_[m];
                receiver.uy = uy_old_[m];
                double link_x = 0.0, link_y = 0.0;
                link_momentum(k,
                              donor,
                              receiver,
                              parameters_.rho1,
                              parameters_.rho2,
                              &link_x,
                              &link_y,
                              &dissipation_[m * kQ + k]);
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

// Collision of g, construction of h, and streaming of both.
void Solver::collide_and_stream() {
    const int nx = parameters_.nx;
    const int ny = parameters_.ny;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int m = node(i, j);
            double eq[kQ], source[kQ], post[kQ], phase[kQ];
            hydrodynamic_equilibrium(pressure_number_[m], ux_[m], uy_[m], eq);
            forcing(ux_[m], uy_[m], ax_[m], ay_[m], source);
            const double tau = dynamic_viscosity(c_[m]) / (rho_[m] * kCs2) + 0.5;
            VelocityGradient gradient;
            gradient_periodic(
                ux_.data(), nx, ny, i, j, GradientStencil::E4, &gradient.dux_dx, &gradient.dux_dy);
            gradient_periodic(
                uy_.data(), nx, ny, i, j, GradientStencil::E4, &gradient.duy_dx, &gradient.duy_dy);
            collide_hybrid(&g_[m * kQ],
                           eq,
                           source,
                           tau,
                           parameters_.tau_bulk,
                           parameters_.hybrid_weight,
                           gradient,
                           post);
            phase_populations(
                c_[m], ux_[m], uy_[m], normal_x_[m], normal_y_[m], parameters_.width, phase);
            for (int k = 0; k < kQ; ++k) {
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
