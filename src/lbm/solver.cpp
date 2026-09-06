#include "lbm/solver.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "lbm/isotropic_gradient.h"

namespace cglbm {
namespace lbm {

Solver::Solver(CaseConfig config)
    : config_(std::move(config)), nx_(config_.nx), ny_(config_.ny) {
    if (nx_ <= 0 || ny_ <= 0) {
        throw std::invalid_argument("lattice must have at least one node on each axis");
    }
    if (!config_.initial_phase) {
        throw std::invalid_argument("case '" + config_.name + "' has no initial phase field");
    }
    // A wall-bounded case reflects at j = 0 and j = ny - 1; with fewer than
    // three rows those are the same node and the bounce-back is meaningless.
    if (config_.boundary == Boundary::WallY && ny_ < 3) {
        throw std::invalid_argument("a wall-bounded case needs at least three rows along y");
    }

    dx_ = config_.units.dx;
    dt_ = config_.units.dt;
    cs2_ = config_.units.cs2();
    cs4_ = config_.units.cs4();
    cs6_ = config_.units.cs6();

    wall_y_ = config_.boundary == Boundary::WallY;
    parallel_ = config_.parallel;
    components_ = {config_.physics.c1 * config_.physics.c1,
                   config_.physics.c2 * config_.physics.c2, config_.physics.p1_inf,
                   config_.physics.p2_inf};

    rho_ = Field(nx_, ny_);
    rho_mdt_ = Field(nx_, ny_);
    u_ = Field(nx_, ny_, 2);
    p_ = Field(nx_, ny_);
    p_mdt_ = Field(nx_, ny_);
    phi_ = Field(nx_, ny_);
    force_ = Field(nx_, ny_, 2);

    f_ = Field(nx_, ny_, kQ);
    g_ = Field(nx_, ny_, kQ);
    f_eq_ = Field(nx_, ny_, kQ);
    omega_1_ = Field(nx_, ny_, kQ);
    omega_2_ = Field(nx_, ny_, kQ);
    omega_3_ = Field(nx_, ny_, kQ);
    source_ = Field(nx_, ny_, kQ);
}

MacroscopicState Solver::state() const {
    return MacroscopicState{&rho_, &u_, &phi_, &p_};
}

void Solver::colour_gradient(int i, int j, double* grad_x, double* grad_y) const {
    // Both stencils already carry the 1 / cs^2 factor of the colour gradient.
    if (wall_y_) {
        // A neighbour beyond a y wall contributes nothing, so the gradient
        // becomes one-sided within `stencil_reach` nodes of it.
        gradient_wall_y(phi_.data(), nx_, ny_, i, j, config_.stencil, grad_x, grad_y);
    } else {
        gradient_periodic(phi_.data(), nx_, ny_, i, j, config_.stencil, grad_x, grad_y);
    }
}

void Solver::equilibrium() {
    const bool parallel = parallel_;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            const double rho_local = rho_(i, j);
            const double u_x = u_(i, j, 0);
            const double u_y = u_(i, j, 1);
            const double p_local = p_(i, j);
            for (int k = 0; k < kQ; k++) {
                // Hermite polynomials of the discrete velocity, up to fourth
                // order: the equilibrium is their projection.
                const double H0 = 1.0;
                const double Hx = kXi[k][0];
                const double Hy = kXi[k][1];

                const double Hxx = Hx * Hx - cs2_;
                const double Hyy = Hy * Hy - cs2_;
                const double Hxy = kXi[k][0] * kXi[k][1];

                const double Hxxy = Hx * Hx * Hy - cs2_ * Hy;
                const double Hyyx = Hy * Hy * Hx - cs2_ * Hx;
                const double Hxxx = std::pow(Hx, 3) - cs2_ * 3. * Hx;
                const double Hyyy = std::pow(Hy, 3) - cs2_ * 3. * Hy;

                const double Hxxyy = Hx * Hx * Hy * Hy - cs2_ * (Hx * Hx + Hy * Hy) + cs4_;

                // E carries the part of the equilibrium that responds to
                // p - rho cs^2, i.e. to the departure from the ideal gas.
                const double E = kW[k] * ((Hxx + Hyy) / (2. * cs4_) - Hxxyy / (4. * cs6_));
                const double term1 =
                    rho_local * kW[k] *
                    (H0 + u_x * Hx / cs2_ + u_y * Hy / cs2_ +
                     0.5 * (u_x * u_x * Hxx + u_x * u_y * Hxy * 2. + u_y * u_y * Hyy) / cs4_);
                f_eq_(i, j, k) =
                    term1 + (p_local - rho_local * cs2_) *
                                (E + kW[k] * (u_x * (Hyyx + Hxxx) + u_y * (Hyyy + Hxxy)) /
                                         (2. * cs6_));
            }
        }
    }
}

void Solver::macroscopic() {
    const bool parallel = parallel_;
    const double gravity = config_.physics.gravity;
    const bool has_gravity = gravity != 0.0;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            double sum_f = 0.0;
            double sum_xi_x = 0.0;
            double sum_xi_y = 0.0;
            for (int k = 0; k < kQ; k++) {
                const double f_local = f_(i, j, k);
                sum_f += f_local;
                sum_xi_x += f_local * kXi[k][0];
                sum_xi_y += f_local * kXi[k][1];
            }
            rho_(i, j) = sum_f;
            force_(i, j, 0) = 0.0;
            force_(i, j, 1) = has_gravity ? -sum_f * gravity : 0.0;
            // Guo's forcing: half the force acts on the velocity of this step.
            u_(i, j, 0) = (sum_xi_x + force_(i, j, 0) * dt_ * 0.5) / sum_f;
            u_(i, j, 1) = (sum_xi_y + force_(i, j, 1) * dt_ * 0.5) / sum_f;
        }
    }
}

void Solver::phase_field() {
    const bool parallel = parallel_;
    const ComponentPair components = components_;
    const bool warn = config_.warn_phase_out_of_range;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            double sum_f = 0.0;
            double sum_g = 0.0;
            for (int k = 0; k < kQ; k++) {
                sum_f += f_(i, j, k);
                sum_g += g_(i, j, k);
            }
            phi_(i, j) = sum_g / sum_f;
            if (warn && std::fabs(phi_(i, j)) > 1.0) {
                std::cout << "Error : phi = " << phi_(i, j) << std::endl;
            }
            const double rho_local = rho_(i, j);
            const double phi_local = sum_g / sum_f;
            // Two-component equation of state: each component keeps its own
            // sound speed and pressure at infinity, which is what decouples the
            // density ratio from the sound-speed ratio.
            //   T. Lafarge et al., Phys. Fluids 33, 082110 (2021) -- see
            //   src/lbm/equation_of_state.h, which also keeps the linear mixing
            //   rule this replaced, for comparison.
            // Qualified: the accessor `pressure()` would otherwise hide the
            // equation of state of the same name.
            p_(i, j) = ::cglbm::lbm::pressure(rho_local, phi_local, components);
        }
    }
}

void Solver::collide() {
    const bool parallel = parallel_;
    const double nu = config_.physics.nu;
    const double nu_b = config_.physics.nu_b;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            const double rho_local = rho_(i, j);
            const double p_local = p_(i, j);
            // Relaxation times from the viscosities, p. 284 of Kruger et al.
            const double tau_nu = rho_local * nu / (p_local * dt_) + 0.5;
            const double tau_b = rho_local * nu_b / (p_local * dt_) + 0.5;
            double sum_nu_neq = 0.0, sum_b_neq = 0.0, sum_xy_neq = 0.0;

            // Project the non-equilibrium part onto the shear, bulk and
            // off-diagonal Hermite moments: f_{k,i}^{r,neq}, k in {nu, b, xy}.
            for (int k = 0; k < kQ; k++) {
                const double xi_x = kXi[k][0], xi_y = kXi[k][1];
                const double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
                const double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - cs2_;
                const double H_xy = xi_x * xi_y;

                const double f_neq = f_(i, j, k) - f_eq_(i, j, k) + 0.5 * source_(i, j, k);
                sum_nu_neq += f_neq * H_nu;
                sum_b_neq += f_neq * H_b;
                sum_xy_neq += f_neq * H_xy;
            }
            // Rebuild the regularised operator Omega^(1) from those moments,
            // relaxing shear and bulk at their own rates.
            for (int k = 0; k < kQ; k++) {
                const double xi_x = kXi[k][0], xi_y = kXi[k][1];
                const double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
                const double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - cs2_;
                const double H_xy = xi_x * xi_y;

                const double f_nu_neq = H_nu / cs4_ * sum_nu_neq;
                const double f_b_neq = H_b / cs4_ * sum_b_neq;
                const double f_xy_neq = H_xy / cs4_ * sum_xy_neq;
                omega_1_(i, j, k) = kW[k] * (1.0 - 1.0 / tau_nu) * (f_nu_neq + f_xy_neq) +
                                    kW[k] * (1.0 - 1.0 / tau_b) * f_b_neq;
            }
        }
    }
}

void Solver::force() {
    const bool parallel = parallel_;
    const bool wall = wall_y_;
    // The wall-bounded cases have always started at j = 1, leaving the source
    // term of the j = 0 row at zero. Preserved here; see docs/numerics.md.
    const int j_start = wall ? 1 : 0;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = j_start; j < ny_; j++) {
            // S_F, S_Sp and S_t are node-local: each is built and consumed
            // within this iteration. They used to be three more lattice-sized
            // arrays, which for the production Rayleigh-Taylor case was about a
            // gigabyte of state that never outlived a single node.
            double S_F[kQ];   ///< body force, Guo's forcing
            double S_Sp[kQ];  ///< third-order moment correction, unresolved on D2Q9
            double S_t[kQ];   ///< temporal correction from d(p - rho cs^2)/dt

            const double u_x = u_(i, j, 0);
            const double u_y = u_(i, j, 1);
            const double F_x = force_(i, j, 0);
            const double F_y = force_(i, j, 1);
            for (int k = 0; k < kQ; k++) {
                const double xik0 = kXi[k][0];
                const double xik1 = kXi[k][1];
                const double Hxxk = xik0 * xik0 - cs2_;
                const double Hyyk = xik1 * xik1 - cs2_;
                const double Hxyk = xik0 * xik1;

                const double term1 = (F_x * xik0 + F_y * xik1) / cs2_;
                const double term23 =
                    (u_x * F_x * Hxxk + u_y * F_y * Hyyk + (u_x * F_y + u_y * F_x) * Hxyk) / cs4_;

                S_F[k] = kW[k] * (term1 + term23);
            }

            // Lattice divergence of (p - rho cs^2) u, the quantity whose
            // third-order moment D2Q9 gets wrong.
            double derive_x = 0., derive_y = 0.;
            for (int k = 0; k < kQ; k++) {
                const int xi_x = static_cast<int>(kXi[k][0]);
                const int xi_y = static_cast<int>(kXi[k][1]);
                const int ind_x = (i + xi_x + nx_) % nx_;  // x is always periodic
                int ind_y;
                if (wall) {
                    // A neighbour beyond a wall contributes nothing. The
                    // original wrote this as `derive_x += 0.`, with no matching
                    // line for derive_y; skipping the node is the same value.
                    if (j == 0 && (k == 4 || k == 7 || k == 8)) {
                        continue;
                    }
                    if (j == ny_ - 1 && (k == 2 || k == 5 || k == 6)) {
                        continue;
                    }
                    ind_y = j + xi_y;
                } else {
                    ind_y = (j + xi_y + ny_) % ny_;
                }
                derive_x += kW[k] * xi_x * (p_(ind_x, ind_y) - rho_(ind_x, ind_y) * cs2_) *
                            u_(ind_x, ind_y, 0);
                derive_y += kW[k] * xi_y * (p_(ind_x, ind_y) - rho_(ind_x, ind_y) * cs2_) *
                            u_(ind_x, ind_y, 1);
            }
            derive_x /= (dt_ * cs2_);
            derive_y /= (dt_ * cs2_);

            for (int k = 0; k < kQ; k++) {
                const double H_nu = (kXi[k][0] * kXi[k][0] - kXi[k][1] * kXi[k][1]) / 2.;
                const double H_b = (kXi[k][0] * kXi[k][0] + kXi[k][1] * kXi[k][1]) / 2. - cs2_;
                S_Sp[k] =
                    kW[k] * (derive_y * (3. * H_nu - H_b) + derive_x * (-3. * H_nu - H_b)) /
                    (2. * cs4_);
            }

            for (int k = 0; k < kQ; k++) {
                const double H_xx = kXi[k][0] * kXi[k][0] - cs2_;
                const double H_yy = kXi[k][1] * kXi[k][1] - cs2_;
                const double H_xxyy = kXi[k][0] * kXi[k][0] * kXi[k][1] * kXi[k][1] -
                                      cs2_ * (kXi[k][0] * kXi[k][0] + kXi[k][1] * kXi[k][1]) +
                                      cs4_;
                const double E = kW[k] * ((H_xx + H_yy) / (2. * cs4_) - H_xxyy / (4. * cs6_));
                S_t[k] = (p_(i, j) - p_mdt_(i, j) - (rho_(i, j) - rho_mdt_(i, j)) * cs2_) * E;

                source_(i, j, k) = S_F[k] + S_Sp[k] + S_t[k];
            }
        }
    }
}

void Solver::collide_surface() {
    const bool parallel = parallel_;
    const double nu = config_.physics.nu;
    const double nu_b = config_.physics.nu_b;
    const double sigma = config_.physics.sigma;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            double Cx = 0.0, Cy = 0.0;
            colour_gradient(i, j, &Cx, &Cy);

            const double norm_C = std::sqrt(Cx * Cx + Cy * Cy);
            for (int k = 0; k < kQ; k++) {
                const double xi_x = kXi[k][0], xi_y = kXi[k][1];
                const double H_nu = 0.5 * (xi_x * xi_x - xi_y * xi_y);
                const double H_b = 0.5 * (xi_x * xi_x + xi_y * xi_y) - cs2_;
                const double H_xy = xi_x * xi_y;
                const double tau_nu = rho_(i, j) * nu / (p_(i, j) * dt_) + 0.5;
                const double tau_b = rho_(i, j) * nu_b / (p_(i, j) * dt_) + 0.5;
                // The operator divides by the norm of the colour gradient,
                // which vanishes away from the interface.
                if (norm_C > kGradientEpsilon) {
                    omega_2_(i, j, k) =
                        sigma * kW[k] / (4 * norm_C * cs4_) *
                        ((2 * Cx * Cy * H_xy + (Cx * Cx - Cy * Cy) * H_nu) / tau_nu -
                         ((Cx * Cx + Cy * Cy) * H_b) / tau_b);
                } else {
                    omega_2_(i, j, k) = 0.0;
                }
            }
        }
    }
}

void Solver::recolor() {
    const bool parallel = parallel_;
    const double ch_width_ope = config_.physics.ch_width_ope;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            double Cx = 0.0, Cy = 0.0;
            colour_gradient(i, j, &Cx, &Cy);

            const double grad_phi_x = Cx / dt_;
            const double grad_phi_y = Cy / dt_;
            const double norm_grad_phi =
                std::sqrt(grad_phi_x * grad_phi_x + grad_phi_y * grad_phi_y);
            if (norm_grad_phi > kGradientEpsilon) {
                // Push the components apart along the gradient, at the rate
                // that holds the interface at ch_width_ope.
                for (int k = 0; k < kQ; k++) {
                    const double xi_x = kXi[k][0], xi_y = kXi[k][1];
                    omega_3_(i, j, k) = kW[k] * p_(i, j) * (1 - phi_(i, j) * phi_(i, j)) /
                                        (2. * ch_width_ope) *
                                        (xi_x * grad_phi_x + xi_y * grad_phi_y) /
                                        (cs2_ * norm_grad_phi);
                }
            } else {
                for (int k = 0; k < kQ; k++) {
                    omega_3_(i, j, k) = 0.0;
                }
            }
        }
    }
}

void Solver::stream() {
    const bool parallel = parallel_;
    const bool wall = wall_y_;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < kQ; k++) {
                const int ip = (i + static_cast<int>(kXi[k][0]) + nx_) % nx_;
                int jp, kp;
                if (wall && j == 0 && (k == 4 || k == 7 || k == 8)) {
                    // Half-way bounce-back off the bottom wall: the population
                    // stays put and comes back along the opposite direction,
                    // which for this velocity order is k - 2.
                    jp = j;
                    kp = k - 2;
                } else if (wall && j == ny_ - 1 && (k == 2 || k == 5 || k == 6)) {
                    jp = j;
                    kp = k + 2;
                } else {
                    jp = wall ? j + static_cast<int>(kXi[k][1])
                              : (j + static_cast<int>(kXi[k][1]) + ny_) % ny_;
                    kp = k;
                }
                f_(ip, jp, kp) = f_eq_(i, j, k) + omega_1_(i, j, k) + omega_2_(i, j, k) +
                                 0.5 * source_(i, j, k);
                g_(ip, jp, kp) = f_(ip, jp, k) * phi_(i, j) + omega_3_(i, j, k);
            }
            // The temporal correction of the next step needs this step's state.
            rho_mdt_(i, j) = rho_(i, j);
            p_mdt_(i, j) = p_(i, j);
        }
    }
}

void Solver::initialize() {
    const Physics& physics = config_.physics;
    const bool has_gravity = physics.gravity != 0.0;

    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            const double phi_local = config_.initial_phase(config_, i, j);
            phi_(i, j) = phi_local;
            u_(i, j, 0) = 0.0;  // the flow starts at rest
            u_(i, j, 1) = 0.0;

            const double rho_local =
                physics.rho1 * (0.5 + 0.5 * phi_local) + physics.rho2 * (0.5 - 0.5 * phi_local);
            rho_(i, j) = rho_local;
            // At t = 0 the previous step is this one.
            rho_mdt_(i, j) = rho_local;

            // The pressure is laid down with the linear mixing rule rather than
            // with the equation of state, which is what makes the Laplace jump
            // exact at t = 0 given `matched_p1_inf`.
            const double p_local = rho_local * ((1 + phi_local) * 0.5 * physics.c1 * physics.c1 +
                                                (1 - phi_local) * 0.5 * physics.c2 * physics.c2) -
                                   (1 + phi_local) * 0.5 * physics.p1_inf -
                                   (1 - phi_local) / 2. * physics.p2_inf;
            p_(i, j) = p_local;
            p_mdt_(i, j) = p_local;

            force_(i, j, 0) = 0.0;
            force_(i, j, 1) = has_gravity ? -rho_local * physics.gravity : 0.0;
        }
    }

    equilibrium();
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < kQ; k++) {
                f_(i, j, k) = f_eq_(i, j, k);
                g_(i, j, k) = f_(i, j, k) * phi_(i, j);
            }
        }
    }
}

void Solver::step() {
    force();
    collide();
    collide_surface();
    recolor();
    stream();
    macroscopic();
    phase_field();
    equilibrium();
}

void Solver::run() {
    CsvWriter writer(config_.output_precision);

    initialize();
    writer.write_grids(0, state());
    if (config_.track_interface) {
        writer.open_interface_track();
    }

    for (int timestep = 1; timestep <= config_.steps; timestep++) {
        step();

        if (config_.track_interface) {
            writer.write_interface(timestep, phi_);
        }
        if (timestep % config_.interval == 0) {
            std::cout << "Step " << timestep << std::endl;
            writer.write_grids(timestep, state());
        }
    }
}

}  // namespace lbm
}  // namespace cglbm
