#include "lbm/two_population_solver.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "lbm/isotropic_gradient.h"
#include "lbm/mrt.h"

namespace cglbm {
namespace lbm {

namespace {

/// Squared lattice speed of each direction, 0, 1 or 2 on D2Q9.
inline double speed_squared(int k) {
    return kXi[k][0] * kXi[k][0] + kXi[k][1] * kXi[k][1];
}

}  // namespace

TwoPopulationSolver::TwoPopulationSolver(CaseConfig config)
    : config_(std::move(config)), nx_(config_.nx), ny_(config_.ny) {
    if (nx_ <= 0 || ny_ <= 0) {
        throw std::invalid_argument("lattice must have at least one node on each axis");
    }
    if (!config_.initial_phase) {
        throw std::invalid_argument("case '" + config_.name + "' has no initial phase field");
    }
    if (config_.boundary == Boundary::WallY && ny_ < 3) {
        throw std::invalid_argument("a wall-bounded case needs at least three rows along y");
    }
    const Physics& physics = config_.physics;
    if (physics.rho1 <= 0.0 || physics.rho2 <= 0.0) {
        throw std::invalid_argument("both densities must be positive");
    }

    dt_ = config_.units.dt;
    cs2_ = config_.units.cs2();
    wall_y_ = config_.boundary == Boundary::WallY;
    parallel_ = config_.parallel;
    mrt_ = config_.collision == Collision::MRT;
    correction_ = config_.third_moment_correction;
    const double rates[3] = {physics.s_e, physics.s_eps, physics.s_q};
    for (double rate : rates) {
        if (mrt_ && !(rate > 0.0 && rate < 2.0)) {
            throw std::invalid_argument("MRT relaxation rates must lie in (0, 2)");
        }
    }

    // The density ratio is carried by the rest weights:
    //     rho1 / rho2 = (1 - alpha_2) / (1 - alpha_1).
    alpha2_ = physics.alpha2;
    alpha1_ = 1.0 - (1.0 - alpha2_) * physics.rho2 / physics.rho1;
    if (!(alpha1_ >= 0.0 && alpha1_ <= 1.0) || !(alpha2_ >= 0.0 && alpha2_ <= 1.0)) {
        throw std::invalid_argument(
            "alpha must lie in [0, 1] for both fluids; lower physics.alpha2 for this density "
            "ratio");
    }

    const double alpha[2] = {alpha1_, alpha2_};
    for (int fluid = 0; fluid < 2; fluid++) {
        // (c_s^k)^2 = 3/5 (1 - alpha_k), Ba et al. after Eq. (6).
        cs_squared_[fluid] = 0.6 * (1.0 - alpha[fluid]);
        phi_rest_[fluid] = alpha[fluid];
        phi_near_[fluid] = (1.0 - alpha[fluid]) / 5.0;
        phi_diag_[fluid] = (1.0 - alpha[fluid]) / 20.0;
    }
    // Bulk dynamic viscosities. A negative nu2 means "the same kinematic
    // viscosity as fluid 1", as it does for `Solver`.
    const double nu2 = physics.nu2 < 0.0 ? physics.nu : physics.nu2;
    mu_[0] = physics.rho1 * physics.nu;
    mu_[1] = physics.rho2 * nu2;

    rho_1_ = Field(nx_, ny_);
    rho_2_ = Field(nx_, ny_);
    rho_ = Field(nx_, ny_);
    p_ = Field(nx_, ny_);
    u_ = Field(nx_, ny_, 2);
    phi_n_ = Field(nx_, ny_);
    force_ = Field(nx_, ny_, 2);
    grad_x_ = Field(nx_, ny_);
    grad_y_ = Field(nx_, ny_);
    normal_x_ = Field(nx_, ny_);
    normal_y_ = Field(nx_, ny_);
    if (correction_) {
        q_x_ = Field(nx_, ny_);
        q_y_ = Field(nx_, ny_);
    }
    f1_ = Field(nx_, ny_, kQ);
    f2_ = Field(nx_, ny_, kQ);
}

MacroscopicState TwoPopulationSolver::state() const {
    return MacroscopicState{&rho_, &u_, &phi_n_, &p_};
}

void TwoPopulationSolver::equilibrium(
    int fluid, double rho_k, double u_x, double u_y, double* out) const {
    // Ba et al. Eq. (14): the standard alpha_k equilibrium with the enhanced
    // first-order term, which repairs the third-order velocity moment when the
    // two fluids' sound speeds differ. The bracket vanishes for a fluid whose
    // sound speed is the lattice one.
    const double enhancement = 0.5 * (3.0 * cs_squared_[fluid] - 1.0);
    const double u_squared = u_x * u_x + u_y * u_y;
    for (int k = 0; k < kQ; k++) {
        const double e_squared = speed_squared(k);
        const double phi_i = k == 0            ? phi_rest_[fluid]
                             : e_squared < 1.5 ? phi_near_[fluid]
                                               : phi_diag_[fluid];
        const double eu = kXi[k][0] * u_x + kXi[k][1] * u_y;
        const double first = 3.0 * eu * (1.0 + enhancement * (3.0 * e_squared - 4.0));
        out[k] = rho_k * phi_i + rho_k * kW[k] * (first + 4.5 * eu * eu - 1.5 * u_squared);
    }
}

void TwoPopulationSolver::densities() {
    const bool parallel = parallel_;
    const double rho1_bulk = config_.physics.rho1;
    const double rho2_bulk = config_.physics.rho2;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            double sum1 = 0.0, sum2 = 0.0;
            for (int k = 0; k < kQ; k++) {
                sum1 += f1_(i, j, k);
                sum2 += f2_(i, j, k);
            }
            rho_1_(i, j) = sum1;
            rho_2_(i, j) = sum2;
            rho_(i, j) = sum1 + sum2;
            // p = sum_k rho_k (c_s^k)^2. The two bulk pressures are equal by
            // construction, which is what keeps the pressure continuous across
            // the interface at any density ratio.
            p_(i, j) = sum1 * cs_squared_[0] + sum2 * cs_squared_[1];
            // Ba et al. Eq. (21). Each fluid is measured against its own bulk
            // density, so the zero contour is where they occupy equal volume.
            const double scaled1 = sum1 / rho1_bulk;
            const double scaled2 = sum2 / rho2_bulk;
            const double total = scaled1 + scaled2;
            phi_n_(i, j) = total > 0.0 ? (scaled1 - scaled2) / total : 0.0;
        }
    }
}

void TwoPopulationSolver::gradient_at(
    const double* field, int i, int j, double* grad_x, double* grad_y) const {
    if (wall_y_) {
        gradient_wall_y(field, nx_, ny_, i, j, config_.stencil, grad_x, grad_y);
    } else {
        gradient_periodic(field, nx_, ny_, i, j, config_.stencil, grad_x, grad_y);
    }
}

void TwoPopulationSolver::update_colour_gradient() {
    const bool parallel = parallel_;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            double gx = 0.0, gy = 0.0;
            gradient_at(phi_n_.data(), i, j, &gx, &gy);
            grad_x_(i, j) = gx;
            grad_y_(i, j) = gy;
            const double norm = std::sqrt(gx * gx + gy * gy);
            if (norm > kInterfaceGradientFloor) {
                normal_x_(i, j) = -gx / norm;
                normal_y_(i, j) = -gy / norm;
            } else {
                normal_x_(i, j) = 0.0;
                normal_y_(i, j) = 0.0;
            }
        }
    }
}

void TwoPopulationSolver::surface_force() {
    const bool parallel = parallel_;
    const double sigma = config_.physics.sigma;
    const double gravity = config_.physics.gravity;
    const bool has_gravity = gravity != 0.0;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            // Ba et al. Eqs. (24)-(26): F_s = -1/2 sigma K grad phi_N, with the
            // curvature written so it stays bounded where the discrete normal
            // is not exactly a unit vector.
            double dnx_dx = 0.0, dnx_dy = 0.0, dny_dx = 0.0, dny_dy = 0.0;
            gradient_at(normal_x_.data(), i, j, &dnx_dx, &dnx_dy);
            gradient_at(normal_y_.data(), i, j, &dny_dx, &dny_dy);
            const double nx = normal_x_(i, j);
            const double ny = normal_y_(i, j);
            const double curvature =
                nx * ny * (dnx_dy + dny_dx) - nx * nx * dny_dy - ny * ny * dnx_dx;

            const double gx = grad_x_(i, j);
            const double gy = grad_y_(i, j);
            const bool on_interface = std::sqrt(gx * gx + gy * gy) > kInterfaceGradientFloor;
            force_(i, j, 0) = on_interface ? -0.5 * sigma * curvature * gx : 0.0;
            force_(i, j, 1) = on_interface ? -0.5 * sigma * curvature * gy : 0.0;
            if (has_gravity) {
                force_(i, j, 1) -= rho_(i, j) * gravity;
            }
        }
    }
}

void TwoPopulationSolver::update_velocity() {
    const bool parallel = parallel_;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            double sum_x = 0.0, sum_y = 0.0;
            for (int k = 0; k < kQ; k++) {
                const double total = f1_(i, j, k) + f2_(i, j, k);
                sum_x += total * kXi[k][0];
                sum_y += total * kXi[k][1];
            }
            // Ba et al. Eq. (29): half the capillary force belongs to the
            // velocity, the other half to Guo's forcing in `collide`.
            const double density = rho_(i, j);
            u_(i, j, 0) = (sum_x + 0.5 * force_(i, j, 0) * dt_) / density;
            u_(i, j, 1) = (sum_y + 0.5 * force_(i, j, 1) * dt_) / density;
        }
    }
}

void TwoPopulationSolver::third_moment_error() {
    const bool parallel = parallel_;
    // 1 - 3 (c_s^k)^2 in lattice units: zero for a fluid at the lattice sound
    // speed, and 1.8 alpha_k - 0.8 in Ba et al.'s notation.
    const double weight1 = 1.0 - cs_squared_[0] / cs2_;
    const double weight2 = 1.0 - cs_squared_[1] / cs2_;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            const double mass = weight1 * rho_1_(i, j) + weight2 * rho_2_(i, j);
            q_x_(i, j) = mass * u_(i, j, 0);
            q_y_(i, j) = mass * u_(i, j, 1);
        }
    }
}

void TwoPopulationSolver::collide_node(int i, int j, double* out1, double* out2) const {
    const double u_x = u_(i, j, 0);
    const double u_y = u_(i, j, 1);
    const double density = rho_(i, j);
    const double fraction = 0.5 * (1.0 + phi_n_(i, j));  // volume fraction of fluid 1
    // Ba et al. Eq. (19): mu = dt (1/s_nu - 1/2) p, with mu interpolated
    // across the interface (their Eq. 22 blends the rate parabolically;
    // the volume fraction is used here, as in `Solver`).
    const double mu = fraction * mu_[0] + (1.0 - fraction) * mu_[1];
    const double tau = mu / (p_(i, j) * dt_) + 0.5;
    const double omega = 1.0 / tau;
    const double guo = (1.0 - 0.5 * omega) * dt_;

    const double f_x = force_(i, j, 0);
    const double f_y = force_(i, j, 1);

    double eq1[kQ], eq2[kQ];
    equilibrium(0, rho_1_(i, j), u_x, u_y, eq1);
    equilibrium(1, rho_2_(i, j), u_x, u_y, eq2);

    // Guo's forcing. It is split between the fluids in proportion to their
    // share of the mass, so the total momentum comes out right.
    double source[kQ];
    for (int k = 0; k < kQ; k++) {
        const double xi_x = kXi[k][0], xi_y = kXi[k][1];
        const double eu = xi_x * u_x + xi_y * u_y;
        source[k] =
            kW[k] *
            (((xi_x - u_x) + eu * xi_x / cs2_) * f_x + ((xi_y - u_y) + eu * xi_y / cs2_) * f_y) /
            cs2_;
    }
    const double share1 = rho_1_(i, j) / density;
    const double s_e = mrt_ ? config_.physics.s_e : omega;

    // Ba et al. Eqs. (17)-(18): the divergence of Q added back into the trace
    // of the stress and its normal difference. Their Eq. (20) is the
    // nine-point isotropic derivative, which is E4 here.
    double correction[kQ] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    if (correction_) {
        double dqx_dx = 0.0, dqx_dy = 0.0, dqy_dx = 0.0, dqy_dy = 0.0;
        if (wall_y_) {
            gradient_wall_y(q_x_.data(), nx_, ny_, i, j, GradientStencil::E4, &dqx_dx, &dqx_dy);
            gradient_wall_y(q_y_.data(), nx_, ny_, i, j, GradientStencil::E4, &dqy_dx, &dqy_dy);
        } else {
            gradient_periodic(q_x_.data(), nx_, ny_, i, j, GradientStencil::E4, &dqx_dx, &dqx_dy);
            gradient_periodic(q_y_.data(), nx_, ny_, i, j, GradientStencil::E4, &dqy_dx, &dqy_dy);
        }
        add_third_moment_source(dqx_dx, dqy_dy, s_e, omega, dt_, correction);
    }

    if (!mrt_) {
        for (int k = 0; k < kQ; k++) {
            const double f1 = f1_(i, j, k);
            const double f2 = f2_(i, j, k);
            out1[k] = f1 + (-omega * (f1 - eq1[k]) + guo * share1 * source[k]);
            out2[k] = f2 + (-omega * (f2 - eq2[k]) + guo * (1.0 - share1) * source[k]);
            if (correction_) {
                out1[k] += share1 * correction[k];
                out2[k] += (1.0 - share1) * correction[k];
            }
        }
        return;
    }

    // MRT on the total, which is all the recolouring reads.
    double rates[kQ];
    mrt_rates(omega, s_e, config_.physics.s_eps, config_.physics.s_q, rates);
    double total[kQ], total_eq[kQ], post[kQ];
    for (int k = 0; k < kQ; k++) {
        total[k] = f1_(i, j, k) + f2_(i, j, k);
        total_eq[k] = eq1[k] + eq2[k];
    }
    mrt_collide(total, total_eq, source, rates, dt_, post);
    for (int k = 0; k < kQ; k++) {
        post[k] += correction[k];
        out1[k] = share1 * post[k];
        out2[k] = post[k] - out1[k];
    }
}

void TwoPopulationSolver::collide() {
    if (correction_) {
        third_moment_error();
    }
    const bool parallel = parallel_;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            double out1[kQ], out2[kQ];
            collide_node(i, j, out1, out2);
            for (int k = 0; k < kQ; k++) {
                f1_(i, j, k) = out1[k];
                f2_(i, j, k) = out2[k];
            }
        }
    }
}

void TwoPopulationSolver::collide_node_for_test(int i, int j, double* out) {
    if (correction_) {
        third_moment_error();
    }
    double out1[kQ], out2[kQ];
    collide_node(i, j, out1, out2);
    for (int k = 0; k < kQ; k++) {
        out[k] = out1[k] + out2[k];
    }
}

double TwoPopulationSolver::rest_weight(int i, int j, int k) const {
    const double e_squared = speed_squared(k);
    const double phi_1 = k == 0 ? phi_rest_[0] : (e_squared < 1.5 ? phi_near_[0] : phi_diag_[0]);
    const double phi_2 = k == 0 ? phi_rest_[1] : (e_squared < 1.5 ? phi_near_[1] : phi_diag_[1]);
    return (rho_1_(i, j) * phi_1 + rho_2_(i, j) * phi_2) / rho_(i, j);
}

void TwoPopulationSolver::recolor() {
    const bool parallel = parallel_;
    const double beta = config_.physics.beta;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            const double density = rho_(i, j);
            const double share1 = rho_1_(i, j) / density;
            const double gx = grad_x_(i, j);
            const double gy = grad_y_(i, j);
            const double norm = std::sqrt(gx * gx + gy * gy);
            // Ba et al. Eq. (30), the Latva-Kokko segregation operator. It is
            // usable here, unlike in `Solver`, because the two populations
            // stream separately and their positivity is what bounds the phase
            // field -- there is no `g = f phi + omega_3` to overshoot.
            const double strength = norm > kGradientEpsilon
                                        ? beta * rho_1_(i, j) * rho_2_(i, j) / (density * norm)
                                        : 0.0;
            for (int k = 0; k < kQ; k++) {
                const double total = f1_(i, j, k) + f2_(i, j, k);
                double push = 0.0;
                if (k != 0 && strength != 0.0) {
                    // cos(phi_i) = (e_i . grad phi_N) / (|e_i| |grad phi_N|)
                    const double speed = std::sqrt(speed_squared(k));
                    push =
                        rest_weight(i, j, k) * strength * (kXi[k][0] * gx + kXi[k][1] * gy) / speed;
                }
                f1_(i, j, k) = share1 * total + push;
                f2_(i, j, k) = total - f1_(i, j, k);
            }
        }
    }
}

void TwoPopulationSolver::stream() {
    Field next1(nx_, ny_, kQ);
    Field next2(nx_, ny_, kQ);
    const bool parallel = parallel_;
    const bool wall = wall_y_;
#pragma omp parallel for collapse(2) if (parallel)
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            for (int k = 0; k < kQ; k++) {
                const int ip = (i + static_cast<int>(kXi[k][0]) + nx_) % nx_;
                int jp, kp;
                if (wall && j == 0 && (k == 4 || k == 7 || k == 8)) {
                    jp = j;
                    kp = k - 2;  // half-way bounce-back, as in `Solver::stream`
                } else if (wall && j == ny_ - 1 && (k == 2 || k == 5 || k == 6)) {
                    jp = j;
                    kp = k + 2;
                } else {
                    jp = wall ? j + static_cast<int>(kXi[k][1])
                              : (j + static_cast<int>(kXi[k][1]) + ny_) % ny_;
                    kp = k;
                }
                next1(ip, jp, kp) = f1_(i, j, k);
                next2(ip, jp, kp) = f2_(i, j, k);
            }
        }
    }
    f1_ = std::move(next1);
    f2_ = std::move(next2);
}

void TwoPopulationSolver::initialize() {
    const Physics& physics = config_.physics;
    for (int i = 0; i < nx_; i++) {
        for (int j = 0; j < ny_; j++) {
            // `initial_phase` gives phi_N = 2c - 1, so the case's prescribed
            // radius is the radius of the region actually full of fluid 1.
            double indicator = config_.initial_phase(config_, i, j);
            if (indicator > 1.0) {
                indicator = 1.0;
            } else if (indicator < -1.0) {
                indicator = -1.0;
            }
            const double fraction = 0.5 * (1.0 + indicator);
            const double rho_1 = fraction * physics.rho1;
            const double rho_2 = (1.0 - fraction) * physics.rho2;
            double u_x = 0.0, u_y = 0.0;
            if (config_.initial_velocity) {
                config_.initial_velocity(config_, i, j, &u_x, &u_y);
            }

            double eq1[kQ], eq2[kQ];
            equilibrium(0, rho_1, u_x, u_y, eq1);
            equilibrium(1, rho_2, u_x, u_y, eq2);
            for (int k = 0; k < kQ; k++) {
                f1_(i, j, k) = eq1[k];
                f2_(i, j, k) = eq2[k];
            }
        }
    }
    densities();
    update_colour_gradient();
    surface_force();
    update_velocity();
}

void TwoPopulationSolver::refresh() {
    densities();
    update_velocity();
}

void TwoPopulationSolver::step() {
    densities();
    update_colour_gradient();
    surface_force();
    update_velocity();
    collide();
    recolor();
    stream();
}

void TwoPopulationSolver::run() {
    CsvWriter writer(config_.output_precision);
    initialize();
    writer.write_grids(0, state());
    for (int timestep = 1; timestep <= config_.steps; timestep++) {
        step();
        if (timestep % config_.interval == 0) {
            std::cout << "Step " << timestep << std::endl;
            refresh();
            writer.write_grids(timestep, state());
        }
    }
}

}  // namespace lbm
}  // namespace cglbm
