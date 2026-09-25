#ifndef CGLBM_LBM_TWO_POPULATION_SOLVER_H
#define CGLBM_LBM_TWO_POPULATION_SOLVER_H

#include "lbm/case_config.h"
#include "lbm/d2q9.h"
#include "lbm/field.h"
#include "lbm/output_writer.h"

/// The classical colour-gradient scheme, carrying one distribution per fluid.
///
/// `Solver` and this one are both colour-gradient models and they differ in
/// where the density ratio lives, which turns out to decide how far it can go.
///
///  - `Solver`, after Lafarge et al., carries the mixture `f` and the colour
///    difference `g`, and gets the density ratio from a two-component equation
///    of state `p = p(rho, phi)`. The whole contrast sits in one density field,
///    and the pressure has to stay smooth across a jump of that size. Measured,
///    that scheme is steady to a ratio of about 500 and diverges at 1000.
///  - This one, after Grunau, Reis & Phillips, Leclaire et al. and Ba et al.,
///    carries `f^1` and `f^2` and streams them separately. The density ratio
///    comes from a free parameter in the *equilibrium*: each fluid keeps its own
///    rest-particle weight `alpha_k`, hence its own sound speed
///
///        (c_s^k)^2 = 3/5 (1 - alpha_k),   p_k = rho_k (c_s^k)^2
///
///    and the ratio is fixed by `rho_1/rho_2 = (1 - alpha_2)/(1 - alpha_1)`.
///    Two things follow that the other model cannot have. The bulk pressures
///    match identically, `rho_1 (c_s^1)^2 = rho_2 (c_s^2)^2`, so the pressure is
///    continuous across the interface however large the ratio; and each
///    `rho_k` has its own smooth profile, so nothing has to resolve a jump of
///    `rho_1/rho_2` in a single field.
///
/// The cost is the classical limitation the other model was written to escape:
/// the density ratio and the sound-speed ratio are no longer independent. At a
/// ratio of 10^4 the heavy fluid's sound speed is 0.0069 in lattice units.
/// For a static or slow flow that is a fair trade; for anything acoustic it is
/// not, and `Solver` remains the right choice below a ratio of a few hundred.
///
/// One step is
///
///     densities  ->  colour gradient  ->  surface_force  ->  velocity
///                ->  collide  ->  recolor  ->  stream
///
/// split that way because the velocity carries half the capillary force
/// (Ba et al. Eq. 29) while the force itself is built from the densities.
///
/// References
///  - D. Grunau, S. Chen, K. Eggert, "A lattice Boltzmann model for multiphase
///    fluid flows", Phys. Fluids A 5, 2557 (1993). The alpha_k equilibrium.
///  - Y. Ba, H. Liu, Q. Li, Q. Kang, J. Sun, Phys. Rev. E 94, 023310 (2016).
///    The form implemented here: Eq. (6) for phi_i^k, Eq. (7) for the density
///    ratio, Eq. (14) for the enhanced equilibrium, Eq. (21) for the phase
///    field, Eqs. (23)-(29) for the continuum-surface-force tension, Eq. (30)
///    for the Latva-Kokko recolouring.
///  - S. Leclaire, M. Reggio, J.-Y. Trepanier, Computers & Fluids 48, 98
///    (2011). Reaches O(10^4) on a static bubble with this model family and an
///    isotropic colour gradient.

namespace cglbm {
namespace lbm {

class TwoPopulationSolver {
public:
    /// Build the lattice described by `config` and allocate its fields.
    ///
    /// Uses `physics.rho1`, `rho2`, `nu`, `nu2`, `sigma`, `radius`, `beta`,
    /// `alpha2` and `gravity`, plus `boundary`, `stencil` and `initial_phase`.
    /// The equation-of-state fields -- `c1`, `c2`, `p1_inf`, `p2_inf` -- and
    /// `initial_state` play no part here: this model has no equation of state
    /// to invert, and its sound speeds come from `alpha2` and the density
    /// ratio. `initial_phase` is read as the *volume fraction* indicator,
    /// `2c - 1`, which is what `phi_N` means in both models.
    ///
    /// Throws `std::invalid_argument` when the configuration cannot be run.
    explicit TwoPopulationSolver(CaseConfig config);

    void initialize();
    void step();
    void run();

    MacroscopicState state() const;

    /// Recompute rho, p, phi_N and u from the distributions.
    ///
    /// `step()` leaves the populations streamed but the macroscopic fields as
    /// they were before the collision, so anything reading `density()` or
    /// `pressure()` between steps wants this first. `run()` calls it before
    /// writing.
    void refresh();

    const CaseConfig& config() const {
        return config_;
    }
    const Field& density() const {
        return rho_;
    }
    const Field& velocity() const {
        return u_;
    }
    /// The phase field, as the bulk-normalised indicator `2c - 1`.
    const Field& phase() const {
        return phi_n_;
    }
    const Field& pressure() const {
        return p_;
    }

    /// One fluid's distribution, `fluid` being 0 or 1.
    ///
    /// Exposed because non-negativity of both is the invariant this model
    /// stands on: it is what bounds the phase field, and it is what the
    /// recolouring has to respect. See `rest_weight`.
    const Field& population(int fluid) const {
        return fluid == 0 ? f1_ : f2_;
    }

    /// Density of one fluid, `fluid` being 0 or 1.
    const Field& component_density(int fluid) const {
        return fluid == 0 ? rho_1_ : rho_2_;
    }

    /// The equilibrium of one fluid, for the unit tests to take moments of.
    void equilibrium_for_test(int fluid, double rho_k, double u_x, double u_y, double* out) const {
        equilibrium(fluid, rho_k, u_x, u_y, out);
    }

    /// Rest-particle weights of the two fluids, as resolved from the densities.
    double alpha1() const {
        return alpha1_;
    }
    double alpha2() const {
        return alpha2_;
    }

private:
    /// rho_k, rho, p and phi_N from the distributions.
    void densities();
    /// u from the distributions and half the body force, Ba et al. Eq. (29).
    void update_velocity();
    void update_colour_gradient();
    void surface_force();
    void collide();
    void recolor();
    void stream();

    /// The mixture's rest weight in direction `k`, `f_i^eq(rho, 0) / rho`.
    ///
    /// This is what the recolouring pushes along, in place of the lattice
    /// weight `w_i` that Latva-Kokko and Rothman use. The two agree exactly
    /// when the fluids share a rest weight, and they must differ when the
    /// fluids do not -- which is the whole of how this model carries a density
    /// ratio. `w_i` there would push a fixed fraction of the *lattice* weight,
    /// while what has to stay positive is a fraction of the *population*, and
    /// in the heavy fluid the non-rest populations are `(1 - alpha_1)/5` of the
    /// density: 1.6e-4 at a ratio of 1000, against `w_i = 1/9`. Latva-Kokko's
    /// form drives `f^2` negative there by a factor of 240, and the run is gone
    /// within ten steps.
    ///
    /// Writing `Phi_i = (rho_1 phi_i^1 + rho_2 phi_i^2) / rho` instead keeps
    /// both requirements the operator has to meet:
    ///
    ///  - each fluid's mass is conserved exactly, because `Phi_i` depends only
    ///    on `|e_i|` and the cosine is odd, so the push sums to zero;
    ///  - both populations stay non-negative for any `beta <= 1`, because near
    ///    equilibrium the push is `beta (rho_1/rho)` of fluid 2's own share.
    ///
    /// Leclaire, Reggio & Trepanier (2012) report adapting Latva-Kokko's
    /// operator "for the case of variable density ratios", and Leclaire et al.
    /// (2011) credit that adaptation with reaching a density ratio of 10^4.
    /// Their text was not available here, so this form is derived from the two
    /// requirements above rather than transcribed; what can be checked is that
    /// it reduces to Latva-Kokko exactly at equal rest weights, and what it
    /// measures is in docs/numerics.md.
    double rest_weight(int i, int j, int k) const;

    /// Equilibrium of fluid `k` (0 or 1) at a node, into `out[kQ]`.
    void equilibrium(int fluid, double rho_k, double u_x, double u_y, double* out) const;

    void gradient_at(const double* field, int i, int j, double* grad_x, double* grad_y) const;

    CaseConfig config_;
    int nx_;
    int ny_;
    double dt_;
    double cs2_;
    bool wall_y_;
    bool parallel_;

    double alpha1_;         ///< rest weight of fluid 1, from the density ratio
    double alpha2_;         ///< rest weight of fluid 2, a free parameter
    double cs_squared_[2];  ///< (c_s^k)^2 = 3/5 (1 - alpha_k)
    double phi_rest_[2];    ///< phi_i^k for the rest direction: alpha_k
    double phi_near_[2];    ///< for the four axial directions: (1 - alpha_k)/5
    double phi_diag_[2];    ///< for the four diagonals: (1 - alpha_k)/20
    double mu_[2];          ///< bulk dynamic viscosity rho_k0 nu_k of each fluid

    Field rho_1_;     ///< density of fluid 1
    Field rho_2_;     ///< density of fluid 2
    Field rho_;       ///< mixture density, rho_1 + rho_2
    Field p_;         ///< mixture pressure, sum of rho_k (c_s^k)^2
    Field u_;         ///< velocity, two components per node
    Field phi_n_;     ///< bulk-normalised phase field, Ba et al. Eq. (21)
    Field force_;     ///< body force: capillary plus gravity
    Field grad_x_;    ///< colour gradient, x
    Field grad_y_;    ///< colour gradient, y
    Field normal_x_;  ///< interface normal, x
    Field normal_y_;  ///< interface normal, y

    Field f1_;  ///< distribution of fluid 1
    Field f2_;  ///< distribution of fluid 2
};

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_TWO_POPULATION_SOLVER_H
