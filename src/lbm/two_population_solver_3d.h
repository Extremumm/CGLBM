#ifndef CGLBM_LBM_TWO_POPULATION_SOLVER_3D_H
#define CGLBM_LBM_TWO_POPULATION_SOLVER_3D_H

#include "lbm/case_config.h"
#include "lbm/d3q19.h"
#include "lbm/field3d.h"
#include "lbm/isotropic_gradient_3d.h"

/// The two-population colour-gradient scheme on D3Q19.
///
/// The three-dimensional counterpart of `TwoPopulationSolver`, and the same
/// model: one distribution per fluid, the density ratio carried by each fluid's
/// rest-particle weight rather than by an equation of state. What changes with
/// the extra dimension is small and worth stating, because it is where a
/// careless port would go wrong.
///
///  - The rest weights. Requiring `sum_i phi_i^k = 1` and a fourth-order
///    isotropic `sum_i phi_i^k e e e e` fixes them, on D3Q19, at
///
///        phi_0 = alpha_k,  phi_axial = (1 - alpha_k)/12,
///        phi_edge = (1 - alpha_k)/24,  (c_s^k)^2 = (1 - alpha_k)/2,
///
///    where the same derivation on D2Q9 gives `(1-alpha)/5`, `(1-alpha)/20` and
///    `3(1-alpha)/5`. The relation between the rest weights and the density
///    ratio, `rho_1/rho_2 = (1 - alpha_2)/(1 - alpha_1)`, is unchanged, and so
///    is its consequence that the two bulk pressures are identically equal.
///
///  - The enhanced equilibrium. Its correction runs over `3|e_i|^2 - (d + 2)`,
///    which is `3|e_i|^2 - 4` in two dimensions and `3|e_i|^2 - 5` here. That
///    follows from the third-order Hermite contraction
///    `H_xxx + H_yyx + H_zzx = e_x(|e|^2 - 5 c_s^2)`, and it is what keeps the
///    correction from disturbing the momentum: `sum_i w_i e e (3|e|^2 - 5) = 0`
///    on D3Q19 exactly as `sum_i w_i e e (3|e|^2 - 4) = 0` on D2Q9.
///
///  - The curvature. `K = -(div n - n_a n_b d_a n_b)`, the surface divergence,
///    which is `-div n` for a unit normal and stays bounded when the discrete
///    normal is not quite one. On a sphere of radius R it is `-2/R`, so
///    Laplace's law reads `2 sigma / R` rather than `sigma / R`.
///
/// Everything else -- the continuum-surface-force tension, the segregation
/// operator and its adaptation to a density ratio, the matched dynamic
/// viscosities -- carries over unchanged. See `two_population_solver.h` for the
/// derivations and docs/report for the proofs.

namespace cglbm {
namespace lbm {

/// The macroscopic fields of a three-dimensional run.
struct MacroscopicState3D {
    const Field3D* density;
    const Field3D* velocity;
    const Field3D* phase;
    const Field3D* pressure;
};

class TwoPopulationSolver3D {
public:
    /// Build the lattice described by `config` and allocate its fields.
    ///
    /// Reads `nx`, `ny`, `nz`, `stencil_3d` and `initial_phase_3d` alongside
    /// the physics; the two-dimensional `stencil` and `initial_phase` are
    /// ignored. Throws `std::invalid_argument` when the case cannot be run.
    explicit TwoPopulationSolver3D(CaseConfig config);

    void initialize();
    void step();
    void run();

    /// Recompute the macroscopic fields from the distributions.
    void refresh();

    MacroscopicState3D state() const;

    const CaseConfig& config() const {
        return config_;
    }
    const Field3D& density() const {
        return rho_;
    }
    const Field3D& velocity() const {
        return u_;
    }
    /// The phase field, as the bulk-normalised indicator `2c - 1`.
    const Field3D& phase() const {
        return phi_n_;
    }
    const Field3D& pressure() const {
        return p_;
    }
    const Field3D& population(int fluid) const {
        return fluid == 0 ? f1_ : f2_;
    }
    const Field3D& component_density(int fluid) const {
        return fluid == 0 ? rho_1_ : rho_2_;
    }

    double alpha1() const {
        return alpha1_;
    }
    double alpha2() const {
        return alpha2_;
    }

    /// The equilibrium of one fluid, for the unit tests to take moments of.
    void equilibrium_for_test(int fluid, double rho_k, const double* u, double* out) const {
        equilibrium(fluid, rho_k, u[0], u[1], u[2], out);
    }

private:
    void densities();
    void update_colour_gradient();
    void surface_force();
    void update_velocity();
    void collide();
    void recolor();
    void stream();

    /// Write the z = nz/2 slice through `CsvWriter`, in the same four files a
    /// two-dimensional run produces, so the existing post-processing reads a
    /// three-dimensional run unchanged.
    void write_midplane(class CsvWriter& writer, int timestep) const;

    void
    equilibrium(int fluid, double rho_k, double u_x, double u_y, double u_z, double* out) const;

    /// The mixture's rest weight in direction `q`, `f_q^eq(rho, 0) / rho`.
    ///
    /// What the segregation operator pushes along, in place of the lattice
    /// weight. See `TwoPopulationSolver::rest_weight`; the argument is
    /// identical here and the numbers only change through `phi_q^k`.
    double rest_weight(int i, int j, int k, int q) const;

    void gradient_at(const double* field,
                     int i,
                     int j,
                     int k,
                     double* grad_x,
                     double* grad_y,
                     double* grad_z) const;

    CaseConfig config_;
    int nx_;
    int ny_;
    int nz_;
    double dt_;
    double cs2_;
    bool wall_y_;
    bool parallel_;

    double alpha1_;
    double alpha2_;
    double cs_squared_[2];  ///< (c_s^k)^2 = (1 - alpha_k)/2 on D3Q19
    double phi_rest_[2];    ///< alpha_k
    double phi_axial_[2];   ///< (1 - alpha_k)/12
    double phi_edge_[2];    ///< (1 - alpha_k)/24
    double mu_[2];          ///< bulk dynamic viscosity of each fluid

    Field3D rho_1_;
    Field3D rho_2_;
    Field3D rho_;
    Field3D p_;
    Field3D u_;
    Field3D phi_n_;
    Field3D force_;
    Field3D grad_x_;
    Field3D grad_y_;
    Field3D grad_z_;
    Field3D normal_x_;
    Field3D normal_y_;
    Field3D normal_z_;

    Field3D f1_;
    Field3D f2_;
};

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_TWO_POPULATION_SOLVER_3D_H
