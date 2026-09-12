#ifndef CGLBM_LBM_TWO_POPULATION_SOLVER_3D_H
#define CGLBM_LBM_TWO_POPULATION_SOLVER_3D_H

#include <memory>
#include <vector>

#include "lbm/case_config.h"
#include "lbm/field3d.h"
#include "lbm/isotropic_gradient_3d.h"
#include "lbm/lattice3d.h"
#include "lbm/quasi_static_mhd_3d.h"

/// The two-population colour-gradient scheme on D3Q19 or D3Q27.
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
///    D3Q27 has a third moving shell and needs a third condition; `lattice3d.h`
///    carries both solutions and the solver reads them from there rather than
///    branching on the lattice.
///
///  - The enhanced equilibrium. *Two* things change here, and the second one is
///    easy to miss. The correction runs over `3|e_i|^2 - (d + 2)`, which is
///    `3|e_i|^2 - 4` in two dimensions and `3|e_i|^2 - 5` here, from the
///    third-order Hermite contraction
///    `H_xxx + H_yyx + H_zzx = e_x(|e|^2 - 5 c_s^2)`. But its *amplitude*
///    changes too:
///
///        lambda = ((c_s^k)^2 - c_s^2) / (3 T),
///        T = sum_i w_i e_x^2 e_y^2 (3|e_i|^2 - (d + 2)),
///
///    and `T` is `2/9` on D2Q9 against `1/9` on D3Q19 -- the (+-1, +-1, 0)
///    shell carries `3|e|^2 - 5 = 1` where in two dimensions it carried
///    `3|e|^2 - 4 = 2`. So `lambda` is `(3 (c_s^k)^2 - 1) / 2` there and
///    `3 (c_s^k)^2 - 1` here, twice as large.
///
///    What that amplitude buys is the shear stress. The correction exists to
///    drag the equilibrium's third moment from the lattice's `c_s^2 = 1/3`
///    towards the fluid's own `(c_s^k)^2`, and for `a != b` the viscous stress
///    depends on the third moment only through `M3_aab`, so `lambda` is fixed
///    by `M3_xxy = (c_s^k)^2 u_y`. Carry the two-dimensional half over unchanged
///    and the correction closes only half the gap, leaving an effective
///    `((c_s^k)^2 + 1/3) / 2` -- the *average* of the lattice sound speed and
///    the fluid's. Since `collide()` sets tau from `(c_s^k)^2`, that is a shear
///    viscosity too large by `1 / (6 (c_s^k)^2) + 1/2`: a factor of 4.7 in the
///    droplet at a density ratio of 10, and 417 at a ratio of 1000. It is
///    invisible to the mass, momentum and second-moment checks, because
///    `lambda` does not enter any of them; `report_equilibrium` in the unit
///    test measures `M3_xxy` for exactly that reason.
///
///    What does *not* change with the dimension is that the correction leaves
///    the momentum alone: `sum_i w_i e e (3|e|^2 - 5) = 0` on D3Q19 exactly as
///    `sum_i w_i e e (3|e|^2 - 4) = 0` on D2Q9. That is the necessary condition,
///    not the sufficient one, and checking it alone is what let the amplitude
///    stay wrong.
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
///
/// Nothing above is written down twice for the two lattices. The rest weights
/// and `(c_s^k)^2` come from the lattice descriptor, and the amplitude
/// `lambda` is *summed over the lattice* at construction rather than tabulated,
/// precisely because it is the quantity a port gets wrong and no conservation
/// check notices. See `lattice3d.h`.
///
/// # The magnetic coupling
///
/// A case that sets `CaseConfig::mhd` gains an inductionless magnetic force.
/// The potential is solved and the current rebuilt once per step, and the
/// Lorentz force joins the surface tension, gravity and the body force in the
/// same field, so it reaches the populations through the same Guo forcing.
///
/// Where it sits in the step matters. The force needs a velocity, and the
/// velocity the solver reports carries half a step of the force already --
/// using it would make the magnetic force depend on itself. So with a field
/// present the velocity is formed twice: once as the bare momentum
/// `sum_q f_q e_q / rho`, which is what the field sees, and then again with the
/// half-step correction once the force is known. Without a field the extra pass
/// is skipped and the time loop is what it always was, to the last bit.
///
/// See `quasi_static_mhd_3d.h` for the model, and note the step bound it
/// carries: magnetic damping alone relaxes the velocity on
/// `rho / (sigma B^2)`.

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

    /// The lattice this run is on.
    const Lattice3D& lattice() const {
        return *lattice_;
    }

    /// `(c_s^k)^2` of fluid `k`, as the lattice and `alpha_k` give it.
    double sound_speed_squared(int fluid) const {
        return cs_squared_[fluid];
    }

    /// The enhanced equilibrium's amplitude for fluid `k`.
    ///
    /// Derived from the lattice at construction, not tabulated. A test that
    /// wants to know whether the third moment comes out right should measure
    /// the moment; this is here so it can also report what produced it.
    double enhancement(int fluid) const {
        return enhancement_[fluid];
    }

    /// The magnetic module, or null when the case has no field.
    const QuasiStaticMhd3D* mhd() const {
        return mhd_.get();
    }

    /// Distance from the domain centre to the interface along +x, +y and +z.
    ///
    /// The `phi_N = 0` crossing, linearly interpolated between the two nodes
    /// that bracket it -- the same estimate `pycglbm` takes of the mid-plane
    /// slice, in three directions instead of one. `nan` for an axis along which
    /// the phase field never changes sign.
    ///
    /// For a droplet these are its three semi-axes, and the mode-2 oscillation
    /// `oscillation_3d` measures is the difference between the polar one and
    /// the equatorial ones. `run()` appends them to `interface.csv` every step
    /// when `config.track_interface` is set.
    void interface_axes(double* radii) const;

    /// The equilibrium of one fluid, for the unit tests to take moments of.
    void equilibrium_for_test(int fluid, double rho_k, const double* u, double* out) const {
        equilibrium(fluid, rho_k, u[0], u[1], u[2], out);
    }

private:
    void densities();
    void update_colour_gradient();
    void surface_force();

    /// The velocity, with or without the half step of the force.
    ///
    /// `with_force` false is the bare momentum `sum_q f_q e_q / rho`, which is
    /// the velocity the magnetic force is evaluated at; true adds Guo's half
    /// step and is what the collision and the output use.
    void update_velocity(bool with_force);
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
    const Lattice3D* lattice_;
    int nx_;
    int ny_;
    int nz_;
    double dt_;
    double cs2_;
    bool wall_y_;
    bool parallel_;
    bool has_body_force_;

    double alpha1_;
    double alpha2_;
    double cs_squared_[2];    ///< (c_s^k)^2, from the lattice and alpha_k
    double enhancement_[2];   ///< lambda, summed over the lattice at construction
    double mu_[2];            ///< bulk dynamic viscosity of each fluid
    std::vector<double> rest_weights_[2];  ///< phi_q^k, one entry per velocity

    std::unique_ptr<QuasiStaticMhd3D> mhd_;

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
