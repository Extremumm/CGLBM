#ifndef CGLBM_LBM_SOLVER_H
#define CGLBM_LBM_SOLVER_H

#include "lbm/case_config.h"
#include "lbm/d2q9.h"
#include "lbm/equation_of_state.h"
#include "lbm/field.h"
#include "lbm/output_writer.h"

/// The improved colour-gradient scheme on D2Q9.
///
/// One time step is
///
///     force  ->  collide  ->  collide_surface  ->  recolor  ->  stream
///            ->  macroscopic  ->  phase_field  ->  surface_force
///            ->  equilibrium
///
/// carrying two distribution functions: `f`, the sum of the two components'
/// populations, and `g`, their difference. `f` transports the mixture; `phi =
/// sum(g) / sum(f)` is the phase field that says which component a node holds.
///
/// The collision is split into three operators applied to the same populations:
///
///  - `omega_1`, the regularised viscous relaxation, on the shear and bulk
///    Hermite moments with separate relaxation times `tau_nu` and `tau_b`;
///  - `omega_2`, the surface-tension operator, which acts along the colour
///    gradient and vanishes wherever the gradient does. A case may instead ask
///    for the continuum-surface-force form of Ba et al., in which case
///    `omega_2` stays zero and `surface_force` builds a body force from an
///    explicit curvature; see `SurfaceTension`;
///  - `omega_3`, the recolouring operator, which pushes the two components
///    apart along that same gradient and holds the interface at `ch_width_ope`.
///
/// The source term `S` gathers the body force (Guo's forcing), a correction for
/// the third-order moment that D2Q9 cannot resolve, and a temporal correction
/// from the change in `p - rho cs^2` since the previous step.
///
/// References
///  - T. Lafarge, P. Boivin, N. Odier, B. Cuenot, "Improved color-gradient
///    method for lattice Boltzmann modeling of two-phase flows", Physics of
///    Fluids 33(8), 082110 (2021), doi:10.1063/5.0061638. The scheme
///    implemented here, including the two-component equation of state and the
///    corrective source terms.
///  - Y. Ba, H. Liu, Q. Li, Q. Kang, J. Sun, "Multiple-relaxation-time
///    color-gradient lattice Boltzmann model for simulating two-phase flows
///    with high density ratio", Phys. Rev. E 94, 023310 (2016),
///    doi:10.1103/PhysRevE.94.023310. The normalised phase field, Eq. (21),
///    and the continuum-surface-force tension, Eqs. (23)-(29).
///  - S. Leclaire, N. Pellerin, M. Reggio, J.-Y. Trepanier, "Enhanced
///    equilibrium distribution functions for simulating immiscible multiphase
///    flows with variable density ratios in a class of lattice Boltzmann
///    models", Int. J. Multiphase Flow 57, 159 (2013). The third-order
///    correction to the equilibrium, which the Hermite equilibrium used here
///    already contains term for term -- see `report_enhanced_equilibrium` in
///    programs/unit_testing/lbm/solver.
///  - T. Kruger et al., "The Lattice Boltzmann Method: Principles and
///    Practice", Springer (2017). Periodic boundaries p. 170, half-way
///    bounce-back p. 176, viscosity in lattice units p. 284.

namespace cglbm {
namespace lbm {

class Solver {
public:
    /// Build the lattice described by `config` and allocate its fields.
    ///
    /// Throws `std::invalid_argument` when the configuration cannot be run --
    /// a non-positive lattice, or no initial phase field.
    explicit Solver(CaseConfig config);

    /// Lay down the initial state: phase field, density, pressure, and the
    /// equilibrium populations they imply.
    void initialize();

    /// Advance one time step.
    void step();

    /// Initialise, then advance `config.steps` steps, writing output as the
    /// configuration asks.
    ///
    /// Progress goes to stdout; the CSV grids go to the working directory.
    void run();

    /// The macroscopic fields, for a writer or a test.
    MacroscopicState state() const;

    const CaseConfig& config() const {
        return config_;
    }

    /// Read-only views of the state, for tests that check an invariant.
    const Field& density() const {
        return rho_;
    }
    const Field& velocity() const {
        return u_;
    }
    const Field& phase() const {
        return phi_;
    }
    const Field& pressure() const {
        return p_;
    }

private:
    // The eight stages of a step, in the order `step()` applies them.
    void force();
    void collide();
    void collide_surface();
    void recolor();
    void stream();
    void macroscopic();
    /// Capillary body force from an explicit curvature, Ba et al. Eq. (24).
    ///
    /// Only runs under `SurfaceTension::ContinuumSurfaceForce`; otherwise
    /// `force_surface_` stays zero and the tension comes from `omega_2_`.
    void surface_force();
    void phase_field();
    void equilibrium();

    /// Colour gradient at (i, j), following the case's boundary along y.
    ///
    /// Taken of whichever field `config.interface_field` selects, which is what
    /// centres the surface-tension and recolouring operators on the interface
    /// rather than on the zero of the raw colour field.
    void colour_gradient(int i, int j, double* grad_x, double* grad_y) const;

    /// Gradient of any flat `nx * ny` lattice field, following the same
    /// boundary rule. `colour_gradient` is this applied to the phase field;
    /// `surface_force` also differentiates the interface normal with it.
    void gradient_at(const double* field, int i, int j, double* grad_x, double* grad_y) const;

    /// Refresh `phi_n_` from `phi_`, when the case asks for a normalised field.
    void update_interface_field();

    CaseConfig config_;
    int nx_;
    int ny_;

    // Lattice constants, resolved once from config_.units.
    double dx_;
    double dt_;
    double cs2_;
    double cs4_;
    double cs6_;

    bool wall_y_;               ///< config_.boundary == Boundary::WallY
    bool normalise_interface_;  ///< config_.interface_field == BulkNormalised
    bool parallel_;             ///< config_.parallel, copied out for the OpenMP if clause
    ComponentPair components_;

    // Macroscopic fields.
    Field rho_;            ///< density
    Field rho_mdt_;        ///< density at the previous step, for the temporal correction
    Field u_;              ///< velocity, two components per node
    Field p_;              ///< pressure
    Field p_mdt_;          ///< pressure at the previous step
    Field phi_;            ///< phase field
    Field phi_n_;          ///< phase field normalised by the bulk densities, Ba Eq. (21)
    Field force_;          ///< external volume force, two components per node
    Field force_surface_;  ///< capillary part of it, Ba et al. Eq. (24)
    Field normal_x_;       ///< interface normal, x, kept flat for the gradient stencil
    Field normal_y_;       ///< interface normal, y
    Field gradient_norm_;  ///< |grad phi_N|, the interface delta the force rides on

    // Populations and collision operators.
    Field f_;        ///< sum of the two components' distributions
    Field g_;        ///< difference of the two components' distributions
    Field f_eq_;     ///< equilibrium of f
    Field omega_1_;  ///< viscous relaxation
    Field omega_2_;  ///< surface tension
    Field omega_3_;  ///< recolouring
    Field source_;   ///< source term S
};

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_SOLVER_H
