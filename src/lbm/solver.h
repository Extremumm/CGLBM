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
///            ->  macroscopic  ->  phase_field  ->  equilibrium
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
///    gradient and vanishes wherever the gradient does;
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

    const CaseConfig& config() const { return config_; }

    /// Read-only views of the state, for tests that check an invariant.
    const Field& density() const { return rho_; }
    const Field& velocity() const { return u_; }
    const Field& phase() const { return phi_; }
    const Field& pressure() const { return p_; }

  private:
    // The eight stages of a step, in the order `step()` applies them.
    void force();
    void collide();
    void collide_surface();
    void recolor();
    void stream();
    void macroscopic();
    void phase_field();
    void equilibrium();

    /// Colour gradient at (i, j), following the case's boundary along y.
    void colour_gradient(int i, int j, double* grad_x, double* grad_y) const;

    CaseConfig config_;
    int nx_;
    int ny_;

    // Lattice constants, resolved once from config_.units.
    double dx_;
    double dt_;
    double cs2_;
    double cs4_;
    double cs6_;

    bool wall_y_;    ///< config_.boundary == Boundary::WallY
    bool parallel_;  ///< config_.parallel, copied out for the OpenMP if clause
    ComponentPair components_;

    // Macroscopic fields.
    Field rho_;      ///< density
    Field rho_mdt_;  ///< density at the previous step, for the temporal correction
    Field u_;        ///< velocity, two components per node
    Field p_;        ///< pressure
    Field p_mdt_;    ///< pressure at the previous step
    Field phi_;      ///< phase field
    Field force_;    ///< external volume force, two components per node

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
