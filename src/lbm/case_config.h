#ifndef CGLBM_LBM_CASE_CONFIG_H
#define CGLBM_LBM_CASE_CONFIG_H

#include <functional>
#include <string>

#include "lbm/d2q9.h"
#include "lbm/equation_of_state.h"
#include "lbm/isotropic_gradient.h"

/// What distinguishes one colour-gradient case from another.
///
/// The five solvers shipped here run the same scheme; they differ only in the
/// lattice, the fluid properties, how the domain is closed along y, and how the
/// phase field starts. Those differences used to be spelled out by copying the
/// whole program -- `capillary` and `gravity_capillary` were 588-line files
/// that disagreed on seven lines. They are gathered here instead, so a case is
/// a value and the scheme has one implementation.
///
/// Everything is a runtime value. The lattice size and the step count used to
/// be `const int` at file scope, which meant a rebuild to change a resolution
/// and a compile-time constant that a test could only duplicate and hope stayed
/// in sync.

namespace cglbm {
namespace lbm {

struct CaseConfig;

/// Initial phase field at node (i, j), given the case it belongs to.
///
/// phi = +1 is component 1 and phi = -1 component 2; the interface is where it
/// crosses zero.
using PhaseFieldInit = std::function<double(const CaseConfig&, int i, int j)>;

/// Which field the colour gradient is taken of.
///
/// The surface-tension and recolouring operators both act along the gradient
/// of a phase field, and both are centred on where that field crosses zero. The
/// question is which field marks the interface.
///
/// `Colour` uses phi = (rho_1 - rho_2) / rho directly, as the model has since
/// the beginning. Its zero contour is the interface only when the two
/// components have the same bulk density: in general phi = 0 sits at
///
///     phi = (rho1 - rho2) / (rho1 + rho2)
///
/// away from the true interface, which is 0.90 at a density ratio of 20 and
/// 0.9998 at 10^4. The operators are then centred inside the light fluid rather
/// than on the interface, and further inside it the larger the ratio.
///
/// `BulkNormalised` divides each component's density by its own bulk value
/// before forming the phase field:
///
///     phi_N = (rho_1/rho1_0 - rho_2/rho2_0) / (rho_1/rho1_0 + rho_2/rho2_0)
///           = (rho2 (1 + phi) - rho1 (1 - phi)) / (rho2 (1 + phi) + rho1 (1 - phi))
///
/// so that phi_N = 0 is the interface at any density ratio, and phi_N = phi
/// when rho1 == rho2.
///
/// **Measured, and it does not help here.** On the Laplace case it degrades the
/// pressure jump monotonically with the density ratio, and it neither closes
/// the gap between the phase and density interfaces nor raises the ratio at
/// which the case diverges:
///
///     ratio            2      5     10     20     50    100
///     colour       0.956  0.925  0.933  0.962  1.016  1.013
///     normalised   0.983  0.905  0.809  0.696  0.554  0.445
///
/// The reason is that Ba et al. use the normalised field only for the interface
/// *normal*, in a continuum-surface-force operator where sigma enters
/// explicitly. The operator here, from Lafarge et al., carries its calibration
/// in the *magnitude* of the gradient: swapping in a field with a different
/// profile changes the surface tension it produces. Transferring Eq. (21) to
/// this operator would mean recalibrating Omega^(2) with it, which has not been
/// done. `Colour` is therefore the default; `BulkNormalised` is kept because it
/// is what the literature specifies and is the starting point for that work.
///
/// Reference
///  - Y. Ba, H. Liu, Q. Li, Q. Kang, J. Sun, "Multiple-relaxation-time
///    color-gradient lattice Boltzmann model for simulating two-phase flows
///    with high density ratio", Phys. Rev. E 94, 023310 (2016),
///    doi:10.1103/PhysRevE.94.023310, Eq. (21). "This definition, however,
///    becomes increasingly incorrect in identifying the interface as the
///    density ratio increases."
///  - S. Leclaire, M. Reggio, J.-Y. Trepanier, "Progress and investigation on
///    lattice Boltzmann modeling of multiple immiscible fluids or components
///    with variable density and viscosity ratios", J. Comput. Phys. 246, 318
///    (2013), which first replaced the gradient of the raw colour field.
enum class InterfaceField {
    Colour,         ///< grad phi, the historical behaviour
    BulkNormalised  ///< grad phi_N, Ba et al. Eq. (21)
};

/// How the density and pressure are laid down at t = 0, given the phase field.
///
/// The phase field starts as a tanh profile, which is the shape the recolouring
/// operator maintains. What density and pressure go with it is a separate
/// choice, and at high density ratio it decides whether the case survives its
/// own first few steps.
///
/// `LinearDensity` interpolates the density linearly in phi and takes the
/// pressure from a linear mixing rule. The time loop, however, recovers the
/// pressure from the two-component equation of state, and the two disagree
/// across the interface by a factor that grows with the density ratio: 3.5 at
/// a ratio of 10, 32 at 100, 316 at 1000. The solver applies its own value at
/// the first step, and that discontinuity is an acoustic blast.
///
/// `EquationOfStateP` keeps the linear density but takes the initial pressure
/// from the equation of state, removing that first-step jump. It does not
/// remove the underlying problem: a density interpolated linearly in phi is not
/// in mechanical equilibrium, and the equation of state reads a pressure of
/// about 0.31 rho1 cs^2 at the interface against bulk pressures of order cs^2.
/// The scheme does relax that away -- measured peak pressure at a ratio of 100
/// falls from 10.44 to 0.43 over 3000 steps -- but the transient reaches Mach
/// 1.4, and past a ratio of about 100 it destroys the run first.
///
/// `MechanicalEquilibrium` instead solves the equation of state for the density
/// that makes the pressure smooth across the interface: for each node, rho is
/// found such that p(rho, phi) equals a target interpolating the two bulk
/// pressures. Both bulks come out exactly rho1 and rho2, so the Laplace jump at
/// t = 0 is unchanged; only the interface profile differs, and it starts where
/// the scheme would have taken it anyway.
enum class InitialState {
    /// Density linear in phi, pressure from the linear mixing rule. The
    /// historical initialisation; kept so earlier runs can be reproduced.
    LinearDensity,
    /// Density linear in phi, pressure from the equation of state.
    EquationOfStateP,
    /// Density solved so the interface starts in mechanical equilibrium.
    MechanicalEquilibrium
};

/// How the domain is closed along y. Both cases are periodic along x.
enum class Boundary {
    /// Periodic on both axes. The droplet cases float in an unbounded fluid.
    PeriodicY,
    /// Resting walls at j = 0 and j = ny - 1, applied as half-way bounce-back.
    WallY
};

/// The two components and the transport coefficients of the mixture.
///
/// Everything is in lattice units. `c_dx` and `c_dt` in the case programs are
/// the only place a physical unit appears; by the time a value reaches here it
/// has been converted.
struct Physics {
    double rho1 = 1.0;  ///< density of component 1 (phi = +1)
    double rho2 = 1.0;  ///< density of component 2 (phi = -1)
    double c1 = 1.0;    ///< sound speed of component 1
    double c2 = 1.0;    ///< sound speed of component 2

    double nu = 0.0;    ///< kinematic shear viscosity
    double nu_b = 0.0;  ///< kinematic bulk viscosity

    double sigma = 0.0;   ///< surface tension
    double radius = 0.0;  ///< prescribed interface radius, for sigma / R

    /// Uniform acceleration along -y. Zero switches the body force off
    /// entirely rather than multiplying by zero, so a case without gravity
    /// carries no force term at all.
    double gravity = 0.0;

    double ch_width_init = 1.1;  ///< interface width used to lay down phi at t = 0
    double ch_width_ope = 1.6;   ///< interface width the recolouring operator maintains

    double p1_inf = 0.0;  ///< pressure at infinity of component 1
    double p2_inf = 0.0;  ///< pressure at infinity of component 2
};

/// Density at which the equation of state gives `target_pressure` for `phi`.
///
/// Used to start an interface in mechanical equilibrium. Solved by bisection on
/// a bracket around the two bulk densities; returns the linear interpolation
/// `rho1 (1+phi)/2 + rho2 (1-phi)/2` unchanged if the bracket holds no root,
/// so a case can never fail to initialise because of this.
double density_at_pressure(
    double phi, double target_pressure, double rho1, double rho2, const ComponentPair& components);

/// Phase field normalised by the two bulk densities, Ba et al. Eq. (21).
///
/// Maps `phi` in [-1, 1] onto [-1, 1] so that the zero contour marks the
/// interface at any density ratio:
///
///     phi_N = (rho2 (1 + phi) - rho1 (1 - phi)) / (rho2 (1 + phi) + rho1 (1 - phi))
///
/// It fixes the two bulks, `phi_N(+-1) = +-1`, is the identity when
/// `rho1 == rho2`, and crosses zero at `phi = (rho1 - rho2) / (rho1 + rho2)`,
/// which is where the two components are present in equal proportion of their
/// own bulk densities. `phi` outside [-1, 1] is clamped; the denominator is
/// positive on that range for any positive densities.
double normalised_phase(double phi, double rho1, double rho2);

/// The `p1_inf` that makes the Laplace jump exact at t = 0.
///
/// Matching the two ideal-gas branches across an interface of radius R while
/// carrying a jump of sigma / R fixes the offset:
///
///     p1_inf = rho1 c1^2 - rho2 c2^2 - sigma / R
///
/// Every shipped case uses it, including those with sigma = 0, where it reduces
/// to matching the two bulk pressures.
double matched_p1_inf(const Physics& physics);

/// One case: a lattice, a fluid, a boundary, an initial state, and how long to
/// run it.
struct CaseConfig {
    /// Name used in the run log. Not a file name -- output names are fixed.
    std::string name;

    int nx = 128;  ///< lattice nodes along x
    int ny = 128;  ///< lattice nodes along y

    int steps = 1000;    ///< time steps to advance
    int interval = 100;  ///< write the CSV grids every this many steps

    LatticeUnits units;
    Physics physics;

    Boundary boundary = Boundary::PeriodicY;
    PhaseFieldInit initial_phase;

    /// How the density and pressure are laid down at t = 0.
    ///
    /// `MechanicalEquilibrium` is what lets a case start at a density ratio
    /// past about a hundred; `LinearDensity` reproduces the earlier runs.
    InitialState initial_state = InitialState::MechanicalEquilibrium;

    /// Which field the colour gradient is taken of.
    ///
    /// `Colour` is the default: the normalised field measures worse with the
    /// surface-tension operator as it is currently calibrated. See
    /// `InterfaceField`.
    InterfaceField interface_field = InterfaceField::Colour;

    /// Isotropy order of the colour gradient.
    ///
    /// The wall-bounded cases default to E4: it reaches one node, so it stays
    /// exact next to a wall, while a wider stencil becomes one-sided over two
    /// nodes there and has not been validated against those cases.
    GradientStencil stencil = GradientStencil::E8;

    /// Append the interface position along the mid-plane to `interface.csv`,
    /// every step. The capillary cases measure an oscillation period from it.
    bool track_interface = false;

    /// Run the lattice loops across OpenMP threads.
    ///
    /// Off by default so that a case keeps serial semantics even in an
    /// OpenMP-enabled build. Every loop but the streaming step writes only its
    /// own node and so gives identical results either way.
    bool parallel = false;

    /// Report a phase field that has left [-1, 1] on stdout.
    bool warn_phase_out_of_range = false;

    /// Significant digits in the CSV output.
    ///
    /// Six is the iostream default and what every run so far has written; it
    /// loses about ten digits of a double. Raise it to 17 for output that
    /// round-trips.
    int output_precision = 6;

    /// Number of nodes, `nx * ny`.
    long node_count() const {
        return static_cast<long>(nx) * static_cast<long>(ny);
    }
};

/// Ready-made initial phase fields.
///
/// A tanh profile of width `ch_width` is the equilibrium shape of the
/// recolouring operator, so starting from one avoids a transient in which the
/// interface finds its own thickness.
///
/// A droplet of `physics.radius` centred in the domain, component 1 inside.
PhaseFieldInit droplet_interface();

/// A flat interface at mid-height, perturbed by one cosine wavelength across
/// the domain with amplitude `amplitude * nx`.
///
/// `inverted` puts the dense component on top, which is what makes the
/// Rayleigh-Taylor case unstable.
PhaseFieldInit cosine_layer(double amplitude, bool inverted);

/// Outcome of reading the command line.
enum class CommandLineResult {
    Run,       ///< carry on and run the case
    Finished,  ///< nothing to do, exit 0 (`--help`)
    Error      ///< bad usage, exit non-zero
};

/// Apply command-line overrides to `config`.
///
/// The historical form, a bare stencil name (`laplace E4`), is still accepted.
/// Beyond it every option is `--key=value`:
///
///     --stencil=E4|E6|E8   --nx=N        --ny=N
///     --steps=N            --interval=N  --precision=N
///     --threads=N          --help
///
/// Diagnostics go to stderr; `program_name` is what they call the program.
CommandLineResult
parse_command_line(CaseConfig& config, int argc, char** argv, const std::string& program_name);

/// The case as `key = value` lines, one per line, for the run log.
///
/// The test suite reads these back rather than duplicating the constants, so a
/// parameter can only be stated once. `pycglbm.testing.parse_key_values`
/// understands the format.
std::string describe(const CaseConfig& config);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_CASE_CONFIG_H
