#ifndef CGLBM_LBM_CASE_CONFIG_H
#define CGLBM_LBM_CASE_CONFIG_H

#include <functional>
#include <string>

#include "lbm/d2q9.h"
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
    long node_count() const { return static_cast<long>(nx) * static_cast<long>(ny); }
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
CommandLineResult parse_command_line(CaseConfig& config,
                                     int argc,
                                     char** argv,
                                     const std::string& program_name);

/// The case as `key = value` lines, one per line, for the run log.
///
/// The test suite reads these back rather than duplicating the constants, so a
/// parameter can only be stated once. `pycglbm.testing.parse_key_values`
/// understands the format.
std::string describe(const CaseConfig& config);

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_CASE_CONFIG_H
