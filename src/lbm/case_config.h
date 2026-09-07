#ifndef CGLBM_LBM_CASE_CONFIG_H
#define CGLBM_LBM_CASE_CONFIG_H

#include <functional>
#include <string>

#include "lbm/d2q9.h"
#include "lbm/equation_of_state.h"
#include "lbm/isotropic_gradient.h"
#include "lbm/isotropic_gradient_3d.h"

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

/// The same, for a three-dimensional case.
using PhaseFieldInit3D = std::function<double(const CaseConfig&, int i, int j, int k)>;

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
/// when rho1 == rho2. Equivalently phi_N = 2c - 1, where c is the volume
/// fraction of component 1: it is the volume-fraction indicator, and `phi_N = 0`
/// is where the two components occupy equal volume.
///
/// Which one to take the gradient of is not free: it decides *where* the
/// surface-tension operator acts. With `Colour` the operator is centred on the
/// zero of phi, five nodes outside the droplet at a ratio of 1000; with
/// `BulkNormalised` it sits on the interface. The integrated tension is the
/// same either way -- measured on one relaxed state, `sum |grad phi| = 2.00000`
/// against `sum |grad phi_N| = 1.99939`.
///
/// Moving it is only safe together with `SurfaceTension::ContinuumSurfaceForce`.
/// The stress form of the operator carries a `1 / tau` (see `SurfaceTension`),
/// and putting it where the heavy fluid is puts it where `tau` is largest: the
/// jump falls to 6 % of `sigma / R` at a ratio of 1000. The two options belong
/// together, and the shipped Laplace case sets both.
///
/// A note on an earlier measurement recorded here, which said the normalised
/// field was simply worse. It was taken with `phi_n_` built once in
/// `initialize()` and never refreshed -- so the gradient was of a field frozen
/// at t = 0 -- and with the initial profile and the tension operator both left
/// as they were. All three are fixed now.
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

/// How the surface tension is applied to the populations.
///
/// Both forms produce the same force in the continuum limit and differ in how
/// they reach the momentum equation.
///
/// `Perturbation` is the operator of Lafarge et al.: a capillary *stress*
/// `sigma |grad phi| (nn - I)` injected into the second-order moment, with no
/// curvature ever formed. Because a post-collision addition reaches the
/// momentum flux multiplied by the relaxation time, the operator carries a
/// `1 / tau` of its own so the two cancel. That cancellation is exact only
/// where `tau` is uniform over the interface, and `tau = rho nu / (p dt) + 1/2`
/// spans the whole density ratio: on the Laplace case at a ratio of 1000 it
/// runs from 5.5 in the light fluid to 5.0e3 in the heavy one, three nodes
/// apart. The stress the light side injects is then damped ~400x faster than
/// it accumulates on the heavy side, and the tension collapses -- measured, the
/// jump falls to 6 % of `sigma / R` once the interface is placed where the mass
/// actually is.
///
/// `ContinuumSurfaceForce` is the operator of Ba et al.: the curvature is
/// formed explicitly and the tension enters as a body force
///
///     F_s = -1/2 sigma K grad phi_N,   K = -div n,   n = -grad phi_N / |grad phi_N|
///
/// which reaches the momentum equation through Guo's forcing, whose factor
/// `1 - 1/(2 tau)` stays in [1/2, 1] however large `tau` grows. It is what lets
/// the density ratio reach 1000 here.
///
/// Reference
///  - Y. Ba, H. Liu, Q. Li, Q. Kang, J. Sun, Phys. Rev. E 94, 023310 (2016),
///    Eqs. (23)-(29): the perturbation operator in continuum-surface-force
///    form, the curvature, and the velocity redefinition that goes with it.
enum class SurfaceTension {
    Perturbation,          ///< capillary stress in Omega^(2), Lafarge et al.
    ContinuumSurfaceForce  ///< body force from an explicit curvature, Ba et al.
};

/// How the two components are pushed back apart after each collision.
///
/// Both forms push colour along the interface normal with the same
/// `(1 - phi^2)` shape, and in this code's variables they differ only in what
/// sets the strength:
///
///     InterfaceWidth:  omega_3 = W_i  p (1 - phi^2) (xi.n) / (2 w cs^2)
///     LatvaKokko:      omega_3 = W_i  beta rho (1 - phi^2) cos(phi_i) / 2
///
/// The difference is `p` against `rho`, and at high density contrast that is
/// the density ratio itself. The pressure is nearly uniform across an interface
/// -- 0.36 against 0.33 on the Laplace case -- while the density spans the whole
/// contrast, so `InterfaceWidth` segregates the heavy side about `rho1/rho2`
/// times more weakly than Latva-Kokko would. Matching the two in the light bulk
/// of the shipped case gives beta = p / (w cs^2 rho) = 0.625, close to the 0.7
/// Ba et al. use; on the heavy side the same expression gives 6.2e-4.
///
/// `cos(phi_i)` is the cosine of the angle between the colour gradient and
/// `xi_i`, so the Latva-Kokko term carries `(xi_i.n) / |xi_i|` where the other
/// carries `(xi_i.n)`. That is the form the literature states, and it weights
/// the diagonals by 1/sqrt(2) relative to the axes.
///
/// **Measured, and `InterfaceWidth` is the default for a reason.** Latva-Kokko
/// is the better operator at a density ratio of 2 (+0.3 % on the Laplace jump
/// against +2.4 %) and comes apart above about 10; by 1000 no beta in
/// [0.01, 0.9] survives. It also leaves the range the equation of state is
/// defined on -- |phi| reaches 1.00058 at a ratio of 2, where every other
/// configuration here holds 1 to within 1e-9 -- because `g` is rebuilt each
/// step as `f * phi + omega_3` rather than streamed as two component
/// populations, so the positivity that bounds phi in the two-population form is
/// not available. See docs/numerics.md.
///
/// Reference
///  - M. Latva-Kokko, D. H. Rothman, "Diffusion properties of gradient-based
///    lattice Boltzmann models of immiscible fluids", Phys. Rev. E 71, 056702
///    (2005). The segregation operator, which removes the lattice pinning of
///    the original colour-gradient recolouring.
///  - Y. Ba, H. Liu, Q. Li, Q. Kang, J. Sun, Phys. Rev. E 94, 023310 (2016),
///    Eq. (30), which is where the form implemented here is taken from.
enum class Recolouring {
    InterfaceWidth,  ///< p (1 - phi^2) / (2 w), Lafarge et al.
    LatvaKokko       ///< beta rho (1 - phi^2) / 2, Latva-Kokko & Rothman
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

    /// Kinematic shear viscosity of component 1 (phi = +1).
    ///
    /// The relaxation time is `tau = rho nu / (p dt) + 1/2`, so with one
    /// viscosity for both components the *dynamic* viscosity `rho nu` -- and
    /// with it tau -- spans the whole density ratio: 5.5 in the light fluid
    /// against 5.0e3 in the heavy one at a ratio of 1000 on the Laplace case.
    /// That is what makes the capillary stress of `SurfaceTension::Perturbation`
    /// fail to cancel its own `1 / tau`, and it is the reason a case at high
    /// density contrast may want to set the two separately. Choosing
    /// `nu2 = nu * rho1 / rho2` matches the dynamic viscosities and makes tau
    /// uniform across the interface.
    double nu = 0.0;
    double nu_b = 0.0;  ///< kinematic bulk viscosity of component 1

    /// Kinematic viscosities of component 2 (phi = -1).
    ///
    /// Negative means "the same as component 1", which is what every case
    /// shipped here uses and what the solver takes as a fast path: the
    /// interpolation below is skipped entirely, so a single-viscosity case is
    /// unaffected to the last bit.
    ///
    /// Where they differ, the kinematic viscosity at a node is interpolated on
    /// the volume fraction of component 1, `c = (1 + phi_N) / 2`:
    ///
    ///     nu(x) = c nu_1 + (1 - c) nu_2
    ///
    /// which is `nu_1` and `nu_2` exactly in the two bulks. Ba et al. Eq. (22)
    /// instead blends the relaxation *rate* parabolically over a band around
    /// the interface; both are smooth and agree in the bulks, and the volume
    /// fraction is used here because it is the same quantity `phi_N` already
    /// measures and because it reduces to the single-viscosity case exactly.
    double nu2 = -1.0;
    double nu_b2 = -1.0;

    double sigma = 0.0;   ///< surface tension
    double radius = 0.0;  ///< prescribed interface radius, for sigma / R

    /// Uniform acceleration along -y. Zero switches the body force off
    /// entirely rather than multiplying by zero, so a case without gravity
    /// carries no force term at all.
    double gravity = 0.0;

    double ch_width_init = 1.1;  ///< interface width used to lay down phi at t = 0
    double ch_width_ope = 1.6;   ///< interface width `Recolouring::InterfaceWidth` maintains

    /// Segregation strength of `Recolouring::LatvaKokko`, in [0, 1].
    ///
    /// Plays the role `ch_width_ope` plays for the other operator: larger beta
    /// is a thinner interface. Ba et al. use 0.7, which they report as an
    /// interface 4 to 5 nodes wide. Ignored by `Recolouring::InterfaceWidth`.
    double beta = 0.7;

    double p1_inf = 0.0;  ///< pressure at infinity of component 1
    double p2_inf = 0.0;  ///< pressure at infinity of component 2

    /// Rest-particle weight of component 2, for `TwoPopulationSolver` only.
    ///
    /// That model gets its density ratio from the equilibrium rather than from
    /// an equation of state: each fluid keeps its own rest weight `alpha_k`,
    /// hence its own sound speed `(c_s^k)^2 = 3/5 (1 - alpha_k)`, and
    ///
    ///     rho1 / rho2 = (1 - alpha_2) / (1 - alpha_1)
    ///
    /// fixes `alpha_1` once `alpha_2` is chosen. Both must lie in [0, 1], so
    /// `alpha_2` is what bounds the reachable density ratio: with the 0.2 Ba
    /// et al. use, `1 - alpha_1 = 0.8 rho2/rho1`, which stays positive at any
    /// ratio. `Solver` ignores this field.
    double alpha2 = 0.2;
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

/// Inverse of :func:`normalised_phase`: the colour field with a given `phi_n`.
///
///     phi = (S phi_N - D) / (S - D phi_N),   S = rho1 + rho2,  D = rho2 - rho1
///
/// The denominator is `rho1 (1 + phi_N) + rho2 (1 - phi_N)`, positive for any
/// `phi_n` in [-1, 1] and any positive densities. Used to lay an interface down
/// in the normalised field, which is the one whose zero contour is the
/// interface; see `CaseConfig::initial_profile_field`.
double phase_from_normalised(double phi_n, double rho1, double rho2);

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
    /// Lattice nodes along z. Read only by the three-dimensional solver; a
    /// two-dimensional case leaves it at 1.
    int nz = 1;

    int steps = 1000;    ///< time steps to advance
    int interval = 100;  ///< write the CSV grids every this many steps

    LatticeUnits units;
    Physics physics;

    Boundary boundary = Boundary::PeriodicY;
    PhaseFieldInit initial_phase;
    /// The same for a three-dimensional case; the solver uses whichever of the
    /// two its dimensionality calls for.
    PhaseFieldInit3D initial_phase_3d;

    /// How the density and pressure are laid down at t = 0.
    ///
    /// `MechanicalEquilibrium` is what lets a case start at a density ratio
    /// past about a hundred; `LinearDensity` reproduces the earlier runs.
    InitialState initial_state = InitialState::MechanicalEquilibrium;

    /// Which field `initial_phase` prescribes its profile in.
    ///
    /// A case says where its interface starts by handing back a tanh profile of
    /// a prescribed radius or height. The question this answers is *of which
    /// field* -- and the two choices do not describe the same droplet.
    ///
    /// The volume fraction of the heavy component is
    ///
    ///     c = rho_1 / rho1 = rho (1 + phi) / (2 rho1)
    ///
    /// and the physical interface is c = 1/2, where the two components occupy
    /// equal volume. In terms of the colour field that sits at
    /// `phi = (rho1 - rho2) / (rho1 + rho2)`, not at `phi = 0`: at a density
    /// ratio of 1000 the half-volume point is `phi = 0.998`. Prescribing the
    /// tanh in `phi` therefore puts the *droplet* nowhere near the prescribed
    /// radius -- measured on the Laplace case, a droplet asked for R = 10 is
    /// born at R = 8.41 at a ratio of 20 and at R = 6.27 at 1000 -- while
    /// `p1_inf` still carries the jump `sigma / R` for the radius that was
    /// asked for. The case starts out of equilibrium by that mismatch, and the
    /// mismatch grows with the density ratio.
    ///
    /// `BulkNormalised` prescribes the profile in `phi_N = 2c - 1` instead and
    /// inverts it through :func:`phase_from_normalised`, so the interface
    /// starts where the case asked for it at any density ratio. It is the
    /// default. `Colour` reproduces the earlier behaviour.
    InterfaceField initial_profile_field = InterfaceField::BulkNormalised;

    /// Which field the colour gradient is taken of. See `InterfaceField`.
    InterfaceField interface_field = InterfaceField::Colour;

    /// How the surface tension reaches the populations. See `SurfaceTension`.
    SurfaceTension surface_tension = SurfaceTension::Perturbation;

    /// How the two components are segregated after collision. See `Recolouring`.
    Recolouring recolouring = Recolouring::InterfaceWidth;

    /// Keep `physics.p1_inf` at the value `matched_p1_inf` gives.
    ///
    /// A case that sets `p1_inf = matched_p1_inf(physics)` has tied it to the
    /// densities, the sound speeds, sigma and the radius. `--rho1`, `--sigma`
    /// and the rest would then leave it stale, which is a silent wrong answer
    /// rather than a visible one: the case would start with the wrong pressure
    /// jump. Setting this makes `parse_command_line` recompute it after the
    /// overrides are in. `TwoPopulationSolver` ignores it -- that model has no
    /// pressure offset to match.
    bool matched_pressure_offset = false;

    /// Isotropy order of the colour gradient.
    ///
    /// The wall-bounded cases default to E4: it reaches one node, so it stays
    /// exact next to a wall, while a wider stencil becomes one-sided over two
    /// nodes there and has not been validated against those cases.
    GradientStencil stencil = GradientStencil::E8;

    /// Isotropy order of the colour gradient in three dimensions.
    ///
    /// `E6` reaches two nodes along each axis and is sixth-order isotropic;
    /// `E4` is the D3Q19 neighbourhood and only fourth-order. In two dimensions
    /// the equivalent step was worth an order of magnitude on the spurious
    /// currents, so `E6` is the default.
    GradientStencil3D stencil_3d = GradientStencil3D::E6;

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

    /// Also write the whole field as `field_<t>.vtk`, every `interval` steps.
    ///
    /// Only the three-dimensional solver reads this, and it is off by default
    /// because of what it costs: a 48^3 run writes 14 MB per dump against
    /// 300 kB for the four CSV slices, and the slices are what the shipped
    /// post-processing reads. It is here because without it a three-dimensional
    /// run has no way at all to show what happens off the mid-plane -- CSV is
    /// not offered for a whole field, at these sizes it would be gigabytes --
    /// so a case that needs ParaView asks for `--vtk`.
    bool write_vtk_field = false;

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
/// `inverted` flips the profile, so that component 2 -- the light one in every
/// case shipped here -- is the one above the interface. It is therefore
/// `inverted = false` that puts the dense component on top and makes the
/// Rayleigh-Taylor case unstable, which is what that case passes.
PhaseFieldInit cosine_layer(double amplitude, bool inverted);

/// A spherical droplet of `physics.radius` centred in the domain, component 1
/// inside, as a profile in the bulk-normalised phase field.
PhaseFieldInit3D droplet_interface_3d();

/// The same droplet, deformed into the second spherical harmonic.
///
/// The surface is `r(theta) = R [1 + eps P_2(cos theta)]` with theta measured
/// from the z axis and `P_2(x) = (3x^2 - 1)/2`, so the droplet starts as a
/// spheroid: prolate along z for `deformation > 0`, oblate for negative. It is
/// the shape whose free oscillation Lamb's mode-2 frequency describes, and
/// releasing it from rest is how `oscillation_3d` measures that frequency.
///
/// `R` is scaled by `(1 + 3 eps^2/5)^(-1/3)` so the spheroid encloses the
/// volume of a sphere of `physics.radius`. That correction is second order in
/// eps -- 0.2 % at eps = 0.1 -- but it is the difference between a droplet that
/// relaxes to the radius the case asked for and one that does not, and the
/// frequency goes as `R^(-3/2)`.
///
/// Reference
///  - H. Lamb, *Hydrodynamics*, 6th ed., Cambridge (1932), art. 275. The
///    frequency of the free oscillations of a liquid globe.
PhaseFieldInit3D oscillating_droplet_3d(double deformation);

/// A flat interface at mid-height, perturbed by one cosine wavelength along
/// both x and z.
///
/// The three-dimensional counterpart of `cosine_layer`, with the same
/// displacement `amplitude * nx` and the same meaning for `inverted`. The
/// perturbation is the product `cos(2 pi x / nx) cos(2 pi z / nz)`, the single
/// mode of a square cell: the interface is displaced upward by the full
/// amplitude at the cell's four corners and at its centre, downward by it at
/// the midpoint of each edge, and not at all along the lines between. With the
/// dense component on top those are the bubbles and the spikes respectively.
PhaseFieldInit3D cosine_layer_3d(double amplitude, bool inverted);

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
