#ifndef CGLBM_LBM_VELOCITY_BASED_H
#define CGLBM_LBM_VELOCITY_BASED_H

/// A velocity-based two-phase scheme, for interfaces that move at large density
/// ratios, with a momentum exchange that conserves momentum exactly.
///
/// The colour-gradient solvers stream populations whose moments are the
/// density rho and the momentum rho u. Both jump by the density ratio across
/// an interface two or three nodes wide. Keeping the populations non-negative
/// in the heavy fluid needs p / rho >~ |u|, which at a density ratio of 1e4
/// holds only below |u| ~ 3e-5. Past that, the heavy side streams populations
/// of about -rho u / 2 into light nodes that hold a mass of order one, and a
/// droplet moving at 1e-3 lattice units per step diverges within a few hundred
/// steps (docs/numerics.md).
///
/// Here the density is never streamed. The hydrodynamic populations g carry
///
///     P = p / (rho cs^2)   and   u,
///
/// both continuous across the interface, and the density follows from the
/// volume fraction c of component 1, rho = rho_2 + c (rho_1 - rho_2). c is
/// carried by a second, memoryless set of populations h whose carrier is the
/// Maxwellian Gamma_i(u) >= 0, so c stays in [0, 1]; the sharpening along the
/// interface normal is the colour-gradient recolouring, and together they
/// solve the conservative Allen-Cahn equation
///
///     dc/dt + div(c u) = div( M (grad c - 2 c (1 - c) / W  n) ),   M = cs^2 / 2,
///
/// whose equilibrium is c = (1 + tanh(x / W)) / 2. The components are
/// incompressible: the lattice sound speed cs is an artificial compressibility
/// in both fluids, and the equation of state of the colour-gradient solvers is
/// not used.
///
/// Streaming u makes the scheme robust but, read as it stands, not
/// conservative: a population leaving node y carries rho(y) g_i per unit of
/// velocity, and arrives at x counted with rho(x). Across the interface the two
/// differ by up to the density ratio. The momentum is instead updated link by
/// link, every link exchanging equal and opposite amounts (link_momentum):
///
///  - the lattice's own non-equilibrium and thermal exchange, weighted by the
///    lighter density of the two ends. A light node sees its heavy neighbour
///    as the velocity-based scheme does; a heavy node sees a light neighbour
///    through the light fluid's momentum only, as a free surface;
///  - the advection, as the mass flux of the phase populations times the mean
///    velocity of the two ends, so that mass and momentum move together (the
///    consistent mass and momentum transport of the phase-field literature).
///    A uniform velocity is then carried by the mass flux alone, whatever the
///    densities on either side of a link.
/// The pressure force, the lattice gradient of p itself (pressure_force), and
/// the capillary stress divergence both sum to zero over the lattice. What a
/// node holds after the exchanges, divided by its new density, is its new
/// velocity; set_velocity hands it back to the populations.
///
/// Two more link terms are diffusions of the velocity: the upwind part of the
/// mass flux the lighter density does not account for, which a node the heavy
/// fluid is leaving needs, and the viscous stress the lighter-density
/// weighting takes away, restored to the harmonic mean of the two dynamic
/// viscosities. Neither goes into the exchange. At tau close to 1/2 the
/// non-equilibrium flips sign every step, and a dissipation added to the
/// velocity within the step, which the populations never see, drives that
/// period-2 mode: the same diffusion that is harmless through the forcing term
/// makes a single fluid diverge in about 400 steps at |u| = 0.1. link_momentum
/// hands their coefficients back, and dissipation_force applies them as a
/// force, half before and half after the step like any other.
///
/// The collision relaxes the trace of the non-equilibrium at its own rate: a
/// heavy fluid with a modest dynamic viscosity has tau close to 1/2, and
/// without a separate bulk relaxation its acoustic modes are undamped. Its
/// shear modes, as nearly undamped, are held by the hybrid regularised
/// collision (collide_hybrid).
///
/// References
///  - A. Fakhari, T. Mitchell, C. Leonardi, D. Bolster, "Improved locality of
///    the phase-field lattice-Boltzmann model for immiscible fluids at high
///    density ratios", Phys. Rev. E 96, 053301 (2017). The velocity-based
///    hydrodynamic equilibrium.
///  - Y. Q. Zu, S. He, "Phase-field-based lattice Boltzmann model for
///    incompressible binary fluid systems with density and viscosity
///    contrasts", Phys. Rev. E 87, 043301 (2013).
///  - P.-H. Chiu, Y.-T. Lin, "A conservative phase field method for solving
///    incompressible two-phase flows", J. Comput. Phys. 230, 185-204 (2011).
///    The conservative Allen-Cahn equation.
///  - Z. Huang, G. Lin, A. M. Ardekani, "Consistent and conservative scheme
///    for incompressible two-phase flows using the conservative Allen-Cahn
///    model", J. Comput. Phys. 420, 109718 (2020); S. Mirjalili, A. Mani,
///    "Consistent, energy-conserving momentum transport for simulations of
///    two-phase flows using the phase field equations", J. Comput. Phys. 426,
///    109918 (2021). Momentum convected by the mass flux of the phase-field
///    equation, including its Allen-Cahn part.
///  - C. Zhan, Z. Chai, B. Shi, "Consistent and conservative phase-field-based
///    lattice Boltzmann method for incompressible two-phase flows", Phys. Rev.
///    E 106, 025319 (2022). The same consistency in a lattice Boltzmann model.
///  - H. Otomo et al., "Lattice Boltzmann models for the hydrodynamic equations
///    in multiphase flow with high density ratio", arXiv:2512.01027 (2025).
///    Why schemes that stream rho u lose accuracy near a large density jump,
///    and stream u instead.
///  - O. Malaspinas, "Increasing stability and accuracy of the lattice
///    Boltzmann scheme: recursivity and regularization", arXiv:1505.06900
///    (2015); J. Jacob, O. Malaspinas, P. Sagaut, "A new hybrid recursive
///    regularised Bhatnagar-Gross-Krook collision model for Lattice Boltzmann
///    method-based large eddy simulation", J. Turbul. 19, 1051-1076 (2018).
///    The regularised and hybrid regularised collisions.

namespace cglbm {
namespace lbm {
namespace velocity_based {

/// D2Q9, in the order of the solvers: rest, the four axes, the four diagonals.
constexpr int kQ = 9;
constexpr int kVelocity[kQ][2] = {
    {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
constexpr double kWeight[kQ] = {4.0 / 9.0,
                                1.0 / 9.0,
                                1.0 / 9.0,
                                1.0 / 9.0,
                                1.0 / 9.0,
                                1.0 / 36.0,
                                1.0 / 36.0,
                                1.0 / 36.0,
                                1.0 / 36.0};
/// The direction opposite each of kVelocity.
constexpr int kOpposite[kQ] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
constexpr double kSoundSpeedSquared = 1.0 / 3.0;

/// Gamma_i(u), the second-order Maxwellian of unit mass. Its moments are 1, u
/// and cs^2 I + u u, and it is non-negative for |u| below about 0.4.
void velocity_equilibrium(double ux, double uy, double* gamma);

/// g_i^eq = w_i P + Gamma_i(u) - w_i, with moments P, u and P cs^2 I + u u.
void hydrodynamic_equilibrium(double pressure_number, double ux, double uy, double* equilibrium);

/// Post-collision phase populations c Gamma_i(u) + theta w_i A (xi_i . n) / cs^2.
///
/// A = M 2 c (1 - c) / W with M = cs^2 / 2 is the sharpening flux along the
/// unit normal n that balances the diffusion of the memoryless transport on
/// the profile c = (1 + tanh(x / W)) / 2. Streamed, they give the new c.
///
/// theta in [0, 1] is the largest weight that keeps every population between
/// 0 and Gamma_i(u), and with them the component-2 populations Gamma_i - h_i,
/// non-negative. Then the new c is a sum of non-negative parts, and so is
/// 1 - c up to the compressibility of the carriers. It is 1 at rest and
/// below |u| ~ 0.03 with W = 1.6; faster, the carrier against the flow drops
/// below the sharpening on the far side of the interface, and a population a
/// thousandth of a unit negative there is a negative mass as large as the
/// light node's own at a density ratio of 1e4.
void phase_populations(double c,
                       double ux,
                       double uy,
                       double normal_x,
                       double normal_y,
                       double width,
                       double* populations);

/// Guo et al. forcing populations for an acceleration a (force per unit mass).
void forcing(double ux, double uy, double ax, double ay, double* source);

/// Regularised collision with separate shear and bulk relaxation.
///
/// The non-equilibrium g - g^eq + S/2 is projected onto the second-order
/// Hermite polynomials; its deviatoric part relaxes with `tau_shear`, its
/// trace with `tau_bulk`, and the result is g^eq + non-equilibrium + S/2.
void collide(const double* populations,
             const double* equilibrium,
             const double* source,
             double tau_shear,
             double tau_bulk,
             double* post_collision);

/// A velocity gradient, du_a / dx_b.
struct VelocityGradient {
    double dux_dx;
    double dux_dy;
    double duy_dx;
    double duy_dy;
};

/// Hybrid regularised collision (Jacob, Malaspinas & Sagaut 2018).
///
/// As collide, but the second-order non-equilibrium that relaxes is
/// `sigma` times that of the populations plus 1 - sigma times its
/// Chapman-Enskog value from the velocity gradient,
///
///     -tau_shear cs^2 (grad u + grad u^T - div u I) - tau_bulk cs^2 div u I.
///
/// The two agree on resolved flow; on the grid-scale modes the populations
/// carry and the finite-difference gradient does not, the blend damps what a
/// shear relaxation time within 1e-4 of 1/2 would leave undamped. sigma = 1
/// is collide.
void collide_hybrid(const double* populations,
                    const double* equilibrium,
                    const double* source,
                    double tau_shear,
                    double tau_bulk,
                    double sigma,
                    const VelocityGradient& gradient,
                    double* post_collision);

/// Pressure force per unit mass, -grad_lat(p) / rho(x), with p = rho cs^2 P.
///
///     a(x) = sum_i w_i xi_i p(x - xi_i) / (cs^2 rho(x))
///
/// is the lattice's own pressure gradient, applied to p where it is continuous
/// rather than to P, which jumps by the density ratio. As a force density
/// rho a it sums to zero over a periodic lattice, and a uniform p exerts no
/// force whatever rho does. `pressure_number` and `rho` are `nx * ny` values
/// indexed `[i * ny + j]`; both axes are periodic.
void pressure_force(const double* pressure_number,
                    const double* rho,
                    int nx,
                    int ny,
                    int i,
                    int j,
                    double* ax,
                    double* ay);

/// One end of a lattice link, as link_momentum reads it after streaming.
struct LinkEnd {
    /// The post-collision hydrodynamic population this end sent along the
    /// link, less its pressure part w_i P.
    double outgoing;
    /// The phase population this end sent along the link.
    double phase;
    /// Density at the new time.
    double rho;
    /// Dynamic viscosity at the new time.
    double mu;
    /// Velocity at the old time, the one the equilibria were built with.
    double ux;
    double uy;
};

/// Momentum that `receiver` gains through its link with `donor`, where donor
/// sits at receiver - xi_k, and the dissipation coefficient of the link. The
/// donor loses the same momentum: swapping the two ends and taking the
/// opposite direction gives exactly the opposite vector and the same
/// coefficient.
///
/// With rho_l = min(donor.rho, receiver.rho), F the velocity-based exchange
/// per unit mass (the two `outgoing` values), F_adv its advective part
/// (the u u terms of both equilibria), K the exchange of the carriers Gamma,
/// K_lin its part linear in u, and M the mass carried by the phase populations,
///
///     J = rho_l (F - F_adv) xi_k                 thermal, viscous, forcing
///       + (M - rho_l (K - K_lin)) u_mid          advection by the mass flux
///
///     D = |M - rho_l K| / 2 + beta 2 w_k / cs^2  (dissipation_force)
///
/// u_mid is the mean of the two velocities. K - K_lin, the carriers'
/// quadratic part, is a diffusion of the velocity along the flow; carried
/// within the step it destabilises a single fluid at |u| = 0.1, and it
/// vanishes for a uniform velocity, which the mass flux still carries exactly.
/// The first term of D turns the centred advection of the excess mass flux
/// M - rho_l K into an upwind one; beta raises the link viscosity
/// rho_l (nu_donor + nu_receiver) / 2 to the harmonic mean of the dynamic
/// viscosities. Where the two densities are equal and the phase populations
/// are those of a single component, M = rho K and D = 0: the exchange is the
/// lattice's own, with its advection written as a mass flux.
void link_momentum(int k,
                   const LinkEnd& donor,
                   const LinkEnd& receiver,
                   double rho1,
                   double rho2,
                   double* jx,
                   double* jy,
                   double* dissipation);

/// Force density sum_k D_k (u(x - xi_k) - u(x)) at node (i, j), from the
/// link coefficients link_momentum returns.
///
/// `coefficients` holds kQ values per node, `[(i * ny + j) * kQ + k]` for the
/// link to (i, j) - xi_k; `ux` and `uy` are `nx * ny` values indexed
/// `[i * ny + j]`, both axes periodic. With the same coefficient at the two
/// ends of every link, the force sums to zero over the lattice.
void dissipation_force(const double* coefficients,
                       const double* ux,
                       const double* uy,
                       int nx,
                       int ny,
                       int i,
                       int j,
                       double* fx,
                       double* fy);

/// Replace the first moment of `populations` by u, leaving the zeroth and
/// second moments unchanged.
void set_velocity(double* populations, double ux, double uy);

}  // namespace velocity_based
}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_VELOCITY_BASED_H
