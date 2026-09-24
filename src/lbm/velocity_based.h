#ifndef CGLBM_LBM_VELOCITY_BASED_H
#define CGLBM_LBM_VELOCITY_BASED_H

/// A velocity-based two-phase scheme, for interfaces that move at large density
/// ratios.
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
/// Streaming P and u gives -grad(p / rho) and div(nu S) per unit mass. The
/// physics wants -grad(p) / rho and div(mu S) / rho, and two corrections make
/// up the difference:
///  - pressure_correction, written on the lattice's own stencil so that the
///    total pressure force is exactly -grad_lat(p) / rho. The finite-difference
///    form P cs^2 grad(rho) / rho cancels two terms of order rho_1 / rho_2
///    between different stencils, and diverges at a density ratio of 1e4;
///  - viscous_correction, nu S . grad(ln rho), bounded across the interface
///    where grad(rho) / rho is not.
/// The collision relaxes the trace of the non-equilibrium at its own rate: a
/// heavy fluid with a modest dynamic viscosity has tau close to 1/2, and
/// without a separate bulk relaxation its acoustic modes are undamped.
///
/// Known limitation: momentum is not conserved exactly. The scheme evolves u,
/// not rho u; the mass moved by the Allen-Cahn flux does not carry its
/// momentum, and the viscous corrections are not in conservative form. On a
/// droplet launched into quiescent fluid, total momentum drifts by 4.4 % over
/// 10 000 steps at a density ratio of 1e4, and by 11.7 % at 100. The
/// colour-gradient solvers are exactly conservative and remain the reference
/// below a density ratio of about 100.
///
/// References
///  - A. Fakhari, T. Mitchell, C. Leonardi, D. Bolster, "Improved locality of
///    the phase-field lattice-Boltzmann model for immiscible fluids at high
///    density ratios", Phys. Rev. E 96, 053301 (2017). The velocity-based
///    hydrodynamic equilibrium and its pressure and viscous corrections.
///  - Y. Q. Zu, S. He, "Phase-field-based lattice Boltzmann model for
///    incompressible binary fluid systems with density and viscosity
///    contrasts", Phys. Rev. E 87, 043301 (2013).
///  - P.-H. Chiu, Y.-T. Lin, "A conservative phase field method for solving
///    incompressible two-phase flows", J. Comput. Phys. 230, 185-204 (2011).
///    The conservative Allen-Cahn equation.

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
constexpr double kSoundSpeedSquared = 1.0 / 3.0;

/// Gamma_i(u), the second-order Maxwellian of unit mass. Its moments are 1, u
/// and cs^2 I + u u, and it is non-negative for |u| below about 0.4.
void velocity_equilibrium(double ux, double uy, double* gamma);

/// g_i^eq = w_i P + Gamma_i(u) - w_i, with moments P, u and P cs^2 I + u u.
void hydrodynamic_equilibrium(double pressure_number, double ux, double uy, double* equilibrium);

/// Post-collision phase populations c Gamma_i(u) + w_i A (xi_i . n) / cs^2.
///
/// A = M 2 c (1 - c) / W with M = cs^2 / 2 is the sharpening flux along the
/// unit normal n that balances the diffusion of the memoryless transport on
/// the profile c = (1 + tanh(x / W)) / 2. Streamed, they give the new c.
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

/// Acceleration that turns the lattice's pressure force into -grad(p) / rho.
///
/// Streaming w_i P gives -cs^2 grad_lat(P) per unit mass, with grad_lat the
/// stencil sum_i w_i xi_i f(x + xi_i) / cs^2. Adding
///
///     a(x) = sum_i w_i xi_i P(x + xi_i) (1 - rho(x + xi_i) / rho(x))
///
/// makes the total exactly -grad_lat(rho cs^2 P) / rho(x): the pressure p is
/// differentiated where it is continuous, and nothing of order rho_1 / rho_2
/// has to cancel. `pressure_number` and `rho` are `nx * ny` values indexed
/// `[i * ny + j]`; both axes are periodic.
void pressure_correction(const double* pressure_number,
                         const double* rho,
                         int nx,
                         int ny,
                         int i,
                         int j,
                         double* ax,
                         double* ay);

/// Acceleration nu S . grad(ln rho), with S = grad u + grad u^T.
///
/// The lattice gives div(nu S) per unit mass for the local nu = mu / rho; the
/// physics wants div(mu S) / rho. For any mu(c) the difference is
/// nu S . grad(ln rho), which stays bounded across the interface.
void viscous_correction(double nu,
                        double dux_dx,
                        double dux_dy,
                        double duy_dx,
                        double duy_dy,
                        double dlnrho_dx,
                        double dlnrho_dy,
                        double* ax,
                        double* ay);

}  // namespace velocity_based
}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_VELOCITY_BASED_H
