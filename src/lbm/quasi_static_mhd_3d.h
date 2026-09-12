#ifndef CGLBM_LBM_QUASI_STATIC_MHD_3D_H
#define CGLBM_LBM_QUASI_STATIC_MHD_3D_H

#include <cstddef>
#include <vector>

#include "lbm/field3d.h"

/// Inductionless magnetohydrodynamics: the low-`Rm` limit, with no induced
/// field ever formed.
///
/// A conducting fluid moving through an imposed field `B_0` separates charge;
/// the resulting current closes through the fluid and pushes back on it. The
/// full description evolves the magnetic field alongside the flow. At the
/// magnetic Reynolds numbers of liquid metals -- `Rm = u L / eta`, below 1e-2
/// in almost every industrial and laboratory flow -- the field the current
/// induces is negligible beside `B_0`, and the induction equation degenerates
/// into an instantaneous constraint. What is left is
///
///     J = sigma (-grad phi + u x B_0),     div J = 0,     F = J x B_0,
///
/// three statements with no time derivative in them: the current adjusts to the
/// velocity within a step. Eliminating `J` between the first two gives the
/// potential equation this class solves,
///
///     div (sigma grad phi) = div (sigma (u x B_0)),
///
/// after which `J` follows from Ohm's law and `F` from the cross product. There
/// is no magnetic field to store, no induction equation to march, and no
/// magnetic diffusion time step to respect -- which is the entire reason the
/// approximation is used, because that time scale is `1e-7` of the flow's.
///
/// # Why the current is built on the faces
///
/// The tempting implementation evaluates `grad phi` at a node, forms
/// `J = sigma(-grad phi + u x B_0)` there, and takes the cross product. It is
/// wrong in a way that gets worse exactly where the method is worth using.
/// `div J = 0` is not a diagnostic of this system; it is the equation that
/// determines `phi`. A node-centred gradient is not the operator the potential
/// was solved against, so the current it produces carries a spurious divergence
/// -- a fictitious charge source -- and the force that divergence generates
/// grows with the Hartmann number, which is to say with the field strength that
/// motivated the calculation.
///
/// So the current lives on the faces, and is built from the same difference the
/// potential equation was discretised with:
///
///     J_f = sigma_f (e_f - (phi_+ - phi_0)),    e_f = (u x B_0) . n_f,
///
/// which makes the discrete `sum_faces J_f` identically the residual of the
/// linear solve, and therefore zero to the tolerance it converged to rather
/// than to the truncation error of some other stencil. Only for the force is
/// the face current averaged back to the node. `charge_imbalance()` reports
/// what is left, and the unit test holds it at the solver tolerance.
///
/// This is the lattice-Boltzmann reading of the conservative scheme of Ni et
/// al., who make the same argument for finite volumes and show that the
/// non-conservative alternative loses accuracy at high Hartmann number.
///
/// # Two ways to solve for the potential
///
/// The potential equation is the whole cost of the coupling, and there are two
/// reasonable ways to discretise it. Both are here, because they fail
/// differently and the choice belongs to the case.
///
/// `PotentialSolver::FiniteVolume` is the one described above: the
/// conservative face operator, solved by a Jacobi-preconditioned conjugate
/// gradient warm-started from the previous step. It is the accurate one. The
/// face current it produces is conservative to the tolerance the solve reached,
/// the conductivity jump is carried exactly by the face coefficient, and the
/// iteration count does not grow when the conductivity ratio does. It is also
/// the one that needs global reductions, three per iteration.
///
/// `PotentialSolver::LatticeBoltzmann` marches the pseudo-transient problem
///
///     d(phi)/dt = div (sigma grad phi) - div (sigma (u x B_0))
///
/// to its steady state on a D3Q7 lattice, which is the Poisson model of Chai
/// and Shi. It is entirely local -- a collision and a stream, the same shape as
/// the flow solver beside it, with the insulating wall falling out of the same
/// half-way bounce-back -- so it carries no reductions at all and parallelises
/// exactly as the rest of the time loop does.
///
/// Four things about it are worth stating plainly, because they are what the
/// choice turns on. All four are measured by
/// `programs/unit_testing/lbm/mhd_potential`; the numbers below are that
/// program's output, not estimates.
///
///  - *It converges, at second order, and it converges slowly.* Against an
///    analytic potential at a uniform conductivity it reaches order 1.98 and
///    2.00 under refinement -- the same order as the finite volume, at about
///    2.5 times its error on the same lattice. The cost is the diffusive one:
///    1152, 4272 and 15856 sweeps from cold at `n` = 16, 32 and 64, which is
///    `O(L^2)` to within a per cent, against a single conjugate-gradient
///    iteration for the same problem. Inside the time loop that gap narrows a
///    long way, because the potential is warm-started and has to move very
///    little; from cold it is three to four orders of magnitude.
///  - Only its *steady state* is wanted, and that is unchanged if `sigma` and
///    the source are scaled together. The march therefore runs at a scaled
///    diffusivity chosen to put the largest relaxation time at
///    `MhdPhysics::lbm_tau_max`, which is a free speed-up of the ratio between
///    the largest conductivity and what a relaxation time of one would have
///    allowed.
///  - The collision is two-relaxation-time, not BGK, with the magic parameter
///    `Lambda = (tau^+ - 1/2)(tau^- - 1/2)` held at 1/4. With BGK the effective
///    position of a discontinuity in `sigma` drifts with the relaxation time,
///    and here `sigma` is discontinuous at every interface and the relaxation
///    time spans the conductivity ratio. Holding `Lambda` fixed is what makes
///    the answer independent of that.
///  - *A large conductivity ratio defeats it.* This is the reason the finite
///    volume is the default, and it is not a matter of tuning. The march runs
///    at a scaled diffusivity, and whatever the scale, the ratio between the
///    fastest and slowest relaxing region is the conductivity ratio itself; the
///    poorly conducting fluid then needs that factor more sweeps than the well
///    conducting one. Measured at a ratio of 10^4 -- a liquid metal beside an
///    electrolyte, which is an ordinary pairing rather than an extreme one --
///    it had not converged after 40000 sweeps, and the current built from the
///    unconverged potential carried a charge imbalance of 45 % of itself. The
///    conjugate gradient solved the same problem in 38 iterations to 1e-12.
///    The elliptic problem's condition number *is* that ratio, and no purely
///    local relaxation escapes it.
///
/// Its steady state is not the finite-volume one -- two discretisations of the
/// same equation, agreeing at second order and not before; measured, 2.6 %
/// apart on a swirling field at `n` = 24. The face current is built the same
/// way from either potential, so `charge_imbalance()` is comparable between
/// them, and it is the honest way to see the difference: against the
/// finite-volume potential it is the linear-solver tolerance -- measured
/// 1.2e-12 of the current itself, periodic, and 2.3e-12 against an insulating
/// wall -- and against the lattice-Boltzmann one it is that scheme's own error,
/// which is small only once the march has actually converged.
///
/// # Boundaries
///
/// The walls are electrically insulating: `J . n = 0`. In flux form that is
/// simply a face the current cannot cross, so the wall face is dropped from the
/// operator and from the right-hand side alike -- no ghost node, no one-sided
/// difference, and the conservation statement above still holds at the wall.
/// A perfectly conducting wall would instead pin `phi`, and is not offered.
///
/// With every boundary insulating or periodic the potential is determined only
/// up to a constant, and the right-hand side is a discrete divergence and so
/// orthogonal to that constant. The solve projects it out each iteration and
/// returns the zero-mean potential.
///
/// # The conductivity of a mixture, and why neither average is right
///
/// The two components carry their own conductivities and the interface here is
/// diffuse, so a face between two nodes of different composition needs a value
/// and there is no single correct one. A laminate has an *anisotropic*
/// effective conductivity: harmonic across the layers, arithmetic along them.
/// A scalar face coefficient cannot carry both, so the choice has to be made
/// against the direction the current actually runs at the interface.
///
///  - Current crossing the interface -- a conducting drop in a field, say --
///    sees the two phases in series, and harmonic is right. Arithmetic would
///    let the current through a resistive shell as though it were not there.
///  - Current running *along* the interface -- which is what a horizontal
///    layer in a vertical field produces -- sees them in parallel, and
///    arithmetic is right. Harmonic then makes the interface band very nearly
///    an insulator: at a conductivity ratio of 10^4 and a face half in each
///    phase, the harmonic value is 2700 times smaller than the arithmetic one,
///    so two or three nodes of interface stop carrying current at all.
///
/// Harmonic is the default, because it is the conservative choice for a current
/// that has to cross and because it is what the reference implementation this
/// module was written against uses. It is the wrong one for a flow like
/// `magnetic_rayleigh_taylor`, and `--conductivity=arithmetic` is there for
/// that. A sharp-interface treatment -- imposing `[phi] = 0` and
/// `[sigma d_n phi] = [sigma (u x B).n]` on the reconstructed interface rather
/// than blending at all -- is what removes the choice, and is not implemented
/// here.
///
/// # The time step
///
/// The force is evaluated at the velocity the populations carry into the step,
/// which makes the magnetic damping explicit. Damping alone relaxes the
/// velocity on `tau_m = rho / (sigma B_0^2)`, so a step much beyond that is
/// unstable however well the potential is solved.
/// :func:`magnetic_damping_time` returns it, and the cases that use a field
/// check their step against it.
///
/// References
///  - M.-J. Ni, R. Munipalli, N. B. Morley, P. Huang, M. A. Abdou, "A current
///    density conservative scheme for incompressible MHD flows at a low
///    magnetic Reynolds number. Part I: On a rectangular collocated grid
///    system", J. Comput. Phys. 227(1), 174 (2007),
///    doi:10.1016/j.jcp.2007.07.025. The conservative face current and why the
///    node-centred one fails at high Hartmann number.
///  - S. Smolentsev, R. Moreau, L. Buhler, C. Mistrangelo, "MHD thermofluid
///    issues of liquid-metal blankets", Fusion Eng. Des. 85, 1196 (2010). The
///    inductionless approximation and its range of validity.
///  - P. A. Davidson, *An Introduction to Magnetohydrodynamics*, Cambridge
///    (2001), ch. 5. The low-`Rm` limit and the damping time `rho / (sigma
///    B^2)`.
///  - Z. Chai, B. Shi, "A novel lattice Boltzmann model for the Poisson
///    equation", Applied Mathematical Modelling 32(10), 2050 (2008),
///    doi:10.1016/j.apm.2007.06.033. The D3Q7 Poisson model.
///  - I. Ginzburg, "Equilibrium-type and link-type lattice Boltzmann models for
///    generic advection and anisotropic-dispersion equation", Adv. Water
///    Resour. 28, 1171 (2005). The two-relaxation-time collision and why the
///    magic parameter is what fixes the placement of a coefficient jump.

namespace cglbm {
namespace lbm {

/// How the potential equation is discretised. See the header.
enum class PotentialSolver {
    /// Conservative face operator, preconditioned conjugate gradient.
    FiniteVolume,
    /// D3Q7 two-relaxation-time pseudo-time march, Chai and Shi.
    LatticeBoltzmann
};

/// Parse "fv"/"finite-volume" or "lbm"/"lattice-boltzmann".
bool potential_solver_from_name(const char* name, PotentialSolver* solver);

/// Name of `solver`, as accepted by :func:`potential_solver_from_name`.
const char* potential_solver_name(PotentialSolver solver);

/// The imposed field, the two conductivities, and how hard to solve.
struct MhdPhysics {
    /// Run the magnetic coupling at all. Off leaves every case untouched.
    bool enabled = false;

    /// The imposed uniform field `B_0`, in lattice units.
    double b[3] = {0.0, 0.0, 0.0};

    /// Electrical conductivity of component 1 (phi = +1).
    double conductivity1 = 1.0;

    /// Electrical conductivity of component 2 (phi = -1).
    double conductivity2 = 1.0;

    /// Average the two conductivities harmonically on a face, not
    /// arithmetically. See the header: series, not parallel.
    bool harmonic_conductivity = true;

    /// Which discretisation of the potential equation to use.
    PotentialSolver solver = PotentialSolver::FiniteVolume;

    /// Relative residual, or relative change per sweep, the solve stops at.
    double tolerance = 1e-10;

    /// Iteration or sweep limit of the potential solve.
    ///
    /// Reaching it is reported, not thrown: a case that drifts out of tolerance
    /// should say so and carry on rather than lose the run. The
    /// lattice-Boltzmann march needs far more of these than the conjugate
    /// gradient does, and its default reflects that.
    int max_iterations = 500;

    /// Largest relaxation time of the lattice-Boltzmann march.
    ///
    /// The steady state does not depend on it -- scaling the conductivity and
    /// the source together leaves the answer alone -- so it is purely how fast
    /// the march runs, and one is the usual compromise between the sweep count
    /// and the accuracy of the most conducting region. Ignored by the
    /// finite-volume solver.
    double lbm_tau_max = 1.0;

    /// Sweeps between two convergence checks of the lattice-Boltzmann march.
    ///
    /// The check costs a pass over the lattice, which is a third of a sweep, so
    /// testing every sweep would be a third of the run spent asking whether it
    /// was over. Ignored by the finite-volume solver.
    int lbm_check_interval = 16;
};

/// The velocity at which magnetic damping acts, `rho / (sigma B^2)`.
///
/// Returns infinity for a vanishing field or conductivity, so a case can
/// compare against it unconditionally.
double magnetic_damping_time(double density, double conductivity, const double* b);

/// The potential, the current and the Lorentz force of one lattice.
class QuasiStaticMhd3D {
public:
    /// `wall_y` makes the two y boundaries insulating walls; otherwise every
    /// axis is periodic. `parallel` runs the loops across OpenMP threads, which
    /// does not change the result: the reductions sum a fixed partition of the
    /// lattice in a fixed order, so the answer does not depend on the thread
    /// count.
    QuasiStaticMhd3D(int nx, int ny, int nz, bool wall_y, MhdPhysics physics, bool parallel);

    /// Solve for the potential and rebuild the current from `velocity` and the
    /// composition `phase`, the bulk-normalised phase field in [-1, 1].
    ///
    /// The potential is kept between calls and used as the starting guess, so
    /// the solve after the first costs a handful of iterations rather than a
    /// hundred.
    void solve(const Field3D& velocity, const Field3D& phase);

    /// Add `J x B_0` to `force`, a `depth == 3` field.
    void add_lorentz_force(Field3D& force) const;

    const Field3D& potential() const {
        return phi_;
    }

    /// The node-centred current, `depth == 3`. Averaged from the faces.
    const Field3D& current() const {
        return current_;
    }

    /// Iterations the last solve took.
    int iterations() const {
        return iterations_;
    }

    /// Relative residual the last solve reached.
    double residual() const {
        return residual_;
    }

    /// Whether the last solve met its tolerance.
    bool converged() const {
        return converged_;
    }

    /// Largest `|sum_faces J_f|` over the lattice after the last solve.
    ///
    /// The discrete charge that the scheme fails to conserve. It is the
    /// residual of the linear solve and nothing else, so it tracks
    /// `tolerance`; anything larger means the face current and the operator
    /// have come apart.
    double charge_imbalance() const;

    /// Apply the potential operator to `in`, writing `A in` into `out`.
    ///
    /// Exposed so a test can manufacture a right-hand side from a potential it
    /// chose, rather than compare against a continuum solution the discrete
    /// operator was never going to reproduce exactly.
    ///
    /// The operator is built from the face conductivities, which `solve()` is
    /// what fills in; call it once with the composition wanted before using
    /// this, or the operator is the zero one.
    void apply_operator(const Field3D& in, Field3D& out) const;

    /// Solve `A phi = rhs` directly, for the same test. Uses the finite-volume
    /// path whatever `MhdPhysics::solver` says, since it is that operator's
    /// inverse that is being asked for.
    void solve_potential(const Field3D& rhs);

    /// The conductivity on the `+axis` face of node (i, j, k).
    double face_conductivity(int i, int j, int k, int axis) const {
        return sigma_face_(i, j, k, axis);
    }

private:
    /// Fill `sigma_face_` from the phase field, and zero the wall faces.
    void update_conductivity(const Field3D& phase);

    /// Fill `drive_` with `sigma_f (u x B_0) . n_f` and `rhs_` with its
    /// negated divergence.
    void update_drive(const Field3D& velocity);

    /// Solve `A phi = rhs_` for the potential, warm-started from `phi_`.
    void solve_flat();

    /// The same, marched to steady state on the D3Q7 lattice.
    void solve_lattice_boltzmann();

    /// Rebuild `current_` from `drive_`, `sigma_face_` and the potential.
    void update_current();

    /// Sum of `a * b` over the lattice, in an order that does not depend on
    /// the thread count.
    double dot(const std::vector<double>& a, const std::vector<double>& b) const;

    /// Subtract the mean of `v`, which is the projection off the nullspace the
    /// all-insulating potential problem has.
    void remove_mean(std::vector<double>& v) const;

    int nx_;
    int ny_;
    int nz_;
    bool wall_y_;
    bool parallel_;
    bool singular_;  ///< every boundary insulating or periodic: phi is defined up to a constant
    MhdPhysics physics_;
    std::size_t nodes_;
    std::size_t chunk_count_;

    Field3D phi_;
    Field3D sigma_face_;
    Field3D sigma_node_;
    Field3D g_;
    Field3D g_next_;
    Field3D phi_previous_;
    Field3D drive_;
    Field3D emf_;
    Field3D face_current_;
    Field3D current_;

    std::vector<double> rhs_;
    std::vector<double> solution_;
    std::vector<double> residual_vector_;
    std::vector<double> direction_;
    std::vector<double> operator_direction_;
    std::vector<double> preconditioned_;
    std::vector<double> diagonal_;
    mutable std::vector<double> partials_;

    int iterations_ = 0;
    double residual_ = 0.0;
    bool converged_ = true;
};

}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_QUASI_STATIC_MHD_3D_H
