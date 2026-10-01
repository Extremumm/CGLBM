#ifndef CGLBM_LBM_VELOCITY_BASED_SOLVER_H
#define CGLBM_LBM_VELOCITY_BASED_SOLVER_H

#include <functional>
#include <vector>

#include "lbm/isotropic_gradient.h"
#include "lbm/velocity_based.h"

/// The time step of the velocity-based two-phase scheme of velocity_based.h,
/// on a lattice periodic along x and, along y, periodic or closed by resting
/// walls.
///
/// A step collides the hydrodynamic populations, builds the phase
/// populations, streams both, updates the momentum link by link, and applies
/// the forces through the forcing term, half before and half after the step:
/// surface tension, pressure, the link dissipation, gravity and an optional
/// body force. The programs under programs/solvers/velocity_based set up the
/// initial state and write the output; this does everything in between.
///
/// Surface tension is the capillary stress of surface_force.h on psi = 2c - 1,
/// with each layer of the interface weighted by layer_weight so that a
/// circular interface carries sigma / R whatever its width.
///
/// Walls (Boundary::WallY) sit halfway between the first and last rows and
/// their ghosts, at j = -1/2 and j = ny - 1/2. Both sets of populations bounce
/// back there: the fluid is at rest on the wall, and no volume crosses it, nor
/// any of P, the pressure number the populations carry. A link through a wall
/// joins a node to its own point image (wall_image), so
/// the momentum exchange is the lattice's own friction and nothing else. The
/// stencil operators read every field mirrored across the wall
/// (velocity_based.h, field_at), which leaves a uniform phase field without a
/// normal at the wall and meets an interface there at a right angle: the walls
/// are neutrally wetting.

namespace cglbm {
namespace lbm {
namespace velocity_based {

/// How the dynamic viscosity is mixed across the diffuse interface, on the
/// volume fraction c of component 1.
enum class InterfaceViscosity {
    /// mu = c mu1 + (1 - c) mu2.
    Arithmetic,
    /// 1 / mu = c / mu1 + (1 - c) / mu2.
    Harmonic,
    /// The harmonic mean for the strain that shears across the interface and
    /// the arithmetic mean for the strain that stretches along it
    /// (collide_laminate); the links take the harmonic mean.
    Laminate
};

struct SolverParameters {
    int nx = 128;
    int ny = 128;
    /// Densities of component 1 (c = 1) and component 2 (c = 0).
    double rho1 = 1.0;
    double rho2 = 1.0;
    /// Dynamic viscosities of the two components.
    double mu1 = 1.0 / 6.0;
    double mu2 = 1.0 / 6.0;
    double surface_tension = 0.0;
    /// Interface width W of the profile c = (1 + tanh(x / W)) / 2.
    double width = 1.6;
    /// How the two viscosities are mixed across the interface.
    InterfaceViscosity interface_viscosity = InterfaceViscosity::Arithmetic;
    /// Relaxation time of the trace of the non-equilibrium.
    double tau_bulk = 1.0;
    /// Weight of the populations' own non-equilibrium in collide_filtered; the
    /// rest is its mean over this step and the previous one. At 1 a droplet
    /// launched at 0.1 lattice units per step diverges within 1000 steps at a
    /// density ratio of 1e4; at 0.7 it holds.
    ///
    /// The solver used to take the rest from the finite-difference velocity
    /// gradient instead (collide_hybrid, at the same 0.7). That misread the
    /// heavy fluid's oscillatory boundary layer wherever it is thinner than
    /// the interface: a capillary wave at 1e4 was damped 2.10 times the exact
    /// rate, and 1.27 times with the mean in time (before fourth_order_phase).
    double filter_weight = 0.7;
    /// Stencil of the gradients: colour field, normals and capillary stress.
    /// With fourth_order_phase the phase populations sharpen along a normal of
    /// their own, on the lattice's stencil whatever this one is.
    ///
    /// E8 is the best of them here. On the static droplet of radius 10 at
    /// 1e4, E10 and E12 raise the currents from 2.8e-6 to 4.6e-6 and 7.6e-6
    /// and the mode-4 deformation from 6e-5 to 1.0e-4 and 1.4e-4.
    GradientStencil stencil = GradientStencil::E8;
    /// Lattice temperature of the phase populations' carrier, which sets the
    /// interface mobility M = phase_temperature / 2 (phase_carrier).
    ///
    /// 0.2 rather than the lattice's own cs^2 = 1/3. It took the capillary
    /// wave at a density ratio of 1000 from 1.084 to 1.049 times the exact
    /// damping rate (with the hybrid collision then in use) and left the
    /// static, moving and fast droplets and the sheared layers where they
    /// were. What it took out was not the mobility but the surface diffusion
    /// fourth_order_phase describes, which is proportional to T: 0.088 of the
    /// exact rate at cs^2, 0.053 at 0.2. 0.2 is the lowest
    /// round value at which the carrier stays non-negative up to |u| = 0.2,
    /// the fastest a droplet has been launched in a research copy (that needs
    /// 0.184; the 0.1 of the long tests needs 0.106). Must lie below 0.6,
    /// where the rest weight is still positive.
    double phase_temperature = 0.2;
    /// Build the phase populations on the fourth-order normal of
    /// interface_normal, within kFourthOrderNormalBand of the interface, and
    /// with the flux of phase_correction_flux. capillary_wave_vb and
    /// oscillation_vb take it as --fourth-order-phase.
    ///
    /// Without them the phase field alone, the fluid held at rest, relaxes a
    /// mode-2 droplet of radius 20 at 2.9e-6 per step with the E8 normal:
    /// surface diffusion, proportional to phase_temperature and to R^-4, three
    /// quarters of it from the normal's truncation error and the rest from the
    /// transport's own fourth-order error. With them it relaxes at 1.1e-7. At
    /// a density ratio of 1e4 the droplet's own viscous damping is 2e-6 per
    /// step and half the shape relaxation adds to it: the droplet is damped
    /// 1.76 times the exact rate, and 1.07 times with them.
    ///
    /// Off by default. The same surface diffusion had been offsetting the
    /// under-resolved boundary layer of the capillary waves on a wavelength of
    /// 64, which come to 0.976, 0.931 and 0.719 of the exact rate at 100, 1000
    /// and 1e4 with it on; on a wavelength of 128 the wave at 1000 comes to
    /// 0.976, and at 1e4 with the layer resolved by a heavier viscosity to
    /// 0.992. The static droplet at 1e4 keeps ringing in its mode 4, which the
    /// surface diffusion had damped, with currents swinging up to 1.5e-5
    /// instead of 3e-6.
    /// See docs/numerics.md, "The phase field's own surface diffusion".
    bool fourth_order_phase = false;
    /// fourth_order_phase carried a further two orders: the phase populations
    /// sharpen along sixth_order_normal, within the same band, and carry the
    /// flux of sixth_order_flux. It implies fourth_order_phase.
    ///
    /// It pays little. The droplet at 1e4 goes from 1.069 to 1.061 times the
    /// exact damping, the waves on a wavelength of 64 stay within 0.003. The
    /// fourth-order correction leaves the phase field's own relaxation of a
    /// mode-2 droplet of radius 20 at 1.1e-7 per step, and this one at
    /// 1.0e-7; at radius 10 it takes it from 3.3e-7 to 1.5e-6, and on an
    /// interface twice as wide (W = 3.2) from 7.5e-8 to 9.3e-8. What the
    /// fourth-order correction leaves is no longer the truncation error of the
    /// lattice operators, which this one takes to eighth order, but the band's
    /// edge and the sharpening's own nonlinearity. See docs/numerics.md,
    /// "The phase field's own surface diffusion".
    bool sixth_order_phase = false;
    /// How the lattice is closed along y: periodic, or by resting walls.
    /// Periodic along x either way.
    Boundary boundary = Boundary::PeriodicY;
    /// Gravitational acceleration. The force per unit volume is
    /// (rho - gravity_reference_density) g, which follows the interface as it
    /// moves. The reference changes only the pressure, by rho_ref g . x, in an
    /// incompressible flow; between walls any value holds, and on a lattice
    /// periodic along g only the mean density leaves no net force.
    double gravity_x = 0.0;
    double gravity_y = 0.0;
    double gravity_reference_density = 0.0;
    /// Run the per-node loops across OpenMP threads.
    ///
    /// Every loop of a step writes only its own node, or in the streaming a
    /// slot that exactly one node writes, and none sums across nodes, so the
    /// fields are the same to the last bit whatever the thread count. Off by
    /// default, as in the colour-gradient solvers; the programs turn it on
    /// with --threads=N.
    bool parallel = false;
};

/// The state of one node, to initialise from.
struct NodeState {
    /// Volume fraction of component 1.
    double c = 0.0;
    double ux = 0.0;
    double uy = 0.0;
    /// Pressure relative to the far field.
    double p = 0.0;
};

class Solver {
public:
    explicit Solver(const SolverParameters& parameters);

    /// Equilibrium populations for the state `state(i, j)` gives every node.
    void initialize(const std::function<NodeState(int i, int j)>& state);

    /// A body force per unit volume, constant in time, applied from the next
    /// step on. Zero unless set.
    void set_body_force(const std::function<void(int i, int j, double* fx, double* fy)>& force);

    void step();

    int nx() const {
        return parameters_.nx;
    }
    int ny() const {
        return parameters_.ny;
    }
    const SolverParameters& parameters() const {
        return parameters_;
    }

    double volume_fraction(int i, int j) const {
        return c_[node(i, j)];
    }
    double density(int i, int j) const {
        return rho_[node(i, j)];
    }
    /// psi = 2c - 1 with c bounded to [0, 1]: +1 in component 1.
    double phase(int i, int j) const {
        return psi_[node(i, j)];
    }
    /// Pressure relative to the far field, p = rho cs^2 P.
    double pressure(int i, int j) const;
    double velocity_x(int i, int j) const {
        return ux_[node(i, j)];
    }
    double velocity_y(int i, int j) const {
        return uy_[node(i, j)];
    }

private:
    int node(int i, int j) const {
        return i * parameters_.ny + j;
    }
    /// The viscosity the links and the isotropic collision take, by
    /// interface_viscosity.
    double dynamic_viscosity(double c) const;
    double arithmetic_viscosity(double c) const;
    double harmonic_viscosity(double c) const;

    /// Link k of node (i, j) reaches through a wall: its donor (i, j) - xi_k
    /// lies outside the lattice.
    bool through_wall(int j, int k) const {
        const int jd = j - kVelocity[k][1];
        return walls_ && (jd < 0 || jd >= parameters_.ny);
    }

    void macroscopic();
    void acceleration();
    void momentum();
    void velocity();
    void save_step();
    void collide_and_stream();

    SolverParameters parameters_;
    bool walls_;

    // populations, kQ per node
    std::vector<double> g_;
    std::vector<double> h_;
    std::vector<double> g_new_;
    std::vector<double> h_new_;
    // link dissipation coefficients, kQ per node
    std::vector<double> dissipation_;
    // the non-equilibrium of the previous collision, xx yy xy per node
    std::vector<double> previous_non_equilibrium_;

    std::vector<double> c_;
    std::vector<double> rho_;
    std::vector<double> psi_;
    std::vector<double> pressure_number_;
    std::vector<double> ux_;
    std::vector<double> uy_;
    std::vector<double> ax_;
    std::vector<double> ay_;
    std::vector<double> lattice_ux_;
    std::vector<double> lattice_uy_;
    std::vector<double> body_x_;
    std::vector<double> body_y_;
    // the previous step, which the momentum exchange reads
    std::vector<double> rho_old_;
    std::vector<double> pressure_number_old_;
    std::vector<double> ux_old_;
    std::vector<double> uy_old_;
    std::vector<double> ax_old_;
    std::vector<double> ay_old_;
    // surface tension
    std::vector<double> grad_x_;
    std::vector<double> grad_y_;
    std::vector<double> normal_x_;
    std::vector<double> normal_y_;
    // the phase populations' normal and correction flux, and the lattice
    // Laplacians of psi and c they are taken from
    std::vector<double> laplacian_psi_;
    std::vector<double> laplacian_c_;
    // and L L psi and L L c, for sixth_order_phase
    std::vector<double> laplacian2_psi_;
    std::vector<double> laplacian2_c_;
    std::vector<double> phase_normal_x_;
    std::vector<double> phase_normal_y_;
    std::vector<double> phase_flux_x_;
    std::vector<double> phase_flux_y_;
    std::vector<double> stress_xx_;
    std::vector<double> stress_xy_;
    std::vector<double> stress_yy_;
};

}  // namespace velocity_based
}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_VELOCITY_BASED_SOLVER_H
