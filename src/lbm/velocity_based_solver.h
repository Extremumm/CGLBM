#ifndef CGLBM_LBM_VELOCITY_BASED_SOLVER_H
#define CGLBM_LBM_VELOCITY_BASED_SOLVER_H

#include <functional>
#include <vector>

#include "lbm/isotropic_gradient.h"
#include "lbm/velocity_based.h"

/// The time step of the velocity-based two-phase scheme of velocity_based.h,
/// on a doubly periodic lattice.
///
/// A step collides the hydrodynamic populations, builds the phase
/// populations, streams both, updates the momentum link by link, and applies
/// the forces through the forcing term, half before and half after the step:
/// surface tension, pressure, the link dissipation, and an optional body
/// force. The programs under programs/solvers/velocity_based set up the
/// initial state and write the output; this does everything in between.
///
/// Surface tension is the capillary stress of surface_force.h on psi = 2c - 1,
/// with each layer of the interface weighted by layer_weight so that a
/// circular interface carries sigma / R whatever its width.

namespace cglbm {
namespace lbm {
namespace velocity_based {

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
    double dynamic_viscosity(double c) const;

    void macroscopic();
    void acceleration();
    void momentum();
    void velocity();
    void save_step();
    void collide_and_stream();

    SolverParameters parameters_;

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
