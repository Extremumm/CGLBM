#ifndef CGLBM_LBM_VELOCITY_BASED_SOLVER_H
#define CGLBM_LBM_VELOCITY_BASED_SOLVER_H

#include <functional>
#include <vector>

#include "lbm/isotropic_gradient.h"

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
    /// Weight of the populations' own non-equilibrium in collide_hybrid. 0.98
    /// lets a droplet launched at 0.05 lattice units per step diverge after
    /// 7700 steps at a density ratio of 1e4; at 0.7 it holds up to 0.2.
    double hybrid_weight = 0.7;
    /// Stencil of the gradients: colour field, normals and capillary stress.
    GradientStencil stencil = GradientStencil::E8;
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
    std::vector<double> stress_xx_;
    std::vector<double> stress_xy_;
    std::vector<double> stress_yy_;
};

}  // namespace velocity_based
}  // namespace lbm
}  // namespace cglbm

#endif  // CGLBM_LBM_VELOCITY_BASED_SOLVER_H
