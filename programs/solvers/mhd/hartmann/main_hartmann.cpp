/// Hartmann flow: a conducting fluid driven along a channel across a magnetic
/// field.
///
/// The one configuration of inductionless MHD with a closed form, and the
/// standard thing to check a coupling against. A uniform force `G` drives the
/// flow along x between insulating walls at `y = +-L`, with the field `B_0`
/// normal to them. The current is `J = sigma u B_0` along z, the Lorentz force
/// `-sigma B_0^2 u` opposes the motion, and the balance
///
///     mu u'' - sigma B_0^2 u = -G,     u(+-L) = 0
///
/// gives
///
///     u(y) = (G / (sigma B_0^2)) [1 - cosh(Ha y/L) / cosh(Ha)],
///     Ha = B_0 L sqrt(sigma / mu).
///
/// What that profile does is worth seeing rather than only reading. The
/// parabola of ordinary channel flow is replaced by a flat core, where the
/// pressure gradient is balanced by the magnetic force alone and the velocity
/// is simply `G / (sigma B_0^2)`, and two boundary layers of thickness `L / Ha`
/// in which viscosity brings it to rest. At the `Ha = 10` of this case the core
/// is flat to a part in 10^4 and the layers are 3.2 nodes thick.
///
/// # What this validates, and what it does not
///
/// Everything downstream of the potential: Ohm's law, the face current, its
/// average back to the node, the cross product, the insulating wall, and the
/// magnitude of the force that reaches the populations. A sign error, a factor
/// of two, or a current averaged the wrong way all show up as the wrong core
/// velocity or the wrong layer thickness, and the two are independent -- the
/// core measures `sigma B_0^2` and the layers measure `sqrt(mu / sigma) / B_0`.
///
/// It does *not* exercise the potential solve, and that is a property of the
/// configuration rather than an omission. With `u = u(y) x^` and `B_0` along y,
/// the motional field `u x B_0` points along z and does not vary along z, so
/// its divergence vanishes identically and `phi` is uniform. That is a check in
/// itself -- the solve must return zero here, and the run reports it -- but the
/// potential is exercised by
/// `programs/unit_testing/lbm/mhd_potential` and by
/// `magnetic_rayleigh_taylor`, not here.
///
/// # Why these numbers
///
/// The relaxation time is set to 0.9330, where `(tau - 1/2)^2 = 3/16`. That is
/// the value at which the half-way bounce-back wall sits on the half node
/// rather than a `tau`-dependent distance from it, so the channel is the width
/// the analytic profile is evaluated on. Left at the viscosity a case would
/// otherwise pick, the wall drifts and the comparison measures that drift as
/// much as it measures the physics.
///
/// The field is then chosen for `Ha = 10`: high enough for a core and two
/// distinguishable layers, low enough that 3.2 nodes resolve a layer. The drive
/// puts the core at Mach 0.03. Both relaxation times of the problem -- the
/// magnetic `rho / (sigma B_0^2)` and the viscous `(L/Ha)^2 / nu` -- are about
/// 62 steps, so 20000 is steady several hundred times over.
///
/// The lattice is D3Q27. Nothing here needs it, the flow being one-dimensional
/// and single-phase, and that is the point: it is the case that would show a
/// mistake in the corner velocities as a broken profile rather than as a subtly
/// different droplet.
///
/// Reference
///  - J. Hartmann, "Hg-dynamics I: theory of the laminar flow of an electrically
///    conductive liquid in a homogeneous magnetic field", Det Kgl. Danske
///    Videnskabernes Selskab, Mat.-fys. Medd. 15(6) (1937).
///  - P. A. Davidson, *An Introduction to Magnetohydrodynamics*, Cambridge
///    (2001), sec. 5.4.

#include <cmath>
#include <iostream>

#include "lbm/case_config.h"
#include "lbm/two_population_solver_3d.h"

namespace {

using cglbm::lbm::CaseConfig;

/// Density of the single fluid.
const double kDensity = 1.0;

/// Kinematic viscosity, chosen for `tau = 0.9330`.
///
/// `tau = nu / (c_s^k)^2 + 1/2` with `(c_s^k)^2 = 9(1 - alpha)/19 = 0.378947`
/// at `alpha = 0.2` on D3Q27, so `nu = 0.378947 * sqrt(3)/4`.
const double kViscosity = 0.378947368421052632 * 0.433012701892219323;

/// Electrical conductivity of the fluid.
const double kConductivity = 1.0;

/// Lattice nodes across the channel. The walls sit half a node outside the
/// first and last, so the half-width is `ny / 2`.
const int kAcross = 64;

/// Hartmann number this case is built for.
const double kHartmann = 10.0;

/// Core velocity the drive is chosen for.
const double kCoreVelocity = 0.02;

double half_width() {
    return 0.5 * kAcross;
}

double field_strength() {
    // Ha = B L sqrt(sigma / mu), with mu = rho nu.
    return kHartmann / (half_width() * std::sqrt(kConductivity / (kDensity * kViscosity)));
}

double drive() {
    // u_core = (G / (sigma B^2)) [1 - 1/cosh(Ha)].
    const double b = field_strength();
    return kCoreVelocity * kConductivity * b * b / (1.0 - 1.0 / std::cosh(kHartmann));
}

CaseConfig hartmann_case() {
    CaseConfig config;
    config.name = "hartmann";
    config.nx = 8;  // the flow does not vary along x or z; these are here to
    config.ny = kAcross;
    config.nz = 8;  // keep the case three-dimensional, not to resolve anything
    config.steps = 20000;
    config.interval = 4000;

    // One fluid. Equal densities make the two populations identical, the colour
    // gradient vanishes, and the surface-tension and segregation operators are
    // inert; what is left is a single-component lattice Boltzmann fluid.
    config.physics.rho1 = kDensity;
    config.physics.rho2 = kDensity;
    config.physics.sigma = 0.0;
    config.physics.nu = kViscosity;
    config.physics.nu_b = kViscosity;
    config.physics.nu2 = kViscosity;
    config.physics.nu_b2 = kViscosity;
    config.physics.alpha2 = 0.2;

    config.physics.body_force[0] = drive();

    config.mhd.enabled = true;
    config.mhd.b[1] = field_strength();  // normal to the walls
    config.mhd.conductivity1 = kConductivity;
    config.mhd.conductivity2 = kConductivity;
    config.mhd.tolerance = 1e-12;
    config.mhd.max_iterations = 500;

    config.boundary = cglbm::lbm::Boundary::WallY;
    config.lattice_3d = cglbm::lbm::Lattice3DKind::D3Q27;
    config.stencil_3d = cglbm::lbm::GradientStencil3D::E4;
    config.initial_phase_3d = [](const CaseConfig&, int, int, int) { return 1.0; };
    config.parallel = true;

    return config;
}

/// The case's derived numbers, for the test to read back rather than restate.
void describe_hartmann(const CaseConfig& config) {
    const double b = config.mhd.b[1];
    const double mu = config.physics.rho1 * config.physics.nu;
    const double sigma = config.mhd.conductivity1;
    const double length = half_width();
    const double hartmann = b * length * std::sqrt(sigma / mu);
    std::cout.precision(17);
    std::cout << "hartmann_number = " << hartmann << "\n"
              << "half_width = " << length << "\n"
              << "layer_thickness = " << length / hartmann << "\n"
              << "core_velocity = "
              << config.physics.body_force[0] / (sigma * b * b) *
                     (1.0 - 1.0 / std::cosh(hartmann))
              << "\n"
              << "magnetic_damping_time = "
              << cglbm::lbm::magnetic_damping_time(config.physics.rho1, sigma, config.mhd.b)
              << std::endl;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = hartmann_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "hartmann")) {
    case cglbm::lbm::CommandLineResult::Finished:
        return 0;
    case cglbm::lbm::CommandLineResult::Error:
        return 2;
    case cglbm::lbm::CommandLineResult::Run:
        break;
    }

    std::cout << cglbm::lbm::describe(config) << std::endl;
    describe_hartmann(config);
    try {
        cglbm::lbm::TwoPopulationSolver3D solver(config);
        solver.run();

        // The potential must be uniform in this configuration; anything else is
        // a drive that should not exist. Reported rather than asserted, because
        // the program's job is to produce the run and the test's is to judge it.
        const cglbm::lbm::QuasiStaticMhd3D* mhd = solver.mhd();
        if (mhd != nullptr) {
            double largest = 0.0;
            const cglbm::lbm::Field3D& potential = mhd->potential();
            for (int i = 0; i < potential.nx(); ++i) {
                for (int j = 0; j < potential.ny(); ++j) {
                    for (int k = 0; k < potential.nz(); ++k) {
                        largest = std::max(largest, std::fabs(potential(i, j, k)));
                    }
                }
            }
            std::cout.precision(17);
            std::cout << "final_max_potential = " << largest << "\n"
                      << "final_charge_imbalance = " << mhd->charge_imbalance() << std::endl;
        }
    } catch (const std::exception& error) {
        std::cerr << "hartmann: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
