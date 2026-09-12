/// Magnetic Rayleigh-Taylor: a single mode, grown against a magnetic field.
///
/// The case that exercises the coupling as a whole -- two fluids of different
/// density *and* different conductivity, a moving interface, and a growth rate
/// that the field is supposed to change by a factor of three. Hartmann flow
/// says the force is the right size in a steady one-dimensional flow; this says
/// it is still the right size when the flow is unsteady, two-phase, and the
/// conductivity jumps by four orders of magnitude across the interface.
///
/// What it does *not* exercise, and this is a property of the configuration
/// rather than an omission, is the electric potential. With a two-dimensional
/// perturbation in the x-y plane and `B_0` along y, the motional field
/// `u x B_0` points along z and nothing in the problem varies along z, so its
/// divergence is identically zero and `phi` is uniform. That is true of the
/// analytic problem, which is why the dispersion relation below has a closed
/// form at all, and the run reports the potential as a check that the drive
/// assembly does not invent structure the physics does not have. The potential
/// solve itself is covered by `programs/unit_testing/lbm/mhd_potential`, and
/// engaged in a running flow by `--bz`, which tilts the field into the plane of
/// the interface and leaves the potential with real work to do -- at the cost
/// of the closed form.
///
/// A heavy layer rests on a light one between insulating walls, the interface
/// carries one cosine along x, and gravity pulls the arrangement apart. The
/// field `B_0` is normal to the interface, along y. As the interface grows the
/// fluid on either side moves across `B_0`, a current closes through both
/// layers, and the Lorentz force opposes the motion: the instability still
/// grows, but slower, and the stronger the field the slower it goes.
///
/// # The reference
///
/// The configuration, the parameters and the dispersion relation are those of
/// `basiliskMHD/QSMHD/mrt`, whose growth rates come from the finite-depth
/// quasi-static relation
///
///     rho1^.5 w^1.5 sqrt(w rho1 + B^2 s1) coth(k h1 / sqrt(1 + s1 B^2/(rho1 w)))
///   + rho2^.5 w^1.5 sqrt(w rho2 + B^2 s2) coth(k h2 / sqrt(1 + s2 B^2/(rho2 w)))
///   - [(rho1 - rho2) g k - gamma k^3]  =  0,
///
/// solved for its positive root. Two properties of it decide how this case is
/// built. It is **inviscid** -- no viscosity appears -- so the comparison is
/// only as good as the flow is inviscid, and it is **dimensionally
/// consistent**, so it may be evaluated directly in lattice units.
///
/// The second point is what makes the case possible, and the first is why it
/// has to be. The physical pairing behind the `mrt` campaign is a liquid metal
/// over a molten salt: a kinematic viscosity of 6.8e-7 m^2/s in a domain of
/// 0.1 m. Matching that at 128 nodes across the domain puts the lattice
/// viscosity at 1.1e-5 and the relaxation time at 0.50003, which no
/// lattice-Boltzmann collision survives. So this case does not reproduce the
/// `mrt` *numbers*; it reproduces the `mrt` *relation*, evaluated at lattice
/// parameters the scheme can actually carry, with the density ratio and the
/// conductivity ratio -- the two that set the character of the problem --
/// carried over unchanged at 3.6501 and 10825.
///
/// What remains of viscosity is a correction of about 7 %: `nu k^2 / omega` at
/// the light fluid's viscosity. That is too large to score an absolute growth
/// rate against and small enough that it very largely cancels from a *ratio* of
/// growth rates, which is why the measurement the tests make is
/// `omega(B) / omega(0)` against the same ratio from the relation. The absolute
/// rate at zero field is reported beside it, and is expected to sit a few per
/// cent low for exactly this reason.
///
/// # The field strengths
///
/// The magnetic suppression is governed by `N = sigma B^2 / (rho omega)`, the
/// ratio of the magnetic damping rate to the growth rate. The `mrt` campaign
/// spans `N` from 0.3 at 0.2 T to 225 at 2 T; this case takes `N = 0`, 2 and
/// 10 in the heavy fluid, which covers the transition from barely affected to
/// strongly suppressed. `--sigma-e1` selects it, and the magnetic damping time
/// `rho / (sigma B^2)` stays above 96 steps at the strongest, so the explicit
/// force is nowhere near its stability bound.
///
/// The conductivity ratio is 10825, which is the one number in this case that
/// the potential solve genuinely has to work for. It is also the reason the
/// finite-volume solve is the default and the lattice-Boltzmann march is not:
/// measured in `programs/unit_testing/lbm/mhd_potential`, the march does not
/// converge at that contrast.

#include <cmath>
#include <iostream>
#include <string>

#include "lbm/case_config.h"
#include "lbm/two_population_solver_3d.h"

namespace {

using cglbm::lbm::CaseConfig;

constexpr double kPi = 3.14159265358979323846;

/// Density ratio of the `mrt` pairing, 6260 / 1715.
const double kDensityRatio = 3.6501457725947524;

/// Conductivity ratio of the same pairing, 8.66e5 / 80.
const double kConductivityRatio = 10825.0;

/// Dynamic viscosity shared by the two fluids, so tau is uniform.
///
/// 0.03 puts tau at 0.579 on D3Q27 at this density ratio, and leaves the
/// viscous correction `nu k^2 / omega` at 6.8 % in the light fluid. Lower would
/// sharpen the comparison and bring tau towards 1/2; this is the compromise.
const double kViscosity = 0.03;

/// Initial displacement of the interface, in nodes.
///
/// Half a node: a 250th of the wavelength, well inside the linear regime the
/// dispersion relation describes, and small enough that the run has four
/// e-folds of growth before the mode leaves it.
const double kSeed = 0.5;

/// A flat interface displaced by one cosine along x, and nothing along z.
///
/// The three-dimensional `cosine_layer_3d` carries a mode along both x and z,
/// whose wavenumber is `sqrt(2)` times either one; the dispersion relation is
/// written for a single wavenumber, so this case takes a mode that has one.
/// The profile is prescribed in the bulk-normalised field, as everywhere else.
cglbm::lbm::PhaseFieldInit3D cosine_layer_x(double displacement) {
    return [displacement](const CaseConfig& config, int i, int j, int) {
        const int y0 = config.ny / 2;
        const double offset = displacement * std::cos(2.0 * kPi * i / config.nx);
        return std::tanh(((j - y0) - offset) / config.physics.ch_width_ope);
    };
}

CaseConfig magnetic_rayleigh_taylor_case() {
    CaseConfig config;
    config.name = "magnetic_rayleigh_taylor";
    config.nx = 128;  // one wavelength
    config.ny = 128;  // wall to wall, interface centred
    config.nz = 4;    // the mode has no structure along z
    config.steps = 20000;
    config.interval = 250;

    config.physics.rho1 = kDensityRatio;  // the heavy fluid, on top
    config.physics.rho2 = 1.0;
    config.physics.gravity = 4.0e-5;
    config.physics.sigma = 1.0e-3;

    config.physics.nu = kViscosity / config.physics.rho1;
    config.physics.nu_b = config.physics.nu;
    config.physics.nu2 = kViscosity / config.physics.rho2;
    config.physics.nu_b2 = config.physics.nu2;

    config.physics.alpha2 = 0.2;
    config.physics.beta = 0.7;
    config.physics.ch_width_init = 1.1;

    // N = sigma B^2 / (rho1 omega) = 10 at the growth rate this case has
    // without a field. `--sigma-e1=0` recovers the hydrodynamic run, which is
    // what the suppression is measured against.
    config.mhd.enabled = true;
    config.mhd.b[1] = 1.0;  // normal to the interface
    config.mhd.conductivity1 = 3.81e-2;
    config.mhd.conductivity2 = config.mhd.conductivity1 / kConductivityRatio;
    config.mhd.tolerance = 1e-10;
    config.mhd.max_iterations = 500;

    // Arithmetic, against the library default. The current here runs *along*
    // the interface -- `u x B_0` points across the mode, in the plane of the
    // layer -- so the two phases are in parallel and arithmetic is the right
    // average; harmonic is right for a current that has to cross. The
    // difference is not academic at a conductivity ratio of 10825: a face half
    // in each phase gets 2.7e3 times more conductivity from the arithmetic
    // average than from the harmonic one, so a harmonic blend turns the two or
    // three nodes of interface into an insulator, exactly where the shear is,
    // and the mode grows very much faster than the relation says it should.
    // See `quasi_static_mhd_3d.h`.
    config.mhd.harmonic_conductivity = false;

    config.boundary = cglbm::lbm::Boundary::WallY;
    config.lattice_3d = cglbm::lbm::Lattice3DKind::D3Q27;
    config.stencil_3d = cglbm::lbm::GradientStencil3D::E4;
    config.initial_phase_3d = cosine_layer_x(kSeed);
    config.parallel = true;

    return config;
}

/// The numbers the growth-rate fit needs, so it restates none of them.
void describe_mode(const CaseConfig& config) {
    const double k = 2.0 * kPi / config.nx;
    std::cout.precision(17);
    std::cout << "wavenumber = " << k << "\n"
              << "layer_depth = " << 0.5 * config.ny << "\n"
              << "seed_amplitude = " << kSeed << "\n"
              << "density_ratio = " << config.physics.rho1 / config.physics.rho2 << "\n"
              << "conductivity_ratio = "
              << (config.mhd.conductivity2 > 0.0
                      ? config.mhd.conductivity1 / config.mhd.conductivity2
                      : 0.0)
              << "\n"
              << "magnetic_damping_time_heavy = "
              << cglbm::lbm::magnetic_damping_time(
                     config.physics.rho1, config.mhd.conductivity1, config.mhd.b)
              << std::endl;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = magnetic_rayleigh_taylor_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "magnetic_rayleigh_taylor")) {
    case cglbm::lbm::CommandLineResult::Finished:
        return 0;
    case cglbm::lbm::CommandLineResult::Error:
        return 2;
    case cglbm::lbm::CommandLineResult::Run:
        break;
    }

    // `--sigma-e1=X` selects the field strength, and `--sigma-e1=0` is the
    // hydrodynamic reference the suppression is measured against. Either way
    // the second conductivity follows the first, so that the pair stays the
    // metal-and-salt contrast the case is about instead of quietly becoming two
    // equally conducting fluids. Saying `--sigma-e2` explicitly overrides that.
    bool second_given = false;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind("--sigma-e2=", 0) == 0) {
            second_given = true;
        }
    }
    if (!second_given) {
        config.mhd.conductivity2 = config.mhd.conductivity1 / kConductivityRatio;
    }

    std::cout << cglbm::lbm::describe(config) << std::endl;
    describe_mode(config);
    try {
        cglbm::lbm::TwoPopulationSolver3D solver(config);
        solver.run();
    } catch (const std::exception& error) {
        std::cerr << "magnetic_rayleigh_taylor: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
