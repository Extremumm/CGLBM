/// A capillary wave on a heavy layer, against the exact viscous normal mode.
///
/// A band of component 1, the heavy one, fills ny/4 < y < 3ny/4 of a doubly
/// periodic box; component 2 fills the rest. The lower interface is displaced
/// by `kAmplitude cos(2 pi x / nx)` and released; the upper one starts flat.
/// Each fluid layer is one wavelength deep, so the wave sees two half-spaces to
/// within `exp(-2 pi)`. Surface tension pulls the interface back, and it rings
/// down at the frequency and the rate of the normal mode of two viscous fluids
/// (`pycglbm.normal_modes.capillary_wave`).
///
/// What the case measures is the dynamics of a moving interface at a large
/// density ratio, which no static case can see. The flow inside the heavy
/// layer is the potential flow `exp(k y)`: it stretches the heavy fluid, so its
/// damping depends on the heavy fluid's extensional viscosity, which
/// `--source-stencil` decides (docs/numerics.md, "The heavy fluid's extensional
/// viscosity").
///
/// Beyond the grids, `mode.csv` records every `kTrackInterval` steps the first
/// Fourier coefficient of the lower interface's height, read off the heavy
/// fluid's volume in each column of the lower half.
///
/// The defaults are the density ratio of 1000 the tests run; `--rho1`, `--nu`
/// and `--nu2` move it. The dynamic viscosities are mixed as rho nu and the
/// tension is the divergence of the capillary stress, as the high-ratio
/// Laplace case runs.

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

#include "lbm/case_config.h"
#include "lbm/output_writer.h"
#include "lbm/solver.h"

namespace {

using cglbm::lbm::CaseConfig;

// Conversion factors from lattice to physical units; everything below is in
// lattice units.
const double c_dx = 1.e-5;                        // m
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

constexpr double kPi = 3.14159265358979323846;

/// Initial displacement of the lower interface, in nodes. Small enough that
/// the wave is linear at a density ratio of 1000, where the light side of a
/// moving interface is what limits the scheme.
constexpr double kAmplitude = 0.3;

/// How often `mode.csv` is written.
constexpr int kTrackInterval = 50;

CaseConfig capillary_wave_case() {
    CaseConfig config;
    config.name = "capillary_wave";
    config.nx = 64;
    config.ny = 128;
    config.steps = 25000;
    config.interval = 5000;

    config.physics.rho1 = 1000.;               // kg m^-3, the heavy layer
    config.physics.rho2 = 1.;                  // kg m^-3
    config.physics.c1 = 347. / (c_dx / c_dt);  // speed of sound, lattice units
    config.physics.c2 = 347. / (c_dx / c_dt);
    config.physics.sigma = 1. / (c_dx * c_dx * c_dx / c_dt / c_dt);
    // Dynamic viscosities of 2 in the heavy layer and 0.05 in the light one,
    // in lattice units: tau is 6.5 and 0.65.
    config.physics.nu = 2.e-3;
    config.physics.nu_b = 2.e-3;
    config.physics.nu2 = 0.05;
    config.physics.nu_b2 = 0.05;
    config.physics.ch_width_init = 1.1 * config.units.dx;
    config.physics.ch_width_ope = 1.6 * config.units.dx;
    // Flat interfaces carry no Laplace jump.
    config.physics.radius = std::numeric_limits<double>::infinity();
    config.physics.p2_inf = 0.;
    config.physics.p1_inf = cglbm::lbm::matched_p1_inf(config.physics);
    config.matched_pressure_offset = true;

    config.boundary = cglbm::lbm::Boundary::PeriodicY;
    config.stencil = cglbm::lbm::GradientStencil::E8;
    config.interface_field = cglbm::lbm::InterfaceField::BulkNormalised;
    config.surface_tension = cglbm::lbm::SurfaceTension::CapillaryStress;
    config.viscosity_mixing = cglbm::lbm::ViscosityMixing::Dynamic;

    config.initial_phase = [](const CaseConfig& c, int i, int j) {
        const double width = c.physics.ch_width_init;
        const double lower = c.ny / 4.0 + kAmplitude * std::cos(2.0 * kPi * i / c.nx);
        const double upper = 3.0 * c.ny / 4.0;
        return j < c.ny / 2 ? std::tanh((j - lower) / width) : std::tanh((upper - j) / width);
    };
    return config;
}

/// First cosine coefficient of the lower interface's height.
///
/// The height of column i is ny/2 minus the heavy fluid's volume in the lower
/// half, the volume fraction being (rho - rho2) / (rho1 - rho2).
double mode_amplitude(const cglbm::lbm::Solver& solver) {
    const CaseConfig& c = solver.config();
    const auto& rho = solver.density();
    double amplitude = 0.0;
    for (int i = 0; i < c.nx; ++i) {
        double volume = 0.0;
        for (int j = 0; j < c.ny / 2; ++j) {
            volume += (rho(i, j) - c.physics.rho2) / (c.physics.rho1 - c.physics.rho2);
        }
        amplitude += 2.0 / c.nx * (c.ny / 2.0 - volume) * std::cos(2.0 * kPi * i / c.nx);
    }
    return amplitude;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = capillary_wave_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "capillary_wave")) {
    case cglbm::lbm::CommandLineResult::Finished:
        return 0;
    case cglbm::lbm::CommandLineResult::Error:
        return 2;
    case cglbm::lbm::CommandLineResult::Run:
        break;
    }

    std::cout << cglbm::lbm::describe(config) << "\namplitude = " << kAmplitude << std::endl;
    try {
        cglbm::lbm::Solver solver(config);
        cglbm::lbm::CsvWriter writer(config.output_precision);
        std::ofstream track("mode.csv");
        if (!track) {
            throw cglbm::lbm::OutputError("cannot open mode.csv");
        }
        track.precision(12);
        track << "timestep,amplitude\n";

        solver.initialize();
        writer.write_grids(0, solver.state());
        for (int timestep = 0; timestep <= config.steps; ++timestep) {
            if (timestep > 0) {
                solver.step();
            }
            if (timestep % kTrackInterval == 0) {
                track << timestep << "," << mode_amplitude(solver) << "\n";
            }
            if (timestep > 0 && timestep % config.interval == 0) {
                std::cout << "Step " << timestep << std::endl;
                writer.write_grids(timestep, solver.state());
            }
        }
        if (!track) {
            throw cglbm::lbm::OutputError("cannot write mode.csv");
        }
    } catch (const std::exception& error) {
        std::cerr << "capillary_wave: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
