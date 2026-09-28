/// Mode-2 oscillation of a two-dimensional droplet, against the exact viscous
/// normal mode.
///
/// A droplet of component 1, the heavy one, is laid down as
/// `r(theta) = R' (1 + kDeformation cos 2 theta)` in a doubly periodic box and
/// released; `R'` is shrunk so the droplet holds the area of a circle of
/// `physics.radius`. Surface tension pulls it back to the circle, and it rings
/// down at the frequency and the rate of the normal mode of two viscous fluids
/// (`pycglbm.normal_modes.droplet_mode`). `oscillation_3d` is the same
/// experiment in three dimensions, scored against Lamb's inviscid frequency.
///
/// Inside the droplet the mode-2 flow is a pure linear strain, which no
/// finite difference of the flow gets wrong: what this case measures is the
/// interface in motion. At a density ratio of 1000 that is the part of the
/// colour-gradient model that fails (docs/numerics.md, "The heavy fluid's
/// extensional viscosity").
///
/// Beyond the grids, `mode.csv` records every `kTrackInterval` steps the
/// droplet's deformation, the second moment of its excess density
///
///     D = sum (rho - rho2) (x^2 - y^2) / sum (rho - rho2) (x^2 + y^2),
///
/// about the domain centre, which is `2 eps` for a small deformation `eps`.
///
/// The defaults are the density ratio of 1000 the tests run, with the dynamic
/// viscosities mixed as rho nu and the tension as the divergence of the
/// capillary stress, as the high-ratio Laplace case runs.

#include <cmath>
#include <fstream>
#include <iostream>

#include "lbm/case_config.h"
#include "lbm/output_writer.h"
#include "lbm/solver.h"

namespace {

using cglbm::lbm::CaseConfig;

// Conversion factors from lattice to physical units; everything below is in
// lattice units.
const double c_dx = 1.e-5;                        // m
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

/// Initial mode-2 deformation, as a fraction of the radius. Small enough to be
/// linear: at 0.1 and a density ratio of 1000 the light side of the interface
/// moves fast enough to break the scheme with `--source-stencil=matched`.
constexpr double kDeformation = 0.03;

/// How often `mode.csv` is written.
constexpr int kTrackInterval = 50;

CaseConfig oscillation_case() {
    CaseConfig config;
    config.name = "oscillation";
    config.nx = 128;
    config.ny = 128;
    config.steps = 24000;
    config.interval = 6000;

    config.physics.rho1 = 1000.;               // kg m^-3, the droplet
    config.physics.rho2 = 1.;                  // kg m^-3
    config.physics.c1 = 347. / (c_dx / c_dt);  // speed of sound, lattice units
    config.physics.c2 = 347. / (c_dx / c_dt);
    config.physics.radius = 20.;
    config.physics.sigma = 1. / (c_dx * c_dx * c_dx / c_dt / c_dt);
    // Dynamic viscosities of 2 in the droplet and 0.05 around it, in lattice
    // units: tau is 6.5 and 0.65.
    config.physics.nu = 2.e-3;
    config.physics.nu_b = 2.e-3;
    config.physics.nu2 = 0.05;
    config.physics.nu_b2 = 0.05;
    config.physics.ch_width_init = 1.1 * config.units.dx;
    config.physics.ch_width_ope = 1.6 * config.units.dx;
    config.physics.p2_inf = 0.;
    config.physics.p1_inf = cglbm::lbm::matched_p1_inf(config.physics);
    config.matched_pressure_offset = true;

    config.boundary = cglbm::lbm::Boundary::PeriodicY;
    config.stencil = cglbm::lbm::GradientStencil::E8;
    config.interface_field = cglbm::lbm::InterfaceField::BulkNormalised;
    config.surface_tension = cglbm::lbm::SurfaceTension::CapillaryStress;
    config.viscosity_mixing = cglbm::lbm::ViscosityMixing::Dynamic;

    config.initial_phase = [](const CaseConfig& c, int i, int j) {
        // r = R' (1 + eps cos 2 theta) encloses pi R'^2 (1 + eps^2 / 2).
        const double scale = 1.0 / std::sqrt(1.0 + 0.5 * kDeformation * kDeformation);
        const double x = i - c.nx / 2;
        const double y = j - c.ny / 2;
        const double distance = std::sqrt(x * x + y * y);
        const double cos_2theta = distance > 0.0 ? (x * x - y * y) / (distance * distance) : 0.0;
        const double radius = c.physics.radius * scale * (1.0 + kDeformation * cos_2theta);
        return -std::tanh((distance - radius) / c.physics.ch_width_init);
    };
    return config;
}

/// The deformation D of the module comment.
double deformation(const cglbm::lbm::Solver& solver) {
    const CaseConfig& c = solver.config();
    const auto& rho = solver.density();
    double difference = 0.0;
    double sum = 0.0;
    for (int i = 0; i < c.nx; ++i) {
        for (int j = 0; j < c.ny; ++j) {
            const double x = i - c.nx / 2;
            const double y = j - c.ny / 2;
            const double excess = rho(i, j) - c.physics.rho2;
            difference += excess * (x * x - y * y);
            sum += excess * (x * x + y * y);
        }
    }
    return difference / sum;
}

}  // namespace

int main(int argc, char** argv) {
    CaseConfig config = oscillation_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "oscillation")) {
    case cglbm::lbm::CommandLineResult::Finished:
        return 0;
    case cglbm::lbm::CommandLineResult::Error:
        return 2;
    case cglbm::lbm::CommandLineResult::Run:
        break;
    }

    std::cout << cglbm::lbm::describe(config) << "\ndeformation = " << kDeformation << std::endl;
    try {
        cglbm::lbm::Solver solver(config);
        cglbm::lbm::CsvWriter writer(config.output_precision);
        std::ofstream track("mode.csv");
        if (!track) {
            throw cglbm::lbm::OutputError("cannot open mode.csv");
        }
        track.precision(12);
        track << "timestep,deformation\n";

        solver.initialize();
        writer.write_grids(0, solver.state());
        for (int timestep = 0; timestep <= config.steps; ++timestep) {
            if (timestep > 0) {
                solver.step();
            }
            if (timestep % kTrackInterval == 0) {
                track << timestep << "," << deformation(solver) << "\n";
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
        std::cerr << "oscillation: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
