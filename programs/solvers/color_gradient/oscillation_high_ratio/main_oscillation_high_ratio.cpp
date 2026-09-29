/// Mode-2 oscillation of a droplet at a density ratio of 1000, in the
/// two-population solver.
///
/// The same benchmark as `oscillation`, run with
/// `cglbm::lbm::TwoPopulationSolver` instead of `cglbm::lbm::Solver`, as
/// `laplace_high_ratio` is to `laplace`: a droplet of radius 20, deformed by
/// 3 % into `cos 2 theta`, released and ringing down, scored against the exact
/// normal mode of two viscous fluids (`pycglbm.normal_modes.droplet_mode`).
/// It is the moving-interface case of Ba et al.'s model, whose static droplet
/// `laplace_high_ratio` holds.
///
/// The model's own parameters follow `laplace_high_ratio` and Ba et al.
/// (2016): alpha_2 = 0.2, beta = 0.7, sigma = 0.1, and their third-moment
/// source term. The viscosities are `oscillation`'s, dynamic viscosities of 2
/// in the droplet and 0.05 around it, which puts tau at 4.7 and 0.60; a BGK
/// collision diverges there within a few hundred steps, so the collision is
/// Ba et al.'s MRT.
///
/// `mode.csv` records every `kTrackInterval` steps the droplet's deformation,
/// the second moment of its excess density
///
///     D = sum (rho - rho2) (x^2 - y^2) / sum (rho - rho2) (x^2 + y^2),
///
/// about the domain centre, which is `2 eps` for a small deformation `eps`.

#include <cmath>
#include <fstream>
#include <iostream>

#include "lbm/case_config.h"
#include "lbm/output_writer.h"
#include "lbm/two_population_solver.h"

namespace {

using cglbm::lbm::CaseConfig;

/// Initial mode-2 deformation, as a fraction of the radius.
constexpr double kDeformation = 0.03;

/// How often `mode.csv` is written.
constexpr int kTrackInterval = 50;

CaseConfig oscillation_high_ratio_case() {
    CaseConfig config;
    config.name = "oscillation_high_ratio";
    config.nx = 128;
    config.ny = 128;
    config.steps = 24000;
    config.interval = 6000;

    config.physics.rho1 = 1000.;  // the droplet
    config.physics.rho2 = 1.;     // the surrounding fluid
    config.physics.radius = 20.;
    config.physics.sigma = 0.1;
    // Dynamic viscosities of 2 in the droplet and 0.05 around it, as in
    // `oscillation`.
    config.physics.nu = 2.e-3;
    config.physics.nu_b = 2.e-3;
    config.physics.nu2 = 0.05;
    config.physics.nu_b2 = 0.05;

    config.physics.alpha2 = 0.2;
    config.physics.beta = 0.7;
    config.physics.ch_width_init = 1.1;

    config.boundary = cglbm::lbm::Boundary::PeriodicY;
    config.stencil = cglbm::lbm::GradientStencil::E8;
    config.collision = cglbm::lbm::Collision::MRT;
    config.third_moment_correction = true;

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
double deformation(const cglbm::lbm::TwoPopulationSolver& solver) {
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
    CaseConfig config = oscillation_high_ratio_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "oscillation_high_ratio")) {
    case cglbm::lbm::CommandLineResult::Finished:
        return 0;
    case cglbm::lbm::CommandLineResult::Error:
        return 2;
    case cglbm::lbm::CommandLineResult::Run:
        break;
    }

    std::cout << cglbm::lbm::describe(config) << "\ndeformation = " << kDeformation << std::endl;
    try {
        cglbm::lbm::TwoPopulationSolver solver(config);
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
            // step() leaves the macroscopic fields as they were before the
            // collision; refresh() brings them up to the populations.
            if (timestep % kTrackInterval == 0) {
                solver.refresh();
                track << timestep << "," << deformation(solver) << "\n";
            }
            if (timestep > 0 && timestep % config.interval == 0) {
                std::cout << "Step " << timestep << std::endl;
                solver.refresh();
                writer.write_grids(timestep, solver.state());
            }
        }
        if (!track) {
            throw cglbm::lbm::OutputError("cannot write mode.csv");
        }
    } catch (const std::exception& error) {
        std::cerr << "oscillation_high_ratio: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
