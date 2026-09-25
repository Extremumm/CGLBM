/// Regenerates the measurements plotted in docs/report/report.tex.
///
/// Build and run it from the repository root:
///
///     g++ -O2 -std=c++17 -ffp-contract=off -fopenmp -I src
///         -o report_data docs/report/tools/report_data.cpp src/lbm/*.cpp src/omp/*.cpp
///     ./report_data <task> docs/report/figures/data
///
/// Tasks: `sweep`, `history`, `profiles`, `currents`. Each writes one CSV. The
/// results are committed under figures/data so the report builds without a
/// three-hour recompute; this is how they were produced.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "lbm/case_config.h"
#include "lbm/isotropic_gradient.h"
#include "lbm/solver.h"
#include "lbm/two_population_solver.h"

using namespace cglbm::lbm;

namespace {

const double c_dx = 1.e-5;
const double c_dt = c_dx / 347. / std::sqrt(3.);

/// Where a monotone profile crosses `level`, linearly interpolated.
double crossing(const std::vector<double>& profile, double level) {
    for (std::size_t k = 0; k + 1 < profile.size(); k++) {
        const double a = profile[k] - level, b = profile[k + 1] - level;
        if (a >= 0.0 && b < 0.0) {
            return double(k) + a / (a - b);
        }
    }
    return NAN;
}

/// The equation-of-state solver on the Laplace droplet.
CaseConfig eos_case(double ratio, bool ba_configuration) {
    CaseConfig config;
    config.name = "eos";
    config.nx = 128;
    config.ny = 128;
    config.interval = 1 << 30;
    config.physics.rho1 = ratio;
    config.physics.rho2 = 1.;
    config.physics.c1 = 347. / (c_dx / c_dt);
    config.physics.c2 = config.physics.c1;
    config.physics.radius = 10.;
    config.physics.sigma = 1. / (c_dx * c_dx * c_dx / c_dt / c_dt);
    config.physics.nu = 1.e-2 / (c_dx * c_dx / c_dt);
    config.physics.nu_b = config.physics.nu;
    config.physics.ch_width_init = 1.1;
    config.physics.ch_width_ope = 1.6;
    config.physics.p2_inf = 0.;
    config.physics.p1_inf = matched_p1_inf(config.physics);
    config.boundary = Boundary::PeriodicY;
    config.initial_phase = droplet_interface();
    config.stencil = GradientStencil::E8;
    if (ba_configuration) {
        config.initial_profile_field = InterfaceField::BulkNormalised;
        config.interface_field = InterfaceField::BulkNormalised;
        config.surface_tension = SurfaceTension::ContinuumSurfaceForce;
    } else {
        config.initial_profile_field = InterfaceField::Colour;
        config.interface_field = InterfaceField::Colour;
        config.surface_tension = SurfaceTension::Perturbation;
    }
    return config;
}

/// The two-population solver on Ba et al.'s droplet, matched dynamic viscosity.
CaseConfig two_population_case(double ratio) {
    const double mu = 0.1667;
    CaseConfig config;
    config.name = "two_population";
    config.nx = 100;
    config.ny = 100;
    config.interval = 1 << 30;
    config.physics.rho1 = ratio;
    config.physics.rho2 = 1.;
    config.physics.radius = 25.;
    config.physics.sigma = 0.1;
    config.physics.nu = mu / ratio;
    config.physics.nu_b = config.physics.nu;
    config.physics.nu2 = mu;
    config.physics.nu_b2 = mu;
    config.physics.alpha2 = 0.2;
    config.physics.beta = 0.7;
    config.physics.ch_width_init = 1.1;
    config.boundary = Boundary::PeriodicY;
    config.initial_phase = droplet_interface();
    config.stencil = GradientStencil::E8;
    return config;
}

/// Pressure jump, interface radius and peak speed of a relaxed droplet.
template <typename S>
bool measure(
    const S& solver, const CaseConfig& config, double* jump, double* radius, double* peak_speed) {
    const Field& p = solver.pressure();
    const Field& phase = solver.phase();
    const Field& u = solver.velocity();
    const int nx = config.nx;
    const double R = config.physics.radius;

    double inside = 0.0, outside = 0.0;
    long n_inside = 0, n_outside = 0;
    *peak_speed = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < nx; j++) {
            if (!std::isfinite(p(i, j)) || !std::isfinite(u(i, j, 0))) {
                return false;
            }
            const double r = std::hypot(i - nx / 2.0, j - nx / 2.0);
            if (r < 0.5 * R) {
                inside += p(i, j);
                n_inside++;
            }
            if (r > 1.5 * R) {
                outside += p(i, j);
                n_outside++;
            }
            *peak_speed = std::max(*peak_speed, std::hypot(u(i, j, 0), u(i, j, 1)));
        }
    }
    std::vector<double> profile;
    for (int i = nx / 2; i < nx; i++) {
        profile.push_back(phase(i, nx / 2));
    }
    // The equation-of-state solver's phase field is the colour field, whose
    // zero is not the interface; its interface is the density midpoint. The
    // two-population solver reports phi_N, whose zero is the interface.
    *radius = crossing(profile, 0.0);
    *jump = inside / n_inside - outside / n_outside;
    return std::isfinite(*jump);
}

/// Interface radius of the equation-of-state solver: the density midpoint.
double density_radius(const Solver& solver, int nx) {
    const Field& rho = solver.density();
    std::vector<double> profile;
    for (int i = nx / 2; i < nx; i++) {
        profile.push_back(rho(i, nx / 2));
    }
    return crossing(profile, 0.5 * (profile.front() + profile.back()));
}

void task_sweep(const std::string& out, double only_ratio) {
    // One density ratio per process when `only_ratio` is given, so the sweep
    // can be run in parallel; the shards are concatenated afterwards.
    const std::string suffix = only_ratio > 0 ? "_" + std::to_string(long(only_ratio)) : "";
    FILE* fp = fopen((out + "/laplace_sweep" + suffix + ".csv").c_str(), "w");
    if (only_ratio <= 0) {
        fprintf(fp, "model,ratio,sigma_ratio,peak_speed,survived\n");
    }
    const double all[] = {2., 5., 10., 20., 50., 100., 200., 500., 1000., 1.e4, 1.e5};
    std::vector<double> ratios;
    for (double candidate : all) {
        if (only_ratio <= 0 || candidate == only_ratio) {
            ratios.push_back(candidate);
        }
    }

    for (double ratio : ratios) {
        for (int variant = 0; variant < 2; variant++) {
            CaseConfig config = eos_case(ratio, variant == 1);
            config.steps = 120000;
            Solver solver(config);
            solver.initialize();
            bool alive = true;
            for (int t = 1; t <= config.steps && alive; t++) {
                solver.step();
                if (t % 2000 == 0) {
                    double j, r, s;
                    alive = measure(solver, config, &j, &r, &s) && s < 1.0;
                }
            }
            double jump = NAN, radius = NAN, speed = NAN;
            if (alive) {
                measure(solver, config, &jump, &radius, &speed);
                radius = density_radius(solver, config.nx);
            }
            fprintf(fp,
                    "%s,%g,%.8g,%.8g,%d\n",
                    variant == 1 ? "eos_csf" : "eos_legacy",
                    ratio,
                    alive ? jump * radius / config.physics.sigma : NAN,
                    alive ? speed : NAN,
                    alive ? 1 : 0);
            fflush(fp);
        }

        CaseConfig config = two_population_case(ratio);
        config.steps = 60000;
        TwoPopulationSolver solver(config);
        solver.initialize();
        bool alive = true;
        for (int t = 1; t <= config.steps && alive; t++) {
            solver.step();
            if (t % 2000 == 0) {
                solver.refresh();
                double j, r, s;
                alive = measure(solver, config, &j, &r, &s) && s < 1.0;
            }
        }
        solver.refresh();
        double jump = NAN, radius = NAN, speed = NAN;
        if (alive) {
            measure(solver, config, &jump, &radius, &speed);
        }
        fprintf(fp,
                "two_population,%g,%.8g,%.8g,%d\n",
                ratio,
                alive ? jump * radius / config.physics.sigma : NAN,
                alive ? speed : NAN,
                alive ? 1 : 0);
        fflush(fp);
    }
    fclose(fp);
}

void task_history(const std::string& out) {
    FILE* fp = fopen((out + "/history_r1000.csv").c_str(), "w");
    fprintf(fp, "model,step,sigma_ratio,peak_speed\n");
    const int report_every = 4000;

    for (int variant = 0; variant < 2; variant++) {
        CaseConfig config = eos_case(1000., variant == 1);
        config.steps = 120000;
        Solver solver(config);
        solver.initialize();
        for (int t = 1; t <= config.steps; t++) {
            solver.step();
            if (t % report_every) {
                continue;
            }
            double j, r, s;
            if (!measure(solver, config, &j, &r, &s) || s >= 1.0) {
                break;
            }
            fprintf(fp,
                    "%s,%d,%.8g,%.8g\n",
                    variant == 1 ? "eos_csf" : "eos_legacy",
                    t,
                    j * density_radius(solver, config.nx) / config.physics.sigma,
                    s);
            fflush(fp);
        }
    }

    CaseConfig config = two_population_case(1000.);
    config.steps = 120000;
    TwoPopulationSolver solver(config);
    solver.initialize();
    for (int t = 1; t <= config.steps; t++) {
        solver.step();
        if (t % report_every) {
            continue;
        }
        solver.refresh();
        double j, r, s;
        if (!measure(solver, config, &j, &r, &s) || s >= 1.0) {
            break;
        }
        fprintf(fp, "two_population,%d,%.8g,%.8g\n", t, j * r / config.physics.sigma, s);
        fflush(fp);
    }
    fclose(fp);
}

void task_profiles(const std::string& out) {
    FILE* fp = fopen((out + "/profiles_r1000.csv").c_str(), "w");
    fprintf(fp, "model,r,phase,rho,p\n");

    CaseConfig eos = eos_case(1000., true);
    eos.steps = 30000;
    Solver eos_solver(eos);
    eos_solver.initialize();
    for (int t = 1; t <= eos.steps; t++) {
        eos_solver.step();
    }
    for (int i = eos.nx / 2; i < eos.nx; i++) {
        fprintf(fp,
                "eos_csf,%d,%.8g,%.8g,%.8g\n",
                i - eos.nx / 2,
                eos_solver.phase()(i, eos.ny / 2),
                eos_solver.density()(i, eos.ny / 2),
                eos_solver.pressure()(i, eos.ny / 2));
    }

    CaseConfig tp = two_population_case(1000.);
    tp.steps = 40000;
    TwoPopulationSolver tp_solver(tp);
    tp_solver.initialize();
    for (int t = 1; t <= tp.steps; t++) {
        tp_solver.step();
    }
    tp_solver.refresh();
    for (int i = tp.nx / 2; i < tp.nx; i++) {
        fprintf(fp,
                "two_population,%d,%.8g,%.8g,%.8g\n",
                i - tp.nx / 2,
                tp_solver.phase()(i, tp.ny / 2),
                tp_solver.density()(i, tp.ny / 2),
                tp_solver.pressure()(i, tp.ny / 2));
    }
    fclose(fp);
}

void task_currents(const std::string& out) {
    CaseConfig tp = two_population_case(1000.);
    tp.steps = 40000;
    TwoPopulationSolver solver(tp);
    solver.initialize();
    for (int t = 1; t <= tp.steps; t++) {
        solver.step();
    }
    solver.refresh();
    FILE* fp = fopen((out + "/currents_r1000.csv").c_str(), "w");
    fprintf(fp, "i,j,ux,uy,phase\n");
    for (int i = 0; i < tp.nx; i++) {
        for (int j = 0; j < tp.ny; j++) {
            fprintf(fp,
                    "%d,%d,%.6g,%.6g,%.6g\n",
                    i,
                    j,
                    solver.velocity()(i, j, 0),
                    solver.velocity()(i, j, 1),
                    solver.phase()(i, j));
        }
    }
    fclose(fp);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr,
                     "usage: report_data <sweep|history|profiles|currents> <out-dir> "
                     "[ratio]\n");
        return 2;
    }
    const std::string task = argv[1];
    const std::string out = argv[2];
    if (task == "sweep") {
        task_sweep(out, argc > 3 ? std::atof(argv[3]) : -1.0);
    } else if (task == "history") {
        task_history(out);
    } else if (task == "profiles") {
        task_profiles(out);
    } else if (task == "currents") {
        task_currents(out);
    } else {
        std::fprintf(stderr, "unknown task '%s'\n", task.c_str());
        return 2;
    }
    return 0;
}
