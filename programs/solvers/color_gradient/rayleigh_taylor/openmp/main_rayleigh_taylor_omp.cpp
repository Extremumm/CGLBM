/// Rayleigh-Taylor instability at production resolution, across OpenMP threads.
///
/// The same case as `rayleigh_taylor`, on a 1024 x 4096 lattice. Every lattice
/// loop is spread over the threads; each writes only its own node, so the
/// result does not depend on how many threads run it.
///
/// The lattice is allocated on the heap rather than in `.bss`, so this program
/// no longer needs the `-mcmodel=medium` exception that the static arrays used
/// to force. It still asks for several gigabytes -- check the available memory
/// before launching it, or lower the resolution with `--nx` and `--ny`.

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "lbm/case_config.h"
#include "lbm/solver.h"
#include "omp/omp_environment.h"

namespace {

using cglbm::lbm::CaseConfig;

const double c_dx = 1.e-5;                        // m
const double c_dt = c_dx / 347. / std::sqrt(3.);  // s

/// Threads to use when `OMP_NUM_THREADS` says nothing, as this case always has.
const int kDefaultThreads = 8;

CaseConfig rayleigh_taylor_omp_case() {
    CaseConfig config;
    config.name = "rayleigh_taylor_omp";
    config.nx = 1024;
    config.ny = 4096;
    config.steps = 2000000;
    config.interval = 10000;

    config.physics.rho1 = 4.;
    config.physics.rho2 = 1.;
    config.physics.c1 = 347. / (c_dx / c_dt);
    config.physics.c2 = 347. / (c_dx / c_dt);
    config.physics.radius = 10.;
    config.physics.sigma = 0. / (c_dx * c_dx * c_dx / c_dt / c_dt);
    config.physics.nu = 1.e-4 / (c_dx * c_dx / c_dt);
    config.physics.nu_b = 1.e-4 / (c_dx * c_dx / c_dt);
    config.physics.gravity = 9.81e2 / (c_dx / c_dt / c_dt);  // m s^-2, lattice units
    config.physics.ch_width_init = 1.1 * config.units.dx;
    config.physics.ch_width_ope = 1.6 * config.units.dx;
    config.physics.p2_inf = 0.;
    config.physics.p1_inf = cglbm::lbm::matched_p1_inf(config.physics);

    config.boundary = cglbm::lbm::Boundary::WallY;
    // A smaller perturbation than the reference case: at this resolution the
    // instability has room to select its own wavelength.
    config.initial_phase = cglbm::lbm::cosine_layer(0.1, /*inverted=*/false);
    config.stencil = cglbm::lbm::GradientStencil::E4;
    config.parallel = true;

    return config;
}

}  // namespace

int main(int argc, char** argv) {
    // Thread count: OMP_NUM_THREADS when it is set, otherwise this case's
    // historical default. Going through src/omp keeps the program linkable in a
    // build configured with WITH_OpenMP=OFF, where it runs serially.
    cglbm::omp::set_thread_count(std::getenv("OMP_NUM_THREADS") ? 0 : kDefaultThreads);

    CaseConfig config = rayleigh_taylor_omp_case();
    switch (cglbm::lbm::parse_command_line(config, argc, argv, "rayleigh_taylor_omp")) {
    case cglbm::lbm::CommandLineResult::Finished:
        return 0;
    case cglbm::lbm::CommandLineResult::Error:
        return 2;
    case cglbm::lbm::CommandLineResult::Run:
        break;
    }

    std::cout << cglbm::omp::describe() << std::endl;
    if (!cglbm::omp::available()) {
        std::cerr << "Warning: built without OpenMP, running serially." << std::endl;
    }
    std::cout << cglbm::lbm::describe(config) << std::endl;

    try {
        cglbm::lbm::Solver solver(config);
        solver.run();
    } catch (const std::exception& error) {
        std::cerr << "rayleigh_taylor_omp: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
