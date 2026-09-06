// Reports what the OpenMP module sees, and checks that a parallel region really
// runs on the requested number of threads. The pytest beside this file parses
// the key = value lines.

#include <cstdlib>
#include <iostream>

#include "omp/omp_environment.h"

int main(int argc, char** argv) {
    const int requested = (argc > 1) ? std::atoi(argv[1]) : 0;

    const int in_force = cglbm::omp::set_thread_count(requested);

    std::cout << "available = " << (cglbm::omp::available() ? 1 : 0) << "\n";
    std::cout << "requested = " << requested << "\n";
    std::cout << "max_threads = " << in_force << "\n";
    std::cout << "describe = " << cglbm::omp::describe() << "\n";

    // Count the threads that actually enter a parallel region, and check that
    // every thread id appears exactly once.
    int observed = 0;
    int id_sum = 0;
#pragma omp parallel reduction(+ : observed, id_sum)
    {
        observed += 1;
        id_sum += cglbm::omp::thread_id();
    }

    std::cout << "observed_threads = " << observed << "\n";
    // 0 + 1 + ... + (n-1)
    std::cout << "expected_id_sum = " << observed * (observed - 1) / 2 << "\n";
    std::cout << "id_sum = " << id_sum << "\n";

    const double start = cglbm::omp::wall_time();
    const double elapsed = cglbm::omp::wall_time() - start;
    std::cout << "wall_time_monotonic = " << (elapsed >= 0.0 ? 1 : 0) << "\n";

    return 0;
}
