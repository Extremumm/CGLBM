#include "omp/omp_environment.h"

#include <chrono>
#include <cstdlib>
#include <string>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cglbm {
namespace omp {

bool available() {
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}

int max_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

int thread_count() {
#ifdef _OPENMP
    return omp_get_num_threads();
#else
    return 1;
#endif
}

int thread_id() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

int set_thread_count(int count) {
#ifdef _OPENMP
    if (count < 1) {
        // Back to the runtime default: OMP_NUM_THREADS when it is set to a
        // usable value, otherwise one thread per available core. Reading the
        // variable is necessary because omp_set_num_threads() elsewhere may
        // already have overridden it.
        int fallback = omp_get_num_procs();
        if (const char* requested = std::getenv("OMP_NUM_THREADS")) {
            const int parsed = std::atoi(requested);
            if (parsed > 0) {
                fallback = parsed;
            }
        }
        omp_set_num_threads(fallback);
    } else {
        omp_set_num_threads(count);
    }
    return omp_get_max_threads();
#else
    (void) count;
    return 1;
#endif
}

double wall_time() {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    const auto now = std::chrono::steady_clock::now().time_since_epoch();
    return std::chrono::duration<double>(now).count();
#endif
}

std::string describe() {
#ifdef _OPENMP
    const unsigned int cores = std::thread::hardware_concurrency();
    return "OpenMP enabled, " + std::to_string(max_threads()) + " threads (of " +
           std::to_string(cores) + " cores)";
#else
    return "OpenMP disabled, serial execution";
#endif
}

}  // namespace omp
}  // namespace cglbm
