#ifndef CGLBM_OMP_OMP_ENVIRONMENT_H
#define CGLBM_OMP_OMP_ENVIRONMENT_H

#include <string>

/// Thread-level parallelism.
///
/// Everything here compiles and behaves sensibly whether or not the translation
/// unit was built with OpenMP: without it the environment simply reports a
/// single thread, so a solver needs no `#ifdef` of its own.
namespace cglbm {
namespace omp {

/// True when this build has OpenMP enabled.
bool available();

/// Threads a parallel region would use, i.e. `omp_get_max_threads()`.
int max_threads();

/// Threads in the enclosing parallel region: 1 outside of one, and always 1
/// without OpenMP.
int thread_count();

/// Index of the calling thread, 0 outside a parallel region.
int thread_id();

/// Request `count` threads for the following parallel regions.
///
/// A `count` below 1 restores the runtime default (`OMP_NUM_THREADS`, else the
/// number of available cores). Returns the value actually in force afterwards,
/// which is 1 without OpenMP.
int set_thread_count(int count);

/// Wall clock seconds, from `omp_get_wtime()` when available.
double wall_time();

/// One line for the run log, e.g. `OpenMP enabled, 8 threads (of 24 cores)`.
std::string describe();

}  // namespace omp
}  // namespace cglbm

#endif  // CGLBM_OMP_OMP_ENVIRONMENT_H
