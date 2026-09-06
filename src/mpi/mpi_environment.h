#ifndef CGLBM_MPI_MPI_ENVIRONMENT_H
#define CGLBM_MPI_MPI_ENVIRONMENT_H

#include <string>

/// Distributed-memory parallelism.
///
/// As in :file:`src/omp/omp_environment.h`, everything here compiles without
/// the dependency: a build with `WITH_MPI=OFF` behaves like a single rank, so a
/// solver written against this interface runs serially with no change.
namespace cglbm {
namespace mpi {

/// Initialise MPI for the lifetime of the object.
///
/// Construct exactly one, at the top of `main`, and let it go out of scope
/// last: the destructor finalises MPI. Copying is forbidden, since finalising
/// twice is an error.
class Environment {
  public:
    Environment(int& argc, char**& argv);
    ~Environment();

    Environment(const Environment&) = delete;
    Environment& operator=(const Environment&) = delete;

    /// True when this build has MPI enabled.
    static bool available();

    /// Rank of the calling process, 0 without MPI.
    static int rank();

    /// Number of ranks in the job, 1 without MPI.
    static int size();

    /// True on the rank that owns standard output and the log files.
    static bool is_root();

    /// Wait for every rank of `MPI_COMM_WORLD`.
    static void barrier();

    /// Sum `value` over all ranks and give the result to all of them.
    static double sum(double value);

    /// Largest `value` over all ranks, given to all of them.
    static double max(double value);

    /// One line for the run log, e.g. `MPI enabled, rank 0 of 4`.
    static std::string describe();
};

}  // namespace mpi
}  // namespace cglbm

#endif  // CGLBM_MPI_MPI_ENVIRONMENT_H
