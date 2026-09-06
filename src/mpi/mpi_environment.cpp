#include "mpi/mpi_environment.h"

#include "mpi/mpi_error.h"

#ifdef CGLBM_WITH_MPI
#include <mpi.h>
#endif

namespace cglbm {
namespace mpi {

Environment::Environment(int& argc, char**& argv) {
#ifdef CGLBM_WITH_MPI
    int already_initialized = 0;
    CGLBM_MPI_CHECK(MPI_Initialized(&already_initialized));
    if (!already_initialized) {
        CGLBM_MPI_CHECK(MPI_Init(&argc, &argv));
    }
#else
    (void) argc;
    (void) argv;
#endif
}

Environment::~Environment() {
#ifdef CGLBM_WITH_MPI
    int already_finalized = 0;
    // The destructor must not throw, so the return codes are not checked here.
    MPI_Finalized(&already_finalized);
    if (!already_finalized) {
        MPI_Finalize();
    }
#endif
}

bool Environment::available() {
#ifdef CGLBM_WITH_MPI
    return true;
#else
    return false;
#endif
}

int Environment::rank() {
#ifdef CGLBM_WITH_MPI
    int value = 0;
    CGLBM_MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &value));
    return value;
#else
    return 0;
#endif
}

int Environment::size() {
#ifdef CGLBM_WITH_MPI
    int value = 1;
    CGLBM_MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &value));
    return value;
#else
    return 1;
#endif
}

bool Environment::is_root() {
    return rank() == 0;
}

void Environment::barrier() {
#ifdef CGLBM_WITH_MPI
    CGLBM_MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
#endif
}

double Environment::sum(double value) {
#ifdef CGLBM_WITH_MPI
    double total = 0.0;
    CGLBM_MPI_CHECK(MPI_Allreduce(&value, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
    return total;
#else
    return value;
#endif
}

double Environment::max(double value) {
#ifdef CGLBM_WITH_MPI
    double largest = value;
    CGLBM_MPI_CHECK(MPI_Allreduce(&value, &largest, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD));
    return largest;
#else
    return value;
#endif
}

std::string Environment::describe() {
#ifdef CGLBM_WITH_MPI
    return "MPI enabled, rank " + std::to_string(rank()) + " of " + std::to_string(size());
#else
    return "MPI disabled, single rank";
#endif
}

}  // namespace mpi
}  // namespace cglbm
