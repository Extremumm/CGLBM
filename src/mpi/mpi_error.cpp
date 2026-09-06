#include "mpi/mpi_error.h"

#include <cstdlib>
#include <iostream>

#ifdef CGLBM_WITH_MPI
#include <mpi.h>
#endif

namespace cglbm {
namespace mpi {

std::string error_string(int errcode) {
#ifdef CGLBM_WITH_MPI
    char buffer[MPI_MAX_ERROR_STRING] = {};
    int length = 0;
    if (MPI_Error_string(errcode, buffer, &length) == MPI_SUCCESS && length > 0) {
        return std::string(buffer, static_cast<std::size_t>(length));
    }
    return "unknown MPI error " + std::to_string(errcode);
#else
    return "MPI error " + std::to_string(errcode) + " (built without MPI)";
#endif
}

void fail(int errcode, const std::string& what, const char* file, int line) {
    std::cerr << "MPI error at " << file << ":" << line << "\n"
              << "  call    : " << what << "\n"
              << "  code    : " << errcode << "\n"
              << "  message : " << error_string(errcode) << std::endl;
#ifdef CGLBM_WITH_MPI
    MPI_Abort(MPI_COMM_WORLD, errcode);
#endif
    std::abort();
}

void check(int errcode, const std::string& what, const char* file, int line) {
#ifdef CGLBM_WITH_MPI
    if (errcode != MPI_SUCCESS) {
        fail(errcode, what, file, line);
    }
#else
    (void) errcode;
    (void) what;
    (void) file;
    (void) line;
#endif
}

}  // namespace mpi
}  // namespace cglbm
