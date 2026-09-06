#ifndef CGLBM_MPI_MPI_ERROR_H
#define CGLBM_MPI_MPI_ERROR_H

#include <string>

namespace cglbm {
namespace mpi {

/// Turn an MPI return code into a message, then abort the job.
///
/// MPI's default error handler already aborts, but it says nothing about which
/// call failed; every call in this module is wrapped so that the file, the line
/// and the operation are reported first.
void fail(int errcode, const std::string& what, const char* file, int line);

/// Check one MPI return code, aborting through :func:`fail` when it is not
/// `MPI_SUCCESS`.
void check(int errcode, const std::string& what, const char* file, int line);

/// The message MPI associates with an error code.
std::string error_string(int errcode);

}  // namespace mpi
}  // namespace cglbm

/// Wrap an MPI call so that a failure reports where it happened.
#define CGLBM_MPI_CHECK(call) ::cglbm::mpi::check((call), #call, __FILE__, __LINE__)

#endif  // CGLBM_MPI_MPI_ERROR_H
