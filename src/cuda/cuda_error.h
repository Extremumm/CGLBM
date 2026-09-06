#ifndef CGLBM_CUDA_CUDA_ERROR_H
#define CGLBM_CUDA_CUDA_ERROR_H

#include <string>

namespace cglbm {
namespace cuda {

/// Turn a CUDA return code into a message, then abort the program.
void fail(int errcode, const std::string& what, const char* file, int line);

/// Check one CUDA return code, aborting through :func:`fail` when it is not
/// `cudaSuccess` (0).
void check(int errcode, const std::string& what, const char* file, int line);

/// The message CUDA associates with an error code.
std::string error_string(int errcode);

}  // namespace cuda
}  // namespace cglbm

/// Wrap a CUDA call so that a failure reports where it happened.
#define CGLBM_CUDA_CHECK(call)                                                                     \
    ::cglbm::cuda::check(static_cast<int>(call), #call, __FILE__, __LINE__)

#endif  // CGLBM_CUDA_CUDA_ERROR_H
