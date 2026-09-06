#include "cuda/cuda_error.h"

#include <cstdlib>
#include <iostream>

#ifdef CGLBM_WITH_CUDA
#include <cuda_runtime.h>
#endif

namespace cglbm {
namespace cuda {

std::string error_string(int errcode) {
#ifdef CGLBM_WITH_CUDA
    const char* str = cudaGetErrorString(static_cast<cudaError_t>(errcode));
    if (str != nullptr) {
        return std::string(str);
    }
    return "unknown CUDA error " + std::to_string(errcode);
#else
    return "CUDA error " + std::to_string(errcode) + " (built without CUDA)";
#endif
}

void fail(int errcode, const std::string& what, const char* file, int line) {
    std::cerr << "CUDA error at " << file << ":" << line << "\n"
              << "  call    : " << what << "\n"
              << "  code    : " << errcode << "\n"
              << "  message : " << error_string(errcode) << std::endl;
    std::abort();
}

void check(int errcode, const std::string& what, const char* file, int line) {
#ifdef CGLBM_WITH_CUDA
    if (errcode != static_cast<int>(cudaSuccess)) {
        fail(errcode, what, file, line);
    }
#else
    (void) errcode;
    (void) what;
    (void) file;
    (void) line;
#endif
}

}  // namespace cuda
}  // namespace cglbm
