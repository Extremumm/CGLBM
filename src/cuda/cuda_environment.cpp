#include "cuda/cuda_environment.h"

#include <chrono>
#include <cstdlib>
#include <string>

#ifdef CGLBM_WITH_CUDA
#include <cuda_runtime.h>
#endif

namespace cglbm {
namespace cuda {

bool available() {
#ifdef CGLBM_WITH_CUDA
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    return (err == cudaSuccess && count > 0);
#else
    return false;
#endif
}

int device_count() {
#ifdef CGLBM_WITH_CUDA
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    if (err == cudaSuccess && count > 0) {
        return count;
    }
    return 0;
#else
    return 0;
#endif
}

int current_device() {
#ifdef CGLBM_WITH_CUDA
    int device = 0;
    if (cudaGetDevice(&device) == cudaSuccess) {
        return device;
    }
    return 0;
#else
    return 0;
#endif
}

bool set_device(int device_id) {
#ifdef CGLBM_WITH_CUDA
    if (device_id < 0 || device_id >= device_count()) {
        return false;
    }
    return (cudaSetDevice(device_id) == cudaSuccess);
#else
    (void) device_id;
    return false;
#endif
}

std::string device_name(int device_id) {
#ifdef CGLBM_WITH_CUDA
    cudaDeviceProp prop = {};
    if (cudaGetDeviceProperties(&prop, device_id) == cudaSuccess) {
        return std::string(prop.name);
    }
    return "";
#else
    (void) device_id;
    return "";
#endif
}

void compute_capability(int& major, int& minor, int device_id) {
#ifdef CGLBM_WITH_CUDA
    cudaDeviceProp prop = {};
    if (cudaGetDeviceProperties(&prop, device_id) == cudaSuccess) {
        major = prop.major;
        minor = prop.minor;
        return;
    }
    major = 0;
    minor = 0;
#else
    (void) device_id;
    major = 0;
    minor = 0;
#endif
}

void device_memory(std::size_t& free_bytes, std::size_t& total_bytes, int device_id) {
#ifdef CGLBM_WITH_CUDA
    int prev_device = current_device();
    if (cudaSetDevice(device_id) == cudaSuccess) {
        cudaMemGetInfo(&free_bytes, &total_bytes);
        cudaSetDevice(prev_device);
        return;
    }
    free_bytes = 0;
    total_bytes = 0;
#else
    (void) device_id;
    free_bytes = 0;
    total_bytes = 0;
#endif
}

void synchronize() {
#ifdef CGLBM_WITH_CUDA
    cudaDeviceSynchronize();
#endif
}

double wall_time() {
    const auto now = std::chrono::steady_clock::now().time_since_epoch();
    return std::chrono::duration<double>(now).count();
}

std::string describe() {
#ifdef CGLBM_WITH_CUDA
    const int count = device_count();
    if (count == 0) {
        return "CUDA enabled at compile time, but no CUDA device detected";
    }
    int major = 0;
    int minor = 0;
    compute_capability(major, minor, 0);
    std::size_t free_bytes = 0;
    std::size_t total_bytes = 0;
    device_memory(free_bytes, total_bytes, 0);
    const double total_mib = static_cast<double>(total_bytes) / (1024.0 * 1024.0);

    return "CUDA enabled, " + std::to_string(count) + " device" + (count > 1 ? "s" : "") + ": " +
           device_name(0) + " (sm_" + std::to_string(major) + std::to_string(minor) + ", " +
           std::to_string(static_cast<int>(total_mib)) + " MiB VRAM)";
#else
    return "CUDA disabled, host execution";
#endif
}

}  // namespace cuda
}  // namespace cglbm
