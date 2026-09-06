#ifndef CGLBM_CUDA_CUDA_ENVIRONMENT_H
#define CGLBM_CUDA_CUDA_ENVIRONMENT_H

#include <cstddef>
#include <string>

/// GPU acceleration via CUDA.
///
/// Everything here compiles and behaves sensibly whether or not the translation
/// unit was built with CUDA: without it the environment simply reports no
/// devices and unavailable, so a caller needs no `#ifdef` of its own.
namespace cglbm {
namespace cuda {

/// True when this build has CUDA enabled and at least one CUDA-capable GPU is available.
bool available();

/// Number of CUDA devices visible to the runtime: 0 without CUDA or when no device is present.
int device_count();

/// Currently selected active device index (0 without CUDA).
int current_device();

/// Select active device index. Returns false if device index is invalid or CUDA is unavailable.
bool set_device(int device_id);

/// Device name string (e.g. "NVIDIA GeForce RTX 3060"), or empty without CUDA.
std::string device_name(int device_id = 0);

/// Compute capability (major, minor), e.g. (8, 6), or (0, 0) without CUDA.
void compute_capability(int& major, int& minor, int device_id = 0);

/// Free and total device global memory in bytes, or (0, 0) without CUDA.
void device_memory(std::size_t& free_bytes, std::size_t& total_bytes, int device_id = 0);

/// Synchronize the current device (safe no-op without CUDA).
void synchronize();

/// Wall clock seconds, using a high-resolution timer.
double wall_time();

/// One line for the run log, e.g. `CUDA enabled, 1 device: NVIDIA GeForce RTX 3060 (sm_86, 12288
/// MiB VRAM)`.
std::string describe();

}  // namespace cuda
}  // namespace cglbm

#endif  // CGLBM_CUDA_CUDA_ENVIRONMENT_H
