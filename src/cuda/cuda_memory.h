#ifndef CGLBM_CUDA_CUDA_MEMORY_H
#define CGLBM_CUDA_CUDA_MEMORY_H

#include <cstddef>
#include <stdexcept>
#include <utility>

#include "cuda/cuda_error.h"

#ifdef CGLBM_WITH_CUDA
#include <cuda_runtime.h>
#endif

namespace cglbm {
namespace cuda {

/// Allocate `count` elements of type T on the active GPU device.
template <typename T> T* allocate_device(std::size_t count) {
#ifdef CGLBM_WITH_CUDA
    T* ptr = nullptr;
    CGLBM_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&ptr), count * sizeof(T)));
    return ptr;
#else
    (void) count;
    throw std::runtime_error("allocate_device called in build without CUDA");
#endif
}

/// Free memory allocated on device with allocate_device.
template <typename T> void free_device(T* ptr) {
#ifdef CGLBM_WITH_CUDA
    if (ptr != nullptr) {
        CGLBM_CUDA_CHECK(cudaFree(ptr));
    }
#else
    (void) ptr;
#endif
}

/// Copy `count` elements of type T from host to device.
template <typename T> void copy_to_device(T* dst_device, const T* src_host, std::size_t count) {
#ifdef CGLBM_WITH_CUDA
    CGLBM_CUDA_CHECK(cudaMemcpy(dst_device, src_host, count * sizeof(T), cudaMemcpyHostToDevice));
#else
    (void) dst_device;
    (void) src_host;
    (void) count;
    throw std::runtime_error("copy_to_device called in build without CUDA");
#endif
}

/// Copy `count` elements of type T from device to host.
template <typename T> void copy_to_host(T* dst_host, const T* src_device, std::size_t count) {
#ifdef CGLBM_WITH_CUDA
    CGLBM_CUDA_CHECK(cudaMemcpy(dst_host, src_device, count * sizeof(T), cudaMemcpyDeviceToHost));
#else
    (void) dst_host;
    (void) src_device;
    (void) count;
    throw std::runtime_error("copy_to_host called in build without CUDA");
#endif
}

/// RAII managed device buffer.
template <typename T> class DeviceBuffer {
public:
    DeviceBuffer() : ptr_(nullptr), size_(0) {}

    explicit DeviceBuffer(std::size_t size) : ptr_(nullptr), size_(size) {
        if (size_ > 0) {
            ptr_ = allocate_device<T>(size_);
        }
    }

    ~DeviceBuffer() {
        if (ptr_ != nullptr) {
            free_device(ptr_);
            ptr_ = nullptr;
        }
    }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    DeviceBuffer(DeviceBuffer&& other) noexcept : ptr_(other.ptr_), size_(other.size_) {
        other.ptr_ = nullptr;
        other.size_ = 0;
    }

    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
        if (this != &other) {
            if (ptr_ != nullptr) {
                free_device(ptr_);
            }
            ptr_ = other.ptr_;
            size_ = other.size_;
            other.ptr_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }

    void allocate(std::size_t size) {
        if (ptr_ != nullptr) {
            free_device(ptr_);
            ptr_ = nullptr;
        }
        size_ = size;
        if (size_ > 0) {
            ptr_ = allocate_device<T>(size_);
        }
    }

    void copy_from_host(const T* host_ptr, std::size_t count = 0) {
        const std::size_t n = (count == 0) ? size_ : count;
        copy_to_device(ptr_, host_ptr, n);
    }

    void copy_to_host(T* host_ptr, std::size_t count = 0) const {
        const std::size_t n = (count == 0) ? size_ : count;
        cglbm::cuda::copy_to_host(host_ptr, ptr_, n);
    }

    void zero() {
#ifdef CGLBM_WITH_CUDA
        if (ptr_ != nullptr && size_ > 0) {
            CGLBM_CUDA_CHECK(cudaMemset(ptr_, 0, size_ * sizeof(T)));
        }
#endif
    }

    T* data() {
        return ptr_;
    }
    const T* data() const {
        return ptr_;
    }
    std::size_t size() const {
        return size_;
    }
    std::size_t bytes() const {
        return size_ * sizeof(T);
    }

private:
    T* ptr_;
    std::size_t size_;
};

}  // namespace cuda
}  // namespace cglbm

#endif  // CGLBM_CUDA_CUDA_MEMORY_H
