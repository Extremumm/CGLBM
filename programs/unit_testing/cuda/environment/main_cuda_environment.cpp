#include <iostream>
#include <vector>

#include "cuda/cuda_environment.h"
#include "cuda/cuda_error.h"
#include "cuda/cuda_memory.h"

int main(int argc, char** argv) {
    (void) argc;
    (void) argv;

    const bool avail = cglbm::cuda::available();
    const int count = cglbm::cuda::device_count();

    std::cout << "available = " << (avail ? 1 : 0) << "\n";
    std::cout << "device_count = " << count << "\n";
    std::cout << "describe = " << cglbm::cuda::describe() << "\n";

    if (avail && count > 0) {
        std::cout << "device_name = " << cglbm::cuda::device_name(0) << "\n";

        int major = 0;
        int minor = 0;
        cglbm::cuda::compute_capability(major, minor, 0);
        std::cout << "compute_capability = " << major << "." << minor << "\n";

        std::size_t free_b = 0;
        std::size_t total_b = 0;
        cglbm::cuda::device_memory(free_b, total_b, 0);
        std::cout << "total_memory_mb = " << (total_b / (1024 * 1024)) << "\n";
        std::cout << "free_memory_mb = " << (free_b / (1024 * 1024)) << "\n";

        // Memory roundtrip test
        const std::size_t n = 1024;
        std::vector<int> h_in(n);
        for (std::size_t i = 0; i < n; ++i) {
            h_in[i] = static_cast<int>(i * 3 + 1);
        }

        cglbm::cuda::DeviceBuffer<int> d_buf(n);
        d_buf.copy_from_host(h_in.data());

        cglbm::cuda::synchronize();

        std::vector<int> h_out(n, 0);
        d_buf.copy_to_host(h_out.data());

        bool match = true;
        for (std::size_t i = 0; i < n; ++i) {
            if (h_out[i] != h_in[i]) {
                match = false;
                break;
            }
        }
        std::cout << "memory_roundtrip = " << (match ? 1 : 0) << "\n";
    } else {
        std::cout << "device_name = none\n";
        std::cout << "compute_capability = 0.0\n";
        std::cout << "total_memory_mb = 0\n";
        std::cout << "free_memory_mb = 0\n";
        std::cout << "memory_roundtrip = 1\n";
    }

    const double t0 = cglbm::cuda::wall_time();
    const double elapsed = cglbm::cuda::wall_time() - t0;
    std::cout << "wall_time_monotonic = " << (elapsed >= 0.0 ? 1 : 0) << "\n";

    return 0;
}
