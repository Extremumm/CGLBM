// Checks the halo exchange, corners included.
//
// Every interior node is filled with a value that identifies its global
// position uniquely. After the exchange, each ghost node must hold the value of
// the global node it shadows -- including the four corners, which D2Q9 reaches
// through its diagonal velocities and which the exchange fills without ever
// sending a diagonal message.
//
//   main_mpi_exchange_ghosts [global_nx] [global_ny] [depth] [periodic_y]

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "mpi/mpi_environment.h"
#include "mpi/mpi_exchange_ghosts.h"
#include "mpi/mpi_topology.h"

namespace {

/// A value that identifies the global node (gx, gy) and the component q.
double node_value(int gx, int gy, int q, int global_ny, int depth) {
    return static_cast<double>((gx * global_ny + gy) * depth + q);
}

/// Wrap a global index onto the lattice.
int wrap(int index, int extent) { return (index % extent + extent) % extent; }

}  // namespace

int main(int argc, char** argv) {
    cglbm::mpi::Environment environment(argc, argv);

    const int global_nx = (argc > 1) ? std::atoi(argv[1]) : 64;
    const int global_ny = (argc > 2) ? std::atoi(argv[2]) : 64;
    const int depth = (argc > 3) ? std::atoi(argv[3]) : 1;
    const bool periodic_y = (argc > 4) ? std::atoi(argv[4]) != 0 : true;
    const bool periodic_x = true;

    cglbm::mpi::CartesianTopology topology(global_nx, global_ny, periodic_x, periodic_y);
    const cglbm::mpi::Decomposition& local = topology.local();

    const int nx = local.nx;
    const int ny = local.ny;
    const int stride_j = depth;
    const int stride_i = (ny + 2) * depth;

    // -1 marks a ghost node nobody has written yet
    std::vector<double> field(static_cast<std::size_t>(local.padded_size()) * depth, -1.0);

    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            const int gx = local.x_offset + i - 1;
            const int gy = local.y_offset + j - 1;
            for (int q = 0; q < depth; ++q) {
                field[i * stride_i + j * stride_j + q] = node_value(gx, gy, q, global_ny, depth);
            }
        }
    }

    cglbm::mpi::exchange_ghosts(topology, field.data(), nx, ny, depth);

    int checked = 0;
    int wrong = 0;
    int untouched = 0;

    for (int i = 0; i <= nx + 1; ++i) {
        for (int j = 0; j <= ny + 1; ++j) {
            const bool is_interior = (i >= 1 && i <= nx && j >= 1 && j <= ny);
            if (is_interior) {
                continue;
            }

            const int gx_raw = local.x_offset + i - 1;
            const int gy_raw = local.y_offset + j - 1;

            // A ghost beyond a non-periodic edge is left to the solver's own
            // boundary condition, so it must still hold the sentinel.
            const bool outside_y = !periodic_y && (gy_raw < 0 || gy_raw >= global_ny);
            if (outside_y) {
                for (int q = 0; q < depth; ++q) {
                    if (field[i * stride_i + j * stride_j + q] != -1.0) {
                        ++wrong;
                    }
                }
                ++untouched;
                continue;
            }

            const int gx = wrap(gx_raw, global_nx);
            const int gy = wrap(gy_raw, global_ny);
            for (int q = 0; q < depth; ++q) {
                const double expected = node_value(gx, gy, q, global_ny, depth);
                const double found = field[i * stride_i + j * stride_j + q];
                if (std::fabs(found - expected) > 1.0e-12) {
                    if (wrong < 5) {
                        std::cerr << "rank " << cglbm::mpi::Environment::rank() << " ghost (" << i
                                  << "," << j << ") q=" << q << " expected " << expected
                                  << " found " << found << std::endl;
                    }
                    ++wrong;
                }
                ++checked;
            }
        }
    }

    const int total_wrong = static_cast<int>(cglbm::mpi::Environment::sum(wrong));
    const int total_checked = static_cast<int>(cglbm::mpi::Environment::sum(checked));
    const int total_untouched = static_cast<int>(cglbm::mpi::Environment::sum(untouched));

    if (cglbm::mpi::Environment::is_root()) {
        std::cout << "ranks = " << cglbm::mpi::Environment::size() << "\n";
        std::cout << "dims = " << topology.dims(0) << "," << topology.dims(1) << "\n";
        std::cout << "depth = " << depth << "\n";
        std::cout << "checked = " << total_checked << "\n";
        std::cout << "untouched = " << total_untouched << "\n";
        std::cout << "wrong = " << total_wrong << "\n";
        std::cout << "result = " << (total_wrong == 0 ? "PASS" : "FAIL") << std::endl;
    }

    return total_wrong == 0 ? 0 : 1;
}
