// Prints the decomposition each rank receives, so that the pytest beside this
// file can check that the blocks tile the global lattice exactly once and that
// the neighbours are consistent.
//
//   main_mpi_topology [global_nx] [global_ny] [periodic_x] [periodic_y]

#include <cstdlib>
#include <iostream>

#include "mpi/mpi_environment.h"
#include "mpi/mpi_topology.h"

int main(int argc, char** argv) {
    cglbm::mpi::Environment environment(argc, argv);

    const int global_nx = (argc > 1) ? std::atoi(argv[1]) : 128;
    const int global_ny = (argc > 2) ? std::atoi(argv[2]) : 128;
    const bool periodic_x = (argc > 3) ? std::atoi(argv[3]) != 0 : true;
    const bool periodic_y = (argc > 4) ? std::atoi(argv[4]) != 0 : true;

    cglbm::mpi::CartesianTopology topology(global_nx, global_ny, periodic_x, periodic_y);
    const cglbm::mpi::Decomposition& local = topology.local();

    // One self-describing line per rank; the order between ranks is not
    // deterministic, so the test sorts them.
    std::cout << "rank = " << cglbm::mpi::Environment::rank()
              << " size = " << cglbm::mpi::Environment::size()
              << " dims = " << topology.dims(0) << "," << topology.dims(1)
              << " coords = " << topology.coords(0) << "," << topology.coords(1)
              << " nx = " << local.nx << " ny = " << local.ny
              << " x_offset = " << local.x_offset << " y_offset = " << local.y_offset
              << " left = " << topology.neighbour(-1, 0)
              << " right = " << topology.neighbour(+1, 0)
              << " below = " << topology.neighbour(0, -1)
              << " above = " << topology.neighbour(0, +1)
              << " lower_left = " << topology.neighbour(-1, -1)
              << " upper_right = " << topology.neighbour(+1, +1) << std::endl;

    if (cglbm::mpi::Environment::is_root()) {
        std::cerr << topology.describe() << std::endl;
    }
    return 0;
}
