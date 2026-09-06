#include "mpi/mpi_topology.h"

#include <stdexcept>
#include <string>

#include "mpi/mpi_environment.h"
#include "mpi/mpi_error.h"

#ifdef CGLBM_WITH_MPI
#include <mpi.h>
#endif

namespace cglbm {
namespace mpi {

int block_size(int total, int parts, int index, int* offset) {
    const int base = total / parts;
    const int remainder = total % parts;
    // the first `remainder` blocks take one extra node
    const int size = base + (index < remainder ? 1 : 0);
    if (offset != nullptr) {
        *offset = index * base + (index < remainder ? index : remainder);
    }
    return size;
}

CartesianTopology::CartesianTopology(int global_nx,
                                     int global_ny,
                                     bool periodic_x,
                                     bool periodic_y,
                                     int dims_x,
                                     int dims_y)
    : global_nx_(global_nx), global_ny_(global_ny) {
    if (global_nx < 1 || global_ny < 1) {
        throw std::invalid_argument("CartesianTopology: the lattice must not be empty");
    }
    periodic_[0] = periodic_x;
    periodic_[1] = periodic_y;

#ifdef CGLBM_WITH_MPI
    const int ranks = Environment::size();

    if (dims_x > 0 && dims_y > 0) {
        if (dims_x * dims_y != ranks) {
            throw std::invalid_argument("CartesianTopology: dims_x * dims_y must equal the "
                                        "number of ranks (" +
                                        std::to_string(ranks) + ")");
        }
        dims_[0] = dims_x;
        dims_[1] = dims_y;
    } else {
        // 0 lets MPI choose a balanced factorisation of the rank count
        dims_[0] = 0;
        dims_[1] = 0;
        CGLBM_MPI_CHECK(MPI_Dims_create(ranks, 2, dims_));
    }

    if (dims_[0] > global_nx || dims_[1] > global_ny) {
        throw std::invalid_argument("CartesianTopology: more ranks than lattice nodes along an "
                                    "axis; use fewer ranks or a larger lattice");
    }

    const int periods[2] = {periodic_x ? 1 : 0, periodic_y ? 1 : 0};
    MPI_Comm* comm = new MPI_Comm(MPI_COMM_NULL);
    // reorder = 1 lets MPI renumber the ranks to match the hardware
    CGLBM_MPI_CHECK(MPI_Cart_create(MPI_COMM_WORLD, 2, dims_, periods, 1, comm));
    communicator_ = comm;

    // the rank may have been renumbered by the reorder above
    int cart_rank = 0;
    CGLBM_MPI_CHECK(MPI_Comm_rank(*comm, &cart_rank));
    CGLBM_MPI_CHECK(MPI_Cart_coords(*comm, cart_rank, 2, coords_));
#else
    (void) dims_x;
    (void) dims_y;
    dims_[0] = 1;
    dims_[1] = 1;
    coords_[0] = 0;
    coords_[1] = 0;
#endif

    local_.nx = block_size(global_nx_, dims_[0], coords_[0], &local_.x_offset);
    local_.ny = block_size(global_ny_, dims_[1], coords_[1], &local_.y_offset);
}

CartesianTopology::~CartesianTopology() {
#ifdef CGLBM_WITH_MPI
    if (communicator_ != nullptr) {
        MPI_Comm* comm = static_cast<MPI_Comm*>(communicator_);
        if (*comm != MPI_COMM_NULL) {
            MPI_Comm_free(comm);
        }
        delete comm;
        communicator_ = nullptr;
    }
#endif
}

int CartesianTopology::dims(int axis) const { return dims_[axis]; }

int CartesianTopology::coords(int axis) const { return coords_[axis]; }

bool CartesianTopology::periodic(int axis) const { return periodic_[axis]; }

int CartesianTopology::neighbour(int step_x, int step_y) const {
#ifdef CGLBM_WITH_MPI
    int target[2] = {coords_[0] + step_x, coords_[1] + step_y};
    for (int axis = 0; axis < 2; ++axis) {
        if (target[axis] < 0 || target[axis] >= dims_[axis]) {
            if (!periodic_[axis]) {
                return MPI_PROC_NULL;
            }
            // wrap: (v % n + n) % n keeps the result non-negative
            target[axis] = (target[axis] % dims_[axis] + dims_[axis]) % dims_[axis];
        }
    }
    int rank = MPI_PROC_NULL;
    CGLBM_MPI_CHECK(MPI_Cart_rank(*static_cast<MPI_Comm*>(communicator_), target, &rank));
    return rank;
#else
    // one rank owning everything: it is its own neighbour where the axis wraps
    for (int axis = 0; axis < 2; ++axis) {
        const int step = (axis == 0) ? step_x : step_y;
        if (step != 0 && !periodic_[axis]) {
            return -1;
        }
    }
    return 0;
#endif
}

std::string CartesianTopology::describe() const {
    return "lattice " + std::to_string(global_nx_) + "x" + std::to_string(global_ny_) +
           " over ranks " + std::to_string(dims_[0]) + "x" + std::to_string(dims_[1]) +
           ", this rank at (" + std::to_string(coords_[0]) + "," + std::to_string(coords_[1]) +
           ") owns " + std::to_string(local_.nx) + "x" + std::to_string(local_.ny) +
           " from (" + std::to_string(local_.x_offset) + "," +
           std::to_string(local_.y_offset) + ")";
}

}  // namespace mpi
}  // namespace cglbm
