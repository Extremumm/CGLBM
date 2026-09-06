#ifndef CGLBM_MPI_MPI_TOPOLOGY_H
#define CGLBM_MPI_MPI_TOPOLOGY_H

#include <string>

namespace cglbm {
namespace mpi {

/// The piece of the global lattice owned by one rank.
///
/// Indices are global lattice nodes: this rank owns the interior nodes
/// `x_offset .. x_offset + nx - 1` and `y_offset .. y_offset + ny - 1`.
struct Decomposition {
    int nx = 0;        ///< interior nodes along x
    int ny = 0;        ///< interior nodes along y
    int x_offset = 0;  ///< global x index of the first interior node
    int y_offset = 0;  ///< global y index of the first interior node

    /// Nodes owned by this rank, ghost layer excluded.
    int size() const {
        return nx * ny;
    }

    /// Extent of the array including one ghost layer on each side.
    int padded_nx() const {
        return nx + 2;
    }
    int padded_ny() const {
        return ny + 2;
    }

    /// Elements of a field stored with one ghost layer.
    int padded_size() const {
        return padded_nx() * padded_ny();
    }
};

/// A two-dimensional Cartesian decomposition of the lattice.
///
/// The ranks are arranged on a `dims(0) x dims(1)` grid; each holds a
/// rectangular block of the global `global_nx x global_ny` lattice, padded with
/// one ghost layer. D2Q9 streams along the diagonals too, so the eight
/// neighbours are all resolved, corners included.
///
/// Without MPI the object still exists: one rank, the whole lattice, and every
/// neighbour is the rank itself.
class CartesianTopology {
public:
    /// Split `global_nx x global_ny` over every rank of `MPI_COMM_WORLD`.
    ///
    /// The rank grid is chosen by MPI (`MPI_Dims_create`) unless `dims_x` and
    /// `dims_y` are both positive, in which case their product must equal the
    /// number of ranks. Periodicity mirrors the boundary conditions of the
    /// solver: the color-gradient cases are periodic along x, and periodic
    /// along y only for `laplace`.
    CartesianTopology(int global_nx,
                      int global_ny,
                      bool periodic_x = true,
                      bool periodic_y = true,
                      int dims_x = 0,
                      int dims_y = 0);
    ~CartesianTopology();

    CartesianTopology(const CartesianTopology&) = delete;
    CartesianTopology& operator=(const CartesianTopology&) = delete;

    /// Ranks along `axis` (0 for x, 1 for y).
    int dims(int axis) const;

    /// Position of this rank on the rank grid, along `axis`.
    int coords(int axis) const;

    /// Whether `axis` wraps around.
    bool periodic(int axis) const;

    /// Rank of the neighbour `step_x` blocks along x and `step_y` along y, each
    /// in {-1, 0, +1}. Returns a negative value when there is no neighbour,
    /// i.e. at a non-periodic edge; `(0, 0)` is this rank.
    int neighbour(int step_x, int step_y) const;

    /// The block owned by this rank.
    const Decomposition& local() const {
        return local_;
    }

    /// Global lattice size along x and y.
    int global_nx() const {
        return global_nx_;
    }
    int global_ny() const {
        return global_ny_;
    }

    /// The Cartesian communicator, as a `void*` so that this header stays free
    /// of `mpi.h`. Null without MPI. Cast it back to `MPI_Comm*` to use it.
    void* communicator() const {
        return communicator_;
    }

    /// One line for the run log.
    std::string describe() const;

private:
    int global_nx_ = 0;
    int global_ny_ = 0;
    int dims_[2] = {1, 1};
    int coords_[2] = {0, 0};
    bool periodic_[2] = {true, true};
    Decomposition local_;
    void* communicator_ = nullptr;
};

/// Split `total` nodes over `parts` blocks, returning the size of block
/// `index` and, through `offset`, where it starts.
///
/// The remainder is spread one node at a time over the first blocks, so no two
/// blocks differ by more than one node.
int block_size(int total, int parts, int index, int* offset);

}  // namespace mpi
}  // namespace cglbm

#endif  // CGLBM_MPI_MPI_TOPOLOGY_H
