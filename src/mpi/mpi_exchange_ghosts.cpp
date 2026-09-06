#include "mpi/mpi_exchange_ghosts.h"

#include "mpi/mpi_error.h"
#include "mpi/mpi_topology.h"

#ifdef CGLBM_WITH_MPI
#include <mpi.h>
#endif

namespace cglbm {
namespace mpi {

namespace {

#ifndef CGLBM_WITH_MPI
/// Serial stand-in for one exchange: a periodic axis copies its own opposite
/// interior layer into the ghost layer, which is what a single periodic rank
/// would have received from itself.
void wrap_periodic(double* field, int nx, int ny, int depth, bool periodic_x, bool periodic_y) {
    const int stride_j = depth;
    const int stride_i = (ny + 2) * depth;

    if (periodic_x) {
        for (int j = 0; j < ny + 2; ++j) {
            for (int q = 0; q < depth; ++q) {
                field[0 * stride_i + j * stride_j + q] = field[nx * stride_i + j * stride_j + q];
                field[(nx + 1) * stride_i + j * stride_j + q] =
                    field[1 * stride_i + j * stride_j + q];
            }
        }
    }
    if (periodic_y) {
        for (int i = 0; i < nx + 2; ++i) {
            for (int q = 0; q < depth; ++q) {
                field[i * stride_i + 0 * stride_j + q] = field[i * stride_i + ny * stride_j + q];
                field[i * stride_i + (ny + 1) * stride_j + q] =
                    field[i * stride_i + 1 * stride_j + q];
            }
        }
    }
}
#endif

}  // namespace

void exchange_ghosts(const CartesianTopology& topology, double* field, int nx, int ny, int depth) {
    if (field == nullptr || nx < 1 || ny < 1 || depth < 1) {
        return;
    }

#ifdef CGLBM_WITH_MPI
    MPI_Comm comm = *static_cast<MPI_Comm*>(topology.communicator());

    const int stride_j = depth;             // one node along y
    const int stride_i = (ny + 2) * depth;  // one column along x
    const int tag_x = 0;
    const int tag_y = 1;

    // --- x direction -------------------------------------------------------
    // A column i = const is contiguous: (ny + 2) * depth values in a row. Only
    // the interior rows are meaningful at this point, but sending the whole
    // column keeps the transfer contiguous and costs two nodes.
    const int left = topology.neighbour(-1, 0);
    const int right = topology.neighbour(+1, 0);
    const int column = (ny + 2) * depth;

    // last interior column -> left ghost column of the right neighbour
    CGLBM_MPI_CHECK(MPI_Sendrecv(field + nx * stride_i,
                                 column,
                                 MPI_DOUBLE,
                                 right,
                                 tag_x,
                                 field + 0 * stride_i,
                                 column,
                                 MPI_DOUBLE,
                                 left,
                                 tag_x,
                                 comm,
                                 MPI_STATUS_IGNORE));
    // first interior column -> right ghost column of the left neighbour
    CGLBM_MPI_CHECK(MPI_Sendrecv(field + 1 * stride_i,
                                 column,
                                 MPI_DOUBLE,
                                 left,
                                 tag_x,
                                 field + (nx + 1) * stride_i,
                                 column,
                                 MPI_DOUBLE,
                                 right,
                                 tag_x,
                                 comm,
                                 MPI_STATUS_IGNORE));

    // --- y direction -------------------------------------------------------
    // A row j = const is strided: (nx + 2) blocks of `depth` values, every
    // stride_i values. The row spans the ghost columns too, so the corners
    // filled above travel on and land in the corner ghosts of the neighbour.
    MPI_Datatype row_type = MPI_DATATYPE_NULL;
    CGLBM_MPI_CHECK(MPI_Type_vector(nx + 2, depth, stride_i, MPI_DOUBLE, &row_type));
    CGLBM_MPI_CHECK(MPI_Type_commit(&row_type));

    const int below = topology.neighbour(0, -1);
    const int above = topology.neighbour(0, +1);

    CGLBM_MPI_CHECK(MPI_Sendrecv(field + ny * stride_j,
                                 1,
                                 row_type,
                                 above,
                                 tag_y,
                                 field + 0 * stride_j,
                                 1,
                                 row_type,
                                 below,
                                 tag_y,
                                 comm,
                                 MPI_STATUS_IGNORE));
    CGLBM_MPI_CHECK(MPI_Sendrecv(field + 1 * stride_j,
                                 1,
                                 row_type,
                                 below,
                                 tag_y,
                                 field + (ny + 1) * stride_j,
                                 1,
                                 row_type,
                                 above,
                                 tag_y,
                                 comm,
                                 MPI_STATUS_IGNORE));

    CGLBM_MPI_CHECK(MPI_Type_free(&row_type));
#else
    wrap_periodic(field, nx, ny, depth, topology.periodic(0), topology.periodic(1));
#endif
}

void exchange_ghosts(const CartesianTopology& topology, double* field, int nx, int ny) {
    exchange_ghosts(topology, field, nx, ny, 1);
}

}  // namespace mpi
}  // namespace cglbm
