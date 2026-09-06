#ifndef CGLBM_MPI_MPI_EXCHANGE_GHOSTS_H
#define CGLBM_MPI_MPI_EXCHANGE_GHOSTS_H

namespace cglbm {
namespace mpi {

class CartesianTopology;

/// Fill the ghost layer of a scalar field from the neighbouring ranks.
///
/// `field` holds `(nx + 2) * (ny + 2)` values, one ghost node on each side, and
/// is indexed `field[i * (ny + 2) + j]` with `i` along x and `j` along y — the
/// same order as the `[Lx][Ly]` arrays of the solvers, so `j` is contiguous.
/// Interior nodes are `i, j` in `1 .. nx` and `1 .. ny`.
///
/// The exchange runs in two passes, x then y, and the second pass sends the
/// ghost columns it has just received along with the interior. The corner ghost
/// nodes, which D2Q9 needs for its diagonal velocities, are therefore filled by
/// the y pass without any diagonal message: what arrives at a corner has
/// travelled through the neighbour that sits between the two ranks.
///
/// A non-periodic edge leaves its ghost layer untouched, for the solver's own
/// boundary condition — bounce-back, say — to fill.
void exchange_ghosts(const CartesianTopology& topology, double* field, int nx, int ny);

/// The same exchange for `depth` values per node, stored contiguously, as in
/// the `[Lx][Ly][Q]` distribution arrays.
///
/// `field` holds `(nx + 2) * (ny + 2) * depth` values, indexed
/// `field[(i * (ny + 2) + j) * depth + q]`.
void exchange_ghosts(const CartesianTopology& topology,
                     double* field,
                     int nx,
                     int ny,
                     int depth);

}  // namespace mpi
}  // namespace cglbm

#endif  // CGLBM_MPI_MPI_EXCHANGE_GHOSTS_H
