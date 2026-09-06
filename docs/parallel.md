# Parallel modules

`src/omp` and `src/mpi` hold the parallelism, kept apart from the physics. Both
are written so that **the same sources build whether or not the dependency is
present**: with `WITH_OpenMP=OFF` the OpenMP module reports one thread, with
`WITH_MPI=OFF` the MPI module behaves as a single rank owning the whole
lattice. A solver written against these interfaces needs no `#ifdef` of its own.

| Option | Effect when OFF |
|---|---|
| `WITH_OpenMP` | `cglbm::omp` reports one thread; `#pragma omp` is ignored by the compiler |
| `WITH_MPI` | `cglbm::mpi` reports rank 0 of 1; the halo exchange wraps periodically in place |

`CGLBM_WITH_MPI` is the guard macro, defined by `cmake/FindDependencies.cmake`
when MPI is both requested and found. OpenMP uses the standard `_OPENMP`.

## `src/omp` — thread-level parallelism

`omp_environment.h` wraps the OpenMP runtime calls a solver needs:

```cpp
#include "omp/omp_environment.h"

cglbm::omp::set_thread_count(0);              // 0: OMP_NUM_THREADS, else all cores
std::cout << cglbm::omp::describe() << "\n";  // "OpenMP enabled, 8 threads (of 24 cores)"

#pragma omp parallel
{
    const int id = cglbm::omp::thread_id();
    const int n = cglbm::omp::thread_count();
}
```

`set_thread_count(n)` with `n >= 1` forces `n` threads; with `n < 1` it restores
the runtime default, which is `OMP_NUM_THREADS` when that is set and one thread
per core otherwise. `wall_time()` gives `omp_get_wtime()`, falling back to
`std::chrono::steady_clock`.

`rayleigh_taylor_omp` uses this module, which is why it links in a build without
OpenMP and why `OMP_NUM_THREADS` now controls it — it previously hard-coded
eight threads.

## `src/mpi` — distributed-memory parallelism

Four units, following the same split as arbalete's `src/mpi`:

| Unit | Role |
|---|---|
| `mpi_environment` | RAII init/finalize, rank, size, barrier, reductions |
| `mpi_topology` | 2D Cartesian decomposition of the lattice, neighbours, local extents |
| `mpi_exchange_ghosts` | halo exchange of a field with one ghost layer |
| `mpi_error` | `CGLBM_MPI_CHECK`, which reports the failing call, file and line |

### Decomposition

```cpp
#include "mpi/mpi_environment.h"
#include "mpi/mpi_topology.h"

int main(int argc, char** argv) {
    cglbm::mpi::Environment environment(argc, argv);   // finalises on scope exit

    // periodic along x, walls along y: the boundary conditions of the
    // capillary and Rayleigh-Taylor cases
    cglbm::mpi::CartesianTopology topology(Lx, Ly, /*periodic_x=*/true, /*periodic_y=*/false);
    const auto& local = topology.local();
    // local.nx, local.ny        interior nodes owned here
    // local.x_offset, y_offset  where they sit in the global lattice
}
```

The rank grid comes from `MPI_Dims_create` unless it is given explicitly. Blocks
are balanced: with a lattice that does not divide evenly, the remainder is
spread one node at a time, so no two blocks differ by more than one node.
`neighbour(dx, dy)` takes steps in `{-1, 0, +1}` and resolves all eight
neighbours; it returns `MPI_PROC_NULL` at a non-periodic edge.

### Halo exchange

A field is stored with one ghost layer, `(nx + 2) * (ny + 2)` values indexed
`field[i * (ny + 2) + j]` — the same order as the solvers' `[Lx][Ly]` arrays, so
`j` is contiguous. Distributions carry a third index:
`field[(i * (ny + 2) + j) * Q + q]`.

```cpp
cglbm::mpi::exchange_ghosts(topology, field, local.nx, local.ny);      // scalar
cglbm::mpi::exchange_ghosts(topology, f, local.nx, local.ny, Q);       // distributions
```

The exchange runs **x first, then y**, and the y pass sends the ghost columns it
has just received along with the interior. The corner ghosts that D2Q9 needs for
its diagonal velocities are therefore filled without any diagonal message: what
reaches a corner has travelled through the rank sitting between the two. A
column is contiguous and goes as-is; a row is strided and travels through an
`MPI_Type_vector`.

Ghosts beyond a non-periodic edge are left untouched, for the solver's own
boundary condition — bounce-back, say — to fill.

### Running

```bash
mpirun -n 4 bin/unit_testing/mpi/topology/mpi_topology_opt 128 128
pytest programs/unit_testing/mpi -m unit_test
```

`pycglbm.testing.run_unit_program(..., nprocs=4)` wraps the launcher for the
test suite, and `run_program(..., nprocs=4)` does the same for a solver.

## Status

The modules are tested on their own — the decomposition is checked to tile the
lattice exactly once, and every ghost node, corners included, is checked against
the global node it shadows, on 1, 2, 4 and 6 ranks with both an even and an
uneven split. **The solvers are not yet distributed**: they still declare the
whole lattice as static `[Lx][Ly][Q]` arrays and run on one rank. Making one of
them use `src/mpi` means replacing those globals with a
`Decomposition`-sized allocation, adding a ghost layer, and calling
`exchange_ghosts` after each `stream()` — the modules are the prerequisite, not
the change itself.
