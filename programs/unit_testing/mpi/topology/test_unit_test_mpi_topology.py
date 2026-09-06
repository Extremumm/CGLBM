"""Unit tests for src/mpi/mpi_topology: the blocks must tile the lattice
exactly once, and the neighbours must agree with the rank grid."""

import pytest

from pycglbm.testing import mpi_launcher, run_unit_program

PROGRAM = "mpi_topology"

requires_mpi = pytest.mark.skipif(mpi_launcher() is None, reason="no MPI launcher available")


def _ranks(nprocs, global_nx=128, global_ny=128, periodic_x=1, periodic_y=1):
    """Run on ``nprocs`` ranks and return one dict of ints per rank, by rank."""
    result = run_unit_program(
        PROGRAM, (global_nx, global_ny, periodic_x, periodic_y), nprocs=nprocs
    )
    ranks = {}
    for line in result.stdout.splitlines():
        if not line.startswith("rank ="):
            continue
        # "rank = 0 size = 4 dims = 2,2 ..." -> {"rank": 0, "dims": (2, 2), ...}
        tokens = line.replace(" = ", "=").split()
        entry = {}
        for token in tokens:
            key, _, value = token.partition("=")
            entry[key] = tuple(int(v) for v in value.split(",")) if "," in value else int(value)
        ranks[entry["rank"]] = entry
    return ranks


@requires_mpi
@pytest.mark.unit_test
@pytest.mark.parametrize("nprocs", [1, 2, 4, 6])
def test_unit_test_mpi_topology_blocks_tile_the_lattice(nprocs):
    global_nx, global_ny = 128, 128
    ranks = _ranks(nprocs, global_nx, global_ny)
    assert len(ranks) == nprocs

    owned = set()
    for entry in ranks.values():
        for i in range(entry["nx"]):
            for j in range(entry["ny"]):
                node = (entry["x_offset"] + i, entry["y_offset"] + j)
                assert node not in owned, f"node {node} owned by two ranks"
                owned.add(node)
    assert len(owned) == global_nx * global_ny


@requires_mpi
@pytest.mark.unit_test
@pytest.mark.parametrize("nprocs", [2, 4, 6])
def test_unit_test_mpi_topology_blocks_are_balanced(nprocs):
    """No two blocks may differ by more than one node along an axis."""
    ranks = _ranks(nprocs, 100, 100)  # 100 divides by neither 3 nor 6
    for axis in ("nx", "ny"):
        sizes = [entry[axis] for entry in ranks.values()]
        assert max(sizes) - min(sizes) <= 1


@requires_mpi
@pytest.mark.unit_test
@pytest.mark.parametrize("nprocs", [1, 4, 6])
def test_unit_test_mpi_topology_neighbours_are_reciprocal(nprocs):
    """The rank to my right must have me on its left."""
    ranks = _ranks(nprocs)
    for rank, entry in ranks.items():
        assert ranks[entry["right"]]["left"] == rank
        assert ranks[entry["above"]]["below"] == rank
        assert ranks[entry["upper_right"]]["lower_left"] == rank


@requires_mpi
@pytest.mark.unit_test
def test_unit_test_mpi_topology_non_periodic_edges_have_no_neighbour():
    """Along a non-periodic axis the outer ranks must report MPI_PROC_NULL."""
    ranks = _ranks(4, 128, 128, periodic_x=1, periodic_y=0)
    dims_y = next(iter(ranks.values()))["dims"][1]
    edge_ranks = [e for e in ranks.values() if e["coords"][1] in (0, dims_y - 1)]
    assert edge_ranks, "expected at least one rank on a y edge"
    for entry in edge_ranks:
        if entry["coords"][1] == 0:
            assert entry["below"] < 0
        if entry["coords"][1] == dims_y - 1:
            assert entry["above"] < 0
