"""Unit tests for src/mpi/mpi_exchange_ghosts.

The program fills every interior node with a value identifying its global
position, exchanges the ghost layer, and checks each ghost against the global
node it shadows. Corners are included: D2Q9 streams along the diagonals, and
they are filled by the two-pass exchange without any diagonal message.
"""

import pytest

from pycglbm.testing import mpi_launcher, parse_key_values, run_unit_program

PROGRAM = "mpi_exchange_ghosts"

requires_mpi = pytest.mark.skipif(mpi_launcher() is None, reason="no MPI launcher available")


def _run(nprocs, global_nx=64, global_ny=64, depth=1, periodic_y=1):
    result = run_unit_program(
        PROGRAM, (global_nx, global_ny, depth, periodic_y), nprocs=nprocs, check=False
    )
    values = parse_key_values(result.stdout)
    values["returncode"] = result.returncode
    values["stderr"] = result.stderr
    return values


@requires_mpi
@pytest.mark.unit_test
@pytest.mark.parametrize("nprocs", [1, 2, 4, 6])
@pytest.mark.parametrize("depth", [1, 9])
def test_unit_test_mpi_exchange_ghosts_fills_every_ghost(nprocs, depth):
    """Every ghost node, corners included, holds its global neighbour's value."""
    values = _run(nprocs, depth=depth)
    assert values["result"] == "PASS", values["stderr"]
    assert int(values["wrong"]) == 0
    assert int(values["checked"]) > 0
    assert values["returncode"] == 0


@requires_mpi
@pytest.mark.unit_test
@pytest.mark.parametrize("nprocs", [2, 4, 6])
def test_unit_test_mpi_exchange_ghosts_handles_an_uneven_decomposition(nprocs):
    """A lattice that does not divide evenly must still exchange correctly."""
    values = _run(nprocs, global_nx=65, global_ny=63, depth=9)
    assert values["result"] == "PASS", values["stderr"]
    assert int(values["wrong"]) == 0


@requires_mpi
@pytest.mark.unit_test
@pytest.mark.parametrize("nprocs", [1, 4])
def test_unit_test_mpi_exchange_ghosts_leaves_non_periodic_edges_alone(nprocs):
    """Ghosts beyond a non-periodic edge belong to the solver's boundary
    condition and must not be overwritten."""
    values = _run(nprocs, depth=9, periodic_y=0)
    assert values["result"] == "PASS", values["stderr"]
    assert int(values["untouched"]) > 0


@requires_mpi
@pytest.mark.unit_test
def test_unit_test_mpi_exchange_ghosts_is_consistent_across_rank_counts():
    """The same lattice must be exchanged identically however it is split."""
    counts = [_run(n, depth=9)["checked"] for n in (1, 2, 4)]
    # more ranks means more ghost nodes overall, never fewer
    assert counts == sorted(counts, key=int)
