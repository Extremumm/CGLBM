"""Unit tests for src/omp: the thread count the module reports is the one a
parallel region actually runs on."""

import pytest

from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "omp_environment"


def _run(*args, env=None):
    return parse_key_values(run_unit_program(PROGRAM, args, env=env).stdout)


@pytest.mark.unit_test
def test_unit_test_omp_environment_reports_a_build_with_openmp():
    values = _run()
    assert values["available"] == "1", "the test build is expected to have OpenMP enabled"
    assert int(values["max_threads"]) >= 1
    assert "OpenMP enabled" in values["describe"]


@pytest.mark.unit_test
@pytest.mark.parametrize("requested", [1, 2, 3])
def test_unit_test_omp_environment_honours_the_requested_thread_count(requested):
    values = _run(requested)
    assert int(values["max_threads"]) == requested
    # what the runtime promised is what the parallel region delivered
    assert int(values["observed_threads"]) == requested


@pytest.mark.unit_test
def test_unit_test_omp_environment_thread_ids_are_unique_and_contiguous():
    values = _run(4)
    # 0 + 1 + ... + (n-1): every id appears exactly once
    assert values["id_sum"] == values["expected_id_sum"]


@pytest.mark.unit_test
def test_unit_test_omp_environment_zero_restores_the_default():
    default = int(_run()["max_threads"])
    assert int(_run(0)["max_threads"]) == default


@pytest.mark.unit_test
def test_unit_test_omp_environment_default_follows_omp_num_threads():
    """Asking for the default must give back what OMP_NUM_THREADS requested."""
    values = _run(0, env={"OMP_NUM_THREADS": "3"})
    assert int(values["max_threads"]) == 3
    assert int(values["observed_threads"]) == 3


@pytest.mark.unit_test
def test_unit_test_omp_environment_explicit_count_overrides_the_environment():
    values = _run(2, env={"OMP_NUM_THREADS": "7"})
    assert int(values["max_threads"]) == 2


@pytest.mark.unit_test
def test_unit_test_omp_environment_wall_time_does_not_go_backwards():
    assert _run()["wall_time_monotonic"] == "1"
