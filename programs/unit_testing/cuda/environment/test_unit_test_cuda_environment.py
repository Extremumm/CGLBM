"""Unit tests for src/cuda: device detection, memory management, and environment queries."""

import pytest
from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "cuda_environment"


def _run(*args, env=None):
    return parse_key_values(run_unit_program(PROGRAM, args, env=env).stdout)


def _require_cuda(values):
    """Skip when the binary was built without CUDA or no GPU is detected."""
    if values["available"] != "1":
        pytest.skip("built without CUDA or no CUDA device available")


@pytest.mark.unit_test
def test_unit_test_cuda_environment_describes_the_build():
    values = _run()
    assert values["available"] in ("0", "1")
    expected = "CUDA enabled" if values["available"] == "1" else "CUDA disabled"
    assert expected in values["describe"]


@pytest.mark.unit_test
def test_unit_test_cuda_environment_device_properties():
    values = _run()
    _require_cuda(values)
    assert int(values["device_count"]) >= 1
    assert len(values["device_name"]) > 0
    assert float(values["compute_capability"]) > 0.0
    assert int(values["total_memory_mb"]) > 0


@pytest.mark.unit_test
def test_unit_test_cuda_environment_memory_roundtrip():
    values = _run()
    _require_cuda(values)
    assert values["memory_roundtrip"] == "1"


@pytest.mark.unit_test
def test_unit_test_cuda_environment_wall_time_is_monotonic():
    assert _run()["wall_time_monotonic"] == "1"
