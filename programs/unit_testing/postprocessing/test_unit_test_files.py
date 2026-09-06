"""Unit tests for the CSV reader, on a synthetic run directory."""

import numpy as np
import pytest
from pycglbm import CaseOutput
from pycglbm.files import load_velocity

LY, LX = 6, 4


@pytest.fixture
def fake_rundir(tmp_path):
    """A run directory holding two complete timesteps, plus one incomplete."""
    for timestep in (0, 100):
        # rows are y, columns are x, exactly as the solvers write them
        density = np.arange(LY * LX).reshape(LY, LX) + timestep
        np.savetxt(tmp_path / f"density_{timestep}.csv", density, delimiter=",")
        np.savetxt(tmp_path / f"phase_{timestep}.csv", np.zeros((LY, LX)), delimiter=",")
        np.savetxt(tmp_path / f"pressure_{timestep}.csv", np.ones((LY, LX)), delimiter=",")
        # velocity interleaves (ux, uy) along the row, giving 2 * LX columns
        velocity = np.zeros((LY, 2 * LX))
        velocity[:, 0::2] = 1.0  # ux
        velocity[:, 1::2] = -2.0  # uy
        np.savetxt(tmp_path / f"velocity_{timestep}.csv", velocity, delimiter=",")
    # a timestep whose velocity file is missing must not be reported
    np.savetxt(tmp_path / "density_200.csv", np.zeros((LY, LX)), delimiter=",")
    return tmp_path


@pytest.mark.unit_test
def test_unit_test_files_timesteps_are_sorted_and_complete(fake_rundir):
    assert CaseOutput(fake_rundir).timesteps == [0, 100]


@pytest.mark.unit_test
def test_unit_test_files_scalar_field_is_indexed_y_then_x(fake_rundir):
    case = CaseOutput(fake_rundir)
    assert case.shape == (LY, LX)
    assert case.density(100)[0, 0] == pytest.approx(100.0)
    assert case.last_timestep == 100


@pytest.mark.unit_test
def test_unit_test_files_velocity_is_folded_into_two_components(fake_rundir):
    velocity = load_velocity(fake_rundir, 0)
    assert velocity.shape == (LY, LX, 2)
    assert np.all(velocity[..., 0] == 1.0)
    assert np.all(velocity[..., 1] == -2.0)


@pytest.mark.unit_test
def test_unit_test_files_missing_output_is_reported(fake_rundir):
    with pytest.raises(FileNotFoundError):
        CaseOutput(fake_rundir).density(999)
