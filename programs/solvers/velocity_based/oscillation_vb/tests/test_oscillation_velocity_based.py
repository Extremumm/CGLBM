"""Mode-2 oscillation of a droplet at a density ratio of 1000 with the
velocity-based solver, against the exact viscous normal mode.

The same case as the colour-gradient `oscillation` (see its tests): a droplet
of radius 20 deformed by 3 % into `cos 2 theta`, released and ringing down,
scored against `pycglbm.normal_modes.droplet_mode`. The colour-gradient solver
damps it eight times too fast with its published source term and not at all
with the matched one, because its interface moves wrongly at this density
ratio; this is the solver whose interface is meant to be right.

The measured values are pinned, so a change in the scheme shows up; the
closeness to the exact mode is asserted separately.
"""

import numpy as np
import pytest
from pycglbm.normal_modes import droplet_mode
from pycglbm.oscillation import fit_damped_oscillation
from pycglbm.testing import artifacts_dir, run_program

#: Steps skipped before the fit: the acoustic start-up of the release.
SETTLING_STEPS = 500

#: Measured decay rate and angular frequency, over the exact mode's.
MEASURED = {"damping": 1.082, "frequency": 0.987}
MEASURED_TOLERANCE = 0.02

#: How close to the exact mode the solver has to be, whatever it measured.
DAMPING_BOUND = 0.2
FREQUENCY_BOUND = 0.03


@pytest.fixture(scope="module")
def droplet_run():
    """One shared run of the case."""
    return run_program(
        "oscillation_vb", artifacts_dir() / "oscillation_vb", args=("E8",), timeout=2400
    )


@pytest.fixture(scope="module")
def result(droplet_run):
    """The fitted mode, the exact one, and the track they come from."""
    exact = droplet_mode(
        2,
        droplet_run.parameter("radius"),
        droplet_run.parameter("rho1"),
        droplet_run.parameter("rho2"),
        droplet_run.parameter("mu1"),
        droplet_run.parameter("mu2"),
        droplet_run.parameter("sigma"),
    )
    timesteps, signal = droplet_run.mode_track()
    keep = timesteps >= SETTLING_STEPS
    fit = fit_damped_oscillation(signal[keep], dt=timesteps[1] - timesteps[0])
    return {
        "damping": fit.decay_rate / exact.decay_rate,
        "frequency": fit.angular_frequency / exact.angular_frequency,
        "fit": fit,
        "signal": signal,
    }


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_vb_stays_bounded(droplet_run, result):
    """No NaN in the track or the last fields."""
    assert np.isfinite(result["signal"]).all()
    for field in droplet_run.fields(droplet_run.last_timestep).values():
        assert np.isfinite(field).all()


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_vb_starts_at_the_laid_down_deformation(result):
    """`D = 2 eps` for a small deformation eps, less what the diffuse profile smears."""
    assert result["signal"][0] == pytest.approx(2 * 0.03, rel=0.05)


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_vb_is_one_mode(result):
    """What is fitted is a single sinusoid in an exponential envelope."""
    assert result["fit"].residual < 0.03 * result["fit"].amplitude


@pytest.mark.long
@pytest.mark.validation
def test_validation_oscillation_vb_against_the_normal_mode(result):
    """Decay rate and frequency close to the exact mode's, and at their measured values."""
    assert result["damping"] == pytest.approx(1.0, abs=DAMPING_BOUND)
    assert result["frequency"] == pytest.approx(1.0, abs=FREQUENCY_BOUND)
    assert result["damping"] == pytest.approx(MEASURED["damping"], abs=MEASURED_TOLERANCE)
    assert result["frequency"] == pytest.approx(MEASURED["frequency"], abs=MEASURED_TOLERANCE / 4)
