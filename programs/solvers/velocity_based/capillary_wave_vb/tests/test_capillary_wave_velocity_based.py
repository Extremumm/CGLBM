"""A capillary wave on a heavy layer with the velocity-based solver, against the
exact viscous normal mode.

The same case as the colour-gradient `capillary_wave` (see its tests for the
geometry and the reference): a wave of wavelength 64 on the lower interface of
a heavy band, released and ringing down, scored against
`pycglbm.normal_modes.capillary_wave`. The colour-gradient solver is off by
1.5 times in damping at a density ratio of 100 and by 2.6 to 3.7 times at 1000,
with the frequency 10 % low; this is the solver that should get it right, and
does to within a few per cent.

The measured values are pinned, so a change in the scheme shows up; the
closeness to the exact mode is asserted separately, with room for the
resolution: the interface is four nodes wide on a wavelength of 64.
"""

import numpy as np
import pytest
from pycglbm.normal_modes import capillary_wave
from pycglbm.oscillation import fit_damped_oscillation
from pycglbm.testing import artifacts_dir, run_program

#: Steps skipped before the fit: the acoustic start-up of the release.
SETTLING_STEPS = 500

#: Dynamic viscosity of the heavy layer and the run length, per density ratio.
HEAVY_VISCOSITY = {"100": 0.5, "1000": 2.0}
STEPS = {"100": 16000, "1000": 25000}

#: Measured decay rate and angular frequency, over the exact mode's.
MEASURED = {
    "100": {"damping": 1.043, "frequency": 0.988},
    "1000": {"damping": 1.084, "frequency": 0.982},
}
MEASURED_TOLERANCE = 0.02

#: How close to the exact mode the solver has to be, whatever it measured.
DAMPING_BOUND = 0.15
FREQUENCY_BOUND = 0.03


@pytest.fixture(scope="module", params=sorted(MEASURED))
def wave_run(request):
    """One shared run per density ratio."""
    ratio = request.param
    run = run_program(
        "capillary_wave_vb",
        artifacts_dir() / f"capillary_wave_vb_{ratio}",
        args=("E8", ratio, str(HEAVY_VISCOSITY[ratio]), str(STEPS[ratio])),
        timeout=2400,
    )
    run.key = ratio
    return run


@pytest.fixture(scope="module")
def result(wave_run):
    """The fitted mode, the exact one, and the track they come from."""
    exact = capillary_wave(
        2 * np.pi / wave_run.parameter("nx"),
        wave_run.parameter("rho1"),
        wave_run.parameter("rho2"),
        wave_run.parameter("mu1"),
        wave_run.parameter("mu2"),
        wave_run.parameter("sigma"),
    )
    timesteps, signal = wave_run.mode_track()
    keep = timesteps >= SETTLING_STEPS
    fit = fit_damped_oscillation(signal[keep], dt=timesteps[1] - timesteps[0])
    return {
        "damping": fit.decay_rate / exact.decay_rate,
        "frequency": fit.angular_frequency / exact.angular_frequency,
        "fit": fit,
        "signal": signal,
        "measured": MEASURED[wave_run.key],
    }


@pytest.mark.long
@pytest.mark.verification
def test_verification_capillary_wave_vb_runs_the_intended_case(wave_run):
    """The arguments reached the solver."""
    assert wave_run.parameter("rho1") / wave_run.parameter("rho2") == float(wave_run.key)
    assert wave_run.parameter("mu1") == HEAVY_VISCOSITY[wave_run.key]
    assert wave_run.parameter("steps", int) == STEPS[wave_run.key]


@pytest.mark.long
@pytest.mark.verification
def test_verification_capillary_wave_vb_stays_bounded(wave_run, result):
    """No NaN in the track or the last fields."""
    assert np.isfinite(result["signal"]).all()
    for field in wave_run.fields(wave_run.last_timestep).values():
        assert np.isfinite(field).all()


@pytest.mark.long
@pytest.mark.verification
def test_verification_capillary_wave_vb_is_one_mode(result):
    """What is fitted is a single decaying sinusoid, not a mixture."""
    assert result["fit"].residual < 0.03 * result["fit"].amplitude


@pytest.mark.long
@pytest.mark.validation
def test_validation_capillary_wave_vb_against_the_normal_mode(result):
    """Decay rate and frequency close to the exact mode's, and at their measured values."""
    assert result["damping"] == pytest.approx(1.0, abs=DAMPING_BOUND)
    assert result["frequency"] == pytest.approx(1.0, abs=FREQUENCY_BOUND)
    assert result["damping"] == pytest.approx(result["measured"]["damping"], abs=MEASURED_TOLERANCE)
    assert result["frequency"] == pytest.approx(
        result["measured"]["frequency"], abs=MEASURED_TOLERANCE / 4
    )
