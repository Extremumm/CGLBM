"""A capillary wave on a heavy layer with the velocity-based solver, against the
exact viscous normal mode.

The same case as the colour-gradient `capillary_wave` (see its tests for the
geometry and the reference): a wave of wavelength 64 on the lower interface of
a heavy band, released and ringing down, scored against
`pycglbm.normal_modes.capillary_wave`. The colour-gradient solver is off by
1.5 times in damping at a density ratio of 100 and by 2.6 to 3.7 times at 1000,
with the frequency 10 % low; this is the solver that should get it right, and
does to within 1 % at both.

The measured values are pinned, so a change in the scheme shows up; the
closeness to the exact mode is asserted separately, with room for the
resolution: the interface is four nodes wide on a wavelength of 64.

A third run starts the same wave at 1000 with an amplitude of 8 nodes instead
of 0.3, for one period. Its interface moves at about 4e-3, thirteen times the
speed at which the colour-gradient solver's heavy populations go negative;
that solver's wave, with the matched stencil, runs at this amplitude only
since its equation of state clamps the phase field, and is damped 1.87 times
the linear rate. This one stays bounded, closer to it. At `k a = 0.79` the
wave is no longer linear, so its frequency is held only to the value it
measured, which is 9 % below the linear mode's.

A fourth runs the wave at a density ratio of 1e4, with the same dynamic
viscosities, for one period. The heavy fluid's oscillatory boundary layer is
then 1.6 nodes thick, and the collision's non-equilibrium, filtered in time
rather than rebuilt from the finite-difference velocity gradient, is what
reads it: the gradient took the damping to 2.10 times the exact rate.

A fifth runs the wave at 1000 on a wavelength of 128, which resolves the heavy
fluid's boundary layer, with the phase populations built to fourth order
(`--fourth-order-phase`, SolverParameters::fourth_order_phase). Before its
normal was confined to the interface, this wave grew a row of cells in the
heavy fluid seven to eleven nodes from it, and its interface's harmonics rose
from 8e-6 at 5000 steps to 4e-3 at 25 000; the linear runs are held to a
clean interface.
"""

import numpy as np
import pytest
from pycglbm.normal_modes import capillary_wave
from pycglbm.oscillation import fit_damped_oscillation
from pycglbm.testing import artifacts_dir, run_program

#: Steps skipped before the fit: the acoustic start-up of the release.
SETTLING_STEPS = 500

#: Density ratio, dynamic viscosity of the heavy layer, run length and initial
#: amplitude (nodes) of each run, and for some the wavelength (64 otherwise) and
#: the fourth-order phase.
CASES = {
    "100": {"ratio": "100", "mu1": 0.5, "steps": 16000, "amplitude": 0.3},
    "1000": {"ratio": "1000", "mu1": 2.0, "steps": 25000, "amplitude": 0.3},
    "1000_large": {"ratio": "1000", "mu1": 2.0, "steps": 12500, "amplitude": 8.0},
    "1e4": {"ratio": "1e4", "mu1": 2.0, "steps": 40000, "amplitude": 0.3},
    "1000_128_fourth_order": {
        "ratio": "1000",
        "mu1": 2.0,
        "steps": 25000,
        "amplitude": 0.3,
        "wavelength": 128,
        "fourth_order_phase": True,
    },
}

#: Measured decay rate and angular frequency, over the exact linear mode's.
MEASURED = {
    "100": {"damping": 1.001, "frequency": 0.991},
    "1000": {"damping": 0.991, "frequency": 0.992},
    "1000_large": {"damping": 1.024, "frequency": 0.909},
    "1e4": {"damping": 1.268, "frequency": 0.992},
    "1000_128_fourth_order": {"damping": 0.976, "frequency": 0.999},
}
MEASURED_TOLERANCE = 0.02

#: How close to the exact mode the solver has to be, whatever it measured. At
#: 1e4 the heavy fluid's oscillatory boundary layer is 1.6 nodes thick.
DAMPING_BOUND = {
    "100": 0.15,
    "1000": 0.15,
    "1000_large": 0.15,
    "1e4": 0.35,
    "1000_128_fourth_order": 0.15,
}
FREQUENCY_BOUND = {
    "100": 0.03,
    "1000": 0.03,
    "1000_large": 0.12,
    "1e4": 0.03,
    "1000_128_fourth_order": 0.03,
}

#: The linear runs whose interface should hold no harmonic above the second,
#: and the bound on them, over the initial amplitude. At 1e4 the third follows
#: the wave's own swing, to 2e-4.
CLEAN_INTERFACE = ("100", "1000", "1000_128_fourth_order")
HARMONIC_BOUND = 1.0e-4


@pytest.fixture(scope="module", params=sorted(CASES))
def wave_run(request):
    """One shared run per case."""
    case = CASES[request.param]
    options = (f"--wavelength={case.get('wavelength', 64)}",)
    if case.get("fourth_order_phase", False):
        options += ("--fourth-order-phase",)
    run = run_program(
        "capillary_wave_vb",
        artifacts_dir() / f"capillary_wave_vb_{request.param}",
        args=(
            "E8",
            case["ratio"],
            str(case["mu1"]),
            str(case["steps"]),
            str(case["amplitude"]),
            *options,
        ),
        timeout=3600,
    )
    run.key = request.param
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
        "key": wave_run.key,
    }


def initial_amplitude(wave_run):
    """The first Fourier coefficient of the lower interface at t = 0."""
    return wave_run.mode_track()[1][0]


def interface_harmonics(wave_run, timestep):
    """Fourier amplitudes of the lower interface's height, as mode.csv reads
    the first: from the heavy fluid's volume in each column of the lower half."""
    fraction = 0.5 * (1.0 + wave_run.phase(timestep))
    ny, nx = fraction.shape
    height = ny / 2 - fraction[: ny // 2].sum(axis=0)
    return np.abs(np.fft.rfft(height)) * 2 / nx


@pytest.mark.long
@pytest.mark.verification
def test_verification_capillary_wave_vb_runs_the_intended_case(wave_run):
    """The arguments reached the solver."""
    case = CASES[wave_run.key]
    assert wave_run.parameter("rho1") / wave_run.parameter("rho2") == float(case["ratio"])
    assert wave_run.parameter("mu1") == case["mu1"]
    assert wave_run.parameter("steps", int) == case["steps"]
    assert wave_run.parameter("amplitude") == case["amplitude"]
    assert wave_run.parameter("nx", int) == case.get("wavelength", 64)
    assert wave_run.parameter("fourth_order_phase", int) == case.get("fourth_order_phase", False)
    assert initial_amplitude(wave_run) == pytest.approx(case["amplitude"], rel=1.0e-3)


@pytest.mark.long
@pytest.mark.verification
def test_verification_capillary_wave_vb_stays_bounded(wave_run, result):
    """No NaN in the track or the last fields, and the phase field in its bounds."""
    assert np.isfinite(result["signal"]).all()
    for field in wave_run.fields(wave_run.last_timestep).values():
        assert np.isfinite(field).all()
    phase = wave_run.phase(wave_run.last_timestep)
    assert phase.min() >= -1.0 - 1.0e-9
    assert phase.max() <= 1.0 + 1.0e-9


@pytest.mark.long
@pytest.mark.verification
def test_verification_capillary_wave_vb_keeps_a_clean_interface(wave_run):
    """The linear runs end with no harmonic above the second standing out."""
    if wave_run.key not in CLEAN_INTERFACE:
        pytest.skip("nonlinear, or at 1e4 where the third harmonic follows the wave")
    harmonics = interface_harmonics(wave_run, wave_run.last_timestep)
    assert harmonics[3:].max() < HARMONIC_BOUND * CASES[wave_run.key]["amplitude"]


@pytest.mark.long
@pytest.mark.verification
def test_verification_capillary_wave_vb_is_one_mode(result):
    """What is fitted is a single decaying sinusoid, not a mixture."""
    assert result["fit"].residual < 0.03 * result["fit"].amplitude


@pytest.mark.long
@pytest.mark.validation
def test_validation_capillary_wave_vb_against_the_normal_mode(result):
    """Decay rate and frequency close to the exact mode's, and at their measured values."""
    assert result["damping"] == pytest.approx(1.0, abs=DAMPING_BOUND[result["key"]])
    assert result["frequency"] == pytest.approx(1.0, abs=FREQUENCY_BOUND[result["key"]])
    assert result["damping"] == pytest.approx(result["measured"]["damping"], abs=MEASURED_TOLERANCE)
    assert result["frequency"] == pytest.approx(
        result["measured"]["frequency"], abs=MEASURED_TOLERANCE / 4
    )
