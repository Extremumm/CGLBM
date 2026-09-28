"""A capillary wave on a heavy layer, against the exact viscous normal mode.

The Laplace cases hold a droplet still; this one moves the interface. The
lower interface of a heavy band is displaced by one cosine wavelength
(lambda = 64) and released, and its first Fourier coefficient rings down. The
reference is the normal mode of two viscous fluids separated by a sharp
interface of tension sigma, `pycglbm.normal_modes.capillary_wave` -- not the
weak-damping rate `2 k^2 (mu1 + mu2) / (rho1 + rho2)`, which misses the
vortical layer at the interface by 15 % here.

The case runs at density ratios of 100 and 1000, each with the published
nine-point source derivative and with `--source-stencil=matched`. What it
measured, as decay rate and frequency over the exact mode's:

                  isotropic            matched
    ratio 100:    1.50 and 0.989       1.16 and 0.980
    ratio 1000:   3.73 and 0.903       2.65 and 0.882

At 100 the error is the heavy fluid's extensional viscosity, which the
potential flow `exp(k y)` inside the layer feels and the matched stencil
corrects. At 1000 it is the interface itself: the light side of a moving
interface is the colour-gradient model's structural limit (docs/numerics.md,
"The heavy fluid's extensional viscosity"), and the frequency is off as well.
The velocity-based solver runs the same case in `capillary_wave_vb`.

These are known gaps, so the values are pinned as measured, with a tolerance
that catches a change in the scheme, and asserted separately are the things
that must hold whatever the gap: the run stays bounded, and the signal is one
decaying mode.
"""

import numpy as np
import pytest
from pycglbm.normal_modes import capillary_wave
from pycglbm.oscillation import fit_damped_oscillation
from pycglbm.testing import artifacts_dir, run_program

#: Steps skipped before the fit: the wave is released with no pressure field
#: under it, and the acoustic start-up that builds one lasts a few hundred.
SETTLING_STEPS = 500

#: Dynamic viscosity of the heavy layer at each density ratio; the light
#: fluid's is 0.05. Both give tau between 0.65 and 6.5.
HEAVY_VISCOSITY = {"100": 0.5, "1000": 2.0}

#: Run length: four periods at 100, two at 1000.
STEPS = {"100": 16000, "1000": 25000}

#: Measured decay rate and angular frequency, over the exact mode's.
MEASURED = {
    ("100", "isotropic"): {"damping": 1.497, "frequency": 0.989},
    ("100", "matched"): {"damping": 1.160, "frequency": 0.980},
    ("1000", "isotropic"): {"damping": 3.727, "frequency": 0.903},
    ("1000", "matched"): {"damping": 2.645, "frequency": 0.882},
}
DAMPING_TOLERANCE = 0.05
FREQUENCY_TOLERANCE = 0.01


def arguments(ratio: str, stencil: str) -> tuple[str, ...]:
    nu = HEAVY_VISCOSITY[ratio] / float(ratio)
    return (
        f"--rho1={ratio}",
        f"--nu={nu!r}",
        f"--nu-b={nu!r}",
        f"--steps={STEPS[ratio]}",
        f"--source-stencil={stencil}",
    )


@pytest.fixture(scope="module", params=sorted(MEASURED))
def wave_run(request):
    """One shared run per density ratio and source stencil."""
    ratio, stencil = request.param
    run = run_program(
        "capillary_wave",
        artifacts_dir() / f"capillary_wave_{ratio}_{stencil}",
        args=arguments(ratio, stencil),
        timeout=2400,
    )
    run.key = request.param
    return run


@pytest.fixture(scope="module")
def result(wave_run):
    """The fitted mode, the exact one, and the track they come from."""
    rho1, rho2 = wave_run.parameter("rho1"), wave_run.parameter("rho2")
    exact = capillary_wave(
        2 * np.pi / wave_run.parameter("nx"),
        rho1,
        rho2,
        rho1 * wave_run.parameter("nu"),
        rho2 * wave_run.parameter("nu2"),
        wave_run.parameter("sigma"),
    )
    timesteps, signal = wave_run.mode_track()
    keep = timesteps >= SETTLING_STEPS
    fit = fit_damped_oscillation(signal[keep], dt=timesteps[1] - timesteps[0])
    return {
        "exact": exact,
        "fit": fit,
        "signal": signal,
        "measured": MEASURED[wave_run.key],
    }


@pytest.mark.long
@pytest.mark.verification
def test_verification_capillary_wave_runs_the_intended_case(wave_run):
    """The options reached the solver: a flag lost on the way is a different case."""
    ratio, stencil = wave_run.key
    assert wave_run.parameter("rho1") / wave_run.parameter("rho2") == float(ratio)
    assert wave_run.config["source_stencil"] == stencil
    assert wave_run.config["viscosity_mixing"] == "dynamic"
    assert wave_run.config["surface_tension"] == "stress"
    assert wave_run.parameter("steps", int) == STEPS[ratio]


@pytest.mark.long
@pytest.mark.verification
def test_verification_capillary_wave_stays_bounded(wave_run, result):
    """No NaN in the track or the last fields, and the phase field within [-1, 1]."""
    assert np.isfinite(result["signal"]).all()
    t = wave_run.last_timestep
    for field in wave_run.fields(t).values():
        assert np.isfinite(field).all()
    phase = wave_run.phase(t)
    assert phase.max() <= 1.0 + 1.0e-6
    assert phase.min() >= -1.0 - 1.0e-6


@pytest.mark.long
@pytest.mark.verification
def test_verification_capillary_wave_is_one_mode(result):
    """What is fitted is a single decaying sinusoid, not a mixture."""
    fit = result["fit"]
    assert fit.residual < 0.03 * fit.amplitude


@pytest.mark.long
@pytest.mark.validation
def test_validation_capillary_wave_against_the_normal_mode(result):
    """Decay rate and frequency over the exact mode's, at their measured values."""
    exact, fit, measured = result["exact"], result["fit"], result["measured"]
    damping = fit.decay_rate / exact.decay_rate
    frequency = fit.angular_frequency / exact.angular_frequency
    assert damping == pytest.approx(measured["damping"], abs=DAMPING_TOLERANCE)
    assert frequency == pytest.approx(measured["frequency"], abs=FREQUENCY_TOLERANCE)
