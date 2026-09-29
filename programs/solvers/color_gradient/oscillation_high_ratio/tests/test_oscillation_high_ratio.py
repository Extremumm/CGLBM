"""Mode-2 oscillation of a droplet at a density ratio of 1000 in the
two-population solver, against the exact viscous normal mode.

The case of `oscillation` -- a droplet of radius 20 deformed by 3 % into
`cos 2 theta`, released and ringing down, with dynamic viscosities of 2 and
0.05 -- in Ba et al.'s model, with their parameters, MRT collision and
third-moment source term. The reference is
`pycglbm.normal_modes.droplet_mode`: decay rate 1.84e-5 per step and a period
of 23 200 steps for sigma = 0.1.

As in `Solver`, the droplet's interior strain is exact and what the case
measures is the interface in motion. With the published nine-point derivative
of the third-moment source it rings down 6.6 times too fast and 10 % low in
frequency. With `--source-stencil=matched` it rings down at 1.45 times the
exact rate, 2 % high in frequency; before the stencil's face values were
limited at the interface it rang down at 0.23 of it. See docs/numerics.md,
"The matched stencil at an interface".

The values are pinned as measured, with a tolerance that catches a change in
the scheme.
"""

import numpy as np
import pytest
from pycglbm.normal_modes import droplet_mode
from pycglbm.oscillation import fit_damped_oscillation
from pycglbm.testing import artifacts_dir, run_program

#: Steps skipped before the fit: the acoustic start-up of the release.
SETTLING_STEPS = 500

#: Measured decay rate and angular frequency, over the exact mode's.
MEASURED = {
    "isotropic": {"damping": 6.64, "frequency": 0.900},
    "matched": {"damping": 1.45, "frequency": 1.023},
}
DAMPING_TOLERANCE = 0.3
FREQUENCY_TOLERANCE = 0.01


@pytest.fixture(scope="module", params=sorted(MEASURED))
def droplet_run(request):
    """One shared run per source stencil."""
    stencil = request.param
    run = run_program(
        "oscillation_high_ratio",
        artifacts_dir() / f"oscillation_high_ratio_{stencil}",
        args=(f"--source-stencil={stencil}",),
        timeout=2400,
    )
    run.key = stencil
    return run


@pytest.fixture(scope="module")
def result(droplet_run):
    """The fitted mode, the exact one, and the track they come from."""
    rho1, rho2 = droplet_run.parameter("rho1"), droplet_run.parameter("rho2")
    exact = droplet_mode(
        2,
        droplet_run.parameter("radius"),
        rho1,
        rho2,
        rho1 * droplet_run.parameter("nu"),
        rho2 * droplet_run.parameter("nu2"),
        droplet_run.parameter("sigma"),
    )
    timesteps, signal = droplet_run.mode_track()
    keep = timesteps >= SETTLING_STEPS
    fit = fit_damped_oscillation(signal[keep], dt=timesteps[1] - timesteps[0])
    return {"exact": exact, "fit": fit, "signal": signal, "measured": MEASURED[droplet_run.key]}


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_high_ratio_runs_the_intended_case(droplet_run):
    """The options reached the solver: a flag lost on the way is a different case."""
    assert droplet_run.parameter("rho1") / droplet_run.parameter("rho2") == 1000.0
    assert droplet_run.config["source_stencil"] == droplet_run.key
    assert droplet_run.config["collision"] == "mrt"
    assert droplet_run.config["third_moment_correction"] == "true"


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_high_ratio_stays_bounded(droplet_run, result):
    """No NaN in the track or the last fields, and the phase field within [-1, 1]."""
    assert np.isfinite(result["signal"]).all()
    t = droplet_run.last_timestep
    for field in droplet_run.fields(t).values():
        assert np.isfinite(field).all()
    phase = droplet_run.phase(t)
    assert phase.max() <= 1.0 + 1.0e-6
    assert phase.min() >= -1.0 - 1.0e-6


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_high_ratio_starts_at_the_laid_down_deformation(result):
    """`D = 2 eps` for a small deformation eps, less what the diffuse profile smears."""
    assert result["signal"][0] == pytest.approx(2 * 0.03, rel=0.05)


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_high_ratio_is_one_mode(result):
    """What is fitted is a single sinusoid in an exponential envelope."""
    fit = result["fit"]
    assert fit.residual < 0.03 * fit.amplitude


@pytest.mark.long
@pytest.mark.validation
def test_validation_oscillation_high_ratio_against_the_normal_mode(result):
    """Decay rate and frequency over the exact mode's, at their measured values."""
    exact, fit, measured = result["exact"], result["fit"], result["measured"]
    damping = fit.decay_rate / exact.decay_rate
    frequency = fit.angular_frequency / exact.angular_frequency
    assert damping == pytest.approx(measured["damping"], abs=DAMPING_TOLERANCE)
    assert frequency == pytest.approx(measured["frequency"], abs=FREQUENCY_TOLERANCE)
