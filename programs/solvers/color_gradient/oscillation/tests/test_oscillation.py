"""Mode-2 oscillation of a droplet at a density ratio of 1000, against the
exact viscous normal mode.

A droplet of radius 20, deformed by 3 % into `cos 2 theta` and released, rings
down. The reference is the normal mode of a liquid cylinder in another viscous
fluid, `pycglbm.normal_modes.droplet_mode`: decay rate 1.91e-5 per step and
angular frequency 4.53e-4 for these parameters, where the weak-damping formula
would say 2.15e-5.

Inside the droplet the mode-2 flow is a pure linear strain, which no finite
difference of it gets wrong, so the case isolates the interface in motion --
the part of the colour-gradient model that fails at this density ratio. The
light side of the interface moves ten to twenty times faster than the flow
around it should, and what becomes of the oscillation depends on how that
motion is damped. With the published nine-point source derivative, whose
excess normal stress damps it, the droplet rings down about eight times too
fast. With `--source-stencil=matched` it rings down at 0.62 of the exact rate.
That stencil's six-point face values used to overshoot the jump of the
corrected quantity at the interface, which cancelled the heavy fluid's
viscous normal stress there, and the droplet did not ring down at all (-0.23);
they are now limited where the density changes (docs/numerics.md, "The
matched stencil at an interface"). What is still missing follows the heavy
fluid's bulk relaxation time, 6.5 here: at 1 the droplet is at 0.95, and the
capillary wave is no longer damped ("The heavy fluid's bulk rate, and the
trace of the correction"). The velocity-based solver runs the same case in
`oscillation_vb`.

Both are known gaps, pinned as measured with a tolerance that catches a change
in the scheme. The frequency, which the interface disturbs much less, is
asserted closer.
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
    "isotropic": {"damping": 7.90, "frequency": 0.963},
    "matched": {"damping": 0.62, "frequency": 1.007},
}
DAMPING_TOLERANCE = 0.3
FREQUENCY_TOLERANCE = 0.01


@pytest.fixture(scope="module", params=sorted(MEASURED))
def droplet_run(request):
    """One shared run per source stencil."""
    stencil = request.param
    run = run_program(
        "oscillation",
        artifacts_dir() / f"oscillation_{stencil}",
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
def test_verification_oscillation_runs_the_intended_case(droplet_run):
    """The options reached the solver: a flag lost on the way is a different case."""
    assert droplet_run.parameter("rho1") / droplet_run.parameter("rho2") == 1000.0
    assert droplet_run.config["source_stencil"] == droplet_run.key
    assert droplet_run.config["viscosity_mixing"] == "dynamic"
    assert droplet_run.config["surface_tension"] == "stress"


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_stays_bounded(droplet_run, result):
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
def test_verification_oscillation_starts_at_the_laid_down_deformation(result):
    """`D = 2 eps` for a small deformation eps, less what the diffuse profile smears."""
    assert result["signal"][0] == pytest.approx(2 * 0.03, rel=0.05)


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_is_one_mode(result):
    """What is fitted is a single sinusoid in an exponential envelope."""
    fit = result["fit"]
    assert fit.residual < 0.03 * fit.amplitude


@pytest.mark.long
@pytest.mark.validation
def test_validation_oscillation_against_the_normal_mode(result):
    """Decay rate and frequency over the exact mode's, at their measured values."""
    exact, fit, measured = result["exact"], result["fit"], result["measured"]
    damping = fit.decay_rate / exact.decay_rate
    frequency = fit.angular_frequency / exact.angular_frequency
    assert damping == pytest.approx(measured["damping"], abs=DAMPING_TOLERANCE)
    assert frequency == pytest.approx(measured["frequency"], abs=FREQUENCY_TOLERANCE)
