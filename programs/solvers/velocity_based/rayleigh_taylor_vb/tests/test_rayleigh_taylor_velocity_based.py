"""The Rayleigh-Taylor instability at large density ratios, between walls,
against the exact growth rate of the viscous normal mode.

A heavy fluid above a light one, their interface displaced by 0.05 nodes in a
cosine of wavelength 64 and released under gravity, in a box four
wavelengths tall closed by resting walls; see
programs/solvers/velocity_based/rayleigh_taylor_vb. The interface's first
Fourier coefficient is fitted with an exponential over the stretch where it is
linear, k a from 0.025 to 0.1, and scored against
`pycglbm.normal_modes.capillary_wave` with gravity: the exact normal mode of
two semi-infinite viscous fluids with surface tension, at an Atwood number of
0.998 and 0.9998, whose growth rate the tension and both viscosities hold 18 %
and 2 % below the inviscid `(A g k)^(1/2)`. The walls are 128 nodes from the
interface, where `tanh(k H)` is 1 to ten digits.

The colour-gradient solvers cannot run this case: their `rayleigh_taylor` is at
a density ratio of 4/1. Beyond the linear stretch the runs continue into the
bubble and spike, which the test holds to finite fields and a bounded phase.

A third run doubles the resolution. The measured rates are pinned, so a
change in the scheme shows up.
"""

import numpy as np
import pytest
from pycglbm.normal_modes import capillary_wave
from pycglbm.testing import artifacts_dir, default_threads, run_program

#: Density ratio, heavy viscosity, steps, amplitude and wavelength of each run.
CASES = {
    "1000": {"ratio": "1000", "mu1": 2.0, "steps": 12000, "amplitude": 0.05, "wavelength": 64},
    "1e4": {"ratio": "1e4", "mu1": 2.0, "steps": 12000, "amplitude": 0.05, "wavelength": 64},
    "1000_128": {
        "ratio": "1000",
        "mu1": 2.0,
        "steps": 9000,
        "amplitude": 0.1,
        "wavelength": 128,
    },
}

#: Measured growth rate over the exact mode's.
MEASURED = {"1000": 0.977, "1e4": 0.98, "1000_128": 0.989}
MEASURED_TOLERANCE = 0.005

#: How close to the exact mode the solver has to be, whatever it measured.
GROWTH_BOUND = 0.03

#: The linear stretch the fit takes, as k a.
LINEAR_WINDOW = (0.025, 0.1)


@pytest.fixture(scope="module", params=sorted(CASES))
def rt_run(request):
    """One shared run per case."""
    case = CASES[request.param]
    run = run_program(
        "rayleigh_taylor_vb",
        artifacts_dir() / f"rayleigh_taylor_vb_{request.param}",
        args=(
            "E8",
            case["ratio"],
            str(case["mu1"]),
            str(case["steps"]),
            str(case["amplitude"]),
            f"--wavelength={case['wavelength']}",
        ),
        timeout=3600,
        threads=default_threads(),
    )
    run.key = request.param
    return run


def exact_growth(run):
    """The growth rate of the exact viscous normal mode; the heavy fluid, on
    top, is fluid 2 of capillary_wave."""
    k = 2.0 * np.pi / run.parameter("nx")
    mode = capillary_wave(
        k,
        run.parameter("rho2"),
        run.parameter("rho1"),
        run.parameter("mu2"),
        run.parameter("mu1"),
        run.parameter("sigma"),
        gravity=run.parameter("gravity"),
    )
    assert mode.angular_frequency == 0.0
    return -mode.decay_rate


@pytest.fixture(scope="module")
def growth(rt_run):
    """The fitted growth rate over the exact one, and the track."""
    track = np.genfromtxt(rt_run.rundir / "mode.csv", delimiter=",", names=True)
    k = 2.0 * np.pi / rt_run.parameter("nx")
    ka = k * track["amplitude"]
    linear = (ka > LINEAR_WINDOW[0]) & (ka < LINEAR_WINDOW[1])
    rate = np.polyfit(track["timestep"][linear], np.log(track["amplitude"][linear]), 1)[0]
    fit = np.polyval(
        np.polyfit(track["timestep"][linear], np.log(track["amplitude"][linear]), 1),
        track["timestep"][linear],
    )
    return {
        "ratio": rate / exact_growth(rt_run),
        "residual": float(np.abs(fit - np.log(track["amplitude"][linear])).max()),
        "points": int(linear.sum()),
        "track": track,
        "key": rt_run.key,
    }


@pytest.mark.long
@pytest.mark.verification
def test_verification_rayleigh_taylor_vb_runs_the_intended_case(rt_run):
    """The arguments reached the solver, the heavy fluid is on top and the
    tension lets the wavelength grow."""
    case = CASES[rt_run.key]
    assert rt_run.parameter("rho1") / rt_run.parameter("rho2") == float(case["ratio"])
    assert rt_run.parameter("mu1") == case["mu1"]
    assert rt_run.parameter("nx", int) == case["wavelength"]
    assert rt_run.parameter("ny", int) == 4 * case["wavelength"]
    phase = rt_run.phase(0)
    assert phase[-1].min() > 0.99 and phase[0].max() < -0.99
    assert exact_growth(rt_run) > 0.0


@pytest.mark.long
@pytest.mark.verification
def test_verification_rayleigh_taylor_vb_stays_bounded(rt_run, growth):
    """No NaN in the track or the last fields, the phase field in its bounds,
    and the instability past its linear stretch by the end."""
    track = growth["track"]
    assert np.isfinite(track["amplitude"]).all()
    for field in rt_run.fields(rt_run.last_timestep).values():
        assert np.isfinite(field).all()
    phase = rt_run.phase(rt_run.last_timestep)
    assert phase.min() >= -1.0 - 1.0e-9
    assert phase.max() <= 1.0 + 1.0e-9
    k = 2.0 * np.pi / rt_run.parameter("nx")
    assert k * track["amplitude"][-1] > 3.0 * LINEAR_WINDOW[1]
    # the spike outruns the bubble, as it does at a large Atwood number
    assert -track["spike"][-1] > track["bubble"][-1] > 0.0


@pytest.mark.long
@pytest.mark.verification
def test_verification_rayleigh_taylor_vb_grows_exponentially(growth):
    """What is fitted is one exponential: enough points, a small residual."""
    assert growth["points"] >= 20
    assert growth["residual"] < 0.02


@pytest.mark.long
@pytest.mark.validation
def test_validation_rayleigh_taylor_vb_against_the_normal_mode(growth):
    """The growth rate close to the exact mode's, and at its measured value."""
    assert growth["ratio"] == pytest.approx(1.0, abs=GROWTH_BOUND)
    assert growth["ratio"] == pytest.approx(MEASURED[growth["key"]], abs=MEASURED_TOLERANCE)
