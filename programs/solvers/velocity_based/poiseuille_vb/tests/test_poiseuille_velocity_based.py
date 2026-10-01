"""Layered Poiseuille flow between resting walls, against its exact profile.

Two fluids side by side in a channel, driven along it by a uniform force, the
heavy one below: the standard wall-bounded test of the phase-field models at
large density ratios (Zu & He 2013, Fakhari et al. 2017, Liang et al. 2018,
Subhedar 2022). The steady profile is a parabola in each fluid, meeting with
the same velocity and the same shear stress, and does not depend on the
densities; see programs/solvers/velocity_based/poiseuille_vb.

What it scores is how the viscous stress crosses a diffuse interface. A shear
stress that crosses it is the same in every layer, so the layers' shear rates
add up in series and the viscosity that carries it is the harmonic mean of
theirs; mixing the viscosity arithmetically on the volume fraction, as the
solver does by default, makes the light side of a 1000-to-1 interface forty
times too viscous at c = 0.04, and the light fluid then flows as if the
interface were a wall five nodes into it: 41 % off. The laminate mixing
(`--viscosity=laminate`, SolverParameters::interface_viscosity) takes the
harmonic mean for the shear across the interface and the arithmetic one for
the stretching along it, which is what the moving-interface benchmarks need.

Four runs:

- one fluid, the walls alone: half-way bounce-back at tau = 1, to 1e-4;
- 1000 apart in density and in viscosity, the same kinematic viscosity, from
  rest, mixed arithmetically and as a laminate;
- the case of Liang et al. (2018), a density ratio of 1000 and a viscosity
  ratio of 100 on a channel of 100 nodes, an interface of their width (W = 2.5
  here is their W = 5 in tanh(2 y / W)), scored with their relative L1 error.
  They report 3.2e-2 for their model, and 0.11 and 0.39, under the same
  conditions, for the models of Ren et al. (2016) and Fakhari et al. (2017).

The measured errors are pinned, so a change in the scheme shows up.
"""

import numpy as np
import pytest
from pycglbm.testing import artifacts_dir, default_threads, run_program

#: Density ratio, viscosity ratio, steps, and options of each run.
CASES = {
    "walls": {"density": "1", "viscosity": "1", "steps": 30000, "options": ()},
    "1000_arithmetic": {
        "density": "1000",
        "viscosity": "1000",
        "steps": 30000,
        "options": ("--viscosity=arithmetic",),
    },
    "1000_laminate": {
        "density": "1000",
        "viscosity": "1000",
        "steps": 30000,
        "options": ("--viscosity=laminate",),
    },
    "liang": {
        "density": "1000",
        "viscosity": "100",
        "steps": 300000,
        "options": ("--ny=100", "--width=2.5", "--viscosity=laminate"),
    },
}

#: Measured relative L2 and L1 errors against the exact profile.
MEASURED = {
    "walls": {"l2": 1.062e-4, "l1": 1.163e-4},
    "1000_arithmetic": {"l2": 0.4105, "l1": 0.3890},
    "1000_laminate": {"l2": 0.0457, "l1": 0.0451},
    "liang": {"l2": 0.0209, "l1": 0.0165},
}
MEASURED_TOLERANCE = 0.05

#: Liang et al.'s own model on their case, relative L1 error.
LIANG_MODEL_L1 = 3.2e-2


@pytest.fixture(scope="module", params=sorted(CASES))
def channel_run(request):
    """One shared run per case."""
    case = CASES[request.param]
    run = run_program(
        "poiseuille_vb",
        artifacts_dir() / f"poiseuille_vb_{request.param}",
        args=("E8", case["density"], case["viscosity"], str(case["steps"]), *case["options"]),
        timeout=3600,
        threads=default_threads(),
    )
    run.key = request.param
    return run


def exact_profile(run, y):
    """The layered Poiseuille profile at height y' from the centre line,
    computed here from the run's parameters, not read from its output."""
    mu1, mu2 = run.parameter("mu1"), run.parameter("mu2")
    drive = run.parameter("drive")
    h = run.parameter("ny") / 2.0
    mu = np.where(y < 0.0, mu1, mu2)
    s = y / h
    return (
        drive
        * h
        * h
        / (2.0 * mu)
        * (-s * s - s * (mu2 - mu1) / (mu1 + mu2) + 2.0 * mu / (mu1 + mu2))
    )


@pytest.fixture(scope="module")
def errors(channel_run):
    """The run's relative L2 and L1 errors, against the exact profile."""
    profile = np.genfromtxt(channel_run.rundir / "profile.csv", delimiter=",", names=True)
    exact = exact_profile(channel_run, profile["y"])
    difference = profile["ux"] - exact
    return {
        "l2": float(np.sqrt((difference**2).sum() / (exact**2).sum())),
        "l1": float(np.abs(difference).sum() / np.abs(exact).sum()),
        "profile": profile,
        "exact": exact,
        "key": channel_run.key,
    }


@pytest.mark.long
@pytest.mark.verification
def test_verification_poiseuille_vb_runs_the_intended_case(channel_run):
    """The arguments reached the solver, and the drive puts the fastest point
    of the exact profile at 0.01."""
    case = CASES[channel_run.key]
    assert channel_run.parameter("rho1") / channel_run.parameter("rho2") == float(case["density"])
    assert channel_run.parameter("mu1") / channel_run.parameter("mu2") == pytest.approx(
        float(case["viscosity"])
    )
    assert channel_run.parameter("steps", int) == case["steps"]
    y = np.linspace(-0.5, 0.5, 20001) * channel_run.parameter("ny")
    assert exact_profile(channel_run, y).max() == pytest.approx(0.01, rel=1e-6)


@pytest.mark.long
@pytest.mark.verification
def test_verification_poiseuille_vb_is_a_channel_flow(channel_run, errors):
    """The flow ends uniform along the channel, with nothing across it, the
    interface on the centre line, and the program's exact profile the test's."""
    t = channel_run.last_timestep
    velocity = channel_run.velocity(t)
    assert np.isfinite(velocity).all()
    assert np.abs(velocity[..., 0] - velocity[:, :1, 0]).max() < 1e-12
    # the artificial compressibility's, over 3e5 steps on the longest run
    assert np.abs(velocity[..., 1]).max() < 1e-5 * 0.01
    profile = errors["profile"]
    np.testing.assert_allclose(profile["exact"], errors["exact"], rtol=1e-9, atol=1e-15)
    # the volume fraction stays symmetric about the centre line
    assert profile["c"] + profile["c"][::-1] == pytest.approx(1.0, abs=1e-6)


@pytest.mark.long
@pytest.mark.validation
def test_validation_poiseuille_vb_against_the_exact_profile(errors):
    """The errors are where they were measured; the laminate runs beat the
    model Liang et al. report on their own case."""
    measured = MEASURED[errors["key"]]
    assert errors["l2"] == pytest.approx(measured["l2"], rel=MEASURED_TOLERANCE)
    assert errors["l1"] == pytest.approx(measured["l1"], rel=MEASURED_TOLERANCE)
    if errors["key"] == "walls":
        assert errors["l2"] < 2e-4
    if errors["key"] in ("1000_laminate", "liang"):
        assert errors["l1"] < LIANG_MODEL_L1 * 1.5
    if errors["key"] == "liang":
        assert errors["l1"] < LIANG_MODEL_L1
