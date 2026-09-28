"""Ba et al.'s static droplet at a density ratio of 1000, as they ran it.

``test_laplace_high_ratio.py`` runs the same droplet with the dynamic
viscosities matched, which keeps tau uniform at 0.85. Ba et al. (2016) matched
the *kinematic* viscosities instead, nu = 0.1667 in both fluids, which puts tau
at 348 in the heavy fluid. With BGK that case reads 0.54 sigma/R after 5e3
steps and falls for the rest of the run, to 0.38 at 4e4. This runs it with the
two things they used to carry it:

    laplace_high_ratio --nu=0.1667 --nu-b=0.1667 --nu2=0.1667 --nu-b2=0.1667
                       --collision=mrt --third-moment-correction

the MRT collision of Lallemand & Luo, which relaxes the moments the
Navier-Stokes equations do not contain at their own rates, and the source term
for the diagonal third moment D2Q9 cannot carry. They report 1.0074 on the
tension and spurious currents of 1.25e-4. See docs/numerics.md, "MRT and the
third-moment source in the two-population solver".
"""

import subprocess

import pytest
from pycglbm.testing import artifacts_dir, find_program, run_program

#: Ba et al.'s kinematic viscosity, the same in both fluids.
NU = 0.1667
STEPS = 120000
ARGS = (
    f"--nu={NU}",
    f"--nu-b={NU}",
    f"--nu2={NU}",
    f"--nu-b2={NU}",
    "--collision=mrt",
    "--third-moment-correction",
    f"--steps={STEPS}",
    "--interval=10000",
)

#: Measured sigma_cal / sigma at the end of the run. It is still creeping up,
#: by about 5e-4 per 1e4 steps: 0.9952 at 6e4 steps, 0.9969 at 8e4, 0.9980 at
#: 1e5, 0.9989 at 1.2e5.
MEASURED_JUMP_RATIO = 0.9989
#: Measured spurious currents at the end of the run, falling monotonically:
#: 2.1e-5 at 6e4 steps, 9.7e-6 at 1e5.
MEASURED_MAX_VELOCITY = 7.86e-6

#: What Ba et al. report for this case.
BA_JUMP_RATIO = 1.0074
BA_MAX_VELOCITY = 1.25e-4


@pytest.fixture(scope="module")
def mrt_run():
    """One shared run of Ba et al.'s case, reused by every test in this module."""
    return run_program(
        "laplace_high_ratio", artifacts_dir() / "laplace_high_ratio_mrt", args=ARGS, timeout=3600
    )


@pytest.fixture(scope="module")
def case(mrt_run):
    radius = mrt_run.parameter("radius")
    return {
        "sigma": mrt_run.parameter("sigma"),
        "steps": mrt_run.parameter("steps", int),
        "inner": 0.5 * radius,
        "outer": 1.5 * radius,
    }


def tension_ratio(run, case, timestep):
    """sigma_cal / sigma from the jump and the phi_N = 0 radius."""
    jump = run.pressure_jump(timestep, inner=case["inner"], outer=case["outer"])
    return jump * run.phase_interface_radius(timestep) / case["sigma"]


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_ratio_mrt_runs_the_intended_case(mrt_run, case):
    """The options reached the solver: this is Ba et al.'s case, not the shipped one."""
    assert mrt_run.config["collision"] == "mrt"
    assert mrt_run.config["third_moment_correction"] == "true"
    assert mrt_run.parameter("nu") == pytest.approx(NU)
    assert mrt_run.parameter("nu2") == pytest.approx(NU)
    assert mrt_run.parameter("rho1") / mrt_run.parameter("rho2") == pytest.approx(1000.0)
    assert case["steps"] == STEPS


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_ratio_mrt_approaches_steadily(mrt_run, case):
    """Over the last quarter the tension moves by under 0.3 %, and only upward.

    Not flat to three digits like the matched-viscosity case: tau = 348 in the
    droplet is its stress relaxation time, and the approach is slow.
    """
    tail = [t for t in mrt_run.timesteps if t >= 3 * case["steps"] // 4]
    ratios = [tension_ratio(mrt_run, case, t) for t in tail]
    assert all(b >= a for a, b in zip(ratios, ratios[1:]))
    assert ratios[-1] - ratios[0] < 3.0e-3


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_high_ratio_mrt_pressure_jump(mrt_run, case):
    """Laplace's law to within Ba et al.'s own 0.74 %."""
    ratio = tension_ratio(mrt_run, case, case["steps"])
    assert ratio == pytest.approx(MEASURED_JUMP_RATIO, abs=2.0e-3)
    assert abs(ratio - 1.0) < abs(BA_JUMP_RATIO - 1.0)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_high_ratio_mrt_spurious_currents(mrt_run, case):
    """Currents at their measured level, an order of magnitude under Ba et al.'s."""
    velocity = mrt_run.velocity(case["steps"])
    speed = (velocity[..., 0] ** 2 + velocity[..., 1] ** 2) ** 0.5
    assert speed.max() == pytest.approx(MEASURED_MAX_VELOCITY, rel=0.2)
    assert speed.max() < BA_MAX_VELOCITY


@pytest.mark.long
@pytest.mark.verification
@pytest.mark.parametrize(
    "argument", ["--collision=trt", "--s-e=2", "--s-eps=0", "--s-q=-1", "--alpha2=1"]
)
def test_verification_laplace_high_ratio_mrt_rejects_a_bad_option(argument, tmp_path):
    """A rate outside (0, 2) or an unknown collision is a usage error, not a run.

    A relaxation rate of 2 or more makes a moment grow instead of decay, and
    one of 0 freezes it; neither is a setting anyone means, so the program says
    so and exits with status 2 before writing anything.
    """
    result = subprocess.run(
        [str(find_program("laplace_high_ratio")), argument, "--steps=1"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 2
    assert result.stderr
