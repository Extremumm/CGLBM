"""Laplace's law in three dimensions, at a density ratio of 1000.

The three-dimensional counterpart of ``laplace_high_ratio``, with the same
model and parameters. What is being tested is different in one respect that
matters: for a sphere the curvature is 2/R, so the law reads

    dp = 2 sigma / R

and an operator carried over from two dimensions without thought would give half
of it. That failure mode is checked more directly by
``programs/unit_testing/lbm/two_population_3d``, which measures the discrete
curvature of an analytic sphere against -2/R; this case is the end-to-end
statement.

The run writes the z = nz/2 slice in the same four CSV files a two-dimensional
run produces, so ``pycglbm`` reads it unchanged. The velocities reported below
are therefore the in-plane ones of that slice, not the full three-dimensional
maximum -- which is the right comparison for a spherically symmetric state and
is stated here so the number is not mistaken for something else.
"""

import pytest
from pycglbm.testing import artifacts_dir, run_program

#: Measured Δp R / (2 sigma) at the end of the run. It is approached from
#: above and still creeping down: 1.0265 here at 1.5e4 steps, 1.0261 at 2e4 and
#: 1.0240 at 3e4. The residual is resolution, not dimension -- R = 10 in two
#: dimensions gives +2.1 % against +0.3 % at R = 25.
MEASURED_JUMP_RATIO = 1.0266
MEASURED_JUMP_TOLERANCE = 0.01

#: In-plane spurious currents on the mid-plane slice, in lattice units.
MEASURED_MAX_VELOCITY = 1.9e-5

#: The droplet's radius from the phi_N = 0 contour of the slice. It starts at
#: exactly the prescribed 10 and relaxes slightly inward.
MEASURED_RADIUS = 9.768


@pytest.fixture(scope="module")
def run_3d():
    """One shared run of the case. It takes about twelve minutes on 24 cores."""
    return run_program("laplace_3d", artifacts_dir() / "laplace_3d", timeout=2700)


@pytest.fixture(scope="module")
def case(run_3d):
    radius = run_3d.parameter("radius")
    return {
        "sigma": run_3d.parameter("sigma"),
        "radius": radius,
        "steps": run_3d.parameter("steps", int),
        "ratio": run_3d.parameter("rho1") / run_3d.parameter("rho2"),
        "inner": 0.5 * radius,
        "outer": 1.8 * radius,
    }


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_3d_is_a_cube_at_the_prescribed_radius(run_3d, case):
    """The lattice is genuinely three-dimensional and the droplet is the right size."""
    assert run_3d.parameter("nz", int) == run_3d.parameter("nx", int)
    assert run_3d.parameter("nz", int) > 1
    assert case["ratio"] == pytest.approx(1000.0)
    assert run_3d.phase_interface_radius(0) == pytest.approx(case["radius"], abs=0.05)


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_3d_starts_without_a_pressure_jump(run_3d, case):
    """The two bulk pressures are equal at t = 0, so the jump has to grow.

    As in two dimensions, this model does not have the answer imposed on it by a
    pressure offset: what it relaxes to is entirely the scheme's doing.
    """
    jump = run_3d.pressure_jump(0, inner=case["inner"], outer=case["outer"])
    assert abs(jump) < 1.0e-9


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_3d_is_close_to_steady(run_3d, case):
    """The jump must have all but settled by the end of the run.

    It approaches its limit from above rather than oscillating, so this asserts
    that the last quarter of the run moves it by well under a per cent -- not
    that it has stopped. Running to 3e4 steps takes it from 1.0265 to 1.0240.
    """
    tail = [t for t in run_3d.timesteps if t >= 3 * case["steps"] // 4]
    jumps = [run_3d.pressure_jump(t, inner=case["inner"], outer=case["outer"]) for t in tail]
    assert max(jumps) - min(jumps) < 5.0e-3 * max(jumps)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_3d_pressure_jump(run_3d, case):
    """Laplace's law for a sphere: dp = 2 sigma / R, to under three per cent."""
    radius = run_3d.phase_interface_radius(case["steps"])
    jump = run_3d.pressure_jump(case["steps"], inner=case["inner"], outer=case["outer"])

    assert radius == pytest.approx(MEASURED_RADIUS, abs=0.1)
    assert jump * radius / (2.0 * case["sigma"]) == pytest.approx(
        MEASURED_JUMP_RATIO, abs=MEASURED_JUMP_TOLERANCE
    )


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_3d_spurious_currents(run_3d, case):
    """Parasitic currents on the mid-plane must stay at their measured level."""
    velocity = run_3d.velocity(case["steps"])
    speed = (velocity[..., 0] ** 2 + velocity[..., 1] ** 2) ** 0.5
    assert speed.max() == pytest.approx(MEASURED_MAX_VELOCITY, rel=0.25)
