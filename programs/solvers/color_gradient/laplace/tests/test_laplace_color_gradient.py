"""Laplace-law benchmark for the color-gradient solver.

A static droplet of radius ``R`` should carry a pressure jump ``dp = sigma / R``
across its interface. The case initialises the pressure field so that the jump
holds exactly at ``t = 0`` (``p1_inf`` carries a ``- sigma / radius`` offset),
which the verification test checks; the validation test then looks at the state
the solver actually relaxes to.
"""

import math

import pytest

from pycglbm.testing import artifacts_dir, run_program

# Compile-time constants of programs/solvers/color_gradient/laplace/main_laplace.cpp.
# Keep in sync with the constants block at the top of that file.
C_DX = 1.0e-5
C_DT = C_DX / 347.0 / math.sqrt(3.0)
SIGMA = 1.0 / (C_DX**3 / C_DT**2)  # surface tension, lattice units
RADIUS = 10.0
NUM_STEPS = 30000

#: Analytic Laplace jump for the prescribed radius.
LAPLACE_JUMP = SIGMA / RADIUS

#: Radii (lattice units, from the domain centre) delimiting bulk and far field.
INNER_RADIUS = 0.5 * RADIUS
OUTER_RADIUS = 3.0 * RADIUS

# The converged solver reaches only a fraction of the analytic jump. This is a
# measured baseline that guards against regressions, NOT a statement that the
# case reproduces the Laplace law -- see docs/numerics.md.
MEASURED_JUMP_RATIO = 0.715
MEASURED_JUMP_TOLERANCE = 0.02


@pytest.fixture(scope="module")
def laplace_run():
    """One shared run of the Laplace case, reused by every test in this module."""
    return run_program("laplace", artifacts_dir() / "laplace", timeout=1800)


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_color_gradient_initial_pressure_jump(laplace_run):
    """At t = 0 the prescribed field must satisfy dp = sigma / R exactly."""
    jump = laplace_run.pressure_jump(0, inner=INNER_RADIUS, outer=OUTER_RADIUS)
    assert jump == pytest.approx(LAPLACE_JUMP, rel=1.0e-3)


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_color_gradient_droplet_is_preserved(laplace_run):
    """The droplet must neither dissolve nor drift in size while relaxing."""
    initial = laplace_run.droplet_radius(0)
    final = laplace_run.droplet_radius(NUM_STEPS)
    assert initial == pytest.approx(RADIUS, abs=0.5)
    assert final == pytest.approx(initial, abs=0.5)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_steady_pressure_jump(laplace_run):
    """The relaxed pressure jump is stationary, at a known fraction of sigma / R."""
    timesteps = [t for t in laplace_run.timesteps if t >= NUM_STEPS // 3]
    jumps = [
        laplace_run.pressure_jump(t, inner=INNER_RADIUS, outer=OUTER_RADIUS) for t in timesteps
    ]

    # stationary: the tail of the run must not drift
    assert max(jumps) - min(jumps) < 1.0e-4 * LAPLACE_JUMP

    ratio = jumps[-1] / LAPLACE_JUMP
    assert ratio == pytest.approx(MEASURED_JUMP_RATIO, abs=MEASURED_JUMP_TOLERANCE)
