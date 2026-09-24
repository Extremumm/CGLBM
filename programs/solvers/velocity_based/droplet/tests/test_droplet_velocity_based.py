"""The velocity-based droplet solver at a density ratio of 1e4.

Two runs of programs/solvers/velocity_based/droplet:

- `droplet E8 1e4 0`: a static droplet, the Laplace benchmark, to compare with
  the colour-gradient `laplace E8 1e4 1`;
- `droplet E8 1e4 0.01`: the same droplet launched at 0.01 lattice units per
  step into a fluid at rest. The colour-gradient solver diverges on this case
  within a few hundred steps, at any speed from 1e-3 up.

The moving run is checked for what the scheme promises -- it stays bounded,
the droplet keeps moving, no velocity appears beyond the flow -- and for what
it does not: total momentum drifts, and the drift is pinned so that it cannot
grow unnoticed. See src/lbm/velocity_based.h and docs/numerics.md.
"""

import math

import numpy as np
import pytest
from pycglbm.testing import artifacts_dir, run_program

# Compile-time constants of main_droplet.cpp.
C_DX = 1.0e-5
C_DT = C_DX / 347.0 / math.sqrt(3.0)
SIGMA = 1.0 / (C_DX**3 / C_DT**2)  # surface tension, lattice units
RADIUS = 10.0
NUM_STEPS = 10000
DENSITY_RATIO = 1.0e4
SPEED = 0.01  # initial speed of the moving droplet

LAPLACE_JUMP = SIGMA / RADIUS
INNER_RADIUS = 0.5 * RADIUS
OUTER_RADIUS = 3.0 * RADIUS

# Measured on these runs, pinned to catch regressions.
MEASURED_JUMP_RATIO = 1.017
MEASURED_MAX_VELOCITY = 5.4e-6
#: Mean speed of the droplet's centre over the run, as a fraction of SPEED.
MEASURED_DROPLET_SPEED = 0.857
#: Relative drift of the total momentum over the run. Not a desired property.
MEASURED_MOMENTUM_DRIFT = 0.044


@pytest.fixture(scope="module")
def static_run():
    return run_program(
        "droplet", artifacts_dir() / "droplet_static", args=("E8", "1e4", "0"), timeout=2400
    )


@pytest.fixture(scope="module")
def moving_run():
    return run_program(
        "droplet", artifacts_dir() / "droplet_moving", args=("E8", "1e4", str(SPEED)), timeout=2400
    )


def _centre_x(case, timestep):
    """x of the droplet's centre, periodic-aware, weighted by the volume fraction."""
    c = 0.5 * (1.0 + case.phase(timestep))
    nx = c.shape[1]
    angle = 2.0 * np.pi * np.arange(nx) / nx
    s = (c * np.sin(angle)[None, :]).sum()
    k = (c * np.cos(angle)[None, :]).sum()
    return (math.atan2(s, k) / (2.0 * np.pi) * nx) % nx


def _momentum(case, timestep):
    return float((case.density(timestep) * case.velocity(timestep)[..., 0]).sum())


@pytest.mark.long
@pytest.mark.verification
def test_verification_droplet_velocity_based_initial_state(static_run):
    jump = static_run.pressure_jump(0, inner=INNER_RADIUS, outer=OUTER_RADIUS)
    assert jump == pytest.approx(LAPLACE_JUMP, rel=1.0e-3)
    assert static_run.density(0).max() == pytest.approx(DENSITY_RATIO, rel=1.0e-4)
    # the phase field is the volume fraction: it crosses zero on the density interface
    assert static_run.phase_interface_radius(0) == pytest.approx(RADIUS, abs=0.05)
    assert static_run.density_interface_radius(0) == pytest.approx(RADIUS, abs=0.05)


@pytest.mark.long
@pytest.mark.verification
def test_verification_droplet_velocity_based_static_is_stationary(static_run):
    tail = [t for t in static_run.timesteps if t >= NUM_STEPS // 3]
    jumps = [static_run.pressure_jump(t, inner=INNER_RADIUS, outer=OUTER_RADIUS) for t in tail]
    radii = [static_run.phase_interface_radius(t) for t in tail]
    assert max(jumps) - min(jumps) < 1.0e-3 * LAPLACE_JUMP
    assert max(radii) - min(radii) < 0.01


@pytest.mark.long
@pytest.mark.validation
def test_validation_droplet_velocity_based_pressure_jump(static_run):
    """Laplace's law; the two percent left is the interface width at R = 10."""
    radius = static_run.phase_interface_radius(NUM_STEPS)
    jump = static_run.pressure_jump(NUM_STEPS, inner=INNER_RADIUS, outer=OUTER_RADIUS)
    ratio = jump / (SIGMA / radius)
    assert ratio == pytest.approx(MEASURED_JUMP_RATIO, abs=0.01)
    assert ratio == pytest.approx(1.0, abs=0.05)


@pytest.mark.long
@pytest.mark.validation
def test_validation_droplet_velocity_based_spurious_currents(static_run):
    velocity = static_run.velocity(NUM_STEPS)
    speed = np.hypot(velocity[..., 0], velocity[..., 1])
    assert speed.max() < 3.0 * MEASURED_MAX_VELOCITY


@pytest.mark.long
@pytest.mark.verification
def test_verification_droplet_velocity_based_moving_stays_bounded(moving_run):
    """The run completes, and the phase field and density stay in range."""
    t = moving_run.last_timestep
    assert t == NUM_STEPS
    for field in moving_run.fields(t).values():
        assert np.isfinite(field).all()
    phase = moving_run.phase(t)
    assert phase.min() >= -1.0 - 1.0e-9
    assert phase.max() <= 1.0 + 1.0e-9
    velocity = moving_run.velocity(t)
    # no velocity beyond the flow the droplet drives
    assert np.hypot(velocity[..., 0], velocity[..., 1]).max() < 1.5 * SPEED


@pytest.mark.long
@pytest.mark.validation
def test_validation_droplet_velocity_based_droplet_keeps_moving(moving_run):
    """The droplet's centre travels at its measured mean speed.

    The interface nodes start at c * SPEED, so the droplet as a whole starts at
    its momentum over its mass, 0.84 SPEED. Dragging the lighter fluid would
    slow it by a fraction of a percent at this density ratio; the momentum drift
    of the scheme instead carries it up to 0.88 SPEED by the end of the run.
    """
    timesteps = moving_run.timesteps
    travelled = 0.0
    previous = _centre_x(moving_run, timesteps[0])
    for t in timesteps[1:]:
        current = _centre_x(moving_run, t)
        step = current - previous
        # unwrap across the periodic boundary
        step -= 128.0 * round(step / 128.0)
        travelled += step
        previous = current
    mean_speed = travelled / (timesteps[-1] - timesteps[0]) / SPEED
    assert mean_speed == pytest.approx(MEASURED_DROPLET_SPEED, abs=0.02)


@pytest.mark.long
@pytest.mark.validation
def test_validation_droplet_velocity_based_momentum_drift(moving_run):
    """Total momentum drifts by the measured amount, and no more.

    The scheme evolves u rather than rho u and is not conservative; this pins
    the drift at a density ratio of 1e4 so that a change which worsens it shows.
    """
    start = _momentum(moving_run, 0)
    end = _momentum(moving_run, NUM_STEPS)
    drift = abs(end - start) / abs(start)
    assert drift < 1.5 * MEASURED_MOMENTUM_DRIFT
