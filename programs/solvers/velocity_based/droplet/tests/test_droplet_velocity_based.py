"""The velocity-based droplet solver at large density ratios.

Three runs of programs/solvers/velocity_based/droplet:

- `droplet E8 1e4 0`: a static droplet, the Laplace benchmark, to compare with
  the colour-gradient `laplace E8 1e4 1`;
- `droplet E8 1e4 0.01`: the same droplet launched at 0.01 lattice units per
  step into a fluid at rest. The colour-gradient solver diverges on this case
  within a few hundred steps, at any speed from 1e-3 up;
- `droplet E8 100 0.01`: the same at a density ratio of 100, where the lighter
  fluid carries a large share of the momentum.

The moving runs are checked for what the scheme promises: they stay bounded,
total momentum is conserved to the precision of the output, and the droplet
slows down as it sets the lighter fluid in motion. See src/lbm/velocity_based.h
and docs/numerics.md.
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
LX = 128
DENSITY_RATIO = 1.0e4
SPEED = 0.01  # initial speed of the moving droplets

LAPLACE_JUMP = SIGMA / RADIUS
INNER_RADIUS = 0.5 * RADIUS
OUTER_RADIUS = 3.0 * RADIUS

# Measured on these runs, pinned to catch regressions.
MEASURED_JUMP_RATIO = 1.030
MEASURED_MAX_VELOCITY = 1.7e-6
#: Mean speed of the droplet's centre over the run, as a fraction of SPEED.
MEASURED_DROPLET_SPEED = {"1e4": 0.844, "100": 0.715}
#: The output is written with ten significant digits; momentum is conserved
#: to rounding, and its sum over the lattice to about that precision.
MOMENTUM_TOLERANCE = 1.0e-9


@pytest.fixture(scope="module")
def static_run():
    return run_program(
        "droplet", artifacts_dir() / "droplet_static", args=("E8", "1e4", "0"), timeout=2400
    )


@pytest.fixture(scope="module", params=["1e4", "100"])
def moving_run(request):
    return request.param, run_program(
        "droplet",
        artifacts_dir() / f"droplet_moving_{request.param}",
        args=("E8", request.param, str(SPEED)),
        timeout=2400,
    )


def _centre_x(case, timestep):
    """x of the droplet's centre, periodic-aware, weighted by the volume fraction."""
    c = 0.5 * (1.0 + case.phase(timestep))
    nx = c.shape[1]
    angle = 2.0 * np.pi * np.arange(nx) / nx
    s = (c * np.sin(angle)[None, :]).sum()
    k = (c * np.cos(angle)[None, :]).sum()
    return (math.atan2(s, k) / (2.0 * np.pi) * nx) % nx


def _speeds(case):
    """Speed of the droplet's centre over each output interval, over SPEED."""
    timesteps = case.timesteps
    speeds = []
    previous = _centre_x(case, timesteps[0])
    for before, after in zip(timesteps[:-1], timesteps[1:]):
        current = _centre_x(case, after)
        step = current - previous
        # unwrap across the periodic boundary
        step -= LX * round(step / LX)
        speeds.append(step / (after - before) / SPEED)
        previous = current
    return speeds


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
    """Laplace's law; the three percent left is the interface width at R = 10."""
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
    _, case = moving_run
    t = case.last_timestep
    assert t == NUM_STEPS
    for field in case.fields(t).values():
        assert np.isfinite(field).all()
    phase = case.phase(t)
    assert phase.min() >= -1.0 - 1.0e-9
    assert phase.max() <= 1.0 + 1.0e-9
    velocity = case.velocity(t)
    # no velocity beyond the flow the droplet drives
    assert np.hypot(velocity[..., 0], velocity[..., 1]).max() < 1.1 * SPEED


@pytest.mark.long
@pytest.mark.verification
def test_verification_droplet_velocity_based_momentum_is_conserved(moving_run):
    """Total momentum is conserved at every output, to the output's precision.

    Every link exchanges equal and opposite momentum, and the pressure and
    capillary forces sum to zero over the lattice. The previous version of the
    scheme, which read the streamed velocity as it stands, drifted by 4.4 % over
    this run at a density ratio of 1e4 and by 11.7 % at 100.
    """
    _, case = moving_run
    start = _momentum(case, 0)
    for t in case.timesteps[1:]:
        assert _momentum(case, t) == pytest.approx(start, rel=MOMENTUM_TOLERANCE)


@pytest.mark.long
@pytest.mark.validation
def test_validation_droplet_velocity_based_droplet_slows_down(moving_run):
    """The droplet travels at its measured mean speed, and slows down.

    The interface nodes start at c * SPEED, so the droplet as a whole starts at
    its momentum over its mass, 0.84 SPEED at a density ratio of 1e4. It can
    only lose speed to the lighter fluid it drags along, which at 1e4 holds a
    fraction of a percent of the momentum and at 100 about a third.
    """
    ratio, case = moving_run
    speeds = _speeds(case)
    assert np.mean(speeds) == pytest.approx(MEASURED_DROPLET_SPEED[ratio], abs=0.01)
    # never faster than over the interval before, and slower at the end
    assert all(later < earlier + 1.0e-3 for earlier, later in zip(speeds, speeds[1:]))
    assert speeds[-1] < speeds[0]
