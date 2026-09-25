"""The velocity-based droplet solver at large density ratios.

Six runs of programs/solvers/velocity_based/droplet:

- `droplet E8 1e4 0`: a static droplet, the Laplace benchmark, to compare with
  the colour-gradient `laplace` at a density ratio of 1e4
  (test_laplace_high_density_ratio.py);
- `droplet E8 1e4 0.01`: the same droplet launched at 0.01 lattice units per
  step into a fluid at rest. The colour-gradient solver diverges on this case
  within a few hundred steps, at any speed from 1e-3 up;
- `droplet E8 100 0.01`: the same at a density ratio of 100, where the lighter
  fluid carries a large share of the momentum;
- `droplet E8 1e4 0.1 <1|10|100>`: launched ten times faster, at viscosity
  ratios of 1, 10 and 100. Before the dissipation went through the forcing
  term these diverged within 1500 steps, and 0.05 after 7700.

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
FAST_SPEED = 0.1  # initial speed of the fast ones

LAPLACE_JUMP = SIGMA / RADIUS
INNER_RADIUS = 0.5 * RADIUS
OUTER_RADIUS = 3.0 * RADIUS

# Measured on these runs, pinned to catch regressions.
MEASURED_JUMP_RATIO = 0.998
MEASURED_MAX_VELOCITY = 2.3e-6
#: Mean speed of the droplet's centre over the run, as a fraction of SPEED.
MEASURED_DROPLET_SPEED = {"1e4": 0.844, "100": 0.714}
#: The same at FAST_SPEED, by viscosity ratio.
MEASURED_FAST_SPEED = {"1": 0.850, "10": 0.850, "100": 0.851}
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


@pytest.fixture(scope="module", params=["1", "10", "100"])
def fast_run(request):
    return request.param, run_program(
        "droplet",
        artifacts_dir() / f"droplet_fast_{request.param}",
        args=("E8", "1e4", str(FAST_SPEED), request.param),
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


def _speeds(case, speed=SPEED):
    """Speed of the droplet's centre over each output interval, over `speed`.

    Unwrapped across the periodic boundary, which a droplet at 0.1 crosses in
    less than an output interval: the displacement is taken nearest to the
    launch speed times the interval.
    """
    timesteps = case.timesteps
    speeds = []
    previous = _centre_x(case, timesteps[0])
    for before, after in zip(timesteps[:-1], timesteps[1:]):
        current = _centre_x(case, after)
        step = current - previous
        expected = 0.8 * speed * (after - before)
        step += LX * round((expected - step) / LX)
        speeds.append(step / (after - before) / speed)
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
    """Laplace's law, to a few parts in 1e3 at R = 10.

    Without layer_weight the capillary stress supports sigma times the mean of
    1/r across the interface, 3 % above sigma / R at this radius.
    """
    radius = static_run.phase_interface_radius(NUM_STEPS)
    jump = static_run.pressure_jump(NUM_STEPS, inner=INNER_RADIUS, outer=OUTER_RADIUS)
    ratio = jump / (SIGMA / radius)
    assert ratio == pytest.approx(MEASURED_JUMP_RATIO, abs=0.002)
    assert ratio == pytest.approx(1.0, abs=0.005)


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


@pytest.mark.long
@pytest.mark.verification
def test_verification_droplet_velocity_based_fast_stays_bounded(fast_run):
    """At 0.1 lattice units per step the run completes, bounded, and conserves
    momentum, whatever the viscosity ratio."""
    _, case = fast_run
    t = case.last_timestep
    assert t == NUM_STEPS
    for field in case.fields(t).values():
        assert np.isfinite(field).all()
    phase = case.phase(t)
    assert phase.min() >= -1.0 - 1.0e-9
    assert phase.max() <= 1.0 + 1.0e-9
    # the droplet is still one droplet: its core is component 1
    assert phase.max() > 0.99
    velocity = case.velocity(t)
    assert np.hypot(velocity[..., 0], velocity[..., 1]).max() < 1.1 * FAST_SPEED
    start = _momentum(case, 0)
    for step in case.timesteps[1:]:
        assert _momentum(case, step) == pytest.approx(start, rel=MOMENTUM_TOLERANCE)


@pytest.mark.long
@pytest.mark.validation
def test_validation_droplet_velocity_based_fast_droplet_slows_down(fast_run):
    """The fast droplet travels at its measured mean speed, and slows down."""
    ratio, case = fast_run
    speeds = _speeds(case, FAST_SPEED)
    assert np.mean(speeds) == pytest.approx(MEASURED_FAST_SPEED[ratio], abs=0.01)
    assert speeds[-1] < speeds[0]
