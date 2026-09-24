"""Laplace-law benchmark at a density ratio of 10^4.

The same program as test_laplace_color_gradient.py, run as `laplace E8 1e4 1`:
a droplet 10^4 times denser than its surroundings, with the same dynamic
viscosity in both fluids. Before the high-density-ratio scheme (volume-fraction
initial state, colour gradient on the normalised colour field, capillary stress
as a body force, per-component viscosities -- see docs/numerics.md) this case
diverged within its first twenty steps.

What is checked:

- the initial state is exact: the density interface on the prescribed radius,
  the Laplace jump already in place, and the phase field phi = 0 contour
  displaced by (W/2) ln(rho_1/rho_2) because phi is a mass fraction;
- the run stays bounded and settles;
- the relaxed jump obeys Laplace's law against the radius the density settles at;
- the spurious currents stay at their measured level.

Only static droplets are claimed at this ratio. A droplet translating faster
than about 1e-4 lattice units per step still breaks the scheme at 10^4; that
limit and its cause are recorded in docs/numerics.md.
"""

import math

import numpy as np
import pytest
from pycglbm.testing import artifacts_dir, run_program

# Compile-time constants of programs/solvers/color_gradient/laplace/main_laplace.cpp.
C_DX = 1.0e-5
C_DT = C_DX / 347.0 / math.sqrt(3.0)
SIGMA = 1.0 / (C_DX**3 / C_DT**2)  # surface tension, lattice units
RADIUS = 10.0  # prescribed initial radius
NUM_STEPS = 30000
INIT_WIDTH = 1.1  # ch_width_init

#: Command-line arguments: stencil, density ratio, dynamic viscosity ratio.
DENSITY_RATIO = 1.0e4
ARGS = ("E8", "1e4", "1")

LAPLACE_JUMP = SIGMA / RADIUS
INNER_RADIUS = 0.5 * RADIUS
OUTER_RADIUS = 3.0 * RADIUS

# Measured on this case, pinned to catch regressions.
#: Delta p / (sigma / R_rho) at the end of the run.
MEASURED_JUMP_RATIO = 1.015
MEASURED_JUMP_TOLERANCE = 0.02
#: Largest spurious velocity at the end of the run, lattice units.
MEASURED_MAX_VELOCITY = 7.2e-4


@pytest.fixture(scope="module")
def high_ratio_run():
    """One shared run of the Laplace case at density ratio 1e4."""
    return run_program("laplace", artifacts_dir() / "laplace_ratio_1e4", args=ARGS, timeout=2400)


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_density_ratio_initial_state(high_ratio_run):
    """The jump holds at t = 0 and the interfaces start where the theory puts them."""
    jump = high_ratio_run.pressure_jump(0, inner=INNER_RADIUS, outer=OUTER_RADIUS)
    assert jump == pytest.approx(LAPLACE_JUMP, rel=1.0e-3)

    density = high_ratio_run.density(0)
    assert density.max() == pytest.approx(DENSITY_RATIO, rel=1.0e-5)
    assert high_ratio_run.density_interface_radius(0) == pytest.approx(RADIUS, abs=0.05)

    # phi is a mass fraction: its zero sits (W/2) ln(rho_1/rho_2) further out
    offset = 0.5 * INIT_WIDTH * math.log(DENSITY_RATIO)
    assert high_ratio_run.phase_interface_radius(0) == pytest.approx(RADIUS + offset, abs=0.1)


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_density_ratio_stays_bounded(high_ratio_run):
    """No NaN, the phase field within [-1, 1], positive pressure everywhere."""
    t = high_ratio_run.last_timestep
    assert t == NUM_STEPS
    for field in high_ratio_run.fields(t).values():
        assert np.isfinite(field).all()
    phase = high_ratio_run.phase(t)
    assert phase.max() <= 1.0 + 1.0e-6
    assert phase.min() >= -1.0 - 1.0e-6
    assert high_ratio_run.pressure(t).min() > 0.0


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_density_ratio_settles(high_ratio_run):
    """Over the last two thirds the jump and the radius barely move."""
    tail = [t for t in high_ratio_run.timesteps if t >= NUM_STEPS // 3]
    jumps = [high_ratio_run.pressure_jump(t, inner=INNER_RADIUS, outer=OUTER_RADIUS) for t in tail]
    radii = [high_ratio_run.density_interface_radius(t) for t in tail]

    assert max(jumps) - min(jumps) < 1.0e-2 * LAPLACE_JUMP
    assert max(radii) - min(radii) < 0.06


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_high_density_ratio_pressure_jump(high_ratio_run):
    """Laplace's law against the radius the density field settles at.

    The remaining two percent is the finite interface width at R = 10, not the
    density ratio: the same scheme gives 1.021 at a density ratio of 20 and
    1.023 at 1000, and the force alone integrates to 1.007 sigma/R at R = 20
    (unit test).
    """
    radius = high_ratio_run.density_interface_radius(NUM_STEPS)
    jump = high_ratio_run.pressure_jump(NUM_STEPS, inner=INNER_RADIUS, outer=OUTER_RADIUS)

    ratio = jump / (SIGMA / radius)
    assert ratio == pytest.approx(MEASURED_JUMP_RATIO, abs=MEASURED_JUMP_TOLERANCE)
    assert ratio == pytest.approx(1.0, abs=0.05)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_high_density_ratio_spurious_currents(high_ratio_run):
    """Parasitic currents stay at their measured level, and well below c_s."""
    velocity = high_ratio_run.velocity(NUM_STEPS)
    speed = np.hypot(velocity[..., 0], velocity[..., 1])
    assert speed.max() < 3.0 * MEASURED_MAX_VELOCITY
