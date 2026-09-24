"""Laplace-law benchmark for the colour-gradient solver.

A static droplet of radius R should carry a pressure jump dp = sigma / R across
its interface. The case initialises the pressure so the jump holds exactly at
t = 0 (p1_inf carries a -sigma/radius offset), which the verification test
checks; the validation tests then look at the state the solver relaxes to.

Which radius to score against is not a formality here. The droplet is
compressible, so it does not stay at the prescribed radius, and the phase field
and the density field sit at *different* radii: phi is a mass fraction, and its
zero contour lies (W/2) ln(rho_1/rho_2) outside the density interface -- see
docs/numerics.md and src/lbm/mixture.h. The surface tension acts on the
density interface, so the jump is scored against the density radius.

This is the shipped case, `laplace` with no arguments: density ratio 20, one
kinematic viscosity for both fluids. test_laplace_high_density_ratio.py runs
the same program at a density ratio of 1e4.
"""

import math

import pytest
from pycglbm.testing import artifacts_dir, run_program

# Compile-time constants of programs/solvers/color_gradient/laplace/main_laplace.cpp.
# Keep in sync with the constants block at the top of that file.
C_DX = 1.0e-5
C_DT = C_DX / 347.0 / math.sqrt(3.0)
SIGMA = 1.0 / (C_DX**3 / C_DT**2)  # surface tension, lattice units
RADIUS = 10.0  # prescribed initial radius
NUM_STEPS = 30000
DENSITY_RATIO = 20.0
INIT_WIDTH = 1.1  # ch_width_init, the width the initial profile is given
WIDTH = 1.6  # ch_width_ope, the width the recolouring maintains

#: Analytic Laplace jump for the initial radius.
LAPLACE_JUMP = SIGMA / RADIUS

#: Radii (lattice units, from the domain centre) delimiting bulk and far field.
INNER_RADIUS = 0.5 * RADIUS
OUTER_RADIUS = 3.0 * RADIUS

# Measured baselines for the current scheme, pinned to catch regressions. The
# history behind them, on this same case at density ratio 20:
#   linear-mixing pressure, E4 gradient      : 0.720 * sigma/R_rho, max|u| 1.44e-3
#   restored EOS, E4 gradient                : 0.962 * sigma/R_rho, max|u| 1.22e-3
#   restored EOS, E8 gradient                : 0.962 * sigma/R_rho, max|u| 1.11e-3
#   high-density-ratio scheme, E4            : 1.017 * sigma/R_rho, max|u| 1.77e-5
#   high-density-ratio scheme, E8 (current)  : 1.021 * sigma/R_rho, max|u| 1.35e-5
# The two percent left over is the finite interface width at R = 10; it falls
# to 0.4 % at R = 20.
MEASURED_JUMP_RATIO = 1.021
MEASURED_JUMP_TOLERANCE = 0.02

#: Offset of the phi = 0 contour from the density interface, (W/2) ln(rho_1/rho_2).
#: This used to be pinned as an unexplained artefact (2.27); it is the mass-
#: fraction nature of phi, and the solver now settles within 2 % of it.
INTERFACE_SPLIT_TOLERANCE = 0.1

#: Spurious currents at steady state, in lattice units.
MEASURED_MAX_VELOCITY = 1.35e-5


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
    # the density interface starts on the prescribed radius, and the phase
    # interface (W/2) ln(rho_1/rho_2) outside it
    assert laplace_run.density_interface_radius(0) == pytest.approx(RADIUS, abs=0.1)
    offset = 0.5 * INIT_WIDTH * math.log(DENSITY_RATIO)
    assert laplace_run.phase_interface_radius(0) == pytest.approx(RADIUS + offset, abs=0.1)


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_color_gradient_reaches_a_steady_state(laplace_run):
    """The tail of the run must be stationary, not still relaxing."""
    tail = [t for t in laplace_run.timesteps if t >= NUM_STEPS // 3]
    jumps = [laplace_run.pressure_jump(t, inner=INNER_RADIUS, outer=OUTER_RADIUS) for t in tail]
    radii = [laplace_run.density_interface_radius(t) for t in tail]

    assert max(jumps) - min(jumps) < 1.0e-4 * LAPLACE_JUMP
    assert max(radii) - min(radii) < 0.01


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_pressure_jump(laplace_run):
    """The relaxed jump, scored against the radius the density field settles at."""
    radius = laplace_run.density_interface_radius(NUM_STEPS)
    jump = laplace_run.pressure_jump(NUM_STEPS, inner=INNER_RADIUS, outer=OUTER_RADIUS)

    ratio = jump / (SIGMA / radius)
    assert ratio == pytest.approx(MEASURED_JUMP_RATIO, abs=MEASURED_JUMP_TOLERANCE)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_interface_split_is_the_mass_fraction_offset(
    laplace_run,
):
    """The phase and density interfaces sit (W/2) ln(rho_1/rho_2) apart.

    phi is a mass fraction and the density follows the volume fraction; a tanh
    profile of width W in one is the same profile in the other, shifted by that
    much (src/lbm/mixture.h). The recolouring maintains the width W, so the
    relaxed split is predicted, not merely pinned.
    """
    phase_radius = laplace_run.phase_interface_radius(NUM_STEPS)
    density_radius = laplace_run.density_interface_radius(NUM_STEPS)

    split = phase_radius - density_radius
    expected = 0.5 * WIDTH * math.log(DENSITY_RATIO)
    assert split == pytest.approx(expected, abs=INTERFACE_SPLIT_TOLERANCE)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_spurious_currents(laplace_run):
    """Parasitic currents at the interface must stay at their measured level."""
    velocity = laplace_run.velocity(NUM_STEPS)
    speed = (velocity[..., 0] ** 2 + velocity[..., 1] ** 2) ** 0.5
    assert speed.max() == pytest.approx(MEASURED_MAX_VELOCITY, rel=0.25)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_higher_isotropy_helps(laplace_run):
    """E8 must not have more parasitic current than the E4 stencil.

    Leclaire et al. (Computers & Fluids 48, 98, 2011) report the isotropic
    gradient cutting spurious currents; this checks the direction of that claim
    on the real case rather than on the stencil alone.
    """
    coarse = run_program("laplace", artifacts_dir() / "laplace_e4", args=("E4",), timeout=1800)

    def peak(case):
        v = case.velocity(NUM_STEPS)
        return float(((v[..., 0] ** 2 + v[..., 1] ** 2) ** 0.5).max())

    assert peak(laplace_run) <= peak(coarse)
