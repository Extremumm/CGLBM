"""Laplace-law benchmark for the colour-gradient solver.

A static droplet of radius R should carry a pressure jump dp = sigma / R across
its interface. The case initialises the pressure so the jump holds exactly at
t = 0 (p1_inf carries a -sigma/radius offset), which the verification test
checks; the validation tests then look at the state the solver relaxes to.

Which radius to score against is not a formality here. The droplet is
compressible, so it does not stay at the prescribed radius, and since the
two-component equation of state was restored the phase field and the density
field settle at *different* radii -- see docs/numerics.md. The tests below
therefore measure the radius rather than assuming it, and say which one they
use.
"""

import pytest
from pycglbm.testing import artifacts_dir, run_program

# The case parameters are not repeated here: the solver reports them as
# `key = value` lines at the top of run.log, and CaseOutput.parameter reads them
# back. A constant duplicated in a test is a constant that drifts.

# Measured baselines for the current scheme (restored equation of state, E8
# colour gradient, interface started in mechanical equilibrium). They are what
# the solver does, pinned to catch regressions. The history behind them, on this
# same case at density ratio 20:
#   linear-mixing pressure, E4 gradient      : 0.720 * sigma/R_rho   (28 % low)
#   restored EOS, E4 gradient                : 0.962 * sigma/R_rho
#   restored EOS, E8 gradient                : 0.962 * sigma/R_rho, fewer currents
#   the above, started in mechanical equilib.: 0.895 * sigma/R_rho
# The last step traded 7 % on this case for a scheme that runs at density ratios
# up to 10^5 instead of diverging above 100 -- see docs/numerics.md. Run with
# `--initial-state=linear` to reproduce the row above it.
MEASURED_JUMP_RATIO = 0.895
MEASURED_JUMP_TOLERANCE = 0.03

#: The phase and density interfaces settle this far apart, in lattice units.
#: An artefact of the restored equation of state, tracked so that a scheme
#: change which removes it shows up here rather than passing unnoticed. It is
#: barely moved by the initialisation (2.270 linear, 2.284 equilibrium), which
#: is part of the evidence that it belongs to the equation of state.
MEASURED_INTERFACE_SPLIT = 2.28
INTERFACE_SPLIT_TOLERANCE = 0.4

#: Spurious currents at steady state, in lattice units.
MEASURED_MAX_VELOCITY = 1.19e-3

#: Radius the density field starts at, in lattice units.
#: The phase field starts exactly on the prescribed radius; the density does
#: not, because mechanical equilibrium -- not a linear interpolation -- is what
#: sets the density profile across the interface.
MEASURED_INITIAL_DENSITY_RADIUS = 8.41


@pytest.fixture(scope="module")
def laplace_run():
    """One shared run of the Laplace case, reused by every test in this module."""
    return run_program("laplace", artifacts_dir() / "laplace", timeout=1800)


@pytest.fixture(scope="module")
def case(laplace_run):
    """The parameters the run reported, and the quantities derived from them."""
    sigma = laplace_run.parameter("sigma")
    radius = laplace_run.parameter("radius")
    return {
        "sigma": sigma,
        "radius": radius,
        "steps": laplace_run.parameter("steps", int),
        # Analytic Laplace jump for the prescribed radius.
        "jump": sigma / radius,
        # Radii, from the domain centre, delimiting bulk and far field.
        "inner": 0.5 * radius,
        "outer": 3.0 * radius,
    }


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_color_gradient_initial_pressure_jump(laplace_run, case):
    """At t = 0 the prescribed field must satisfy dp = sigma / R exactly."""
    jump = laplace_run.pressure_jump(0, inner=case["inner"], outer=case["outer"])
    assert jump == pytest.approx(case["jump"], rel=1.0e-3)
    # The phase field starts exactly on the prescribed radius.
    assert laplace_run.phase_interface_radius(0) == pytest.approx(case["radius"], abs=0.1)
    # The density does not: its profile comes from solving the equation of
    # state for mechanical equilibrium, which is not symmetric about phi = 0.
    assert laplace_run.density_interface_radius(0) == pytest.approx(
        MEASURED_INITIAL_DENSITY_RADIUS, abs=0.1
    )


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_color_gradient_reaches_a_steady_state(laplace_run, case):
    """The tail of the run must be stationary, not still relaxing."""
    tail = [t for t in laplace_run.timesteps if t >= case["steps"] // 3]
    jumps = [laplace_run.pressure_jump(t, inner=case["inner"], outer=case["outer"]) for t in tail]
    radii = [laplace_run.density_interface_radius(t) for t in tail]

    assert max(jumps) - min(jumps) < 1.0e-4 * case["jump"]
    assert max(radii) - min(radii) < 0.01


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_pressure_jump(laplace_run, case):
    """The relaxed jump, scored against the radius the density field settles at."""
    radius = laplace_run.density_interface_radius(case["steps"])
    jump = laplace_run.pressure_jump(case["steps"], inner=case["inner"], outer=case["outer"])

    ratio = jump / (case["sigma"] / radius)
    assert ratio == pytest.approx(MEASURED_JUMP_RATIO, abs=MEASURED_JUMP_TOLERANCE)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_interfaces_stay_split(laplace_run, case):
    """The phase and density interfaces sit a known distance apart.

    They coincided under the old linear-mixing pressure and separate under the
    restored equation of state. This is an open issue, not a desired property:
    the assertion pins the size of the gap so that any scheme change which
    closes it -- or widens it -- is visible immediately.
    """
    phase_radius = laplace_run.phase_interface_radius(case["steps"])
    density_radius = laplace_run.density_interface_radius(case["steps"])

    split = phase_radius - density_radius
    assert split == pytest.approx(MEASURED_INTERFACE_SPLIT, abs=INTERFACE_SPLIT_TOLERANCE)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_spurious_currents(laplace_run, case):
    """Parasitic currents at the interface must stay at their measured level."""
    velocity = laplace_run.velocity(case["steps"])
    speed = (velocity[..., 0] ** 2 + velocity[..., 1] ** 2) ** 0.5
    assert speed.max() == pytest.approx(MEASURED_MAX_VELOCITY, rel=0.15)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_higher_isotropy_helps(laplace_run, case):
    """E8 must not have more parasitic current than the E4 stencil.

    Leclaire et al. (Computers & Fluids 48, 98, 2011) report the isotropic
    gradient cutting spurious currents; this checks the direction of that claim
    on the real case rather than on the stencil alone.
    """
    coarse = run_program("laplace", artifacts_dir() / "laplace_e4", args=("E4",), timeout=1800)

    def peak(run):
        v = run.velocity(case["steps"])
        return float(((v[..., 0] ** 2 + v[..., 1] ** 2) ** 0.5).max())

    assert peak(laplace_run) <= peak(coarse)
