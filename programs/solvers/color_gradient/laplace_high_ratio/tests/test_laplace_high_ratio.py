"""Laplace-law benchmark at a density ratio of 1000.

The same physical test as ``laplace``, run with the two-population solver. The
point of having both is that this one can be here at all: ``laplace`` runs at a
density ratio of 20, and its scheme -- which puts the whole density ratio into a
single field through an equation of state -- diverges at 1000.

The parameters follow Ba et al. (2016) Table I, so the numbers below can be read
against theirs directly. They report a 0.74 % error on the interfacial tension
at this density ratio, with spurious currents of 1.25e-4.
"""

import pytest
from pycglbm.testing import artifacts_dir, run_program

#: Measured relative error on the interfacial tension, sigma_cal / sigma_th.
#: Reached from below and flat to three digits over the last 10^4 steps.
#: Ba et al. report 1.0074 for the same case.
MEASURED_JUMP_RATIO = 1.003
MEASURED_JUMP_TOLERANCE = 0.01

#: Spurious currents at the end of the run, in lattice units. They fall
#: monotonically -- 1.8e-3, 8.4e-4, 3.7e-4, 1.7e-4, ... -- so this is an upper
#: bound on a decaying quantity rather than a stationary value.
#: Ba et al. report 1.25e-4 for the same case.
MEASURED_MAX_VELOCITY = 4.34e-5

#: The droplet's radius, from the phi_N = 0 contour. It starts at exactly the
#: prescribed 25 and relaxes slightly inward as the interface takes its profile.
MEASURED_RADIUS = 24.895


@pytest.fixture(scope="module")
def high_ratio_run():
    """One shared run of the case, reused by every test in this module."""
    return run_program("laplace_high_ratio", artifacts_dir() / "laplace_high_ratio", timeout=1800)


@pytest.fixture(scope="module")
def case(high_ratio_run):
    sigma = high_ratio_run.parameter("sigma")
    radius = high_ratio_run.parameter("radius")
    return {
        "sigma": sigma,
        "radius": radius,
        "steps": high_ratio_run.parameter("steps", int),
        "ratio": high_ratio_run.parameter("rho1") / high_ratio_run.parameter("rho2"),
        "inner": 0.5 * radius,
        "outer": 1.5 * radius,
    }


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_ratio_starts_at_the_prescribed_radius(high_ratio_run, case):
    """At t = 0 the droplet is exactly the size the case asked for.

    This model has no equation of state to invert: the initial profile is the
    volume fraction, so ``rho_1 = c rho1`` and ``rho_2 = (1 - c) rho2`` place
    the interface where the case put it, at any density ratio.
    """
    assert case["ratio"] == pytest.approx(1000.0)
    assert high_ratio_run.phase_interface_radius(0) == pytest.approx(case["radius"], abs=0.05)


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_ratio_starts_without_a_pressure_jump(high_ratio_run, case):
    """The two bulk pressures are equal at t = 0, so the jump has to grow.

    ``rho_k (c_s^k)^2`` is the same for both fluids by construction -- that is
    what ``(1 - alpha_2) / (1 - alpha_1) = rho1 / rho2`` buys -- so unlike
    ``laplace`` this case does not start with the answer imposed on it. What it
    relaxes to is entirely the scheme's doing.
    """
    jump = high_ratio_run.pressure_jump(0, inner=case["inner"], outer=case["outer"])
    assert abs(jump) < 1.0e-9


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_ratio_reaches_a_steady_state(high_ratio_run, case):
    """The jump must have settled, not still be climbing.

    Worth checking rather than assuming: at a density ratio of 1000 a run of
    3 x 10^4 steps covers 0.6 viscous relaxation times in the *other* solver,
    and the high-ratio claims this repository used to make were transients that
    had not yet diverged. Matching the dynamic viscosities is what makes tau
    uniform here and the run settle inside its length.
    """
    tail = [t for t in high_ratio_run.timesteps if t >= 3 * case["steps"] // 4]
    jumps = [
        high_ratio_run.pressure_jump(t, inner=case["inner"], outer=case["outer"]) for t in tail
    ]
    assert max(jumps) - min(jumps) < 1.0e-3 * max(jumps)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_high_ratio_pressure_jump(high_ratio_run, case):
    """Laplace's law at a density ratio of 1000, to under a per cent."""
    radius = high_ratio_run.phase_interface_radius(case["steps"])
    jump = high_ratio_run.pressure_jump(case["steps"], inner=case["inner"], outer=case["outer"])

    assert radius == pytest.approx(MEASURED_RADIUS, abs=0.1)
    assert jump * radius / case["sigma"] == pytest.approx(
        MEASURED_JUMP_RATIO, abs=MEASURED_JUMP_TOLERANCE
    )


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_high_ratio_spurious_currents(high_ratio_run, case):
    """Parasitic currents must stay at their measured level."""
    velocity = high_ratio_run.velocity(case["steps"])
    speed = (velocity[..., 0] ** 2 + velocity[..., 1] ** 2) ** 0.5
    assert speed.max() == pytest.approx(MEASURED_MAX_VELOCITY, rel=0.2)
