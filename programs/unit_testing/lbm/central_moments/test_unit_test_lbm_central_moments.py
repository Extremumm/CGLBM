"""Unit tests for src/lbm/central_moments.

The central-moment collision of Saito et al. (2023) stands on three facts
checked here: the transforms to and from central moments are exact on D2Q9, the
generalized equilibrium built from its central moments is the Hermite form
derived for D2Q9 (at a heavy fluid's pressure too, where p - rho cs^2 is about
-rho cs^2), and the collision conserves mass and momentum, relaxes the stress
at the rates it is given and sets the higher moments to equilibrium.
"""

import pytest
from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "lbm_central_moments"


@pytest.fixture(scope="module")
def values():
    return parse_key_values(run_unit_program(PROGRAM).stdout)


@pytest.mark.unit_test
def test_unit_test_lbm_central_moments_round_trip(values):
    """Populations to central moments and back, to rounding."""
    assert float(values["round_trip_error"]) < 1e-13


@pytest.mark.unit_test
def test_unit_test_lbm_central_moments_generalized_equilibrium(values):
    """The equilibrium built from rho, p, p and p cs^2 is the Hermite form, and
    its central moments do not depend on the velocity."""
    assert float(values["hermite_equilibrium_error"]) < 1e-14
    assert float(values["equilibrium_moment_error"]) < 1e-14


@pytest.mark.unit_test
def test_unit_test_lbm_central_moments_solver_equilibrium(values):
    """Solver's equilibrium agrees up to the second order and differs in k_21,
    k_12 and k_22 by exactly the velocity terms the generalized one removes."""
    assert float(values["solver_equilibrium_low_order_error"]) < 1e-14
    assert float(values["solver_equilibrium_k22_error"]) < 1e-14


@pytest.mark.unit_test
def test_unit_test_lbm_central_moments_collision(values):
    """Mass and momentum kept, shear and trace relaxed at their rates, the third
    and fourth order at equilibrium; at unit rates, the equilibrium itself."""
    assert float(values["collision_conservation_error"]) < 1e-14
    assert float(values["collision_relaxation_error"]) < 1e-14
    assert float(values["collision_unit_rate_error"]) < 1e-14


@pytest.mark.unit_test
def test_unit_test_lbm_central_moments_trace_at_unit_rate_is_stable(values):
    """A uniform fluid at tau = 100 with a heavy fluid's non-ideal pressure:
    with the trace at rate 1, as Solver runs it, noise of 1e-6 stays small for
    400 steps; with the trace at the slow rate 1/tau it grows past 1."""
    assert float(values["uniform_speed_unit_bulk"]) < 1e-5
    assert float(values["uniform_speed_slow_bulk"]) >= 1.0
