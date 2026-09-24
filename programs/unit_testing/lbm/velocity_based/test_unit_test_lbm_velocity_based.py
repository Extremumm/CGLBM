"""Unit tests for src/lbm/velocity_based.

The scheme stands on a few exact properties, checked here without running a
solver: the moments of its equilibria and forcing, the conservation of P and u
by the collision, the positivity of the phase-field carrier, a pressure
correction that leaves no force from a uniform pressure across a density jump
of 1e4, and a phase transport that keeps its profile, its mass and its bounds.
"""

import pytest
from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "lbm_velocity_based"


@pytest.fixture(scope="module")
def values():
    return parse_key_values(run_unit_program(PROGRAM).stdout)


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_moments(values):
    """g_eq has moments P, u, P cs^2 I + u u; the forcing has 0, a, u a + a u."""
    assert float(values["equilibrium_moment_error"]) < 1e-15
    assert float(values["forcing_moment_error"]) < 1e-15
    assert float(values["phase_moment_error"]) < 1e-15


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_collision_conserves(values):
    """The collision keeps P, and u up to the half-step force."""
    assert float(values["collision_conservation_error"]) < 1e-15
    assert float(values["unit_tau_error"]) < 1e-15


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_carrier_is_positive(values):
    """Gamma_i(u) >= 0 up to |u| = 0.1, which is what bounds the phase field."""
    for speed in (1, 5, 10):
        assert float(values[f"carrier_min_{speed}"]) > 0.0


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_uniform_pressure_exerts_no_force(values):
    """Across a density jump of 1e4, a uniform p gives no force, to rounding.

    P = p / (rho cs^2) itself jumps by 1e4 there; the correction is written on
    the lattice's own stencil so that the two cancel exactly rather than to the
    accuracy of a finite difference.
    """
    assert float(values["uniform_pressure_residual"]) < 1e-12
    assert float(values["uniform_density_correction"]) == 0.0


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_phase_transport(values):
    """The phase field keeps its tanh profile, its total and its bounds.

    The one percent of deviation is the difference between the discrete
    steady profile and the continuum tanh it is initialised with.
    """
    assert float(values["phase_mass_error"]) < 1e-12
    assert float(values["phase_profile_deviation"]) < 0.02
    assert float(values["phase_min"]) >= 0.0
    assert float(values["phase_max"]) <= 1.0
