"""Unit tests for the shared colour-gradient kernel, src/lbm/solver.h.

The solvers under ``programs/solvers`` run one physical case each for tens of
thousands of steps; they are the validation suite and they are slow. These tests
drive the same kernel on a 32x32 lattice for twenty steps and check the
invariants it must hold at every step, whatever the case: mass is conserved, the
phase field stays inside the range the equation of state is defined on, and a
state with no interface and no gravity does not start moving on its own.
"""

import pytest
from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "lbm_solver"

STENCILS = ("E4", "E6", "E8")

#: Relative mass drift allowed over the run. Measured drift is around 4e-15,
#: i.e. accumulated rounding; anything above this is a leak in the scheme.
MASS_DRIFT_TOLERANCE = 1.0e-12

#: How far outside [-1, 1] the phase field may stray. The recolouring step
#: overshoots by a rounding error, which the equation of state clamps; a real
#: excursion means populations are being mixed up.
PHASE_TOLERANCE = 1.0e-9

#: A state with no interface must stay exactly at rest. Measured peak is 3e-16.
REST_SPEED_TOLERANCE = 1.0e-12


@pytest.fixture(scope="module")
def solver_report():
    """One run of the unit program per stencil, shared by the tests below."""
    return {s: parse_key_values(run_unit_program(PROGRAM, (s,)).stdout) for s in STENCILS}


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_lbm_solver_reports_its_stencil(solver_report, stencil):
    assert solver_report[stencil]["stencil"] == stencil


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
@pytest.mark.parametrize("boundary", ["periodic", "wall"])
def test_unit_test_lbm_solver_conserves_mass(solver_report, stencil, boundary):
    """Streaming moves populations; none of the three operators creates them."""
    values = solver_report[stencil]
    assert float(values[f"{boundary}_mass_drift"]) < MASS_DRIFT_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_lbm_solver_keeps_the_phase_field_bounded(solver_report, stencil):
    """phi must stay in [-1, 1] on a periodic domain."""
    values = solver_report[stencil]
    assert float(values["periodic_max_abs_phase"]) <= 1.0 + PHASE_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.xfail(
    reason=(
        "Known defect in the bounce-back streaming step: g is rebuilt from "
        "f[ip][jp][k] while f was written to f[ip][jp][kp], so on a reflected "
        "direction it reads an unrelated population and phi leaves [-1, 1]. "
        "Measured peak is 1.29. See docs/numerics.md."
    ),
    strict=True,
)
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_lbm_solver_keeps_the_phase_field_bounded_at_a_wall(solver_report, stencil):
    """phi must stay in [-1, 1] against a bounce-back wall too."""
    values = solver_report[stencil]
    assert float(values["wall_max_abs_phase"]) <= 1.0 + PHASE_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
@pytest.mark.parametrize("boundary", ["rest_periodic", "rest_wall"])
def test_unit_test_lbm_solver_leaves_a_uniform_state_at_rest(solver_report, stencil, boundary):
    """With no interface and no gravity, nothing should drive a flow."""
    values = solver_report[stencil]
    assert float(values[f"{boundary}_max_speed"]) < REST_SPEED_TOLERANCE
    assert float(values[f"{boundary}_mass_drift"]) < MASS_DRIFT_TOLERANCE
