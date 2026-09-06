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
@pytest.mark.parametrize("boundary", ["periodic", "wall"])
def test_unit_test_lbm_solver_keeps_the_phase_field_bounded(solver_report, stencil, boundary):
    """phi must stay in [-1, 1], against a bounce-back wall as well as away from one.

    The wall case is the one with history: the streaming step used to rebuild g
    from a population other than the one it had just written, and phi reached
    1.29 there while the periodic case stayed within 1e-13 of 1.
    """
    values = solver_report[stencil]
    assert float(values[f"{boundary}_max_abs_phase"]) <= 1.0 + PHASE_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
@pytest.mark.parametrize("boundary", ["rest_periodic", "rest_wall"])
def test_unit_test_lbm_solver_leaves_a_uniform_state_at_rest(solver_report, stencil, boundary):
    """With no interface and no gravity, nothing should drive a flow."""
    values = solver_report[stencil]
    assert float(values[f"{boundary}_max_speed"]) < REST_SPEED_TOLERANCE
    assert float(values[f"{boundary}_mass_drift"]) < MASS_DRIFT_TOLERANCE


#: Density ratios the normalised phase field is reported at.
RATIOS = (1, 20, 1000, 100000)


@pytest.mark.unit_test
@pytest.mark.parametrize("ratio", RATIOS)
def test_unit_test_lbm_solver_normalised_phase_fixes_the_bulks(solver_report, ratio):
    """phi_N must leave both bulks where they are, at any density ratio."""
    values = solver_report["E8"]
    assert float(values[f"phin_r{ratio}_at_plus_one"]) == pytest.approx(1.0, abs=1e-12)
    assert float(values[f"phin_r{ratio}_at_minus_one"]) == pytest.approx(-1.0, abs=1e-12)


@pytest.mark.unit_test
@pytest.mark.parametrize("ratio", RATIOS)
def test_unit_test_lbm_solver_normalised_phase_marks_the_interface(solver_report, ratio):
    """phi_N must cross zero where the components balance, not where phi does.

    The two coincide only at a density ratio of 1. This is the whole point of
    the normalisation: the raw colour field's zero drifts into the light fluid
    as the ratio grows, taking the surface-tension and recolouring operators
    with it.
    """
    values = solver_report["E8"]
    assert float(values[f"phin_r{ratio}_at_interface"]) == pytest.approx(0.0, abs=1e-12)

    offset = abs(float(values[f"phin_r{ratio}_at_zero"]))
    if ratio == 1:
        assert offset == pytest.approx(0.0, abs=1e-12)
    else:
        # (rho1 - rho2) / (rho1 + rho2), i.e. how far phi = 0 sits from the
        # interface once measured in the normalised field.
        assert offset == pytest.approx((ratio - 1) / (ratio + 1), rel=1e-9)


@pytest.mark.unit_test
def test_unit_test_lbm_solver_normalised_phase_is_the_identity_at_equal_densities(solver_report):
    """With rho1 == rho2 the normalisation must change nothing."""
    assert float(solver_report["E8"]["phin_identity_error"]) < 1.0e-14


@pytest.mark.unit_test
@pytest.mark.parametrize("ratio", RATIOS)
def test_unit_test_lbm_solver_normalised_phase_clamps_overshoot(solver_report, ratio):
    """A phase field past 1 must saturate, not run off."""
    assert float(solver_report["E8"][f"phin_r{ratio}_clamped_above"]) == pytest.approx(1.0)


@pytest.mark.unit_test
def test_unit_test_lbm_solver_equilibrium_carries_the_enhanced_third_order_term(solver_report):
    """The equilibrium here must equal Ba et al. Eq. (14), term for term.

    Leclaire et al. (2013) fixed the third-order velocity moment of the
    colour-gradient equilibrium, and Ba et al. (2016) build their
    high-density-ratio model on that correction. The equilibrium in
    ``Solver::equilibrium`` is a Hermite expansion and never names their
    parameter ``alpha``, so the two look unrelated; they are the same
    expression. Pinning it here means a change to the equilibrium that quietly
    drops the correction fails a fast test rather than a validation run.
    """
    difference = float(solver_report["E8"]["enhanced_equilibrium_difference"])
    assert difference < 1.0e-12


#: Density ratios the solver must survive, as decimal exponents.
#: The scheme diverged before step 200 above a ratio of about 100 while the
#: interface was started out of mechanical equilibrium. See docs/numerics.md.
HIGH_RATIOS = (3, 5)

#: Peak spurious velocity tolerated on the short high-ratio runs, in lattice
#: units. Measured: 3.3e-3 at 10^3 and 5.4e-3 at 10^5. The lattice sound speed
#: is 0.577, so these are Mach 1e-2; the transient that used to break these runs
#: reached Mach 1.4.
HIGH_RATIO_MAX_SPEED = 2.0e-2


@pytest.mark.unit_test
@pytest.mark.parametrize("exponent", HIGH_RATIOS)
def test_unit_test_lbm_solver_survives_high_density_ratio(solver_report, exponent):
    """A droplet at a density ratio of 10^n must still be there after 300 steps."""
    values = solver_report["E8"]
    assert int(values[f"ratio_1e{exponent}_finite"]) == 1


@pytest.mark.unit_test
@pytest.mark.parametrize("exponent", HIGH_RATIOS)
def test_unit_test_lbm_solver_high_density_ratio_stays_subsonic(solver_report, exponent):
    """The start must not launch an acoustic transient, which is what used to
    destroy these runs before they could relax."""
    values = solver_report["E8"]
    assert float(values[f"ratio_1e{exponent}_max_speed"]) < HIGH_RATIO_MAX_SPEED


@pytest.mark.unit_test
@pytest.mark.parametrize("exponent", HIGH_RATIOS)
def test_unit_test_lbm_solver_high_density_ratio_keeps_phase_bounded(solver_report, exponent):
    """phi must stay in [-1, 1] at high density ratio too."""
    values = solver_report["E8"]
    assert float(values[f"ratio_1e{exponent}_max_abs_phase"]) <= 1.0 + PHASE_TOLERANCE
