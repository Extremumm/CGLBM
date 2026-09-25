"""Unit tests for the two-population colour-gradient solver.

``src/lbm/solver.h`` and ``src/lbm/two_population_solver.h`` are both
colour-gradient models, and they differ in where the density ratio lives. This
one carries a distribution per fluid and puts the ratio in the equilibrium's
rest weight, which changes which invariants matter: each fluid's mass is
conserved separately rather than only the total, and both distributions must
stay non-negative, because that -- not a clamp -- is what bounds the phase
field.

The equilibrium moments are checked directly. They are the model: the second
moment has to come out as ``rho_k (c_s^k)^2 + rho_k u u`` with each fluid's own
sound speed, since that is how the density ratio reaches the pressure.
"""

import pytest
from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "two_population"

STENCILS = ("E4", "E6", "E8")

#: Density ratios the solver is exercised at, as they appear in the report.
RATIOS = ("r20", "r1000", "r100000")

#: Relative error allowed in the equilibrium's moments. Measured: 2e-16.
MOMENT_TOLERANCE = 1.0e-12

#: Relative mass drift allowed per fluid. Measured: below 1e-11 at every ratio.
MASS_DRIFT_TOLERANCE = 1.0e-9

#: How far outside [-1, 1] the phase field may stray. Measured: exactly 1.
PHASE_TOLERANCE = 1.0e-12

#: A uniform fluid must stay exactly at rest. Measured: 0.
REST_SPEED_TOLERANCE = 1.0e-12


@pytest.fixture(scope="module")
def report():
    """One run of the unit program per stencil, shared by the tests below."""
    return {s: parse_key_values(run_unit_program(PROGRAM, (s,)).stdout) for s in STENCILS}


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_equilibrium_has_the_right_moments(report, stencil):
    """Mass, momentum and momentum flux, at density ratios from 1 to 1e5.

    The momentum flux is the one that matters: it must be
    ``rho_k (c_s^k)^2 delta + rho_k u u``, each fluid with its own sound speed.
    That is the whole mechanism by which this model carries a density ratio, and
    an equilibrium that got it wrong would still conserve mass and momentum.
    """
    values = report[stencil]
    assert float(values["equilibrium_mass_error"]) < MOMENT_TOLERANCE
    assert float(values["equilibrium_momentum_error"]) < MOMENT_TOLERANCE
    assert float(values["equilibrium_stress_error"]) < MOMENT_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_rest_weights_encode_the_density_ratio(report, stencil):
    """(1 - alpha_2) / (1 - alpha_1) must be rho1 / rho2, Ba et al. Eq. (7)."""
    assert float(report[stencil]["alpha_density_ratio_error"]) < 1.0e-9


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
@pytest.mark.parametrize("ratio", RATIOS)
def test_unit_test_two_population_conserves_each_fluid_separately(report, stencil, ratio):
    """Neither fluid may leak into the other.

    The recolouring redistributes colour between the two distributions every
    step, so this is a stronger statement than total mass conservation and it is
    the one that can actually fail.
    """
    values = report[stencil]
    assert float(values[f"{ratio}_mass1_drift"]) < MASS_DRIFT_TOLERANCE
    assert float(values[f"{ratio}_mass2_drift"]) < MASS_DRIFT_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
@pytest.mark.parametrize("ratio", RATIOS)
def test_unit_test_two_population_keeps_both_distributions_non_negative(report, stencil, ratio):
    """Neither distribution may go negative, at any density ratio.

    This is the invariant the recolouring had to be adapted for. Latva-Kokko and
    Rothman push along the lattice weight ``w_i``; with the fluids on different
    rest weights, what has to stay positive is a fraction of the *population*,
    and in the heavy fluid the non-rest populations are ``(1 - alpha_1)/5`` of
    the density -- 1.6e-4 at a ratio of 1000 against ``w_i = 1/9``. Pushing by
    ``w_i`` there drives the light fluid's distribution negative by a factor of
    240 and the run is gone within ten steps. See ``rest_weight``.
    """
    assert float(report[stencil][f"{ratio}_min_population"]) >= 0.0


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
@pytest.mark.parametrize("ratio", RATIOS)
def test_unit_test_two_population_keeps_the_phase_field_bounded(report, stencil, ratio):
    """phi_N stays in [-1, 1] without being clamped there.

    It is ``(s1 - s2) / (s1 + s2)`` of two non-negative densities, so this
    follows from the test above rather than from any explicit bound -- which is
    why it holds to the last bit here and does not in the other solver.
    """
    assert float(report[stencil][f"{ratio}_max_abs_phase"]) <= 1.0 + PHASE_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_leaves_a_uniform_fluid_at_rest(report, stencil):
    """No interface and no gravity: nothing may start moving."""
    assert float(report[stencil]["rest_max_speed"]) < REST_SPEED_TOLERANCE
