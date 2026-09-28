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

#: Density ratios the solver is exercised at, as they appear in the report. The
#: ``mrt_`` runs use the MRT collision and the third-moment source, with Ba et
#: al.'s equal kinematic viscosities: tau reaches 348 in the droplet at 1000.
RATIOS = ("r20", "r1000", "r100000", "mrt_r1000", "mrt_r100000")

#: Absolute error allowed in the MRT kernel and the third-moment source, whose
#: inputs are of order 0.1 and 1e-4. Measured: 5e-16 and below.
KERNEL_TOLERANCE = 1.0e-14

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


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_mrt_is_bgk_at_one_rate(report, stencil):
    """With every rate equal, the MRT collision is BGK with Guo's forcing.

    That is what makes the MRT an extension rather than a different model: the
    moment basis only changes *which* moments relax at which rate.
    """
    assert float(report[stencil]["mrt_bgk_difference"]) < KERNEL_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_mrt_relaxes_each_moment_at_its_own_rate(report, stencil):
    """Each moment of the result is ``m - s (m - m_eq) + (1 - s/2) m_S``.

    And the basis is orthogonal, so its transpose over the row norms inverts it.
    """
    values = report[stencil]
    assert float(values["mrt_moment_error"]) < KERNEL_TOLERANCE
    assert float(values["mrt_orthogonality_error"]) == 0.0


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_third_moment_source_is_ba_eqs_17_18(report, stencil):
    """The source carries no mass, no momentum and no shear.

    It adds ``(1 - s_e/2) div Q`` to the trace of the second moment and
    ``(1 - s_nu/2) (d_x Q_x - d_y Q_y)`` to its normal difference -- the two
    places D2Q9's diagonal third moment reaches the viscous stress.
    """
    values = report[stencil]
    assert float(values["third_moment_mass"]) < KERNEL_TOLERANCE
    assert float(values["third_moment_momentum"]) < KERNEL_TOLERANCE
    assert float(values["third_moment_shear"]) < KERNEL_TOLERANCE
    assert float(values["third_moment_trace_error"]) < 1.0e-12
    assert float(values["third_moment_normal_error"]) < 1.0e-12


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_mrt_conserves_mass_and_momentum(report, stencil):
    """Node by node, on a droplet at 1000 with tau = 348 in the heavy fluid.

    The collision must leave the mass alone and give the momentum exactly the
    force, ``2 rho u - sum f e`` with Guo's half-force velocity, whatever the
    rates and with the third-moment source on.
    """
    values = report[stencil]
    assert float(values["mrt_solver_mass_error"]) < 1.0e-13
    assert float(values["mrt_solver_momentum_error"]) < 1.0e-13


@pytest.mark.unit_test
def test_unit_test_two_population_third_moment_source_repairs_the_normal_stress(report):
    """A Taylor-Green vortex in one fluid decays at ``2 nu k^2`` only with the source.

    The vortex's strain is purely normal, so its decay measures the normal
    viscous stress alone. D2Q9 gives it ``(1 - c_k^2) / (2 c_k^2)`` times its
    value for a fluid of sound speed ``c_k``: 0.5417 at ``c_k^2 = 0.48`` (a
    density ratio of 1 with alpha = 0.2) and 9.917 at 0.048 (a ratio of 10),
    against 0.5425 and 9.919 measured. The source brings the first to 0.995.

    With the nine-point isotropic stencil it cannot cancel the error exactly,
    because it is a finite difference of a quantity ``1 / c_k^2`` times larger
    than the stress it leaves behind, and the streaming that made the error is
    not that stencil. The difference is of order ``k^2 / c_k^2``: 1.139 at a
    ratio of 10 on this 32^2 lattice, 1.035 on 64^2, 17.5 at a ratio of 1000.
    The stencil that matches the streaming (``--source-stencil=matched``)
    removes it: 0.997 at a ratio of 10 and 0.998 at 1000. Pinned as measured,
    so a change to either operator shows up here.
    """
    values = report["E8"]
    assert float(values["tg_r1_uncorrected"]) == pytest.approx(0.5425, abs=2e-3)
    assert float(values["tg_r10_uncorrected"]) == pytest.approx(9.919, abs=0.02)
    assert float(values["tg_r1_corrected"]) == pytest.approx(0.9950, abs=2e-3)
    assert float(values["tg_r10_corrected"]) == pytest.approx(1.1393, abs=2e-3)
    assert float(values["tg_matched_r10"]) == pytest.approx(1.0, abs=5e-3)
    assert float(values["tg_matched_r1000"]) == pytest.approx(1.0, abs=5e-3)


@pytest.mark.unit_test
def test_unit_test_two_population_shear_stress_is_right_at_any_density_ratio(report):
    """A shear wave in the heavy fluid decays at ``nu k^2`` without any source.

    The off-diagonal third moment is the one the enhanced equilibrium repairs,
    so the shear viscosity is right unaided even at a density ratio of 1000.
    """
    assert float(report["E8"]["tg_r1000_shear"]) == pytest.approx(1.0, abs=5e-3)
