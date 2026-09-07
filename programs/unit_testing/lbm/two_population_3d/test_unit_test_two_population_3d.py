"""Unit tests for the three-dimensional colour-gradient solver.

The model is the one in ``two_population_solver.h`` with an extra dimension, so
most of what has to hold is the same: each fluid's mass conserved separately,
both distributions non-negative, a uniform fluid staying at rest. Three things
are new in three dimensions, and each is checked directly because each is a
place a port goes wrong quietly rather than loudly.

The gradient stencils claim an isotropy order, and the shell weights were solved
for it rather than taken from a table. The rest weights ``phi_i^k`` change --
``(1-alpha)/12`` and ``(1-alpha)/24`` rather than ``(1-alpha)/5`` and
``(1-alpha)/20`` -- and with them the sound speed, ``(c_s^k)^2 = (1-alpha)/2``,
which is how the density ratio reaches the pressure. And the curvature of a
sphere is ``-2/R``, not ``-1/R``: Laplace's law reads ``2 sigma / R`` here, and
getting it wrong would surface only as a factor of two in a validation run.

Two more were added with the cases that move. ``oscillation_3d`` starts from a
spheroid and measures its semi-axes back, and either half of that could hide an
error in the other, so the shape and the measurement are checked against the
analytic surface here. ``rayleigh_taylor_3d`` uses the wall and the body force
together in a flow that has no closed form, so the wall and the body force are
separated here: mass may not cross a wall, and a fluid at rest between two of
them must go hydrostatic.
"""

import pytest
from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "two_population_3d"

STENCILS = ("E4", "E6")

#: Density ratios the solver is exercised at, as they appear in the report.
RATIOS = ("r20", "r1000")

#: Relative error allowed in the equilibrium's moments. Measured: 1e-16.
MOMENT_TOLERANCE = 1.0e-12

#: Relative mass drift allowed per fluid over thirty steps. Measured: 6e-14.
MASS_DRIFT_TOLERANCE = 1.0e-9

#: How far outside [-1, 1] the phase field may stray. Measured: it does not.
PHASE_TOLERANCE = 1.0e-12

#: A uniform fluid must stay exactly at rest. Measured: 0.
REST_SPEED_TOLERANCE = 1.0e-12

#: Relative error allowed on the curvature of a resolved sphere. Measured: 0.06 %.
CURVATURE_TOLERANCE = 0.02
#: Relative error allowed on the semi-axes of an analytic spheroid.
#:
#: Measured: 0.13 % on the equator and 0.18 % at the pole, both positive. The
#: bias is the linear interpolation between the two nodes bracketing the
#: crossing: tanh is convex on the inside of the interface, so the straight line
#: between the samples crosses zero slightly late.
AXIS_TOLERANCE = 0.005

#: Relative error allowed on `dp/dy = -rho g`. Measured: 1.0 % after 4000 steps.
#:
#: It is a residual, not a floor. The whole pressure varies by 0.8 % across the
#: box, so the balance is being read out of the fifth digit of the pressure, and
#: the run still carries motion at 4e-6 -- enough to move that digit by about a
#: per cent. What the test is for is the failure that would matter: a force that
#: reached the momentum equation with the wrong sign or magnitude, or a wall
#: that pushed back on the column.
HYDROSTATIC_TOLERANCE = 0.03


@pytest.fixture(scope="module")
def report():
    """One run of the unit program per stencil, shared by the tests below."""
    return {s: parse_key_values(run_unit_program(PROGRAM, (s,)).stdout) for s in STENCILS}


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_stencils_are_normalised(report, stencil):
    """Both stencils must differentiate exactly: sum W c_x^2 = 1, no cross term."""
    values = report[stencil]
    for name in STENCILS:
        assert float(values[f"stencil_{name}_second_moment"]) == pytest.approx(1.0, abs=1e-12)
        assert float(values[f"stencil_{name}_off_diagonal"]) < 1.0e-14


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_stencils_have_the_isotropy_they_claim(report, stencil):
    """E4 is fourth-order isotropic and E6 is additionally sixth-order.

    The weights are 1/6, 1/12 for E4 and 2/15, 1/15, 1/60, 1/120 for E6, and
    they are the solution of exactly these conditions -- so this is a check of
    the derivation, not a tolerance on a measurement. E4's sixth-order error is
    asserted to be *large*: it is a fourth-order stencil and claiming otherwise
    would be the mistake.
    """
    values = report[stencil]
    assert float(values["stencil_E4_fourth_error"]) < 1.0e-14
    assert float(values["stencil_E4_sixth_error"]) > 0.1
    assert float(values["stencil_E6_fourth_error"]) < 1.0e-14
    assert float(values["stencil_E6_sixth_error"]) < 1.0e-14


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_equilibrium_has_the_right_moments(report, stencil):
    """Mass, momentum and momentum flux, at density ratios from 1 to 1e5."""
    values = report[stencil]
    assert float(values["equilibrium_mass_error"]) < MOMENT_TOLERANCE
    assert float(values["equilibrium_momentum_error"]) < MOMENT_TOLERANCE
    assert float(values["equilibrium_stress_error"]) < MOMENT_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_rest_weights_give_the_d3q19_sound_speed(report, stencil):
    """sum_q phi_q^k = 1 and sum_q phi_q^k e_x^2 = (1 - alpha_k) / 2.

    The second is the D3Q19 relation, against 3(1 - alpha)/5 on D2Q9. Both
    follow from requiring a fourth-order isotropic fourth moment of phi, and
    both are what tie the density ratio to the rest weight.
    """
    assert float(report[stencil]["rest_weight_error"]) < MOMENT_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_rest_weights_encode_the_density_ratio(report, stencil):
    """(1 - alpha_2) / (1 - alpha_1) must be rho1 / rho2, as in two dimensions."""
    assert float(report[stencil]["alpha_density_ratio_error"]) < 1.0e-9


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_sphere_curvature_is_two_over_r(report, stencil):
    """The discrete curvature of a resolved sphere must be -2/R.

    This is the one place the extra dimension changes the physics rather than
    the bookkeeping: Laplace's law is 2 sigma / R for a sphere and sigma / R for
    a cylinder, and a curvature operator ported without thought would give the
    latter.
    """
    values = report[stencil]
    measured = float(values["sphere_curvature"])
    exact = float(values["sphere_curvature_exact"])
    assert measured == pytest.approx(exact, rel=CURVATURE_TOLERANCE)


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
@pytest.mark.parametrize("ratio", RATIOS)
def test_unit_test_two_population_3d_conserves_each_fluid_separately(report, stencil, ratio):
    """Neither fluid may leak into the other through the recolouring."""
    values = report[stencil]
    assert float(values[f"{ratio}_mass1_drift"]) < MASS_DRIFT_TOLERANCE
    assert float(values[f"{ratio}_mass2_drift"]) < MASS_DRIFT_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
@pytest.mark.parametrize("ratio", RATIOS)
def test_unit_test_two_population_3d_keeps_both_distributions_non_negative(report, stencil, ratio):
    """The adapted segregation operator must hold in three dimensions too.

    It is the same argument: the push is proportional to the mixture's own rest
    weight, so near equilibrium it is ``beta rho_1 / rho`` of the other fluid's
    share whatever the density ratio. Only the numbers in ``phi_q^k`` change.
    """
    assert float(report[stencil][f"{ratio}_min_population"]) >= 0.0


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
@pytest.mark.parametrize("ratio", RATIOS)
def test_unit_test_two_population_3d_keeps_the_phase_field_bounded(report, stencil, ratio):
    """phi_N stays in [-1, 1], as a consequence of the test above."""
    assert float(report[stencil][f"{ratio}_max_abs_phase"]) <= 1.0 + PHASE_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_leaves_a_uniform_fluid_at_rest(report, stencil):
    """No interface and no gravity: nothing may start moving."""
    assert float(report[stencil]["rest_max_speed"]) < REST_SPEED_TOLERANCE



@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_lays_down_the_spheroid_it_was_asked_for(report, stencil):
    """The mode-2 shape and the semi-axes that measure it must agree.

    `oscillation_3d` prescribes `r(theta) = R' [1 + eps P_2(cos theta)]` and
    reads the oscillation off the difference between the polar and equatorial
    semi-axes. If the deformation came out at the wrong amplitude, or the
    measurement read it back wrong, the case would still oscillate and would
    still be scored against Lamb -- so both are checked against the analytic
    surface, which fixes the poles at `R'(1 + eps)` and the equator at
    `R'(1 - eps/2)`.
    """
    values = report[stencil]
    assert float(values["spheroid_rz"]) == pytest.approx(
        float(values["spheroid_polar_exact"]), rel=AXIS_TOLERANCE
    )
    for axis in ("spheroid_rx", "spheroid_ry"):
        assert float(values[axis]) == pytest.approx(
            float(values["spheroid_equatorial_exact"]), rel=AXIS_TOLERANCE
        )
    assert float(values["spheroid_equal_volume_radius"]) == pytest.approx(
        float(values["spheroid_equal_volume_exact"]), rel=AXIS_TOLERANCE
    )


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_spheroid_is_axisymmetric(report, stencil):
    """The two equatorial semi-axes must be equal to the last bit.

    The deformation is about z, so x and y are interchangeable. Nothing in the
    lattice enforces that -- D3Q19 treats the three axes alike, and this is the
    statement that the solver does too.
    """
    values = report[stencil]
    assert float(values["spheroid_rx"]) == float(values["spheroid_ry"])


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_walls_hold_the_fluid(report, stencil):
    """Nothing may cross a wall, and the fluid must come to rest against it."""
    values = report[stencil]
    assert float(values["wall_mass_drift"]) < MASS_DRIFT_TOLERANCE
    assert float(values["wall_max_speed"]) < 1.0e-4


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_gravity_reaches_hydrostatic_balance(report, stencil):
    """dp/dy must settle at -rho g, which is what the body force is for."""
    assert float(report[stencil]["wall_hydrostatic_error"]) < HYDROSTATIC_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_wall_layer_conserves_both_components(report, stencil):
    """The Rayleigh-Taylor arrangement, with the walls and the force together."""
    values = report[stencil]
    assert float(values["layer_mass1_drift"]) < MASS_DRIFT_TOLERANCE
    assert float(values["layer_mass2_drift"]) < MASS_DRIFT_TOLERANCE
    assert float(values["layer_min_population"]) >= 0.0
    assert float(values["layer_max_abs_phase"]) <= 1.0 + PHASE_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_two_population_3d_treats_x_and_z_alike(report, stencil):
    """A layer symmetric under swapping x and z must stay so.

    The perturbation is `cos(2 pi x / nx) cos(2 pi z / nz)` on a square
    cross-section, which is unchanged by the swap. Everything downstream of it
    -- the gradient stencil, the streaming, the wall -- has to be too, and this
    is the one invariant a two-dimensional run could not have checked.
    """
    assert float(report[stencil]["layer_xz_asymmetry"]) < 1.0e-12
