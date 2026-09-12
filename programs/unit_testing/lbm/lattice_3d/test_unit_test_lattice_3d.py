"""Unit tests for the two three-dimensional lattices and the equilibrium on each.

``two_population_3d`` covers the solver on D3Q19. What is here is what exists
only because there are now two lattices, and it is built around one specific
way a second lattice goes wrong.

Four numbers change between D3Q19 and D3Q27: the rest weights of the moving
shells, the sound speed ``(c_s^k)^2``, and the amplitude ``lambda`` of the
enhanced equilibrium. The first three are visible in a table and a mistake in
them shows up immediately. ``lambda`` is derived, enters no conservation law,
and a wrong value leaves mass, momentum and the second moment all exactly right
while the shear viscosity is wrong by a factor that grows with the density
ratio. That is the bug the D2Q9-to-D3Q19 port shipped, and the correct D3Q27
value happens to equal the incorrect D3Q19 one -- so the same slip made in the
other direction would look, to every other check, entirely right.

The solver therefore never stores it: it sums the two lattice moments that fix
it and divides. These tests measure the moment that fixing was for,
``M3_xxy = rho_k (c_s^k)^2 u_y``, on both lattices at four density ratios.
"""

import pytest
from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "lattice_3d"

LATTICES = ("D3Q19", "D3Q27")

#: Relative error allowed in a lattice moment. Measured: 3e-16.
MOMENT_TOLERANCE = 1.0e-12

#: Relative mass drift allowed per fluid over six hundred steps. Measured: 3e-14.
MASS_DRIFT_TOLERANCE = 1.0e-9

#: How far outside [-1, 1] the phase field may stray. Measured: it does not.
PHASE_TOLERANCE = 1.0e-9

#: Relative error allowed on the predicted positivity bound. Measured: 0.08 %.
#:
#: The prediction drops the `(c_s^k)^2` terms in the amplitude, so it is the
#: high-ratio limit rather than the exact bound; at a ratio of 1000 the two
#: differ in the fourth digit.
BOUND_TOLERANCE = 0.01


@pytest.fixture(scope="module")
def report():
    """One run of the program, which takes about a minute."""
    return parse_key_values(run_unit_program(PROGRAM, timeout=600).stdout)


@pytest.mark.unit_test
@pytest.mark.parametrize("lattice", LATTICES)
def test_unit_test_lattice_weights_satisfy_their_conditions(report, lattice):
    """``sum w = 1``, ``sum w e e = I/3``, and fourth-order isotropy."""
    tag = f"lattice_{lattice}"
    assert float(report[f"{tag}_weight_sum_error"]) < MOMENT_TOLERANCE
    assert float(report[f"{tag}_weight_cs2_error"]) < MOMENT_TOLERANCE
    assert float(report[f"{tag}_weight_off_diagonal"]) < MOMENT_TOLERANCE
    assert float(report[f"{tag}_weight_fourth_error"]) < MOMENT_TOLERANCE


@pytest.mark.unit_test
def test_unit_test_only_d3q27_closes_the_sixth_order_relation(report):
    """``sum w e_x^4 e_y^2 = 3 sum w e_x^2 e_y^2 e_z^2``.

    This is the condition that fixed D3Q27's third shell against the other two,
    and it is the whole of what the corner velocities buy at the level of the
    weights. D3Q19 has no third shell and cannot satisfy it -- it misses by
    1/9 -- which is the concrete sense in which the smaller lattice is less
    isotropic. Asserting that D3Q19 *fails* it is not pedantry: if both passed,
    the relation would be measuring nothing.
    """
    assert float(report["lattice_D3Q27_weight_sixth_partial_error"]) < MOMENT_TOLERANCE
    assert float(report["lattice_D3Q19_weight_sixth_partial_error"]) > 0.1


@pytest.mark.unit_test
@pytest.mark.parametrize("lattice", LATTICES)
@pytest.mark.parametrize("alpha", ("02", "08"))
def test_unit_test_rest_weights_satisfy_their_conditions(report, lattice, alpha):
    """The rest weights, at two values of alpha well apart.

    ``sum phi = 1`` fixes the scale, ``sum phi e e`` must give the sound speed
    the solver derives its relaxation time from, and the fourth-order relation
    is what set the shells against each other. Two values of alpha rather than
    one, because the standard weights satisfy all of this at a single alpha by
    construction and would pass without the relations holding in general.
    """
    tag = f"lattice_{lattice}_alpha{alpha}"
    assert float(report[f"{tag}_rest_sum_error"]) < MOMENT_TOLERANCE
    assert float(report[f"{tag}_rest_cs2_error"]) < MOMENT_TOLERANCE
    assert float(report[f"{tag}_rest_off_diagonal"]) < MOMENT_TOLERANCE
    assert float(report[f"{tag}_rest_fourth_error"]) < MOMENT_TOLERANCE
    assert float(report[f"{tag}_cs2"]) > 0.0


@pytest.mark.unit_test
@pytest.mark.parametrize("lattice", LATTICES)
def test_unit_test_equilibrium_has_the_right_moments(report, lattice):
    """Mass, momentum, stress -- and the third moment, which is the point.

    The first three would pass with the amplitude set to anything at all,
    including zero. ``M3_xxy`` is the only moment that sees it, and it is
    checked at four density ratios spanning 1 to 1e5 on both lattices.
    """
    tag = f"equilibrium_{lattice}"
    assert float(report[f"{tag}_mass_error"]) < MOMENT_TOLERANCE
    assert float(report[f"{tag}_momentum_error"]) < MOMENT_TOLERANCE
    assert float(report[f"{tag}_stress_error"]) < MOMENT_TOLERANCE
    assert float(report[f"{tag}_third_moment_error"]) < MOMENT_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("lattice", LATTICES)
def test_unit_test_derived_amplitude_matches_its_closed_form(report, lattice):
    """The summed amplitude equals the algebra, on each lattice separately.

    ``3 (c_s^k)^2 - 1`` on D3Q19 and half of that on D3Q27. The solver reaches
    both by summing the same two moments over whichever velocity set it holds,
    so this is the check that the summation reproduces the hand derivation --
    and, because the two closed forms differ by exactly two, that the lattice
    the solver used is the lattice it was asked for.
    """
    assert float(report[f"equilibrium_{lattice}_amplitude_error"]) < MOMENT_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("lattice", LATTICES)
def test_unit_test_density_ratio_reaches_the_pressure(report, lattice):
    """``rho1/rho2 = (1 - alpha_2)/(1 - alpha_1)``, and the bulk pressures match.

    The coefficient in ``(c_s^k)^2`` differs between the lattices -- 1/2 against
    9/19 -- and it has to cancel between the two fluids, or the model would not
    have a single pressure and the density ratio would not be the one asked for.
    """
    tag = f"equilibrium_{lattice}"
    assert float(report[f"{tag}_density_ratio_error"]) < 1.0e-9
    assert float(report[f"{tag}_pressure_match_error"]) < 1.0e-9


@pytest.mark.unit_test
def test_unit_test_d3q27_leaves_more_room_before_the_equilibrium_goes_negative(report):
    """``(c_s^k)^2/3`` on D3Q19 against ``(c_s^k)^2/2`` on D3Q27.

    The two-population model buys its density ratio with the heavy fluid's
    sound speed, which leaves that fluid with very small moving weights; past
    the velocity measured here its equilibrium has a negative component. It is
    a real bound -- ``two_population_3d`` records the D3Q19 case crossing it
    during the transient at a ratio of 1000 -- and half again as much room is
    the practical reason to prefer the larger lattice, over and above the
    isotropy.

    Both bounds are minimised over the directions rather than assumed, so this
    also checks *which* direction binds: on D3Q27 the axial shell does, at
    ``1/2``, even though the corner shell has by far the smallest weight.
    """
    for lattice, expected in (("D3Q19", 1.0 / 3.0), ("D3Q27", 0.5)):
        tag = f"positivity_{lattice}"
        assert float(report[f"{tag}_bound_error"]) < BOUND_TOLERANCE
        assert float(report[f"{tag}_bound_over_cs2"]) == pytest.approx(expected, rel=BOUND_TOLERANCE)

    assert float(report["positivity_D3Q27_bound_over_cs2"]) > float(
        report["positivity_D3Q19_bound_over_cs2"]
    )


@pytest.mark.unit_test
@pytest.mark.parametrize("lattice", LATTICES)
def test_unit_test_solver_runs_on_both_lattices_at_a_ratio_of_1000(report, lattice):
    """Six hundred steps of a droplet, per fluid mass and phase field.

    What is *not* asserted is the spurious velocity of one lattice against the
    other. Six hundred steps is the transient and not the relaxed state, and
    measured there D3Q27 is marginally the worse of the two; the comparison the
    corner velocities were added for is the converged one, which belongs to
    ``laplace_3d``.
    """
    tag = f"droplet_{lattice}"
    assert float(report[f"{tag}_mass1_drift"]) < MASS_DRIFT_TOLERANCE
    assert float(report[f"{tag}_mass2_drift"]) < MASS_DRIFT_TOLERANCE
    assert float(report[f"{tag}_max_abs_phase"]) <= 1.0 + PHASE_TOLERANCE
    assert float(report[f"{tag}_spurious_velocity"]) < 0.05
