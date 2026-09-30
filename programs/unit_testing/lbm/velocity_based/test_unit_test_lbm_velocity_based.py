"""Unit tests for src/lbm/velocity_based.

The scheme stands on a few exact properties, checked here without running a
solver: the moments of its equilibria and forcing, the conservation of P and u
by the collisions, the positivity of the phase-field carrier, a pressure force
that leaves no force from a uniform pressure across a density jump of 1e4, a
link momentum exchange that is equal and opposite at the two ends of every
link, a dissipation force that sums to zero, and a phase transport that keeps
its profile, its mass and its bounds, at speeds up to 0.2, and with its
fourth-order corrections the shape of a curved interface.
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
def test_unit_test_lbm_velocity_based_hybrid_collision(values):
    """The hybrid collision is the regularised one at sigma = 1, conserves P
    and u, and does not depend on sigma for a resolved non-equilibrium."""
    assert float(values["hybrid_unit_weight_error"]) == 0.0
    assert float(values["hybrid_conservation_error"]) < 1e-15
    assert float(values["hybrid_resolved_error"]) < 1e-18


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_filtered_collision(values):
    """The filtered collision is the regularised one at sigma = 1, conserves P
    and u, leaves a non-equilibrium that has not changed since the step before
    alone whatever sigma, and at sigma = 0 lets nothing of one that has flipped
    sign through."""
    assert float(values["filtered_unit_weight_error"]) == 0.0
    assert float(values["filtered_steady_error"]) < 1e-18
    assert float(values["filtered_flip_error"]) < 1e-18
    assert float(values["filtered_conservation_error"]) < 1e-15


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_carrier_is_positive(values):
    """Gamma_i(u) >= 0 up to |u| = 0.1, which is what bounds the phase field."""
    for speed in (1, 5, 10):
        assert float(values[f"carrier_min_{speed}"]) > 0.0


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_phase_carrier_at_lower_temperature(values):
    """The phase carrier at T = 0.2 has the moments 1, u and T I + u u exactly.

    The second-order Hermite form on D2Q9 has them only at cs^2; away from it
    the second moment is off by O(u^2), which moved a sheared interface by
    1e-6 nodes per step, so the carrier puts it back. Its populations stay
    non-negative up to |u| = 0.1, it is velocity_equilibrium itself at cs^2,
    the phase populations built on it carry c and c u + A n at mobility T / 2,
    and the link exchange read with it is still equal and opposite.
    """
    assert float(values["phase_carrier_moment_error"]) < 1e-15
    assert float(values["phase_carrier_min_10"]) > 0.0
    assert float(values["phase_carrier_cs2_identity"]) == 0.0
    assert float(values["phase_carrier_populations_error"]) < 1e-15
    assert float(values["phase_carrier_link_antisymmetry"]) < 1e-14


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_uniform_pressure_exerts_no_force(values):
    """Across a density jump of 1e4, a uniform p gives no force, to rounding.

    P = p / (rho cs^2) itself jumps by 1e4 there; the force is the lattice
    gradient of p, so nothing of the order of the jump has to cancel.
    """
    assert float(values["uniform_pressure_force"]) < 1e-15


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_pressure_force_conserves(values):
    """As a force density, the pressure force sums to zero over the lattice."""
    assert float(values["pressure_force_total"]) < 1e-13


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_link_momentum_conserves(values):
    """What one end of a link gains, the other loses; over a periodic lattice
    the exchanges sum to zero. Both ends see the same dissipation coefficient,
    and it is never negative."""
    assert float(values["link_antisymmetry"]) < 1e-14
    assert float(values["lattice_momentum_total"]) < 1e-14
    assert float(values["link_dissipation_asymmetry"]) < 1e-12
    assert float(values["link_dissipation_min"]) >= 0.0


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_dissipation_force(values):
    """With one coefficient per link, the dissipation force sums to zero over
    the lattice, and a uniform velocity feels none."""
    assert float(values["dissipation_force_total"]) < 1e-13
    assert float(values["dissipation_force_uniform"]) == 0.0


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_link_momentum_uniform_velocity(values):
    """A uniform velocity is carried by the phase populations' mass flux
    alone, whatever the densities on the two sides: no spurious exchange."""
    assert float(values["link_uniform_velocity_error"]) < 1e-12


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_link_momentum_single_component(values):
    """Within one component the exchange is the lattice's own, with its
    advection written as a mass flux; no dissipation."""
    assert float(values["link_single_component_error"]) < 1e-14
    assert float(values["link_single_component_dissipation"]) < 1e-12


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_set_velocity(values):
    """set_velocity changes the first moment and leaves the others alone."""
    assert float(values["set_velocity_error"]) < 1e-15


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


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_phase_limiter(values):
    """At |u| = 0.2 the unlimited sharpening would make populations negative;
    limited, each stays between 0 and its carrier, and c is kept."""
    assert float(values["phase_unlimited_min"]) < -1e-4
    assert float(values["phase_limited_min"]) > -1e-17
    assert float(values["phase_limited_complement_min"]) > -1e-17
    assert float(values["phase_limited_mass_error"]) < 1e-15


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_fourth_order_operators(values):
    """On a quartic the lattice Laplacian is del^2 f + del^4 f / 12, the
    interface normal is grad f / |grad f| and the correction flux is
    -(T/24) grad del^2 f, to rounding; the plain E4 normal is 2 % off."""
    assert float(values["lattice_laplacian_error"]) < 1e-13
    assert float(values["interface_normal_error"]) < 1e-13
    assert float(values["correction_flux_error"]) < 1e-13
    assert float(values["plain_normal_error"]) > 1e-2


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_fourth_order_phase_keeps_the_shape(values):
    """A mode-2 droplet of radius 10, its phase field alone with the fluid at
    rest: over 2000 steps the E8 normal loses 9 % of the deformation to
    surface diffusion, the fourth-order phase, its normal confined to the
    interface as the solver takes it, 0.1 %, conserving c and staying
    non-negative."""
    plain = float(values["plain_phase_shape_kept"])
    corrected = float(values["fourth_order_phase_shape_kept"])
    assert plain < 0.95
    assert abs(corrected - 1.0) < 0.1 * (1.0 - plain)
    assert float(values["fourth_order_phase_mass_error"]) < 1e-12
    assert float(values["fourth_order_phase_min"]) >= 0.0


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_sixth_order_operators(values):
    """On a plane wave, halving the wavenumber divides the fourth-order
    normal's angle and transport residual by 16 and the sixth-order ones by
    64: the sixth-order normal is exact to sixth order, and the sixth-order
    flux leaves the transport an error of eighth order (sixth relative to the
    diffusion)."""
    assert float(values["normal4_convergence"]) == pytest.approx(16.0, rel=0.05)
    assert float(values["normal6_convergence"]) == pytest.approx(64.0, rel=0.05)
    assert float(values["transport4_convergence"]) == pytest.approx(16.0, rel=0.05)
    assert float(values["transport6_convergence"]) == pytest.approx(64.0, rel=0.05)
    assert float(values["normal6_angle"]) < 0.1 * float(values["normal4_angle"])
    assert float(values["transport6_error"]) < 0.1 * float(values["transport4_error"])


@pytest.mark.unit_test
def test_unit_test_lbm_velocity_based_sixth_order_phase_keeps_the_shape(values):
    """The phase-only droplet of the fourth-order test with the sixth-order
    phase: it keeps its deformation to 0.3 %, a little less well than the
    fourth order's 0.1 %, conserving c and staying non-negative."""
    plain = float(values["plain_phase_shape_kept"])
    sixth = float(values["sixth_order_phase_shape_kept"])
    assert abs(sixth - 1.0) < 0.1 * (1.0 - plain)
    assert float(values["sixth_order_phase_mass_error"]) < 1e-12
    assert float(values["sixth_order_phase_min"]) >= 0.0
