"""Unit tests for src/lbm/mixture and src/lbm/surface_force.

The high-density-ratio scheme rests on three facts checked here without running
a solver: the initial state built from a volume fraction is in equilibrium with
the equation of state, the phase field phi is a mass fraction whose zero sits
(W/2) ln(rho_1/rho_2) off the density interface, and the surface force -- the
divergence of the capillary stress -- integrates to the Laplace jump across a
circular interface and conserves momentum.
"""

import math

import pytest
from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "lbm_mixture"

#: Density ratios spanning the shipped case (20) to the high-ratio one (1e4).
DENSITY_RATIOS = ("1", "20", "1000", "10000")

STENCILS = ("E4", "E6", "E8")

#: Interface width of the solvers, as in the program.
WIDTH = 1.6


def _run(density_ratio="20", stencil="E8"):
    return parse_key_values(run_unit_program(PROGRAM, (density_ratio, stencil)).stdout)


@pytest.mark.unit_test
@pytest.mark.parametrize("density_ratio", DENSITY_RATIOS)
def test_unit_test_lbm_mixture_initial_state_is_in_equilibrium(density_ratio):
    """alpha -> (rho, phi) -> p gives back the pressure and the volume fraction.

    The relative pressure error grows with the density ratio because the heavy
    branch computes p as rho c^2 - p_inf, a difference of two numbers about
    3000 times larger than p at 1e4; it stays far below anything physical.
    """
    values = _run(density_ratio)
    assert float(values["round_trip_pressure_error"]) < 1e-8
    assert float(values["round_trip_alpha_error"]) < 1e-12
    assert float(values["round_trip_psi_error"]) < 1e-10


@pytest.mark.unit_test
@pytest.mark.parametrize("density_ratio", DENSITY_RATIOS)
def test_unit_test_lbm_mixture_bulks_are_pure_components(density_ratio):
    values = _run(density_ratio)
    assert float(values["pure1_phi"]) == 1.0
    assert float(values["pure2_phi"]) == -1.0
    assert float(values["pure1_rho_error"]) == 0.0
    assert float(values["pure2_rho_error"]) == 0.0


@pytest.mark.unit_test
@pytest.mark.parametrize("density_ratio", DENSITY_RATIOS)
def test_unit_test_lbm_mixture_normalised_phase_stays_bounded(density_ratio):
    """phi may overshoot +-1 by a rounding error; psi must stay within [-1, 1]."""
    assert float(_run(density_ratio)["psi_overshoot"]) <= 0.0


@pytest.mark.unit_test
@pytest.mark.parametrize("density_ratio", DENSITY_RATIOS)
def test_unit_test_lbm_mixture_phase_zero_sits_off_the_density_interface(density_ratio):
    """The phi = 0 contour is displaced by (W/2) ln(rho_1/rho_2).

    This is the "split" between the phase and density interfaces of the
    Laplace case: 2.4 lattice units at a density ratio of 20, 7.4 at 1e4. It is
    a property of phi being a mass fraction, not an error of the solver.
    """
    values = _run(density_ratio)
    offset = float(values["phase_zero_offset"])
    assert offset == pytest.approx(float(values["phase_zero_offset_expected"]), abs=1e-9)
    # and close to the closed form with the nominal ratio; the program takes
    # rho_k at the ambient pressure, where the Laplace offset of p1_inf makes
    # rho_1 slightly lighter than nominal (19.92 instead of 20)
    ratio = float(density_ratio)
    if ratio >= 20:
        assert offset == pytest.approx(0.5 * WIDTH * math.log(ratio), rel=5e-3)


@pytest.mark.unit_test
def test_unit_test_lbm_mixture_single_viscosity_is_reproduced():
    """nu_1 == nu_2 must give that viscosity back, whatever phi."""
    assert float(_run()["viscosity_single_error"]) == 0.0


@pytest.mark.unit_test
@pytest.mark.parametrize("density_ratio", DENSITY_RATIOS)
def test_unit_test_lbm_mixture_viscosity_mixes_by_volume(density_ratio):
    """rho (Y_1 nu_1 + Y_2 nu_2) equals alpha_1 mu_1 + alpha_2 mu_2."""
    assert float(_run(density_ratio)["viscosity_mixing_error"]) < 1e-12


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_lbm_mixture_surface_force_gives_the_laplace_jump(stencil):
    """Summed across the interface the force is sigma/R, pointing inward.

    What is left is the finite interface width at R = 20: the force is spread
    over radii R +- W, where 1/r is not exactly 1/R.
    """
    values = _run(stencil=stencil)
    assert float(values["force_jump_over_laplace"]) == pytest.approx(1.0, abs=1.5e-2)
    assert int(values["force_outward_nodes"]) == 0
    assert int(values["force_inward_nodes"]) > 0


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_lbm_mixture_surface_force_conserves_momentum(stencil):
    """An asymmetric closed interface exerts no net force on the fluid.

    The force is the divergence of a stress, and a centred divergence sums to
    zero over a periodic lattice. The curvature form sigma/2 kappa grad(psi)
    does not have this property; see src/lbm/surface_force.h.
    """
    assert float(_run(stencil=stencil)["ellipse_net_force"]) < 1e-13


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_lbm_mixture_flat_interface_between_walls_feels_no_force(stencil):
    assert float(_run(stencil=stencil)["flat_force"]) < 1e-14
