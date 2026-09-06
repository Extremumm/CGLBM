"""Unit tests for src/lbm/isotropic_gradient.

Each stencil is checked against the property it is built to satisfy: exactness
on a linear field, and isotropy of the lattice tensors up to its stated order.
The last test is the one that matters physically -- the colour gradient enters
the surface-tension operator divided by its own norm, so any angular error
becomes an anisotropic surface tension.
"""

import pytest

from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "lbm_gradient"

STENCILS = ("E4", "E6", "E8")

#: Isotropy order each stencil is constructed to reach.
ISOTROPY_ORDER = {"E4": 4, "E6": 6, "E8": 8}

#: Neighbours in each stencil, and how far it reaches along an axis.
GEOMETRY = {"E4": (8, 1), "E6": (12, 2), "E8": (24, 2)}


def _run(stencil):
    return parse_key_values(run_unit_program(PROGRAM, (stencil,)).stdout)


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_lbm_gradient_geometry(stencil):
    values = _run(stencil)
    points, reach = GEOMETRY[stencil]
    assert values["stencil"] == stencil
    assert int(values["points"]) == points
    assert int(values["reach"]) == reach


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_lbm_gradient_is_normalised(stencil):
    """sum W c_a c_b must be the identity, and the odd moments must vanish."""
    values = _run(stencil)
    assert float(values["moment_xx"]) == pytest.approx(1.0, abs=1e-12)
    assert float(values["moment_yy"]) == pytest.approx(1.0, abs=1e-12)
    assert float(values["moment_xy"]) == pytest.approx(0.0, abs=1e-12)
    assert float(values["moment_x"]) == pytest.approx(0.0, abs=1e-12)
    assert float(values["moment_y"]) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_lbm_gradient_is_exact_on_a_linear_field(stencil):
    """Normalisation makes every stencil differentiate a linear field exactly."""
    assert float(_run(stencil)["linear_error"]) < 1e-12


@pytest.mark.unit_test
@pytest.mark.parametrize("stencil", STENCILS)
def test_unit_test_lbm_gradient_reaches_its_isotropy_order(stencil):
    """Isotropic up to its own order, and measurably not beyond it."""
    values = _run(stencil)
    order = ISOTROPY_ORDER[stencil]

    for satisfied in range(4, order + 1, 2):
        defect = float(values[f"isotropy_defect_{satisfied}"])
        assert defect < 1e-12, f"{stencil} should be isotropic at rank {satisfied}"

    if order < 8:
        # the first order it does not reach must be visibly broken, otherwise
        # the stencil is not what it claims to be
        defect = float(values[f"isotropy_defect_{order + 2}"])
        assert defect > 1e-3, f"{stencil} should not be isotropic at rank {order + 2}"


@pytest.mark.unit_test
def test_unit_test_lbm_gradient_anisotropy_falls_with_order():
    """On a radial tanh interface, a higher-order stencil points more truly.

    This is the quantity Leclaire et al. (Computers & Fluids 48, 98, 2011)
    report improving by about an order of magnitude, and the reason the wider
    stencils exist at all.
    """
    errors = {s: float(_run(s)["radial_angle_error_deg"]) for s in STENCILS}

    assert errors["E8"] < errors["E6"] < errors["E4"]
    # E4 is already good to a fraction of a degree; E8 must be an order better
    assert errors["E4"] / errors["E8"] > 5.0
    assert errors["E8"] < 0.05
