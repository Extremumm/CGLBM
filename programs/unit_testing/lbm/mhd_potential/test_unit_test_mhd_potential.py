"""Unit tests for the inductionless MHD potential solve and the current it makes.

The magnetic coupling is one elliptic solve and a cross product, so almost
everything that can go wrong with it goes wrong in the solve. These tests are
arranged by what kind of statement each one is, because they are not all the
same kind and conflating them is how an elliptic solver gets believed further
than it deserves.

**Exact statements.** That the linear solve recovers a potential it was handed
the residual of, and that the face current carries no net charge. Neither
involves truncation error: the first is the conjugate-gradient tolerance and
the second *is* the linear residual, because the current is built from the same
difference the operator was. Both are asserted near machine precision, and both
are run with a conductivity that jumps by 1e4 across an interface.

**A truncation statement.** That the discretisation is second order against an
analytic potential. This is the only test that says the operator solves the
right equation rather than merely solving some equation well.

**A measurement, not an assertion.** Whether the lattice-Boltzmann Poisson
march converges, and what it costs. It is a genuinely different discretisation
with a genuinely different failure mode, and the interesting numbers are its
own -- so they are recorded here with loose bounds, and the one place a hard
assertion is made is where the two paths must agree: both are second order
against the analytic solution.
"""

import pytest
from pycglbm.testing import parse_key_values, run_unit_program

PROGRAM = "mhd_potential"

#: Relative error allowed when recovering a manufactured potential.
#:
#: Measured: 4e-14 at a uniform conductivity and 3e-13 across a jump of 1e4.
#: The solve is asked for 1e-12 on the residual, and the error in the solution
#: trails it by the condition number, which is what the second number is.
SOLVE_TOLERANCE = 1.0e-9

#: Charge imbalance allowed, as a fraction of the current carrying it.
#:
#: Measured: 1.2e-12 periodic and 2.3e-12 against an insulating wall, at a
#: conductivity ratio of 1e4. This is the property the face-based current
#: exists for, so the bound is tight on purpose: a node-centred gradient would
#: put it at the truncation error of the lattice instead, which is 1e-3 here.
CHARGE_TOLERANCE = 1.0e-9

#: How far the measured order of accuracy may sit from two. Measured: 0.004.
ORDER_TOLERANCE = 0.1


@pytest.fixture(scope="module")
def report():
    """One run of the program. It takes about two minutes.

    Most of that is the lattice-Boltzmann march, which is asked to converge
    from cold at three resolutions and, once, on a problem it cannot converge
    on at all.
    """
    return parse_key_values(run_unit_program(PROGRAM, timeout=1800).stdout)


@pytest.mark.unit_test
@pytest.mark.parametrize("conductivity", ("uniform", "jump"))
def test_unit_test_linear_solve_recovers_a_manufactured_potential(report, conductivity):
    """Apply the operator to a chosen potential, solve, get it back.

    No truncation error enters this: the right-hand side was made by the same
    operator the solve inverts, so the only error is the one the conjugate
    gradient stopped at. The ``jump`` case runs it across a conductivity
    contrast of 1e4, which is a liquid metal beside an electrolyte and is what
    the Jacobi preconditioner has to survive -- it costs 41 iterations against
    the uniform case's 5.
    """
    tag = f"manufactured_{conductivity}"
    assert int(report[f"{tag}_converged"]) == 1
    assert float(report[f"{tag}_error"]) < SOLVE_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("boundary", ("periodic", "wall"))
def test_unit_test_finite_volume_current_conserves_charge(report, boundary):
    """``sum_faces J = 0`` at every node, to the solver tolerance.

    This is the one property that distinguishes this arrangement from the
    obvious one. ``div J = 0`` is not a diagnostic of the inductionless system;
    it is the equation that determines the potential, so a current evaluated
    with any other stencil carries a fictitious charge source, and the spurious
    force that source makes grows with the field strength that motivated the
    calculation in the first place.

    The ``wall`` case additionally says the insulating boundary holds: the
    current may not cross it, and dropping the wall face from the operator and
    the right-hand side alike makes that exact rather than approximate.
    """
    tag = f"charge_{boundary}_fv"
    assert int(report[f"{tag}_converged"]) == 1
    assert float(report[f"{tag}_current"]) > 0.0
    assert float(report[f"{tag}_relative"]) < CHARGE_TOLERANCE


@pytest.mark.unit_test
@pytest.mark.parametrize("solver", ("fv", "lbm"))
def test_unit_test_potential_is_second_order_accurate(report, solver):
    """Both discretisations, against the analytic potential of a known drive.

    With a uniform conductivity and ``B`` along z, a velocity
    ``(U sin k y, U sin k x, 0)`` drives ``phi = -(U B / k)[cos k x - cos k y]``
    exactly. Refining the lattice refines ``k`` with it, so the error should
    fall as ``n^-2``, and it does for both paths -- 2.00 for the finite volume
    and 1.99 for the lattice-Boltzmann march. This is what says the operator is
    of the right equation; everything else in this file would pass for an
    operator that was merely self-consistent.
    """
    for n in (32, 64):
        tag = f"analytic_{solver}_n{n}"
        assert int(report[f"{tag}_converged"]) == 1
        assert float(report[f"{tag}_order"]) == pytest.approx(2.0, abs=ORDER_TOLERANCE)


@pytest.mark.unit_test
def test_unit_test_lattice_boltzmann_march_converges_and_what_it_costs(report):
    """It converges at a uniform conductivity, diffusively.

    Recorded rather than pinned: 1152, 4272 and 15856 sweeps from cold at
    n = 16, 32 and 64, which is ``O(L^2)`` to within a per cent, against a
    single conjugate-gradient iteration for the same problem. The bound here is
    deliberately loose -- it is checking that the march terminates and that its
    cost scales the way a diffusive relaxation must, not reproducing a
    particular sweep count.
    """
    sweeps = [int(report[f"analytic_lbm_n{n}_iterations"]) for n in (16, 32, 64)]
    assert all(int(report[f"analytic_lbm_n{n}_converged"]) == 1 for n in (16, 32, 64))
    for coarse, fine in zip(sweeps, sweeps[1:]):
        assert 3.0 < fine / coarse < 5.0  # four, for a diffusive march


@pytest.mark.unit_test
def test_unit_test_a_large_conductivity_ratio_defeats_the_lattice_boltzmann_march(report):
    """And this is why the finite volume is the default.

    The march runs at a scaled diffusivity, but no scale helps: the ratio
    between the fastest and the slowest relaxing region is the conductivity
    ratio itself, so the poorly conducting fluid needs that factor more sweeps.
    At a ratio of 1e4 it has not converged after 40000 sweeps, and the current
    built from the unconverged potential carries a charge imbalance of about
    half of itself -- while the conjugate gradient solved the same problem in
    38 iterations to 1e-12.

    This is asserted, not merely recorded, because it is a real limit of the
    method and a future change that quietly made this pass would mean either
    the march had been fixed -- worth knowing -- or the test had stopped
    exercising the contrast it was written for.
    """
    assert int(report["charge_periodic_lbm_converged"]) == 0
    assert float(report["charge_periodic_lbm_relative"]) > 0.01
    assert int(report["charge_periodic_fv_converged"]) == 1
    assert float(report["charge_periodic_fv_relative"]) < CHARGE_TOLERANCE


@pytest.mark.unit_test
def test_unit_test_the_two_discretisations_agree_to_their_own_accuracy(report):
    """2.6 % apart on a swirling field at n = 24, and both converged.

    They cannot agree to round-off -- they are two discretisations of the same
    equation, second order and no better -- so what is checked is that they
    land within the truncation error each was separately measured to have, and
    not, say, a factor of two apart because one of them dropped a term.
    """
    assert int(report["agreement_fv_converged"]) == 1
    assert int(report["agreement_lbm_converged"]) == 1
    assert float(report["agreement_difference"]) < 0.1
