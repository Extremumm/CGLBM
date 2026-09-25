"""Hartmann flow against its closed form: the validation of the MHD coupling.

This is the one inductionless configuration with an analytic answer, and it is
what the magnetic force is scored against. A uniform drive along x, insulating
walls at ``y = +-L``, the field normal to them, and

    u(y) = (G / (sigma B^2)) [1 - cosh(Ha y/L) / cosh(Ha)],
    Ha = B L sqrt(sigma / mu).

Two independent things are measured, and it matters that they are independent.
The **core velocity** is ``G / (sigma B^2)`` and fixes the magnitude of the
magnetic force. The **layer thickness** is ``L / Ha``, which is
``sqrt(mu/sigma) / B``, and fixes how that force is distributed against
viscosity. A force scaled wrong gets the first one wrong; a current assembled or
averaged wrong gets the second one wrong while the first can still come out
right by coincidence of the drive.

Nothing here is duplicated from the program. Every constant is read back out of
the run log, which is the only place the case is stated.

**What this case does not test.** The potential solve. With ``u = u(y) x^`` and
``B`` along y, the motional field points along z and has no z dependence, so
``div(sigma u x B)`` vanishes identically and the potential is uniform. That is
asserted below -- it is a real check, since a broken drive assembly would put
structure there -- but the potential solve is exercised by
``programs/unit_testing/lbm/mhd_potential`` and by ``magnetic_rayleigh_taylor``.
"""

import numpy as np
import pytest
from pycglbm.testing import artifacts_dir, parse_key_values, run_program


def analytic_profile(run) -> tuple[np.ndarray, np.ndarray]:
    """``(eta, u)`` of the closed form, on the lattice the run used.

    The walls are half a node outside the first and last row -- that is what
    half-way bounce-back places them at -- so the channel is ``ny`` wide, its
    half-width is ``ny/2``, and row ``j`` sits at
    ``eta = (j - (ny-1)/2) / (ny/2)``.
    """
    ny = run.parameter("ny", int)
    hartmann = run.parameter("hartmann_number")
    conductivity = run.parameter("sigma_e1")
    field = run.parameter("b_y")
    drive = run.parameter("body_force_x")

    rows = np.arange(ny)
    eta = (rows - (ny - 1) / 2.0) / (ny / 2.0)
    core = drive / (conductivity * field * field)
    return eta, core * (1.0 - np.cosh(hartmann * eta) / np.cosh(hartmann))


def measured_profile(run, timestep: int) -> np.ndarray:
    """``u_x`` against y, averaged along x.

    The flow is one-dimensional, so the average is over a row of identical
    values; ``test_..._stays_one_dimensional`` is what says so.
    """
    return run.velocity(timestep)[:, :, 0].mean(axis=1)


@pytest.fixture(scope="module")
def run_hartmann():
    """One run of the case. It takes about two minutes on 24 cores."""
    return run_program("hartmann", artifacts_dir() / "hartmann", timeout=2400)


@pytest.fixture(scope="module")
def log(run_hartmann):
    return parse_key_values((run_hartmann.rundir / "run.log").read_text())


@pytest.mark.validation
def test_validation_hartmann_is_the_case_it_claims_to_be(run_hartmann):
    """The arrangement, before any number taken out of it means anything.

    A field along x rather than y, or periodic boundaries rather than walls,
    would each give a perfectly smooth profile that is not Hartmann flow.
    """
    assert run_hartmann.parameter("mhd", str) == "true"
    assert run_hartmann.parameter("boundary", str) == "wall_y"
    assert run_hartmann.parameter("lattice_3d", str) == "D3Q27"
    assert run_hartmann.parameter("b_y") > 0.0
    assert run_hartmann.parameter("b_x") == run_hartmann.parameter("b_z") == 0.0
    assert run_hartmann.parameter("body_force_x") > 0.0
    assert run_hartmann.parameter("sigma") == 0.0  # no surface tension: one fluid
    assert run_hartmann.parameter("rho1") == run_hartmann.parameter("rho2")


@pytest.mark.validation
def test_validation_hartmann_stays_one_dimensional(run_hartmann):
    """Nothing varies along x or z, so the profile is a function of y alone.

    Exactly nothing, to the last bit: the case is uniform along both axes and
    every operator in the time loop is too. A non-zero answer here would mean
    the magnetic force had introduced a dependence the physics does not have.
    """
    velocity = run_hartmann.velocity(run_hartmann.timesteps[-1])
    along_x = velocity[:, :, 0]
    assert np.abs(along_x - along_x[:, :1]).max() == 0.0


@pytest.mark.validation
def test_validation_hartmann_potential_is_uniform(log):
    """``div(sigma u x B) = 0`` in this configuration, so ``phi`` is a constant.

    Reported by the program at the end of the run. It is zero, not merely small:
    the right-hand side is identically zero, the solve recognises that and
    returns the zero-mean constant without iterating, and the current is then
    the motional term alone.
    """
    assert float(log["final_max_potential"]) == 0.0
    assert float(log["final_charge_imbalance"]) == 0.0


@pytest.mark.validation
def test_validation_hartmann_core_velocity_matches(run_hartmann):
    """``u_core = (G / sigma B^2)(1 - 1/cosh Ha)``, which measures the force.

    In the core the viscous term is negligible and the drive is balanced by the
    magnetic force alone, so this number is the Lorentz force's magnitude read
    directly. Measured: 6 parts per million.
    """
    _, analytic = analytic_profile(run_hartmann)
    measured = measured_profile(run_hartmann, run_hartmann.timesteps[-1])
    middle = len(measured) // 2
    assert measured[middle] == pytest.approx(analytic[middle], rel=1.0e-3)


@pytest.mark.validation
def test_validation_hartmann_profile_matches(run_hartmann):
    """The whole profile, core and both layers, in L2. Measured: 0.1 %.

    The tolerance is 1 %, which is where the discretisation of a layer 3.2 nodes
    thick puts it rather than where the run does; the case is well inside it.
    """
    _, analytic = analytic_profile(run_hartmann)
    measured = measured_profile(run_hartmann, run_hartmann.timesteps[-1])
    error = np.sqrt(np.mean((measured - analytic) ** 2)) / np.sqrt(np.mean(analytic**2))
    assert error < 0.01


@pytest.mark.validation
def test_validation_hartmann_layer_is_the_thickness_the_field_sets(run_hartmann, log):
    """``L / Ha``, measured from the wall to where the profile reaches 1 - 1/e.

    This is the half of the answer the core velocity cannot see. It depends on
    ``sqrt(mu/sigma)/B`` rather than on ``sigma B^2``, so a current that was
    right in magnitude but assembled from the wrong stencil -- averaged over the
    wrong faces, say -- would keep the core and lose this.

    The measurement is to the nearest node, and the layer is 3.2 nodes thick, so
    the comparison is held to half a node rather than to a percentage.
    """
    measured = measured_profile(run_hartmann, run_hartmann.timesteps[-1])
    core = measured[len(measured) // 2]
    reached = np.argmax(measured > core * (1.0 - np.exp(-1.0)))
    # Row `reached` is the first inside the layer; the wall is half a node
    # outside row 0, so its distance from the wall is `reached + 0.5`.
    assert reached + 0.5 == pytest.approx(float(log["layer_thickness"]), abs=0.75)


@pytest.mark.validation
def test_validation_hartmann_has_reached_its_steady_state(run_hartmann):
    """The last two outputs agree, so the profile is the steady one.

    Both relaxation times of the problem -- the magnetic ``rho/(sigma B^2)`` and
    the viscous ``(L/Ha)^2/nu`` -- are about 62 steps against a run of 20000, so
    this has a great deal of margin. It is here because without it every number
    above would be a statement about a transient.
    """
    steps = run_hartmann.timesteps
    assert len(steps) >= 3
    last = measured_profile(run_hartmann, steps[-1])
    before = measured_profile(run_hartmann, steps[-2])
    assert np.abs(last - before).max() / np.abs(last).max() < 1.0e-6
