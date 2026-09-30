"""A droplet carried by a uniform flow: the Galilean-invariance test.

The shipped Laplace case (a droplet of radius 10 at a density ratio of 20, the
same kinematic viscosity in both fluids) with the whole box started at
U = 0.01 along x:

    laplace --translate=0.01 --surface-tension=stress --viscosity-mixing=dynamic
            --source-stencil=isotropic|matched --steps=3000 --interval=500

In a periodic box nothing should happen but the translation: the droplet moves
at U and the fluid around it too. With the capillary stress the total momentum
is kept to rounding, so what can go wrong is that it is shared out wrongly.

With the nine-point source it is. The droplet falls behind the flow, at 0.92 U,
and the fluid around it runs ahead to keep the total: the deviatoric part of
S_Sp, which cancels D2Q9's third-moment error, is taken on a stencil the
streaming did not make the error with, and across the moving interface the
difference is a stress dipole that holds the heavy side back and pushes the
light side on. It does not depend on the viscosity, since the error stress and
the drag that balances it both scale with it. With the matched stencil the
droplet moves at U to within 0.1 %. docs/numerics.md, "Generalized equilibria
in central moments", has the measurements, and what is left at a lower
viscosity.

The case's default tension, the continuum surface force, does not conserve
momentum on a moving droplet (the box loses 1.5 % of it in 2000 steps), which is
why this runs the capillary stress.
"""

import numpy as np
import pytest
from pycglbm import CaseOutput
from pycglbm.testing import artifacts_dir, run_program

SPEED = 0.01
STEPS = 3000
#: The droplet's speed is measured between these two outputs, after the
#: start-up has passed.
FIRST, LAST = 1000, 3000

#: The droplet's speed over U, measured, per source stencil.
MEASURED_SPEED = {"isotropic": 0.919, "matched": 0.9997}
SPEED_TOLERANCE = 0.003


@pytest.fixture(scope="module", params=sorted(MEASURED_SPEED))
def translating_run(request) -> CaseOutput:
    stencil = request.param
    run = run_program(
        "laplace",
        artifacts_dir() / f"laplace_translating_{stencil}",
        args=(
            f"--translate={SPEED!r}",
            "--surface-tension=stress",
            "--viscosity-mixing=dynamic",
            f"--source-stencil={stencil}",
            f"--steps={STEPS}",
            "--interval=500",
        ),
        timeout=1800,
    )
    run.stencil = stencil
    return run


def centroid_x(run: CaseOutput, timestep: int) -> float:
    """x of the droplet's centre, from the phase of the first Fourier moment of
    its volume fraction: exact for a shape symmetric about its centre, and
    indifferent to the periodic wrap."""
    rho = run.density(timestep)
    rho1, rho2 = run.parameter("rho1"), run.parameter("rho2")
    fraction = (rho - rho2) / (rho1 - rho2)
    nx = rho.shape[1]
    angle = 2.0 * np.pi * np.arange(nx) / nx
    moment = (fraction.sum(axis=0) * np.exp(1j * angle)).sum()
    return (np.angle(moment) / (2.0 * np.pi) * nx) % nx


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_translating_runs_the_intended_case(translating_run):
    assert translating_run.parameter("translate") == SPEED
    assert translating_run.config["surface_tension"] == "stress"
    assert translating_run.config["source_stencil"] == translating_run.stencil


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_translating_keeps_momentum(translating_run):
    """The capillary stress is a divergence: the box keeps its momentum, to
    the six digits the fields are written with."""
    rho = translating_run.density(LAST)
    u = translating_run.velocity(LAST)
    mean_velocity = (rho * u[..., 0]).sum() / rho.sum()
    assert mean_velocity == pytest.approx(SPEED, rel=1e-5)
    assert np.abs((rho * u[..., 1]).sum()) / rho.sum() < 1e-7


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_translating_droplet_speed(translating_run):
    """The droplet keeps up with the flow with the matched stencil, and falls
    8 % behind with the nine-point one, as measured."""
    nx = translating_run.shape[1]
    travelled = centroid_x(translating_run, LAST) - centroid_x(translating_run, FIRST)
    # the droplet moves 20 nodes of the 128: unwrap by the expected distance
    expected = SPEED * (LAST - FIRST)
    travelled -= nx * np.round((travelled - expected) / nx)
    speed = travelled / expected
    assert speed == pytest.approx(MEASURED_SPEED[translating_run.stencil], abs=SPEED_TOLERANCE)
    if translating_run.stencil == "matched":
        assert speed == pytest.approx(1.0, abs=SPEED_TOLERANCE)
