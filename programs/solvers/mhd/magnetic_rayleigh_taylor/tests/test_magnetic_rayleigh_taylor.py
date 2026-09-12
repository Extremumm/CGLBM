"""Magnetic Rayleigh-Taylor: the growth rate, and what the field does to it.

Three runs of one mode at three conductivities, which at a fixed field is three
field strengths in disguise. The measurement is the growth rate of the seeded
cosine, and it is scored against the finite-depth quasi-static dispersion
relation used by the ``mrt`` campaign in ``basiliskMHD/QSMHD`` -- the same
relation, restated in ``dispersion.py``, evaluated at this case's lattice
parameters.

**What is asserted is the ratio, not the rate.** The relation is inviscid, and
this flow is not: ``nu k^2 / omega`` is 6.9 % at the light fluid's viscosity,
and lowering the viscosity far enough to make an absolute comparison sharp puts
the relaxation time on top of 1/2. So the absolute rate is measured and
reported, and what is *scored* is ``omega(B) / omega(0)``, from which the
viscous correction very largely cancels because it is the same correction in
both. The relation predicts that ratio to fall to 0.737 and then 0.342 as the
interaction parameter ``N = sigma B^2 / (rho omega)`` goes from 0 to 2 to 10,
and the measurement gives 0.776 and 0.374 -- errors of 5.3 % and 9.2 % on a
quantity that changes by a factor of two across the sweep, which is far outside
anything the viscosity ambiguity could account for.

**Why this case and not the physical one.** The ``mrt`` campaign is a liquid
metal over a molten salt: 6.8e-7 m^2/s of kinematic viscosity across 0.1 m.
Matching that at 128 nodes puts the lattice viscosity at 1.1e-5 and the
relaxation time at 0.50003, which no lattice-Boltzmann collision survives. What
is carried over instead is the configuration and the two ratios that set the
character of the problem -- density 3.6501 and conductivity 10825 -- with the
rest chosen so the scheme can run it.

**What this case caught.** With the library's default harmonic conductivity
blend, the measured suppression at N = 10 was 0.567 against a predicted 0.342 --
66 % too little braking -- and the absolute growth rate came out *above* the
inviscid relation, which is not physically available. The cause is the blend: a
face half in each phase gets a harmonic conductivity 2700 times smaller than the
arithmetic one at this contrast, so two or three nodes of interface stop
carrying current exactly where the shear is. The current in this flow runs
*along* the interface and the two phases are therefore in parallel, so
arithmetic is the correct average, and the case sets it. The same measurement
with arithmetic gives 0.374 against 0.342. The failure is reproducible with
``--conductivity=harmonic``.
"""

import numpy as np
import pytest
from pycglbm.testing import artifacts_dir, run_program

from dispersion import fit_growth_rate, growth_rate, mode_amplitude

#: Conductivity of the heavy fluid, and the interaction parameter it gives.
#:
#: ``N = sigma B^2 / (rho1 omega_0)`` at the zero-field growth rate, with
#: ``B = 1``. Zero is the hydrodynamic reference the other two are measured
#: against.
FIELDS = {"n0": 0.0, "n2": 7.62e-3, "n10": 3.81e-2}

#: Upper end of the fit window, in nodes. A sixteenth of the wavelength, so the
#: mode is still linear; past this the spike and the bubble stop being mirror
#: images and the single cosine stops describing the interface.
FIT_UPPER = 8.0

#: Tolerance on the suppression ratio.
#:
#: Measured +5.3 % at N = 2 and +9.2 % at N = 10 with the arithmetic blend the
#: case sets. What is
#: left is the diffuse interface -- the conductivity, the density and the
#: velocity all vary over two or three nodes where the relation has them jump --
#: and the residue of viscosity that the ratio does not cancel.
RATIO_TOLERANCE = 0.15


def theoretical_rate(run, conductivity: float) -> float:
    """The relation at this run's own parameters, read back from its log."""
    return growth_rate(
        rho1=run.parameter("rho1"),
        rho2=run.parameter("rho2"),
        gravity=run.parameter("gravity"),
        gamma=run.parameter("sigma"),
        k=run.parameter("wavenumber"),
        depth=run.parameter("layer_depth"),
        field=run.parameter("b_y"),
        conductivity1=conductivity,
        conductivity2=conductivity / run.parameter("conductivity_ratio")
        if conductivity > 0.0
        else 0.0,
    )


def measured_rate(run) -> tuple[float, int]:
    times = run.timesteps
    amplitudes = [mode_amplitude(run.phase(t)) for t in times]
    return fit_growth_rate(times, amplitudes, run.parameter("seed_amplitude"), FIT_UPPER)


@pytest.fixture(scope="module", params=sorted(FIELDS))
def run_mode(request):
    """One run per field strength. Each takes about an hour on 24 cores."""
    name = request.param
    conductivity = FIELDS[name]
    return name, conductivity, run_program(
        "magnetic_rayleigh_taylor",
        artifacts_dir() / f"magnetic_rayleigh_taylor_{name}",
        args=(f"--sigma-e1={conductivity}",),
        timeout=7200,
    )


@pytest.fixture(scope="module")
def rates():
    """Every field strength's measured and predicted rate, in one pass.

    Separate from ``run_mode`` because the suppression ratio needs two runs at
    once and a per-run fixture cannot give it that.
    """
    measured, predicted = {}, {}
    for name, conductivity in FIELDS.items():
        run = run_program(
            "magnetic_rayleigh_taylor",
            artifacts_dir() / f"magnetic_rayleigh_taylor_{name}",
            args=(f"--sigma-e1={conductivity}",),
            timeout=7200,
        )
        rate, count = measured_rate(run)
        assert count >= 4, f"{name}: only {count} samples inside the fit window"
        measured[name] = rate
        predicted[name] = theoretical_rate(run, conductivity)
    return measured, predicted


@pytest.mark.long
@pytest.mark.validation
def test_validation_magnetic_rayleigh_taylor_is_the_case_it_claims_to_be(run_mode):
    """The arrangement, the field, and the conductivity contrast.

    The dense fluid on top and the field normal to the interface are what make
    this the configuration the dispersion relation was written for; the
    conductivity ratio is what makes the potential problem a hard one.
    """
    _, _, run = run_mode
    assert run.parameter("mhd", str) == "true"
    assert run.parameter("boundary", str) == "wall_y"
    assert run.parameter("lattice_3d", str) == "D3Q27"
    assert run.parameter("rho1") > run.parameter("rho2")
    assert run.parameter("gravity") > 0.0
    assert run.parameter("b_y") > 0.0
    assert run.parameter("b_x") == run.parameter("b_z") == 0.0
    assert run.parameter("conductivity_ratio") == pytest.approx(10825.0, rel=1e-6)

    phase = run.phase(0)
    assert phase[-1, :].min() > 0.9  # the dense component at the top
    assert phase[0, :].max() < -0.9


@pytest.mark.long
@pytest.mark.validation
def test_validation_magnetic_rayleigh_taylor_potential_is_uniform(run_mode):
    """``div(sigma u x B) = 0`` for a two-dimensional mode and a normal field.

    The motional field points along z and nothing varies along z, so the
    potential is a constant -- which is a property of the configuration, and the
    reason it has a closed form at all. Asserted because a drive assembled wrong
    would put structure there, and because every number below is read from a
    current the potential contributes nothing to.

    The potential solve itself is covered by
    ``programs/unit_testing/lbm/mhd_potential``.
    """
    _, _, run = run_mode
    log = (run.rundir / "run.log").read_text()
    assert "charge imbalance 0\n" in log or "charge imbalance 0 " in log


@pytest.mark.long
@pytest.mark.validation
def test_validation_magnetic_rayleigh_taylor_grows(rates):
    """Every field strength still goes unstable, and the field slows it.

    Monotone in the field, which is the qualitative statement the whole case is
    for: a magnetic field normal to a Rayleigh-Taylor interface does not stop
    the instability, it brakes it.
    """
    measured, _ = rates
    assert measured["n0"] > measured["n2"] > measured["n10"] > 0.0


@pytest.mark.long
@pytest.mark.validation
def test_validation_magnetic_rayleigh_taylor_suppression_matches_theory(rates):
    """``omega(B)/omega(0)`` against the relation. Predicted 0.737 and 0.342.

    This is the scored measurement. The viscous correction the relation omits is
    the same in numerator and denominator and largely cancels; what is left is
    the magnetic braking, which is a factor of two between the weakest and the
    strongest field and cannot be confused with it. Measured 0.374 against 0.342
    at N = 10.
    """
    measured, predicted = rates
    for name in ("n2", "n10"):
        assert measured[name] / measured["n0"] == pytest.approx(
            predicted[name] / predicted["n0"], rel=RATIO_TOLERANCE
        )


@pytest.mark.long
@pytest.mark.validation
def test_validation_magnetic_rayleigh_taylor_absolute_rate_is_low_by_viscosity(rates):
    """Bounded on the side the physics says it must fall.

    The relation is inviscid, so the simulation must grow *slower*, never
    faster. The deficit is roughly ``2 nu k^2`` in absolute terms -- 1.4e-4 here
    -- which is 14 % of the zero-field rate and a larger fraction of the damped
    ones; measured 0.856, 0.901 and 0.934 of the relation at N = 0, 2 and 10.

    A rate *above* the inviscid one is the failure this catches, and it is not
    hypothetical: the harmonic conductivity blend produced exactly that at
    N = 10, at 1.42 times the relation, before the case was changed to the
    average its own geometry calls for.
    """
    measured, predicted = rates
    for name in FIELDS:
        assert 0.5 < measured[name] / predicted[name] <= 1.0
