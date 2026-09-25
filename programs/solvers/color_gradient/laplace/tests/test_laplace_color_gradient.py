"""Laplace-law benchmark for the colour-gradient solver.

A static droplet of radius R should carry a pressure jump dp = sigma / R across
its interface. The case initialises the pressure so the jump holds exactly at
t = 0 (p1_inf carries a -sigma/radius offset), which the verification test
checks; the validation tests then look at the state the solver relaxes to.

Which radius to score against is not a formality here. The droplet is
compressible, so it does not stay at the prescribed radius, and the phase field
and the density field sit at *different* radii: phi is a mass fraction, and its
zero contour lies (W/2) ln(rho_1/rho_2) outside the density interface -- see
docs/numerics.md and src/lbm/mixture.h. The surface tension acts on the
density interface, so the jump is scored against the density radius.

This is the shipped case, `laplace` with no arguments: density ratio 20, one
kinematic viscosity for both fluids. test_laplace_high_density_ratio.py runs
the same program at a density ratio of 1e4.
"""

import math

import pytest
from pycglbm.testing import artifacts_dir, run_program

# The case parameters are not repeated here: the solver reports them as
# `key = value` lines at the top of run.log, and CaseOutput.parameter reads them
# back. A constant duplicated in a test is a constant that drifts.

# Measured baselines for the current scheme (restored equation of state, E8
# colour gradient, interface started in mechanical equilibrium, located and
# driven through the bulk-normalised phase field, tension applied as a
# continuum surface force). They are what the solver does, pinned to catch
# regressions. The history behind them, on this same case at density ratio 20:
#   linear-mixing pressure, E4 gradient      : 0.720 * sigma/R_rho   (28 % low)
#   restored EOS, E4 gradient                : 0.962 * sigma/R_rho
#   restored EOS, E8 gradient                : 0.962 * sigma/R_rho, fewer currents
#   the above, started in mechanical equilib.: 0.895 * sigma/R_rho
#   the above, Ba et al. phi_N + CSF tension : 1.021 * sigma/R_rho
# The last step is what carries the density ratio to about 500 within a couple
# of per cent rather than a factor of two, and it cut the spurious currents by
# 72x on this case. (500, not 1000: at 1000 the case diverges after 8.7e4
# steps, which shorter runs had hidden. See docs/numerics.md.)
# `--initial-profile=colour --interface-field=colour
# --surface-tension=perturbation` reproduces the row above it.
MEASURED_JUMP_RATIO = 1.021
MEASURED_JUMP_TOLERANCE = 0.03

#: How far the phase interface may sit from (W/2) ln(rho_1/rho_2) outside the
#: density interface, in lattice units. The split used to be pinned as a
#: measured artefact (2.27 linear, 2.28 equilibrium, 2.35 with the profile laid
#: down in phi_N); it is the mass-fraction nature of phi, 2.40 at this density
#: ratio with W = 1.6. The density radius is the physical one:
#: rho = (rho1 + rho2) / 2 is exactly where the two components occupy equal
#: volume, which is the surface the tension acts on. See docs/numerics.md.
INTERFACE_SPLIT_TOLERANCE = 0.1

#: Spurious currents at steady state, in lattice units.
#: The continuum-surface-force tension of Ba et al. dropped this from 1.19e-3,
#: a factor of 72, which is the clearest single sign that the tension is now
#: applied on the interface rather than beside it.
MEASURED_MAX_VELOCITY = 1.66e-5

#: Radius the density field starts at, in lattice units.
#: The case asks for 10 and now gets it. It used to start at 8.41, because the
#: tanh was laid down in the colour field, whose zero contour is not the
#: interface once the two bulk densities differ -- Ba et al. Eq. (21).
MEASURED_INITIAL_DENSITY_RADIUS = 10.06


@pytest.fixture(scope="module")
def laplace_run():
    """One shared run of the Laplace case, reused by every test in this module."""
    return run_program("laplace", artifacts_dir() / "laplace", timeout=1800)


@pytest.fixture(scope="module")
def case(laplace_run):
    """The parameters the run reported, and the quantities derived from them."""
    sigma = laplace_run.parameter("sigma")
    radius = laplace_run.parameter("radius")
    return {
        "sigma": sigma,
        "radius": radius,
        "steps": laplace_run.parameter("steps", int),
        # Analytic Laplace jump for the prescribed radius.
        "jump": sigma / radius,
        "density_ratio": laplace_run.parameter("rho1") / laplace_run.parameter("rho2"),
        "width": laplace_run.parameter("ch_width_ope"),
        "width_init": laplace_run.parameter("ch_width_init"),
        # Radii, from the domain centre, delimiting bulk and far field.
        "inner": 0.5 * radius,
        "outer": 3.0 * radius,
    }


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_color_gradient_initial_pressure_jump(laplace_run, case):
    """At t = 0 the prescribed field must satisfy dp = sigma / R exactly."""
    jump = laplace_run.pressure_jump(0, inner=case["inner"], outer=case["outer"])
    assert jump == pytest.approx(case["jump"], rel=1.0e-3)
    # The *density* field starts on the prescribed radius, because the profile
    # is laid down in phi_N and rho = (rho1 + rho2) / 2 is exactly phi_N = 0.
    # This is the point of `--initial-profile=normalised`: `sigma / radius` in
    # `p1_inf` and the droplet it is applied to now mean the same radius.
    assert laplace_run.density_interface_radius(0) == pytest.approx(
        MEASURED_INITIAL_DENSITY_RADIUS, abs=0.1
    )
    # The colour field's zero contour does not, and cannot: phi is a mass
    # fraction, and a tanh profile of width W in phi_N is the same profile in
    # phi shifted (W/2) ln(rho1/rho2) outward.
    offset = 0.5 * case["width_init"] * math.log(case["density_ratio"])
    assert laplace_run.phase_interface_radius(0) == pytest.approx(
        laplace_run.density_interface_radius(0) + offset, abs=0.1
    )


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_color_gradient_reaches_a_steady_state(laplace_run, case):
    """The tail of the run must be stationary, not still relaxing."""
    tail = [t for t in laplace_run.timesteps if t >= case["steps"] // 3]
    jumps = [laplace_run.pressure_jump(t, inner=case["inner"], outer=case["outer"]) for t in tail]
    radii = [laplace_run.density_interface_radius(t) for t in tail]

    assert max(jumps) - min(jumps) < 1.0e-4 * case["jump"]
    assert max(radii) - min(radii) < 0.01


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_pressure_jump(laplace_run, case):
    """The relaxed jump, scored against the radius the density field settles at."""
    radius = laplace_run.density_interface_radius(case["steps"])
    jump = laplace_run.pressure_jump(case["steps"], inner=case["inner"], outer=case["outer"])

    ratio = jump / (case["sigma"] / radius)
    assert ratio == pytest.approx(MEASURED_JUMP_RATIO, abs=MEASURED_JUMP_TOLERANCE)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_interface_split_is_the_mass_fraction_offset(
    laplace_run, case
):
    """The phase and density interfaces sit (W/2) ln(rho_1/rho_2) apart.

    phi is a mass fraction and the density follows the volume fraction; a tanh
    profile of width W in one is the same profile in the other, shifted by that
    much (src/lbm/mixture.h). The recolouring maintains the width W, so the
    relaxed split is predicted, not merely pinned. Under the old linear-mixing
    pressure the two coincided, because that pressure made phi a volume
    fraction.
    """
    phase_radius = laplace_run.phase_interface_radius(case["steps"])
    density_radius = laplace_run.density_interface_radius(case["steps"])

    split = phase_radius - density_radius
    expected = 0.5 * case["width"] * math.log(case["density_ratio"])
    assert split == pytest.approx(expected, abs=INTERFACE_SPLIT_TOLERANCE)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_spurious_currents(laplace_run, case):
    """Parasitic currents at the interface must stay at their measured level."""
    velocity = laplace_run.velocity(case["steps"])
    speed = (velocity[..., 0] ** 2 + velocity[..., 1] ** 2) ** 0.5
    assert speed.max() == pytest.approx(MEASURED_MAX_VELOCITY, rel=0.15)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_color_gradient_higher_isotropy_helps(laplace_run, case):
    """E8 must not have more parasitic current than the E4 stencil.

    Leclaire et al. (Computers & Fluids 48, 98, 2011) report the isotropic
    gradient cutting spurious currents; this checks the direction of that claim
    on the real case rather than on the stencil alone.
    """
    coarse = run_program("laplace", artifacts_dir() / "laplace_e4", args=("E4",), timeout=1800)

    def peak(run):
        v = run.velocity(case["steps"])
        return float(((v[..., 0] ** 2 + v[..., 1] ** 2) ** 0.5).max())

    assert peak(laplace_run) <= peak(coarse)
