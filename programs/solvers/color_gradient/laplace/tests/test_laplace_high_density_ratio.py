"""Laplace-law benchmark at a density ratio of 10^4.

The same program as test_laplace_color_gradient.py, with a droplet 10^4 times
denser than its surroundings and the same dynamic viscosity in both fluids:

    laplace --rho1=1e4 --nu=<nu/1e4> --nu-b=<nu/1e4> --nu2=<nu> --nu-b2=<nu>
            --viscosity-mixing=dynamic --surface-tension=stress

Everything else is the shipped configuration: the equilibrium start and the
colour gradient on the bulk-normalised phase field. The two options are what
carry it from about 500 to 10^4 (docs/numerics.md, "Beyond 500"): the
viscosities mixed as rho nu, which keeps tau uniform across the interface, and
the tension as the divergence of the capillary stress, which conserves momentum.
Interpolating nu instead leaves tau ~ 1e4 on the interface nodes, and the jump
and the currents drift off within 10^4 steps.

What is checked:

- the initial state: the density interface on the prescribed radius, the
  Laplace jump already in place, and the phase field phi = 0 contour displaced
  by (W/2) ln(rho_1/rho_2) because phi is a mass fraction;
- the run stays bounded and settles;
- the relaxed jump obeys Laplace's law against the radius the density settles at;
- the spurious currents stay at their measured level.

30 000 steps is 5.5e3 relaxation times at tau = 5.5. The same case run to 1.2e5
steps does not diverge but is still converging towards the 1e3 state (the jump
reaches 1.025, the currents fall to 4.5e-5), which docs/numerics.md records;
this test pins the state at 3e4 steps.

Only static droplets are claimed at this ratio. A droplet translating faster
than about 1e-4 lattice units per step still breaks the colour-gradient scheme
at 10^4; the velocity-based droplet solver is the one for that.
"""

import math

import numpy as np
import pytest
from pycglbm.testing import artifacts_dir, run_program

C_DX = 1.0e-5
C_DT = C_DX / 347.0 / math.sqrt(3.0)
#: The shipped kinematic viscosity, in lattice units, given to the light fluid.
NU_LIGHT = 1.0e-2 / (C_DX**2 / C_DT)

DENSITY_RATIO = 1.0e4
#: Same dynamic viscosity in both fluids: nu_heavy = nu_light / ratio.
NU_HEAVY = NU_LIGHT / DENSITY_RATIO
ARGS = (
    "--rho1=1e4",
    f"--nu={NU_HEAVY!r}",
    f"--nu-b={NU_HEAVY!r}",
    f"--nu2={NU_LIGHT!r}",
    f"--nu-b2={NU_LIGHT!r}",
    "--viscosity-mixing=dynamic",
    "--surface-tension=stress",
)

# Measured on this case, pinned to catch regressions.
#: Delta p / (sigma / R_rho) at the end of the run.
MEASURED_JUMP_RATIO = 1.018
MEASURED_JUMP_TOLERANCE = 0.02
#: Largest spurious velocity at the end of the run, lattice units.
MEASURED_MAX_VELOCITY = 6.5e-4


@pytest.fixture(scope="module")
def high_ratio_run():
    """One shared run of the Laplace case at density ratio 1e4."""
    return run_program("laplace", artifacts_dir() / "laplace_ratio_1e4", args=ARGS, timeout=2400)


@pytest.fixture(scope="module")
def case(high_ratio_run):
    """The parameters the run reported, and the quantities derived from them."""
    sigma = high_ratio_run.parameter("sigma")
    radius = high_ratio_run.parameter("radius")
    return {
        "sigma": sigma,
        "radius": radius,
        "steps": high_ratio_run.parameter("steps", int),
        "width_init": high_ratio_run.parameter("ch_width_init"),
        # Analytic Laplace jump for the prescribed radius.
        "jump": sigma / radius,
        # Radii, from the domain centre, delimiting bulk and far field.
        "inner": 0.5 * radius,
        "outer": 3.0 * radius,
    }


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_density_ratio_runs_the_intended_case(high_ratio_run):
    """The options reached the solver: a flag lost on the way is a different case."""
    assert high_ratio_run.parameter("rho1") / high_ratio_run.parameter("rho2") == DENSITY_RATIO
    assert high_ratio_run.config["viscosity_mixing"] == "dynamic"
    assert high_ratio_run.config["surface_tension"] == "stress"
    assert high_ratio_run.config["interface_field"] == "normalised"


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_density_ratio_initial_state(high_ratio_run, case):
    """The jump holds at t = 0 and the interfaces start where the theory puts them."""
    jump = high_ratio_run.pressure_jump(0, inner=case["inner"], outer=case["outer"])
    assert jump == pytest.approx(case["jump"], rel=1.0e-3)

    density = high_ratio_run.density(0)
    assert density.max() == pytest.approx(DENSITY_RATIO, rel=1.0e-5)
    density_radius = high_ratio_run.density_interface_radius(0)
    assert density_radius == pytest.approx(case["radius"], abs=0.1)

    # phi is a mass fraction: its zero sits (W/2) ln(rho_1/rho_2) further out
    offset = 0.5 * case["width_init"] * math.log(DENSITY_RATIO)
    assert high_ratio_run.phase_interface_radius(0) == pytest.approx(
        density_radius + offset, abs=0.1
    )


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_density_ratio_stays_bounded(high_ratio_run, case):
    """No NaN, the phase field within [-1, 1], positive pressure everywhere."""
    t = high_ratio_run.last_timestep
    assert t == case["steps"]
    for field in high_ratio_run.fields(t).values():
        assert np.isfinite(field).all()
    phase = high_ratio_run.phase(t)
    assert phase.max() <= 1.0 + 1.0e-6
    assert phase.min() >= -1.0 - 1.0e-6
    assert high_ratio_run.pressure(t).min() > 0.0


@pytest.mark.long
@pytest.mark.verification
def test_verification_laplace_high_density_ratio_settles(high_ratio_run, case):
    """Over the last two thirds the jump and the radius barely move."""
    tail = [t for t in high_ratio_run.timesteps if t >= case["steps"] // 3]
    jumps = [
        high_ratio_run.pressure_jump(t, inner=case["inner"], outer=case["outer"]) for t in tail
    ]
    radii = [high_ratio_run.density_interface_radius(t) for t in tail]

    assert max(jumps) - min(jumps) < 1.0e-2 * case["jump"]
    assert max(radii) - min(radii) < 0.06


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_high_density_ratio_pressure_jump(high_ratio_run, case):
    """Laplace's law against the radius the density field settles at.

    The remaining couple of percent is the finite interface width at R = 10,
    not the density ratio: the shipped case gives 1.021 at a density ratio of
    20, and the force alone integrates to 1.007 sigma/R at R = 20 (unit test of
    lbm/mixture).
    """
    radius = high_ratio_run.density_interface_radius(case["steps"])
    jump = high_ratio_run.pressure_jump(case["steps"], inner=case["inner"], outer=case["outer"])

    ratio = jump / (case["sigma"] / radius)
    assert ratio == pytest.approx(MEASURED_JUMP_RATIO, abs=MEASURED_JUMP_TOLERANCE)
    assert ratio == pytest.approx(1.0, abs=0.05)


@pytest.mark.long
@pytest.mark.validation
def test_validation_laplace_high_density_ratio_spurious_currents(high_ratio_run, case):
    """Parasitic currents stay at their measured level, and well below c_s."""
    velocity = high_ratio_run.velocity(case["steps"])
    speed = np.hypot(velocity[..., 0], velocity[..., 1])
    assert speed.max() < 3.0 * MEASURED_MAX_VELOCITY
