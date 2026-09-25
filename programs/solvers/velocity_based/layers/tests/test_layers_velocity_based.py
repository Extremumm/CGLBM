"""The velocity-based scheme on two sheared layers 1e4 apart in density.

One run of programs/solvers/velocity_based/layers, `layers E8 1e4 0.01`: a
Kolmogorov flow u = U sin(k y) across a heavy layer (0 < y < Ly/2) and a light
one, with equal dynamic viscosities, started in its steady state. The
interfaces sit where the shear stress is largest. The heavy layer is held there
by its inertia; the light one relaxes within a few thousand steps to whatever
the interface lets through, and a slip there shows as an offset of the whole
light layer.

It pins the accuracy of the viscous stress across the interface, which the
lighter-density weighting of the link exchange halves and the link viscosity
restores; see src/lbm/velocity_based.h and docs/numerics.md.
"""

import numpy as np
import pytest
from pycglbm.testing import artifacts_dir, run_program

NUM_STEPS = 10000
AMPLITUDE = 0.01

# Measured on this run, pinned to catch regressions: largest error of the
# velocity in each layer, as a fraction of AMPLITUDE.
MEASURED_LIGHT_ERROR = 0.0193
MEASURED_HEAVY_ERROR = 0.0025
#: The output is written with ten significant digits.
MOMENTUM_TOLERANCE = 1.0e-9


@pytest.fixture(scope="module")
def layers_run():
    return run_program(
        "layers", artifacts_dir() / "layers", args=("E8", "1e4", str(AMPLITUDE)), timeout=600
    )


def _interfaces(phase):
    """Where the phase field crosses zero, going up the column, periodically."""
    crossings = []
    n = len(phase)
    for j in range(n):
        below, above = phase[j], phase[(j + 1) % n]
        if (below > 0.0) != (above > 0.0):
            crossings.append((j + below / (below - above)) % n)
    return sorted(c if c < n - 0.5 else c - n for c in crossings)


def _errors(case, timestep):
    """Largest |u - U sin(k y)| / U in the light and the heavy layer."""
    velocity = case.velocity(timestep)[..., 0]
    ny = velocity.shape[0]
    exact = AMPLITUDE * np.sin(2.0 * np.pi * np.arange(ny) / ny)[:, None]
    error = np.abs(velocity - exact) / AMPLITUDE
    light = case.phase(timestep) < 0.0
    return error[light].max(), error[~light].max()


@pytest.mark.long
@pytest.mark.verification
def test_verification_layers_velocity_based_stays_a_layered_flow(layers_run):
    """The run completes, the flow stays uniform along x and across y, and the
    interfaces stay where they were."""
    t = layers_run.last_timestep
    assert t == NUM_STEPS
    velocity = layers_run.velocity(t)
    assert np.isfinite(velocity).all()
    assert np.abs(velocity[..., 0] - velocity[:, :1, 0]).max() < 1e-12
    # only the artificial compressibility moves anything across the layers
    assert np.abs(velocity[..., 1]).max() < 1e-5 * AMPLITUDE
    assert _interfaces(layers_run.phase(t)[:, 0]) == pytest.approx([0.0, 64.0], abs=1e-3)


@pytest.mark.long
@pytest.mark.verification
def test_verification_layers_velocity_based_momentum_is_conserved(layers_run):
    """The body force sums to zero over the period, and so does every link
    exchange: the total momentum does not change."""
    start = float((layers_run.density(0) * layers_run.velocity(0)[..., 0]).sum())
    for t in layers_run.timesteps[1:]:
        now = float((layers_run.density(t) * layers_run.velocity(t)[..., 0]).sum())
        assert now == pytest.approx(start, rel=MOMENTUM_TOLERANCE)


@pytest.mark.long
@pytest.mark.validation
def test_validation_layers_velocity_based_stress_crosses_the_interface(layers_run):
    """Both layers stay within 2 % of the steady flow.

    The light layer drifts by a nearly uniform 1.9 % of U within the first few
    thousand steps and stays there: a slip across the diffuse interface. The
    link viscosity decides it: adding 5 % of the lattice's own viscosity to
    every link turns it into 2.2 % the other way, adding 25 % into 16 %.
    """
    light, heavy = _errors(layers_run, NUM_STEPS)
    assert light == pytest.approx(MEASURED_LIGHT_ERROR, abs=0.002)
    assert heavy == pytest.approx(MEASURED_HEAVY_ERROR, abs=0.0005)
    assert light < 0.025
    # steady: the last half of the run adds little
    light_half, _ = _errors(layers_run, NUM_STEPS // 2)
    assert abs(light - light_half) < 0.002
