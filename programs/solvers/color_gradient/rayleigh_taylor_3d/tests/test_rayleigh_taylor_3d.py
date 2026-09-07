"""Rayleigh-Taylor in three dimensions: the wall and the body force, together.

``laplace_3d`` and ``oscillation_3d`` are both triply periodic and neither has a
body force, so this case is what exercises the two remaining paths through the
three-dimensional solver. It is also the only three-dimensional case with no
closed form to be scored against, and the tests here are written to that:

**What is asserted.** That the case is set up as the unstable arrangement and
not its mirror image -- a sign slip in the initial layer gives a perfectly
stable interface that would sit there looking plausible for the whole run; that
the spike falls and the bubble rises; that the mirror symmetry the single mode
starts with survives; and that the phase field stays inside [-1, 1] throughout,
which is the statement that the segregation operator still holds where the
interface is being stretched rather than sitting still.

**What is not.** The inviscid growth rate ``sqrt(A g k)``. At this resolution
the viscous correction ``-nu k^2`` is a quarter of it and the interface is 1.6
nodes wide against a wavelength of 32, so agreement or disagreement with it
would say more about the resolution than about the scheme. The quantitative
statements live in ``oscillation_3d``, and the invariants that hold whatever the
flow does -- mass across a wall, hydrostatic balance, x against z -- are checked
exactly in ``programs/unit_testing/lbm/two_population_3d``.
"""

import numpy as np
import pytest
from pycglbm.testing import artifacts_dir, run_program


def interface_height(phase: np.ndarray) -> np.ndarray:
    """Height of the ``phi_N = 0`` contour in each column of a ``[y, x]`` slice.

    The dense component is on top, so each column runs from -1 at the bottom to
    +1 at the top and crosses once. The crossing is interpolated linearly
    between the two nodes that bracket it, and a column that somehow never
    crosses gives ``nan`` rather than a plausible number.
    """
    ny, nx = phase.shape
    heights = np.full(nx, np.nan)
    for column in range(nx):
        profile = phase[:, column]
        below = np.nonzero((profile[:-1] < 0.0) & (profile[1:] >= 0.0))[0]
        if below.size:
            j = below[0]
            heights[column] = j + (-profile[j]) / (profile[j + 1] - profile[j])
    return heights


@pytest.fixture(scope="module")
def run_3d():
    """One shared run of the case. It takes about eight minutes on 24 cores."""
    return run_program("rayleigh_taylor_3d", artifacts_dir() / "rayleigh_taylor_3d", timeout=2700)


@pytest.fixture(scope="module")
def heights(run_3d):
    """The interface height along the z = nz/2 slice, at every output step."""
    return {t: interface_height(run_3d.phase(t)) for t in run_3d.timesteps}


@pytest.mark.long
@pytest.mark.verification
def test_verification_rayleigh_taylor_3d_is_the_unstable_arrangement(run_3d):
    """Dense on top, walls along y, gravity on, no surface tension.

    Every one of these is a way the case could quietly become a different one.
    The arrangement especially: `cosine_layer_3d(A, inverted=true)` is the same
    interface with the dense component underneath, which is stable, and a run of
    it would look like a case that simply has not gone unstable yet.
    """
    assert run_3d.parameter("nz", int) == run_3d.parameter("nx", int) > 1
    assert run_3d.parameter("boundary", str) == "wall_y"
    assert run_3d.parameter("gravity") > 0.0
    assert run_3d.parameter("sigma") == 0.0
    assert run_3d.parameter("rho1") > run_3d.parameter("rho2")

    phase = run_3d.phase(0)
    assert phase[-1, :].min() > 0.9  # component 1, the dense one, at the top
    assert phase[0, :].max() < -0.9  # component 2 at the bottom


@pytest.mark.long
@pytest.mark.verification
def test_verification_rayleigh_taylor_3d_starts_on_the_prescribed_mode(run_3d, heights):
    """The initial interface is `y0 - A nx cos(2 pi x / nx)` along this slice.

    The slice is at z = nz/2, where `cos(2 pi z / nz)` is -1, so the mode shows
    up along it inverted: the spike is at x = 0 and the bubble at x = nx/2.
    """
    nx = run_3d.parameter("nx", int)
    ny = run_3d.parameter("ny", int)
    start = heights[0]
    assert np.all(np.isfinite(start))

    amplitude = 0.5 * (start.max() - start.min())
    x = np.arange(nx)
    expected = ny // 2 - amplitude * np.cos(2.0 * np.pi * x / nx)
    assert np.max(np.abs(start - expected)) < 0.2
    assert start[0] == pytest.approx(start.min(), abs=0.05)
    assert start[nx // 2] == pytest.approx(start.max(), abs=0.05)


@pytest.mark.long
@pytest.mark.verification
def test_verification_rayleigh_taylor_3d_keeps_its_mirror_symmetry(run_3d, heights):
    """`x -> nx - x` leaves the initial mode alone, so it must leave the flow alone.

    Not to round-off -- the run is long and the flow is unstable, so any
    difference is amplified along with everything else -- but the interface must
    stay symmetric to well within a node for the whole run.
    """
    nx = run_3d.parameter("nx", int)
    mirror = np.concatenate([[0], np.arange(nx - 1, 0, -1)])
    for timestep, height in heights.items():
        assert np.nanmax(np.abs(height - height[mirror])) < 0.05, f"at step {timestep}"


@pytest.mark.long
@pytest.mark.verification
def test_verification_rayleigh_taylor_3d_stays_inside_the_phase_range(run_3d):
    """phi_N must not leave [-1, 1] while the interface is being stretched.

    The bound comes from both distributions staying non-negative, which the unit
    tests check on a droplet sitting still. This is the same statement where the
    interface is being pulled apart, which is where it would fail first.
    """
    for timestep in run_3d.timesteps:
        phase = run_3d.phase(timestep)
        assert np.all(np.isfinite(phase)), f"at step {timestep}"
        assert np.abs(phase).max() <= 1.0 + 1.0e-9, f"at step {timestep}"


@pytest.mark.long
@pytest.mark.validation
def test_validation_rayleigh_taylor_3d_spike_falls_and_bubble_rises(run_3d, heights):
    """The instability has to run: the two fronts must separate, monotonically.

    This is the whole qualitative content of the case, and it is what a wall
    that leaked, a force with the wrong sign, or an interface pinned to the
    lattice would each break in a different visible way.
    """
    timesteps = sorted(heights)
    spikes = np.array([np.nanmin(heights[t]) for t in timesteps])
    bubbles = np.array([np.nanmax(heights[t]) for t in timesteps])

    assert np.all(np.diff(spikes) < 0.0)
    assert np.all(np.diff(bubbles) > 0.0)
    assert bubbles[-1] - spikes[-1] > 4.0 * (bubbles[0] - spikes[0])
