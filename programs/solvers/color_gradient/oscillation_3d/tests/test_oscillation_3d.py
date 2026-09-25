"""Free oscillation of a deformed droplet, against Lamb's mode-2 frequency.

``laplace_3d`` is the three-dimensional scheme with nothing moving. This is the
same scheme in motion, and it is scored against a formula that only exists in
three dimensions:

    omega_2^2 = 24 sigma / { R^3 (3 rho_in + 2 rho_out) }

The two-dimensional analogue, ``6 sigma / [R^3 (rho_in + rho_out)]``, differs in
the numerator, in how the two densities are weighted and in the mode shape
behind both, so no two-dimensional run could have been scored against this. It
is the dynamic counterpart of ``2 sigma / R`` against ``sigma / R``.

**What is asserted, and what is only recorded.** The frequency comes out 18 %
below Lamb's at the shipped resolution, and that gap is the finding rather than
a failure: it closes as `1/R` as the droplet is better resolved -- 0.776 at
R = 8, 0.822 at R = 10, 0.855 at R = 13, with `gap * R` flat at 1.8 -- while the
static surface tension at these exact parameters is right to 2.4 % and removing
the density ratio makes the agreement worse rather than better.

The viscosity is a genuine part of the gap here, unlike in the version of this
docstring that preceded the enhanced-equilibrium fix: the case is pinned at
mu = 0.06 by the D3Q19 positivity bound (see the case header), the damping ratio
is 0.35, and doubling mu costs 10.3 points. Roughly a fifth of the trend across
the three radii is that damping falling with R; the rest is resolution.

So the tests below assert the measured value with a tolerance that would catch a
change in the scheme, and assert separately the things that must be exactly true
whatever the resolution: that the two equatorial semi-axes stay equal, that the
droplet keeps its volume, and that what is being fitted really is a single
decaying mode.
"""

import numpy as np
import pytest
from pycglbm.oscillation import fit_damped_oscillation, lamb_frequency
from pycglbm.testing import artifacts_dir, run_program

#: Steps skipped before the fit starts.
#:
#: The droplet is released from rest with no pressure jump inside it, so the
#: Laplace pressure has to build first -- an acoustic transient over
#: `R / c_s = 50` steps, on top of the mode-2 oscillation the case is about.
SETTLING_STEPS = 50

#: Measured `omega_0 / omega_Lamb` at R = 10 in 48^3. See the module docstring.
MEASURED_LAMB_RATIO = 0.822

#: Measured period, in steps. Lamb's is 725.5 for the same numbers.
MEASURED_PERIOD = 940.9

#: Measured `alpha / omega_0`. The case runs at mu = 0.06 because the D3Q19
#: positivity bound puts a floor under it -- see the case header -- so this is a
#: heavily damped oscillation and the correction from the observed frequency to
#: the undamped one is 6.6 %, not the 0.8 % the pre-fix case claimed.
MEASURED_DAMPING_RATIO = 0.346

#: Measured `r_z - r_x` at t = 0, which the initial shape fixes at `1.5 eps R'`.
MEASURED_INITIAL_DEFORMATION = 1.489

#: Measured equal-volume radius, `(2 r_x + r_z) / 3`, averaged over the run.
#: The droplet is laid down at 10 and settles slightly inward as the Laplace
#: pressure compresses it, exactly as in `laplace_3d`, which relaxes to 9.78.
MEASURED_MEAN_RADIUS = 9.690


@pytest.fixture(scope="module")
def run_3d():
    """One shared run of the case. It takes about five minutes on 24 cores."""
    return run_program("oscillation_3d", artifacts_dir() / "oscillation_3d", timeout=2700)


@pytest.fixture(scope="module")
def track(run_3d):
    """The droplet's semi-axes, and the mode-2 signal read off them."""
    axes = run_3d.droplet_axes()
    return {
        "timestep": axes[:, 0],
        "rx": axes[:, 1],
        "ry": axes[:, 2],
        "rz": axes[:, 3],
        "signal": axes[:, 3] - axes[:, 1],
        "radius": (2.0 * axes[:, 1] + axes[:, 3]) / 3.0,
    }


@pytest.fixture(scope="module")
def fit(track):
    return fit_damped_oscillation(track["signal"][SETTLING_STEPS:])


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_3d_starts_from_the_shape_it_was_asked_for(run_3d, track):
    """A spheroid of the prescribed deformation, on a genuinely cubic lattice."""
    assert run_3d.parameter("nz", int) == run_3d.parameter("nx", int) > 1
    assert run_3d.parameter("track_interface", str) == "true"
    assert track["signal"][0] == pytest.approx(MEASURED_INITIAL_DEFORMATION, rel=0.01)
    assert track["radius"][0] == pytest.approx(run_3d.parameter("radius"), abs=0.1)


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_3d_stays_axisymmetric(track):
    """r_x and r_y must agree to the last bit, for every step of the run.

    The deformation is about z, so the two equatorial semi-axes are the same
    quantity measured along two different lattice directions. Nothing in D3Q19
    or in the solver treats x and y differently, and this is the statement that
    nothing does: not "close", but equal.
    """
    assert np.array_equal(track["rx"], track["ry"])


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_3d_keeps_the_droplet_whole(run_3d, track):
    """The equal-volume radius must hold through the oscillation.

    `(2 r_x + r_z) / 3` is `R` exactly for a mode-2 spheroid, so it is the
    radius that should not move while the axes swing by a node and a half
    either way. A droplet losing volume to the surrounding fluid, or drifting
    off centre, would show here first.
    """
    assert np.all(np.isfinite(track["radius"]))
    assert track["radius"].mean() == pytest.approx(MEASURED_MEAN_RADIUS, abs=0.15)
    assert track["radius"].std() < 0.1


@pytest.mark.long
@pytest.mark.verification
def test_verification_oscillation_3d_rings_down_as_a_single_mode(track, fit):
    """What is fitted has to be one decaying sinusoid, or the frequency is a fiction.

    The residual is asserted against the signal's own spread rather than an
    absolute number: a second mode, a drifting centre or a droplet coming apart
    would all leave the fit unable to follow the track, and that is what would
    make the frequency below meaningless.
    """
    assert fit.residual < 0.06 * track["signal"][SETTLING_STEPS:].std()
    assert 0.0 < fit.decay_rate
    assert fit.damping_ratio == pytest.approx(MEASURED_DAMPING_RATIO, rel=0.2)


@pytest.mark.long
@pytest.mark.validation
def test_validation_oscillation_3d_frequency_against_lamb(run_3d, track, fit):
    """The measured frequency, against `24 sigma / [R^3 (3 rho_in + 2 rho_out)]`.

    Scored at the radius the case asked for, which is also the radius the
    surface tension was measured at in `laplace_3d`. The ratio is 0.822 and the
    tolerance is set to catch a change in the scheme rather than to declare
    agreement: the remaining gap closes as `1/R`, about four fifths of it
    resolution and the rest the viscous damping the positivity bound forces on
    this case. The module docstring has the three radii that show it.
    """
    lamb = lamb_frequency(
        mode=2,
        sigma=run_3d.parameter("sigma"),
        radius=run_3d.parameter("radius"),
        rho_in=run_3d.parameter("rho1"),
        rho_out=run_3d.parameter("rho2"),
    )
    assert fit.period == pytest.approx(MEASURED_PERIOD, rel=0.05)
    assert fit.undamped_frequency / lamb == pytest.approx(MEASURED_LAMB_RATIO, rel=0.05)


@pytest.mark.long
@pytest.mark.validation
def test_validation_oscillation_3d_is_far_from_a_standing_still_droplet(track, fit):
    """The oscillation has to be the thing being measured, not the noise.

    `laplace_3d` leaves spurious currents at 7e-5 in a droplet that is not
    moving. Here the semi-axes swing by more than a node, which is four orders
    of magnitude above the interface's own jitter -- stated so that a future
    change quietly killing the oscillation cannot pass as agreement.

    The peak-to-peak swing is 1.99 nodes: the droplet is released at +1.49 and
    the first return reaches -0.50. It was over 2.0 while the case ran at
    mu = 0.02, and the shortfall is the extra damping the D3Q19 positivity
    bound costs -- the release amplitude is set by the initial shape and has not
    moved. The bound below is 1.5, which a droplet that had stopped ringing
    could not reach.
    """
    assert fit.amplitude > 1.0
    assert track["signal"].max() - track["signal"].min() > 1.5
