"""Unit tests for the exact viscous normal modes the oscillation benchmarks use."""

import numpy as np
import pytest
from pycglbm import CaseOutput
from pycglbm.normal_modes import bessel_i, bessel_k, capillary_wave, droplet_mode

#: The tension every benchmark runs with, in lattice units.
SIGMA = 0.27683423442876637


@pytest.mark.unit_test
@pytest.mark.parametrize("z", [0.07 + 0.07j, 1.3 + 1.4j, 6.7 + 6.7j, 30.6 + 25.8j])
@pytest.mark.parametrize("order", [1, 2, 3])
def test_unit_test_normal_modes_bessel_wronskian(z, order):
    """``I_n K_n' - I_n' K_n = -1/z`` holds for the quadratures, to rounding.

    It ties the two independent integrals together, over the range of arguments
    the droplet mode evaluates them at.
    """
    i_n, k_n = bessel_i(order, z), bessel_k(order, z)
    i_prime = 0.5 * (bessel_i(order - 1, z) + bessel_i(order + 1, z))
    k_prime = -0.5 * (bessel_k(order - 1, z) + bessel_k(order + 1, z))
    wronskian = i_n * k_prime - i_prime * k_n
    assert abs(wronskian * z + 1.0) < 1.0e-10


@pytest.mark.unit_test
def test_unit_test_normal_modes_bessel_known_values():
    """Two tabulated values: ``I_0(1) = 1.2660658778``, ``K_0(1) = 0.4210244382``."""
    assert bessel_i(0, 1.0 + 0.0j).real == pytest.approx(1.2660658777520084, rel=1.0e-13)
    assert bessel_k(0, 1.0 + 0.0j).real == pytest.approx(0.42102443824070834, rel=1.0e-13)


@pytest.mark.unit_test
@pytest.mark.parametrize("nu", [1.0e-3, 1.0e-4])
def test_unit_test_normal_modes_capillary_wave_free_surface_limit(nu):
    """With nothing above it, a slightly viscous wave decays at ``2 nu k^2``.

    The correction is of order ``(nu k^2 / omega)^(1/2)``: 1.7 % at 1e-3, 0.5 %
    at 1e-4 for this wavelength.
    """
    k = 2 * np.pi / 64
    mode = capillary_wave(k, 1.0, 1.0e-6, nu, 1.0e-12, SIGMA)
    assert mode.decay_rate / (2 * nu * k * k) == pytest.approx(1.0, abs=2.0e-2)
    assert mode.angular_frequency == pytest.approx(np.sqrt(SIGMA * k**3), rel=1.0e-3)


@pytest.mark.unit_test
@pytest.mark.parametrize("nu", [1.0e-3, 1.0e-4])
def test_unit_test_normal_modes_droplet_free_drop_limit(nu):
    """A slightly viscous droplet in vacuum decays at ``2 n (n - 1) nu / R^2``."""
    mode = droplet_mode(2, 20.0, 1.0, 1.0e-6, nu, 1.0e-10, SIGMA)
    assert mode.decay_rate / (4 * nu / 400.0) == pytest.approx(1.0, abs=1.0e-2)
    assert mode.angular_frequency == pytest.approx(np.sqrt(6 * SIGMA / 8000.0), rel=1.0e-3)


@pytest.mark.unit_test
def test_unit_test_normal_modes_benchmark_references():
    """The references the benchmarks are scored against, as scipy computes them.

    These come from an independent solution with scipy's Bessel functions and
    root finder; the weak-damping formulas give 9.15e-5 -> 1.05e-4 and
    1.91e-5 -> 2.15e-5 instead, which is why they are not used.
    """
    wave = capillary_wave(2 * np.pi / 64, 100.0, 1.0, 0.5, 0.05, SIGMA)
    assert wave.decay_rate == pytest.approx(9.1519e-5, rel=1.0e-4)
    assert wave.angular_frequency == pytest.approx(1.5884e-3, rel=1.0e-4)
    droplet = droplet_mode(2, 20.0, 1000.0, 1.0, 2.0, 0.05, SIGMA)
    assert droplet.decay_rate == pytest.approx(1.91115e-5, rel=1.0e-4)
    assert droplet.angular_frequency == pytest.approx(4.5304e-4, rel=1.0e-4)


def _lamb_free_surface(k, nu, gravity, sigma):
    """The root of Lamb's exact free-surface relation (art. 349),
    ``(s + 2 nu k^2)^2 + g k + sigma k^3 = 4 nu^2 k^3 (k^2 + s / nu)^(1/2)``,
    by Newton's method from the inviscid root."""
    s = complex(-2.0 * nu * k * k, np.sqrt(gravity * k + sigma * k**3))
    for _ in range(50):
        m = np.sqrt(k * k + s / nu)
        f = (s + 2 * nu * k * k) ** 2 + gravity * k + sigma * k**3 - 4 * nu * nu * k**3 * m
        df = 2 * (s + 2 * nu * k * k) - 2 * nu * k**3 / m
        s -= f / df
    return s


@pytest.mark.unit_test
@pytest.mark.parametrize("sigma", [0.0, SIGMA])
def test_unit_test_normal_modes_gravity_wave_against_lamb(sigma):
    """With gravity, a free surface's viscous wave is Lamb's exact root.

    Gravity enters the normal-stress balance beside the tension, as
    ``sigma k^2 + (rho1 - rho2) g``; the free surface is fluid 1 under a
    fluid a millionth as dense.
    """
    k, nu, gravity = 2 * np.pi / 64, 1.0e-3, 1.0e-4
    mode = capillary_wave(k, 1.0, 1.0e-9, nu, 1.0e-15, sigma, gravity=gravity)
    exact = _lamb_free_surface(k, nu, gravity, sigma)
    assert mode.decay_rate == pytest.approx(-exact.real, rel=1.0e-6)
    assert mode.angular_frequency == pytest.approx(exact.imag, rel=1.0e-6)


@pytest.mark.unit_test
def test_unit_test_normal_modes_rayleigh_taylor_limits():
    """The heavier fluid on top: the mode grows, at ``(A g k)^(1/2)`` without
    viscosity or tension, at ``(A g k - sigma k^3 / (rho1 + rho2))^(1/2)`` with
    tension alone, and slower with viscosity; below the capillary cutoff the
    interface is stable again."""
    k, gravity = 2 * np.pi / 64, 1.0e-5
    atwood = 999.0 / 1001.0
    inviscid = capillary_wave(k, 1.0, 1000.0, 1.0e-7, 1.0e-7, 0.0, gravity=gravity)
    assert inviscid.angular_frequency == 0.0
    assert -inviscid.decay_rate == pytest.approx(np.sqrt(atwood * gravity * k), rel=1.0e-4)
    tension = capillary_wave(k, 1.0, 1000.0, 1.0e-7, 1.0e-7, SIGMA, gravity=gravity)
    expected = np.sqrt(atwood * gravity * k - SIGMA * k**3 / 1001.0)
    assert -tension.decay_rate == pytest.approx(expected, rel=1.0e-4)
    viscous = capillary_wave(k, 1.0, 1000.0, 0.05, 2.0, SIGMA, gravity=gravity)
    assert 0.0 < -viscous.decay_rate < -tension.decay_rate
    short = capillary_wave(2 * np.pi / 16, 1.0, 1000.0, 0.05, 2.0, SIGMA, gravity=gravity)
    assert short.decay_rate > 0.0 and short.angular_frequency > 0.0


@pytest.mark.unit_test
def test_unit_test_normal_modes_reject_mode_below_two():
    """Modes 0 and 1 are not restored by tension."""
    with pytest.raises(ValueError):
        droplet_mode(1, 20.0, 1000.0, 1.0, 2.0, 0.05, SIGMA)


@pytest.mark.unit_test
def test_unit_test_files_mode_track_reads_nan(tmp_path):
    """``mode.csv`` comes back as timesteps and signal, a diverged run's nan included."""
    (tmp_path / "mode.csv").write_text("timestep,amplitude\n0,0.3\n50,0.29\n100,nan\n")
    timesteps, signal = CaseOutput(tmp_path).mode_track()
    assert timesteps.tolist() == [0.0, 50.0, 100.0]
    assert signal[:2].tolist() == [0.3, 0.29]
    assert np.isnan(signal[2])
