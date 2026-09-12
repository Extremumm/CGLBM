"""The quasi-static MHD Rayleigh-Taylor dispersion relation, and a fit for it.

The relation is the finite-depth inviscid one used by the ``mrt`` campaign in
``basiliskMHD/QSMHD``:

    rho1^.5 w^1.5 sqrt(w rho1 + B^2 s1) coth(k h1 / sqrt(1 + s1 B^2/(rho1 w)))
  + rho2^.5 w^1.5 sqrt(w rho2 + B^2 s2) coth(k h2 / sqrt(1 + s2 B^2/(rho2 w)))
  - [(rho1 - rho2) g k - gamma k^3]  =  0,

whose positive root is the growth rate. Two things about it shape everything
downstream. It carries no viscosity, so a simulation can only match it to the
extent that the flow is inviscid; and it is dimensionally consistent, so it may
be evaluated in lattice units directly, which is what is done here.

Setting ``B = 0`` reduces the two magnetic terms to ``rho_k coth(k h_k)`` times
``w^2`` and recovers the hydrodynamic rate, so one function covers both.

It is restated here rather than imported because the campaign it comes from is a
separate project on a separate solver, and a test that silently depended on a
path outside this repository would be a test that breaks for the wrong reason.
"""

from __future__ import annotations

import math

import numpy as np


def _coth(x: float) -> float:
    """``coth`` without the overflow at large argument or the pole at zero."""
    if x < 1.0e-12:
        return 1.0 / (x + 1.0e-300)
    if x > 25.0:
        return 1.0
    return math.cosh(x) / math.sinh(x)


def driving_term(rho1: float, rho2: float, gravity: float, gamma: float, k: float) -> float:
    """Gravity pulling the arrangement apart, less surface tension holding it."""
    return (rho1 - rho2) * gravity * k - gamma * k**3


def _phase_term(rho: float, conductivity: float, depth: float, k: float, field: float,
                omega: float) -> float:
    """One layer's contribution to the relation."""
    alpha = conductivity * field**2 / (rho * omega)
    return (
        rho**0.5
        * omega**1.5
        * math.sqrt(omega * rho + field**2 * conductivity)
        * _coth(k * depth / math.sqrt(1.0 + alpha))
    )


def growth_rate(
    *,
    rho1: float,
    rho2: float,
    gravity: float,
    gamma: float,
    k: float,
    depth: float,
    field: float = 0.0,
    conductivity1: float = 0.0,
    conductivity2: float = 0.0,
) -> float:
    """The positive root, or zero when the arrangement is stable.

    Both layers are of the same depth here, which is what a centred interface in
    a channel gives.
    """
    drive = driving_term(rho1, rho2, gravity, gamma, k)
    if drive <= 0.0:
        return 0.0

    def residual(omega: float) -> float:
        return (
            _phase_term(rho1, conductivity1, depth, k, field, omega)
            + _phase_term(rho2, conductivity2, depth, k, field, omega)
            - drive
        )

    low, high = 1.0e-14, 1.0e-6
    for _ in range(200):
        if residual(high) > 0.0:
            break
        high *= 2.0
    else:
        return 0.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        if residual(middle) > 0.0:
            high = middle
        else:
            low = middle
    return 0.5 * (low + high)


def interface_height(phase: np.ndarray) -> np.ndarray:
    """Height of the ``phi_N = 0`` contour in each column of a ``[y, x]`` slice.

    The heavy component is on top, so each column runs from -1 at the bottom to
    +1 at the top and crosses once. The crossing is interpolated linearly
    between the two nodes bracketing it; a column that never crosses gives
    ``nan`` rather than a plausible number.
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


def mode_amplitude(phase: np.ndarray) -> float:
    """Amplitude of the single cosine the interface was seeded with.

    Taken as the first Fourier coefficient of the interface height rather than
    as its peak-to-trough range, so that the number is unaffected by the spike
    and bubble becoming asymmetric -- which they do well before the mode leaves
    the linear regime, and which would otherwise be read as growth.
    """
    heights = interface_height(phase)
    if np.isnan(heights).any():
        return float("nan")
    nx = heights.size
    return float(2.0 * np.abs(np.fft.rfft(heights - heights.mean())[1]) / nx)


def fit_growth_rate(times, amplitudes, seed: float, upper: float):
    """Least-squares slope of ``log a`` against ``t``, over the linear window.

    The window runs from three times the seed amplitude to ``upper``, which is
    where the amplitude is still small against the wavelength.

    Three times, rather than just clear of the seed, because the case starts
    from rest with a sharp interface and the velocity field takes time to
    establish. Measured on the zero-field run, the apparent rate climbs through
    5.7e-4, 5.5e-4, 8.7e-4 and 9.1e-4 over the first 2000 steps before settling;
    a window opened at the seed reads the average of the startup and the growth,
    which is 23 % low. By three times the seed the transient has passed.

    Returns ``(rate, count)`` with the number of samples the fit used, so that a
    caller can tell a fit from an accident.
    """
    times = np.asarray(times, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)
    inside = np.isfinite(amplitudes) & (amplitudes > 3.0 * seed) & (amplitudes < upper)
    if inside.sum() < 4:
        return float("nan"), int(inside.sum())
    slope, _ = np.polyfit(times[inside], np.log(amplitudes[inside]), 1)
    return float(slope), int(inside.sum())
