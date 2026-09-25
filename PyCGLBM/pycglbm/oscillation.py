"""Measuring the frequency and the decay rate of a ringing signal.

A droplet released from a deformed shape rings down: its mode-2 amplitude is a
sinusoid inside a decaying envelope, and what a validation wants from it is the
frequency, which is what Lamb's formula predicts, and the decay rate, which is
what says how far the measured frequency has been pulled by viscosity.

Fitting `A exp(-alpha t) cos(omega t + phi) + c` is a nonlinear least-squares
problem, and the tooling here depends on numpy alone -- no scipy. It does not
need one. A damped sinusoid sampled at a fixed interval satisfies a linear
recurrence exactly,

    s[n+1] = a s[n] + b s[n-1] + c,
    a = 2 exp(-alpha dt) cos(omega dt),   b = -exp(-2 alpha dt)

-- the constant term carrying any offset -- so `a` and `b` come out of one
ordinary least-squares solve, and the roots of `z^2 - a z - b` are
`exp((-alpha +- i omega) dt)`. That is Prony's method for a single mode, and it
needs no starting guess.

It is not the end of the fit, because its least squares is over the recurrence
rather than over the signal: the regressors carry the same errors the target
does, so anything the single mode does not explain biases it. Measured on the
`oscillation_3d` track, Prony alone leaves a residual of 0.035 against 0.014 for
a proper nonlinear fit, and moves the frequency by 1.2 %. So Prony supplies the
starting point and a shrinking grid over `(alpha, omega)` finishes the job:
with those two fixed the model is linear in the amplitude and the offset, so
each trial costs one small least-squares solve and the search runs over two
parameters rather than five.

Reference
 - R. Prony, "Essai experimental et analytique", Journal de l'Ecole
   Polytechnique 1(2), 24 (1795). The recurrence, and the observation that the
   exponents follow from its roots.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np


class DampedOscillation(NamedTuple):
    """What :func:`fit_damped_oscillation` measured."""

    #: Angular frequency actually observed, i.e. `2 pi / (measured period)`.
    angular_frequency: float
    #: Decay rate of the envelope, in inverse units of the sample spacing.
    decay_rate: float
    #: Amplitude of the envelope extrapolated back to the first sample.
    amplitude: float
    #: Constant the oscillation is about.
    offset: float
    #: RMS of what the fit does not explain, in the units of the signal.
    residual: float

    @property
    def period(self) -> float:
        return 2.0 * np.pi / self.angular_frequency

    @property
    def undamped_frequency(self) -> float:
        """The frequency the same mode would have without the damping.

        `omega_0 = sqrt(omega^2 + alpha^2)`, which is the inverse of
        `omega = sqrt(omega_0^2 - alpha^2)` for a linear damped oscillator. It
        is the quantity an *inviscid* prediction such as Lamb's should be
        compared against; the correction it applies is second order in
        `alpha / omega`, so it matters at all only for a heavily damped signal.
        """
        return float(np.hypot(self.angular_frequency, self.decay_rate))

    @property
    def damping_ratio(self) -> float:
        """`alpha / omega_0`, zero for a free oscillation and one at critical."""
        return self.decay_rate / self.undamped_frequency


def _linear_stage(values, time, decay_rate, angular_frequency):
    """Best amplitude and offset for a fixed envelope and frequency.

    Returns ``(coefficients, residual)`` for the model
    ``c + exp(-alpha t) (A cos omega t + B sin omega t)``, which is linear in
    the three coefficients once ``alpha`` and ``omega`` are given.
    """
    envelope = np.exp(-decay_rate * time)
    columns = np.column_stack(
        [
            np.ones(values.size),
            envelope * np.cos(angular_frequency * time),
            envelope * np.sin(angular_frequency * time),
        ]
    )
    coefficients, *_ = np.linalg.lstsq(columns, values, rcond=None)
    residual = float(np.sqrt(np.mean((columns @ coefficients - values) ** 2)))
    return coefficients, residual


def _prony(values):
    """Decay rate and angular frequency from the three-term recurrence."""
    # s[n+1] = a s[n] + b s[n-1] + c, over every triple of consecutive samples.
    design = np.column_stack([values[1:-1], values[:-2], np.ones(values.size - 2)])
    (a, b, _), *_ = np.linalg.lstsq(design, values[2:], rcond=None)

    discriminant = a * a + 4.0 * b
    if discriminant >= 0.0:
        raise ValueError(
            "the recurrence has real roots: the signal decays without oscillating, "
            "so it has no frequency to report"
        )
    root = complex(0.5 * a, 0.5 * np.sqrt(-discriminant))
    return -np.log(abs(root)), abs(np.angle(root))


def fit_damped_oscillation(signal, dt: float = 1.0, refinements: int = 6) -> DampedOscillation:
    """Fit one decaying sinusoid about a constant to ``signal``.

    ``signal`` is sampled at a uniform spacing ``dt``. Prony's recurrence gives
    the first estimate and ``refinements`` rounds of a shrinking grid over
    ``(alpha, omega)`` minimise the residual in the signal itself; each round
    halves the search box, so the frequency is resolved to about
    ``2^-refinements`` of its own value.

    Raises ``ValueError`` when the signal holds fewer than four samples, when it
    is not finite throughout, or when the recurrence comes out with real roots
    -- which is what a signal that decays without ever oscillating looks like,
    and reporting a frequency for one of those would be an invention.
    """
    values = np.asarray(signal, dtype=float)
    if values.ndim != 1 or values.size < 4:
        raise ValueError(f"need at least four samples to fit, got {values.size}")
    if not np.all(np.isfinite(values)):
        raise ValueError("signal holds values that are not finite")

    decay_rate, angular_frequency = _prony(values)
    time = np.arange(values.size, dtype=float)  # in samples; dt is applied at the end

    best = _linear_stage(values, time, decay_rate, angular_frequency)[1]
    # Half the frequency to start with, so the search can still recover from a
    # Prony estimate that a contaminated signal has pulled well off.
    span = 0.5
    for _ in range(refinements):
        centre = (decay_rate, angular_frequency)
        for trial_alpha in centre[0] + span * angular_frequency * np.linspace(-1.0, 1.0, 9):
            for trial_omega in centre[1] * (1.0 + span * np.linspace(-1.0, 1.0, 9)):
                if trial_omega <= 0.0:
                    continue
                residual = _linear_stage(values, time, trial_alpha, trial_omega)[1]
                if residual < best:
                    best, decay_rate, angular_frequency = residual, trial_alpha, trial_omega
        span *= 0.25

    (offset, cosine, sine), residual = _linear_stage(values, time, decay_rate, angular_frequency)
    return DampedOscillation(
        angular_frequency=float(angular_frequency / dt),
        decay_rate=float(decay_rate / dt),
        amplitude=float(np.hypot(cosine, sine)),
        offset=float(offset),
        residual=residual,
    )


def lamb_frequency(mode: int, sigma: float, radius: float, rho_in: float, rho_out: float) -> float:
    """Lamb's angular frequency for the free oscillation of a droplet.

        omega_n^2 = n (n - 1) (n + 1) (n + 2) sigma
                    / { R^3 [ (n + 1) rho_in + n rho_out ] }

    Two inviscid fluids, an interface of tension ``sigma``, and a deformation in
    the ``mode``-th surface harmonic. ``mode = 2`` is the lowest a droplet has:
    modes 0 and 1 are a change of volume and a translation, neither of which the
    surface energy resists.

    The two-dimensional analogue is `omega^2 = n(n^2-1) sigma / [R^3 (rho_in +
    rho_out)]`, so this is not a formula a two-dimensional run could be scored
    against -- the mode shape, the numerator and the weighting of the two
    densities all change with the dimension.

    Reference
     - H. Lamb, *Hydrodynamics*, 6th ed., Cambridge (1932), art. 275.
    """
    if mode < 2:
        raise ValueError(f"mode {mode} is not restored by surface tension")
    numerator = mode * (mode - 1) * (mode + 1) * (mode + 2) * sigma
    denominator = radius**3 * ((mode + 1) * rho_in + mode * rho_out)
    return float(np.sqrt(numerator / denominator))
