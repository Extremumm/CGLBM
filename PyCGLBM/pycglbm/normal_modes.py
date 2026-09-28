"""Exact linear normal modes of two viscous fluids and the interface between them.

A capillary wave or a droplet oscillation is usually scored against the
inviscid frequency and a weak-damping rate, `2 nu k^2` or `2 n (n - 1) nu / R^2`.
Those assume the potential flow on each side is the whole story. With two
fluids it is not: the two potential flows slip past each other at the
interface, and the vortical layer that takes up the slip dissipates too. For
the cases the benchmarks run -- a heavy fluid under a light one of comparable
dynamic viscosity -- the weak-damping formulas are 12 to 15 % off, which is as
large as the errors being measured.

So the reference here is the normal mode itself. Each fluid is incompressible
and obeys the linearised Navier-Stokes equations; the velocity is split as

    u = grad phi + curl(psi e_z),   lap phi = 0,   (lap - q^2) psi = 0,
    q^2 = s rho / mu,               p = -rho s phi,

with every field proportional to `exp(s t)`. Continuity of both velocity
components and of the shear stress across the interface, and the jump of the
normal stress that the tension balances, give four homogeneous equations in the
four amplitudes; the complex root `s = -decay_rate + i angular_frequency` of
their determinant is the mode. It is found by a secant iteration started from
the weak-damping estimate.

For the droplet the vortical parts are the modified Bessel functions `I_n`
inside and `K_n` outside, of complex argument. The package depends on numpy
alone, so they are computed here from their integral representations,

    I_n(z) = (1/pi) int_0^pi exp(z cos t) cos(n t) dt,
    K_n(z) = int_0^inf exp(-z cosh t) cosh(n t) dt,        Re z > 0,

by the trapezoidal rule, which converges exponentially fast for both: the first
integrand is periodic and analytic, the second decays double-exponentially.

References
 - H. Lamb, *Hydrodynamics*, 6th ed., Cambridge (1932), arts. 275 and 355.
 - S. Chandrasekhar, *Hydrodynamic and Hydromagnetic Stability*, Oxford
   (1961), ch. X: the dispersion relation of a viscous interface.
 - A. Prosperetti, "Free oscillations of drops and bubbles: the initial-value
   problem", J. Fluid Mech. 100, 333 (1980): why the normal mode, and not the
   weak-damping rate, is the right reference at these viscosities.
"""

from __future__ import annotations

from typing import Callable, NamedTuple

import numpy as np


class NormalMode(NamedTuple):
    """A decaying oscillation `exp(-decay_rate t) cos(angular_frequency t)`."""

    decay_rate: float
    angular_frequency: float


def _secant(function: Callable[[complex], complex], start: complex) -> complex:
    """A root of `function` near `start`, by the complex secant method."""
    previous, current = start, start * (1.0 + 1.0e-3)
    f_previous, f_current = function(previous), function(current)
    for _ in range(100):
        if f_current == f_previous:
            break
        step = f_current * (current - previous) / (f_current - f_previous)
        previous, f_previous = current, f_current
        current = current - step
        f_current = function(current)
        if abs(step) < 1.0e-14 * abs(current):
            return current
    if abs(current - previous) < 1.0e-10 * abs(current):
        return current
    raise ValueError(f"the mode near s = {start} did not converge")


def _weak_start(decay_rate: float, angular_frequency: float) -> complex:
    return complex(-decay_rate, angular_frequency)


def capillary_wave(
    wavenumber: float, rho1: float, rho2: float, mu1: float, mu2: float, sigma: float
) -> NormalMode:
    """The capillary wave on a flat interface between two semi-infinite fluids.

    Fluid 1 lies below the interface, fluid 2 above; gravity is absent. The
    inviscid limit is `omega^2 = sigma k^3 / (rho1 + rho2)`, and with one fluid
    only the weak-damping limit is `2 nu k^2`.
    """
    k = float(wavenumber)
    omega_0 = np.sqrt(sigma * k**3 / (rho1 + rho2))
    start = _weak_start(2.0 * k * k * (mu1 + mu2) / (rho1 + rho2), omega_0)
    # Fixed row scales, so that the determinant stays analytic in s.
    scale_shear = 1.0 / (max(mu1, mu2) * k * k)
    scale_normal = 1.0 / ((rho1 + rho2) * omega_0)

    def determinant(s: complex) -> complex:
        m1 = np.sqrt(k * k + s * rho1 / mu1)
        m2 = np.sqrt(k * k + s * rho2 / mu2)
        ik = 1j * k
        # Unknowns: phi_1 = A e^{ky}, psi_1 = B e^{m1 y}, phi_2 = C e^{-ky},
        # psi_2 = D e^{-m2 y}; u_x = d_x phi + d_y psi, u_y = d_y phi - d_x psi.
        matrix = np.array(
            [
                [ik / k, m1 / k, -ik / k, m2 / k],
                [1.0, -1j, 1.0, 1j],
                [
                    scale_shear * mu1 * 2j * k * k,
                    scale_shear * mu1 * (m1 * m1 + k * k),
                    scale_shear * mu2 * 2j * k * k,
                    -scale_shear * mu2 * (m2 * m2 + k * k),
                ],
                [
                    scale_normal * (-(rho1 * s + 2.0 * mu1 * k * k) - sigma * k**3 / s),
                    scale_normal * (2.0 * mu1 * ik * m1 + sigma * k * k * ik / s),
                    scale_normal * (rho2 * s + 2.0 * mu2 * k * k),
                    scale_normal * 2.0 * mu2 * ik * m2,
                ],
            ],
            dtype=complex,
        )
        return complex(np.linalg.det(matrix))

    s = _secant(determinant, start)
    return NormalMode(decay_rate=-s.real, angular_frequency=s.imag)


def bessel_i(order: int, z: complex, points: int = 256) -> complex:
    """The modified Bessel function of the first kind, `I_n(z)`, integer n."""
    theta = np.linspace(0.0, np.pi, points + 1)
    values = np.exp(z * np.cos(theta)) * np.cos(order * theta)
    return complex((values.sum() - 0.5 * (values[0] + values[-1])) / points)


def bessel_k(order: int, z: complex, step: float = 0.01) -> complex:
    """The modified Bessel function of the second kind, `K_n(z)`, `Re z > 0`."""
    if z.real <= 0.0:
        raise ValueError(f"K_n needs Re z > 0, got {z}")
    # Past T the integrand is below exp(-700) of its value at the origin.
    end = np.arccosh(max(1.0, 700.0 / z.real + 1.0)) + 1.0
    t = np.arange(0.0, end + step, step)
    values = np.exp(-z * np.cosh(t)) * np.cosh(order * t)
    return complex(step * (values.sum() - 0.5 * values[0]))


def _bessel_pair(function, order: int, z: complex) -> tuple[complex, complex]:
    """`F_n(z)` and `dF_n/dz`, from `F_{n-1}` and `F_{n+1}`."""
    value = function(order, z)
    below = function(abs(order - 1), z)
    above = function(order + 1, z)
    if function is bessel_k:
        return value, -0.5 * (below + above)
    return value, 0.5 * (below + above)


def droplet_mode(
    mode: int,
    radius: float,
    rho_in: float,
    rho_out: float,
    mu_in: float,
    mu_out: float,
    sigma: float,
) -> NormalMode:
    """The `mode`-th oscillation of a two-dimensional droplet (a liquid cylinder).

    The interface is `r = R + a cos(n theta)`. The inviscid limit is
    `omega^2 = n (n^2 - 1) sigma / [R^3 (rho_in + rho_out)]`, and the
    weak-damping limit of a droplet in vacuum is `2 n (n - 1) nu / R^2`.
    """
    n = int(mode)
    if n < 2:
        raise ValueError(f"mode {n} is not restored by surface tension")
    r = float(radius)
    omega_0 = np.sqrt(n * (n * n - 1) * sigma / (r**3 * (rho_in + rho_out)))
    start = _weak_start(
        2.0 * n * ((n - 1) * mu_in + (n + 1) * mu_out) / ((rho_in + rho_out) * r * r), omega_0
    )
    scale_shear = r * r / max(mu_in, mu_out)
    scale_normal = r / ((rho_in + rho_out) * omega_0)

    def columns(inside: bool, s: complex) -> dict[str, tuple[complex, complex]]:
        """Per unknown (potential, stream function): u_r, u_theta and their r-derivatives."""
        if inside:
            f, fp, fpp = r**n, n * r ** (n - 1), n * (n - 1) * r ** (n - 2)
            q = np.sqrt(s * rho_in / mu_in)
            value, derivative = _bessel_pair(bessel_i, n, q * r)
        else:
            f, fp, fpp = r**-n, -n * r ** (-n - 1), n * (n + 1) * r ** (-n - 2)
            q = np.sqrt(s * rho_out / mu_out)
            value, derivative = _bessel_pair(bessel_k, n, q * r)
        # The stream function is normalised to one at r = R.
        g, gp = 1.0, q * derivative / value
        gpp = q * q * g + n * n * g / r**2 - gp / r
        # phi = F cos(n theta), psi = G sin(n theta):
        # u_r = F' + (n/r) G, u_theta = -(n/r) F - G'.
        return {
            "phi": (f, fp, fpp),
            "ur": (fp, n / r * g),
            "ut": (-n / r * f, -gp),
            "urp": (fpp, -n / r**2 * g + n / r * gp),
            "utp": (n / r**2 * f - n / r * fp, -gpp),
        }

    def determinant(s: complex) -> complex:
        inner, outer = columns(True, s), columns(False, s)
        rows = []
        for key in ("ur", "ut"):
            rows.append([inner[key][0], inner[key][1], -outer[key][0], -outer[key][1]])

        def shear(side, mu, k):
            return mu * (side["utp"][k] - side["ut"][k] / r - n / r * side["ur"][k])

        rows.append(
            [
                scale_shear * shear(inner, mu_in, 0),
                scale_shear * shear(inner, mu_in, 1),
                -scale_shear * shear(outer, mu_out, 0),
                -scale_shear * shear(outer, mu_out, 1),
            ]
        )

        def normal(side, rho, mu, k):
            potential = side["phi"][0] if k == 0 else 0.0
            return rho * s * potential + 2.0 * mu * side["urp"][k]

        tension = sigma * (n * n - 1) / r**2 / s
        rows.append(
            [
                scale_normal * (-normal(inner, rho_in, mu_in, 0) - tension * inner["ur"][0]),
                scale_normal * (-normal(inner, rho_in, mu_in, 1) - tension * inner["ur"][1]),
                scale_normal * normal(outer, rho_out, mu_out, 0),
                scale_normal * normal(outer, rho_out, mu_out, 1),
            ]
        )
        # Column scales: the potentials to their velocity at r = R.
        matrix = np.array(rows, dtype=complex)
        matrix[:, 0] *= r ** (1 - n)
        matrix[:, 2] *= r ** (1 + n)
        return complex(np.linalg.det(matrix))

    s = _secant(determinant, start)
    return NormalMode(decay_rate=-s.real, angular_frequency=s.imag)
