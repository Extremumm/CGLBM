# Two-population solvers

`TwoPopulationSolver` (`src/lbm/two_population_solver.h`, D2Q9) and
`TwoPopulationSolver3D` (`src/lbm/two_population_solver_3d.h`, D3Q19 or
D3Q27) are the classical colour-gradient model of Grunau et al. (1993) in the
form of Ba et al. (2016): one distribution per fluid, `f1` and `f2`, streamed
separately, with the density ratio carried by each fluid's rest-particle
weight `α_k`. There is no equation of state to invert.

Programs: `laplace_high_ratio` (2D); `laplace_3d`, `oscillation_3d`,
`rayleigh_taylor_3d`, `hartmann`, `magnetic_rayleigh_taylor` (3D).

## What they read from `CaseConfig`

`physics.rho1`, `rho2`, `nu`, `nu2`, `sigma`, `radius`, `beta`, `alpha2`,
`gravity`; `boundary`, `parallel`, `output_precision`, `steps`, `interval`.
The 2D solver reads `nx`, `ny`, `stencil` and `initial_phase`; the 3D solver
reads `nx`, `ny`, `nz`, `stencil_3d`, `initial_phase_3d`, `lattice_3d`,
`physics.body_force`, `mhd`, `track_interface` and `write_vtk_field`.
`c1`, `c2`, `p1_inf`, `p2_inf`, `initial_state`, `initial_profile_field`,
`interface_field`, `surface_tension`, `viscosity_mixing`, `recolouring` and
`nu_b`/`nu_b2` play no part. The initial phase is always read as the
volume-fraction indicator φ_N = 2c − 1.

## Construction

Both constructors throw `std::invalid_argument` when:

- a lattice dimension is not positive;
- the initial phase function is empty;
- the boundary is `WallY` and `ny < 3` (2D) or `ny < 2·reach + 1` with
  `reach = stencil_reach_3d(stencil_3d)` (3D);
- `rho1` or `rho2` is not positive;
- `α₁` or `α₂` falls outside [0, 1].

The rest weights follow from the density ratio:

```
α₂ = physics.alpha2
α₁ = 1 − (1 − α₂) ρ₂/ρ₁                  so that ρ₁/ρ₂ = (1 − α₂)/(1 − α₁)
```

| | D2Q9 (`TwoPopulationSolver`) | D3Q19 | D3Q27 |
|---|---|---|---|
| `(c_s^k)²` | `3(1 − α_k)/5` | `(1 − α_k)/2` | `9(1 − α_k)/19` |
| rest weight `φ₀^k` | `α_k` | `α_k` | `α_k` |
| axial `φ^k` | `(1 − α_k)/5` | `(1 − α_k)/12` | `2(1 − α_k)/19` |
| diagonal/edge `φ^k` | `(1 − α_k)/20` | `(1 − α_k)/24` | `(1 − α_k)/38` |
| corner `φ^k` | — | — | `(1 − α_k)/152` |

The bulk pressures `ρ_k(c_s^k)²` are then equal, so the pressure is continuous
across the interface at any ratio. The bulk dynamic viscosities are
`μ₁ = ρ₁ν₁` and `μ₂ = ρ₂ν₂` (`ν₂ = ν₁` when `nu2` is negative).

## Equilibrium

2D (Ba et al. Eq. 14):

```
f_k^eq = ρ_k φ_k + ρ_k w_k [3(ξ·u)(1 + λ(3|ξ|² − 4)) + 4.5(ξ·u)² − 1.5u²],   λ = (3(c_s^k)² − 1)/2
```

### The 3D equilibrium

```
f_q^eq = ρ_k φ_q + ρ_k w_q [3(e·u)(1 + λ_k(3|e|² − 5)) + 4.5(e·u)² − 1.5u²]
```

`λ_k` is not tabulated. The constructor sums, over the lattice it was given,

```
S = Σ_q w_q e_x² e_y²,     T = Σ_q w_q e_x² e_y² (3|e_q|² − 5),
λ_k = ((c_s^k)² − 3S) / (3T)
```

which is `3(c_s^k)² − 1` on D3Q19 and half that on D3Q27. It throws if
`T = 0`. `enhancement(k)` returns `λ_k`; `equilibrium_for_test` exposes the
equilibrium to the unit tests.

## Time step

2D and 3D:

```
step():  densities → update_colour_gradient → surface_force → update_velocity
         → collide → recolor → stream
```

In 3D with a magnetic field, `update_velocity(false)` runs before
`surface_force` and `update_velocity(true)` after it (see below).

| Stage | What it does |
|---|---|
| `densities()` | `ρ_k = Σf_k`, `ρ = ρ₁ + ρ₂`, `p = ρ₁(c_s¹)² + ρ₂(c_s²)²`, `φ_N = (ρ₁/ρ₁⁰ − ρ₂/ρ₂⁰)/(ρ₁/ρ₁⁰ + ρ₂/ρ₂⁰)` (0 where the denominator vanishes), ρ_k⁰ the bulk densities `physics.rho1`, `rho2` |
| `update_colour_gradient()` | `C = ∇φ_N` with the configured stencil (one-sided at walls); `n = −C/|C|` where `|C| > kInterfaceGradientFloor`, else 0 |
| `surface_force()` | curvature `K`, then `F = −½σK C` on interface nodes (`|C| > kInterfaceGradientFloor`), 0 elsewhere; minus `ρg` along y; 3D also adds `body_force` and, with MHD, the Lorentz force |
| `update_velocity()` | `u = (Σ(f1 + f2)ξ + F dt/2)/ρ` (Ba et al. Eq. 29) |
| `collide()` | BGK for each fluid at one rate, with Guo forcing split by mass share |
| `recolor()` | Latva-Kokko segregation with the mixture rest weight |
| `stream()` | stream both populations into new arrays |

**Curvature.** 2D: `K = n_xn_y(∂_yn_x + ∂_xn_y) − n_x²∂_yn_y − n_y²∂_xn_x`. 3D:
`K = −(∇·n − n_an_b∂_bn_a)`, the surface divergence, which is `−2/R` on a
sphere, so Laplace's law reads `2σ/R`.

**Collision.** With `c = (1 + φ_N)/2`:

```
μ = cμ₁ + (1 − c)μ₂,    τ = μ/(p dt) + ½,    ω = 1/τ
S_q = w_q [((e − u) + (e·u)e/c_s²)·F] / c_s²                     (Guo)
f1_q += −ω(f1_q − f1_q^eq) + (1 − ω/2) dt · (ρ₁/ρ) S_q
f2_q += −ω(f2_q − f2_q^eq) + (1 − ω/2) dt · (ρ₂/ρ) S_q
```

`c_s² = 1/3` here is the lattice's, from `config.units`.

**Recolouring** (Ba et al. Eq. 30, adapted). Where `|C| > kGradientEpsilon`:

```
Φ_q = (ρ₁φ_q¹ + ρ₂φ_q²)/ρ                                        rest_weight()
f1_q ← (ρ₁/ρ)(f1_q + f2_q) + Φ_q · βρ₁ρ₂/(ρ|C|) · (e_q·C)/|e_q|  (q ≠ 0)
f2_q ← (f1_q + f2_q) − f1_q
```

The push uses the mixture rest weight `Φ_q` rather than the lattice weight
`w_q`; this keeps both populations non-negative for `β ≤ 1` at any density
ratio, and reduces to Latva-Kokko & Rothman exactly when the fluids share a
rest weight.

**Streaming and walls.** x (and z in 3D) are periodic. With `WallY`, a
population whose destination row is outside `[0, ny − 1]` is written back in
the reversed direction (`k ∓ 2` in 2D, `lattice.opposite[q]` in 3D) at the same
row but with its tangential step kept: it lands at `(i + ξ_x, j[, k + ξ_z])`.
For flow that does not vary along the wall this is the half-way bounce-back;
for a diagonal population where it does, the reflected population is one node
along the wall from where the textbook scheme would put it. `Solver` streams
the same way.

## Initialisation and output

`initialize()` sets `c = (1 + clamp(φ_N))/2` at every node, `ρ₁ = cρ₁⁰`,
`ρ₂ = (1 − c)ρ₂⁰`, the fluid at rest, and each population to its equilibrium;
then it runs `densities`, `update_colour_gradient`, `surface_force` and
`update_velocity`.

`step()` leaves the populations streamed but the macroscopic fields as they
were before the collision. `refresh()` recomputes ρ_k, ρ, p, φ_N and u from the
populations; `run()` calls it before writing.

`state()` points at `rho_`, `u_`, `phi_n_` and `p_`: the phase written to disk
is φ_N, not a colour field. Accessors: `density()`, `velocity()`, `phase()`
(φ_N), `pressure()`, `population(k)`, `component_density(k)`, `alpha1()`,
`alpha2()`; in 3D also `lattice()`, `sound_speed_squared(k)`,
`enhancement(k)`, `mhd()` and `interface_axes(radii)`.

`TwoPopulationSolver::run()`:

```
initialize(); write_grids(0)
for t = 1 .. steps: step(); if t % interval == 0: print "Step t"; refresh(); write_grids(t)
```

`TwoPopulationSolver3D::run()` writes the `z = nz/2` slice through the same
`CsvWriter` (`write_midplane`), so the 2D post-processing reads it unchanged:

```
initialize(); write_midplane(0); [write_vtk_3d("field", 0)]
if track_interface: open_droplet_track(); write_axes(0, interface_axes())
for t = 1 .. steps:
    step()
    if track_interface: densities(); write_axes(t, interface_axes())
    if t % interval == 0:
        print "Step t" [+ "  potential: N iterations, residual R, charge imbalance Q"
                        + "  [did not reach tolerance]" when the solve did not converge]
        refresh(); write_midplane(t); [write_vtk_3d("field", t)]
```

`interface_axes(radii)` walks from the centre `(nx/2, ny/2, nz/2)` along +x,
+y and +z and returns the first place φ_N falls from ≥ 0 to < 0, linearly
interpolated between the two nodes; NaN for an axis with no crossing.

## The magnetic coupling (3D)

With `config.mhd.enabled`, the constructor creates a `QuasiStaticMhd3D`
(`mhd()`), with insulating walls when the boundary is `WallY`. Each step then:

1. `update_velocity(false)`: `u = Σ(f1 + f2)e/ρ`, the bare momentum, without the
   half force — the magnetic force must not depend on itself;
2. `surface_force()` builds capillary, gravity and body forces, then calls
   `mhd_->solve(u, φ_N)` and `mhd_->add_lorentz_force(force_)`;
3. `update_velocity(true)` adds the half force of the full `F`.

Without a field, step 1 is skipped. See [Magnetohydrodynamics](mhd.md).
