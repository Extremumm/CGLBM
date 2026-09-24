# Numerical method and code map

The solver implements the improved color-gradient (CG) lattice Boltzmann model
of Lafarge et al. (2021) on the **D2Q9** lattice. Two immiscible components —
conventionally "red" (1) and "blue" (2) — are tracked by a single set of
distribution functions plus a phase field, rather than by two separate
populations.

## Lattice

Discrete velocities, with $c = \Delta x/\Delta t$:

$$
\vec{\xi}_i = c \begin{cases}
(0,0) & i=0\\
(\pm1,0),\,(0,\pm1) & i=1..4\\
(\pm1,\pm1),\,(\pm1,\mp1) & i=5..8
\end{cases}
\qquad
w_i = \begin{cases} 4/9 & i=0\\ 1/9 & i=1..4\\ 1/36 & i=5..8\end{cases}
$$

Speed of sound $c_s = \Delta x / (\sqrt{3}\,\Delta t)$. All cases use
$\Delta x = \Delta t = 1$ (lattice units).

## Fields

| Symbol | Array | Meaning |
|---|---|---|
| $f_i$ | `f[Lx][Ly][Q]` | sum (total) distribution function |
| $g_i$ | `g[Lx][Ly][Q]` | difference (color) distribution function |
| $f_i^{eq}$ | `f_eq[Lx][Ly][Q]` | equilibrium distribution |
| $\rho$ | `rho`, `rho_mdt` | mixture density, current and previous step |
| $p$ | `p`, `p_mdt` | pressure, current and previous step |
| $\vec u$ | `u[Lx][Ly][2]` | velocity |
| $\phi$ | `phi[Lx][Ly]` | phase field, $\phi = +1$ in fluid 1, $-1$ in fluid 2; a **mass**-fraction difference $Y_1 - Y_2$ |
| $\psi$ | `psi[Lx][Ly]` | `laplace` only: normalised colour field, the **volume**-fraction difference $\alpha_1 - \alpha_2$ |
| $\hat n$ | `normal_x`, `normal_y` | `laplace` only: interface normal $\nabla\psi/\lvert\nabla\psi\rvert$ |
| $\Omega^{(1,2,3)}$ | `omega_1/2/3` | the three collision contributions (`laplace` has no `omega_2`) |
| $S_i$, $\vec F$ | `S`, `F` | forcing term and body force (in `laplace`, the surface force) |

## Collision

The collision operator is split into three successive contributions, applied in
this order:

1. **$\Omega^{(1)}$ — regularized relaxation** (`collide`). The non-equilibrium
   part $f_i - f_i^{eq} + \tfrac12 S_i$ is projected onto the second-order
   Hermite polynomials $H_\nu$, $H_b$, $H_{xy}$ and rebuilt from those three
   moments, so higher-order (ghost) moments are discarded. The deviatoric and
   trace parts relax at separate rates set by the shear and bulk viscosities,
   $\tau_\nu = \rho\nu/(p\,\Delta t) + 1/2$ and
   $\tau_b = \rho\nu_b/(p\,\Delta t) + 1/2$ — note both are built on the local
   pressure from the equation of state, not on $\rho c_s^2$. In `laplace` each
   component has its own viscosity and $\nu = Y_1\nu_1 + Y_2\nu_2$, so that
   $\rho\nu = \alpha_1\mu_1 + \alpha_2\mu_2$.
2. **$\Omega^{(2)}$ — surface tension** (`collide_surface`). A perturbation
   built on the color gradient $\nabla\phi$ that introduces the surface tension
   $\sigma$ and reproduces the Laplace jump. `laplace` replaces it by a body
   force, the divergence of the capillary stress built on $\nabla\psi$, which
   enters through `S_F` (`calSurfaceForce`); see
   [the high-density-ratio scheme](#the-high-density-ratio-scheme).
3. **$\Omega^{(3)}$ — recoloring** (`recolor`). Redistributes the components
   along the color gradient to keep the interface sharp. `epsilon` guards the
   division by the gradient norm; `ch_width_ope` sets the operating interface
   thickness (`ch_width_init` is used only when initialising). In `laplace` the
   direction is the normal of $\psi$ rather than of $\phi$.

The forcing term assembled by `force()` is the sum of three parts:

$$S_i = S_i^{F} + S_i^{Sp} + S_i^{t}$$

- `S_F` — the body force: gravity `a_g` in the gravity cases, the surface force
  in `laplace`, zero otherwise;
- `S_Sp` — a spatial correction for the third-order moment, which the D2Q9
  lattice resolves improperly;
- `S_t` — a temporal correction $\propto (p - p^{-\Delta t}) - (\rho -
  \rho^{-\Delta t})c_s^2$, which is why `p_mdt` and `rho_mdt` are carried from
  the previous step (saved at the end of `stream`).

## Equation of state

Each component is an ideal gas with its own sound speed `c1`, `c2` and its own
pressure at infinity `p1_inf`, `p2_inf`. The mixture pressure is recovered in
`calPhaseField` by solving the quadratic mixing relation — the same expression
implemented independently in `postprocessing/pressure_jump_lafarge.py`, which is
therefore an analytic check on the initial pressure field.

For the Laplace benchmark the initial `p1_inf` is offset by $-\sigma/R$ so that
the initial condition already satisfies the pressure jump.

The same relation defines each component's density at a given pressure,
$\rho_k(p) = (p + p_{k,\infty})/c_k^2$, and so the volume fractions
$\alpha_k = \rho Y_k/\rho_k(p)$ with $Y_1 = (1+\phi)/2$. `src/lbm/mixture.h`
computes them, builds an initial state that is in equilibrium with the equation
of state, and gives the normalised colour field $\psi = \alpha_1 - \alpha_2$.

## Time loop

`runSimulation()` executes, per step:

```
force()            →  S_i from the body force F
collide()          →  Ω(1)  relaxation towards f_eq
collide_surface()  →  Ω(2)  surface-tension perturbation      (not in laplace)
recolor()          →  Ω(3)  interface sharpening
stream()           →  propagation along ξ_i
calMacroscopic()   →  ρ, u from moments of f (force-corrected; laplace: ρ and momentum)
calPhaseField()    →  φ and p from the equation of state      (laplace: and ψ)
calSurfaceForce()  →  laplace: normals, capillary stress T, F = ∇·T
calVelocity()      →  laplace: u, with the half-force correction
calEquilibrium()   →  f_eq for the next step
```

In `laplace` the velocity waits for the surface force, which needs $\psi$,
which needs the pressure.

Output is written every `interval` steps by `outputDataCSV` (and optionally
`outputVTK`).

## Boundary conditions

Neighbour lookups in `collide_surface` and `recolor` (used for the color
gradient) are periodic in both directions in every case: `(i + ξ_x + Lx) % Lx`,
`(j + ξ_y + Ly) % Ly`. The propagation step differs per case:

| Program | x | y (in `stream`) |
|---|---|---|
| `laplace` | periodic | periodic |
| `capillary`, `gravity_capillary` | periodic | half-way bounce-back |
| `rayleigh_taylor`, `..._omp` | periodic | half-way bounce-back |

Half-way bounce-back is implemented inline: at `j == 0` the populations moving
downwards (`k = 4, 7, 8`) stay on the node and are reflected into `k-2`, and at
`j == Ly-1` the upward ones (`k = 2, 5, 6`) are reflected into `k+2`. This
places the resting wall halfway between the last fluid node and the ghost node.
- `applyAbsorbingBoundary()` damps acoustic waves on the outermost `numBoundary`
  nodes. It is defined in every case but **commented out** of the time loop;
  enable it if the initial pressure transient reflects back into the domain.

## Case-specific validation

| Case | Analytic reference |
|---|---|
| `laplace` | $\Delta p = \sigma / R$ across the droplet interface, at density ratios 20 and $10^4$ |
| `capillary` | $T_{theo} = 2\pi\sqrt{(\rho_1+\rho_2)\,r^3/(6\sigma)}$, printed at startup and compared against `interface.csv` |
| `gravity_capillary` | Balance of the Laplace jump against the hydrostatic head |
| `rayleigh_taylor` | Linear growth rate of the instability, $\sigma = 0$ |

## The equation of state and the density ratio

The colour-gradient method is classically limited to modest density ratios
because the density contrast is carried by the rest-particle weight of the
equilibrium, which ties it to the ratio of the components' sound speeds.
Lafarge et al. remove that coupling by giving each component its own ideal-gas
branch — its own `c_k` and its own `p_k_inf` — and mixing them through the phase
field. `src/lbm/equation_of_state.h` implements it, and `calPhaseField` calls it.

That equation of state used to be computed and then **overwritten** by a linear
mixing rule, `p = rho cs2 - (1+phi)/2 p1_inf - (1-phi)/2 p2_inf`, which reads
like a debugging substitution left in place. Because the cases set
`c1 = c2 = cs`, the two forms agree in the bulk of each fluid and differ only
across the interface — which is exactly where the pressure jump is set. The cost
was a Laplace jump that degraded monotonically with the density ratio:

| density ratio | Δp / (σ/R) with the linear mixing rule |
|---|---|
| 1 | 1.007 |
| 2 | 0.986 |
| 5 | 0.922 |
| 10 | 0.832 |
| 20 | 0.720 |
| 50 | −3.07 (sign inverted) |
| 100 | diverged at step 2000 |

At ratio 1 the scheme reproduces Laplace's law to 0.7 %, so σ is calibrated
correctly and Ω⁽²⁾ works; the error tracks the density contrast alone.

## Isotropy of the colour gradient

The colour gradient is differentiated twice per step and divided by its own
norm to give the interface normal (in Ω⁽²⁾, or in the capillary stress of
`laplace`), so gradient anisotropy becomes an anisotropic surface tension.
`src/lbm/isotropic_gradient.h` offers three stencils, selectable as the first
argument of every solver:

| Stencil | Neighbours | Reach | Isotropic to | Angular error on a tanh interface |
|---|---|---|---|---|
| `E4` | 8 | 1 | 4th order | 0.145° |
| `E6` | 12 | 2 | 6th order | 0.038° |
| `E8` | 24 | 2 | 8th order | 0.012° |

All three are exact for a linear field (their normalisation is
`Σ W c_α c_β = δ_αβ`); they differ in the higher lattice tensors. The angular
error is the largest angle between −∇φ and the outward radius over the interface
band of a radial tanh profile, and is measured by
`programs/unit_testing/lbm/gradient`. Going from E4 to E8 reduces it twelvefold,
which is the effect Leclaire et al. report.

`laplace` defaults to E8. The other cases default to E4 because they bounce back
on y: a stencil reaching two nodes becomes one-sided against the wall, which has
not been validated here.

## The high-density-ratio scheme

With the scheme described so far, the Laplace case diverged at a density ratio
of 1000 within its first twenty steps, whatever the viscosities. `laplace` now
carries four changes that take it to $10^4$. Each one lives in `src/lbm` with
its reference and fixes a failure that was measured on this case.

The diagnostic numbers in the subsections come from a parameterised research
copy of this program, which reproduces the old program's output digit for digit
when run with the old algorithm. The "new scheme" figures in the
[results](#results) tables are the output of `laplace` itself.

### The phase field is a mass fraction

$\phi = (\rho_1 - \rho_2)/(\rho_1 + \rho_2)$ is built on the partial densities
carried by $f$ and $g$, so $Y_1 = (1+\phi)/2$ is a mass fraction. The equation of
state is the pressure-equilibrium mixture $1/\rho = Y_1/\rho_1(p) + Y_2/\rho_2(p)$,
and the density at fixed pressure is linear in the *volume* fraction
$\alpha_1 = \rho Y_1/\rho_1(p)$, not in $\phi$. The two fractions differ by a
logistic shift, $\operatorname{logit} Y_1 = \operatorname{logit}\alpha_1 + \ln(\rho_1/\rho_2)$,
so the tanh profile of width $W$ that the recolouring maintains in $\phi$ is the
same profile in $\alpha_1$, displaced by

$$\delta = \tfrac{W}{2}\ln\frac{\rho_1}{\rho_2}.$$

This is the "two interfaces separate" issue that used to be listed here as open.
It is not an artefact: it is where the phase field has to be. With $W = 1.6$,
$\delta$ is 2.39 lattice units at a density ratio of 20 (2.27 was measured
under the old initial state; 2.35 now) and 7.37 at $10^4$. The consequence is
the one that matters: anything meant to act on the interface has to be built on
the volume fraction. On $\phi$ it acts 7.4 nodes out in the light fluid at
$10^4$.

### 1. An initial state in equilibrium with the equation of state

The cases initialised $\rho$ linearly in $\phi$ and the pressure with the linear
mixing rule. Given the mass-fraction nature of $\phi$, that state is far from
the equation of state across the interface. The interface pressure comes out at
6.6 times the ambient one at a density ratio of 20 and 309 times at 1000, and the
run diverges within ten steps. `laplace` now gives the tanh profile to the
volume fraction, steps the pressure by $\sigma/R$ across it, and gives each
component the density its own branch has at that pressure
(`mixture_from_volume_fraction`). The equation of state returns that pressure
at every node. The density interface starts on the prescribed radius and the
$\phi = 0$ contour starts $\delta$ outside it.

### 2. The colour gradient on the normalised colour field

The interface normal is taken on
$\psi = (\alpha_1 - \alpha_2)/(\alpha_1 + \alpha_2)$ (`normalised_phase`). That
is the normalised density $\rho^N$ that Leclaire et al. (2011) and Ba et al.
(2016) differentiate instead of $\phi$, with the reference densities replaced by
$\rho_k(p)$. $\psi$ is centred on the density interface. On the $10^4$ case,
keeping the surface tension on $\nabla\phi$ makes the jump come out at
$\sigma/R_\phi$ with $R_\phi = R + 6.3$, i.e. 0.62 σ/R against the droplet the
density actually describes.

### 3. Surface tension as a body force

$\Omega^{(2)}$ writes the capillary stress into the non-equilibrium populations
scaled by $1/\tau$, and streaming carries that stress between neighbours before
it is relaxed. The two cancel only where $\tau$ is uniform: the leading error
is proportional to $\tau - 1$ times the Laplacian of the capillary stress, and
across an interface where $\tau$ changes it does not integrate away:

| case | Δp / (σ/R) with Ω⁽²⁾ |
|---|---|
| density ratio 1, droplet 100× more viscous | 0.86 |
| density ratio 20, one kinematic viscosity (τ: 5.5 → 100), Ω⁽²⁾ on ∇ψ | 0.71 |
| density ratio $10^4$, surroundings 10× more viscous | 1.03 |

With the force below, the first two cases give 1.03 and 1.02.

The shipped case escaped this only because $\Omega^{(2)}$ acted on $\phi$, 2.4
nodes out in the light fluid, where the density -- and $\tau$ with it -- has
nearly dropped to its outer value. `laplace` applies the same capillary stress
directly, as the divergence of the stress tensor of Lafaurie et al. (1994):

$$\vec F = \nabla\cdot\mathsf T,\qquad \mathsf T = \frac{\sigma}{2}\left(\lvert\nabla\psi\rvert\,\mathsf I - \frac{\nabla\psi\,\nabla\psi}{\lvert\nabla\psi\rvert}\right),$$

added through `S_F` and the half-force velocity correction like gravity
(`src/lbm/surface_force.h`). It does not involve $\tau$. In the continuum it
equals the continuum-surface force $\tfrac{\sigma}{2}\kappa\nabla\psi$ of
Brackbill et al. (1992) and Lishchuk et al. (2003), but only the divergence form
conserves momentum on the lattice: a centred divergence sums to zero over a
periodic lattice, whatever the shape of the interface. The curvature form was
tried first. It gives the same static jumps, but on a droplet translating at
0.01 lattice units per step at density ratio 100 it drained 36 % of the total
momentum in 3000 steps, where the divergence form conserves it to the last
digit. On its own, over a radius-20 circle, the force integrates to 1.007 σ/R,
and its net force on an off-centre ellipse is zero to rounding
(`programs/unit_testing/lbm/mixture`).

### 4. One viscosity per component

With one kinematic viscosity the droplet's dynamic viscosity scales with its
density, and $\tau = \mu/p + 1/2$ with it: τ ≈ 5×10⁴ in the droplet at $10^4$.
There the jump drops to about 0.1 σ/R and the currents rise to 2×10⁻².
`laplace` gives each component its own $\nu_k$ and mixes them as
$\nu = Y_1\nu_1 + Y_2\nu_2$, i.e. $\rho\nu = \alpha_1\mu_1 + \alpha_2\mu_2$
(`mixture_kinematic_viscosity`). With $\nu_1 = \nu_2$ that is the old rule
exactly. The density and viscosity ratios are the second and third arguments of
the program.

### Results

On the shipped case (128², R = 10, density ratio 20, one kinematic viscosity,
30 000 steps), scored against the radius the density settles at:

| | Δp / (σ/R(ρ)) | R(ρ) | R(φ = 0) − R(ρ) | max &#124;u&#124; |
|---|---|---|---|---|
| linear mixing, E4 | 0.720 | 10.06 | 0.00 | 1.44 × 10⁻³ |
| restored EOS, E4 | 0.962 | 9.42 | 2.27 | 1.22 × 10⁻³ |
| restored EOS, E8 | 0.962 | 9.42 | 2.27 | 1.11 × 10⁻³ |
| high-density-ratio scheme, E4 | 1.017 | 9.91 | 2.35 | 1.77 × 10⁻⁵ |
| high-density-ratio scheme, E8 | 1.021 | 9.91 | 2.35 | 1.35 × 10⁻⁵ |

The parasitic currents fall by a factor of 80, and the state is stationary to
the six digits written. The 2 % left on the jump is the finite interface width
at R = 10, not the density ratio: the same case at R = 20 gives 1.004 (research
copy).

Against the density ratio, with the same dynamic viscosity in both fluids
(`laplace E8 <ratio> 1`, 30 000 steps):

| ρ₁/ρ₂ | before: Δp / (σ/R) | now: Δp / (σ/R) | now: R(ρ) | now: R(φ = 0) − R(ρ) | now: max &#124;u&#124; |
|---|---|---|---|---|---|
| 1 | 1.005 | 1.023 | 9.90 | — | 1.78 × 10⁻⁵ |
| 20 | 0.816 | 1.021 | 9.91 | 2.35 | 1.86 × 10⁻⁵ |
| 100 | diverged before step 10 000 | 1.022 | 9.94 | 3.47 | 1.87 × 10⁻⁵ |
| 1000 | diverged within 20 steps | 1.023 | 9.97 | 5.02 | 2.05 × 10⁻⁵ |
| 10⁴ | diverged within 20 steps | 1.015 | 9.93 | 6.52 | 7.18 × 10⁻⁴ |

"Before" is the previous scheme given the same two viscosities (research copy);
with its own single kinematic viscosity it was still unsteady at step 6000 at
100 (0.91 σ/R there), and diverged within 30 steps at 1000. At ratio 1 the
density is uniform and the radius is the phase interface's.

Up to a ratio of 100 the runs are stationary over their last 20 000 steps to the
digits written. At 1000 and $10^4$ the jump still moves by 0.16 % and 0.3 % of
σ/R over that span, and the radius by 0.01 and 0.03. The jump no longer depends
on the density ratio; the currents rise only at $10^4$. The offset
$R(\phi = 0) - R(\rho)$ follows $\tfrac{W}{2}\ln(\rho_1/\rho_2)$ (2.39, 3.68,
5.53, 7.37) but falls increasingly short as the ratio grows, by 2 % at 20 and
12 % at $10^4$, as the $\phi = 0$ contour moves into the far tail of the
profile.

### What is not solved

- **Moving interfaces at high density ratio.** A droplet translating through a
  periodic box with velocity $U$ (the whole domain initialised at $U$, equal
  viscosities, 3000 steps, research copy) stays bounded only at small $U$ as
  the ratio grows, and is accurate only at moderate ratios:

  | ρ₁/ρ₂ | U | outcome | droplet speed | max &#124;u&#124; |
  |---|---|---|---|---|
  | 100 | 10⁻² | bounded | 0.96 U | 2.6 × 10⁻² |
  | 1000 | 10⁻³ | bounded | 1.00 U | 1.6 × 10⁻² |
  | 1000 | 10⁻² | diverges | — | — |
  | $10^4$ | 10⁻⁴ | bounded | 1.03 U | 3.2 × 10⁻² |
  | $10^4$ | 10⁻³ | diverges | — | — |

  Total momentum is conserved to the last digit in every bounded run, but at
  $10^4$ the spurious currents are 300 times the droplet's speed.

  The mechanism is in the colour transport, not in the surface tension; the
  perturbation operator fails at the same velocities. In the heavy fluid the
  moving populations are $\approx p/3 \pm \rho u/2$, strongly negative once
  $\rho u > p$. Colour streamed as $\phi f_i$ with negative weights leaves the
  convex hull, and $\phi$ overshoots 1 at the advancing front. There the
  equation of state is badly conditioned: in the heavy half of the interface the
  pressure depends on $1-\phi \sim \rho_2/\rho_1$, so an error of $10^{-3}$ in
  $\phi$ turns the pressure negative, and $\tau$ with it. A colour transport
  that is not tied to the populations of $f$ is the natural next step, e.g. a
  separate phase-field distribution with its own equilibrium, or the
  velocity-based equilibrium of Subhedar (2022). Flooring the pressure keeps a
  1-D front alive but not a droplet.
- **A viscous heavy fluid at high density ratio.** Keep τ in the droplet of
  order 1–10. At $10^4$ with $\mu_1/\mu_2 = 20$ (τ ≈ 100 in the droplet), the
  jump wanders between 0.90 and 1.01 σ/R over 30 000 steps and the currents
  reach 2×10⁻² (research copy).
- **Only `laplace` carries the scheme.** `capillary`, `gravity_capillary` and
  the three `rayleigh_taylor` programs still use $\Omega^{(2)}$ on $\nabla\phi$,
  the density linear in $\phi$ and one viscosity. They run at a density ratio
  of 4, where the difference is small. Porting them is mechanical with
  `src/lbm/mixture.h` and `src/lbm/surface_force.h` (the force takes
  `Boundary::WallY`), but it changes their results and none of them has a
  validation test yet.

## Status of the modular library

`src/core/constants.h`, `src/lbm/lattice_boltzmann.h` and `src/main_cglbm.cpp`
are an in-progress refactor that factors the shared algorithm out of the
programs. The implementation units behind `lattice_boltzmann.h` do not exist
yet, so the `cglbm` library is header-only and `src/main_cglbm.cpp` is not built
into a target — it would compile but not link. The programs under
`programs/solvers` are the working code.
