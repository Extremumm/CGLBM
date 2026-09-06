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
| $\phi$ | `phi[Lx][Ly]` | phase field, $\phi = +1$ in fluid 1, $-1$ in fluid 2 |
| $\Omega^{(1,2,3)}$ | `omega_1/2/3` | the three collision contributions |
| $S_i$, $\vec F$ | `S`, `F` | forcing term and external body force |

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
   pressure from the equation of state, not on $\rho c_s^2$.
2. **$\Omega^{(2)}$ — surface tension** (`collide_surface`). A perturbation
   built on the color gradient $\nabla\phi$ that introduces the surface tension
   $\sigma$ and reproduces the Laplace jump.
3. **$\Omega^{(3)}$ — recoloring** (`recolor`). Redistributes the components
   along the color gradient to keep the interface sharp. `epsilon` guards the
   division by the gradient norm; `ch_width_ope` sets the operating interface
   thickness (`ch_width_init` is used only when initialising).

The forcing term assembled by `force()` is the sum of three parts:

$$S_i = S_i^{F} + S_i^{Sp} + S_i^{t}$$

- `S_F` — the external body force (gravity `a_g`, zero in the non-gravity cases);
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

## Time loop

`runSimulation()` executes, per step:

```
force()            →  S_i from the external body force F
collide()          →  Ω(1)  relaxation towards f_eq
collide_surface()  →  Ω(2)  surface-tension perturbation
recolor()          →  Ω(3)  interface sharpening
stream()           →  propagation along ξ_i
calMacroscopic()   →  ρ, u from moments of f (force-corrected)
calPhaseField()    →  φ and p from the equation of state
calEquilibrium()   →  f_eq for the next step
```

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
| `laplace` | $\Delta p = \sigma / R$ across the droplet interface (see the open point below) |
| `capillary` | $T_{theo} = 2\pi\sqrt{(\rho_1+\rho_2)\,r^3/(6\sigma)}$, printed at startup and compared against `interface.csv` |
| `gravity_capillary` | Balance of the Laplace jump against the hydrostatic head |
| `rayleigh_taylor` | Linear growth rate of the instability, $\sigma = 0$ |

## Open point: the Laplace pressure jump

The `laplace` program does **not** reproduce the Laplace law at steady state.
Measured with `pycglbm` on the shipped configuration (128², R = 10, ρ₁/ρ₂ = 20,
σ = 0.2768 in lattice units):

| Timestep | Δp | Δp / (σ/R) | droplet radius |
|---|---|---|---|
| 0 | 0.027684 | 1.0000 | 9.85 |
| 1000 | 0.021959 | 0.7932 | 10.17 |
| 10000 | 0.019806 | 0.7154 | 10.17 |
| 30000 | 0.019806 | 0.7154 | 10.17 |

The jump is exact at `t = 0` because `p1_inf` carries the `- sigma / radius`
offset, so the equation of state and the initial condition are consistent. The
solver then relaxes to a **stationary** state — the jump is constant to six
digits from step 10000 on — holding only 72 % of the analytic jump, while the
droplet radius, the interface profile and the phase bounds stay clean. So this
is not a convergence, dissolution or sharpness problem: the surface tension
effectively realised by $\Omega^{(2)}$ is smaller than the prescribed $\sigma$
for this parameter set.

The parasitic currents at the interface are small and stationary
(max |u| = 1.44 × 10⁻³, i.e. Ma ≈ 2.5 × 10⁻³, unchanged between step 10000 and
30000), so the deficit is not being sustained by an unresolved flow either.

Worth checking, in order: the calibration between $\sigma$ and the amplitude of
the $\Omega^{(2)}$ perturbation (a factor involving `ch_width_ope` or the
relaxation time is the usual suspect), then the finite interface thickness,
which makes the effective radius ambiguous by about one lattice unit — too small
an effect to explain 28 %.

`programs/solvers/color_gradient/laplace/tests/test_laplace_color_gradient.py`
pins the measured ratio so that a change of behaviour is caught; it deliberately
does not assert agreement with $\sigma/R$.

## Status of the modular library

`src/core/constants.h`, `src/lbm/lattice_boltzmann.h` and `src/main_cglbm.cpp`
are an in-progress refactor that factors the shared algorithm out of the
programs. The implementation units behind `lattice_boltzmann.h` do not exist
yet, so the `cglbm` library is header-only and `src/main_cglbm.cpp` is not built
into a target — it would compile but not link. The programs under
`programs/solvers` are the working code.
