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

All fields are members of `cglbm::lbm::Solver` and hold `cglbm::lbm::Field`,
one heap allocation laid out as `[nx][ny][depth]` — the memory order the old
`double f[Lx][Ly][Q]` arrays had, without fixing the lattice size at compile
time.

| Symbol | Member | Meaning |
|---|---|---|
| $f_i$ | `f_` | sum (total) distribution function |
| $g_i$ | `g_` | difference (color) distribution function |
| $f_i^{eq}$ | `f_eq_` | equilibrium distribution |
| $\rho$ | `rho_`, `rho_mdt_` | mixture density, current and previous step |
| $p$ | `p_`, `p_mdt_` | pressure, current and previous step |
| $\vec u$ | `u_` | velocity, two components per node |
| $\phi$ | `phi_` | phase field, $\phi = +1$ in fluid 1, $-1$ in fluid 2 |
| $\Omega^{(1,2,3)}$ | `omega_1_/2_/3_` | the three collision contributions |
| $S_i$, $\vec F$ | `source_`, `force_` | forcing term and external body force |

$S^F$, $S^{Sp}$ and $S^t$ are not fields: each is built and consumed within one
node of `force()`, so they are local arrays of nine doubles. As lattice-sized
arrays they cost about a gigabyte on the production Rayleigh-Taylor case.

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

The forcing term assembled by `Solver::force()` is the sum of three parts:

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
`Solver::phase_field()` by solving the quadratic mixing relation, implemented in
`src/lbm/equation_of_state.cpp`.

For the Laplace benchmark the initial `p1_inf` is offset by $-\sigma/R$ so that
the initial condition already satisfies the pressure jump.

## Time loop

`Solver::step()` executes, per step:

```
force()            →  S_i from the external body force F
collide()          →  Ω(1)  relaxation towards f_eq
collide_surface()  →  Ω(2)  surface-tension perturbation
recolor()          →  Ω(3)  interface sharpening
stream()           →  propagation along ξ_i
macroscopic()      →  ρ, u from moments of f (force-corrected)
phase_field()      →  φ and p from the equation of state
equilibrium()      →  f_eq for the next step
```

`Solver::run()` initialises, then loops `config.steps` times, writing the CSV
grids every `config.interval` steps through `cglbm::lbm::CsvWriter`.
`cglbm::lbm::write_vtk` produces a ParaView/VisIt file from the same state.

## Boundary conditions

`CaseConfig::boundary` selects one of two policies, applied consistently by the
colour gradient, by `force()` and by `stream()`:

| Program | `boundary` | x | y |
|---|---|---|---|
| `laplace` | `PeriodicY` | periodic | periodic |
| `capillary`, `gravity_capillary` | `WallY` | periodic | half-way bounce-back |
| `rayleigh_taylor`, `..._omp` | `WallY` | periodic | half-way bounce-back |

Under `WallY` the colour gradient goes through `gradient_wall_y`, which drops
any neighbour beyond a wall and so becomes one-sided within `stencil_reach`
nodes of it; `force()` drops the same neighbours from its lattice divergence.

Half-way bounce-back is implemented inline: at `j == 0` the populations moving
downwards (`k = 4, 7, 8`) stay on the node and are reflected into `k-2`, and at
`j == Ly-1` the upward ones (`k = 2, 5, 6`) are reflected into `k+2`. This
places the resting wall halfway between the last fluid node and the ghost node.

An `applyAbsorbingBoundary()` that damped acoustic waves on the outermost nodes
was defined in every solver and commented out of every time loop. It was removed
with the rest of the duplication; recover it from the history if the initial
pressure transient ever needs damping.

## Case-specific validation

| Case | Analytic reference |
|---|---|
| `laplace` | $\Delta p = \sigma / R$ across the droplet interface (see the open point below) |
| `capillary` | $T_{theo} = 2\pi\sqrt{(\rho_1+\rho_2)\,r^3/(6\sigma)}$, printed at startup and compared against `interface.csv` |
| `gravity_capillary` | Balance of the Laplace jump against the hydrostatic head |
| `rayleigh_taylor` | Linear growth rate of the instability, $\sigma = 0$ |

## The equation of state and the density ratio

The colour-gradient method is classically limited to modest density ratios
because the density contrast is carried by the rest-particle weight of the
equilibrium, which ties it to the ratio of the components' sound speeds.
Lafarge et al. remove that coupling by giving each component its own ideal-gas
branch — its own `c_k` and its own `p_k_inf` — and mixing them through the phase
field. `src/lbm/equation_of_state.h` implements it, and `Solver::phase_field()` calls it.

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

The colour gradient is differentiated twice per step, and Ω⁽²⁾ divides it by its
own norm, so gradient anisotropy becomes an anisotropic surface tension.
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

## Where the Laplace case stands now

With the equation of state restored and the E8 gradient, on the shipped case
(128², prescribed R = 10, density ratio 20):

| | Δp | R(ρ) | Δp / (σ/R(ρ)) | max &#124;u&#124; |
|---|---|---|---|---|
| linear mixing, E4 | 0.019806 | 10.06 | 0.720 | 1.44 × 10⁻³ |
| restored EOS, E4 | 0.028662 | 9.42 | 0.962 | 1.22 × 10⁻³ |
| restored EOS, E8 | 0.028294 | 9.42 | 0.962 | 1.11 × 10⁻³ |

The jump error falls from 28 % to under 4 %, and the parasitic currents by 23 %.

### Open issue: the two interfaces separate

The radius above is the one the **density** field settles at. It has to be
stated, because the restored equation of state introduces an artefact that the
linear mixing rule did not have: the phase field and the density field come to
rest at different radii.

| | R(φ = 0) | R(ρ midpoint) | gap |
|---|---|---|---|
| linear mixing | 10.06 | 10.07 | 0.00 |
| restored EOS | 11.69 | 9.42 | 2.27 |

Both start at 10.0, separate within the first ~3000 steps, and then hold to six
digits — it is a converged state, not a drift, and the total mass of each colour
is conserved exactly throughout. The consequence is that "the radius of the
droplet" is ambiguous to about two lattice units, and the Laplace ratio reads
0.962 against the density interface and 1.194 against the phase interface.

`programs/solvers/color_gradient/laplace/tests/` pins the gap so that a scheme
change closing it is visible. Closing it is the natural next piece of work, and
the equilibrium is the place to look: the enhanced equilibria of Leclaire et al.
(2013) and the third-order Hermite equilibrium of Ba et al. (2016) both target
the interfacial momentum error that a truncated equilibrium leaves behind at
high density contrast. See [`references.md`](references.md).

## Structure of the solver

The scheme has one implementation, `cglbm::lbm::Solver` in `src/lbm/solver.cpp`.
A program under `programs/solvers` is a *case definition*: it fills a
`cglbm::lbm::CaseConfig` — lattice, `Physics`, boundary, initial phase field,
stencil — and hands it to the solver.

| Header | Holds |
|---|---|
| `src/lbm/d2q9.h` | the velocity set, the weights, $c_s$ and its powers |
| `src/lbm/field.h` | `Field`, the heap-allocated lattice storage |
| `src/lbm/case_config.h` | `CaseConfig`, `Physics`, the boundary and initial-condition choices, the command line |
| `src/lbm/solver.h` | the scheme itself |
| `src/lbm/output_writer.h` | the CSV and VTK writers |
| `src/lbm/equation_of_state.h` | the two-component pressure |
| `src/lbm/isotropic_gradient.h` | the E4/E6/E8 colour-gradient stencils |

Before this, the five solvers were 566 to 610 line programs holding a copy each
of the same scheme; `capillary` and `gravity_capillary` differed on seven lines
out of 588. The case definitions are now 76 to 94 lines.

Every parameter is a runtime value. The lattice size, the step count and the
output interval are `--nx`, `--ny`, `--steps` and `--interval`, so a resolution
study no longer means editing a `const int` and rebuilding, and a test can ask a
run what it was configured with instead of repeating the constant.

### Reproducibility

`cmake/CompilerFlags.cmake` passes `-ffp-contract=off`. Fusing `a*b+c` into an
FMA changes the result, and the compiler does it inconsistently: without the
flag the same source gives different last bits at `-O2` and at `-O3`, and moving
an expression into a function can move the answer. With it, a run is reproducible
across optimisation levels, and a refactor can be checked against a recorded run
bit for bit — which is how this one was checked, over 82 output files of four
cases.

## Known defects

Two defects in the wall-bounded cases (`capillary`, `gravity_capillary`,
`rayleigh_taylor`, `rayleigh_taylor_omp`). `laplace` is periodic and is affected
by neither.

### The colour distribution is rebuilt from the wrong population at a wall

`Solver::stream()` writes the reflected population to `f_(ip, jp, kp)` and then
reads `f_(ip, jp, k)` to build `g`:

```cpp
f_(ip, jp, kp) = f_eq_(i, j, k) + omega_1_(i, j, k) + omega_2_(i, j, k) + 0.5 * source_(i, j, k);
g_(ip, jp, kp) = f_(ip, jp, k) * phi_(i, j) + omega_3_(i, j, k);
```

On an interior node `kp == k` and the two agree. On a bounce-back direction —
`j == 0` with `k` in {4, 7, 8}, `j == ny-1` with `k` in {2, 5, 6} — they do not,
so `g` is built from an unrelated direction of `f`. Since $\phi = \sum g / \sum f$,
the phase field leaves $[-1, 1]$: `programs/unit_testing/lbm/solver` measures a
peak of 1.29 against a wall, where the periodic case stays within $10^{-13}$ of 1.

### The source term is never computed on the `j = 0` row

`Solver::force()` runs `for (int j = j_start; j < ny_; j++)` with `j_start == 1`
under `WallY`. Row `j = 0` therefore keeps a source term of zero for the whole
run, while row `j = ny-1` gets a full one — the two walls are not treated alike.
It also makes the `j == 0` branch inside the neighbour loop unreachable.

Both are preserved as they were so that the extraction of the shared kernel
could be verified bit for bit against the previous code. They are fixed in the
commits that follow this one, each on its own, with the change in the results
reported.
