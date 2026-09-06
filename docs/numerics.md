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

## How far the density ratio goes

The colour-gradient family is classically limited to modest density ratios, and
this code has been pushed against that limit deliberately. Where it stands, on
the shipped Laplace case (128², prescribed R = 10, E8 gradient, 10⁴ steps), is:

| ratio | before | after | Δp / (σ/R) after | max &#124;u&#124; after |
|---|---|---|---|---|
| 2 | ok | ok | 0.954 | 1.2 × 10⁻³ |
| 10 | ok | ok | 0.894 | 1.2 × 10⁻³ |
| 20 | ok | ok | 0.895 | 1.2 × 10⁻³ |
| 100 | ok | ok | 0.898 | 1.2 × 10⁻³ |
| 500 | **diverged, step < 200** | ok | 0.820 | 2.4 × 10⁻³ |
| 10³ | **diverged, step < 200** | ok | 0.741 | 4.4 × 10⁻³ |
| 10⁴ | **diverged, step < 200** | ok | 0.580 | 9.7 × 10⁻³ |
| 10⁵ | **diverged, step < 200** | ok | 0.515 | 1.2 × 10⁻² |

Read that as two separate claims. **Stability** now reaches 10⁵: the case runs
and the phase field stays inside [−1, 1] to the last digit. **Accuracy** is a
different matter — the pressure jump holds to about 10 % up to a ratio of 100,
then degrades to roughly a factor of two at 10⁵. Nothing here certifies
Laplace's law at 10⁵; it says the scheme survives there and by how much it is
wrong.

### What was breaking

Not the collision operator, and not the colour gradient. The initial condition.

`initialize()` laid the density down as a linear interpolation in φ and the
pressure from a linear mixing rule. The time loop, however, recovers the
pressure from the two-component equation of state, and a density linear in φ is
not in mechanical equilibrium under that equation of state. Two things followed:

1. The initial pressure and the pressure the solver computes at step 1 disagreed
   across the interface by a factor growing with the density ratio — 3.5 at a
   ratio of 10, 32 at 100, 316 at 1000.
2. Even with that fixed, a linearly interpolated density makes the equation of
   state read about 0.31 ρ₁c_s² at the interface, against bulk pressures of
   order c_s². At a ratio of 100 that is a peak pressure of 10.4 against bulks
   of 0.33 and 0.36.

The scheme *does* relax that away — the measured peak falls from 10.44 to 0.43
over 3000 steps at a ratio of 100 — but the relaxation is an acoustic blast. The
velocity reached 0.80 in lattice units by step 4, against a sound speed of
0.577: **Mach 1.4**. Below a ratio of about 100 the scheme rides it out; above,
it does not.

### The fix

`InitialState::MechanicalEquilibrium` solves the equation of state for the
density that makes the pressure smooth across the interface, instead of reading
off whatever pressure a linear density implies. For each node, given φ, it
bisects for the ρ satisfying

    p(ρ, φ) = p_light + (p_heavy − p_light) (1 + φ)/2

Both bulks come out exactly ρ₁ and ρ₂, so the Laplace jump at t = 0 is still
exactly σ/R. Only the interface profile changes — and it starts where the scheme
was going to take it anyway. At a ratio of 100 the initial transient drops from
Mach 1.4 to Mach 0.035, a factor of 40, and the pressure never collapses:

| | linear start | equilibrium start |
|---|---|---|
| peak &#124;u&#124; over the first 30 steps | 0.80 | 0.020 |
| min pressure at step 30 | −2.6 × 10⁻⁴ | 0.28 |
| max τ_ν at step 30 | 4.8 × 10⁸ | 5.9 × 10² |

`--initial-state=linear` reproduces the old behaviour, and
`--initial-state=eos` is the intermediate step: linear density, but the pressure
the solver will actually compute.

The cost, on the shipped Laplace case at ratio 20, is Δp/(σ/R) falling from
0.962 to 0.895, and the density interface starting at R = 8.41 rather than at
the prescribed 10 — mechanical equilibrium, not a linear interpolation, is what
sets the density profile, and it is not symmetric about φ = 0. The gap between
the phase and density interfaces is essentially untouched (2.270 → 2.284), which
is more evidence that it belongs to the equation of state rather than to the
initialisation.

### What did not work

Ba et al. (2016) identify the phase field used to locate the interface as the
thing that fails at high density ratio: φ = (ρ₁ − ρ₂)/ρ has its zero contour at
the interface only when the two bulk densities are equal. In general φ = 0 sits
at (ρ₁ − ρ₂)/(ρ₁ + ρ₂), which is 0.90 at a ratio of 20 and 0.9998 at 10⁴ — well
inside the light fluid. Their Eq. (21) normalises each component by its own bulk
density first, and is available here as
`InterfaceField::BulkNormalised` (`--interface-field=normalised`).

Measured on this case, it is worse:

| ratio | 2 | 5 | 10 | 20 | 50 | 100 |
|---|---|---|---|---|---|---|
| `colour` | 0.956 | 0.925 | 0.933 | 0.962 | 1.016 | 1.013 |
| `normalised` | 0.983 | 0.905 | 0.809 | 0.696 | 0.554 | 0.445 |

It also leaves the interface gap where it was (3.310 → 3.342 at a ratio of 100)
and does not change the ratio at which the old initialisation diverged. The
reason is that Ba et al. use the normalised field for the interface *normal*, in
a continuum-surface-force operator where σ appears explicitly. Ω⁽²⁾ here, from
Lafarge et al., carries its calibration in the *magnitude* of the gradient, so
substituting a field with a different profile changes the surface tension it
produces. Transferring Eq. (21) means recalibrating Ω⁽²⁾ against it. The option
is kept because it is what the literature specifies and is where that work would
start.

### Where the literature is

No colour-gradient model in the literature is validated at 10⁵.

| Work | Reached | How |
|---|---|---|
| Leclaire et al. 2011 | O(10⁴), static Laplace, 0.5 % error | isotropic colour gradient + Latva-Kokko recolouring |
| Leclaire et al. 2013 | 10³, dynamic | enhanced equilibrium distributions |
| Ba et al. 2016 | 10³, dynamic, high Re | MRT collision, CSF perturbation, normalised phase field |
| Saito et al. 2023 | 10 accurately, 10³ marginally | sixth-order Hermite equilibria, central moments |

Ratios of 10⁵ and beyond are reported by other families — chemical-potential
pseudopotential models (> 6.5 × 10⁴), phase-field Allen–Cahn models, entropic
MRT — not by colour-gradient ones. The accuracy figures in the first table are
consistent with that: quantitative to O(10²–10³), qualitative beyond.

The next ingredients to try, in the order the literature suggests, are the
enhanced equilibrium of Leclaire et al. (2013), then an MRT collision. Both
target the same thing: the interfacial momentum error a truncated equilibrium
leaves behind at high density contrast — which is also the most likely cause of
the interface gap documented above.

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
| the above, equilibrium start | 0.029726 | 8.33 | 0.895 | 1.19 × 10⁻³ |

The jump error fell from 28 % to under 4 % with the equation of state restored.
The last row is the current default: it gives back 7 % on this case in exchange
for a scheme that runs to a density ratio of 10⁵ rather than diverging above
100 — see [How far the density ratio goes](#how-far-the-density-ratio-goes).
`--initial-state=linear` reproduces the row above it.

### Open issue: the two interfaces separate

The radius above is the one the **density** field settles at. It has to be
stated, because the restored equation of state introduces an artefact that the
linear mixing rule did not have: the phase field and the density field come to
rest at different radii.

| | R(φ = 0) | R(ρ midpoint) | gap |
|---|---|---|---|
| linear mixing | 10.06 | 10.07 | 0.00 |
| restored EOS | 11.69 | 9.42 | 2.27 |
| restored EOS, equilibrium start | 10.62 | 8.33 | 2.28 |

The initialisation barely moves it, which is part of why it is attributed to the
equation of state. Under the equilibrium start the gap is already 1.59 at t = 0,
because mechanical equilibrium is what sets the density profile and it is not
symmetric about φ = 0.

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

### Fixed: the normalised phase field was never refreshed after t = 0

`Solver::update_interface_field()` builds `phi_n_`, the bulk-normalised phase
field of Ba et al. Eq. (21). It was called once, at the end of `initialize()`,
and never again — so `--interface-field=normalised` took its gradient of the
field as it stood at t = 0, for the whole run.

The option was therefore never actually tested, and the negative result recorded
against it in this file was not evidence about the normalised field. Re-measured
with the refresh in place, the option on its own moves the Laplace jump by 0.2 %
of itself. `update_interface_field()` now runs in `step()`, after
`phase_field()`, from the phase field that step just computed.

### Fixed: the colour distribution was rebuilt from the wrong population at a wall

`Solver::stream()` wrote the reflected population to `f_(ip, jp, kp)` and then
read `f_(ip, jp, k)` to build `g`:

```cpp
f_(ip, jp, kp) = f_eq_(i, j, k) + omega_1_(i, j, k) + omega_2_(i, j, k) + 0.5 * source_(i, j, k);
g_(ip, jp, kp) = f_(ip, jp, k) * phi_(i, j) + omega_3_(i, j, k);   // k, not kp
```

On an interior node `kp == k` and the two agree. On a bounce-back direction —
`j == 0` with `k` in {4, 7, 8}, `j == ny-1` with `k` in {2, 5, 6} — they do not,
so `g` was built from an unrelated direction of `f`. Since
$\phi = \sum g / \sum f$, the phase field left $[-1, 1]$:
`programs/unit_testing/lbm/solver` measured a peak of **1.29** against a wall,
where the periodic case stays within $10^{-13}$ of 1. With `kp` on both lines
the wall case matches the periodic one to the same $10^{-13}$.

The fix changes the wall-bounded cases. After 200 steps of the shipped
configurations, relative to the peak of each field:

| Case | ρ | φ | p | u |
|---|---|---|---|---|
| `capillary` | 6.5 % | 4.1 % | 4.8 % | 14.1 % |
| `gravity_capillary` | 11.1 % | 6.0 % | 10.0 % | 11.5 % |
| `rayleigh_taylor` | 1.3 % | 3.4 % | 1.2 % | 5.4 % |

`laplace` is periodic and its output is unchanged, byte for byte.

### Fixed: the source term was never computed on the `j = 0` row

`Solver::force()` ran `for (int j = j_start; j < ny_; j++)` with `j_start == 1`
under `WallY`. Row `j = 0` therefore kept a source term of zero for the whole
run — no body force, no third-order correction, no temporal correction — while
row `j = ny-1` got a full one. The two walls were not treated alike. It also
made the `j == 0` branch inside the neighbour loop unreachable, though that
branch is exactly what makes starting at `j = 0` safe: it drops the directions
that would read below the wall.

The loop now starts at 0 on both boundaries. The effect is much smaller than
the streaming defect above, because in three of the four cases the interface
starts far from the wall and the bottom row is quiescent early on. After 200
steps, relative to the peak of each field:

| Case | ρ | φ | p | u |
|---|---|---|---|---|
| `capillary` | 0.09 % | 0.03 % | 0.32 % | 0.41 % |
| `gravity_capillary` | 0 | 0 | 0 | 1e-8 |
| `rayleigh_taylor` | 2e-7 | 1e-8 | 2e-6 | 1e-6 |

The figures below 1e-6 are at the resolution of the six-digit CSV output rather
than a measurement of the change. `laplace` is periodic and unchanged.

Both were preserved through the extraction of the shared kernel, so that it
could be verified bit for bit against the previous code, and fixed afterwards
one commit at a time. Any measurement taken before those commits carries them.
