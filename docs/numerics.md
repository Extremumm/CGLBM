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
colour gradient    →  ∇φ, evaluated once and reused by the three operators
surface_force()    →  F_s, the capillary body force, when the case asks for it
equilibrium()      →  f_eq for the next step
```

The colour gradient is stored rather than recomputed: `collide_surface()`,
`recolor()` and `surface_force()` all want the gradient of the phase field as it
stood at the end of the previous step, and `phase_field()` is the only thing
that changes it. Evaluating it once is the same number to the last bit and, with
the E8 stencil's 24 neighbours, about a fifth of the run.

A case chooses between two ways of applying the surface tension, and only one
of `collide_surface()` and `surface_force()` does anything in a given run.
`SurfaceTension::Perturbation` injects a capillary stress into Ω(2);
`SurfaceTension::ContinuumSurfaceForce` forms the curvature explicitly and adds
a body force, which is what `force()` and `macroscopic()` then pick up along
with gravity. See [How far the density ratio goes](#how-far-the-density-ratio-goes).

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
this code has been pushed against that limit deliberately. Two rounds of work
are recorded below: getting the scheme to *survive* past a ratio of about 100,
which was an initial-condition problem, and getting it to be *right* there,
which took the two papers this section is named for.

**How long the run is matters more than anything else here, and getting that
wrong is how this file came to contain claims it should not have.** A short run
at a high density ratio does not diverge; it has simply not got there yet. Every
figure below is from 1.2 × 10⁵ steps on the shipped Laplace case (128²,
prescribed R = 10, E8 gradient, ν = 1.664), reported as Δp / (σ/R) against the
prescribed radius, with the whole time series checked rather than the last
value:

| ratio | Ω⁽²⁾ stress, colour φ | φ_N interface + CSF tension | steady? |
|---|---|---|---|
| 20 | 1.140 | **1.021** | yes, to 6 s.f. |
| 100 | 1.205 | **1.023** | yes, constant from 4.8 × 10⁴ steps |
| 200 | — | **1.021** | yes |
| 500 | — | **0.995** | still rising slowly, no divergence |
| 10³ | **diverges, ≈ 6 × 10⁴** | **diverges, 8.7 × 10⁴** | no |
| 10⁴ | diverges | diverges | no |
| 10⁵ | diverges | diverges | no |

So, plainly:

- **to a density ratio of about 500 the scheme is quantitative and steady** — a
  couple of per cent on Laplace's law, converged and holding;
- **at 10³ and above it diverges.** Not slowly, and not only in one
  configuration: the continuum-surface-force operator dies at step 8.7 × 10⁴ and
  the older stress operator by 6 × 10⁴.

The 10³ run *looks* healthy for its first 5 × 10⁴ steps — it reads 0.96, then
0.89, then 0.85, then collapses. Measured at 3 × 10⁴ steps it reports −7.6 %,
which is why this file previously claimed a density ratio of 10³ and, before
that, of 10⁵. Both claims were runs that had not been given time to fail. The
shipped Laplace case at a ratio of 20 is genuinely converged, and its steady-state
test checks exactly that; nothing above a few hundred had ever been given the
same check.

The trade between the two operators is unchanged and still worth making: at 100
the CSF column reads 1.023 against 1.205, and it cut the spurious currents by a
factor of 72.

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

### Reading Ba et al. (2016) and Leclaire et al. (2013)

Both papers target exactly this problem, and between them they name four
ingredients. Two of them turned out to be here already, and saying which is
half the work — the equilibrium used here comes from a Hermite expansion and
never mentions their parameter α, so nothing about it looks like their
equations until the algebra is done.

**Already present: the enhanced equilibrium.** Leclaire et al. (2013) add a
term to the colour-gradient equilibrium to repair its third-order velocity
moment, which a two-component lattice gets wrong whenever the components'
sound speeds differ; Ba et al. restate it as their Eq. (14) and build on it.
Writing their free parameter through $(c_s^k)^2 = \tfrac{3}{5}(1 - \alpha_k)$
and $p_k = \rho_k (c_s^k)^2$, their extra term is

$$\rho_k W_i\,(3\,\mathbf e_i\!\cdot\!\mathbf u)\,\tfrac12\bigl(3 (c_s^k)^2 - 1\bigr)\bigl(3|\mathbf e_i|^2 - 4\bigr)
 = (p_k - \rho_k c_s^2)\,W_i (\mathbf e_i\!\cdot\!\mathbf u)\,4.5\,\bigl(3|\mathbf e_i|^2 - 4\bigr)$$

and the third-order Hermite term of `Solver::equilibrium` reduces to the same
expression, because $H_{xxx} + H_{yyx} = \xi_x(|\xi|^2 - 4c_s^2)$ and
$1/(2c_s^6) = 13.5$. They agree to $1.1 \times 10^{-14}$ over a sweep of sound
speed, density and velocity; `report_enhanced_equilibrium` in
`programs/unit_testing/lbm/solver` pins it so that an edit which quietly drops
the correction fails a fast test.

**Already present: the Chapman–Enskog source correction.** Ba et al. Eqs.
(17)–(18) add $C_1 \propto \partial_x Q_x + \partial_y Q_y$ to the energy mode
and $C_7 \propto \partial_x Q_x - \partial_y Q_y$ to the diagonal stress, with
$Q_\alpha = (1.8\alpha_k - 0.8)\rho_k u_\alpha$. That prefactor is
$1 - 3(c_s^k)^2$, so $Q_\alpha = -3(p_k - \rho_k c_s^2)u_\alpha$ — and
$(p - \rho c_s^2)\mathbf u$ is precisely the quantity `S_Sp` in
`Solver::force()` differentiates, on the same two moments, with the same
nine-point isotropic stencil (their Eq. 20 is $\sum_i W_i \psi\, e_{i\alpha} / c_s^2$).
The $(1 - s/2)$ prefactors appear here as the Guo half-step: `f - f_eq + S/2` in
`collide()` together with `+ S/2` in `stream()` gives the source its
$1 - 1/(2\tau)$.

**Not present, and needed: locating the interface.** This is Ba et al. Eq. (21),
and it is the subject of the next section.

**Not present, and needed: how the tension is applied.** Ba et al. Eqs.
(23)–(29), the section after that.

**Not present, and now added: two viscosities.** Ba et al. Eq. (22) interpolates
the relaxation parameter across the interface because their $s_\nu$ is defined
per fluid. This code had no per-fluid viscosity at all: `Physics::nu` was one
number for the mixture, so the viscosity *ratio* was not a parameter a case
could set. `Physics::nu2` and `nu_b2` add it, interpolated on the volume
fraction; see [Two viscosities](#two-viscosities-ba-et-al-eq-22).

One ingredient of Ba et al. is still not implemented. Their MRT collision is
there for stability; the collision here is regularised, which discards the ghost
moments outright rather than relaxing them at a chosen rate.

### Locating the interface: Ba et al. Eq. (21)

The interface is where the two components occupy equal volume. In terms of the
volume fraction of the heavy component,

$$c = \frac{\rho_1}{\rho_{1,0}} = \frac{\rho\,(1 + \varphi)}{2\rho_{1,0}},$$

that is $c = 1/2$ — and in terms of the colour field it is
$\varphi = (\rho_1 - \rho_2)/(\rho_1 + \rho_2)$, **not** $\varphi = 0$. At a
density ratio of 20 the half-volume point sits at φ = 0.905, at 1000 at
φ = 0.998. The colour field's zero contour is by then deep inside the light
fluid. Ba et al.'s normalised field

$$\varphi_N = \frac{\rho_1/\rho_{1,0} - \rho_2/\rho_{2,0}}{\rho_1/\rho_{1,0} + \rho_2/\rho_{2,0}} = 2c - 1$$

is exactly the volume-fraction indicator, so $\varphi_N = 0$ is the interface at
any density ratio. Note also that $\rho = (\rho_{1,0} + \rho_{2,0})/2$ happens
iff $c = 1/2$: **the density interface radius is the physical one**, which is
why the tests here score Laplace's law against it.

Two things followed, and the first was a plain bug.

**The prescribed radius was not the radius the case got.** `droplet_interface()`
lays a tanh down and the case sets `p1_inf` to carry a jump of σ/R for that same
R. But the tanh was laid down in φ, so the *droplet* — the region actually full
of heavy fluid — came out much smaller: measured, a droplet asked for at R = 10
was born at R = 8.41 at a ratio of 20 and at **R = 6.27 at 1000**. The case
therefore started out of equilibrium by the difference between σ/10 and σ/6.27,
and the mismatch grew with the density ratio. `CaseConfig::initial_profile_field`
now prescribes the profile in φ_N and inverts it through
`phase_from_normalised`; the droplet comes out at 10.06 at a ratio of 20 and at
9.95 at 1000.

**Ω⁽²⁾ was reading a stale field.** `phi_n_` was built once in `initialize()`
and never refreshed, so `--interface-field=normalised` had always taken its
gradient of the field as it was at t = 0. That is fixed — `update_interface_field()`
now runs every step, after `phase_field()`. It does not change the conclusion
recorded here previously (the measured jump moved by 0.2 % of itself), but the
earlier negative result was reached with only half the change in place: the
gradient swapped, the initial condition not, and the tension operator not.

### Applying the tension: Ba et al. Eqs. (23)–(29)

Moving the gradient onto φ_N *alone* makes things far worse — the jump collapses
to 6 % of σ/R at a ratio of 1000. The integrated tension is not the problem:
measured on the same relaxed state, $\sum|\nabla\varphi| = 2.00000$ and
$\sum|\nabla\varphi_N| = 1.99939$. What changes is *where* the operator acts,
and Ω⁽²⁾ is sensitive to that in a way that is easy to miss.

A perturbation added post-collision reaches the momentum flux multiplied by the
relaxation time — the injected stress accumulates over successive steps as
$\sum_n (1 - 1/\tau)^n = \tau$ — so Ω⁽²⁾ carries a $1/\tau$ of its own to
cancel it. That cancellation needs τ to be roughly uniform over the region where
the operator acts. It is not: $\tau = \rho\nu/(p\,\mathrm dt) + 1/2$ runs from
5.5 in the light fluid to 5.0 × 10³ in the heavy one on the Laplace case at a
ratio of 1000, three nodes apart. Placed where the mass actually is — deep on
the heavy side — the operator injects an amplitude 10³ times smaller, into
populations that stream within a few nodes into fluid that damps them 400 times
faster than they accumulate. The tension does not survive the trip.

Ba et al. apply the tension as a body force instead:

$$\mathbf F_s = -\tfrac12 \sigma K \nabla\varphi_N, \qquad
  K = -\nabla\!\cdot\!\mathbf n, \qquad
  \mathbf n = -\nabla\varphi_N / |\nabla\varphi_N|$$

with the curvature formed explicitly (their Eq. 26, which is $-\nabla\cdot\mathbf n$
rewritten to stay bounded when the discrete normal is not exactly a unit
vector). The ½ makes $|\nabla\varphi_N|/2$ an interface delta function, since
φ_N runs from −1 to +1. This reaches the momentum equation through Guo's
forcing, whose prefactor $1 - 1/(2\tau)$ stays in [½, 1] however large τ grows,
and it needs no separate velocity redefinition: `macroscopic()` already forms
$\rho\mathbf u = \sum_i \xi_i f_i + \mathbf F\,\mathrm dt/2$, which is their
Eq. (29). `SurfaceTension::ContinuumSurfaceForce` selects it;
`Solver::surface_force()` is the implementation.

An explicit curvature brings a hazard the stress form does not have. The
curvature is built from the **unit** normal, so it does not shrink with the
gradient it came from: in the bulk, where that gradient is nothing but
round-off, the curvature is a large random number and the force
σK|∇φ_N|/2 fails to vanish. It bites in proportion to the density ratio,
because φ_N amplifies noise in φ by dφ_N/dφ = ρ₁/ρ₂ at φ = +1. Measured in the
heavy bulk of the droplet at a ratio of 10⁴, |∇φ_N| starts at 7 × 10⁻¹² and
grows steadily; it crosses the old `kGradientEpsilon` of 10⁻¹⁰ at step 9000 and
the run is destroyed by step 19000. `kInterfaceGradientFloor` (10⁻⁶) is the
threshold that keeps the force on the interface — five orders below the physical
gradient of 0.66 and three above the noise at a ratio of 10⁵. Ba et al. state
the same rule in words after their Eq. (24): the tension "is applied only at the
lattice sites where two fluids coexist".

That the force is applied correctly can be checked directly rather than
inferred. On a relaxed droplet at a ratio of 1000, integrating the solver's own
capillary force along the outward normal gives 0.027755 against a measured
pressure jump of 0.025991 — and the part of the integral outside the interface,
0.026097, matches the jump to 0.4 %. The scheme balances the force it applies;
the residual is the diffuse interface, over which the local curvature 1/r varies
by a third.

### What it costs and what it buys

Measured on the shipped Laplace case at a density ratio of 20, switching the
tension to the CSF form and locating the interface with φ_N:

| | before | after |
|---|---|---|
| Δp / (σ/R_ρ) | 0.895 | **1.021** |
| max &#124;u&#124;, spurious currents | 1.19 × 10⁻³ | **1.66 × 10⁻⁵** |
| R_ρ at t = 0, asked for 10 | 8.41 | **10.06** |

The 72-fold drop in spurious currents is the clearest single sign that the
tension is now applied on the interface rather than beside it, and it puts this
case in the range Ba et al. report (10⁻⁵ for ratios below 100).

The remaining error is not uniform in the parameters, and it is worth being
plain about which parts are physics and which are the measurement. At a ratio of
20 the error falls from +2.1 % at R = 10 to +0.2 % at R = 25, so most of it at
R = 10 is the interface width against the radius. At a ratio of 1000 it is about
−6 % at both radii, and it depends on the interface width — +10 % at
`ch_width_ope` = 1.2, −7 % at 1.6, −29 % at 2.4 — because on a curved diffuse
interface the local curvature varies across the profile and "the radius" stops
being a single number. No width was tuned: 1.6 is the value the case already
had.

What did change is that the error is now *bounded*. The old perturbation
operator, scored the same way, gives −8.5 % at (ratio 20, R = 10), +0.3 % at
(20, 25), −1.4 % at (1000, 10) and **+117 %** at (1000, R = 25, ν = 1.664): it
swings with viscosity and radius. The CSF form gives +2.09 % and +2.10 % at a
ratio of 20 for viscosities a decade apart — the same number.

### The recolouring: why `p` and not `rho`

The last ingredient of Ba et al. not shared with this code was the recolouring.
They use the Latva-Kokko segregation operator, their Eq. (30), with β = 0.7,
and the error at a ratio of 1000 is strongly sensitive to the width the operator
here holds — +10 % at `ch_width_ope` = 1.2, −7 % at 1.6, −29 % at 2.4, while at
a ratio of 20 the same sweep barely moves. That pointed straight at it.

Written in this code's variables the two operators turn out to be the same
object with a different prefactor. From Eq. (30),

$$g_i'' = \varphi f_i' + 2\beta W_i \frac{\rho_R \rho_B}{\rho}\cos\phi_i,
  \qquad \frac{\rho_R\rho_B}{\rho} = \frac{\rho\,(1-\varphi^2)}{4},$$

so Latva-Kokko is $W_i\,\beta\rho(1-\varphi^2)\cos\phi_i/2$ where
`Solver::recolor()` has $W_i\,p\,(1-\varphi^2)(\xi\cdot\hat n)/(2wc_s^2)$.
The difference is **p against ρ**. Matching them in the light bulk of the
shipped case gives β = p/(w c_s² ρ) = 0.625, close to Ba et al.'s 0.7 — so at a
density ratio of 1 they are the same operator, and the question is what happens
as the two bulks separate.

`Recolouring::LatvaKokko` implements it. Measured against the operator already
here, everything else held at the configuration above:

| ratio | `width` | LK β = 0.3 | LK β = 0.62 | LK β = 0.7 | LK β = 0.9 |
|---|---|---|---|---|---|
| 1 | +2.6 % | +3.1 % | +3.6 % | +3.0 % | +2.1 % |
| 2 | +2.4 % | +3.5 % | +0.6 % | **+0.3 %** | **+0.3 %** |
| 10 | +2.1 % | −3.0 % | −21.6 % | −7.8 % | diverges |
| 10³ | −7.6 % | diverges | diverges | diverges | diverges |

At a ratio of 1 the two agree, which is the check that the implementation is
faithful rather than the conclusion. At 2 Latva-Kokko is the better operator.
From 10 upwards it comes apart, and by 10³ no β in [0.01, 0.9] survives.

The reason is the prefactor. The segregation flux is proportional to ρ under
Latva-Kokko and to p under the operator here, and **p is continuous across an
interface where ρ is not** — 0.36 against 0.33 on the Laplace case, while ρ goes
from 10³ to 1. Relative to the local `f`, the strength is a uniform β/2 under
Latva-Kokko and 0.31/ρ under the other: no single β is both strong enough on the
light side to hold the interface and weak enough on the heavy side not to tear
it.

What it tears is the bound on φ, and that is visible long before the divergence.
Every other configuration here holds |φ| ≤ 1 to 10⁻⁹. Latva-Kokko does not: at a
density ratio of **2**, where it gives the best Laplace jump of any operator in
this code, it already runs to |φ| = 1.00058. The equation of state is only
defined on [−1, 1], so that overshoot is what grows into the divergence as the
ratio rises. `test_unit_test_lbm_solver_latva_kokko_overshoots_the_phase_field`
pins it.

The overshoot is structural rather than a matter of tuning. In the two-population
form the bound comes from positivity of `f_R` and `f_B`, which stream separately.
Here `g` is not streamed: it is rebuilt every step as `f * phi + omega_3`, so
there is no positivity left to enforce |φ| ≤ 1, and a segregation flux strong
enough to hold a sharp interface is strong enough to push past it.

This is worth stating plainly because it inverts the expectation. The
pressure-weighted recolouring is not an approximation to Latva-Kokko that this
code happens to use; at high density contrast it is the reason the scheme works
at all, and Latva-Kokko is the one that cannot be carried over. The option is
kept, defaulting off, because it is the better operator below a ratio of about
2 and because the negative result should stay reproducible:

    laplace_opt --recolouring=latva-kokko --beta=0.7

That leaves MRT as the only untried ingredient of Ba et al., and it is the one
aimed at stability — which is what the CSF configuration lacks above 10³.

### Two viscosities: Ba et al. Eq. (22)

The relaxation time is

$$\tau = \frac{\rho\,\nu}{p\,\mathrm dt} + \frac12,$$

so with a single kinematic viscosity for both components the *dynamic*
viscosity ρν — and with it τ — spans the whole density ratio. On the Laplace
case at a ratio of 1000 that is τ = 5.5 in the light fluid against 5.0 × 10³ in
the heavy one, three nodes apart, which is what breaks the stress form of the
tension operator (above). It is also, on its own, a missing capability: with one
`nu`, the viscosity *ratio* of the two fluids is not something a case could set,
and it was pinned at exactly the density ratio.

`Physics::nu2` and `nu_b2` give component 2 its own viscosities, negative
meaning "the same as component 1". Where they differ the kinematic viscosity at
a node is interpolated on the volume fraction of component 1, c = (1 + φ_N)/2:

    nu(x) = c nu_1 + (1 - c) nu_2

which is exactly ν₁ and ν₂ in the two bulks. Ba et al. instead blend the
relaxation *rate* parabolically over a band around the interface; both are
smooth and agree in the bulks. The volume fraction is used here because it is
the quantity φ_N already measures, and because when the two viscosities are
equal the solver skips the interpolation entirely — so every case shipped here
is unaffected to the last bit, which is checked.

Setting ν₂ = ν₁ ρ₁/ρ₂ matches the dynamic viscosities and makes τ uniform across
the interface. **That does not by itself rescue the density ratio.** Measured at
10⁴, with τ = 5.49 on both sides instead of 5.5 against 5 × 10⁴, the run still
diverges. The τ span is what makes Ω⁽²⁾'s `1/τ` cancellation fail; it is not
what limits the CSF configuration.

### Trying for 10⁴, and how the earlier claims went wrong

Leclaire et al. (2011) is the only colour-gradient work reporting O(10⁴): a
steady bubble, 0.5 % on Laplace's law, across "a sharp interface, with a
thickness of about 5-6 lattice units". Their ingredients are the isotropic
colour gradient — which is already here, `src/lbm/isotropic_gradient.h` cites
that paper — and a recolouring operator in the Latva-Kokko family, which
[does not transfer](#the-recolouring-why-p-and-not-rho). Everything else in the
recent literature stops at 10³: Ba et al. (2016), Leclaire et al. (2013), Saito
et al. (2023), and Subhedar (2022), whose survey puts colour-gradient accuracy
as "still limited to a density ratio of 1000".

Before any of that could be tested here, a plainer problem had to be dealt with:
none of the runs behind the earlier claims were long enough to mean anything.
The viscous relaxation time is τ = ρν/(p dt) + 1/2, so with the shipped
viscosity it grows with the density ratio, and the number of relaxation times a
run of 3 × 10⁴ steps actually covers is:

| ratio | τ in the heavy fluid | relaxation times in 3 × 10⁴ steps |
|---|---|---|
| 20 | 1.0 × 10² | 299 |
| 10³ | 5.0 × 10³ | 6.0 |
| 10⁴ | 5.0 × 10⁴ | **0.60** |
| 10⁵ | 5.0 × 10⁵ | **0.06** |

**Every number this file quoted above a density ratio of a few hundred was a
transient.** The runs had not relaxed; neither the divergences nor the pressure
jumps meant what they appeared to. Run to 1.2 × 10⁵ steps instead, the density
ratio of 10³ diverges in both configurations, and so does everything above it.
The earlier claim that the scheme "runs at 10⁵" was a run given six per cent of
one relaxation time before being declared a success.

Matching the dynamic viscosities with `nu2 = nu rho1/rho2` fixes the budget — τ
is then 5.5 everywhere and 3 × 10⁴ steps is 5.5 × 10³ relaxation times — and is
the reason the two-viscosity option above is worth having even for a case that
does not care about the viscosity ratio physically. With that in place, the
question can finally be asked properly. Measured, R = 10 in 128², CSF tension,
φ_N interface, matched viscosities:

| interface width | 10³ | 10⁴ | 10⁵ |
|---|---|---|---|
| 1.6 | −4.4 % | diverges | diverges |
| 2.0 | −12.7 % | +6.8 % | diverges |
| 2.4 | −21.2 % | +4.8 % | −13.8 % |

The +6.8 % and +4.8 % look like success. They are not: **run the same two cases
to 6 × 10⁴ steps and both are `nan`.** At 3 × 10⁴ steps they had not yet
diverged, which is all that column says. The same caution applies to the 10⁵
entry, and at R = 25 the same settings mostly diverge outright, reading −21 %
and −36 % where they do not.

Nor is 10³ steady under this configuration: the first column's −4.4 % becomes
+13.5 % at 6 × 10⁴ steps and −20.7 % at 1.2 × 10⁵. Matching the viscosities buys
a τ small enough to relax the *viscous* modes; it does not make the interface
settle, and something slower is still moving.

Two other candidates were tried and ruled out. Matching the dynamic viscosities
so τ is uniform does not help (above). Raising `kInterfaceGradientFloor` from
10⁻⁶ to 10⁻³, in case the bulk curvature noise was still driving it, changes the
trajectory at 10³ not at all — same values step for step, and it dies slightly
sooner.

What the failure looks like, traced at 10⁴, is a pressure collapse rather than
anything to do with the phase field: |φ| stays inside 1 to 5 × 10⁻¹² throughout
while the minimum pressure falls 0.27 → 0.24 → 0.18 → 0.07 and then goes
negative. Since τ = ρν/(p dt) + ½, a pressure heading for zero takes the
dissipation with it, and that is the runaway. The density changes by about 5.6×
per lattice node at 10³ with the shipped width and 10× at 10⁴, which is the
scale of the problem the pressure has to stay smooth across.

So: **10⁴ was not reached, and neither was 10³.** Getting there inside this
model would mean resolving that density profile, and the model puts the whole
density ratio into a single field through the equation of state. The literature
that reaches 10⁴ does not: Leclaire et al. carry ρ_R and ρ_B as separate smooth
fields, and their density ratio lives in the equilibrium's rest weight rather
than in an equation of state. That is a different model, not a missing term.

### Where the literature is

| Work | Reached | How |
|---|---|---|
| Leclaire et al. 2011 | O(10⁴), static Laplace, 0.5 % error | isotropic colour gradient + Latva-Kokko recolouring |
| Leclaire et al. 2013 | 10³, dynamic | enhanced equilibrium distributions |
| Ba et al. 2016 | 10³, dynamic, high Re; 0.74 % on σ, u_max 1.3 × 10⁻⁴ | MRT collision, CSF perturbation, normalised phase field |
| Saito et al. 2023 | 10 accurately, 10³ marginally | sixth-order Hermite equilibria, central moments |

Ratios of 10⁵ and beyond are reported by other families — chemical-potential
pseudopotential models (> 6.5 × 10⁴), phase-field Allen–Cahn models, entropic
MRT — not by colour-gradient ones.

Ba et al.'s static droplet is the closest comparison: R = 25 in 100², ν = 0.1667,
σ = 0.1, interface 4–5 nodes, 0.74 % error on σ at a ratio of 1000. This code at
the same R and ν gives **+0.20 % at a ratio of 20**, which matches them; at 1000
it reads −4.94 % after 3 × 10⁴ steps, but that is not a converged state and the
run does not survive a long one, so there is no comparable figure to give.

The gap is not in the ingredients: all four are implemented, one of them was
already here, and the fourth — the recolouring — turned out to be the one that
must *not* be carried over. What remains untried is MRT.

## A second model, for density ratios past 10³

The section above ends by saying that reaching 10⁴ would mean a different model
rather than a missing term. `src/lbm/two_population_solver.h` is that model.
Both are colour-gradient methods; they differ in where the density ratio lives,
and that turns out to decide how far it can go.

`Solver` carries the mixture `f` and the colour difference `g`, and gets the
ratio from a two-component equation of state p = p(ρ, φ). The whole contrast
sits in one density field. `TwoPopulationSolver` carries `f¹` and `f²` and
streams them separately; the ratio comes from a free parameter in the
*equilibrium*, each fluid keeping its own rest-particle weight and hence its own
sound speed:

$$(c_s^k)^2 = \tfrac{3}{5}(1-\alpha_k), \qquad p_k = \rho_k (c_s^k)^2,
  \qquad \frac{\rho_1}{\rho_2} = \frac{1-\alpha_2}{1-\alpha_1}.$$

Two things follow that the first model cannot have. The bulk pressures are
**identically equal** — substituting the third relation into the second gives
ρ₁(c_s¹)² = ρ₂(c_s²)² — so the pressure is continuous across the interface
however large the ratio. And each ρ_k has its own smooth profile, so nothing has
to resolve a jump of ρ₁/ρ₂ within a single field.

The price is the classical limitation the first model was written to escape: the
density ratio and the sound-speed ratio are no longer independent. At a ratio of
10⁴ the heavy fluid's sound speed is 0.0069 in lattice units. For a static or
slow flow that is a fair trade; for anything acoustic it is not.

### Adapting the recolouring to a density ratio

Latva-Kokko and Rothman's segregation operator, as Ba et al. write it in their
Eq. (30), pushes colour along the **lattice weight** `w_i`:

$$f_i^{1\ddagger} = \frac{\rho_1}{\rho} f_i^\dagger
  + \beta\,w_i\,\frac{\rho_1\rho_2}{\rho}\cos\phi_i .$$

That cannot work once the two fluids sit on different rest weights, and the
arithmetic is worth doing because it is not obvious from the formula. What has
to stay non-negative is fluid 2's share of the population, `(ρ₂/ρ) f_i†`. At a
density ratio of 1000, at the point where the two occupy equal volume, the
non-rest populations are

    f_i^1 = rho_1 (1 - alpha_1)/5 = 500 x 1.6e-4 = 0.08
    f_i^2 = rho_2 (1 - alpha_2)/5 = 0.5 x 0.16   = 0.08

so `f_i† = 0.16` and fluid 2's share of it is `(0.5/500.5) × 0.16 = 1.6e-4`.
The push is `β w_i ρ₁ρ₂/ρ = 0.7 × (1/9) × 0.4995 = 0.039`, **240 times larger**.
The light fluid's distribution goes negative on the first step and the run is
gone within ten. Measured, before the fix: fine at a ratio of 2, −69 % at 5,
and `nan` from 30 upwards.

Leclaire, Reggio and Trépanier (2012) report adapting this operator "for the
case of variable density ratios", and Leclaire et al. (2011) credit that
adaptation with reaching 10⁴. Their text was not available here, so the form
used is derived from the two requirements the operator has to meet rather than
transcribed. Replacing `w_i` with the mixture's own rest weight

$$\Phi_i = \frac{\rho_1\varphi_i^1 + \rho_2\varphi_i^2}{\rho}
  = \frac{f_i^{\mathrm{eq}}(\rho, 0)}{\rho}$$

gives both:

- **each fluid's mass is conserved exactly**, because `Φ_i` depends only on
  `|e_i|` and the cosine is odd, so the push sums to zero over the directions;
- **both distributions stay non-negative for any β ≤ 1**, because near
  equilibrium the push is `β(ρ₁/ρ)` of fluid 2's own share, and ρ₁ ≤ ρ.

It reduces to Latva-Kokko exactly when the two fluids share a rest weight, which
is the case their form was written for. Measured, `min f` over a fifty-step
droplet is **exactly 0** at density ratios of 20, 10³ and 10⁵, and |φ_N| is
exactly 1 — the phase field is bounded because the populations are, not because
anything clamps it.

### The viscosity has to be matched, not the kinematic viscosity

With the recolouring fixed the model runs at every ratio, and at 10³ it is still
63 % low. The reason is the same τ that limits the first model, arriving by a
different route. Since τ = μ/(p dt) + ½ and the two bulk pressures are equal,
giving both fluids the same *kinematic* viscosity — as Ba et al. do, ν = 0.1667
— puts τ at 348 in the heavy fluid against 0.85 in the light one. A BGK
collision at τ = 348 leaves a viscous stress large enough to absorb most of the
capillary force. Measured on a relaxed droplet at a ratio of 1000: the force
integral is 0.00413, correct to within 3 % of σ/R, while the pressure jump is
0.00176 — **43 % of the force actually applied**. The rest is viscous stress
against the residual spurious current.

Ba et al. carry that with an MRT collision, which relaxes the ghost moments at
their own rate and leaves the shear mode at τ = 348. The alternative, taken
here, is to match the **dynamic** viscosities, which makes τ uniform at 0.85 and
needs no MRT. It is a different physical case — a viscosity ratio of ρ₁/ρ₂
rather than 1 — and the shipped case says so in its own comments.

### What it measures

Ba et al.'s static droplet, R = 25 in 100², σ = 0.1, α₂ = 0.2, β = 0.7,
periodic, matched dynamic viscosity μ = 0.1667, run to 1.2 × 10⁵ steps:

| ratio | σ_cal / σ | steady? | max &#124;u&#124; | Ba et al. |
|---|---|---|---|---|
| 100 | **1.0024** | constant from 6 × 10⁴ steps | 2.9 × 10⁻⁵ | 1.0069, 6.8 × 10⁻⁵ |
| 10³ | **1.0030** | constant from 6 × 10⁴ steps | 4.0 × 10⁻⁵ | 1.0074, 1.25 × 10⁻⁴ |

At 10² and 10³ this is better than the paper on both the tension and the
spurious currents, and unlike anything in the first model it is *converged* —
the last six reports of the run agree to every digit printed.

The `laplace_high_ratio` program is that case, and its tests pin it.

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
| the above, φ_N interface + CSF tension | 0.028415 | 9.95 | **1.021** | **1.66 × 10⁻⁵** |

The jump error fell from 28 % to under 4 % with the equation of state restored,
went back up to 10 % when the interface was started in mechanical equilibrium —
the price of running past a density ratio of 100 at all — and is now 2 %. The
last row is the current default. It is also where the spurious currents fell by
a factor of 72, which is the part that matters most: those currents are what
limit every other use of the scheme.
`--initial-profile=colour --interface-field=colour --surface-tension=perturbation`
reproduces the row above it, and `--initial-state=linear` on top of that the row
above that.

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
| the above, φ_N interface + CSF | 12.30 | 9.95 | 2.35 |

Nothing moves it, which is why it is attributed to the equation of state.

It is, however, no longer an open question *which* of the two is the interface.
The density midpoint is: ρ = (ρ₁ + ρ₂)/2 holds exactly where the two components
occupy equal volume, whereas the colour field's zero contour sits at
φ = (ρ₁ − ρ₂)/(ρ₁ + ρ₂) — 0.905 here, 0.998 at a density ratio of 1000 — which
is inside the light fluid and further inside it the larger the ratio. So the gap
is not two candidate radii of one interface; it is the interface, and a contour
of φ that is not on it. See
[Locating the interface](#locating-the-interface-ba-et-al-eq-21).

The gap opens within the first ~3000 steps and then holds to six digits — a
converged state, not a drift — and the total mass of each colour is conserved
exactly throughout.

`programs/solvers/color_gradient/laplace/tests/` pins the gap so that a scheme
change moving it is visible. The equilibrium is no longer the place to look for
it: the enhanced equilibrium of Leclaire et al. (2013) and the third-order term
of Ba et al. (2016) are the same expression as the one already in
`Solver::equilibrium`, checked to 10⁻¹⁴ — see
[Reading Ba et al. (2016) and Leclaire et al. (2013)](#reading-ba-et-al-2016-and-leclaire-et-al-2013).

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
of itself, so the fix does not by itself change any conclusion; what changed the
conclusion was applying the whole of Ba et al.'s prescription rather than a
third of it. `update_interface_field()` now runs in `step()`, after
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
