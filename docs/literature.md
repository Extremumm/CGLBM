# High density ratios: what the literature does, and what limits it

This is a review of how lattice Boltzmann two-phase models reach density ratios
of 10³ and beyond, written to decide what to try next in this code. It asks
three questions of every paper: how far the density ratio went, **for a static
or a moving interface**, and which ingredient carried it there. The last
section says which of those ingredients this code has, and what the
measurements here add to the picture.

A note on sources. Ba et al. (2016) was read in full, from the author's
accepted manuscript, and so was Saito et al. (2023), from its arXiv
version. Every other entry is from its abstract, from the publisher's
summary, or from how Ba et al. and others describe it, and it is quoted only
for what those say. Full references are in
[`references.md`](references.md).

## Three families, sorted by where the density ratio lives

What decides how far a model goes is less its collision operator or its
surface-tension term than **what its populations carry**. Three answers exist.

**1. Colour-gradient models whose populations carry the density.** The fluid
is `f = f_R + f_B` (or `f` and a colour difference `g`), `sum f = rho`,
`sum f e = rho u`, and the density contrast is either in the equilibrium's rest
weight — Grunau, Chen & Eggert (1993), Reis & Phillips (2007), Leclaire et al.
(2011, 2012, 2013), Ba et al. (2016), Wen et al. (2019), Saito et al. (2023);
`TwoPopulationSolver` here — or in a two-component equation of state — Lafarge
et al. (2021); `Solver` here. The heavy fluid is then a lattice gas whose
kinetic temperature `p / rho` is the density ratio times smaller than the
light fluid's.

**2. Chromodynamic models built from a continuum description.** Lishchuk,
Halliday & Care (2008) map a single inhomogeneous, essentially incompressible
fluid onto a multicomponent lattice Boltzmann method and report "correct static
and dynamic operation up to a fluid density contrast ratio of more than 500".
Burgin, Spendlove, Xu & Halliday (2019) and Spendlove et al. (2020) develop the
same line with an MRT collision and test the kinematic condition (the fluids
do not interpenetrate) and the continuity of traction directly.

**3. Pressure- or velocity-based hydrodynamics with a separate interface.**
The populations carry the pressure and the velocity, not the density; the
density follows from an order parameter transported by its own equation,
usually the conservative Allen–Cahn equation. Inamuro et al. (2004, with a
pressure projection), Lee & Lin (2005), Zu & He (2013), Fakhari et al. (2017),
Hajabdollahi, Premnath & Welch (2021, central moments), Reis (2022, the one-fluid
model), Otomo et al. (2025). Two recent papers come to this family *from* the
colour-gradient one: Subhedar (2022) keeps the colour-gradient segregation but
uses Zu & He's velocity-based equilibrium, so that the interface mobility no
longer depends on the density ratio; Haghani et al. (2024) show that "the CG
method is not fluid invariant" and turn it into a phase-field equation coupled
to a separate hydrodynamic solver. And the group behind `Solver`'s own scheme
has made the same move: Gregorczyk, Zhao & Boivin (2025) solve the mean mass
and momentum with a lattice Boltzmann scheme and the order parameter with a
finite-volume one, following Shao & Shu (2015) and Reis (2022), and simulate
liquid jets through dripping, sinuous and atomisation regimes.

## What each reached

| Work | Family | Static | Moving interface | Carried by |
|---|---|---|---|---|
| Reis & Phillips 2007 | 1, rest weight | — | 18.5, droplet coalescence (per Ba et al.) | alpha_k equilibrium on D2Q9 |
| Lishchuk, Halliday & Care 2008 | 2 | > 500 | > 500 | continuum mapping |
| Leclaire, Reggio & Trépanier 2011 | 1, rest weight | O(10⁴), 0.5 % on Laplace | — | isotropic colour gradient |
| Leclaire, Reggio & Trépanier 2012 | 1, rest weight | — | — | Latva-Kokko recolouring adapted to unequal rest weights |
| Huang, Huang, Lu & Sukop 2013 | 1, rest weight | — | layered channel flow to 8 | a source term cancelling an unwanted term in the momentum equation |
| Leclaire, Pellerin, Reggio & Trépanier 2013 | 1, rest weight | — | layered Couette flow to 1000 | enhanced equilibrium |
| Leclaire, Reggio & Trépanier 2013 (JCP) | 1, rest weight | "high density ratio flows for both steady and unsteady cases" (per Ba et al.) | | MRT + enhanced equilibrium |
| **Ba et al. 2016** | 1, rest weight | **1000**, 0.74 % on σ, u_max 1.25 × 10⁻⁴ | layered channel 1000; Rayleigh–Taylor 3; splashing **100**, Re ≤ 500 | MRT, third-order equilibrium, Chapman–Enskog source term, normalised phase field, CSF tension |
| Wen, Li, Yu & Luo 2019 | 1, 3D | — | droplet impact at "a relatively large density ratio" | error terms of the 3D model removed |
| Burgin et al. 2019; Spendlove et al. 2020 | 2 | — | kinematics verified "for a range of density contrast" | MRT; analysis of the segregation |
| Lafarge et al. 2021 | 1, equation of state | "validation tests up to density ratios of 1000" (abstract) | | two-component EOS, temporal correction |
| Saito et al. 2023 | 1, rest weight | — | accurate at 10, "still limited" at 1000 | sixth-order Hermite equilibria, central moments |
| Subhedar 2022 | 3 (CG segregation) | — | stationary and translating drop, layered Poiseuille, Rayleigh–Taylor | velocity-based equilibrium |
| Inamuro et al. 2004; Lee & Lin 2005 | 3 | — | rising bubbles; splashing, both at 1000 (per Ba et al.) | pressure projection; stable discretisation |
| Zu & He 2013 | 3 | — | "moderate density ratios": moving droplet, Rayleigh–Taylor, layered Poiseuille | velocity-based equilibrium |
| Fakhari et al. 2017 | 3 | — | "large density and viscosity contrasts": Rayleigh–Taylor, Taylor bubble | velocity-based equilibrium, conservative Allen–Cahn, local stress |
| Liang et al. 2018 | 3 | yes | **splashing at 1000**, Re 20–500 | conservative Allen–Cahn, forcing distribution |
| Gregorczyk, Zhao & Boivin 2025 | 3 | — | liquid-jet atomisation | pressure-based LBM + finite-volume order parameter |

Two things stand out. **No colour-gradient model with density-carrying
populations reports a moving interface far past 100**: Ba et al.'s splashing is
at 100, Saito et al. call 1000 "still limited", and the one result past 10³
(Leclaire et al. 2011) is static. Lishchuk et al.'s "more than 500", static and
dynamic, is the furthest any chromodynamic model goes, and it is built from a
continuum description of an incompressible fluid rather than on the rest-weight
or equation-of-state equilibria. The same splashing test is run at 1000 by Liang
et al. (2018) in family 3. And **every model reported at 10³ with a
moving interface stops carrying the density in its populations**.

## The ingredients, and which are here

| Ingredient | Source | `Solver` | `TwoPopulationSolver` |
|---|---|---|---|
| Isotropic colour gradient | Leclaire et al. 2011 | E4, E6, E8 | E4, E6, E8 |
| Interface at the half-volume contour, φ_N | Ba et al. Eq. (21); Leclaire et al. 2013 | yes | yes |
| Tension as a body force | Brackbill et al.; Lishchuk et al. 2003; Ba et al. Eqs. (23)–(29) | CSF or capillary stress | CSF |
| Enhanced (third-order) equilibrium | Leclaire et al. 2013; Ba et al. Eq. (14) | yes, the Hermite term | yes |
| Source term for the diagonal third moment | Huang et al. 2013; Ba et al. Eqs. (17)–(18) | yes, `S_Sp` | **added here**, `--third-moment-correction` |
| MRT collision | Ba et al. Eqs. (11)–(13); Leclaire et al. 2013; Spendlove et al. 2020 | regularised (ghost moments discarded) | **added here**, `--collision=mrt` |
| Recolouring for unequal rest weights | Leclaire et al. 2012 | not applicable | yes, `rest_weight` |
| Viscosity mixed as ρν | — | `--viscosity-mixing=dynamic` | on the volume fraction |
| Generalized equilibria, central-moment collision | Saito et al. 2023, Eqs. (25), (56)–(63) | **added here**, `--collision=central` (D2Q9 form derived here) | no |
| Velocity-based equilibrium | Zu & He 2013; Fakhari et al. 2017; Subhedar 2022 | no | no — the `droplet` solver |
| Interface mobility set on its own | Subhedar 2022 | no | no — `phase_temperature` in the `droplet` solver |
| Walls and gravity, for layered Poiseuille flow and Rayleigh–Taylor | Zu & He 2013; Fakhari et al. 2017; Liang et al. 2018 | yes | no — **added here** to the `droplet` solver: `poiseuille_vb`, `rayleigh_taylor_vb` |
| A viscosity that carries shear across the interface | Liang et al. 2018 (a step at φ = 1/2) | no | no — **added here** to the `droplet` solver, `--viscosity=laminate`: harmonic for the shear across, arithmetic for the stretching along |

Before this review, MRT was the one ingredient of Ba et al. neither
colour-gradient solver had, and the two-population solver also lacked the
third-moment source term. Both are now in it; see
[`numerics.md`](numerics.md#mrt-and-the-third-moment-source-in-the-two-population-solver).

## What limits the density ratio

**Static interfaces.** The literature's account — where the interface is, how
the tension is applied, how uniform τ is — is borne out here and is no longer
the limit: `Solver` holds Laplace's law at 10⁴, 10⁵ and 10⁶ (see
[`numerics.md`](numerics.md#static-droplets-up-to-a-density-ratio-of-a-million)), and the
two-population solver now runs Ba et al.'s own benchmark, with their equal
kinematic viscosities, rather than only a matched-μ variant of it.

**Moving interfaces.** Two structural limits follow from the density being in
the populations, and together they explain the table above.

1. *Positivity of the mass flux.* A heavy fluid moving at `u` has populations
   of order `p / c_s^2 ± rho_1 u` in the directions along and against the
   motion, and `p` is fixed by the light fluid. Behind the interface the
   heavy side streams about `-rho_1 u / 2` into light nodes that hold a mass
   of order `rho_2`. Keeping that non-negative needs `u rho_1 / rho_2` of
   order one — measured here in the earlier monolithic `laplace` program:
   bounded at `u = 10^-3` and 1000, diverging at `10^-2`; bounded at `10^-4`
   and 10⁴, diverging at `10^-3` (`numerics.md`,
   [What is not solved](numerics.md#what-is-not-solved)). The rest-weight
   model has the same bound, since the heavy fluid's moving populations there
   are `(1 - alpha_1)/5 = 0.16 rho_2 / rho_1` of its density.
2. *The heavy fluid is a cold lattice gas.* D2Q9 has `e_x^3 = e_x`, so the
   diagonal third moment of any population set is the momentum itself:
   `sum f e_x^3 = rho u_x`, where the Navier–Stokes equations need `3 p u_x`.
   For a fluid at `p = rho c_s^2` the two agree. For the heavy fluid
   `p / (rho c_s^2)` is the density ratio times smaller, and the defect is
   the whole of `rho_1 u`. It reaches the **normal** viscous stress only — the
   off-diagonal third moment is what the enhanced equilibrium repairs — and
   every model in family 1 cancels it with a finite difference: Huang et al.'s
   source term, Ba et al.'s `C`, `S_Sp` here. Ba et al.'s verdict on Huang et
   al.'s term is that "the evaluation of the source term may lead to
   considerable numerical errors at higher density ratios". The measurements
   here say how large: the finite difference is not the stencil the streaming
   used, the two differ at order `k^2`, and the difference multiplies a
   quantity `rho_1 / rho_2` times larger than the stress it corrects. On a
   Taylor–Green vortex in the heavy fluid — pure normal strain — `Solver`'s
   decay rate is **24.8 times** the Navier–Stokes one at a density ratio of
   1000 on a 32² lattice, and the two-population solver's is 17.5 times with
   Ba et al.'s source term and about 760 times without it; a shear wave in the
   same fluid is right to 0.5 %. The heavy fluid resists extension far more
   than it should, at every wavelength a droplet cares about. That part has a
   cure, found here: the source is split half before and half after the
   collision, and a split source cancels what the streaming did only if its
   derivative has the Fourier symbol `2i tan(k/2)` — the Cayley transform of
   the streaming's shift — rather than `ik`. A reach-three axis stencil that
   matches it through `k⁵` brings the vortex to 0.998 at 1000 and 1.010 at
   10⁴ in both solvers (`--source-stencil=matched`). What it does not cure is
   the interface: at 10³ a moving one is wrong whichever stencil the source
   takes, which is the first limit again. Burgin et al.
   (2019) reach the same place from the other end: the reported restrictions
   on the density contrast "in rapid flow ... originate in the effect on the
   model's kinematics of the terms ... which correct its dynamics for large
   density differences".

**What Saito et al. (2023) find with the same scheme.** Theirs is the most
recent colour-gradient model of family 1 tested at 10³ with a moving
interface, and it carries the same correction as `Solver`,
`Q = −3∇·[(p − ρc_s²)u]` (their Eq. (65)), taken by a finite difference.
Its equilibria are exact to sixth order in the velocity and its collision is
in central moments. At 10³ a rising bubble stays stable for 728 000 steps.
Its centre of mass follows the reference, but its rise velocity does not, and
they conclude that "further improvements are required to solve it more
accurately within the framework of the CG model". The error they point at is
the second limit: "the gradient computation of the correction term Q by
finite differences introduces numerical errors that distort the droplet",
and "a higher-order lattice (e.g., D3Q39 lattice) and the corresponding
third-order equilibrium may be needed to solve it completely". The matched
stencil here makes the same repair inside the heavy fluid without a new
lattice. The new lattice itself does not work, tried here in two dimensions:
on the 17-velocity lattice of Shan, Yuan & Chen (2006) the same equilibrium
carries the third moment exactly, and is linearly unstable for any fluid
cooler than the lattice, by 5 % a step at the heavy fluid's temperature
([numerics.md](numerics.md#a-higher-order-lattice-for-the-heavy-fluid)). Neither touches the first limit, which no lattice can lift. On any
lattice whose moving velocities have integer components, non-negative
populations satisfy `Σ f e_x² ≥ Σ f |e_x| ≥ |Σ f e_x|`, so a fluid carried
by them has `T + u_x² ≥ |u_x|`, with `T = p/ρ` its kinetic temperature: it
cannot move, in lattice units, much faster than it is hot. The heavy fluid's
`T` is the light fluid's divided by the density ratio. The one kinetic
formulation without that bound is Particles on Demand (Dorschner, Bösch &
Karlin 2018; Kallikounis, Dorschner & Karlin 2022). It moves the lattice to
each node's own velocity and temperature, so every fluid is carried at the
lattice's own temperature. It has been shown on compressible single-phase
flows with near-vacuum regions, not on two immiscible fluids. It also
replaces exact streaming with an interpolation or a finite-volume step, so
it would be a different code rather than a change to this one.

Neither limit is a missing term. Both come from asking D2Q9 populations to
carry a fluid whose kinetic temperature is far below the lattice's. Family 3
does not: its populations carry the velocity, or the momentum, at the lattice
temperature `c_s^2` whatever the density, with the pressure as a separate
isotropic part; the density comes from a bounded order parameter, and no
population has to carry a mass flux across the interface. That is why the literature past
10³ with moving interfaces is all family 3, and it is where Lafarge's own group
has gone. In this code it is the velocity-based
`droplet` solver, which moves a droplet at `u = 0.1` at a density ratio of 10⁴.

## What follows for this code

- The **static** density ratio of the colour-gradient solvers is limited by
  nothing the literature points at any more. `Solver` is measured to 10⁶.
- The **dynamic** density ratio of the colour-gradient solvers is limited by
  the two structural effects above, not by any ingredient the literature
  offers. The two-population solver's MRT and source term make it run Ba et
  al.'s benchmark correctly; they do not lift either limit, and the
  measurements say so: its mode-2 droplet at 1000 (`oscillation_high_ratio`)
  rings down 6.6 times too fast with their source as published, and 1.45 times
  with the matched stencil below. The heavy fluid's extensional viscosity is fixed in
  its interior by taking the source on the stencil the streaming uses
  (`--source-stencil=matched`, off by default); against exact normal modes it
  brings a capillary wave at a ratio of 100 from 1.50 to 1.19 times its
  damping, but at 10³ the moving interface is wrong with either stencil —
  3.7 and 2.8 times for the wave, and a mode-2 droplet damped eight times too
  fast or at 0.62 of the rate — where the velocity-based solver is within 5 %.
  The 0.62 follows the heavy fluid's bulk relaxation time (0.95 at τ_b = 1,
  where the wave is no longer damped); taking the trace of the correction in
  its continuum form, (p − ρc²)∇·u, brings the droplet to 1.26 and the wave
  to 2.4, and diverges on static droplets at 10⁴, because only the published
  time difference telescopes over the steps a slow bulk rate holds it
  ([numerics.md](numerics.md#the-heavy-fluids-bulk-rate-and-the-trace-of-the-correction)).
  For the wave that is resolution rather than a wall, below the positivity
  bound: at twice the resolution it is at 2.9 with the nine-point source and
  1.6 with the matched one, and a wider interface no longer helps. The
  droplet was not damped at all with the matched source until its face values
  were limited where the density jumps: the six-point stencil overshoots the
  jump of the quantity it differentiates, and the overshoot, by every
  measurement, cancelled the heavy fluid's viscous normal stress at the
  interface. The limit also carries the wave further past the bound: at 1000
  it diverged only from an amplitude of 8 nodes (from 4 before), and not for
  its negative populations but because φ drifted past 1 on the heavy side,
  where the equation of state's extrapolated pressure fell below zero. With
  φ clamped there the waves at 8 and 10 run to the end, damped 1.9 and 2.0
  times the linear rate; the velocity-based solver rings the one at 8 down at
  1.02.
- Flows with moving interfaces at 10³ and above belong to the velocity-based
  solver, which is the colour-gradient segregation on a family-3
  hydrodynamics — the model Subhedar (2022) describes. As Subhedar also
  argues, its interface mobility is set on its own rather than inherited from
  the lattice: the phase populations are built on a carrier at a lattice
  temperature of 0.2 instead of `c_s²`, which took the capillary wave at
  1000 from 1.084 to 1.049 times the exact damping rate. Filtering the
  collision's non-equilibrium in time, rather than rebuilding part of it from
  the finite-difference velocity gradient, took it to 0.991, and the same
  wave at $10^4$ from 2.10 to 1.27
  ([`numerics.md`](numerics.md#oscillations-against-exact-normal-modes)).
  What the lower temperature took out turns out to be a spurious surface
  diffusion of the phase field, proportional to it: the truncation errors of
  the interface normal and of the memoryless transport move a curved
  interface at `C ∂_s²κ`. Taken out to fourth order (`fourth_order_phase`,
  opt-in), the droplet at $10^4$ goes from 1.76 to 1.07 times the exact rate.
  The fourth-order normal has to be kept to the interface itself: taken in
  its far tails as well, it let a capillary wave on a wavelength of 128 grow
  cells in the heavy fluid.
- The velocity-based solver now runs the wall-bounded benchmarks of family 3:
  resting walls and gravity, with layered Poiseuille flow and the
  Rayleigh–Taylor instability scored against exact solutions. The
  Rayleigh–Taylor mode at a density ratio of 1000 grows at 0.98 of the exact
  viscous rate, 0.99 on twice the resolution. The layered flow showed what
  none of the periodic cases had: the shear stress that crosses a diffuse
  interface sees the two fluids in series, and mixing the viscosity
  arithmetically on the volume fraction, as most of family 3 does, makes the
  light side of a 1000-to-1 interface forty times too viscous a few nodes
  into it. The light layer was 41 % off its profile. Liang et al. (2018) take
  a viscosity that jumps at φ = 1/2 for this test, which they warn may not
  hold up under large changes of topology; the laminate mixing here takes the
  harmonic mean for the shear across the interface and the arithmetic one for
  the stretching along it, and brings the layered flow to 4.6 % (on Liang et
  al.'s own case 0.017 in their L1 error, against their 0.032). It lowers the
  damping of the capillary waves, whose boundary layer the interface does not
  resolve, by 2 to 17 %, and is not the default
  ([numerics.md](numerics.md#walls-gravity-and-the-viscosity-across-an-interface)).
