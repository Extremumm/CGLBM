# Numerical kernels

Stateless functions the solvers are built from. All are in `cglbm::lbm`.

## Equation of state

`src/lbm/equation_of_state.h` — the two-component equation of state of
Lafarge et al. (2021).

```cpp
struct ComponentPair {
    double c1_squared, c2_squared;   // squared sound speeds of components 1 and 2
    double p1_inf, p2_inf;           // pressures at infinity
};
double pressure(double rho, double phi, const ComponentPair&);
double pressure_linear_mixing(double rho, double phi, double cs_squared, const ComponentPair&);
```

`pressure` returns

```
p = ½ [ ρĉ² − p1_inf − p2_inf + √( (p2_inf − p1_inf + ρc̄²)² + ρ²(1 − φ²)c₁²c₂² ) ]
ĉ² = (c₁² + c₂²)/2 + φ(c₁² − c₂²)/2
c̄² = (c₁² − c₂²)/2 + φ(c₁² + c₂²)/2
```

which is `ρc₁² − p1_inf` at φ = +1 and `ρc₂² − p2_inf` at φ = −1. `1 − φ²` is
clamped at zero, so a φ that overshoots ±1 by a rounding error does not produce
a NaN.

`pressure_linear_mixing` returns `ρc_s² − (1 + φ)p1_inf/2 − (1 − φ)p2_inf/2`.
It is not the model's equation of state; it is kept so that the old behaviour
can be reproduced for comparison.

## Gradient stencils

`src/lbm/isotropic_gradient.h` (2D) and `src/lbm/isotropic_gradient_3d.h` (3D).
The gradient of a scalar field is

```
∇φ(x) = Σ_l W(|c_l|²) c_l φ(x + c_l)
```

with the shell weights normalised so that `Σ_l W c_a c_b = δ_ab`, which makes
every stencil exact for a linear field.

| Stencil | Shells `|c|²` : W | Points | Reach | Isotropic to |
|---|---|---|---|---|
| 2D `E4` | 1 : 1/3, 2 : 1/12 | 8 | 1 | 4th order |
| 2D `E6` | 1 : 4/15, 2 : 1/10, 4 : 1/120 | 12 | 2 | 6th order |
| 2D `E8` | 1 : 4/21, 2 : 4/45, 4 : 1/60, 5 : 2/315, 8 : 1/5040 | 24 | 2 | 8th order |
| 3D `E4` | 1 : 1/6, 2 : 1/12 | 18 | 1 | 4th order |
| 3D `E6` | 1 : 2/15, 2 : 1/15, 3 : 1/60, 4 : 1/120 | 32 | 2 | 6th order |

2D E4 is exactly `w_k/c_s²` over the D2Q9 neighbours. 3D E4 is the D3Q19
neighbourhood.

2D functions:

| Function | Behaviour |
|---|---|
| `stencil_points(stencil, &count)` | the static array of `StencilPoint {cx, cy, weight}`; `count` may be null |
| `stencil_reach(stencil)` | 1 for E4, 2 otherwise |
| `stencil_name(stencil)`, `stencil_from_name(name, &stencil)` | `"E4"`/`"E6"`/`"E8"`; parsing accepts upper or lower case, returns false and leaves `stencil` untouched otherwise |
| `gradient_periodic(field, nx, ny, i, j, stencil, &gx, &gy)` | both axes wrap |
| `gradient_wall_y(field, nx, ny, i, j, stencil, &gx, &gy)` | x wraps; a point with `j + cy` outside `[0, ny − 1]` is dropped, so the gradient becomes one-sided within `reach` nodes of a wall. For E4 this is term for term the original treatment |
| `gradient(field, nx, ny, i, j, stencil, boundary, &gx, &gy)` | `gradient_wall_y` for `Boundary::WallY`, `gradient_periodic` otherwise |

`field` is `nx·ny` values indexed `[i·ny + j]`. The 3D functions are the same
with a `k` index and `nz`: `stencil_points_3d`, `stencil_reach_3d`,
`stencil_name_3d`, `stencil_from_name_3d` (E4, E6), `gradient_periodic_3d`,
`gradient_wall_y_3d` (x and z wrap), field indexed `[(i·ny + j)·nz + k]`.

`Boundary` is declared here: `PeriodicY` (both axes periodic) or `WallY`
(resting walls at `j = 0` and `j = ny − 1`, x periodic).

## Mixture quantities

`src/lbm/mixture.h`. φ is a mass-fraction difference: `Y₁ = (1 + φ)/2`. The
equation of state is the pressure-equilibrium mixture
`1/ρ = Y₁/ρ₁(p) + Y₂/ρ₂(p)` with `ρ_k(p) = (p + p_k,inf)/c_k²`, so component k
fills the volume fraction `α_k = ρY_k/ρ_k(p)`. Mass and volume fractions differ
by a logistic shift, `logit Y₁ = logit α₁ + ln(ρ₁/ρ₂)`: a tanh profile of width
W in one is the same profile in the other, moved by `(W/2)ln(ρ₁/ρ₂)`.

| Function | Returns |
|---|---|
| `component1_density(p, components)` | `(p + p1_inf)/c₁²` |
| `component2_density(p, components)` | `(p + p2_inf)/c₂²` |
| `volume_fraction(ρ, φ, p, components)` | `α₁ = ρ(1 + φ)/2 / ρ₁(p)`, not clamped |
| `normalised_phase(φ, p, components)` | `ψ = (α₁ − α₂)/(α₁ + α₂)` with φ clamped to [−1, 1] and each specific volume `c_k²/max(p + p_k,inf, 10⁻³⁰⁰)`; ρ cancels, so ψ ∈ [−1, 1] |
| `mixture_from_volume_fraction(α₁, p, components)` | `MixtureState {ρ, φ}` with `ρ = α₁ρ₁(p) + α₂ρ₂(p)` and `φ = (α₁ρ₁(p) − α₂ρ₂(p))/ρ`; `pressure(ρ, φ)` returns `p` for it |
| `mixture_kinematic_viscosity(φ, ν₁, ν₂)` | `Y₁ν₁ + Y₂ν₂` with φ clamped, i.e. `ρν = α₁μ₁ + α₂μ₂` |

`normalised_phase(φ, p, components)` uses the pressure-dependent reference
densities `ρ_k(p)`; `normalised_phase(φ, rho1, rho2)` in `case_config.h` uses
the constant bulk densities. They agree where the pressure is the bulk one.
`Solver` uses the latter for φ_N and `mixture_kinematic_viscosity` for
`ViscosityMixing::Dynamic`; the rest of this header is used by the unit tests.

## Capillary stress

`src/lbm/surface_force.h`.

| Name | Definition |
|---|---|
| `kInterfaceGradientThreshold` | 1e-10: below this `|∇ψ|` there is no interface |
| `unit_normal(gx, gy, &nx, &ny)` | `∇ψ/|∇ψ|`, or (0, 0) at or below the threshold |
| `capillary_stress(σ, gx, gy, &Txx, &Txy, &Tyy)` | `T = σ/2 (|g| I − g⊗g/|g|)`, zero at or below the threshold |
| `layer_weight(ψ, ∇·n, W)` | `1/J` with `J = 1 + d ∇·n`, `d = −W atanh(ψ)` (ψ bounded to ±(1 − 10⁻¹²)), J clamped to [¼, 4] |
| `surface_force(Txx, Txy, Tyy, nx, ny, i, j, stencil, boundary, &fx, &fy)` | `F_a = ∂_b T_ab`, the three components differentiated with `gradient(…, stencil, boundary)` |

The stress integrates to σ across a flat interface, because ψ changes by 2 and
`|∇ψ|/2` is a delta function. In the continuum `∇·T = ½σκ∇ψ` with
`κ = −∇·(∇ψ/|∇ψ|)`. On the lattice, a centred divergence sums to zero over a
periodic domain, so the force exerts no net force on the fluid for any
interface shape. `layer_weight` divides each layer's stress by the Jacobian
between it and ψ = 0 (Kublik & Tsai 2016), so that a circle of radius R
carries σ/R whatever the interface width; only the velocity-based solver uses
it.
