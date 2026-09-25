# Colour-gradient solver — `src/lbm/solver.h`

`cglbm::lbm::Solver` is the equation-of-state colour-gradient model of
Lafarge et al. (2021), with the options of Ba et al. (2016) and the capillary
stress. It carries two D2Q9 distributions per node:

- `f`, the mixture: `Σf = ρ`, `Σf ξ = ρu − F/2`;
- `g`, the colour difference: `Σg = ρφ`, with φ = (ρ₁ − ρ₂)/ρ.

The density ratio comes from the two-component equation of state
`p = p(ρ, φ)` ([Numerical kernels](numerical-kernels.md#equation-of-state)),
not from the equilibrium, so it is independent of the sound-speed ratio.

Programs: `laplace`, `capillary`, `gravity_capillary`, `rayleigh_taylor`,
`rayleigh_taylor_omp`.

## Interface

```cpp
class Solver {
public:
    explicit Solver(CaseConfig config);
    void initialize();                   // t = 0 state
    void step();                         // one time step
    void run();                          // initialize, step `steps` times, write output
    MacroscopicState state() const;      // pointers to rho, u, phi, p
    const CaseConfig& config() const;
    const Field& density() const;        // rho
    const Field& velocity() const;       // u, depth 2
    const Field& phase() const;          // phi, the colour field (not phi_N)
    const Field& pressure() const;       // p
};
```

The constructor throws `std::invalid_argument` if `nx` or `ny` is not
positive, if `initial_phase` is empty, or if `boundary == WallY` with
`ny < 3`. It resolves the constant parts of the configuration once:

| Member | Value |
|---|---|
| `dx_`, `dt_`, `cs2_`, `cs4_`, `cs6_` | from `config.units` |
| `nu2_`, `nu_b2_` | `physics.nu2`, `nu_b2`, or `nu`, `nu_b` when negative |
| `uniform_viscosity_` | `nu2_ == nu && nu_b2_ == nu_b` |
| `dynamic_mixing_` | `viscosity_mixing == Dynamic` |
| `wall_y_` | `boundary == WallY` |
| `normalise_interface_` | `interface_field == BulkNormalised` |
| `components_` | `{c1², c2², p1_inf, p2_inf}` |

## Fields

| Member | depth | Contents |
|---|---|---|
| `rho_`, `rho_mdt_` | 1 | density, and its value at the previous step |
| `p_`, `p_mdt_` | 1 | pressure, and its value at the previous step |
| `u_` | 2 | velocity |
| `phi_` | 1 | colour field φ |
| `phi_n_` | 1 | φ_N = `normalised_phase(φ, rho1, rho2)`; maintained only when `interface_field == BulkNormalised` |
| `grad_phi_` | 2 | gradient of the interface field (φ or φ_N), evaluated once per step |
| `force_` | 2 | total body force of the step: capillary plus gravity |
| `force_surface_` | 2 | capillary body force (`ContinuumSurfaceForce`, `CapillaryStress`) |
| `normal_x_`, `normal_y_`, `gradient_norm_` | 1 | interface normal and `|∇φ_N|` (`ContinuumSurfaceForce`) |
| `stress_xx_`, `stress_xy_`, `stress_yy_` | 1 | capillary stress; allocated only for `CapillaryStress` |
| `f_`, `g_`, `f_eq_` | 9 | populations and the equilibrium of `f` |
| `omega_1_`, `omega_2_`, `omega_3_` | 9 | viscous relaxation, surface tension (Ω⁽²⁾), recolouring |
| `source_` | 9 | forcing term S |

## Time step

```
step():  force → collide → collide_surface → recolor → stream
         → macroscopic → phase_field → update_interface_field
         → update_colour_gradient → surface_force → equilibrium
```

Every stage loops over all nodes, under
`#pragma omp parallel for collapse(2) if (config.parallel)`. Each writes only
its own node except `stream`, which writes the neighbour it streams to. In what
follows, `H` are the Hermite polynomials of the velocity `ξ_k`:
`H_xx = ξ_x² − c_s²`, `H_yy = ξ_y² − c_s²`, `H_xy = ξ_xξ_y`,
`H_ν = (ξ_x² − ξ_y²)/2`, `H_b = (ξ_x² + ξ_y²)/2 − c_s²`, and
`E_k = w_k((H_xx + H_yy)/(2c_s⁴) − H_xxyy/(4c_s⁶))` with
`H_xxyy = ξ_x²ξ_y² − c_s²(ξ_x² + ξ_y²) + c_s⁴`.

### `equilibrium()`

```
f_eq_k = ρ w_k [1 + ξ·u/c_s² + (u_x²H_xx + 2u_xu_yH_xy + u_y²H_yy)/(2c_s⁴)]
       + (p − ρc_s²) [E_k + w_k (u_x(H_yyx + H_xxx) + u_y(H_yyy + H_xxy))/(2c_s⁶)]
```

with `H_xxx = ξ_x³ − 3c_s²ξ_x`, `H_xxy = ξ_x²ξ_y − c_s²ξ_y`, and likewise. The
second line carries the departure from the ideal gas; its third-order term is
the enhanced equilibrium of Leclaire et al. (2013) and Ba et al. (2016)
Eq. (14), identical to 10⁻¹⁴ (checked by the `lbm/solver` unit test).

### `force()`

Builds the source `S_k = S_F + S_Sp + S_t` at every node, all rows included.

- `S_F = w_k [(F·ξ)/c_s² + (u_xF_xH_xx + u_yF_yH_yy + (u_xF_y + u_yF_x)H_xy)/c_s⁴]`,
  Guo's forcing of `force_`.
- `S_Sp = w_k (d_y(3H_ν − H_b) + d_x(−3H_ν − H_b)) / (2c_s⁴)`, where
  `d_x = Σ_k w_k ξ_x (p − ρc_s²)u_x |_(x+ξ_k) / (dt c_s²)` and
  `d_y` likewise with `ξ_y` and `u_y`: the lattice derivative of
  `(p − ρc_s²)u`, which corrects the third-order moment D2Q9 cannot carry. With
  walls, neighbours across a wall are skipped.
- `S_t = (p − p_mdt − (ρ − ρ_mdt)c_s²) E_k`, the temporal correction.

### `collide()`

Regularised collision with separate shear and bulk rates. At each node:

```
ν, ν_b = viscosity_at(i, j)
τ_ν = ρν/(p dt) + ½,   τ_b = ρν_b/(p dt) + ½
f_neq = f − f_eq + S/2
Σ_ν = Σ f_neq H_ν,   Σ_b = Σ f_neq H_b,   Σ_xy = Σ f_neq H_xy
Ω¹_k = w_k (1 − 1/τ_ν)(H_νΣ_ν + H_xyΣ_xy)/c_s⁴ + w_k (1 − 1/τ_b) H_bΣ_b/c_s⁴
```

The relaxation times use the local pressure, not `ρc_s²`.

`viscosity_at(i, j)` returns:

| Case | ν |
|---|---|
| `uniform_viscosity_` | `physics.nu`, exactly |
| `ViscosityMixing::Dynamic` | `mixture_kinematic_viscosity(φ, ν₁, ν₂) = Y₁ν₁ + Y₂ν₂`, `Y₁ = (1 + φ)/2`, i.e. `ρν = cμ₁ + (1 − c)μ₂` |
| `ViscosityMixing::Kinematic` | `cν₁ + (1 − c)ν₂`, `c = (1 + normalised_phase(φ, rho1, rho2))/2` |

and the same for `ν_b`.

### `collide_surface()`

Only for `SurfaceTension::Perturbation`; the body-force forms leave `omega_2_`
at zero. With `C = grad_phi_` at the node:

```
Ω²_k = σ w_k / (4|C| c_s⁴) · [ (2C_xC_y H_xy + (C_x² − C_y²) H_ν)/τ_ν − (C_x² + C_y²) H_b/τ_b ]
```

and zero where `|C| ≤ kGradientEpsilon`. The `1/τ` factors cancel the relaxation
the stress receives before it reaches the momentum equation; that cancellation
is exact only where τ is uniform.

### `recolor()`

With `C = grad_phi_ / dt` and `|C| > kGradientEpsilon`:

| `recolouring` | Ω³_k |
|---|---|
| `InterfaceWidth` | `w_k p (1 − φ²)/(2W) · (ξ_k·C)/(c_s²|C|)`, `W = ch_width_ope` |
| `LatvaKokko` | `w_k · β ρ(1 − φ²)/2 · (ξ_k·C)/(|ξ_k||C|)` for k ≠ 0, and 0 for k = 0 |

and zero elsewhere. φ here is the colour field; `C` is the gradient of the
interface field.

### `stream()`

```
f(x + ξ_k, k') = f_eq_k + Ω¹_k + Ω²_k + S_k/2
g(x + ξ_k, k') = f(x + ξ_k, k') · φ(x) + Ω³_k
```

with `k' = k` except at a wall: at `j = 0`, directions 4, 7, 8 stay at the node
and come back as `k − 2`; at `j = ny − 1`, directions 2, 5, 6 come back as
`k + 2` (half-way bounce-back). x is always periodic; y is periodic under
`PeriodicY`. The same loop stores `rho_mdt_ = rho_` and `p_mdt_ = p_`.

### `macroscopic()`

```
ρ = Σf
F = force_surface_ + (0, −ρg)          (gravity only when g ≠ 0)
u = (Σf ξ + F dt/2) / ρ
```

`force_` is set to that `F`, so the next `force()` uses the same force as the
velocity, as Guo's scheme requires. The capillary force is the one built at the
end of the previous step.

### `phase_field()`

`φ = Σg / Σf`, then `p = pressure(ρ, φ, components_)`. With
`warn_phase_out_of_range`, a node with `|φ| > 1` prints `Error : phi = <φ>`.

### `update_interface_field()` and `update_colour_gradient()`

When `interface_field == BulkNormalised`, `phi_n_ = normalised_phase(φ, rho1,
rho2)` at every node. Then `grad_phi_` is the gradient of `phi_n_` (or of
`phi_` for `Colour`), with the configured stencil, periodic or one-sided at the
walls (`gradient_wall_y`). The three consumers — `collide_surface`, `recolor`,
`surface_force` — read this stored gradient.

### `surface_force()`

| `surface_tension` | `force_surface_` |
|---|---|
| `Perturbation` | unchanged (zero) |
| `ContinuumSurfaceForce` | `n = −C/|C|` where `|C| > kInterfaceGradientFloor`, else 0; `K = n_xn_y(∂_yn_x + ∂_xn_y) − n_x²∂_yn_y − n_y²∂_xn_x` (Ba et al. Eq. 26, the gradients with the same stencil); `F = ½σK|C| n` on interface nodes, 0 elsewhere |
| `CapillaryStress` | `T = σ/2 (|C| I − C⊗C/|C|)` at every node, zero where `|C| ≤ kInterfaceGradientThreshold` = 1e-10 (`capillary_stress`); then `F = ∇·T` with the same stencil and boundary (`surface_force` of `surface_force.h`) |

`C` is `grad_phi_`. The two body-force forms are the same force in the
continuum. Only the divergence form sums to zero over a periodic lattice for
any interface shape, so only it conserves momentum; the `lbm/solver` unit test
checks that on a tilted, off-centre ellipse.

## `initialize()`

For every node:

1. `φ = initial_phase(config, i, j)`; when `initial_profile_field ==
   BulkNormalised` the value is φ_N and is converted with
   `phase_from_normalised(φ_N, rho1, rho2)`.
2. `u = 0`.
3. ρ and p according to `initial_state`:

   | `initial_state` | ρ | p |
   |---|---|---|
   | `LinearDensity` | `ρ₁(1 + φ)/2 + ρ₂(1 − φ)/2` | linear mixing, `ρ((1 + φ)c₁² + (1 − φ)c₂²)/2 − (1 + φ)p1_inf/2 − (1 − φ)p2_inf/2` |
   | `EquationOfStateP` | same | `pressure(ρ, φ)` |
   | `MechanicalEquilibrium` | `density_at_pressure(φ, target)` with `target = p_light + (p_heavy − p_light)(1 + φ)/2`, `p_heavy = ρ₁c₁² − p1_inf`, `p_light = ρ₂c₂² − p2_inf` | `pressure(ρ, φ)` |

4. `rho_mdt_ = ρ`, `p_mdt_ = p`, `force_ = (0, −ρg)`.

Then `update_interface_field`, `update_colour_gradient`, `surface_force` and
`equilibrium`, and finally `f = f_eq`, `g = f φ`. The capillary force computed
here enters from the first step's `macroscopic()`; the first `force()` sees
gravity only.

## `run()`

```
CsvWriter writer(output_precision)
initialize(); writer.write_grids(0, state())
if track_interface: writer.open_interface_track()
for t = 1 .. steps:
    step()
    if track_interface: writer.write_interface(t, phi_)
    if t % interval == 0: print "Step t"; writer.write_grids(t, state())
```

`state()` points at `rho_`, `u_`, `phi_` and `p_`, so `phase_<t>.csv` holds the
colour field φ, not φ_N. `run()` does not print the configuration; the program
prints `describe(config)` first. See [Output](output.md).

## What each option costs

| Option | Extra work per step |
|---|---|
| `BulkNormalised` | one pass for `phi_n_` |
| `ContinuumSurfaceForce` | one pass for the normals, one with two gradients |
| `CapillaryStress` | one pass for the stress, one with three gradients; three extra fields |
| `Dynamic` / `Kinematic` mixing | none when the viscosities are equal (fast path) |
| `track_interface` | a scan of the column `i = nx/2` every step |
