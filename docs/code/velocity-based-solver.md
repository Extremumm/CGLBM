# Velocity-based solver

`src/lbm/velocity_based.h` holds the kernels and `src/lbm/velocity_based_solver.h`
the time step of a two-phase scheme for interfaces that *move* at large density
ratios. Everything is in `cglbm::lbm::velocity_based`. The scheme never streams
the density:

- the hydrodynamic populations `g` carry `P = p/(ρc_s²)` (zeroth moment) and
  `u` (first moment), both continuous across the interface;
- the phase populations `h` carry the volume fraction `c` of component 1;
- `ρ = ρ₂ + c̄(ρ₁ − ρ₂)`, with `c̄ = clamp(c, 0, 1)`.

Both components are incompressible; the lattice sound speed is an artificial
compressibility. The domain is periodic in x and y. The solver is serial.

Programs: `droplet`, `layers`. Background and measurements:
[`../numerics.md`](../numerics.md#the-velocity-based-droplet-solver).

## Kernels — `velocity_based.h`

Constants: `kQ = 9`, `kVelocity` and `kWeight` in the D2Q9 order of `d2q9.h`,
`kOpposite = {0, 3, 4, 1, 2, 7, 8, 5, 6}`, `kSoundSpeedSquared = 1/3`.

| Function | Computes |
|---|---|
| `velocity_equilibrium(ux, uy, gamma)` | `Γ_k(u) = w_k(1 + ξ·u/c_s² + (ξ·u)²/(2c_s⁴) − u²/(2c_s²))`, the unit-mass Maxwellian; moments 1, u, `c_s²I + uu` |
| `hydrodynamic_equilibrium(P, ux, uy, eq)` | `g_k^eq = Γ_k(u) + w_k(P − 1)`; moments P, u, `Pc_s²I + uu` |
| `phase_populations(c, ux, uy, nx, ny, W, h)` | `h_k = cΓ_k(u) + θ·flux_k` with `flux_k = w_k A (ξ_k·n)/c_s²`, `A = M·2c̄(1 − c̄)/W`, `M = c_s²/2`. θ ∈ [0, 1] is the largest value keeping `0 ≤ c̄Γ_k + θ·flux_k ≤ Γ_k` for every k |
| `forcing(ux, uy, ax, ay, S)` | Guo's source for an acceleration `a`: `S_k = w_k[ξ·a/c_s² + (H_xx·2u_xa_x + H_yy·2u_ya_y + 2ξ_xξ_y(u_xa_y + u_ya_x))/(2c_s⁴)]` |
| `collide(g, eq, S, τ, τ_b, out)` | regularised collision: the second moments `Π` of `g − g^eq + S/2` are relaxed, the trace at `1 − 1/τ_b`, the deviator and `Π_xy` at `1 − 1/τ`, and rebuilt as `out_k = g_k^eq + w_k/(2c_s⁴)[(ξ_x² − c_s²)Q_xx + (ξ_y² − c_s²)Q_yy + 2ξ_xξ_yQ_xy] + S_k/2` |
| `collide_hybrid(g, eq, S, τ, τ_b, σ, ∇u, out)` | as `collide`, but relaxing `σΠ + (1 − σ)Π^CE` with `Π^CE = −τc_s²(∇u + ∇uᵀ − (∇·u)I) − τ_bc_s²(∇·u)I`; σ = 1 is `collide` |
| `pressure_force(P, ρ, nx, ny, i, j, ax, ay)` | `a(x) = Σ_{k≥1} w_k ξ_k ρ(x − ξ_k)P(x − ξ_k) / ρ(x)`, i.e. `Σ w_k ξ_k p(x − ξ_k)/(c_s²ρ(x))`: the lattice gradient of `p` itself, so a uniform `p` exerts no force whatever ρ does, and `ρa` sums to zero over the lattice |
| `link_momentum(k, donor, receiver, ρ₁, ρ₂, jx, jy, D)` | momentum the receiver gains from the donor at `receiver − ξ_k`, and the link's dissipation coefficient; see below |
| `dissipation_force(D, ux, uy, nx, ny, i, j, fx, fy)` | `f(x) = Σ_{k≥1} D_k(x)(u(x − ξ_k) − u(x))`, with `D` stored `[(i·ny + j)·kQ + k]` |
| `set_velocity(g, ux, uy)` | adds `w_k ξ_k·(u − Σgξ)/c_s²` to each moving population: the first moment becomes u, the zeroth and second are unchanged |

### `link_momentum`

`LinkEnd` describes one end after streaming: `outgoing` (the post-collision
population it sent along the link, less its pressure part `w_kP`), `phase`
(the phase population it sent), `rho` and `mu` at the new time, `ux`, `uy` at
the old time. With `ρ_l = min(ρ_donor, ρ_receiver)`:

```
lattice  = ρ_l · (outgoing_d + outgoing_r − A_k(u_d) − A_k(u_r))      A_k: the uu part of Γ_k
K        = Γ_k(u_d) − Γ_k̄(u_r)                                          volume carried
K_lin    = w_k ξ_k·(u_d + u_r)/c_s²                                     its linear part
M        = (ρ₁ − ρ₂)(phase_d − phase_r) + ρ₂K                           mass carried
excess   = M − ρ_lK
J        = lattice · ξ_k + (ρ_lK_lin + excess) · (u_d + u_r)/2
β        = max(0, 2μ_dμ_r/(μ_d + μ_r) − ρ_l(μ_d/ρ_d + μ_r/ρ_r)/2)
D        = |excess|/2 + β · 2w_k/c_s²
```

Swapping the ends and taking the opposite direction gives `−J` and the same `D`,
so the exchange conserves momentum link by link. The quadratic part of `K` is
left out of the advection: it is a diffusion along the flow, and carried within
the step it destabilises a single fluid at `|u| = 0.1`. `D` is not added to the
momentum here; `dissipation_force` applies it through the forcing term.

## `Solver` — `velocity_based_solver.h`

```cpp
struct SolverParameters {
    int nx = 128, ny = 128;
    double rho1 = 1.0, rho2 = 1.0;             // densities of c = 1 and c = 0
    double mu1 = 1.0/6.0, mu2 = 1.0/6.0;       // dynamic viscosities
    double surface_tension = 0.0;
    double width = 1.6;                        // W of c = (1 + tanh(x/W))/2
    double tau_bulk = 1.0;                     // relaxation of the trace
    double hybrid_weight = 0.7;                // sigma of collide_hybrid
    GradientStencil stencil = GradientStencil::E8;
};
struct NodeState { double c = 0, ux = 0, uy = 0, p = 0; };  // p relative to the far field

class Solver {
public:
    explicit Solver(const SolverParameters&);
    void initialize(const std::function<NodeState(int i, int j)>& state);
    void set_body_force(const std::function<void(int i, int j, double* fx, double* fy)>&);
    void step();
    int nx() const; int ny() const; const SolverParameters& parameters() const;
    double volume_fraction(int i, int j) const;   // c, unbounded
    double density(int i, int j) const;           // rho
    double phase(int i, int j) const;             // psi = 2 c_bounded - 1
    double pressure(int i, int j) const;          // P rho cs^2
    double velocity_x(int i, int j) const;        // macroscopic u, with the half acceleration
    double velocity_y(int i, int j) const;
};
```

Storage is `std::vector<double>`, node `m = i·ny + j`, populations at `m·kQ + k`.
The constructor only allocates; nothing is validated.

### `initialize(state)`

For each node: `ρ = ρ₂ + c(ρ₁ − ρ₂)` (from the given `c`, not clamped),
`g = g^eq(p/(ρc_s²), u)`, `h = cΓ(u)`, `D = 0`, `u` and the lattice velocity set
to the given `u`. Then `macroscopic()` and `acceleration()`.

### `set_body_force(force)`

Stores a force *per unit volume* at every node, constant in time. It enters
the next `acceleration()`, divided by ρ.

### `step()`

```
save_step → collide_and_stream → macroscopic → momentum → acceleration → velocity
```

| Stage | What it does |
|---|---|
| `save_step()` | copies ρ, P, u and a to their `_old` arrays |
| `collide_and_stream()` | at each node: `g^eq(P, u)`, `S = forcing(u, a)`, `τ = μ(c)/(ρc_s²) + ½` with `μ(c) = c̄μ₁ + (1 − c̄)μ₂`, `∇u` with the **E4** stencil, `collide_hybrid(g, g^eq, S, τ, tau_bulk, hybrid_weight, ∇u)`; `h = phase_populations(c, u, n, width)`; both streamed periodically to `x + ξ_k` |
| `macroscopic()` | `c = Σh`, `ρ = ρ₂ + c̄(ρ₁ − ρ₂)`, `ψ = 2c̄ − 1`, `P = Σg` |
| `momentum()` | per node, `J = ρ_old(u_old + a_old/2) + Σ_{k≥1} link_momentum(k, donor = x − ξ_k, receiver = x)`, storing each link's `D`; lattice velocity `J/ρ` |
| `acceleration()` | `a = (∇·T + dissipation_force + body)/ρ + pressure_force`, see below |
| `velocity()` | `set_velocity(g, u_lattice)`; macroscopic `u = u_lattice + a/2` |

`acceleration()` in detail:

1. `∇ψ` with the configured stencil, and the unit normal `n = ∇ψ/|∇ψ|`
   (`unit_normal`, zero below 1e-10);
2. `∇·n` with the same stencil; the stress
   `T = capillary_stress(σ · layer_weight(ψ, ∇·n, W), ∇ψ)`: each layer of the
   diffuse interface is weighted by the Jacobian between it and ψ = 0, so a
   circle carries σ/R whatever its width;
3. `F_s = surface_force(T)` (periodic), `F_d = dissipation_force(D, u_lattice)`,
   `a_p = pressure_force(P, ρ)`, and `a = (F_s + F_d + body)/ρ + a_p`.

Every force is applied through the forcing term: half in the post-collision
velocity `u + a/2` of the next step's `J`, half in the macroscopic velocity.
Applied within the step instead, the dissipation drives a period-2 mode at
τ close to ½.

Output is the program's job: `droplet` and `layers` write the four CSV grids
themselves (see [Programs](programs.md#velocity-based-programs)), with `phase`
= ψ and `pressure` relative to the far field.
