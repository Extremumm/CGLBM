# Magnetohydrodynamics — `src/lbm/quasi_static_mhd_3d.h`

`QuasiStaticMhd3D` computes the Lorentz force of an imposed uniform field `B₀`
on a conducting two-fluid flow in the inductionless (low magnetic Reynolds
number) limit:

```
J = σ(−∇φ + u × B₀),    ∇·J = 0,    F = J × B₀
⇒  ∇·(σ∇φ) = ∇·(σ(u × B₀))
```

Here φ is the electric potential, not the phase field. There is no magnetic
field to store and no induction equation to march; the potential is solved
once per time step. `TwoPopulationSolver3D` owns one when
`config.mhd.enabled` is set (see
[Two-population solvers](two-population-solvers.md#the-magnetic-coupling-3d)).

Programs: `hartmann`, `magnetic_rayleigh_taylor`.

## `MhdPhysics`

| Member | Default | Flag | Meaning |
|---|---|---|---|
| `enabled` | false | `--mhd`, `--no-mhd`, any of `--bx/--by/--bz` | run the coupling |
| `b[3]` | 0, 0, 0 | `--bx`, `--by`, `--bz` | `B₀`, lattice units |
| `conductivity1`, `conductivity2` | 1, 1 | `--sigma-e1`, `--sigma-e2` | electrical conductivity of component 1 (φ_N = +1) and 2 |
| `harmonic_conductivity` | true | `--conductivity=harmonic\|arithmetic` | how a face between two compositions averages the two |
| `solver` | `FiniteVolume` | `--potential-solver=fv\|lbm` | discretisation of the potential equation |
| `tolerance` | 1e-10 | `--mhd-tolerance` | relative residual (FV), or relative change between two checks (LBM), to stop at |
| `max_iterations` | 500 | `--mhd-iterations` | iteration or sweep limit; reaching it is reported, not thrown |
| `lbm_tau_max` | 1.0 | none | largest relaxation time of the LBM march (floored at ½ + 10⁻⁶) |
| `lbm_check_interval` | 16 | none | sweeps between convergence checks of the LBM march |

`magnetic_damping_time(ρ, σ, b)` returns `ρ/(σ|b|²)`, or infinity when
`σ|b|²` is not positive. A time step much longer than this is unstable, since
the force is explicit; the MHD programs print it.

`potential_solver_from_name` accepts `fv`, `finite-volume`, `lbm`,
`lattice-boltzmann`; `potential_solver_name` returns `fv` or `lbm`.

## Interface

```cpp
QuasiStaticMhd3D(int nx, int ny, int nz, bool wall_y, MhdPhysics physics, bool parallel);
void solve(const Field3D& velocity, const Field3D& phase);   // phase: phi_N in [-1, 1]
void add_lorentz_force(Field3D& force) const;                // force += J x B0, depth 3
const Field3D& potential() const;
const Field3D& current() const;       // node-centred J, depth 3
int iterations() const;               // of the last solve
double residual() const;              // relative residual reached
bool converged() const;
double charge_imbalance() const;      // max |sum_faces J_f| over the lattice
void apply_operator(const Field3D& in, Field3D& out) const;  // for the tests
void solve_potential(const Field3D& rhs);                    // FV solve of A phi = rhs, for the tests
double face_conductivity(int i, int j, int k, int axis) const;
```

The constructor throws `std::invalid_argument` if a dimension is not positive,
a conductivity is negative, or `max_iterations ≤ 0`.

## `solve(velocity, phase)`

1. **`update_conductivity`.** On the `+axis` face of every node, with
   `c = clamp((1 + φ_N)/2, 0, 1)` at each end and `f` the mean of the two:
   harmonic `σ₁σ₂/(fσ₂ + (1 − f)σ₁)`, or arithmetic `fσ₁ + (1 − f)σ₂`. With
   `wall_y`, the face that wraps around y is set to zero: that is the whole of
   the insulating condition `J·n = 0`.
2. **`update_drive`.** The EMF `u × B₀` at the nodes; on each face
   `drive = σ_f · ½(emf_here + emf_ahead)·n_f`; the right-hand side is
   `rhs = −Σ_axis (drive(x) − drive(x − e_axis))`.
3. **The potential.**
   - `FiniteVolume` (`solve_flat`): conjugate gradient on
     `(Aφ)(x) = Σ_faces σ_f (φ(x) − φ(neighbour))`, Jacobi-preconditioned by the
     row sum, warm-started from the previous potential. Every boundary is
     periodic or insulating, so A has the constant in its nullspace; the mean
     is projected out of the right-hand side, the iterate, the residual and the
     preconditioned residual, and the zero-mean potential is returned. It stops
     when `|r| ≤ tolerance · |rhs|` or after `max_iterations`. A zero
     right-hand side returns φ = 0 with no iteration.
   - `LatticeBoltzmann` (`solve_lattice_boltzmann`): a pseudo-time march of
     `∂φ/∂t = ∇·(σ∇φ) − ∇·(σ(u × B₀))` on D3Q7 (rest weight ¼, `c_s² = 1/4`)
     with a two-relaxation-time collision at magic parameter
     `Λ = (τ⁺ − ½)(τ⁻ − ½) = ¼`, the diffusivity scaled so that the largest
     relaxation time is `lbm_tau_max`. Every `lbm_check_interval` sweeps it
     compares φ with its value at the previous check and stops when
     `max|Δφ| / max|φ| ≤ tolerance`, which is what `residual()` then reports;
     the mean is removed at the end.
4. **`update_current`.** On each face,
   `J_f = drive_f − σ_f(φ(x + e_axis) − φ(x))`, the same difference the operator
   used, so `Σ_faces J_f` at a node is the linear residual; then the node
   current is the mean of the two faces along each axis.

`add_lorentz_force` adds `J × B₀` at every node.

Reductions (`dot`) sum fixed chunks of 4096 nodes in a fixed order, so the
result does not depend on the thread count.

## Choosing the potential solver

| | `FiniteVolume` | `LatticeBoltzmann` |
|---|---|---|
| Accuracy | conservative to the linear tolerance | second order, about 2.5 times the FV error on the same lattice |
| Cost per solve | a few CG iterations when warm-started; three global reductions per iteration | local sweeps, no reductions; from cold `O(L²)` sweeps |
| Large conductivity ratio | unaffected (38 iterations at 10⁴) | does not converge (not after 40 000 sweeps at 10⁴) |

These figures are measured by `programs/unit_testing/lbm/mhd_potential`. The
default is `FiniteVolume`.

## Face conductivity

Harmonic is right for current crossing the interface (series), arithmetic for
current running along it (parallel). At a ratio of 10⁴ and `f = ½`, harmonic is
about 2700 times smaller than arithmetic, which makes the interface band an
insulator for current along it. `magnetic_rayleigh_taylor`, where the current
runs along the interface, sets `harmonic_conductivity = false`.
