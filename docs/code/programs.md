# Programs

Each directory under `programs/solvers` with a `main_<name>.cpp` builds the
executables `<name>_opt` and `<name>_dbg` into the same path under `bin/`. A
program is a case definition: it fills a configuration, hands it to a solver
and returns. Values below are in lattice units; the programs derive them from
physical ones with `c_dx = 1e-5 m` and `c_dt = c_dx / (347√3) s` where they say
so.

All programs write into the current working directory
([Output](output.md)). Run them through `utilities/run_case.sh <name>`, which
makes `artifacts/<name>` and logs to `run.log` there, or through
`pycglbm.testing.run_program`.

## Colour-gradient programs on `CaseConfig`

These accept every option of [the command line](case-configuration.md#command-line--parse_command_line),
print `describe(config)` (plus the lines listed) before running, and return 0,
1 (a solver exception, printed as `<name>: <message>`) or 2 (a bad argument).

### `laplace` — `programs/solvers/color_gradient/laplace`

A static droplet; Laplace's law Δp = σ/R. Solver: `Solver`.

| Setting | Value |
|---|---|
| lattice, steps, interval | 128 × 128, 30 000, 1000 |
| ρ₁, ρ₂ | 20, 1 |
| c₁ = c₂ | 1/√3 (347 m/s) |
| σ, R | 0.27683 (1/(c_dx³/c_dt²)), 10 |
| ν = ν_b | 1.66383 (10⁻² m²/s), both components |
| W_init, W | 1.1, 1.6 |
| p2_inf, p1_inf | 0, `matched_p1_inf`; `matched_pressure_offset = true` |
| boundary, stencil | `PeriodicY`, E8 |
| initial phase | `droplet_interface()` |
| options | `interface_field = BulkNormalised`, `surface_tension = ContinuumSurfaceForce`, `warn_phase_out_of_range = true` |

Tests: `tests/test_laplace_color_gradient.py` (the shipped case and the same
with `E4`: initial jump, steady state, jump ≈ 1.021 σ/R, interface split
against (W/2) ln(ρ₁/ρ₂), spurious currents, E8 ≤ E4) and
`tests/test_laplace_high_density_ratio.py`, which runs

```
laplace --rho1=1e4 --nu=1.6638e-4 --nu-b=1.6638e-4 --nu2=1.6638 --nu-b2=1.6638 \
        --viscosity-mixing=dynamic --surface-tension=stress
```

(the exact viscosities are computed in the test) and checks the flags reached
the solver, the initial state, boundedness, settling, the jump ≈ 1.018 σ/R and
the currents.

### `capillary` — `color_gradient/capillary`

A perturbed flat interface between walls, oscillating under surface tension.
Solver: `Solver`.

| Setting | Value |
|---|---|
| lattice, steps, interval | 128 × 128, 10 000, 100 |
| ρ₁, ρ₂; c₁ = c₂ | 4, 1; 1/√3 |
| σ, R | 0.0055367 (0.02 of the laplace unit), 10 (R only enters `p1_inf`) |
| ν = ν_b | 0.016638 (10⁻⁴ m²/s) |
| boundary, stencil | `WallY`, E4 |
| initial phase | `cosine_layer(0.2, inverted = true)`: component 1 below |
| options | `matched_pressure_offset`, `track_interface` (writes `interface.csv`); scheme options at their defaults (`Colour`, `Perturbation`) |

### `gravity_capillary` — `color_gradient/gravity_capillary`

As `capillary`, with gravity 2.7157e-10 (9.81 m/s²), 50 000 steps, interval
1000, and `cosine_layer(0.2, inverted = false)`: component 1 above.

### `rayleigh_taylor` — `color_gradient/rayleigh_taylor/serial`

Rayleigh–Taylor instability, σ = 0. Solver: `Solver`, serial.

| Setting | Value |
|---|---|
| lattice, steps, interval | 128 × 1028, 5 000 000, 10 000 |
| ρ₁, ρ₂; gravity | 4, 1; 2.7157e-8 (9.81·10² m/s²) |
| ν = ν_b | 0.016638 |
| boundary, stencil | `WallY`, E4 |
| initial phase | `cosine_layer(0.2, inverted = false)`: heavy above |

`matched_pressure_offset` is not set.

### `rayleigh_taylor_omp` — `color_gradient/rayleigh_taylor/openmp`

The same case at 1024 × 4096, 2 000 000 steps, amplitude 0.1, with
`parallel = true`. It allocates several GB; lower it with `--nx`/`--ny`.

### `rayleigh_taylor_cuda` — `color_gradient/rayleigh_taylor/cuda`

A standalone CUDA port of the Rayleigh–Taylor case, built only with
`WITH_CUDA` and a CUDA toolkit. It does not use `Solver` or `CaseConfig`: the
scheme is re-implemented as kernels, with its own constants (256 × 1024,
50 000 steps, interval 5000, ρ₁ = 4, ρ₂ = 1). Command line:
`rayleigh_taylor_cuda [E4|E6|E8] [steps]`. The stencil argument is parsed and
printed, but the kernels use the E4 stencil whatever it says. It exits with 1
when no CUDA device is available.

### `laplace_high_ratio` — `color_gradient/laplace_high_ratio`

Laplace's law at ρ₁/ρ₂ = 1000 with `TwoPopulationSolver`; the static droplet
benchmark of Ba et al. (2016).

| Setting | Value |
|---|---|
| lattice, steps, interval | 100 × 100, 40 000, 5000 |
| ρ₁, ρ₂; σ; R | 1000, 1; 0.1; 25 |
| μ | 0.1667 in both: `nu = μ/ρ₁`, `nu2 = μ/ρ₂` (and the bulk values) |
| α₂, β | 0.2, 0.7 |
| boundary, stencil, initial phase | `PeriodicY`, E8, `droplet_interface()` |

### `laplace_3d` — `color_gradient/laplace_3d`

Laplace's law in 3D, Δp = 2σ/R, at ρ₁/ρ₂ = 1000. Solver:
`TwoPopulationSolver3D` on D3Q19.

| Setting | Value |
|---|---|
| lattice, steps, interval | 48³, 15 000, 2500 |
| ρ₁, ρ₂; σ; R; μ | 1000, 1; 0.1; 10; 0.1667 |
| α₂, β | 0.2, 0.7 |
| boundary, stencil, initial phase | periodic, 3D E6, `droplet_interface_3d()` |
| parallel | true |

### `oscillation_3d` — `color_gradient/oscillation_3d`

A droplet released from a mode-2 spheroid; its frequency against Lamb's.
Solver: `TwoPopulationSolver3D`, D3Q19.

| Setting | Value |
|---|---|
| lattice, steps, interval | 48³, 4000, 1000 |
| ρ₁, ρ₂; σ; R; μ | 10, 1; 0.1; 10; 0.06 |
| initial phase | `oscillating_droplet_3d(0.1)` |
| stencil, options | 3D E6; `track_interface` (the semi-axes every step), `parallel` |

### `rayleigh_taylor_3d` — `color_gradient/rayleigh_taylor_3d`

Rayleigh–Taylor in 3D, one square-cell mode, σ = 0. Solver:
`TwoPopulationSolver3D`, D3Q19.

| Setting | Value |
|---|---|
| lattice, steps, interval | 32 × 128 × 32, 3000, 250 |
| ρ₁, ρ₂; gravity; μ | 3, 1; 10⁻⁴; 0.04 |
| boundary, stencil | `WallY`, 3D E4 |
| initial phase | `cosine_layer_3d(0.047, inverted = false)`: heavy above |

## MHD programs

### `hartmann` — `programs/solvers/mhd/hartmann`

Hartmann flow between insulating walls, against its closed form. Solver:
`TwoPopulationSolver3D` on D3Q27, one fluid (ρ₁ = ρ₂ = 1, the phase field +1
everywhere).

| Setting | Value |
|---|---|
| lattice, steps, interval | 8 × 64 × 8, 20 000, 4000 |
| ν | 0.164089 (both components, shear and bulk) |
| σ_e | 1 (both) |
| Ha | 10: `B_y = Ha / (L √(σ_e/μ))`, `L = 32` the half-width |
| drive | `body_force[0] = u_core σ_e B² / (1 − 1/cosh Ha)` with `u_core = 0.02` |
| MHD | enabled, `b = (0, B_y, 0)`, tolerance 10⁻¹², 500 iterations |
| boundary, stencil | `WallY`, 3D E4; parallel |

Before running it prints `hartmann_number`, `half_width`, `layer_thickness`,
`core_velocity` and `magnetic_damping_time`; after, `final_max_potential` and
`final_charge_imbalance`.

### `magnetic_rayleigh_taylor` — `mhd/magnetic_rayleigh_taylor`

Rayleigh–Taylor under a field normal to the interface; the growth rate against
the quasi-static MHD dispersion relation. Solver: `TwoPopulationSolver3D` on
D3Q27.

| Setting | Value |
|---|---|
| lattice, steps, interval | 128 × 128 × 4, 20 000, 250 |
| ρ₁/ρ₂ | 3.6501457725947524 (heavy above); μ = 0.03 in both |
| gravity, σ | 4·10⁻⁵, 10⁻³ |
| MHD | enabled, `b = (0, 1, 0)`, σ_e1 = 3.81·10⁻², σ_e2 = σ_e1/10 825, **arithmetic** face average, tolerance 10⁻¹⁰ |
| initial phase | `tanh(((j − ny/2) − 0.5 cos(2πi/nx))/W)`, one wavelength along x |
| boundary, stencil | `WallY`, 3D E4; parallel |

`--sigma-e1=X` sets the field strength; unless `--sigma-e2` is also given, the
program resets σ_e2 to σ_e1/10 825 after parsing, so the conductivity contrast
is kept. `--sigma-e1=0` is the hydrodynamic reference. It prints
`wavenumber`, `layer_depth`, `seed_amplitude`, `density_ratio`,
`conductivity_ratio` and `magnetic_damping_time_heavy`. The test runs it at
σ_e1 = 0, 7.62·10⁻³ and 3.81·10⁻², about an hour each on 24 cores.

## Velocity-based programs

`droplet` and `layers` use [`velocity_based::Solver`](velocity-based-solver.md)
directly. They take positional arguments, not the `--key=value` options, print
their parameters as `key = value` lines, and write the four CSV grids
themselves (precision 10) at t = 0 and every 1000 steps over 10 000 steps. The
domain is periodic. The light fluid has ρ₂ = 1 and `ν₂ = 0.166383`
(10⁻³ m²/s, τ ≈ 1); `μ₁ = viscosity_ratio · μ₂`. σ = 0.27683 and W = 1.6.

### `droplet` — `programs/solvers/velocity_based/droplet`

```
droplet [E4|E6|E8] [density_ratio] [velocity] [viscosity_ratio]
```

Defaults E8, 10⁴, 0, 1. A droplet of radius 10 in 128 × 128:
`c = (1 − tanh((r − 10)/W))/2`, `u_x = velocity · c` (the surroundings at rest),
`p = (σ/R)c` (the Laplace jump in place). Velocity 0 is the static Laplace
benchmark.

Tests (`tests/test_droplet_velocity_based.py`): the static droplet at 10⁴ (jump
0.998 σ/R, currents 2.3·10⁻⁶); a droplet launched at 0.01 at ratios 10⁴ and 100
(bounded, momentum conserved to 10⁻⁹ of the initial momentum, slowing down);
launched at 0.1 at viscosity ratios 1, 10, 100.

### `layers` — `programs/solvers/velocity_based/layers`

```
layers [E4|E6|E8] [density_ratio] [amplitude] [viscosity_ratio]
```

Defaults E8, 10⁴, 0.01, 1. A Kolmogorov flow across two layers in 8 × 128:
component 1 fills `0 < y < 64` (a tanh profile of width W about each
interface), the body force per unit volume `G = U μ₂ k² sin(ky)`,
`k = 2π/128`, drives the steady flow `u_x = U (μ₂/μ) sin(ky)`, and both layers
start in it.

Tests (`tests/test_layers_velocity_based.py`): stays layered, conserves
momentum, and the light layer's offset from the analytic profile is 1.9 % of
U, the heavy one's 0.25 %.
