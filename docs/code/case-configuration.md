# Case configuration — `src/lbm/case_config.h`

A colour-gradient case is a `CaseConfig` value: the lattice, the fluids, the
boundary, the initial state and the options of the scheme. A program builds one
in a `*_case()` function, lets `parse_command_line` override it, prints
`describe(config)` and hands it to `Solver`, `TwoPopulationSolver` or
`TwoPopulationSolver3D`. The velocity-based solver does not use it (see
[Velocity-based solver](velocity-based-solver.md)).

## `CaseConfig`

| Member | Type | Default | Read by | Meaning |
|---|---|---|---|---|
| `name` | `std::string` | `""` | run log | case name; output file names do not depend on it |
| `nx`, `ny` | `int` | 128, 128 | all | nodes along x and y |
| `nz` | `int` | 1 | 3D | nodes along z |
| `steps` | `int` | 1000 | `run()` | time steps to advance |
| `interval` | `int` | 100 | `run()` | write the CSV grids every `interval` steps |
| `units` | `LatticeUnits` | dx = dt = 1 | all | lattice spacing and time step |
| `physics` | `Physics` | see below | all | fluid properties |
| `boundary` | `Boundary` | `PeriodicY` | all | `PeriodicY`, or `WallY`: resting walls at `j = 0` and `j = ny − 1` |
| `initial_phase` | `PhaseFieldInit` | empty | 2D | initial phase field at `(i, j)`; required |
| `initial_phase_3d` | `PhaseFieldInit3D` | empty | 3D | initial phase field at `(i, j, k)`; required by the 3D solver |
| `initial_state` | `InitialState` | `MechanicalEquilibrium` | `Solver` | how ρ and p are laid down at t = 0 |
| `initial_profile_field` | `InterfaceField` | `BulkNormalised` | `Solver` | which field `initial_phase` prescribes |
| `interface_field` | `InterfaceField` | `Colour` | `Solver` | which field the colour gradient is taken of |
| `surface_tension` | `SurfaceTension` | `Perturbation` | `Solver` | how the tension reaches the populations |
| `viscosity_mixing` | `ViscosityMixing` | `Kinematic` | `Solver` | how two viscosities mix across the interface |
| `recolouring` | `Recolouring` | `InterfaceWidth` | `Solver` | segregation operator |
| `matched_pressure_offset` | `bool` | false | `parse_command_line` | recompute `physics.p1_inf = matched_p1_inf(physics)` after the overrides |
| `stencil` | `GradientStencil` | `E8` | 2D | isotropy of the 2D gradient stencil |
| `lattice_3d` | `Lattice3DKind` | `D3Q19` | 3D | velocity set |
| `mhd` | `MhdPhysics` | disabled | 3D | imposed magnetic field and conductivities, see [MHD](mhd.md) |
| `stencil_3d` | `GradientStencil3D` | `E6` | 3D | isotropy of the 3D gradient stencil |
| `track_interface` | `bool` | false | `Solver`, 3D | write `interface.csv` every step |
| `parallel` | `bool` | false | all three | run the lattice loops across OpenMP threads |
| `warn_phase_out_of_range` | `bool` | false | `Solver` | print `Error : phi = <value>` when `|φ| > 1` |
| `write_vtk_field` | `bool` | false | 3D | also write `field_<t>.vtk` every `interval` steps |
| `output_precision` | `int` | 6 | all | significant digits in the CSV output |

`node_count()` returns `nx * ny` as a `long`.

`PhaseFieldInit` is `std::function<double(const CaseConfig&, int i, int j)>`
and `PhaseFieldInit3D` is `std::function<double(const CaseConfig&, int i,
int j, int k)>`. Both return a value in [−1, 1], +1 in component 1.

"`Solver`" in the table means the option is read only by the equation-of-state
solver; the two-population solvers ignore it.

## `Physics`

| Member | Default | Meaning |
|---|---|---|
| `rho1`, `rho2` | 1, 1 | bulk densities of component 1 (φ = +1) and 2 (φ = −1) |
| `c1`, `c2` | 1, 1 | sound speeds of the two components (`Solver` only) |
| `nu`, `nu_b` | 0, 0 | kinematic shear and bulk viscosity of component 1 |
| `nu2`, `nu_b2` | −1, −1 | the same for component 2; negative means "equal to component 1" |
| `sigma` | 0 | surface tension |
| `radius` | 0 | prescribed interface radius; `droplet_interface*` use it, and `matched_p1_inf` divides σ by it |
| `gravity` | 0 | acceleration along −y; the force density is `−ρ g` |
| `body_force[3]` | 0, 0, 0 | uniform force density, added as is (`TwoPopulationSolver3D` only) |
| `ch_width_init` | 1.1 | interface width of the droplet profiles at t = 0 |
| `ch_width_ope` | 1.6 | width the `InterfaceWidth` recolouring maintains; also the width of the `cosine_layer*` profiles |
| `beta` | 0.7 | segregation strength of `LatvaKokko` and of the two-population recolouring, in (0, 1] |
| `p1_inf`, `p2_inf` | 0, 0 | pressures at infinity of the two branches of the equation of state (`Solver` only) |
| `alpha2` | 0.2 | rest weight of component 2 (two-population solvers only) |

## Option enums

| Enum | Values | Command line | What it changes |
|---|---|---|---|
| `InterfaceField` | `Colour`, `BulkNormalised` | `colour`/`color`, `normalised`/`normalized` | `Colour` is φ = (ρ₁ − ρ₂)/ρ, a mass-fraction difference. `BulkNormalised` is φ_N = `normalised_phase(φ, rho1, rho2)`, the volume-fraction difference 2c − 1, whose zero is the interface at any density ratio |
| `InitialState` | `LinearDensity`, `EquationOfStateP`, `MechanicalEquilibrium` | `linear`, `eos`, `equilibrium` | see [`Solver::initialize`](colour-gradient-solver.md#initialize) |
| `SurfaceTension` | `Perturbation`, `ContinuumSurfaceForce`, `CapillaryStress` | `perturbation`, `csf`, `stress` | Ω⁽²⁾ in the collision; a body force `½σK|∇φ_N| n` from an explicit curvature; or the body force `∇·T` of the capillary stress |
| `ViscosityMixing` | `Kinematic`, `Dynamic` | `kinematic`, `dynamic` | `ν = cν₁ + (1 − c)ν₂`, or `ρν = cμ₁ + (1 − c)μ₂` |
| `Recolouring` | `InterfaceWidth`, `LatvaKokko` | `width`, `latva-kokko` | segregation rate `p(1 − φ²)/(2W)`, or `βρ(1 − φ²)/2` |
| `Boundary` (in `isotropic_gradient.h`) | `PeriodicY`, `WallY` | none | printed as `periodic_y` / `wall_y` |
| `GradientStencil` (in `isotropic_gradient.h`) | `E4`, `E6`, `E8` | `E4`/`e4`, `E6`/`e6`, `E8`/`e8` | 2D stencil |
| `GradientStencil3D` (in `isotropic_gradient_3d.h`) | `E4`, `E6` | `E4`/`e4`, `E6`/`e6` | 3D stencil |
| `Lattice3DKind` (in `lattice3d.h`) | `D3Q19`, `D3Q27` | `D3Q19`/`d3q19`, `D3Q27`/`d3q27` | 3D velocity set |
| `PotentialSolver` (in `quasi_static_mhd_3d.h`) | `FiniteVolume`, `LatticeBoltzmann` | `fv`/`finite-volume`, `lbm`/`lattice-boltzmann` | MHD potential solve |

The names are matched exactly as listed: `E8` and `e8` are accepted, `E8 ` or
`D3q19` are not.

## Free functions

| Function | Definition |
|---|---|
| `normalised_phase(phi, rho1, rho2)` | clamps φ to [−1, 1], then `(ρ₂(1 + φ) − ρ₁(1 − φ)) / (ρ₂(1 + φ) + ρ₁(1 − φ))`. Maps ±1 to ±1, is the identity when ρ₁ = ρ₂, and crosses zero at φ = (ρ₁ − ρ₂)/(ρ₁ + ρ₂). Not to be confused with the overload `normalised_phase(phi, p, components)` of `mixture.h`, which uses the pressure-dependent reference densities |
| `phase_from_normalised(phi_n, rho1, rho2)` | the inverse: clamps φ_N to [−1, 1], then `(Sφ_N − D)/(S − Dφ_N)` with `S = ρ₁ + ρ₂`, `D = ρ₂ − ρ₁` |
| `matched_p1_inf(physics)` | `ρ₁c₁² − ρ₂c₂² − σ/R`: the offset that puts component 1 at density ρ₁ under the pressure `p₂ + σ/R` while component 2 sits at ρ₂ under `p₂ = ρ₂c₂²` (with `p2_inf = 0`) |
| `density_at_pressure(phi, target, rho1, rho2, components)` | the ρ for which `pressure(ρ, φ, components) = target`, by 100 bisection steps on [¼ min(ρ₁, ρ₂), 4 max(ρ₁, ρ₂)]; returns `ρ₁(1 + φ)/2 + ρ₂(1 − φ)/2` when that bracket holds no sign change |

## Initial-phase factories

Each returns a function suitable for `initial_phase` or `initial_phase_3d`.
Distances are in nodes from the centre `(nx/2, ny/2[, nz/2])`, with integer
division.

| Factory | Profile |
|---|---|
| `droplet_interface()` | `−tanh((r − R)/W_init)`, `r` the distance to the centre, `R = physics.radius · units.dx`, `W_init = ch_width_init`; +1 inside |
| `cosine_layer(A, inverted)` | `tanh(((j − ny/2) − A·nx·cos(2πi/nx)) / W_ope)`, `W_ope = ch_width_ope`; negated when `inverted`. Not inverted, component 1 is above the interface (the unstable arrangement when it is the heavier one) |
| `droplet_interface_3d()` | the sphere, `−tanh((r − R)/W_init)` |
| `oscillating_droplet_3d(ε)` | `−tanh((r − R'(1 + ε P₂(cos θ)))/W_init)`, θ from the z axis, `P₂(x) = (3x² − 1)/2`, `R' = R/(1 + 3ε²/5)^{1/3}` so the spheroid holds the volume of a sphere of radius R; `cos θ = 1` at the centre |
| `cosine_layer_3d(A, inverted)` | `tanh(((j − ny/2) − A·nx·cos(2πi/nx)·cos(2πk/nz)) / W_ope)`, negated when `inverted` |

How the value is read depends on the solver: `Solver` reads it as φ_N when
`initial_profile_field == BulkNormalised` (the default) and as φ otherwise; the
two-population solvers always read it as φ_N = 2c − 1.

## Command line — `parse_command_line`

```cpp
enum class CommandLineResult { Run, Finished, Error };
CommandLineResult parse_command_line(CaseConfig& config, int argc, char** argv,
                                     const std::string& program_name);
```

`Finished` is returned for `--help`/`-h` after printing the usage; `Error` for
an unknown argument or an invalid value, with a message on stderr. After a
successful parse, a positive `--threads` is passed to
`omp::set_thread_count`, and `p1_inf` is recomputed when
`matched_pressure_offset` is set.

Every option is `--key=value`, except the switches `--mhd`, `--no-mhd`,
`--vtk` and `--help`. A bare `E4`, `E6` or `E8` anywhere on the line is the
historical form of `--stencil`. Later arguments override earlier ones.

| Flag | Sets | Accepts |
|---|---|---|
| `--stencil=S` or bare `S` | `stencil`; and `stencil_3d` when S is E4 or E6 | E4, E6, E8 |
| `--interface-field=F` | `interface_field` | colour, color, normalised, normalized |
| `--initial-profile=F` | `initial_profile_field` | same |
| `--initial-state=S` | `initial_state` | equilibrium, eos, linear |
| `--surface-tension=S` | `surface_tension` | perturbation, csf, stress |
| `--viscosity-mixing=M` | `viscosity_mixing` | kinematic, dynamic |
| `--recolouring=R`, `--recoloring=R` | `recolouring` | width, latva-kokko |
| `--beta=X` | `physics.beta` | 0 < X ≤ 1 (read with `atof`) |
| `--rho1=X`, `--rho2=X` | `physics.rho1`, `rho2` | X > 0 |
| `--sigma=X` | `physics.sigma` | X ≥ 0 |
| `--radius=X` | `physics.radius` | X > 0 |
| `--width=X` | `physics.ch_width_ope` | X ≥ 0 |
| `--width-init=X` | `physics.ch_width_init` | X ≥ 0 |
| `--nu=X`, `--nu-b=X` | `physics.nu`, `nu_b` | X ≥ 0 |
| `--nu2=X`, `--nu-b2=X` | `physics.nu2`, `nu_b2` | X ≥ 0 |
| `--nx=N`, `--ny=N`, `--nz=N` | `nx`, `ny`, `nz` | N > 0 |
| `--steps=N`, `--interval=N` | `steps`, `interval` | N > 0 |
| `--precision=N` | `output_precision` | N > 0 |
| `--threads=N` | sets `parallel = true`; thread count applied after parsing | N > 0 |
| `--lattice=L` | `lattice_3d` | D3Q19, d3q19, D3Q27, d3q27 |
| `--bx=X`, `--by=X`, `--bz=X` | `mhd.b[0..2]`, and `mhd.enabled = true` | finite |
| `--mhd`, `--no-mhd` | `mhd.enabled` | switch |
| `--sigma-e1=X`, `--sigma-e2=X` | `mhd.conductivity1`, `conductivity2` | X ≥ 0 |
| `--conductivity=A` | `mhd.harmonic_conductivity` | harmonic, arithmetic |
| `--potential-solver=P` | `mhd.solver` | fv, finite-volume, lbm, lattice-boltzmann |
| `--mhd-tolerance=X` | `mhd.tolerance` | X > 0 |
| `--mhd-iterations=N` | `mhd.max_iterations` | N > 0 |
| `--drive-x=X`, `--drive-y=X`, `--drive-z=X` | `physics.body_force[0..2]` | finite |
| `--vtk` | `write_vtk_field = true` | switch |
| `--help`, `-h` | prints the usage, returns `Finished` | |

Real values are parsed with `std::stod` and must be consumed entirely; integers
are parsed with `std::stoi`, which stops at the first non-digit (`--nx=64abc`
reads 64). There is no flag for `gravity`, `alpha2`, `c1`, `c2`, `p1_inf` or
`p2_inf`: those belong to the case.

## Run log — `describe`

`describe(config)` returns one `key = value` line per setting, doubles printed
with 17 significant digits. Every program prints it at the top of its output,
which `pycglbm` reads back from `run.log` (`CaseOutput.parameter`), so the
tests never restate a constant.

Keys, in order: `case`, `nx`, `ny`, `nz`, `steps`, `interval`, `boundary`
(`periodic_y`/`wall_y`), `stencil`, `stencil_3d`, `lattice_3d`,
`initial_state` (`equilibrium`/`eos`/`linear`), `initial_profile`,
`interface_field` (`colour`/`normalised`), `surface_tension`
(`perturbation`/`csf`/`stress`), `viscosity_mixing` (`kinematic`/`dynamic`),
`recolouring` (`width`/`latva-kokko`), `dx`, `dt`, `rho1`, `rho2`, `c1`, `c2`,
`nu`, `nu_b`, `nu2`, `nu_b2` (resolved: equal to `nu`, `nu_b` when negative),
`sigma`, `radius`, `gravity`, `ch_width_init`, `ch_width_ope`, `beta`,
`alpha2`, `p1_inf`, `p2_inf`, `body_force_x`, `body_force_y`, `body_force_z`,
`mhd` (`true`/`false`), `b_x`, `b_y`, `b_z`, `sigma_e1`, `sigma_e2`,
`conductivity_average` (`harmonic`/`arithmetic`), `potential_solver`
(`fv`/`lbm`), `mhd_tolerance`, `mhd_iterations`, `output_precision`,
`track_interface`, `vtk`, `parallel`.
