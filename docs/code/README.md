# Code reference

This is the reference for the source code: what each file declares, what each
function computes, how the data is laid out in memory, and in what order a time
step does its work. It describes the code as it is. For the physics and the
measurements behind the design, see [`../numerics.md`](../numerics.md) and
[`../report`](../report/); for the OpenMP, MPI and CUDA modules, see
[`../parallel.md`](../parallel.md).

| Page | Covers |
|---|---|
| [Storage and lattices](storage-and-lattices.md) | `Field`, `Field3D`, D2Q9, D3Q19, D3Q27, `Lattice3D` |
| [Case configuration](case-configuration.md) | `CaseConfig`, `Physics`, the option enums, the initial-phase factories, every command-line flag, the run-log keys |
| [Colour-gradient solver](colour-gradient-solver.md) | `cglbm::lbm::Solver`, the equation-of-state model, stage by stage |
| [Two-population solvers](two-population-solvers.md) | `TwoPopulationSolver` (D2Q9) and `TwoPopulationSolver3D` (D3Q19, D3Q27) |
| [Magnetohydrodynamics](mhd.md) | `QuasiStaticMhd3D`, `MhdPhysics`, the potential solvers |
| [Velocity-based solver](velocity-based-solver.md) | `velocity_based::Solver` and the kernels of `velocity_based.h` |
| [Numerical kernels](numerical-kernels.md) | equation of state, gradient stencils, mixture quantities, capillary stress |
| [Output](output.md) | `CsvWriter`, the VTK writers, the file formats, `run.log` |
| [Programs](programs.md) | every program under `programs/solvers`: the case, the command line, the output, the tests |
| [Build and tests](build-and-tests.md) | CMake targets and options, the pytest suite, the unit-test programs, `pycglbm` |

## Layout

```
src/
  main_cglbm.cpp      placeholder that names the library target; not compiled
  lbm/                the schemes (namespace cglbm::lbm)
  omp/                OpenMP wrappers (namespace cglbm::omp)
  mpi/                MPI environment, decomposition, halo exchange (cglbm::mpi)
  cuda/               CUDA device query, errors, device memory (cglbm::cuda)
programs/
  solvers/<group>/<case>/main_<case>.cpp     one executable per case
  solvers/<group>/<case>/tests/test_*.py     its long tests
  unit_testing/<module>/<name>/main_*.cpp    one executable per unit test
  unit_testing/<module>/<name>/test_*.py     the pytest that runs it
PyCGLBM/pycglbm/      reading, plotting and testing helpers (Python)
cmake/                options, compiler flags, per-file exceptions, CTest
docs/                 this reference, numerics.md, parallel.md, references.md, report/
utilities/            run_case.sh, clean_artifacts.sh
```

Every `.cpp` under `src/` (recursively, except files whose name starts with `_`
and `src/main_cglbm.cpp`) is compiled into one static library,
`libcglbm_opt.a` / `libcglbm_dbg.a`, and every program links it. See
[Build and tests](build-and-tests.md).

## Modules of `src/lbm`

| File | Declares | Used by |
|---|---|---|
| `field.h`, `field3d.h` | `Field`, `Field3D`: lattice-shaped `double` storage | every solver |
| `d2q9.h` | `kQ`, `kXi`, `kW`, `LatticeUnits`, `kGradientEpsilon`, `kInterfaceGradientFloor` | 2D solvers |
| `d3q19.h`, `d3q27.h`, `lattice3d.h/.cpp` | 3D velocity sets and the `Lattice3D` descriptor | `TwoPopulationSolver3D` |
| `equation_of_state.h/.cpp` | `ComponentPair`, `pressure`, `pressure_linear_mixing` | `Solver`, `case_config` |
| `isotropic_gradient.h/.cpp` | `GradientStencil`, `Boundary`, the 2D gradient stencils | every 2D solver |
| `isotropic_gradient_3d.h/.cpp` | `GradientStencil3D`, the 3D gradient stencils | `TwoPopulationSolver3D` |
| `mixture.h/.cpp` | volume fractions, `normalised_phase(phi, p, components)`, `mixture_from_volume_fraction`, `mixture_kinematic_viscosity` | `Solver` (dynamic viscosity mixing), unit tests |
| `surface_force.h/.cpp` | `capillary_stress`, `surface_force`, `unit_normal`, `layer_weight` | `Solver` (capillary stress), `velocity_based::Solver` |
| `case_config.h/.cpp` | `CaseConfig`, `Physics`, option enums, initial-phase factories, `parse_command_line`, `describe` | every colour-gradient program |
| `solver.h/.cpp` | `Solver`: the equation-of-state colour-gradient model | `laplace`, `capillary`, `gravity_capillary`, `rayleigh_taylor`, `rayleigh_taylor_omp` |
| `two_population_solver.h/.cpp` | `TwoPopulationSolver` | `laplace_high_ratio` |
| `two_population_solver_3d.h/.cpp` | `TwoPopulationSolver3D`, `MacroscopicState3D`, `write_vtk_3d` | the 3D and MHD cases |
| `quasi_static_mhd_3d.h/.cpp` | `QuasiStaticMhd3D`, `MhdPhysics`, `PotentialSolver`, `magnetic_damping_time` | `TwoPopulationSolver3D` |
| `velocity_based.h/.cpp` | the kernels of the velocity-based scheme | `velocity_based::Solver` |
| `velocity_based_solver.h/.cpp` | `velocity_based::Solver`, `SolverParameters`, `NodeState` | `droplet`, `layers` |
| `output_writer.h/.cpp` | `MacroscopicState`, `CsvWriter`, `OutputError`, `write_vtk` | the solvers' `run()` |

Dependencies run one way: `field*.h` and the lattice headers depend on nothing;
the kernels (`equation_of_state`, `isotropic_gradient*`, `mixture`,
`surface_force`) depend only on them; `case_config` depends on the kernels, the
3D lattices and `quasi_static_mhd_3d.h` (for `MhdPhysics`); the solvers depend on
all of the above and on `output_writer`. `velocity_based*` depends on
`isotropic_gradient`, `mixture` and `surface_force`, and not on `case_config`.

## The four solvers

| | `Solver` | `TwoPopulationSolver` | `TwoPopulationSolver3D` | `velocity_based::Solver` |
|---|---|---|---|---|
| Header | `solver.h` | `two_population_solver.h` | `two_population_solver_3d.h` | `velocity_based_solver.h` |
| Lattice | D2Q9 | D2Q9 | D3Q19 or D3Q27 | D2Q9 |
| Populations | `f` (mixture), `g` (colour difference) | `f1`, `f2`, one per fluid | `f1`, `f2` | `g` (carries `P = p/(ρc_s²)` and `u`), `h` (carries the volume fraction) |
| Density ratio from | two-component equation of state | rest weights `α_k` | rest weights `α_k` | `ρ = ρ₂ + c(ρ₁ − ρ₂)`; incompressible components |
| Configured by | `CaseConfig` | `CaseConfig` | `CaseConfig` | `SolverParameters` |
| Boundaries | periodic, or walls at `j = 0` and `j = ny − 1` | same | periodic, or walls in y | periodic only |
| Writes its own output | yes, `run()` | yes, `run()` | yes, `run()` (mid-plane slice) | no; the program writes it |
| OpenMP | `config.parallel` | `config.parallel` | `config.parallel` | no |

## Conventions

**Namespaces.** Library code is in `cglbm::lbm`, `cglbm::omp`, `cglbm::mpi` or
`cglbm::cuda`. The velocity-based scheme is in `cglbm::lbm::velocity_based`, so
its `Solver` and `kQ` do not collide with `cglbm::lbm::Solver` and
`cglbm::lbm::kQ`.

**Units.** Everything inside the library is in lattice units, `dx = dt = 1`
(`LatticeUnits` keeps them explicit). A case program converts physical values
with its own `c_dx` and `c_dt` before handing them to the solver; this is the
only place a physical unit appears.

**Memory layout.** A 2D field stores element `(i, j, k)` at
`(i * ny + j) * depth + k`: `i` along x, `j` along y, `k` the component or the
velocity index. A 3D field stores `(i, j, k, q)` at
`((i * ny + j) * nz + k) * depth + q`. A scalar field's `data()` pointer is
what the gradient stencils take, indexed `[i * ny + j]` (2D) or
`[(i * ny + j) * nz + k]` (3D). The velocity-based solver uses plain
`std::vector<double>` with the same `i * ny + j` node index, times `kQ` for
populations.

**Phase field sign.** Component 1 is `phi = +1` (or `c = 1`), component 2 is
`phi = −1` (`c = 0`). Every shipped case puts the denser fluid in component 1.

**Errors.** Constructors throw `std::invalid_argument` for a configuration that
cannot run; the writers throw `OutputError` (a `std::runtime_error`) when a
file cannot be opened or written. The programs built on `CaseConfig` catch
`std::exception` around the solver, print it to stderr prefixed with the
program name, and return 1; they return 2 for a command-line error and 0 for
`--help`. `droplet` and `layers` return 2 for a bad argument and do not catch.

**Naming.** Types are `CamelCase`, functions and variables `snake_case`,
compile-time constants `kCamelCase`, private members end in `_`. See
[`../../CONTRIBUTING.md`](../../CONTRIBUTING.md).
