# Build and tests

## Targets

`CMakeLists.txt` globs `programs/**/main_*.cpp` (and `main_*.cu` when CUDA is
enabled and found) plus `src/main_cglbm.cpp`. Each file names a target after
its stem without `main_`:

- `src/main_cglbm.cpp` names the library. Its sources are every file matching
  `[!_]*.cpp` (and `.cu` with CUDA) under `src/`, recursively, except
  `main_cglbm.cpp` itself. It builds as `cglbm_opt` / `cglbm_dbg`
  (`libcglbm_opt.a`, `libcglbm_dbg.a`), with `src/` as its public include
  directory, so library headers are included as `"lbm/solver.h"`.
- Every other `main_<name>` builds the executables `<name>_opt` and
  `<name>_dbg` from all `[!_]*.cpp` files in its directory, linked to the
  library, and placed in `BIN_DIR/<path under programs>`.

Aggregate targets: `opt` and `dbg` build every optimised or debug target;
`<name>_all` builds both variants of one program. A new program needs only a
new directory and a reconfigure (`cmake --preset gnu`).

## Options — `cmake/Options.cmake`

| Option | Default | Meaning |
|---|---|---|
| `RELEASE` | ON | build the `_opt` variants |
| `DEBUG` | ON | build the `_dbg` variants |
| `ARCH` | `$ENV{ARCH}` | passed as `-march=` (GCC, Clang) or `-x` (Intel) |
| `WITH_OpenMP` | ON | OpenMP (`_OPENMP`) |
| `WITH_MPI` | ON | MPI (`CGLBM_WITH_MPI`) |
| `WITH_CUDA` | ON | CUDA, when a toolkit is found |
| `WITH_IPO` | OFF | link-time optimisation |
| `WITH_Python` | ON | register the pytest suite with CTest |
| `EXCEPT` | ON | apply the per-file flags of `cmake/Exceptions.cmake` |
| `BIN_DIR` | `<source>/bin` | where the executables go |
| `ARTIFACTS_DIR` | `<source>/artifacts` | where the tests write their runs |

Each defaults to the environment variable of the same name when it is set.
Release flags are `-O3 -fno-fast-math` (GCC adds `-fno-finite-math-only`); see
`cmake/CompilerFlags.cmake`.

Presets (`CMakePresets.json` → `cmake/presets/`):

| Preset | Compiler | Build directory | `ARCH` |
|---|---|---|---|
| `gnu` | `g++` | `build-gnu` | `native` |
| `clang` | `clang++` | `build-clang` | `native` |

Both inherit `base` (OpenMP, CUDA, Python, `EXCEPT` on, IPO off). Their build
presets use 9 jobs; their test presets run CTest with 4 jobs and fail when no
test is registered.

```bash
cmake --preset gnu && cmake --build --preset gnu            # everything
cmake --build --preset gnu --target laplace_opt             # one program
cmake -S . -B build-serial -DWITH_MPI=OFF -DWITH_OpenMP=OFF # the serial configuration CI checks
```

## Continuous integration — `.github/workflows`

| Workflow | Job | Runs |
|---|---|---|
| `ci.yml` | `gcc · MPI + OpenMP`, `clang · MPI + OpenMP` | configure with the preset, build everything, `ctest --preset` |
| `ci.yml` | `gcc · no MPI, no OpenMP` | build `build-serial`, `pytest -m unit_test` |
| `ci.yml` | `Formatting and linters` | `pre-commit` on every file: whitespace, YAML/JSON, clang-format, ruff, ruff-format, cmake-format |
| `validation.yml` | `Laplace benchmark` | weekly and on demand: builds `laplace_opt`, `droplet_opt`, `layers_opt`, runs their long tests with `--runlong`, and writes the measured values to the job summary |

## The pytest suite

Every test is a `test_*.py` under `programs/` (`pyproject.toml` sets
`testpaths = ["programs"]`). `conftest.py` puts `PyCGLBM` on the path when
`pycglbm` is not installed, and defines:

| Marker / option | Meaning |
|---|---|
| `unit_test` | fast test of one component, no simulation |
| `verification` | an internal consistency property of a run |
| `validation` | a run against a reference or a measured value |
| `long` | runs a full simulation; skipped unless `--runlong` is given |
| `--allow_plot` | keep the interactive matplotlib backend (otherwise `Agg`) |

Every solver test is `long`, so a plain `pytest` runs only the unit tests and
the quick ones.

```bash
python -m pytest -m unit_test                                   # the unit tests
python -m pytest programs/solvers/color_gradient/laplace --runlong -v
```

With `WITH_Python`, `cmake/Tests.cmake` registers one CTest test per pytest
file and marker, named `cglbm-<marker>-<dotted path>`, run from the source
directory with `CGLBM_BIN_DIR` and `CGLBM_ARTIFACTS_DIR` set. pytest's exit code
5 (nothing collected for that marker) counts as a pass.

Environment variables read by the Python helpers:

| Variable | Default | Used for |
|---|---|---|
| `CGLBM_BIN_DIR` | `<root>/bin` | where `find_program` looks |
| `CGLBM_ARTIFACTS_DIR` | `<root>/artifacts` | where `artifacts_dir()` points |
| `CGLBM_MPI_LAUNCHER` | `mpirun`, else `mpiexec` on `PATH` | the MPI launcher |
| `CGLBM_MPI_LAUNCHER_FLAGS` | `--oversubscribe` | extra launcher flags, space-separated |

## Unit-test programs — `programs/unit_testing`

Each is a small executable that prints `key = value` lines; the pytest beside
it runs it through `run_unit_program` and asserts on the values.

| Program | Checks |
|---|---|
| `lbm/gradient` (`lbm_gradient`) | `Σ W c_a c_b = δ_ab`, the isotropy defect of the rank-4, 6 and 8 tensors, and the angular error on a radial tanh profile, per stencil |
| `lbm/mixture` (`lbm_mixture`) | the volume-fraction initial state against the equation of state, the (W/2) ln(ρ₁/ρ₂) offset, the mixture viscosity, and the capillary force: σ/R on a circle, zero net force on a closed interface, zero on a flat one |
| `lbm/solver` (`lbm_solver`) | `Solver` on 32 × 32 for 20 steps: mass conservation and φ ∈ [−1, 1] in every configuration (periodic, wall, CSF, Latva-Kokko, two viscosities, capillary stress, dynamic mixing), a uniform state at rest, φ_N at its fixed points, the enhanced equilibrium, a zero net force for the stress form on a tilted ellipse, survival of 300 steps at 10³ and 10⁵ |
| `lbm/two_population` | the equilibrium's moments, the rest weights against the density ratio, each fluid's mass, non-negative populations, φ_N ∈ [−1, 1], rest |
| `lbm/two_population_3d` | the 3D stencils' isotropy, the D3Q19 rest weights, the sphere's curvature −2/R, the spheroid of `oscillation_3d`, a wall and a hydrostatic column |
| `lbm/lattice_3d` | D3Q19 and D3Q27 against their weight conditions, the equilibrium's moments including `M3_xxy`, the derived amplitude λ, the positivity bound |
| `lbm/mhd_potential` | the potential solve against a manufactured solution (uniform and 10⁴ conductivity), charge conservation, the insulating wall, both discretisations |
| `lbm/velocity_based` (`lbm_velocity_based`) | the equilibria, forcing and collisions, the carrier's positivity, the pressure force, the link exchange, the dissipation force, the hybrid collision, `set_velocity`, the phase transport and its limiter |
| `mpi/topology`, `mpi/exchange_ghosts` | the block decomposition and the halo exchange, corners included; skipped without an MPI launcher |
| `omp/environment` | thread count and a parallel region |
| `cuda/environment` | device query; skipped without CUDA |
| `postprocessing/test_unit_test_files.py` | `pycglbm.files` on synthetic output (no executable) |

## `pycglbm`

`PyCGLBM/pycglbm`, installable with `pip install -e PyCGLBM`. It depends on
numpy; matplotlib only for plotting.

### `pycglbm.files`

| Name | Behaviour |
|---|---|
| `field_path(rundir, field, t)` | `rundir/<field>_<t>.csv` |
| `load_field(rundir, field, t)` | an `[y, x]` array; `FileNotFoundError` if missing |
| `load_velocity(rundir, t)` | a `[y, x, 2]` array |
| `CaseOutput(rundir)` | one run directory; `NotADirectoryError` if it does not exist |
| `.config` | the `key = value` lines of `run.log` as a dict (empty without a log) |
| `.parameter(key, cast=float)` | one of them, converted; `KeyError` naming the log when missing |
| `.timesteps`, `.last_timestep` | sorted timesteps for which all four files exist |
| `.shape` | `(ny, nx)` |
| `.density(t)`, `.phase(t)`, `.pressure(t)`, `.velocity(t)`, `.fields(t)` | the arrays |
| `.phase_interface_radius(t)` | first zero crossing of φ along the +x ray from `(nx//2, ny//2)`, linearly interpolated; NaN if none |
| `.density_interface_radius(t)` | first crossing of the midpoint between the density at the centre and at the edge, along the same ray |
| `.droplet_radius(t)` | `√(area(φ > 0)/π)` |
| `.pressure_jump(t, inner, outer)` | mean p within `inner` of `(nx/2, ny/2)` minus mean p beyond `outer` |
| `.droplet_axes()` | the 3D track of `interface.csv` as an `[n, 4]` array; `ValueError` for the 2D track |

### `pycglbm.testing`

| Name | Behaviour |
|---|---|
| `bin_dir()`, `artifacts_dir()` | from the environment variables above |
| `find_program(name, variant="opt")` | the first `<name>_<variant>` under `bin_dir()`; `FileNotFoundError` with the build command otherwise |
| `run_program(name, rundir, variant="opt", clean=True, timeout=None, env=None, args=(), nprocs=None)` | empties and creates `rundir`, runs the program there with stdout and stderr to `run.log`, raises on a non-zero exit, returns a `CaseOutput` |
| `run_unit_program(name, args=(), nprocs=None, variant="opt", timeout=120, check=True, env=None)` | runs a unit-test program and returns the `CompletedProcess` with its captured output |
| `parse_key_values(output)` | the `key = value` lines as a dict; lines whose key contains a space are ignored |
| `mpi_launcher()`, `launch_command(executable, args, nprocs)` | the argv, wrapped in the launcher when `nprocs` is given |

### `pycglbm.oscillation`

`fit_damped_oscillation(signal, dt=1.0, refinements=6)` fits
`A e^(−αt) cos(ωt + φ) + c` by Prony's method followed by a shrinking grid
search over (α, ω), and returns a `DampedOscillation` (`angular_frequency`,
`decay_rate`, `amplitude`, `offset`, `residual`, and the properties `period`,
`undamped_frequency = √(ω² + α²)`, `damping_ratio`). It raises `ValueError` for
fewer than four samples, non-finite values, or a signal that does not
oscillate. `lamb_frequency(mode, σ, R, ρ_in, ρ_out)` returns Lamb's
3D frequency `√(n(n − 1)(n + 1)(n + 2)σ / (R³((n + 1)ρ_in + nρ_out)))`.

### Command line and plotting

```
pycglbm describe <rundir>
pycglbm plot <rundir> [-t T] [-o file.png]     # four-panel snapshot, the last timestep by default
pycglbm diff <rundir> <first> <second> [-o file.png]
```

`pycglbm.plotter.plot_fields(case, t)` and `plot_difference(case, t1, t2)`
return a matplotlib `(figure, axes)` pair.

## Utilities

| Script | Behaviour |
|---|---|
| `utilities/run_case.sh <program> [opt\|dbg] [rundir]` | finds `<program>_<variant>` under `CGLBM_BIN_DIR` (default `bin/`), creates `rundir` (default `artifacts/<program>`), runs the program there and tees its output to `run.log`. It passes no arguments to the program; to give flags, run the executable yourself inside a run directory |
| `utilities/clean_artifacts.sh` | asks, then deletes everything under `CGLBM_ARTIFACTS_DIR` (default `artifacts/`) except `.gitkeep` |
