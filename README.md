# CGLBM

Color-gradient lattice Boltzmann method in the D2Q9 scheme, for two-phase
flows. The algorithm is the improved color-gradient model of T. Lafarge,
P. Boivin, N. Odier and B. Cuenot, *Improved color-gradient method for lattice
Boltzmann modeling of two-phase flows*, Physics of Fluids 33 (8), 082110, 2021
([10.1063/5.0061638](https://doi.org/10.1063/5.0061638), `hal-03324224`).

The project is laid out as follows:
a library under `src`, the executables and their tests under `programs`, the
build configuration under `cmake`, the Python tooling in `PyCGLBM`.

## Installation

### System dependencies

- a C++17 compiler (GCC ≥ 9 or Clang ≥ 10),
- CMake ≥ 3.23 (3.16 is enough without the presets),
- OpenMP, for `src/omp` and the `rayleigh_taylor_omp` program,
- an MPI implementation (OpenMPI, MPICH), for `src/mpi`,
- CUDA Toolkit ≥ 11.0 (optional), for GPU acceleration with `src/cuda` and `rayleigh_taylor_cuda`,
- Python ≥ 3.9 for the post-processing and the test suite.

```bash
pip install -r requirements.txt   # numpy, matplotlib, pytest
pip install -e PyCGLBM            # the pycglbm package and its CLI
```

### Compilation

The presets in `cmake/presets` configure a compiler each. Both an optimized and
a debug version of every program are built, suffixed `_opt` and `_dbg`:

```bash
cmake --preset gnu          # or: --preset clang
cmake --build --preset gnu
```

Useful targets:

| Target | Builds |
|---|---|
| `opt` | every optimized program |
| `dbg` | every debug program |
| `laplace_all` | `laplace_opt` and `laplace_dbg` |
| `laplace_opt` | that one program |
| `rayleigh_taylor_cuda_opt` | CUDA-accelerated Rayleigh-Taylor simulation |

Configuration is driven by the options printed at configure time — `RELEASE`,
`DEBUG`, `ARCH`, `WITH_OpenMP`, `WITH_MPI`, `WITH_CUDA`, `WITH_IPO`, `WITH_Python`, `EXCEPT`,
`BIN_DIR`, `ARTIFACTS_DIR`. `WITH_OpenMP=OFF`, `WITH_MPI=OFF`, and `WITH_CUDA=OFF` still build
everything: the parallel modules fall back to one thread, one rank, or host execution. Set them on the command line (`-DDEBUG=OFF`), through the
environment, or in a `CMakeUserPresets.json` — see
`cmake/CMakeUserPresets.example.json`.

To change the compiler flags, edit `cmake/CompilerFlags.cmake`; for a flag that
applies to a single program, `cmake/Exceptions.cmake`.

## User's guide

### Run your first case

Each program writes its output to the current working directory under fixed
names, so every run needs a directory of its own. `utilities/run_case.sh` takes
care of that:

```bash
utilities/run_case.sh laplace          # runs bin/.../laplace_opt in artifacts/laplace
pycglbm describe artifacts/laplace     # what the run produced
pycglbm plot artifacts/laplace -o snapshot.png
```

Equivalently, by hand:

```bash
mkdir -p artifacts/laplace && cd artifacts/laplace
../../bin/solvers/color_gradient/laplace/laplace_opt | tee run.log
```

Every solver takes its case as built-in defaults and accepts overrides on the
command line, so a resolution or a run length is a flag rather than a rebuild:

```
  --stencil=E4|E6|E8  isotropy order of the colour gradient
  --initial-profile=colour|normalised
                      field the initial tanh profile is prescribed in
  --interface-field=colour|normalised
                      field the colour gradient is taken of
  --surface-tension=perturbation|csf
                      capillary stress in Omega^(2), or a body force from an
                      explicit curvature (Ba et al. 2016)
  --recolouring=width|latva-kokko
                      segregation strength from the pressure and an interface
                      width, or from the density and beta
  --beta=X            segregation strength of `latva-kokko`, in (0, 1]
  --nu=X, --nu-b=X    kinematic shear and bulk viscosity of component 1
  --nu2=X, --nu-b2=X  the same for component 2; unset means equal to
                      component 1's
  --initial-state=equilibrium|eos|linear
                      how rho and p are laid down at t = 0
  --nx=N, --ny=N      lattice size
  --steps=N           number of time steps
  --interval=N        write the CSV grids every N steps
  --precision=N       significant digits in the CSV output (default 6)
  --threads=N         OpenMP threads; implies parallel execution
  --help              the same list
```

```bash
utilities/run_case.sh laplace                    # the case as shipped
bin/solvers/color_gradient/laplace/laplace_opt E4        # the original stencil
bin/solvers/color_gradient/laplace/laplace_opt --nx=256 --steps=5000
```

The bare stencil name is still accepted as the first argument, as before.

The run log opens with the full configuration as `key = value` lines, which is
what `pycglbm` reads back:

```python
CaseOutput("artifacts/laplace").parameter("sigma")
```

The physical constants of a case — densities, viscosities, surface tension, the
unit conversion factors `c_dx` and `c_dt` — are in the `*_case()` function at
the top of its `main_*.cpp`. Set the thread count of the OpenMP program with
`OMP_NUM_THREADS` or `--threads`; without either it falls back to eight.

### About the generated files

Every program writes four ASCII CSV grids per `interval` steps, named after the
timestep:

| File | Shape (rows × columns) | Contents |
|---|---|---|
| `density_<t>.csv` | Ly × Lx | mixture density ρ |
| `velocity_<t>.csv` | Ly × 2·Lx | velocity, (uₓ, u_y) interleaved per node |
| `phase_<t>.csv` | Ly × Lx | phase field φ ∈ [−1, 1] |
| `pressure_<t>.csv` | Ly × Lx | pressure p |

One row per lattice row `j` (y), one column per node `i` (x), so `numpy` reads
them as `array[y, x]`. `capillary` and `gravity_capillary` also append the
interface position to `interface.csv`. `cglbm::lbm::write_vtk` produces a
ParaView/VisIt file from the same state; no case calls it by default.

Six significant digits is the default and loses about ten digits of a double.
Pass `--precision=17` for output that round-trips.

### Programs

| Program | Lattice | Steps | ρ₁/ρ₂ | Gravity | Purpose |
|---|---|---|---|---|---|
| `laplace` | 128×128 | 3×10⁴ | 20/1 | no | Laplace law Δp = σ/R across a static droplet |
| `capillary` | 128×128 | 10⁴ | 4/1 | no | oscillation period of a perturbed droplet |
| `gravity_capillary` | 128×128 | 5×10⁴ | 4/1 | yes | droplet under gravity and surface tension |
| `rayleigh_taylor` | 128×1028 | 5×10⁶ | 4/1 | yes | Rayleigh–Taylor instability, σ = 0, serial |
| `rayleigh_taylor_omp` | 1024×4096 | 2×10⁶ | 4/1 | yes | same case, OpenMP, production resolution |
| `laplace_high_ratio` | 100×100 | 4×10⁴ | 1000/1 | no | Laplace law at a density ratio of 1000, two-population solver |

> `rayleigh_taylor_omp` allocates several GB of lattice at its production
> resolution. Check the available memory before launching it, or lower it with
> `--nx` and `--ny`.

### Organization

#### The `bin` folder

Contains all the executables, following the directory structure of `programs`,
each in an optimized (`_opt`) and a debug (`_dbg`) version.

#### The `programs` folder

Contains all the programs

- `programs/initial_conditions` scripts producing or checking initial conditions
- `programs/solvers` the solvers, grouped by method, with their associated pytests
- `programs/unit_testing` tests for small portions of the tooling

#### The `src` folder

Contains all the sources for the library

- `src/lbm` the schemes: the D2Q9 lattice, the two solvers and their time loops, the
  two-component equation of state, the isotropic colour-gradient stencils, the
  lattice storage and the output writers — each carrying its reference
- `src/omp` thread-level parallelism: thread count, ids, timing
- `src/mpi` distributed memory: environment, Cartesian decomposition, halo
  exchange, error checking
- `src/cuda` GPU support: device query, error checking, device memory

The time loop lives in `cglbm::lbm::Solver` (`src/lbm/solver.cpp`), and a
program under `programs/solvers` is a case definition handed to it. See
[`docs/numerics.md`](docs/numerics.md) and [`docs/parallel.md`](docs/parallel.md).

#### The `artifacts` folder

Where runs and tests write their output. Git-ignored; clean it with
`utilities/clean_artifacts.sh`.

### Structure of a program

A program is a directory under `programs` holding one `main_<name>.cpp`. CMake
discovers it automatically and builds the target `<name>_opt` / `<name>_dbg`
from every source in that directory, so adding a case means adding a directory —
no CMake edit. Its tests are `test_*.py` files placed beside it, usually in a
`tests` subdirectory.

A solver program fills a `cglbm::lbm::CaseConfig` and hands it to
`cglbm::lbm::Solver`; the scheme itself is in the library and is not repeated
per case:

```cpp
CaseConfig config;
config.nx = 128;
config.ny = 128;
config.steps = 10000;
config.physics.rho1 = 4.;
config.physics.rho2 = 1.;
config.boundary = Boundary::WallY;
config.initial_phase = cosine_layer(0.2, /*inverted=*/true);
```

## Testing

```bash
pytest -m unit_test              # fast, no simulation
pytest --runlong                 # includes the tests that run a full case
ctest --preset gnu               # the same suite, through CTest
```

The `src/omp` and `src/mpi` tests are part of `-m unit_test`: they run their
programs directly, the MPI ones on 1 to 6 ranks through `mpirun`. They skip
themselves when the dependency is missing — no launcher, or a build configured
with `WITH_MPI=OFF` / `WITH_OpenMP=OFF` — so the same suite runs against any
configuration.

### Continuous integration

| Workflow | Trigger | What it does |
|---|---|---|
| `.github/workflows/ci.yml` | push to `main`, pull requests, manual | builds with GCC and with Clang (MPI + OpenMP) and runs the CTest suite; builds again with both backends off and runs the tests that still apply; runs the formatters |
| `.github/workflows/validation.yml` | Mondays 04:00 UTC, manual | runs the long Laplace tests and writes the measured jump, interface radii and spurious currents into the job summary |

The formatters are the ones in `.pre-commit-config.yaml`, so
`pre-commit run --all-files` locally is the same check CI runs.

Tests are named `test_{marker}_{name}.py`, and each function carries the
matching marker: `unit_test`, `validation` or `verification`. Tests that launch
a solver are additionally marked `long` and skipped unless `--runlong` is given.
`pycglbm.testing.run_program` runs a program in its own directory and returns
the resulting `CaseOutput`.

> **Known gap.** The `laplace` case relaxes to a stationary pressure jump of
> about 1.02 σ/R instead of σ/R, measured against the radius the density field
> settles at, which is the surface where the two components occupy equal volume.
> The jump is exact at t = 0 by construction. The validation test pins this
> measured value to catch regressions; it does not certify the Laplace law.
> See [`docs/numerics.md`](docs/numerics.md).

### Density ratio

There are two colour-gradient solvers here, and which one you want depends
entirely on the density ratio.

**`cglbm::lbm::Solver`** (`laplace`, `capillary`, `gravity_capillary`,
`rayleigh_taylor`) is the improved method of Lafarge et al.: the density ratio
comes from a two-component equation of state, so it is independent of the
sound-speed ratio and the model is fully compressible. It is quantitative — a
couple of per cent on the Laplace jump — and steady up to a density ratio of
about **500**. Two things carry it there, and both come from
[Ba et al. (2016)](docs/references.md):

- the interface is located by the bulk-normalised phase field φ_N = 2c − 1,
  whose zero contour is the surface where the two components occupy equal
  volume, rather than by the colour field, whose zero sits at φ = 0.998 at a
  ratio of 10³ — well inside the light fluid;
- the surface tension is applied as a body force built from an explicit
  curvature, rather than as a capillary stress inside the collision, which
  needs the relaxation time to be uniform across the interface and it is not.

On the shipped `laplace` case that moved the jump from 0.895 to 1.021 σ/R and
cut the spurious currents by a factor of 72, to 1.7 × 10⁻⁵. The enhanced
equilibrium of Leclaire et al. (2013) turned out to be already present — it is
the same expression as the third-order Hermite term the scheme already had,
agreeing to 10⁻¹⁴.

**At 10³ and above it diverges.** Be careful reading short runs here: a density
ratio of 10³ looks healthy for its first 5 × 10⁴ steps and dies at 8.7 × 10⁴,
which is how this repository came to claim 10³ and, before that, 10⁵. Everything
quoted now is from 1.2 × 10⁵ steps with the whole time series checked.

**`cglbm::lbm::TwoPopulationSolver`** (`laplace_high_ratio`) is the classical
model of Grunau, Reis & Phillips, Leclaire et al. and Ba et al.: one
distribution per fluid, and the density ratio carried by the equilibrium's
rest-particle weight. The two bulk pressures then match identically, so the
pressure is continuous across the interface at any ratio, and each fluid's
density has its own smooth profile. On Ba et al.'s own static-droplet
benchmark it gives

| density ratio | σ measured / σ | spurious currents | Ba et al. |
|---|---|---|---|
| 100 | 1.0024 | 2.9 × 10⁻⁵ | 1.0069, 6.8 × 10⁻⁵ |
| 1000 | 1.0030 | 4.0 × 10⁻⁵ | 1.0074, 1.25 × 10⁻⁴ |

converged, with the last half of the run flat to every digit. The price is the
limitation the first model exists to avoid: the density ratio and the
sound-speed ratio are tied together, so at 10³ the heavy fluid's sound speed is
0.022 in lattice units. Use it for static or slow flows at high contrast, and
`Solver` for anything acoustic below a few hundred.

The measurements, the derivations and what is still missing are in
[`docs/numerics.md`](docs/numerics.md#how-far-the-density-ratio-goes).

## PyCGLBM

PyCGLBM is the Python library for CGLBM: reading a run directory, plotting it,
and the helpers used by the test suite. See [`PyCGLBM/README.md`](PyCGLBM/README.md).

```python
from pycglbm import CaseOutput

case = CaseOutput("artifacts/laplace")
case.pressure_jump(30000, inner=5, outer=30)
```

## Documentation

- [`docs/report/report.tex`](docs/report/) — the reference document: user
  guide, command-line and configuration reference, the derivation of both
  schemes with proofs, and the validation results with figures. `make` in that
  directory builds `report.pdf`
- [`docs/numerics.md`](docs/numerics.md) — the discretisation, the collision
  operators, the time loop, and where each step lives in the code
- [`docs/parallel.md`](docs/parallel.md) — the OpenMP and MPI modules
- [`docs/references.md`](docs/references.md) — bibliography and unit conversion
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — coding standards and workflow
- [`AUTHORS.md`](AUTHORS.md) — credits

## License

GNU General Public License v3 — see [LICENSE](LICENSE).
