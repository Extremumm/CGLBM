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

Configuration is driven by the options printed at configure time — `RELEASE`,
`DEBUG`, `ARCH`, `WITH_OpenMP`, `WITH_MPI`, `WITH_IPO`, `WITH_Python`, `EXCEPT`,
`BIN_DIR`, `ARTIFACTS_DIR`. `WITH_OpenMP=OFF` and `WITH_MPI=OFF` still build
everything: the parallel modules fall back to one thread and one rank. Set them on the command line (`-DDEBUG=OFF`), through the
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

To change a configuration — lattice size, viscosities, surface tension, number
of steps — edit the constants block at the top of the program's `main_*.cpp`
and rebuild. Set the thread count of the OpenMP program with `OMP_NUM_THREADS`;
without it that case falls back to eight threads.

Every solver takes one optional argument, the isotropy order of the colour
gradient — `E4`, `E6` or `E8`:

```bash
utilities/run_case.sh laplace            # the program's own default
../../bin/solvers/color_gradient/laplace/laplace_opt E4   # the original stencil
```

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
interface position to `interface.csv`. A VTK writer (`outputVTK`) exists in
every program for ParaView/VisIt, commented out of the time loop by default.

### Programs

| Program | Lattice | Steps | ρ₁/ρ₂ | Gravity | Purpose |
|---|---|---|---|---|---|
| `laplace` | 128×128 | 3×10⁴ | 20/1 | no | Laplace law Δp = σ/R across a static droplet |
| `capillary` | 128×128 | 10⁴ | 4/1 | no | oscillation period of a perturbed droplet |
| `gravity_capillary` | 128×128 | 5×10⁴ | 4/1 | yes | droplet under gravity and surface tension |
| `rayleigh_taylor` | 128×1028 | 5×10⁶ | 4/1 | yes | Rayleigh–Taylor instability, σ = 0, serial |
| `rayleigh_taylor_omp` | 1024×4096 | 2×10⁶ | 4/1 | yes | same case, OpenMP, production resolution |

> `rayleigh_taylor_omp` declares about 2 GB of static lattice arrays, which is
> why `cmake/Exceptions.cmake` builds it with `-mcmodel=medium` on x86-64.
> Check the available memory before launching it.

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

- `src/core` the constants
- `src/lbm` the scheme itself: the two-component equation of state and the
  isotropic colour-gradient stencils, each carrying its reference
- `src/omp` thread-level parallelism: thread count, ids, timing
- `src/mpi` distributed memory: environment, Cartesian decomposition, halo
  exchange, error checking

The time loop still lives inside each program's `main_*.cpp`, and
`src/main_cglbm.cpp` is the entry point waiting for it to move into the library.
`src/lbm`, `src/omp` and `src/mpi` are complete and tested; see
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
> about 0.72 σ/R instead of σ/R, while the droplet radius and the interface stay
> stable. The jump is exact at t = 0 by construction. The validation test pins
> this measured value to catch regressions; it does not certify the Laplace law.
> See [`docs/numerics.md`](docs/numerics.md).

## PyCGLBM

PyCGLBM is the Python library for CGLBM: reading a run directory, plotting it,
and the helpers used by the test suite. See [`PyCGLBM/README.md`](PyCGLBM/README.md).

```python
from pycglbm import CaseOutput

case = CaseOutput("artifacts/laplace")
case.pressure_jump(30000, inner=5, outer=30)
```

## Documentation

- [`docs/numerics.md`](docs/numerics.md) — the discretisation, the collision
  operators, the time loop, and where each step lives in the code
- [`docs/parallel.md`](docs/parallel.md) — the OpenMP and MPI modules
- [`docs/references.md`](docs/references.md) — bibliography and unit conversion
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — coding standards and workflow
- [`AUTHORS.md`](AUTHORS.md) — credits

## License

GNU General Public License v3 — see [LICENSE](LICENSE).
