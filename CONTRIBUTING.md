# Contributing

## General remarks

CGLBM is a research code: a change is worth making when it can be reproduced and
checked. Every physical change should come with the run that shows its effect,
and ideally with a test under `programs`.

### Managing branches

Work on a branch, keep `main` buildable. A branch is ready when
`cmake --build --preset gnu` succeeds for both the `opt` and the `dbg` variants
and `pytest` is green.

## Coding standards

### Linters

The formatting is enforced by [pre-commit](https://pre-commit.com):

```bash
pip install pre-commit
pre-commit install
pre-commit run --all-files
```

- C++ is formatted by `clang-format` (`.clang-format`: LLVM style, 4 spaces,
  100 columns).
- Python is linted and formatted by `ruff` (configured in `pyproject.toml`).
- CMake is formatted by `cmake-format` (`cmake-format.py`).

### Naming style

| Kind | Convention | Example |
|---|---|---|
| program entry point | `main_<target>.cpp` | `main_gravity_capillary.cpp` |
| library header | lower snake case, prefixed by its module | `src/mpi/mpi_topology.h` |
| library symbol | `cglbm::<module>` namespace | `cglbm::mpi::CartesianTopology` |
| function | lower camel case, as in the existing solvers | `calMacroscopic`, `collide_surface` |
| lattice constant | as in the paper's notation | `cs2`, `w`, `xi`, `p1_inf` |
| test file | `test_{marker}_{name}.py` | `test_unit_test_files.py` |

Physical quantities stay in lattice units inside the solvers; the conversion
factors `c_dx` and `c_dt` at the top of each program are the only place where
physical units appear. Document the unit of any new constant in a trailing
comment, as the existing ones do.

### Comments

Comment the physics, not the syntax: which term of which equation a block
implements, and why a correction is there. Reference the paper or the Krüger
book with a page number when a formula is not obvious — the existing sources do
this and it is the fastest way back into the algorithm.

## Adding to the library

`src/<module>` holds the library, one directory per concern (`core`, `lbm`,
`omp`, `mpi`). Sources there are compiled into `libcglbm_opt` / `libcglbm_dbg`
and linked into every program, so anything added must build in every supported
configuration.

An optional dependency must not become mandatory: guard it behind its macro
(`_OPENMP`, `CGLBM_WITH_MPI`) and give the disabled branch behaviour that still
makes sense — one thread, one rank — so that a solver written against the
interface needs no `#ifdef`. `WITH_OpenMP=OFF WITH_MPI=OFF` must still build and
pass the suite:

```bash
cmake -S . -B build-serial -DWITH_MPI=OFF -DWITH_OpenMP=OFF
cmake --build build-serial -j
```

Keep `mpi.h` out of the public headers — `CartesianTopology::communicator()`
returns a `void*` for that reason — so a program that only needs the
decomposition does not inherit the MPI include path.

## Adding a program

1. Create `programs/<group>/<name>/main_<name>.cpp`. CMake globs
   `programs/*/main_*.cpp`, so the targets `<name>_opt` and `<name>_dbg` appear
   with no CMake edit; every other source in that directory is compiled into it.
2. Reconfigure (`cmake --preset gnu`) — the new program is listed at configure
   time.
3. If the program needs a flag of its own, add it to `cmake/Exceptions.cmake`
   rather than widening the global flags.
4. Add a `tests/test_<name>.py` beside it.

## Testing

All tests live in `programs`, next to the code they exercise. `pytest` is the
framework, and `pycglbm.testing` provides the helpers.

### Test structure

```
programs/solvers/color_gradient/laplace/
├── main_laplace.cpp
└── tests/
    └── test_laplace_color_gradient.py
```

### Naming convention

1. **Test files** must start with `test_` and may include the marker:
   `test_{marker}_{name}.py` or `test_{name}.py`. The name should carry enough
   of the path to be unique.
2. **Test functions** must start with `test_`, include the marker, and match
   their file: `test_{marker}_{testfile_suffix}_{function_suffix}`.
3. **Markers**: every function is decorated with `@pytest.mark.{marker}`, one of
   `unit_test`, `validation` or `verification`. The decorator and the name must
   agree.
4. A test that launches a solver is additionally marked `@pytest.mark.long` and
   only runs under `--runlong`.

### Running

```bash
pytest -m unit_test          # fast checks only
pytest --runlong             # everything, simulations included
ctest --preset gnu           # through CTest: one entry per file and marker
```

CTest registers each test file once per marker; a file with no test for a marker
exits 5, which the runner maps to success.

### Writing a test that needs several ranks

`pycglbm.testing.run_unit_program(name, args, nprocs=4)` wraps the program in
`mpirun`; `parse_key_values` reads the `key = value` lines it prints. Skip the
test when there is no launcher:

```python
requires_mpi = pytest.mark.skipif(mpi_launcher() is None, reason="no MPI launcher available")
```

Parametrise over rank counts that exercise different rank grids — 1, 2, 4 and 6
give 1x1, 2x1, 2x2 and 3x2 — and include a lattice size that does not divide
evenly, since the remainder handling is where a decomposition usually breaks.

### Writing a validation test

A validation test compares a run against something external: an analytic
solution, a reference file, or a previously measured value. When the code does
not reproduce the analytic result, say so in the test: pin the measured value,
name it as measured in a comment, and record the gap in `docs/numerics.md`. A
test that quietly loosens its tolerance until it passes hides the physics.

Use `pycglbm.testing.run_program`, which gives each run its own directory under
`artifacts`, and share one run between the tests of a module with a
module-scoped fixture — the cases take minutes.
