# Output — `src/lbm/output_writer.h`

Every program writes into its current working directory under fixed names, so
each run needs a directory of its own (`utilities/run_case.sh` makes one under
`artifacts/`). A failure to open or write a file throws `OutputError`, a
`std::runtime_error`, rather than leaving a silently empty run.

## Files a run produces

| File | Written by | When | Contents |
|---|---|---|---|
| `density_<t>.csv` | every solver | t = 0 and every `interval` steps | ρ, `ny` rows × `nx` columns |
| `velocity_<t>.csv` | every solver | same | `ny` rows × `2·nx` columns: `u_x,u_y` interleaved per node |
| `phase_<t>.csv` | every solver | same | φ (`Solver`), φ_N (two-population solvers), ψ = 2c − 1 (velocity-based) |
| `pressure_<t>.csv` | every solver | same | p (velocity-based: relative to the far field) |
| `interface.csv` | `Solver` and `TwoPopulationSolver3D` with `track_interface` | every step | see below |
| `field_<t>.vtk` | `TwoPopulationSolver3D` with `write_vtk_field` (`--vtk`) | t = 0 and every `interval` steps | the whole 3D field |
| `run.log` | the shell or `pycglbm.testing.run_program` | — | the program's stdout: `describe(config)`, then `Step <t>` lines |

Row `j` of a CSV file is lattice row `j` (y), column `i` is node `i` (x), so
`numpy.loadtxt(..., delimiter=",")` gives `array[y, x]`. A 3D run writes the
`z = nz/2` slice into the same four files, velocity with its x and y
components only. Values are written with `output_precision` significant digits
(6 by default; `--precision=17` round-trips a double). `droplet` and `layers`
write their own CSV files with precision 10.

`interface.csv` holds one of two tracks, identified by its header:

| Header | Written by | Lines |
|---|---|---|
| `Timestep, phi, y` | `CsvWriter::write_interface` (2D) | one per node of the column `i = nx/2` whose φ is strictly inside (−1, 1), every step |
| `Timestep, rx, ry, rz` | `CsvWriter::write_axes` (3D) | one per step: the distances from the centre to the φ_N = 0 crossing along +x, +y, +z; `nan` where there is none |

## `MacroscopicState`

```cpp
struct MacroscopicState {
    const Field* density;    // one value per node
    const Field* velocity;   // depth 2
    const Field* phase;
    const Field* pressure;
};
```

The 2D solvers' `state()` returns one; `TwoPopulationSolver3D::state()` returns
the 3D analogue `MacroscopicState3D` (in `two_population_solver_3d.h`).

## `CsvWriter`

```cpp
class CsvWriter {
public:
    explicit CsvWriter(int precision);
    void write_grids(int timestep, const MacroscopicState&) const;  // the four CSV files
    void open_interface_track();                  // interface.csv, header "Timestep, phi, y"
    void write_interface(int timestep, const Field& phase);
    void open_droplet_track();                    // interface.csv, header "Timestep, rx, ry, rz"
    void write_axes(int timestep, const double* radii);
};
```

`write_grids` opens all four files, writes them, and checks each on close.
The track methods write through one `std::ofstream` member kept open for the
run and check the stream after every line. No run opens both tracks.

## VTK

| Function | File | Format |
|---|---|---|
| `write_vtk(prefix, t, state)` | `<prefix><t>.vtk` (no separator) | legacy ASCII `STRUCTURED_POINTS`, `nx × ny × 1`; `density`, `velocity` (z = 0), `phase_field`, `pressure` as point data |
| `write_vtk_3d(prefix, t, ρ, u, φ, p)` | `<prefix>_<t>.vtk` | legacy ASCII `STRUCTURED_POINTS`, `nx × ny × nz`, scalars and the three-component velocity as `double`; declared in `output_writer.h`, defined in `two_population_solver_3d.cpp` |

No 2D program calls `write_vtk`. `TwoPopulationSolver3D::run()` calls
`write_vtk_3d("field", …)` when `write_vtk_field` is set, giving
`field_<t>.vtk`.

## Reading it back

`pycglbm.CaseOutput(rundir)` reads all of the above; see
[Build and tests](build-and-tests.md#pycglbm).
