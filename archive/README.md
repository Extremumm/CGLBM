# Archive

Superseded material, kept for provenance. Nothing here is built or maintained.

- **`CG_D2Q9.h`** — the first draft of a combined constants/declarations header.
  It does not compile (several `const` objects are declared without an
  initialiser, and `f` is declared twice). It was replaced by
  `src/core/constants.h` + `src/lbm/lattice_boltzmann.h`, and by the
  self-contained programs under `programs/solvers`.
