# Report

`report.tex` documents both colour-gradient solvers: the algorithms, the
derivations and proofs the implementations rest on, the measured results, and
the limits of each.

```bash
cd docs/report && make      # report.pdf
```

The figures are built from `figures/data/*.csv` by `figures/make_figures.py`.
Those CSVs are committed, and `make data` is what produced them — it rebuilds
`tools/report_data.cpp` against the library and re-runs the solvers, which takes
hours, so the results are kept rather than regenerated on every build.

| file | produced by | what it holds |
|---|---|---|
| `laplace_sweep.csv` | `report_data sweep` | relaxed Laplace jump and peak spurious current for both solvers, density ratios 2 to 10⁵, 1.2 × 10⁵ steps |
| `history_r1000.csv` | `report_data history` | the same measured every 4000 steps at a density ratio of 1000, which is where the short-run trap shows |
| `profiles_r1000.csv` | `report_data profiles` | phase, density and pressure along a ray through the interface |
| `currents_r1000.csv` | `report_data currents` | the relaxed velocity field of the two-population solver |

`fig_recolouring.pdf` needs no data: it plots the positivity criterion of the
segregation operator, which is analytic.
