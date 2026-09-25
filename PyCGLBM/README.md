# PyCGLBM

PyCGLBM is the Python library for CGLBM. It reads the CSV output of a run,
provides the figures used for post-processing, and supplies the helpers the
pytest suite uses to launch solvers.

## Install

```bash
pip install -e PyCGLBM
```

## Library

```python
from pycglbm import CaseOutput

case = CaseOutput("artifacts/laplace")
case.timesteps  # every timestep with a complete set of fields
case.pressure(30000)  # [y, x] array
case.velocity(30000)  # [y, x, 2] array
case.droplet_radius(30000)
case.pressure_jump(30000, inner=5, outer=30)
case.droplet_axes()  # [step, 4] track of a three-dimensional droplet
```

`pycglbm.plotter` builds the four-panel figures (`plot_fields`,
`plot_difference`), and `pycglbm.testing` provides `run_program`, which runs a
solver inside a dedicated directory and hands back its `CaseOutput`.

`pycglbm.oscillation` measures a ringing signal:

```python
from pycglbm import fit_damped_oscillation, lamb_frequency

axes = case.droplet_axes()
fit = fit_damped_oscillation(axes[50:, 3] - axes[50:, 1])
fit.period, fit.undamped_frequency, fit.damping_ratio
lamb_frequency(mode=2, sigma=0.1, radius=10.0, rho_in=10.0, rho_out=1.0)
```

`fit_damped_oscillation` fits `A exp(-alpha t) cos(omega t + phi) + c` with
numpy alone -- Prony's recurrence for the starting point, then a shrinking grid
over the two nonlinear parameters -- and reports what it could not explain, so a
signal that is not one decaying mode can be told apart from one that is.

## Command line

```bash
pycglbm describe artifacts/laplace
pycglbm plot artifacts/laplace -t 30000 -o snapshot.png
pycglbm diff artifacts/laplace 10000 30000
```

## Examples

`code_examples/` holds the original interactive scripts (`plot.py`,
`animation.py`, `difference.py`). They now load their data through `pycglbm`
and are run from inside a run directory.
