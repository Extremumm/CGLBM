"""PyCGLBM -- Python tooling for the CGLBM color-gradient lattice Boltzmann solver.

The solvers under ``programs/`` write plain CSV grids into their run directory.
This package turns such a directory into arrays and figures:

    >>> from pycglbm import CaseOutput
    >>> case = CaseOutput("artifacts/laplace")
    >>> case.timesteps[:3]
    [0, 1000, 2000]
    >>> case.pressure(30000).shape
    (128, 128)
"""

from pycglbm.files import CaseOutput, load_field, load_velocity
from pycglbm.oscillation import fit_damped_oscillation, lamb_frequency
from pycglbm.version import __version__

__all__ = [
    "CaseOutput",
    "load_field",
    "load_velocity",
    "fit_damped_oscillation",
    "lamb_frequency",
    "__version__",
]
