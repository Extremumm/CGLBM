"""Reading the CSV output of a CGLBM run.

Every solver writes, once per ``interval`` steps, four files named after the
timestep: ``density_<t>.csv``, ``velocity_<t>.csv``, ``phase_<t>.csv`` and
``pressure_<t>.csv``. Each file holds one row per lattice row ``j`` (y) and one
column per node ``i`` (x), so the arrays returned here are indexed ``[y, x]``.
The velocity file stores the two components interleaved along the row, giving
``2 * Lx`` columns; :func:`load_velocity` folds them back into ``[y, x, 2]``.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

#: Scalar fields written by every solver.
SCALAR_FIELDS = ("density", "phase", "pressure")

#: All fields written by every solver.
FIELDS = (*SCALAR_FIELDS, "velocity")

_TIMESTEP_RE = re.compile(r"^(?P<field>[a-z]+)_(?P<timestep>\d+)\.csv$")


def field_path(rundir: Path | str, field: str, timestep: int) -> Path:
    """Return the path of one field file, without checking that it exists."""
    return Path(rundir) / f"{field}_{timestep}.csv"


def load_field(rundir: Path | str, field: str, timestep: int) -> np.ndarray:
    """Load one scalar field as an ``[y, x]`` array."""
    path = field_path(rundir, field, timestep)
    if not path.is_file():
        raise FileNotFoundError(f"No {field} output for timestep {timestep}: {path}")
    return np.loadtxt(path, delimiter=",")


def load_velocity(rundir: Path | str, timestep: int) -> np.ndarray:
    """Load the velocity field as a ``[y, x, 2]`` array."""
    raw = np.loadtxt(field_path(rundir, "velocity", timestep), delimiter=",")
    ny, ncols = raw.shape
    if ncols % 2:
        raise ValueError(f"Velocity file has an odd number of columns ({ncols})")
    return raw.reshape(ny, ncols // 2, 2)


class CaseOutput:
    """The output of a single run, i.e. one directory of CSV files."""

    def __init__(self, rundir: Path | str) -> None:
        self.rundir = Path(rundir)
        if not self.rundir.is_dir():
            raise NotADirectoryError(f"No such run directory: {self.rundir}")

    def __repr__(self) -> str:
        return f"CaseOutput({str(self.rundir)!r}, {len(self.timesteps)} timesteps)"

    @property
    def timesteps(self) -> list[int]:
        """Every timestep for which all four field files are present, sorted."""
        found: dict[int, set[str]] = {}
        for path in self.rundir.glob("*.csv"):
            match = _TIMESTEP_RE.match(path.name)
            if match and match["field"] in FIELDS:
                found.setdefault(int(match["timestep"]), set()).add(match["field"])
        return sorted(t for t, fields in found.items() if fields >= set(FIELDS))

    @property
    def last_timestep(self) -> int:
        timesteps = self.timesteps
        if not timesteps:
            raise FileNotFoundError(f"No complete output found in {self.rundir}")
        return timesteps[-1]

    @property
    def shape(self) -> tuple[int, int]:
        """Lattice shape ``(Ly, Lx)`` as stored, i.e. ``(rows, columns)``."""
        return self.density(self.timesteps[0]).shape

    def density(self, timestep: int) -> np.ndarray:
        return load_field(self.rundir, "density", timestep)

    def phase(self, timestep: int) -> np.ndarray:
        return load_field(self.rundir, "phase", timestep)

    def pressure(self, timestep: int) -> np.ndarray:
        return load_field(self.rundir, "pressure", timestep)

    def velocity(self, timestep: int) -> np.ndarray:
        return load_velocity(self.rundir, timestep)

    def fields(self, timestep: int) -> dict[str, np.ndarray]:
        """All four fields at one timestep, keyed by field name."""
        return {
            "density": self.density(timestep),
            "velocity": self.velocity(timestep),
            "phase": self.phase(timestep),
            "pressure": self.pressure(timestep),
        }

    def droplet_radius(self, timestep: int) -> float:
        """Effective radius of the ``phi > 0`` region, from its area."""
        area = float((self.phase(timestep) > 0).sum())
        return float(np.sqrt(area / np.pi))

    def pressure_jump(self, timestep: int, inner: float, outer: float) -> float:
        """Mean pressure within ``inner`` of the centre, minus the mean beyond ``outer``.

        Both radii are in lattice units and measured from the domain centre,
        which is where every droplet case initialises its droplet.
        """
        pressure = self.pressure(timestep)
        ny, nx = pressure.shape
        y, x = np.mgrid[0:ny, 0:nx]
        radius = np.hypot(x - nx / 2, y - ny / 2)
        return float(pressure[radius < inner].mean() - pressure[radius > outer].mean())
