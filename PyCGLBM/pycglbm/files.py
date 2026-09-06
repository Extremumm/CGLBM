"""Reading the CSV output of a CGLBM run.

Every solver writes, once per ``interval`` steps, four files named after the
timestep: ``density_<t>.csv``, ``velocity_<t>.csv``, ``phase_<t>.csv`` and
``pressure_<t>.csv``. Each file holds one row per lattice row ``j`` (y) and one
column per node ``i`` (x), so the arrays returned here are indexed ``[y, x]``.
The velocity file stores the two components interleaved along the row, giving
``2 * Lx`` columns; :func:`load_velocity` folds them back into ``[y, x, 2]``.

A run also leaves a ``run.log`` whose header is the case configuration, one
``key = value`` per line. :attr:`CaseOutput.config` reads it back, so a test
states a parameter by asking the run rather than by repeating a constant that
lives in the solver.
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
    def config(self) -> dict[str, str]:
        """The case configuration the solver reported, read back from ``run.log``.

        Returns an empty mapping when there is no log, which is what a run
        launched by hand without redirecting its output leaves behind.
        """
        logfile = self.rundir / "run.log"
        if not logfile.is_file():
            return {}
        values: dict[str, str] = {}
        for line in logfile.read_text(errors="replace").splitlines():
            key, separator, value = line.partition("=")
            if separator and key.strip() and " " not in key.strip():
                values[key.strip()] = value.strip()
        return values

    def parameter(self, key: str, cast=float):
        """One reported parameter, converted with ``cast``.

        Raises ``KeyError`` naming the run when the solver did not report it,
        which is the useful failure: it means the log is from an older binary,
        not that the value is zero.
        """
        config = self.config
        if key not in config:
            raise KeyError(f"{key!r} not reported in {self.rundir / 'run.log'}")
        return cast(config[key])

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

    def _centreline(self, field: np.ndarray) -> np.ndarray:
        """The outward ray from the domain centre along +x, as a 1-D profile."""
        ny, nx = field.shape
        return field[ny // 2, nx // 2 :]

    @staticmethod
    def _first_crossing(profile: np.ndarray, level: float) -> float:
        """Where ``profile`` first falls through ``level``, linearly interpolated.

        Returns NaN when it never does.
        """
        shifted = profile - level
        for k in range(len(shifted) - 1):
            if shifted[k] >= 0.0 > shifted[k + 1]:
                return float(k + shifted[k] / (shifted[k] - shifted[k + 1]))
        return float("nan")

    def phase_interface_radius(self, timestep: int) -> float:
        """Radius where the phase field crosses zero, along the +x centreline."""
        return self._first_crossing(self._centreline(self.phase(timestep)), 0.0)

    def density_interface_radius(self, timestep: int) -> float:
        """Radius where the density crosses halfway between its two bulk values.

        For a droplet this is the interface as the *density* field sees it. It
        need not coincide with :meth:`phase_interface_radius`: the two can
        settle apart, and which one Laplace's law should be scored against
        matters when they do.
        """
        profile = self._centreline(self.density(timestep))
        midpoint = 0.5 * (float(profile[0]) + float(profile[-1]))
        return self._first_crossing(profile, midpoint)

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
