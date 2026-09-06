"""Helpers for the pytest suite living next to the programs."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from pycglbm.files import CaseOutput

#: Repository root, i.e. the directory holding ``programs`` and ``src``.
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def bin_dir() -> Path:
    """Where CMake put the executables (``BIN_DIR``, default ``<root>/bin``)."""
    return Path(os.environ.get("CGLBM_BIN_DIR", PROJECT_ROOT / "bin"))


def artifacts_dir() -> Path:
    """Where runs are written (``ARTIFACTS_DIR``, default ``<root>/artifacts``)."""
    return Path(os.environ.get("CGLBM_ARTIFACTS_DIR", PROJECT_ROOT / "artifacts"))


def find_program(name: str, variant: str = "opt") -> Path:
    """Locate the executable ``<name>_<variant>`` anywhere under :func:`bin_dir`.

    ``bin/`` mirrors the ``programs/`` tree, so the program is found by name
    rather than by hard-coding its group.
    """
    matches = sorted(bin_dir().rglob(f"{name}_{variant}"))
    if not matches:
        raise FileNotFoundError(
            f"No executable {name}_{variant} under {bin_dir()}. Build it first, "
            f"e.g. `cmake --build --preset gnu --target {name}_{variant}`."
        )
    return matches[0]


def run_program(
    name: str,
    rundir: Path | str,
    variant: str = "opt",
    clean: bool = True,
    timeout: float | None = None,
    env: dict[str, str] | None = None,
) -> CaseOutput:
    """Run one solver inside ``rundir`` and return its output.

    The programs write their CSV files into the current working directory with
    fixed names, so each run needs a directory of its own.
    """
    executable = find_program(name, variant)
    rundir = Path(rundir)
    if clean and rundir.exists():
        shutil.rmtree(rundir)
    rundir.mkdir(parents=True, exist_ok=True)

    with (rundir / "run.log").open("w") as logfile:
        subprocess.run(
            [str(executable)],
            cwd=rundir,
            stdout=logfile,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=timeout,
            env={**os.environ, **(env or {})},
        )
    return CaseOutput(rundir)
