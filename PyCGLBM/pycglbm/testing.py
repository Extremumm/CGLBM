"""Helpers for the pytest suite living next to the programs."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from pycglbm.files import CaseOutput

#: Default flags handed to the MPI launcher. Oversubscription keeps a test that
#: asks for more ranks than the machine has cores from failing to start.
MPI_LAUNCHER_FLAGS = ("--oversubscribe",)

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


def mpi_launcher() -> str | None:
    """The MPI launcher to use, or None when MPI is unavailable.

    ``CGLBM_MPI_LAUNCHER`` overrides the lookup; CMake sets it from
    ``MPIEXEC_EXECUTABLE`` when it configures the test suite.
    """
    from shutil import which

    launcher = os.environ.get("CGLBM_MPI_LAUNCHER") or which("mpirun") or which("mpiexec")
    return launcher or None


def launch_command(executable: Path, args=(), nprocs: int | None = None) -> list[str]:
    """Build the argv for one run, wrapping it in the MPI launcher if asked."""
    command = [str(executable), *(str(a) for a in args)]
    if nprocs is None:
        return command

    launcher = mpi_launcher()
    if launcher is None:
        raise RuntimeError("No MPI launcher found; cannot run on several ranks.")
    flags = os.environ.get("CGLBM_MPI_LAUNCHER_FLAGS")
    extra = flags.split() if flags is not None else list(MPI_LAUNCHER_FLAGS)
    return [launcher, "-n", str(nprocs), *extra, *command]


def run_unit_program(
    name: str,
    args=(),
    nprocs: int | None = None,
    variant: str = "opt",
    timeout: float | None = 120.0,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    """Run a unit-test program and capture its output.

    Unlike :func:`run_program` these write nothing to disk: they report through
    stdout, which :func:`parse_key_values` turns into a dictionary.
    """
    command = launch_command(find_program(name, variant), args, nprocs)
    return subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=check,
        env={**os.environ, **(env or {})},
    )


def parse_key_values(output: str) -> dict[str, str]:
    """Read the ``key = value`` lines a unit-test program prints.

    Lines without a ``=`` are ignored, so a program may print anything else
    alongside. A repeated key keeps its last value.
    """
    values: dict[str, str] = {}
    for line in output.splitlines():
        key, separator, value = line.partition("=")
        if separator and key.strip() and " " not in key.strip():
            values[key.strip()] = value.strip()
    return values


def run_program(
    name: str,
    rundir: Path | str,
    variant: str = "opt",
    clean: bool = True,
    timeout: float | None = None,
    env: dict[str, str] | None = None,
    args=(),
    nprocs: int | None = None,
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
            launch_command(executable, args, nprocs),
            cwd=rundir,
            stdout=logfile,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=timeout,
            env={**os.environ, **(env or {})},
        )
    return CaseOutput(rundir)
