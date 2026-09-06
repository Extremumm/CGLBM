"""Command line interface: ``pycglbm <command> ...``."""

from __future__ import annotations

import argparse
from pathlib import Path

from pycglbm.files import CaseOutput
from pycglbm.version import __version__


def _resolve(case: CaseOutput, timestep: int | None) -> int:
    return case.last_timestep if timestep is None else timestep


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pycglbm", description="Inspect a CGLBM run directory.")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    describe = subparsers.add_parser("describe", help="summarise a run directory")
    describe.add_argument("rundir", type=Path)

    plot = subparsers.add_parser("plot", help="four-panel snapshot of one timestep")
    plot.add_argument("rundir", type=Path)
    plot.add_argument("-t", "--timestep", type=int, default=None, help="default: the last one")
    plot.add_argument("-o", "--output", type=Path, default=None, help="save instead of showing")

    diff = subparsers.add_parser("diff", help="difference between two timesteps")
    diff.add_argument("rundir", type=Path)
    diff.add_argument("first", type=int)
    diff.add_argument("second", type=int)
    diff.add_argument("-o", "--output", type=Path, default=None, help="save instead of showing")

    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    case = CaseOutput(args.rundir)

    if args.command == "describe":
        timesteps = case.timesteps
        print(f"run directory : {case.rundir}")
        print(f"lattice       : {case.shape[1]} x {case.shape[0]} (Lx x Ly)")
        print(f"timesteps     : {len(timesteps)}")
        if timesteps:
            print(f"range         : {timesteps[0]} .. {timesteps[-1]}")
            phase = case.phase(timesteps[-1])
            print(f"phase range   : [{phase.min():.3f}, {phase.max():.3f}]")
        return 0

    # plotting commands need matplotlib, which stays an optional import
    import matplotlib.pyplot as plt

    from pycglbm import plotter

    if args.command == "plot":
        figure, _ = plotter.plot_fields(case, _resolve(case, args.timestep))
    else:
        figure, _ = plotter.plot_difference(case, args.first, args.second)

    if args.output:
        figure.savefig(args.output, dpi=120, bbox_inches="tight")
        print(f"wrote {args.output}")
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
