"""Figures built from a :class:`pycglbm.files.CaseOutput`."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from pycglbm.files import CaseOutput

#: How many nodes to skip between two velocity arrows.
QUIVER_STRIDE = 4


def _quiver(axis, velocity: np.ndarray, stride: int = QUIVER_STRIDE) -> None:
    ny, nx = velocity.shape[:2]
    y, x = np.mgrid[0:ny, 0:nx]
    axis.quiver(
        x[::stride, ::stride],
        y[::stride, ::stride],
        velocity[::stride, ::stride, 0],
        velocity[::stride, ::stride, 1],
        color="black",
        angles="xy",
    )
    axis.set_aspect("equal")


def plot_fields(case: CaseOutput, timestep: int, stride: int = QUIVER_STRIDE):
    """Four-panel snapshot: density, velocity, phase field and pressure."""
    fields = case.fields(timestep)
    figure, axes = plt.subplots(1, 4, figsize=(20, 5))

    for axis, (name, cmap) in zip(
        (axes[0], axes[2], axes[3]),
        (("density", "viridis"), ("phase", "coolwarm"), ("pressure", "viridis")),
    ):
        image = axis.imshow(fields[name], origin="lower", cmap=cmap)
        figure.colorbar(image, ax=axis)
        axis.set_title(name.capitalize())

    _quiver(axes[1], fields["velocity"], stride)
    axes[1].set_title("Velocity Field")

    figure.suptitle(f"{case.rundir.name} - timestep {timestep}")
    return figure, axes


def plot_difference(case: CaseOutput, first: int, second: int, stride: int = QUIVER_STRIDE):
    """The same four panels, showing ``second - first``."""
    before, after = case.fields(first), case.fields(second)
    figure, axes = plt.subplots(1, 4, figsize=(20, 5))

    for axis, (name, cmap) in zip(
        (axes[0], axes[2], axes[3]),
        (("density", "viridis"), ("phase", "coolwarm"), ("pressure", "viridis")),
    ):
        image = axis.imshow(after[name] - before[name], origin="lower", cmap=cmap)
        figure.colorbar(image, ax=axis)
        axis.set_title(f"{name.capitalize()} Difference")

    _quiver(axes[1], after["velocity"] - before["velocity"], stride)
    axes[1].set_title("Velocity Field Difference")

    figure.suptitle(f"{case.rundir.name} - timestep {second} minus {first}")
    return figure, axes
