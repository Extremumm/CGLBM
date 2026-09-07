"""Regenerate the figures of docs/report/report.tex from the data under data/.

    python3 docs/report/figures/make_figures.py

The CSVs are produced by docs/report/tools/report_data.cpp; the last figure is
analytic and computed here. Everything is written as PDF beside this script.
"""

from __future__ import annotations

import csv
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"

#: Written into every PDF instead of the current time, so that regenerating the
#: figures from unchanged data produces byte-identical files and committing them
#: does not churn.
SAVE = {"metadata": {"CreationDate": None}}

#: One colour and label per model, used by every figure.
MODELS = {
    "eos_legacy": ("#B3462F", r"EOS solver, $\Omega^{(2)}$ stress"),
    "eos_csf": ("#C89632", r"EOS solver, $\varphi_N$ + CSF"),
    "two_population": ("#2B6CB0", "Two-population solver"),
}

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "axes.axisbelow": True,
        "figure.constrained_layout.use": True,
        "legend.frameon": False,
    }
)


def read(name):
    """Read one CSV into a dict of columns, numbers where possible."""
    with open(DATA / name, newline="") as handle:
        rows = list(csv.DictReader(handle))
    columns = {key: [] for key in rows[0]}
    for row in rows:
        for key, value in row.items():
            try:
                columns[key].append(float(value))
            except ValueError:
                columns[key].append(value)
    return {key: np.array(value) for key, value in columns.items()}


def figure_sweep():
    """Laplace's law and spurious currents against the density ratio."""
    data = read("laplace_sweep.csv")
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(5.4, 5.0), sharex=True)
    diverges_at = np.inf

    for model, (colour, label) in MODELS.items():
        mask = data["model"] == model
        ratio, value = data["ratio"][mask], data["sigma_ratio"][mask]
        speed, alive = data["peak_speed"][mask], data["survived"][mask]
        ok = alive == 1
        top.plot(ratio[ok], value[ok], "o-", color=colour, label=label, markersize=4, lw=1.3)
        bottom.plot(ratio[ok], speed[ok], "o-", color=colour, markersize=4, lw=1.3)
        if (~ok).any():
            diverges_at = min(diverges_at, ratio[~ok].min())

    # Both configurations of the equation-of-state solver stop at the same
    # density ratio, so shade the region rather than overplot two markers.
    for axis in (top, bottom):
        axis.axvspan(diverges_at * 0.72, data["ratio"].max() * 1.6, color="0.88", zorder=0)
    top.annotate(
        "equation-of-state\nsolver diverges",
        xy=(diverges_at * 0.95, 0.30),
        ha="left",
        va="center",
        fontsize=7.5,
        color="0.35",
    )
    top.axhline(1.0, color="0.35", lw=0.8, ls="--", zorder=1)
    top.set_ylabel(r"$\sigma_{\mathrm{measured}} \, / \, \sigma$")
    top.set_ylim(0.0, 1.30)
    top.set_xlim(1.7, 1.6e5)
    top.legend(loc="lower left", fontsize=8)
    top.set_title("Laplace's law against density ratio", fontsize=9.5, loc="left")

    bottom.set_yscale("log")
    bottom.set_xscale("log")
    bottom.set_xlabel(r"density ratio $\rho_1/\rho_2$")
    bottom.set_ylabel(r"spurious currents $\max|\mathbf{u}|$")
    bottom.set_title("Spurious currents", fontsize=9.5, loc="left")
    fig.savefig(HERE / "fig_sweep.pdf", **SAVE)
    plt.close(fig)


def figure_history():
    """The measured tension against time, at a density ratio of 1000."""
    data = read("history_r1000.csv")
    fig, axis = plt.subplots(figsize=(5.4, 3.0))
    for model, (colour, label) in MODELS.items():
        mask = data["model"] == model
        step, value = data["step"][mask], data["sigma_ratio"][mask]
        if step.size == 0:
            continue
        low, high = -0.35, 1.6
        axis.plot(step / 1e3, np.clip(value, low, high), "-", color=colour, label=label, lw=1.4)
        # mark where the run stopped, if it stopped before the end
        if step.max() < 118000:
            end = float(np.clip(value[-1], low + 0.05, high - 0.05))
            axis.plot([step.max() / 1e3], [end], "x", color=colour, markersize=8, mew=1.8)
            axis.annotate(
                f"gone by {step.max() / 1e3:.0f}k",
                xy=(step.max() / 1e3, end),
                xytext=(8, 4),
                textcoords="offset points",
                ha="left",
                fontsize=7,
                color=colour,
            )
    axis.axhline(1.0, color="0.35", lw=0.8, ls="--", zorder=0)
    axis.set_xlabel(r"time step / $10^3$")
    axis.set_ylabel(r"$\sigma_{\mathrm{measured}} \, / \, \sigma$")
    axis.set_ylim(-0.35, 1.6)
    axis.legend(loc="lower left", fontsize=8)
    axis.set_title(
        r"Density ratio $10^3$: what a short run would have reported", fontsize=9.5, loc="left"
    )
    fig.savefig(HERE / "fig_history.pdf", **SAVE)
    plt.close(fig)


def figure_profiles():
    """Phase, density and pressure across the interface at a ratio of 1000."""
    data = read("profiles_r1000.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))
    for model in ("eos_csf", "two_population"):
        colour, label = MODELS[model]
        mask = data["model"] == model
        radius = data["r"][mask]
        # Centre on the *physical* interface, which is the density midpoint in
        # both models (Proposition 3). It is not the phase field's zero contour
        # for the equation-of-state model, and using that would misalign the
        # two profiles by several nodes.
        density = data["rho"][mask]
        midpoint = 0.5 * (density[0] + density[-1])
        crossing = np.interp(-midpoint, -density, radius)
        offset = radius - crossing
        axes[0].plot(offset, data["phase"][mask], "-", color=colour, label=label, lw=1.4)
        axes[1].plot(offset, density, "-", color=colour, lw=1.4)
        axes[2].plot(offset, data["p"][mask] / data["p"][mask][-1], "-", color=colour, lw=1.4)

    for axis in axes:
        axis.set_xlim(-8, 8)
        axis.set_xlabel("distance from interface")
    axes[0].set_ylabel("phase field")
    axes[1].set_ylabel(r"density $\rho$")
    axes[1].set_yscale("log")
    axes[2].set_ylabel(r"pressure $p \, / \, p_\infty$")
    axes[0].legend(loc="lower left", fontsize=6.5)
    axes[2].set_title("pressure", fontsize=8, loc="left")
    fig.savefig(HERE / "fig_profiles.pdf", **SAVE)
    plt.close(fig)


def figure_currents():
    """Spurious current field around the relaxed droplet."""
    data = read("currents_r1000.csv")
    size = int(np.sqrt(data["i"].size))
    grid = lambda key: data[key].reshape(size, size).T  # noqa: E731
    speed = np.hypot(grid("ux"), grid("uy"))

    fig, axis = plt.subplots(figsize=(4.0, 3.4))
    image = axis.imshow(speed, origin="lower", cmap="magma")
    axis.contour(grid("phase"), levels=[0.0], colors="w", linewidths=1.0)
    step = 5
    coords = np.arange(0, size, step)
    axis.quiver(
        *np.meshgrid(coords, coords),
        grid("ux")[::step, ::step],
        grid("uy")[::step, ::step],
        color="w",
        alpha=0.65,
        scale=6e-4,
        width=0.003,
    )
    axis.set_xticks([])
    axis.set_yticks([])
    axis.grid(False)
    axis.set_title(
        r"Two-population solver, $\rho_1/\rho_2 = 10^3$"
        "\n"
        rf"peak $|\mathbf{{u}}| = {speed.max():.1e}$",
        fontsize=8.5,
        loc="left",
    )
    fig.colorbar(image, ax=axis, label=r"$|\mathbf{u}|$", shrink=0.85)
    fig.savefig(HERE / "fig_currents.pdf", **SAVE)
    plt.close(fig)


def figure_recolouring():
    """Why the segregation operator has to be rewritten for a density ratio.

    Analytic. At the point where the two fluids occupy equal volume, compare the
    colour the operator wants to move with the colour actually available to
    move: anything above one takes a distribution negative.
    """
    ratios = np.logspace(0, 5, 400)
    alpha2, beta = 0.2, 0.7
    alpha1 = 1.0 - (1.0 - alpha2) / ratios

    rho1, rho2 = ratios / 2.0, 0.5
    rho = rho1 + rho2
    # non-rest populations of an axial direction, and the total
    f1 = rho1 * (1.0 - alpha1) / 5.0
    f2 = rho2 * (1.0 - alpha2) / 5.0
    total = f1 + f2

    available = (rho2 / rho) * total  # fluid 2's share, which must stay positive
    push_w = beta * (1.0 / 9.0) * rho1 * rho2 / rho  # Latva-Kokko, lattice weight
    push_phi = beta * (total / rho) * rho1 * rho2 / rho  # this work, mixture rest weight

    fig, axis = plt.subplots(figsize=(5.4, 3.0))
    axis.plot(ratios, push_w / available, color="#B3462F", lw=1.5, label=r"push $\propto w_i$")
    axis.plot(
        ratios,
        push_phi / available,
        color="#2B6CB0",
        lw=1.5,
        label=r"push $\propto \Phi_i$ (this work)",
    )
    axis.axhline(1.0, color="0.2", lw=1.0, ls="--")
    crossing = ratios[np.argmin(np.abs(push_w / available - 1.0))]
    axis.axvline(crossing, color="#B3462F", lw=0.8, ls=":")
    axis.annotate(
        rf"negative beyond $\rho_1/\rho_2 \approx {crossing:.1f}$",
        xy=(crossing, 30),
        xytext=(6, 0),
        textcoords="offset points",
        fontsize=7.5,
        color="#B3462F",
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel(r"density ratio $\rho_1/\rho_2$")
    axis.set_ylabel("colour moved / colour available")
    axis.set_ylim(1e-2, 1e5)
    axis.legend(loc="upper left", fontsize=8)
    axis.set_title("Why the recolouring needs adapting", fontsize=9.5, loc="left")
    fig.savefig(HERE / "fig_recolouring.pdf", **SAVE)
    plt.close(fig)
    return crossing


def main():
    figure_sweep()
    figure_history()
    figure_profiles()
    figure_currents()
    crossing = figure_recolouring()
    print(f"figures written to {HERE}")
    print(f"  Latva-Kokko positivity is lost above a density ratio of {crossing:.2f}")


if __name__ == "__main__":
    main()
