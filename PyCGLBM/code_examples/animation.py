"""Animate a run, with a box to jump to a given timestep.

Run it from inside a run directory.
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, TextBox

from pycglbm import CaseOutput
from pycglbm.plotter import plot_fields

case = CaseOutput(".")
timesteps = case.timesteps

figure, axes = plot_fields(case, timesteps[0])


def update_plot(frame):
    """Redraw the four panels in place for the timestep at index ``frame``."""
    timestep = timesteps[frame % len(timesteps)]
    fields = case.fields(timestep)
    for axis, name in zip((axes[0], axes[2], axes[3]), ("density", "phase", "pressure")):
        axis.images[0].set_data(fields[name])
    figure.suptitle(f"{case.rundir.name} - timestep {timestep}")


animation = FuncAnimation(figure, update_plot, frames=len(timesteps), interval=200)

axtbox = plt.axes([0.1, 0.01, 0.05, 0.05])
text_box = TextBox(axtbox, "Timestep", initial=str(timesteps[0]))
axbutton = plt.axes([0.18, 0.01, 0.03, 0.05])
button = Button(axbutton, "Plot")


def plot_specific_timestep(event):
    plot_fields(case, int(text_box.text))
    plt.show()


button.on_clicked(plot_specific_timestep)
plt.show()
