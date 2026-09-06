"""Difference between two timesteps. Run it from inside a run directory."""

import matplotlib.pyplot as plt
from pycglbm import CaseOutput
from pycglbm.plotter import plot_difference

case = CaseOutput(".")
first = int(input("Enter the first time step: "))
second = int(input("Enter the second time step: "))

plot_difference(case, first, second)
plt.show()
