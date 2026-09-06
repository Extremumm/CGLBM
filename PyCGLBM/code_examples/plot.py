"""Four-panel snapshot of one timestep. Run it from inside a run directory."""

import matplotlib.pyplot as plt
from pycglbm import CaseOutput
from pycglbm.plotter import plot_fields

case = CaseOutput(".")
timestep = int(input(f"Enter the time step {case.timesteps[:1]}..{case.timesteps[-1:]}: "))

plot_fields(case, timestep)

# pressure along the horizontal centreline, used to read off the Laplace jump
pressure = case.pressure(timestep)
row = pressure.shape[0] // 2
print("p at x = 64, 74, 84, 94:", [round(pressure[row, x], 6) for x in (64, 74, 84, 94)])

plt.show()
