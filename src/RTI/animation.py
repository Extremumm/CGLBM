import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox

def plot_data(timestep):
    density = np.loadtxt(f'density_{timestep}.csv', delimiter=',')
    velocity = np.loadtxt(f'velocity_{timestep}.csv', delimiter=',').reshape(density.shape[0], density.shape[1], 2)
    phase = np.loadtxt(f'phase_{timestep}.csv', delimiter=',')
    pressure = np.loadtxt(f'pressure_{timestep}.csv', delimiter=',')

    fig, axes = plt.subplots(1, 4, figsize=(20, 5))

    im1 = axes[0].imshow(density, origin='lower', cmap='viridis')
    fig.colorbar(im1, ax=axes[0])
    axes[0].set_title('Density')

    nx = density.shape[1]
    ny = density.shape[0]
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    Y, X = np.meshgrid(y, x)

    axes[1].quiver(X[::5, ::5], Y[::5, ::5], velocity[::5, ::5, 0], velocity[::5, ::5, 1],
                   color='black', angles='xy', scale_units='xy', scale=1)
    axes[1].set_title('Velocity Field')

    im3 = axes[2].imshow(phase, origin='lower', cmap='coolwarm')
    fig.colorbar(im3, ax=axes[2])
    axes[2].set_title('Phase Field')

    im4 = axes[3].imshow(pressure, origin='lower', cmap='viridis')
    fig.colorbar(im4, ax=axes[3])
    axes[3].set_title('Pressure')
    fig.suptitle(f'Time Step: {timestep}', fontsize=16)
    plt.show()

def animate_plot():
    # last_timestep = int(input("Enter the last time step: "))
    last_timestep = 5000000
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    pause = False

    def update_plot(timestep):
        if pause:
            return

        density = np.loadtxt(f'density_{timestep}.csv', delimiter=',')
        velocity = np.loadtxt(f'velocity_{timestep}.csv', delimiter=',').reshape(density.shape[0], density.shape[1], 2)
        phase = np.loadtxt(f'phase_{timestep}.csv', delimiter=',')
        pressure = np.loadtxt(f'pressure_{timestep}.csv', delimiter=',')

        im1.set_data(density)
        im1.set_clim(vmin=density.min(), vmax=density.max())
        im1.colorbar.update_normal(im1)

        im2.set_UVC(velocity[::5, ::5, 0], velocity[::5, ::5, 1])

        im3.set_data(phase)
        im3.set_clim(vmin=phase.min(), vmax=phase.max())
        im3.colorbar.update_normal(im3)

        im4.set_data(pressure)
        im4.set_clim(vmin=pressure.min(), vmax=pressure.max())
        im4.colorbar.update_normal(im4)
        fig.suptitle(f'Time Step: {timestep}', fontsize=16)

    density = np.loadtxt(f'density_0.csv', delimiter=',')
    velocity = np.loadtxt(f'velocity_0.csv', delimiter=',').reshape(density.shape[0], density.shape[1], 2)
    phase = np.loadtxt(f'phase_0.csv', delimiter=',')
    pressure = np.loadtxt(f'pressure_0.csv', delimiter=',')
    
    im1 = axes[0].imshow(density, origin='lower', cmap='viridis')
    fig.colorbar(im1, ax=axes[0])
    axes[0].set_title('Density')

    nx = density.shape[1]
    ny = density.shape[0]
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    Y, X = np.meshgrid(y, x)

    im2 = axes[1].quiver(X[::5, ::5], Y[::5, ::5], velocity[::5, ::5, 0], velocity[::5, ::5, 1],
                         color='black', angles='xy', scale_units='xy', scale=1)
    axes[1].set_title('Velocity Field')

    im3 = axes[2].imshow(phase, origin='lower', cmap='coolwarm')
    fig.colorbar(im3, ax=axes[2])
    axes[2].set_title('Phase Field')

    im4 = axes[3].imshow(pressure, origin='lower', cmap='viridis')
    fig.colorbar(im4, ax=axes[3])
    axes[3].set_title('Pressure')

    axpause = plt.axes([0.22, 0.01, 0.03, 0.075])
    bpause = Button(axpause, 'Pause')

    def pause_animation(event):
        nonlocal pause
        pause = not pause
        if pause:
            bpause.label.set_text('Resume')
        else:
            bpause.label.set_text('Pause')

    bpause.on_clicked(pause_animation)

    axtbox = plt.axes([0.1, 0.01, 0.05, 0.075])
    text_box = TextBox(axtbox, 'Timestep', initial='0')

    axbutton = plt.axes([0.18, 0.01, 0.03, 0.075])
    bplot = Button(axbutton, 'Plot')

    def plot_specific_timestep(event):
        timestep = int(text_box.text)
        plot_data(timestep)

    bplot.on_clicked(plot_specific_timestep)

    animation = FuncAnimation(fig, update_plot, frames=range(0, last_timestep, 10000), interval=200)
    # animation.save('cap.mp4', writer='ffmpeg')
    plt.show()

animate_plot()