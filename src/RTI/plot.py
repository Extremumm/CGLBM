import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

def plot_data(timestep):
    # timestep = int(input("Enter the last time step: "))
    # Load the data
    density = np.loadtxt(f'density_{timestep}.csv', delimiter=',')
    velocity = np.loadtxt(f'velocity_{timestep}.csv', delimiter=',').reshape(density.shape[0], density.shape[1], 2) 
    phase = np.loadtxt(f'phase_{timestep}.csv', delimiter=',')
    pressure = np.loadtxt(f'pressure_{timestep}.csv', delimiter=',')
    # Create a figure with subplots
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))

    # Density plot
    im1 = axes[0].imshow(density, origin='lower', cmap='viridis')
    fig.colorbar(im1, ax=axes[0])
    axes[0].set_title('Density')
    print(density.shape)
    # Velocity field plot (with quivers)
    nx = density.shape[1]
    ny = density.shape[0]
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    Y, X = np.meshgrid(y, x)

    # Adjust the density of arrows for better visualization (e.g., every 5th point)
    axes[1].quiver(X[::4, ::4], Y[::4, ::4], velocity[::4, ::4, 0], velocity[::4, ::4, 1], 
                  color='black', angles='xy', scale_units='xy', scale=1)
    axes[1].set_title('Velocity Field')

    # Phase field plot
    im3 = axes[2].imshow(phase, origin='lower', cmap='coolwarm')
    fig.colorbar(im3, ax=axes[2])
    axes[2].set_title('Phase Field')

    # Pressure plot
    im4 = axes[3].imshow(pressure, origin='lower', cmap='viridis')
    fig.colorbar(im4, ax=axes[3])
    axes[3].set_title('Pressure')
    fig.suptitle(f'Time Step: {timestep}', fontsize=16)
    plt.show()

def animate_plot():
    # Get the last time step entered by the user
    last_timestep = int(input("Enter the last time step: "))

    # Create a figure with subplots
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))

    def update_plot(timestep):
        # Load the data for the current timestep
        density = np.loadtxt(f'density_{timestep}.csv', delimiter=',')
        velocity = np.loadtxt(f'velocity_{timestep}.csv', delimiter=',').reshape(density.shape[0], density.shape[1], 2) 
        phase = np.loadtxt(f'phase_{timestep}.csv', delimiter=',')
        pressure = np.loadtxt(f'pressure_{timestep}.csv', delimiter=',')

        # Update the plots for the current timestep
        im1.set_data(density)
        im2.set_UVC(velocity[::5, ::5, 0], velocity[::5, ::5, 1])
        im3.set_data(phase)
        im4.set_data(pressure)

    # Load the data for the initial timestep
    density = np.loadtxt('density_0.csv', delimiter=',')
    velocity = np.loadtxt('velocity_0.csv', delimiter=',').reshape(density.shape[0], density.shape[1], 2) 
    phase = np.loadtxt('phase_0.csv', delimiter=',')
    pressure = np.loadtxt('pressure_0.csv', delimiter=',')

    # Density plot
    im1 = axes[0].imshow(density, origin='lower', cmap='viridis')
    fig.colorbar(im1, ax=axes[0])
    axes[0].set_title('Density')

    # Velocity field plot (with quivers)
    nx = density.shape[1]
    ny = density.shape[0]
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    Y, X = np.meshgrid(y, x)

    # Adjust the density of arrows for better visualization (e.g., every 5th point)
    im2 = axes[1].quiver(X[::5, ::5], Y[::5, ::5], velocity[::5, ::5, 0], velocity[::5, ::5, 1], 
                  color='black', angles='xy', scale_units='xy', scale=1)
    axes[1].set_title('Velocity Field')

    # Phase field plot
    im3 = axes[2].imshow(phase, origin='lower', cmap='coolwarm')
    fig.colorbar(im3, ax=axes[2])
    axes[2].set_title('Phase Field')

    # Pressure plot
    im4 = axes[3].imshow(pressure, origin='lower', cmap='viridis')
    fig.colorbar(im4, ax=axes[3])
    axes[3].set_title('Pressure')

    # Create the animation
    animation = FuncAnimation(fig, update_plot, frames=range(0, last_timestep), interval=100)
    plt.show()

# Example usage
t = input("Enter the time step: ")
plot_data(t)
p = np.loadtxt(f'pressure_{t}.csv', delimiter=',')
p0 = p[64, 64]
p1 = p[74, 64]
p2 = p[84, 64]
p3 = p[94, 64]
print(p0, p1, p2, p3)
#animate_plot()