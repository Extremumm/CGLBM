import numpy as np
import matplotlib.pyplot as plt


def plot_differences(time_step_1, time_step_2):
    density_1 = np.loadtxt(f'density_{time_step_1}.csv', delimiter=',')
    density_2 = np.loadtxt(f'density_{time_step_2}.csv', delimiter=',')
    velocity_1 = np.loadtxt(f'velocity_{time_step_1}.csv', delimiter=',').reshape(density_1.shape[0], density_1.shape[1], 2)
    velocity_2 = np.loadtxt(f'velocity_{time_step_2}.csv', delimiter=',').reshape(density_2.shape[0], density_2.shape[1], 2)
    phase_1 = np.loadtxt(f'phase_{time_step_1}.csv', delimiter=',')
    phase_2 = np.loadtxt(f'phase_{time_step_2}.csv', delimiter=',')
    pressure_1 = np.loadtxt(f'pressure_{time_step_1}.csv', delimiter=',')
    pressure_2 = np.loadtxt(f'pressure_{time_step_2}.csv', delimiter=',')

    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    im1 = axes[0].imshow(density_2 - density_1, origin='lower', cmap='viridis')
    fig.colorbar(im1, ax=axes[0])
    axes[0].set_title('Density Difference')

    nx = density_1.shape[1]
    ny = density_1.shape[0]
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    Y, X = np.meshgrid(y, x)

    axes[1].quiver(X[::5, ::5], Y[::5, ::5], velocity_2[::5, ::5, 0] - velocity_1[::5, ::5, 0], velocity_2[::5, ::5, 1] - velocity_1[::5, ::5, 1],
                    color='black', angles='xy', scale_units='xy', scale=1)
    axes[1].set_title('Velocity Field Difference')
    im3 = axes[2].imshow(phase_2 - phase_1, origin='lower', cmap='coolwarm')
    fig.colorbar(im3, ax=axes[2])
    axes[2].set_title('Phase Field Difference')
    im4 = axes[3].imshow(pressure_2 - pressure_1, origin='lower', cmap='viridis')
    fig.colorbar(im4, ax=axes[3])
    axes[3].set_title('Pressure Difference')
    plt.show()

start = int(input("Enter the first time step: "))
end = int(input("Enter the second time step: "))
plot_differences(start, end)