from math import sqrt, tanh, cos, pi
import numpy as np
Lx = Ly = 128
x0 = Lx / 2
y0 = Ly / 2
radius = 10
dx = 1
r =  radius * dx 
ch_width_ope = 1.0  # Assuming a value for ch_width_ope
rho1 = 1
rho2 = 1
c1 = c2 = 1/ sqrt(3)
c_dx = 1.e-5
c_dt = c_dx/347./sqrt(3.)
sigma = 0/(c_dx * c_dx * c_dx / c_dt / c_dt)
p1_inf = rho1 * c1 * c1 - rho2 * c2 * c2; 
p2_inf = 0.
p= np.zeros([Lx, Ly])
for i in range(Lx):
    for j in range(Ly):
        distance = sqrt((i - x0) * (i - x0) + (j - y0) * (j - y0))
        phi_local = -tanh((distance - r) / ch_width_ope)
        rho_local = rho1 * (0.5 + 0.5*phi_local) + rho2 * (0.5 - 0.5*phi_local)
        p_local = (1/2) * (rho_local * c1 * c1 - p1_inf - p2_inf + sqrt((p2_inf - p1_inf + rho_local * c1 * c1 * c2 * c2) ** 2 + rho_local ** 2 * (1 - phi_local ** 2) * c1 * c1 * c2 * c2))
        p[i][j] = p_local
import matplotlib.pyplot as plt

plt.imshow(p, cmap='jet')
plt.title('Pression')
plt.colorbar()
plt.show()
