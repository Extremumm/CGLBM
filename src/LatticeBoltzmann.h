#ifndef LATTICE_BOLTZMANN_H
#define LATTICE_BOLTZMANN_H

#include "Constants.h"

// Define your lattice data structures
// Discrete velocities:
// $$
// \vec{\xi}_i=\frac{\Delta x}{\Delta t} \begin{cases}(0,0) & i=0 \\ ( \pm 1,0) \text { or }(0, \pm 1) & i=1,2,3,4 \\ ( \pm 1, \pm 1) \text { or }( \pm 1, \mp 1) & i=5,6,7,8\end{cases}
// $$
// Speed of sound:
// $$
// c_s=\frac{\Delta x}{\sqrt{3} \Delta t} .
// $$
// Weight $w_i= \begin{cases}\frac{4}{9} , & i=0 \\ \frac{1}{9} , & i=1,2,3,4 \\ \frac{1}{36}, & i=5,6,7,8\end{cases}$
extern double rho[Lx][Ly]; // Density
extern double rho_mdt[Lx][Ly]; // Density at the previous time step
extern double u[Lx][Ly][2]; // Velocity components
extern double p[Lx][Ly]; // Pressure
extern double p_mdt[Lx][Ly]; // Pressure at the previous time step
extern double phi[Lx][Ly]; // Phase field function


extern double omega_1[Lx][Ly][Q]; // Collision term
extern double omega_2[Lx][Ly][Q]; // Collision term related to surface tension
extern double omega_3[Lx][Ly][Q]; // Recoloring term
extern double f[Lx][Ly][Q]; // Sum of distribution functions
extern double g[Lx][Ly][Q]; // Difference of distribution functions
extern double f_eq[Lx][Ly][Q]; // Sum of distribution functions at equilibrium
//double g_eq[Lx][Ly][Q]; // Difference of distribution functions at equilibrium

void runSimulation();
void collide();
void force();
void stream();
void recolor();

#endif
