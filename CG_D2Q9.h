#ifndef CG_D2Q9_H
#define CG_D2Q9_H

#include <cmath>

// Define your lattice dimensions and parameters
const int Lx = 100; // Number of lattice nodes in the x-direction
const int Ly = 100; // Number of lattice nodes in the y-direction
const int Q = 9; // Number of discrete velocities

// Define your simulation parameters without dimensions
const double dx; // Lattice spacing
const double dt; // Time step
const int numSteps; // Number of simulation steps

// Define boundary conditions
//The pressure at infinity is used at the equation of state to calculate the pressure and equilibrium distribution function
double p1_inf; // Pressure at infinity for component 1
double p2_inf; // Pressure at infinity for component 2

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
double f[Lx][Ly][Q]; // Distribution functions
double rho[Lx][Ly]; // Density
double rho_mdt[Lx][Ly]; // Density at the previous time step
double u[Lx][Ly][2]; // Velocity components
double p[Lx][Ly]; // Pressure
double p_mdt[Lx][Ly]; // Pressure at the previous time step

double omega_1[Lx][Ly][Q]; // Collision term
double omega_2[Lx][Ly][Q]; // Collision term related to surface tension
double omega_3[Lx][Ly][Q]; // Recoloring term
//discrete velocities set
const double xi[Q][2] = {
    {0, 0},
    {1, 0}, {0, 1}, {-1, 0}, {0, -1},
    {1, 1}, {-1, 1}, {-1, -1}, {1, -1}
};

// viscosities to calculate Relaxation times in collide step
const double nu; // shear viscosity
const double nu_b; //  bulk viscosity

const double sigma; //surface tension

double S[Lx][Ly][Q]; // Force term
double F[Lx][Ly][2]; // External force

//weights
const double w1 = 4./9.;
const double w2 = 1./9.;
const double w3 = 1./36.;
const double w[Q] = {w1, w2, w2, w2, w2, w3, w3, w3, w3};

// Speed of sound and related constants
const double cs = 1.0 * dx / sqrt(3.0) / dt; // Speed of sound in the lattice
const double cs2 = cs * cs;
const double cs4 = cs2 * cs2;
const double cs6 = cs4 * cs2;
const double c1 ; // Speed of sound for component 1
const double c2 ; // Speed of sound for component 2
const double c1_squared = c1 * c1;
const double c2_squared = c2 * c2;

double f[Lx][Ly][Q]; // Sum of distribution functions
double g[Lx][Ly][Q]; // Difference of distribution functions
double f_eq[Lx][Ly][Q]; // Sum of distribution functions at equilibrium

//double g_eq[Lx][Ly][Q]; // Difference of distribution functions at equilibrium

double phi[Lx][Ly]; // Phase field function

// Function declarations
void initializeFields();
void collideAndStream();
void calculateEquilibrium();
void collide();

#endif // CG_D2Q9_H
