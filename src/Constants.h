#ifndef CONSTANTS_H
#define CONSTANTS_H

const int Lx = 256; // Number of lattice nodes in the x-direction
const int Ly = Lx; // Number of lattice nodes in the y-direction
const int Q = 9; // Number of discrete velocities

// Define your simulation parameters without dimensions
const double dx = 1.; // Lattice spacing
const double dt = 1.; // Time step
const int numSteps = 1000; // Number of simulation steps

// Define boundary conditions
//The pressure at infinity is used at the equation of state to calculate the pressure and equilibrium distribution function
const double p1_inf = 0.; // Pressure at infinity for component 1
const double p2_inf = 0.; // Pressure at infinity for component 2

// viscosities to calculate Relaxation times in collide step
const double nu    = 0.1; // shear viscosity
const double nu_b  = 0.1; // bulk viscosity

const double sigma = 0.1; // surface tension

const double c1    = 0.1; // Speed of sound for component 1
const double c2    = 0.1; // Speed of sound for component 2

const double ch_width = 1.6*dx; // Characteristic width of the interface

//discrete velocities set
const double xi[Q][2] = {
    {0, 0},
    {1, 0}, {0, 1}, {-1, 0}, {0, -1},
    {1, 1}, {-1, 1}, {-1, -1}, {1, -1}
};

//weights
const double w1 = 4./9.;
const double w2 = 1./9.;
const double w3 = 1./36.;
const double w[Q] = {w1, w2, w2, w2, w2, w3, w3, w3, w3};

// Speed of sound and related constants
const double cs = dx * 0.5773502691896258 / dt; // Speed of sound in the lattice with 0.5773502691896258 = 1/sqrt(3.0)
const double cs2 = cs * cs;
const double cs4 = cs2 * cs2;
const double cs6 = cs4 * cs2;
const double c1_squared = c1 * c1;
const double c2_squared = c2 * c2;

#endif
