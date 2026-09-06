# References

## Method

1. **T. Lafarge, P. Boivin, N. Odier, B. Cuenot.**
   *Improved color-gradient method for lattice Boltzmann modeling of two-phase
   flows.* Physics of Fluids **33** (8), 082110, 2021.
   DOI: [10.1063/5.0061638](https://doi.org/10.1063/5.0061638) · HAL: `hal-03324224`

   The algorithm implemented here: the three-part collision operator, the
   two-component equation of state, and the recoloring step.

## Lattice Boltzmann background

2. **T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E. M. Viggen.**
   *The Lattice Boltzmann Method: Principles and Practice.* Graduate Texts in
   Physics, Springer, 2017.

   Referenced throughout the sources: periodic boundary conditions (p. 170),
   conversion between lattice and physical units for the viscosity (p. 284).

## Unit conversion

The programs work in lattice units and convert with two factors defined at the
top of each `main_*.cpp`:

- `c_dx = 1e-5` m — length of one lattice spacing,
- `c_dt = c_dx / (347 · √3)` s — duration of one time step, fixed by matching
  the lattice sound speed $c_s = \Delta x/(\sqrt3\,\Delta t)$ to 347 m/s in air.

A physical quantity of dimension $L^a T^b$ is divided by `c_dx^a · c_dt^b` to
obtain its lattice-unit value — the pattern visible in the definitions of `nu`,
`sigma` and `a_g`. `programs/initial_conditions/laplace/init_laplace.py`
evaluates the same equation of state in Python, and is the reference for the
initial pressure field of the `laplace` program.
