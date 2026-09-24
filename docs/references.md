# References

## Method

1. **T. Lafarge, P. Boivin, N. Odier, B. Cuenot.**
   *Improved color-gradient method for lattice Boltzmann modeling of two-phase
   flows.* Physics of Fluids **33** (8), 082110, 2021.
   DOI: [10.1063/5.0061638](https://doi.org/10.1063/5.0061638) · HAL: `hal-03324224`

   The algorithm implemented here: the three-part collision operator, the
   two-component equation of state, and the recoloring step.

## High density ratio in colour-gradient models

3. **S. Leclaire, M. Reggio, J.-Y. Trépanier.** *Isotropic color gradient for
   simulating very high-density ratios with a two-phase flow lattice Boltzmann
   model.* Computers & Fluids **48**(1), 98–112, 2011.
   DOI: [10.1016/j.compfluid.2011.04.001](https://doi.org/10.1016/j.compfluid.2011.04.001)

   The isotropic gradient stencils of `src/lbm/isotropic_gradient.h`, and the
   colour gradient taken on the normalised density rather than on the phase
   field — `normalised_phase` in `src/lbm/mixture.h`.

4. **M. Sbragaglia, R. Benzi, L. Biferale, S. Succi, K. Sugiyama, F. Toschi.**
   *Generalized lattice Boltzmann method with multirange pseudopotential.*
   Physical Review E **75**, 026702, 2007.

   Construction of the isotropic shell weights; source of the E4 and E6 sets.
   The E8 weights are tabulated in
   [arXiv:2505.23647](https://arxiv.org/abs/2505.23647), appendix B.

5. **Y. Ba, H. Liu, Q. Li, Q. Kang, J. Sun.** *Multiple-relaxation-time
   color-gradient lattice Boltzmann model for simulating two-phase flows with
   high density ratio.* Physical Review E **94**, 023310, 2016.
   DOI: [10.1103/PhysRevE.94.023310](https://doi.org/10.1103/PhysRevE.94.023310)

   MRT collision plus a third-order Hermite equilibrium and a Chapman–Enskog
   source term; validated to density ratio 1000. Its colour gradient on the
   normalised density is used by `laplace`; the collision and equilibrium are
   not implemented here.

6. **S. Leclaire, M. Reggio, J.-Y. Trépanier.** *Enhanced equilibrium
   distribution functions for simulating immiscible multiphase flows with
   variable density ratios.* International Journal of Multiphase Flow, 2013.
   DOI: [10.1016/j.ijmultiphaseflow.2013.05.009](https://doi.org/10.1016/j.ijmultiphaseflow.2013.05.009)

7. **A. Subhedar.** *Color-gradient lattice Boltzmann model for immiscible
   fluids with density contrast.* Physical Review E **106**, 045308, 2022.
   DOI: [10.1103/PhysRevE.106.045308](https://doi.org/10.1103/PhysRevE.106.045308)

   Velocity-based equilibrium; interface mobility independent of the density
   ratio. Not implemented here — a candidate for the moving-interface limit
   recorded in [`numerics.md`](numerics.md#what-is-not-solved).

8. **S. Saito, N. Takada, S. Baba, S. Someya, H. Ito.** *Generalized equilibria
   for color-gradient lattice Boltzmann model based on higher-order Hermite
   polynomials.* Physical Review E **108**, 065305, 2023.
   [arXiv:2309.07801](https://arxiv.org/abs/2309.07801)

   Also the clearest recent survey of this model family.

9. **M. Latva-Kokko, D. H. Rothman.** *Diffusion properties of gradient-based
   lattice Boltzmann models of immiscible fluids.* Physical Review E **71**,
   056702, 2005.
   DOI: [10.1103/PhysRevE.71.056702](https://doi.org/10.1103/PhysRevE.71.056702)

   The recolouring lineage that `recolor()` sits in.

## Surface tension as a force

10. **B. Lafaurie, C. Nardone, R. Scardovelli, S. Zaleski, G. Zanetti.**
    *Modelling merging and fragmentation in multiphase flows with SURFER.*
    Journal of Computational Physics **113**, 134–147, 1994.

    The capillary stress tensor
    $\tfrac{\sigma}{2}(\lvert\nabla c\rvert\mathsf I - \nabla c\nabla c/\lvert\nabla c\rvert)$,
    whose divergence is the surface force. `src/lbm/surface_force.h`, used by
    `laplace` in place of the perturbation operator Ω⁽²⁾.

11. **J. U. Brackbill, D. B. Kothe, C. Zemach.** *A continuum method for
    modeling surface tension.* Journal of Computational Physics **100**,
    335–354, 1992. And **S. V. Lishchuk, C. M. Care, I. Halliday.** *Lattice
    Boltzmann algorithm for surface tension with greatly reduced
    microcurrents.* Physical Review E **67**, 036701, 2003.

    The same force in curvature form, $\sigma\kappa\nabla c$, the second on the
    colour field of a colour-gradient model. Not used: it does not conserve
    momentum on the lattice (see [`numerics.md`](numerics.md#3-surface-tension-as-a-body-force)).

12. **Z. Guo, C. Zheng, B. Shi.** *Discrete lattice effects on the forcing term
    in the lattice Boltzmann method.* Physical Review E **65**, 046308, 2002.

    The forcing scheme that `S_F` and the half-force velocity correction
    implement, and through which the surface force enters.

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
