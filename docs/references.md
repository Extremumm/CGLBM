# References

## Method

1. **T. Lafarge, P. Boivin, N. Odier, B. Cuenot.**
   *Improved color-gradient method for lattice Boltzmann modeling of two-phase
   flows.* Physics of Fluids **33** (8), 082110, 2021.
   DOI: [10.1063/5.0061638](https://doi.org/10.1063/5.0061638) · HAL: `hal-03324224`

   The algorithm implemented here: the three-part collision operator, the
   two-component equation of state, and the recoloring step.

## High density ratio in colour-gradient models

3. **D. Grunau, S. Chen, K. Eggert.** *A lattice Boltzmann model for multiphase
   fluid flows.* Physics of Fluids A **5**, 2557, 1993.
   DOI: [10.1063/1.858769](https://doi.org/10.1063/1.858769)

   The rest-particle weight `alpha_k` that gives each fluid its own sound speed
   and carries the density ratio in the equilibrium. Implemented in
   `TwoPopulationSolver`; it is what lets that model hold a density ratio of
   1000 where the equation-of-state model diverges.

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
   source term; validated to density ratio 1000, at 0.74 % error on σ and
   u_max = 1.3 × 10⁻⁴ for a static droplet of R = 25 in 100².

   The main source for the high-density-ratio work here. On the phase field
   used to locate the interface:

   > Usually, ρᴺ is defined by ρᴺ = (ρᴿ − ρᴮ)/(ρᴿ + ρᴮ). This definition,
   > however, becomes increasingly incorrect in identifying the interface as
   > the density ratio increases.

   What is implemented, and where:

   | Their equation | Here |
   |---|---|
   | (14) enhanced equilibrium | already present as the third-order Hermite term of `Solver::equilibrium()`; the two agree to 10⁻¹⁴ |
   | (17)–(18) source correction | already present as `S_Sp` in `Solver::force()`, same two moments, same isotropic derivative; in `TwoPopulationSolver`, `--third-moment-correction` (`add_third_moment_source`, `src/lbm/mrt.h`) |
   | (20) nine-point derivative | `S_Sp` uses it |
   | (21) normalised phase field | `normalised_phase()`, `InterfaceField::BulkNormalised`, and `CaseConfig::initial_profile_field` |
   | (23)–(26) CSF perturbation | `Solver::surface_force()`, `SurfaceTension::ContinuumSurfaceForce` |
   | (29) velocity redefinition | `Solver::macroscopic()` already forms ρu = Σξf + F dt/2 |
   | (30) Latva-Kokko recolouring | `Recolouring::LatvaKokko`, off by default — measured better below a density ratio of 2 and unusable above 10 |
   | (11)–(13) MRT collision | `--collision=mrt` in `TwoPopulationSolver` (`src/lbm/mrt.h`); `Solver`'s collision is regularised |
   | (22) relaxation interpolation | `Physics::nu2` / `nu_b2`, interpolated on the volume fraction |

   See [`numerics.md`](numerics.md#reading-ba-et-al-2016-and-leclaire-et-al-2013).

6. **S. Leclaire, N. Pellerin, M. Reggio, J.-Y. Trépanier.** *Enhanced
   equilibrium distribution functions for simulating immiscible multiphase
   flows with variable density ratios in a class of lattice Boltzmann models.*
   International Journal of Multiphase Flow **57**, 159–168, 2013.
   DOI: [10.1016/j.ijmultiphaseflow.2013.05.009](https://doi.org/10.1016/j.ijmultiphaseflow.2013.05.009)

   The third-order correction to the colour-gradient equilibrium, which repairs
   the third-order velocity moment when the components' sound speeds differ.
   Reached density ratio 1000 on two-layered Couette flow.

   **Already present here**, though not by that name: the equilibrium in
   `Solver::equilibrium()` is a Hermite expansion and its third-order term is
   the same expression, restated by Ba et al. as their Eq. (14).
   `report_enhanced_equilibrium` in `programs/unit_testing/lbm/solver` checks
   the two against each other over a sweep of sound speed, density and
   velocity; they agree to 1.1 × 10⁻¹⁴.

7. **A. Subhedar.** *Color-gradient lattice Boltzmann model for immiscible
   fluids with density contrast.* Physical Review E **106**, 045308, 2022.
   DOI: [10.1103/PhysRevE.106.045308](https://doi.org/10.1103/PhysRevE.106.045308)

   Velocity-based equilibrium; interface mobility independent of the density
   ratio. The velocity-based `droplet` solver is this kind of model, a
   colour-gradient segregation on a velocity-based hydrodynamics, and it sets
   its mobility on its own too, through the phase carrier's temperature
   (`SolverParameters::phase_temperature`; see
   [`numerics.md`](numerics.md#oscillations-against-exact-normal-modes)).

8. **S. Saito, N. Takada, S. Baba, S. Someya, H. Ito.** *Generalized equilibria
   for color-gradient lattice Boltzmann model based on higher-order Hermite
   polynomials.* Physical Review E **108**, 065305, 2023.
   [arXiv:2309.07801](https://arxiv.org/abs/2309.07801)

   Also the clearest recent survey of this model family. Read in full: its
   correction term is the one `Solver` takes as `S_Sp`, and at a density
   ratio of 1000 a rising bubble stays stable but rises at the wrong speed,
   which the authors put down to the finite-difference correction (see
   [`literature.md`](literature.md#what-limits-the-density-ratio)).

9. **S. Leclaire, M. Reggio, J.-Y. Trépanier.** *Progress and investigation on
   lattice Boltzmann modeling of multiple immiscible fluids or components with
   variable density and viscosity ratios.* Journal of Computational Physics
   **246**, 318–342, 2013.
   DOI: [10.1016/j.jcp.2013.03.039](https://doi.org/10.1016/j.jcp.2013.03.039)

   The first replacement of the raw colour gradient by a density-weighted one,
   which Ba et al. Eq. (21) simplifies.

10. **A. Q. Zhang, et al.** *A quantitative comparison of physical accuracy and
    numerical stability of lattice Boltzmann colour gradient and pseudopotential
    multicomponent models for microfluidic applications.*
    [arXiv:2110.05197](https://arxiv.org/abs/2110.05197)

    Places the colour-gradient family against Shan–Chen: a wider accessible
    range of density ratio, viscosity ratio and surface tension, and numerical
    stability at O(1000).

11. **M. Latva-Kokko, D. H. Rothman.** *Diffusion properties of gradient-based
    lattice Boltzmann models of immiscible fluids.* Physical Review E **71**,
    056702, 2005.
    DOI: [10.1103/PhysRevE.71.056702](https://doi.org/10.1103/PhysRevE.71.056702)

    The segregation operator every modern colour-gradient model recolours with,
    including Ba et al. Eq. (30).

    **Implemented twice here, to different ends.**

    In `Solver` it is `Recolouring::LatvaKokko`, and off by default: in that
    code's variables it is the operator already in `Solver::recolor()` with
    `beta * rho` in place of `p / (w cs^2)`, and that difference decides the
    scheme at high density contrast, because the pressure is continuous across
    an interface and the density is not. Measured better at a density ratio of
    2 (+0.3 % against +2.4 %) and unusable from 10 upwards. See
    [`numerics.md`](numerics.md#the-recolouring-why-p-and-not-rho).

    In `TwoPopulationSolver` it is the only recolouring, and it had to be
    adapted: pushing along the lattice weight `w_i` drives the light fluid's
    distribution negative by a factor of 240 at a density ratio of 1000, because
    the two fluids sit on different rest weights. Pushing along the mixture's
    own rest weight instead conserves each fluid's mass exactly and keeps both
    distributions non-negative for any `beta <= 1`, and reduces to Latva-Kokko
    when the rest weights agree. See
    [`numerics.md`](numerics.md#adapting-the-recolouring-to-a-density-ratio).

## Surface tension as a force

12. **B. Lafaurie, C. Nardone, R. Scardovelli, S. Zaleski, G. Zanetti.**
    *Modelling merging and fragmentation in multiphase flows with SURFER.*
    Journal of Computational Physics **113**, 134–147, 1994.

    The capillary stress tensor
    $\tfrac{\sigma}{2}(\lvert\nabla c\rvert\mathsf I - \nabla c\nabla c/\lvert\nabla c\rvert)$,
    whose divergence is the surface force. `src/lbm/surface_force.h`, used by
    `SurfaceTension::CapillaryStress` in place of the perturbation operator Ω⁽²⁾,
    and by the velocity-based solver.

13. **J. U. Brackbill, D. B. Kothe, C. Zemach.** *A continuum method for
    modeling surface tension.* Journal of Computational Physics **100**,
    335–354, 1992. And **S. V. Lishchuk, C. M. Care, I. Halliday.** *Lattice
    Boltzmann algorithm for surface tension with greatly reduced
    microcurrents.* Physical Review E **67**, 036701, 2003.

    The same force in curvature form, $\sigma\kappa\nabla c$, the second on the
    colour field of a colour-gradient model. `SurfaceTension::ContinuumSurfaceForce`
    is this form, after Ba et al.; it does not conserve momentum on the lattice
    (see [`numerics.md`](numerics.md#the-capillary-stress)).

14. **Z. Guo, C. Zheng, B. Shi.** *Discrete lattice effects on the forcing term
    in the lattice Boltzmann method.* Physical Review E **65**, 046308, 2002.

    The forcing scheme that `S_F` and the half-force velocity correction
    implement, and through which the surface force enters.

## Velocity-based scheme

15. **A. Fakhari, T. Mitchell, C. Leonardi, D. Bolster.** *Improved locality
    of the phase-field lattice-Boltzmann model for immiscible fluids at high
    density ratios.* Physical Review E **96**, 053301, 2017.

    The velocity-based hydrodynamic equilibrium that `src/lbm/velocity_based.h`
    follows. Its pressure correction is replaced there by the lattice gradient
    of $p$ itself, and its viscous correction by the link momentum exchange,
    which is what lets it run at $10^4$ and conserve momentum.

16. **Y. Q. Zu, S. He.** *Phase-field-based lattice Boltzmann model for
    incompressible binary fluid systems with density and viscosity
    contrasts.* Physical Review E **87**, 043301, 2013.

17. **P.-H. Chiu, Y.-T. Lin.** *A conservative phase field method for solving
    incompressible two-phase flows.* Journal of Computational Physics **230**,
    185–204, 2011.

    The conservative Allen–Cahn equation the phase populations solve.

18. **Z. Huang, G. Lin, A. M. Ardekani.** *Consistent and conservative scheme
    for incompressible two-phase flows using the conservative Allen–Cahn
    model.* Journal of Computational Physics **420**, 109718, 2020.

    The consistency conditions: the momentum is convected by the mass flux of
    the phase-field equation, Allen–Cahn part included. `link_momentum`
    carries the advected momentum with the mass flux of the phase
    populations, less a part quadratic in the velocity that vanishes for a
    uniform flow; see
    [`numerics.md`](numerics.md#momentum-conservation).

19. **S. Mirjalili, A. Mani.** *Consistent, energy-conserving momentum
    transport for simulations of two-phase flows using the phase field
    equations.* Journal of Computational Physics **426**, 109918, 2021.
    arXiv:1912.10096.

    The same principle for the conservative phase-field equation, with the
    momentum flux built from the discrete mass flux and a centred,
    kinetic-energy-conserving link velocity; and the observation that
    without it high density ratios need special treatment.

20. **M. Raessi, H. Pitsch.** *Consistent mass and momentum transport for
    simulating incompressible flows with large density ratios using the level
    set method.* Computers and Fluids **63**, 70–81, 2012.

21. **C. Zhan, Z. Chai, B. Shi.** *Consistent and conservative phase-field-based
    lattice Boltzmann method for incompressible two-phase flows.* Physical
    Review E **106**, 025319, 2022. arXiv:2111.00847.

    Consistent mass and momentum transport in a lattice Boltzmann model that
    streams $\rho\vec u$, run up to a density ratio of about 1000.

22. **H. Otomo, C. Sun, T. Inamuro, W. Li, M. Dressler, H. Chen, Y. Li,
    R. Zhang.** *Lattice Boltzmann models for the hydrodynamic equations in
    multiphase flow with high density ratio.* arXiv:2512.01027, 2025.

    Why models that stream $\rho\vec u$ lose accuracy next to a large density
    jump (truncation errors in $\vec u\,\partial^n\rho$ that spoil the
    momentum the dense fluid passes to the light one), and why models that
    stream $p/\rho$ depend on the absolute pressure; hence the pressure force
    on $p$ and the streaming of $\vec u$ here.

23. **O. Malaspinas.** *Increasing stability and accuracy of the lattice
    Boltzmann scheme: recursivity and regularization.* arXiv:1505.06900, 2015.

24. **J. Jacob, O. Malaspinas, P. Sagaut.** *A new hybrid recursive
    regularised Bhatnagar–Gross–Krook collision model for Lattice Boltzmann
    method-based large eddy simulation.* Journal of Turbulence **19**,
    1051–1076, 2018.

    The regularised and hybrid regularised collisions of `collide` and
    `collide_hybrid`.

25. **C. Kublik, R. Tsai.** *Integration over curves and surfaces defined by
    the closest point mapping.* Research in the Mathematical Sciences **3**, 3,
    2016. arXiv:1504.05478.

    Proposition 2: in 2D, the Jacobian $1 + \eta\kappa_\eta$ between a curve
    and its level set at distance $\eta$. `layer_weight` divides the capillary
    stress of each layer of a diffuse interface by it, which carries every
    layer's tension onto $\psi = 0$; see
    [`numerics.md`](numerics.md#surface-tension-at-the-interface-not-across-it).

26. **C. Kublik, N. M. Tanushev, R. Tsai.** *An implicit interface boundary
    integral method for Poisson's equation on arbitrary domains.* Journal of
    Computational Physics **247**, 279–311, 2013.

    The same Jacobian, used to write integrals over a curve as integrals over a
    band around it.

## High density ratios: the wider literature

These are the sources of [`literature.md`](literature.md). Except where an
entry above says otherwise, each is cited for what its abstract or Ba et al.
(2016) say of it; Ba et al. were read in full.

27. **P. Lallemand, L.-S. Luo.** *Theory of the lattice Boltzmann method:
    Dispersion, dissipation, isotropy, Galilean invariance, and stability.*
    Physical Review E **61**, 6546–6562, 2000.
    DOI: [10.1103/PhysRevE.61.6546](https://doi.org/10.1103/PhysRevE.61.6546)

    The moment basis of `src/lbm/mrt.h`, `--collision=mrt`.

28. **H. Huang, J.-J. Huang, X.-Y. Lu, M. C. Sukop.** *On simulations of
    high-density ratio flows using color-gradient multiphase lattice Boltzmann
    models.* International Journal of Modern Physics C **24**, 1350021, 2013.
    DOI: [10.1142/S0129183113500216](https://doi.org/10.1142/S0129183113500216)

    The unwanted term in the momentum equation of the rest-weight model, found
    on layered channel flow, and a source term cancelling it, validated to a
    density ratio of 8. Its evaluation "may lead to considerable numerical
    errors at higher density ratios" (Ba et al.'s summary), which the
    Taylor–Green measurements in `numerics.md` quantify.

29. **S. Leclaire, M. Reggio, J.-Y. Trépanier.** *Numerical evaluation of two
    recoloring operators for an immiscible two-phase flow lattice Boltzmann
    model.* Applied Mathematical Modelling **36**, 2237–2252, 2012.
    DOI: [10.1016/j.apm.2011.08.027](https://doi.org/10.1016/j.apm.2011.08.027)

    Latva-Kokko's recolouring adapted to variable density ratios.

30. **S. V. Lishchuk, I. Halliday, C. M. Care.** *Multicomponent lattice
    Boltzmann method for fluids with a density contrast.* Physical Review E
    **77**, 036702, 2008.
    DOI: [10.1103/PhysRevE.77.036702](https://doi.org/10.1103/PhysRevE.77.036702)

31. **K. Burgin, J. Spendlove, X. Xu, I. Halliday.** *Kinematics of
    chromodynamic multicomponent lattice Boltzmann simulation with a large
    density contrast.* Physical Review E **100**, 043310, 2019.
    DOI: [10.1103/PhysRevE.100.043310](https://doi.org/10.1103/PhysRevE.100.043310)

    That the restrictions on density contrast in rapid flow come from the terms
    correcting the dynamics for large density differences.

32. **J. Spendlove, X. Xu, O. J. Halliday, T. Schenkel, I. Halliday.**
    *Chromodynamic multirelaxation-time lattice Boltzmann scheme for fluids with
    density difference.* Physical Review E **102**, 013309, 2020.
    DOI: [10.1103/PhysRevE.102.013309](https://doi.org/10.1103/PhysRevE.102.013309)

33. **Z. X. Wen, Q. Li, Y. Yu, K. H. Luo.** *Improved three-dimensional
    color-gradient lattice Boltzmann model for immiscible two-phase flows.*
    Physical Review E **100**, 023301, 2019.
    DOI: [10.1103/PhysRevE.100.023301](https://doi.org/10.1103/PhysRevE.100.023301)

34. **R. Haghani, H. Erfani, J. E. McClure, E. G. Flekkøy, C. F. Berg.**
    *Color-gradient-based phase-field equation for multiphase flow.* Physical
    Review E **109**, 035301, 2024.
    DOI: [10.1103/PhysRevE.109.035301](https://doi.org/10.1103/PhysRevE.109.035301)

35. **T. Gregorczyk, S. Zhao, P. Boivin.** *Hybrid lattice Boltzmann method for
    multiphase flows.* Physical Review E **112**, 025311, 2025.
    DOI: [10.1103/gzwh-g7hz](https://doi.org/10.1103/gzwh-g7hz)

    Lafarge et al.'s group, moving to a pressure-based lattice Boltzmann scheme
    for the mixture and a finite-volume order parameter.

36. **T. Reis.** *A lattice Boltzmann formulation of the one-fluid model for
    multiphase flow.* Journal of Computational Physics **453**, 110962, 2022.
    DOI: [10.1016/j.jcp.2022.110962](https://doi.org/10.1016/j.jcp.2022.110962)

37. **J. Shao, C. Shu.** *A hybrid phase field multiple relaxation time lattice
    Boltzmann method for the incompressible multiphase flow with large density
    contrast.* International Journal for Numerical Methods in Fluids **77**,
    526–543, 2015. DOI: [10.1002/fld.3995](https://doi.org/10.1002/fld.3995)

38. **F. Hajabdollahi, K. N. Premnath, S. Welch.** *Central moment lattice
    Boltzmann method using a pressure-based formulation for multiphase flows at
    high density ratios and including effects of surface tension and Marangoni
    stresses.* Journal of Computational Physics **425**, 109893, 2021.
    DOI: [10.1016/j.jcp.2020.109893](https://doi.org/10.1016/j.jcp.2020.109893)

39. **T. Inamuro, T. Ogata, S. Tajima, N. Konishi.** *A lattice Boltzmann method
    for incompressible two-phase flows with large density differences.* Journal
    of Computational Physics **198**, 628–644, 2004.
    DOI: [10.1016/j.jcp.2004.01.019](https://doi.org/10.1016/j.jcp.2004.01.019)

40. **T. Lee, C.-L. Lin.** *A stable discretization of the lattice Boltzmann
    equation for simulation of incompressible two-phase flows at high density
    ratio.* Journal of Computational Physics **206**, 16–47, 2005.
    DOI: [10.1016/j.jcp.2004.12.001](https://doi.org/10.1016/j.jcp.2004.12.001)

41. **H. Liang, J. Xu, J. Chen, H. Wang, Z. Chai, B. Shi.** *Phase-field-based
    lattice Boltzmann modeling of large-density-ratio two-phase flows.* Physical
    Review E **97**, 033309, 2018.
    DOI: [10.1103/PhysRevE.97.033309](https://doi.org/10.1103/PhysRevE.97.033309)

    Droplet splashing at a density ratio of 1000, Re 20–500: the test Ba et
    al. run with a colour-gradient model at 100.

42. **B. Dorschner, F. Bösch, I. V. Karlin.** *Particles on demand for kinetic
    theory.* Physical Review Letters **121**, 130602, 2018. And
    **N. G. Kallikounis, B. Dorschner, I. V. Karlin.** *Particles on demand for
    flows with strong discontinuities.* Physical Review E **106**, 015301,
    2022. [arXiv:2203.05477](https://arxiv.org/abs/2203.05477)

    Populations whose lattice follows each node's velocity and temperature,
    which lifts the bound a fixed lattice puts on a cold fluid's speed
    (`T + u² ≥ |u|`, [`literature.md`](literature.md#what-limits-the-density-ratio)).
    Shown on compressible single-phase flows; not implemented here.

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
