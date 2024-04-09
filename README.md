# CutFEM for phase transitions in solids

_Modelling of propagation of phase boundaries in non-linear solids in 2D._

## Summary

The code represents the implementation of the cut-finite-element method[^1] (CutFEM) to solve 2D problems with time-dependent interfaces. It solves PDEs describing the mechanical deformation of a solid and, optionally, the diffusion of the gaseous-type chemical reactant. The PDEs are defined on domains with non-stationary phase boundaries. The kinetics of the latter depends on the solution of the PDEs. 

The code was originally written for Ref.[^2]. It subsequently evolved over the years, acquiring numerous improvements. The main highlights of the code are:
- the mechanical problem is fully non-linear (i.e. large deformations, arbitrary rheology), the total Lagrange formulation is used;
- hyperelastic and visco-elastic constitutive models are implemented, as well as linear elasticity;
- the topology change of the phase boundaries is handled automatically;
- the code has minimalist design, but is capable of parallel computations.

The code is written in MATLAB.

## Getting started

The user can run the code by executing the main file and providing the folder name where the results should be stored: `mainPhaseTran('tmp1');`

After the calculations are finished (the default example takes approximately 100 seconds on a laptop with 11th Gen Intel Core i5), the user can run the postprocessing script to create a video file visualising the solution: `animateFront;`

Parallel computations can be optionally switched on by uncommenting the corresponding lines in `calcMechResPhaseTran.m`

The same file incorporates the implementation of the boundary conditions (BCs) applied externally to the computational domain (not to be confused with the interface conditions). BCs are imposed strongly, replacing the equations for the predefined sets of boundary nodes.

Alternative constitutive models can be switched on in `constLaw.m`

For understanding the problem parameters and the setup of the example, the user is referred to Ref.[^3].

## Developer and acknowledgements

The code was developed by Dr. M. Poluektov during the research fellowship at the University of Warwick. 

The author would like to acknowledge the help of Prof. G. Kreiss (Uppsala University) in understanding the CutFEM approach and the help of Prof. A.B. Freidin (St. Petersburg Polytechnic University) in understanding the physics of phase transitions in solids. 

The author is also grateful to Dr. Ł. Figiel (University of Warwick) for discussions on computational mechanics and to Dr. A. Morozov (Technische Universität Berlin) for discussions on stability of phase boundaries. 

[^1]: E. Burman _et al._, _Internat. J. Numer. Methods Engrg._ 104(7):472-501, 2015, [link](https://onlinelibrary.wiley.com/doi/10.1002/nme.4823).
[^2]: M. Poluektov and Ł. Figiel, _Comput. Mech._ 63:885-911, 2019, [link](https://link.springer.com/article/10.1007/s00466-018-1628-z).
[^3]: A. Morozov _et al._, _Eur. J. Mech. A/Solids_ 104:105211, 2024, [link](https://www.sciencedirect.com/science/article/pii/S0997753823003030).