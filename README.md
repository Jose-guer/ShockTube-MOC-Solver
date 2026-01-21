# ShockTube-MOC-Solver
This code implements a method of characteristics solver for standard shocktube problems in MATLAB. NASA9 polynomials are used to obtained thermodyanmic properties such as the specific heats ratio (gamma) and molecular weights used by the solver.

## System requirements
No specific version of MATLAB is required to run the code.

## Description of Code
- ShockTube_MOC_Solve.m - This is the actual MOC solver which takes in the shocktube geometry and initial conditions.
- ShockTube_MOC_Plot.m - This is the wrapper function used to run ShockTube_MOC_Solve.m and print post-incident and post-refelcted states as well as test time duration.

## Mixtures
-Air
-Argon
-Helium
-Hydrogen
-Nitrogen
-Helium-Argon mixtures

Note that for mixtures, the variable "x" which is the mole fraction of each consituent needs to be updated in the mixture thermodat file.
