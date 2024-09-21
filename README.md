# BioFOAM
#### The BioFOAM solver was developed based on the micro-continuum theory and the Darcy-Brinkman-Stokes (DBS) equation to simulate pore-scale biofilm growth in heterogeneous porous media. This code is developed based on the pioneer work conducted by Cyprien Soulaine from the package named porousMedia4Foam (https://github.com/csoulain/porousMedia4Foam).

### 1. Check the environment

The BioFOAM solver is built up in OpenFOAM Version 7. The required environment is GNU Linux. Before installing the code, please source the OpenFOAM environment:
```
source ~/etc/bashrc
```

### 2. Install the solver

1) enter the correct path ```cd BioFOAM/solver/```
2) run  ```wmake```
3) if you want to clear the solver, run ```wclean```

### Please note that the solver requires debugging process when it is first installed in a new environment.

### 3. Run a case

1) enter the correct path ```cd BioFOAM/cases/SinglePore/```
2) run  ```blockMesh```
3) run ```setField```
3) run ```BioFOAM```
