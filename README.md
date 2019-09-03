# CAETÊ

## __AUTHORS__:

 - Bianca Rius
 - David Lapola
 - Helena Alves
 - João Paulo Darela Filho
 - Put you name here! (labterra@unicamp.br)

## SOFTWARE REQUIREMENTS

To run the model in development mode:

- Python3 + numpy + netCDF4 + (libs for plotting things)
- gfortran >= 5

## FORTRAN FILES

        global.f90
        funcs.f90
        productivity.f90
        budget.f90
        caete.f90

- global.f90 <- contains several modules.
  - types: Real, integer and logical standard types defined in gfortran+amd64.
  - global_par: global parameters.
  - photo_par: constants and parameters for the photosynthesis models.
  - soil_dec: new decomposition module, contains the subroutine carbon3 and others.

- funcs.f90 <- contains two modules.
  - photo: definition of several functions related with ecophysiological/ecosystem function and model management.
  - water: definition of several functions related with the hydrologic cycle.

- productivity.f90 <- Contains one module with one subroutine.
  - productivity: definition of subrotine prod that calculates the primary productivity using some functions of photo(func.f90) module.

- budget.f90 <- Contains one module (budget) with one subroutine.
  - budget: definition of subroutine daily_budget which calculates (a) daily productivity using prod subroutine. and (b) daily water budget (hydrology) using funtions of the water(funcs.f90) module.

- caete.f90 <- Contains one module (caete) with one subroutine.
  - caete: contains the definition of caete_dyn that execute the model trought a timespan.

## PYTHON FILES

- caete.py: main python script with functions that define the execution of the model. Currently the only classes and funcions related with model application in individual gridcells are implemented. (contém gambiarras)

- plsgen.py: python script that creates the table of CAETÊ Plant Life Strategies.

- test_mod.py: Script that executes the model to some gridcells. Contains the base structure for model execution to be implemented in caete.py

## Other python files
 
 - create_input.py: how to use the caete.data_in class for generate/organize the environmental input data for CAETÊ
 
 - write_output.py: old version- Made to save the outputs of CAETE-PVM
 
 - caete_driver.py old version - Made to drive the execution of CAETE-PVM runs and save the models outputs

## MODEL INPUTS

GAMBIARRA

## Compile the module

Compile caete_module running the script `build.sh`:

`$ ./build.sh`

See the `build.sh` script.

All the fortran files described are compiled using the f2py tool that is part of numpy. Details can be found here: [f2py user guide](https://docs.scipy.org/doc/numpy-1.11.0/f2py/index.html).


The resulting module (caete_module) is 'importable' in python environments.
