#!/bin/bash

## apply the numpy.f2py tool to CAETE fortan code

## Curious person: - How to compile the CAETE fortran code?
## Spider man:  - This way:

## Spider man: - First of all, generate the file with modules and procedures interfaces
python -m numpy.f2py -h caete_module.pyf global.f90 utils.f90  soil_pools.f90 soil_dec.f90\
 funcs.f90 productivity.f90 budget.f90 caete.f90 -m caete_module --overwrite-signature

## Spider man: - Finally, you compile the model (as a shared object or dll)
python -m numpy.f2py caete_module.pyf -c global.f90 utils.f90 soil_pools.f90 soil_dec.f90\
 funcs.f90 productivity.f90 budget.f90 caete.f90
