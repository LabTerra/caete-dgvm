#!/bin/sh

rm -rf __pycache__/
rm -rf caete_module.pyf
rm -rf caete_module.cpython*
rm -rf pls_attrs.csv
rm -rf pools_din*

# COMPILE FORTRAn MODULES OF INTEREST

#gfortran -Wall -g -S global.f90 utils.f90 soil_dec.f90 funcs.f90\
# productivity.f90 budget.f90 caete.f90
