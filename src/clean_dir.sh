#!/bin/sh

rm -rf __pycache__/
#rm -rf exec_log.txt
#rm -rf caete_module.pyf
#rm -rf caete_module.cpython-35m-i386-linux-gnu.so
#rm -rf caete_module.cpython-35m-x86_64-linux-gnu.so
#rm -rf pls.bin
#rm -rf pls.txt
rm -rf pls_attrs.csv
rm -rf *.mod
rm -rf pools_din*

# COMPILE FORTRAn MODULES OF INTEREST

gfortran -Wall -g -S global.f90 soil_dec.f90 funcs.f90 productivity.f90 budget.f90 caete.f90

rm -rf *.s