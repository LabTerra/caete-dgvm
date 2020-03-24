FC = gfortran
PYEXEC = python3.8
PYFILE = create_pls_table.py
SO_BUILD = ./build.sh
CLEAN_DIR = ./clean_dir.sh

##SO
MFLAG = -m
HFLAG = -h
CFLAG = -c
OVRTFLAG = --overwrite-signature
F2PY = python3.8 -m numpy.f2py
MODNAME = caete_module
INTERFACES = $(MODNAME).pyf

FCFLAGS = -g -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fbounds-check -Wconversion -pedantic
FLFLAGS = -c -g -Wall -fcheck=all -Wextra -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fbounds-check -Wconversion -pedantic

sources = global.f90 utils.f90 funcs.f90 soil_pools.f90 soil_dec.f90 productivity.f90 budget.f90 caete.f90

objects = global.o utils.o funcs.o soil_pools.o soil_dec.o productivity.o budget.o test_caete.o

test: $(objects)
	$(SO_BUILD)
	$(PYEXEC) $(PYFILE)
	$(FC) -o test $(FCFLAGS) $(objects)

test1: $(objects)
	$(FC) -o test $(FCFLAGS) $(objects)

interfaces: $(sources)
	$(F2PY) $(HFLAG) $(INTERFACES) $(sources) $(MFLAG) $(MODNAME) $(OVRTFLAG)

so: $(sources) interfaces
	$(F2PY) $(INTERFACES) $(CFLAG) $(sources)

global.o: global.f90
	$(FC) $(FLFLAGS) global.f90

types.mod: global.o global.f90
	$(FC) $(FLFLAGS) global.f90

global_par.mod: global.o global.f90
	$(FC) $(FLFLAGS) global.f90

photo_par.mod: global.o global.f90
	$(FC) $(FLFLAGS) global.f90

utils.o: utils.f90
	$(FC) $(FLFLAGS) utils.f90

utils.mod: utils.o utils.f90
	$(FC) $(FLFLAGS) utils.f90

funcs.o: funcs.f90
	$(FC) $(FLFLAGS) funcs.f90

photo.mod: funcs.o funcs.f90
	$(FC) $(FLFLAGS) funcs.f90

water.mod: funcs.o funcs.f90
	$(FC) $(FLFLAGS) funcs.f90

soil_pools.o: soil_pools.f90
	$(FC) $(FLFLAGS) soil_pools.f90

soil_pools.mod: soil_pools.o soil_pools.f90
	$(FC) $(FLFLAGS) soil_pools.f90

soil_dec.o: soil_dec.f90
	$(FC) $(FLFLAGS) soil_dec.f90

soil_dec.mod: soil_dec.o soil_dec.f90
	$(FC) $(FLFLAGS) soil_dec.f90

productivity.o: productivity.f90
	$(FC) $(FLFLAGS) productivity.f90

productivity.mod: productivity.o productivity.f90
	$(FC) $(FLFLAGS) productivity.f90

budget.o: budget.f90
	$(FC) $(FLFLAGS) budget.f90

budget.mod: budget.o budget.f90
	$(FC) $(FLFLAGS) budget.f90

test_caete.o: test_caete.f90
	$(FC) $(FLFLAGS) test_caete.f90

clean:
	rm -rf test $(objects) *.mod *.s *.so *.pyf pls_ex.txt
	$(CLEAN_DIR)

clean_python:
	rm -rf __pycache__
