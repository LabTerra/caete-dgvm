FC = gfortran
FCFLAGS = -Wall -Wextra -g -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -Wconversion -pedantic
FLFLAGS = -c -g -Wall -Wextra -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -Wconversion -pedantic

sources = global.f90 utils.f90 funcs.f90 soil_dec.f90 test_carbon3.f90

objects = global.o utils.o funcs.o soil_dec.o test_carbon3.o


test: $(objects)
	$(FC) -o test $(FCFLAGS) $(objects)

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

soil_dec.o: soil_dec.f90
	$(FC) $(FLFLAGS) soil_dec.f90

soil_dec.mod: soil_dec.o soil_dec.f90
	$(FC) $(FLFLAGS) soil_dec.f90

test_carbon3.o: test_carbon3.f90
	$(FC) $(FLFLAGS) test_carbon3.f90

clean:
	rm -rf test $(objects) *.mod *.s *.so *.pyf
