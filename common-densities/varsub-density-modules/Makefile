objects= precision.o constants.o memcheck_c.o mpi_const.o gauss.o spline.o clebsch.o besselfu.o hdf_tool.o parallel.o pwbook.o meshpoints.o meshscatt.o amplitudes.o CompDens.o  

include make.inc

LIBDENS=libdensity.a

LIB= $(LIBDENS) $(LIBSTD) $(LIBDENS)

CFLAGS=$(COPT) 
FFLAGS=$(FOPT) 

# Reset the default goal.
.DEFAULT_GOAL := testcompdens

testcompdens:  $(LIBDENS) testcompdens.o 
	$(FC) -o testcompdens $(LIB) testcompdens.o $(LIB)

produce2Ndensityformath:  $(LIBDENS) produce2Ndensityformath.o 
	$(FC) -o produce2Ndensityformath $(LIB) produce2Ndensityformath.o $(LIB)

$(LIBDENS):  $(objects)
	$(AR) $(ARFLAGS) $(LIBDENS) $(objects)
	$(RANLIB) $(LIBDENS)


.PHONY : clean
clean: 
	rm -f  *.mod *.o *.a  testcompdens produce2Ndensityformath

.SUFFIXES:
.SUFFIXES: .F90 .f90 .F .f .c .o

.F90.o          :
	$(F90) $(FFLAGS) -c $<

.f90.o          :
	$(F90) $(FFLAGS) -c $<

.F.o          :
	$(FC) $(FFLAGS) -c $<

.f.o          :
	$(FC) $(FFLAGS) -c $<

.c.o          :
	$(CC) $(CFLAGS) -c $<
