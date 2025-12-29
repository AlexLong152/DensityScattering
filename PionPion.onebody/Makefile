SHELL   := /bin/sh

# ─── compilers ────────────────────────────────────────────────────────────────
FC      := mpif90

# ─── where is your HDF5? ──────────────────────────────────────────────────────
HDF5_DIR := /usr/local/hdf5/openmpi

# ─── compiler flags ───────────────────────────────────────────────────────────
FFLAGS  := -O -g -fbacktrace -Wall -fcheck=all -c -Wno-tabs -ffixed-form -ffixed-line-length-none \
           -ffpe-summary=none -mcmodel=large \
           -I$(HDF5_DIR)/include \
           -I../common-densities/density-modules \
           -I../common-densities

# FFLAGS = -O1 -fopenmp -g -Wall -Wno-tabs -ffixed-line-length-0 \
#   -ffpe-summary=none -fPIC -cpp -fcheck=all -fbacktrace \
#   -finit-integer=-1000 -finit-real=nan -fallow-argument-mismatch \
#   -I/usr/local/hdf5/openmpi/include\
#   -I$(HDF5_DIR)/include \
#   -I../common-densities/density-modules \
#   -I../common-densities

# ─── linker flags & libraries ─────────────────────────────────────────────────
LDFLAGS := -Wl,-rpath,$(HDF5_DIR)/lib \
           -L../common-densities/density-modules -ldensity \
           -L$(HDF5_DIR)/lib \
             -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
           -llapack -lgomp -ldl -lz -lsz

# ─── object files ─────────────────────────────────────────────────────────────
OBJS    := \
    ../onebodyvia1Ndensity/calcAmpVars.o \
    ../onebodyvia1Ndensity/calculateAs.OQ3poles.o \
    ../onebodyvia1Ndensity/calculateAs.OQ3nopoles.o \
    ../onebodyvia1Ndensity/calculateAs.Delta.o \
    ../onebodyvia1Ndensity/construct1NAmps.o \
    ../onebodyvia1Ndensity/read1Ndensity.o \
    ../onebodyvia1Ndensity/trafo1Namps-to-helicity.o \
    ../onebodyvia1Ndensity/readinput-onebody.o \
    ../common-densities/gauleg.o \
    ../common-densities/gridGauleg.o \
    ../common-densities/makedensityfilename.o \
    ../common-densities/makeoutputfilename.o \
    ../common-densities/outputtomath.o \
    ../common-densities/photonmomenta.o \
    ../common-densities/readinput-densities.o \
    ../common-densities/outputroutine.o \
    main.onebodyvia1Ndensity.o \
    utility_suite.o \
    parseFile.o \
    PionScat.o

# ─── include‐file dependencies ─────────────────────────────────────────────────
INCS    := \
    ../common-densities/constants.def \
    ../common-densities/params.def \
    ../common-densities/calctype.def

# ─── final executable ─────────────────────────────────────────────────────────
COMMAND := onebodyvia1Ndensity

.PHONY: all clean test testMat testGetMat testGetCS

all: $(COMMAND)

$(COMMAND): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)
	chmod +x $@

# Define test objects excluding the main program object
TEST_OBJS := $(filter-out main.onebodyvia1Ndensity.o, $(OBJS))


test: testGetMat.o $(TEST_OBJS)
	$(FC) -o testGetMat testGetMat.o $(TEST_OBJS) $(LDFLAGS)
	chmod +x testGetMat

testGetMat: testGetMat.o $(TEST_OBJS)
	$(FC) -o testGetMat testGetMat.o $(TEST_OBJS) $(LDFLAGS)
	chmod +x testGetMat

testGetCS: testGetCS.o $(TEST_OBJS)
	$(FC) -o testGetCS testGetCS.o $(TEST_OBJS) $(LDFLAGS)
	chmod +x testGetCS

# Special rule for PionScat that depends on parseFile module
PionScat.o: PionScat.f $(INCS) parseFile.o
	$(FC) $(FFLAGS) $< -o $@

# Special rule for parseFile module (no special dependencies)
parseFile.o: parseFile.f $(INCS)
	$(FC) $(FFLAGS) $< -o $@

# Special rule for main program that depends on PionScat module
main.onebodyvia1Ndensity.o: main.onebodyvia1Ndensity.f $(INCS) PionScat.o
	$(FC) $(FFLAGS) $< -o $@

# Special rule for testGetMat that depends on PionScat module
testGetMat.o: testGetMat.f PionScat.o
	$(FC) $(FFLAGS) $< -o $@

# Special rule for testGetCS that depends on PionScat module
testGetCS.o: testGetCS.f PionScat.o
	$(FC) $(FFLAGS) $< -o $@

# compile any .f → .o if sources or .def files change
%.o: %.f $(INCS)
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(COMMAND) testCrossSec testCrossSec.o testGetMat testGetMat.o testGetCS testGetCS.o *.mod
