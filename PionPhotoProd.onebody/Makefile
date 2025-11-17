SHELL   := /bin/sh

# ─── compilers ────────────────────────────────────────────────────────────────
FC      := mpif90

# ─── where is your HDF5? ──────────────────────────────────────────────────────
HDF5_DIR := /usr/local/hdf5/openmpi

# ─── compiler flags ───────────────────────────────────────────────────────────
FFLAGS  := -O -g -c -Wall -Wno-tabs -ffixed-form -ffixed-line-length-none \
           -ffpe-summary=none -mcmodel=large \
           -I$(HDF5_DIR)/include \
           -I../common-densities/density-modules \
           -I../common-densities

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
    ../common-densities/gauleg.o \
    ../common-densities/gridGauleg.o \
    ../common-densities/makedensityfilename.o \
    ../common-densities/makeoutputfilename.o \
    ../common-densities/outputtomath.o \
    ../common-densities/photonmomenta.o \
    ../common-densities/readinput-densities.o \
    ../common-densities/outputroutine.o \
    main.onebodyvia1Ndensity.o \
    poles.o \
    readinput-onebody.o \
    said_subs.o \
    subroutines.o \
    utility_suite.o

# ─── include‐file dependencies ─────────────────────────────────────────────────
INCS    := \
    ../common-densities/constants.def \
    ../common-densities/params.def \
    ../common-densities/calctype.def

# ─── final executable ─────────────────────────────────────────────────────────
COMMAND := onebodyvia1Ndensity

.PHONY: all clean test check

all: $(COMMAND)

$(COMMAND): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)
	chmod +x $@

# Define test objects excluding the main program object
TEST_OBJS := $(filter-out main.onebodyvia1Ndensity.o, $(OBJS))

# test: testCrossSec.o $(TEST_OBJS)
# 	$(FC) -o testCrossSec testCrossSec.o $(TEST_OBJS) $(LDFLAGS)
# 	chmod +x testCrossSec

# Define check objects excluding the main program object
CHECK_OBJS := $(filter-out main.onebodyvia1Ndensity.o, $(OBJS))

test: testing.o $(CHECK_OBJS)
	$(FC) -o TestCrossSec testing.o $(CHECK_OBJS) $(LDFLAGS)
	chmod +x TestCrossSec

# compile any .f → .o if sources or .def files change
%.o: %.f $(INCS)
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -f $(OBJS) testing.o $(COMMAND) testCrossSec testCrossSec.o TestCrossSec
