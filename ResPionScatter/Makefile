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

.PHONY: all clean test

all: $(COMMAND)

$(COMMAND): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)
	chmod +x $@

# Define test objects excluding the main program object
TEST_OBJS := $(filter-out main.onebodyvia1Ndensity.o, $(OBJS))

test: testCrossSec.o $(TEST_OBJS)
	$(FC) -o testCrossSec testCrossSec.o $(TEST_OBJS) $(LDFLAGS)
	chmod +x testCrossSec

# Special rule for main program that depends on PionScat module
main.onebodyvia1Ndensity.o: main.onebodyvia1Ndensity.f $(INCS) PionScat.o
	$(FC) $(FFLAGS) $< -o $@

# Special rule for parseFile that depends on PionScat module
parseFile.o: parseFile.f $(INCS) PionScat.o
	$(FC) $(FFLAGS) $< -o $@

# compile any .f → .o if sources or .def files change
%.o: %.f $(INCS)
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(COMMAND) testCrossSec testCrossSec.o
