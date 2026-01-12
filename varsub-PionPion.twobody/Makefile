# Use the first available shell: zsh, bash, or sh
SHELL := $(shell which zsh 2>/dev/null || which bash 2>/dev/null || echo /bin/sh)

# Executable name
COMMAND = twobodyvia2Ndensity.PionPion
TEST_COMMAND = testIsospin

# Directory paths
THIS_DIR   := .
COMMON_DIR := ../common-densities
VAR_MODS   := $(COMMON_DIR)/varsub-density-modules
TWOBODY_VS := ../varsub-twobodyvia2Ndensity

# Add source search paths
VPATH := $(THIS_DIR) $(COMMON_DIR) $(TWOBODY_VS) tests

# Object files in this folder
OBJS_THIS := \
  varsub-calculateQs.PionPion.o \
  varsub-readinput.twobody.PionPion.o \
  varsub-2Bkernel.PionPion.o \
  varsub-2Bspinsym.PionPion.o \
  varsub-2Bspinasy.PionPion.o \
  varsub-usesymmetries.PionPion.o \
  varsub-IsospinMap.PionPion.o

# Object files from varsub-twobodyvia2Ndensity
OBJS_TWOBODY := \
  varsub-calculate2BI2.o \
  varsub-finalstatesums.twobodyvia2Ndensity.o \
  varsub-LebedevLaikov.o \
  varsub-main.twobodyvia2Ndensity.o \
  varsub-read2Ndensity.o \
  varsub-setquads.o \
  varsub-spinstructures.o

# Object files from common-densities
OBJS_COMMON := \
  andreasquads.o \
  sphereCartConvert.o \
  gauleg.o \
  getylm.o \
  gridGauleg.o \
  makedensityfilename.o \
  makeoutputfilename.o \
  outputtomath.o \
  photonmomenta.o \
  readinput-densities.o \
  outputroutine.o

# All objects
OBJS := $(OBJS_THIS) $(OBJS_TWOBODY) $(OBJS_COMMON)

# Test objects
TEST_OBJS := testIsospin.o varsub-IsospinMap.PionPion.o

# Compiler and flags
FC = mpif90

# FFLAGS = -O0 -fopenmp -g -Wall -Wno-tabs -ffixed-line-length-0 \
#   -ffpe-summary=None -fPIC -cpp -fcheck=all -fbacktrace \
#   -finit-integer=-1000 -finit-real=inf -fallow-argument-mismatch\
#   -I$(VAR_MODS) -I$(COMMON_DIR) -I/usr/local/hdf5/openmpi/include

FFLAGS  = -O3 -fopenmp -Wall -Wno-tabs -ffixed-line-length-0 \
  -ffpe-summary=none -fPIC -cpp -fbacktrace \
  -finit-integer=-1000 -finit-real=inf -fallow-argument-mismatch \
  -I$(VAR_MODS) -I$(COMMON_DIR) -I/usr/local/hdf5/openmpi/include

debug_FLAGS  = -O0 -fopenmp -Wall -Wno-tabs -ffixed-line-length-0 \
  -ffpe-summary=none -fPIC -cpp -fbacktrace \
  -finit-integer=-1000 -finit-real=inf -fallow-argument-mismatch \
  -I$(VAR_MODS) -I$(COMMON_DIR) -I/usr/local/hdf5/openmpi/include

LDFLAGS = \
  -L$(VAR_MODS) -ldensity \
  -L/usr/local/hdf5/openmpi/lib \
  -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
  -llapack -lgomp -ldl -lz -lsz

.PHONY: all clean fast test 

all: $(COMMAND)
	chmod +x $(COMMAND)


$(COMMAND): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(TEST_COMMAND): $(TEST_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

test: $(TEST_COMMAND)
	chmod +x $(TEST_COMMAND)
	./$(TEST_COMMAND)

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

# Build an optimized version with FAST_FLAGS (same executable name)
debug: clean
	$(MAKE) -B FFLAGS="$(debug_FLAGS)" all

clean:
	find . \( -name '*.o' -o -name '*.mod' \) -delete
	rm -f $(COMMAND) $(TEST_COMMAND)
