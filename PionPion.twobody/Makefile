# Use the first available shell: zsh, bash, or sh
SHELL := $(shell which zsh 2>/dev/null || which bash 2>/dev/null || echo /bin/sh)

# Executable name
COMMAND = twobodyvia2Ndensity.PionPion

# Directory paths
THIS_DIR   := .
COMMON_DIR := ../common-densities
VAR_MODS   := $(COMMON_DIR)/density-modules
TWOBODY_VS := ../twobodyvia2Ndensity

# Add source search paths
VPATH := $(THIS_DIR) $(COMMON_DIR) $(TWOBODY_VS)

# Object files in this folder
OBJS_THIS := \
  calculateQs.PionPion.o \
  readinput.twobody.PionPion.o \
  2Bkernel.PionPion.o \
  2Bspinsym.PionPion.o \
  2Bspinasy.PionPion.o \
  usesymmetries.PionPion.o

# Object files from twobodyvia2Ndensity
OBJS_TWOBODY := \
  calculate2BI2.o \
  finalstatesums.twobodyvia2Ndensity.o \
  LebedevLaikov.o \
  main.twobodyvia2Ndensity.o \
  read2Ndensity.o \
  setquads.o \
  spinstructures.o

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

# Compiler and flags
FC = mpif90
FFLAGS = -O3 -fopenmp -g -c -Wall -Wno-tabs -ffixed-line-length-0 \
  -ffpe-summary=none -mcmodel=large -fPIC -cpp -fcheck=all \
  -finit-integer=-1000 -finit-real=inf -fallow-argument-mismatch \
  -I$(VAR_MODS) -I$(COMMON_DIR) -I/usr/local/hdf5/openmpi/include

LDFLAGS = \
  -L$(VAR_MODS) -ldensity \
  -L/usr/local/hdf5/openmpi/lib \
  -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
  -llapack -lgomp -ldl -lz -lsz

.PHONY: all clean

all: $(COMMAND)
	chmod +x $(COMMAND)

$(COMMAND): $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	find . -name '*.o' -o -name '*.mod' -delete
	rm -f $(COMMAND)
