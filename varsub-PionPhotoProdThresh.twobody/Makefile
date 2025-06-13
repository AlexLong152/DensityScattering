# Use the first available shell: zsh, bash, or sh
SHELL := $(shell which zsh 2>/dev/null || which bash 2>/dev/null || echo /bin/sh)

# Executable name
COMMAND = twobodyvia2Ndensity.PionPhotoProdThresh

# Directory paths
THIS_DIR   := .
COMMON_DIR := ../common-densities
VAR_MODS   := $(COMMON_DIR)/varsub-density-modules
TWOBODY_VS := ../varsub-twobodyvia2Ndensity

# Add source search paths
VPATH := $(THIS_DIR) $(COMMON_DIR) $(TWOBODY_VS)

# Object files in this folder
OBJS_THIS := \
  varsub-calculateQs.PionPhotoProdThresh.o \
  varsub-readinput.twobody.PionPhotoProdThresh.o \
  varsub-2Bkernel.PionPhotoProdThresh.o \
  varsub-2Bspinsym.PionPhotoProdThresh.o \
  varsub-2Bspinasy.PionPhotoProdThresh.o \
  varsub-usesymmetries.PionPhotoProdThresh.o

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
