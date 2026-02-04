# Use the first available shell: zsh, bash, or sh
SHELL := $(shell which zsh 2>/dev/null || which bash 2>/dev/null || echo /bin/sh)

# Executable target name
COMMAND = onebodyvia1Ndensity

# Source directories
SRC_DIR      = ../onebodyvia1Ndensity
THIS_DIR     = .#PionPhotoProd.onebody
COMMON_DIR   = ../common-densities
MOD_DIR      = $(COMMON_DIR)/density-modules

# Add search paths for source files
VPATH := $(SRC_DIR) $(THIS_DIR) $(COMMON_DIR)

# Object files from onebodyvia1Ndensity
OBJS_ONEBODY = \
  $(SRC_DIR)/calcAmpVars.o \
  $(SRC_DIR)/calculateAs.OQ3poles.o \
  $(SRC_DIR)/calculateAs.OQ3nopoles.o \
  $(SRC_DIR)/calculateAs.Delta.o \
  $(SRC_DIR)/construct1NAmps.o \
  $(SRC_DIR)/read1Ndensity.o \
  $(SRC_DIR)/trafo1Namps-to-helicity.o

# Object files from common-densities
OBJS_COMMON = \
  $(COMMON_DIR)/gauleg.o \
  $(COMMON_DIR)/gridGauleg.o \
  $(COMMON_DIR)/makedensityfilename.o \
  $(COMMON_DIR)/makeoutputfilename.o \
  $(COMMON_DIR)/outputtomath.o \
  $(COMMON_DIR)/photonmomenta.o \
  $(COMMON_DIR)/readinput-densities.o \
  $(COMMON_DIR)/outputroutine.o

# Object files from PionPhotoProd.onebody
OBJS_THIS = \
  $(THIS_DIR)/main.onebodyvia1Ndensity.o \
  $(THIS_DIR)/readinput-onebody.o \
	$(THIS_DIR)/outputKernel.o

OBJS = $(OBJS_ONEBODY) $(OBJS_COMMON) $(OBJS_THIS)

# Compiler and flags
FC = mpif90
FFLAGS = -O -g -c -Wall -Wno-tabs -ffixed-line-length-132 \
  -ffpe-summary=none -mcmodel=large -fPIC -cpp -fcheck=all \
  -finit-integer=-1000 -finit-real=inf -fallow-argument-mismatch \
  -ffixed-form \
  -I$(MOD_DIR) -I$(COMMON_DIR) -I/usr/local/hdf5/openmpi/include

LDFLAGS = \
  -L$(MOD_DIR) -ldensity \
  -L/usr/local/hdf5/openmpi/lib \
  -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
  -llapack -lgomp -ldl -lz -lsz

.PHONY: all clean

all: $(COMMAND)
	   chmod +x $(COMMAND)

$(COMMAND): $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)
	chmod +x $@

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	find . -name '*.o' -o -name '*.mod' -delete
	rm -f $(COMMAND)
