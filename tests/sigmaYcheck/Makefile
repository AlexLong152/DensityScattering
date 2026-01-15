# Makefile for sigmaYCheck test program
# Compares sigma_y implementations between singlesigmasym and IsospinMap

# Executable name
COMMAND = sigmaYCheck

# Directory paths
THIS_DIR   := .
TWOBODY_VS := ../../varsub-twobodyvia2Ndensity
PIPI_DIR   := ../../varsub-PionPion.twobody

# Add source search paths
VPATH := $(THIS_DIR) $(TWOBODY_VS) $(PIPI_DIR)

# Object files
OBJS := \
  sigmaYCheck.o \
  varsub-spinstructures.o \
  varsub-IsospinMap.PionPion.o

# Compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall -Wno-tabs -ffixed-line-length-0

.PHONY: all run clean

all: $(COMMAND)
	chmod +x $(COMMAND)

$(COMMAND): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

run: $(COMMAND)
	./$(COMMAND)

clean:
	rm -f $(COMMAND) *.o *.mod
