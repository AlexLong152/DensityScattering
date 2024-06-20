
SHELL=/bin/sh
# hgrie Nov 2023: uses user-specific path to hd5 library
#	hgrie May 2018: rewritten for hdf5 and parallel processing (mpi)
# This makefile is used as follows ([option] last):
# make -f makefile		-> compile with FFLAGS (usually non-optimized)
# make -f makefile global	-> compile with O3FFLAGS (optimized)
# make -f makefile findincl	-> List include files found in sources
# make -f makefile showincl	-> Write dependencies of OBJS on include files
# make -f makefile tar		-> tar and compress
# make -f makefile untar	-> uncompress and extract files from tar
#
# Edit: OBJS, INCS, INPUTS, COMMAND, INCLUDE DEPENDENCIES, tar(?).
#
# Comments:
#   It will compile the fortran modules indicated by SRCS into the object
# modules indicated by OBJS and then produce an executable indicated by
# COMMAND.  Note that objects are not recompiled unless they are out of date
# with respect to their dependencies (sources and includes).
# If the SRCS are in different directories, the objects will
# be placed in same directory as source and includes will be referenced
# correctly.  Includes may not be referenced correctly with global;
# it not, change FFLAGS to O3FFLAGS and make with no option.
#
#   makefile represents the name of this makefile.  If it is actually named 
# makefile, then the -f makefile is unnecessary.  Edit the list of .o files
# (OBJS) below to represent your files.  The .f files (SRCS) can be listed 
# separately, but it is easier to determine them directly from the OBJS.
# The include files and dependence of OBJS on include files can be 
# determined using options findincl and showincl, and the results can be
# pasted into this file (and edited if needed).  Edit INPUTS to list other
# files needed to run COMMAND (*.dat, etc); INPUTS is only used with tar 
# option but it doesnt hurt to document the files anyway.
# Note that a line starting with a space follows every rule (called a 
# separator), and \ can be used to continue lines.
#  
#   The tar option creates the file $(COMMAND).tar containing this makefile,
# SRCS, INCS (include files), and INPUTS. (any input or other files that
# are needed to run code).  Tar preserves dates and is a good way to move
# code between systems.  You may wish to edit the shell commands for the
# tar option directly to include *.f, etc., as needed.
#
#   This makefile is a modified version of the standard makefile.
#   Options void, clean, and clobber which appear in the usual makefile 
# have been removed here.

# OBJS defined below, modify as needed.  OBJS are the source files (.f)
# with .f -> .o .   Directories are relative to make directory.
OBJS=\
varsub-calculateQs.PionPhotoProdThresh.o \
varsub-readinput.twobody.PionPhotoProdThresh.o \
varsub-2Bkernel.PionPhotoProdThresh.o \
varsub-2Bspinsym.PionPhotoProdThresh.o \
varsub-2Bspinasy.PionPhotoProdThresh.o \
varsub-usesymmetries.PionPhotoProdThresh.o \
../varsub-twobodyvia2Ndensity/varsub-calculate2BI2.o \
../varsub-twobodyvia2Ndensity/varsub-finalstatesums.twobodyvia2Ndensity.o \
../varsub-twobodyvia2Ndensity/varsub-LebedevLaikov.o \
../varsub-twobodyvia2Ndensity/varsub-main.twobodyvia2Ndensity.o \
../varsub-twobodyvia2Ndensity/varsub-read2Ndensity.o \
../varsub-twobodyvia2Ndensity/varsub-setquads.o \
../varsub-twobodyvia2Ndensity/varsub-spinstructures.o \
../common-densities/andreasquads.o \
../common-densities/sphereCartConvert.o \
../common-densities/gauleg.o \
../common-densities/getylm.o \
../common-densities/gridGauleg.o \
../common-densities/makedensityfilename.o \
../common-densities/makeoutputfilename.o \
../common-densities/outputtomath.o \
../common-densities/photonmomenta.o \
../common-densities/readinput-densities.o \
../common-densities/outputroutine.o 

#../twobody/calculate2BI2.o \

# Source files (.f) = Object files (.o) with .o -> .f
SRCS = $(OBJS:.o=.f)

# objects, source files and directory to deal with .mod files 
MODS =

MODSRC = $(MODS:.o=.f)

MODDIR = ../common-densities/varsub-density-modules/

# Include files as generated with findincl option and modified as needed;
# any directory reference should be relative to make directory, not the
# directory of .f file which can be different.
# Include files assumed to end with .inc or .def.
# NO continuation (\) ON LAST LINE.
#
INCS=\
../common-densities/constants.def \
../common-densities/params.def \
../common-densities/calctype.def 

# Define other files needed with COMMAND for taring and documentation.
#
INPUTS= 

# Compiler parameters
#
# CFT		Name of the fortran compiler to use [f77]
# FFLAGS	Flags to the fortran compiler [use -O1 when debugging]
# O3FFLAGS      Flags to the fortran compiler with -O3 option
# LDFLAGS	Flags to the loader
# LIBS		A list of libraries 
# COMMAND	Name of the executable
#
CFT=  mpif90
#FFLAGS=  -static -O -g -C -extend_source
#O3FFLAGS=  -static -O3 -extend_source
# FFLAGS and O3FLAGS input via include.  File .FFLAGS is machine dependent.
#include .FFLAGS
#include $(HOME)/.FFLAGS
# Ensure that .mod files are FIRST read from dir common-densities/density-modules/,
# and then from common-densities/. Necessary since CompDens needs newer version of
#	spline.f than the onebody and twobody parts of the code....
# hgrie May 2019: included new flag -mcmodel=large to allow for larger arrays
# FFLAGS= -O -fopenmp -g -c -Wall -Wno-tabs -ffixed-line-length-132  \
# 	-ffpe-summary=none  -mcmodel=large \
# 	-I../common-densities/density-modules -I../common-densities 
#Along June 2024 optimized compiling with -O3 and changed fixed line length to infinite
# -O3 change my break debugging with GDB debugger needs to be done with -O0
FFLAGS= -O3 -fopenmp -g -c -Wall -Wno-tabs -ffixed-line-length-0  \
	-ffpe-summary=none  -mcmodel=large \
	-I../common-densities/varsub-density-modules -I../common-densities 

UNAME := $(shell echo $USER)

ifeq ($(USER),hgrie)
	FFLAGS += -I /usr/include/hdf5/openmpi/
else ifeq ($(USER),alex) # your username
	FFLAGS += -I /usr/include/ #whatever your path is
else ifeq ($(USER),along) # your username
	FFLAGS += -I /usr/include/ #whatever your path is
endif

# library for density
     LIBSTD1 = -L../common-densities/varsub-density-modules -ldensity
     LIBSTD2 = -llapack
     LIBSTD3 = 
     LIBSTD4 = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
     LIBSTD5 = -llapack -lgomp -ldl -lz -lsz
     LIBSTD6 = 
     LIBSTD7 = 
     LIBSTD8 = 

LDFLAGS= $(LIBSTD1) $(LIBSTD2) $(LIBSTD3) $(LIBSTD4) $(LIBSTD5) $(LIBSTD6) $(LIBSTD7) $(LIBSTD8)

ifeq ($(USER),hgrie)
	LDFLAGS += -L /usr/lib/x86_64-linux-gnu/hdf5/openmpi/ 
else ifeq ($(USER),alex) # your username
#LDFLAGS += -L /usr/lib/  -I /usr/include/ #whatever your path is
	LDFLAGS += -L /usr/lib/openmpi # -I /usr/include/ #whatever your path is
else ifeq ($(USER),along) # your username
#	LDFLAGS += -L /usr/lib/  -I /usr/include/ #whatever your path is
#   LDFLAGS += -L /usr/lib/  -I /usr/include/ #whatever your path is
endif

LIBS=  
COMMAND= twobodyvia2Ndensity.PionPhotoProdThresh

#
# Lines from here on are the actual rules that make uses to build $(COMMAND).
# Definition of command must precede include dependencies!

all: $(MODS) MODFILES $(COMMAND)

#Compile the module files separately, tell the Makefile where to find them
#
#
#MODFILES: $(MODS)
#	gfortran -I$(MODDIR) -c $(MODSRC)

$(COMMAND):	$(OBJS)
	$(CFT) -o $(COMMAND) $(OBJS) $(MODS) $(LDFLAGS) $(LIBS) 
	chmod go+x $(COMMAND)

#
# INCLUDE DEPENDENCIES: Dependencies on include files, Edit as needed.  
# Can be generated via showincl option; the directory of include file
# must be relative to make directory (check it).  May need to break long 
# lines with \ .

#main.o readinput.o setquads.o: params.def
#main.o readinput.o : calctype.def
#main.o setquads.o calcmomenta.o calcmatrix.o polnamplitudes.o: constants.def

#
#
# Lines from here on down should rarely need to be changed.
#
# Rule for dependence of .o on .f files.
# Within subshell, cd to the directory of .f file and compile it.
# This leaves the .o file in the .f files directory and insures that
# include files are referenced from the .f files directory and not the
# directory make is run from (if different).
#
.f.o:
	( cd $(@D) ; $(CFT) -c $(FFLAGS) $(<F) )
# Simply moving the .o file does not solve the include reference problem:
#	@ if [ "$(*F)" != "$*" ]; then mv $(@F) $@ ; fi
 
# global should fail if .f file in different directory has include reference
# that is not also valid from this directory.
global:
	$(CFT) $(O3FFLAGS) -o $(COMMAND) $(SRCS) $(LDFLAGS) $(LIBS) 
	chmod go+x $(COMMAND)
 
# tar -- makefile is name of this makefile.
# choice 1) tar everything of use in THIS directory.
# choice 2) tar files related to COMMAND, from whatever directory.
tar:
#	tar cf - makefile .FFLAGS *.f *.inc *.def make* input.dat > $(COMMAND).tar
	tar cf - makefile .FFLAGS $(SRCS) $(INCS) $(INPUTS) > $(COMMAND).tar
	compress $(COMMAND).tar
 
untar:
	uncompress $(COMMAND).tar
	tar xf $(COMMAND).tar
 
# If you want your .o files sorted:
sortobjs:
	@ echo $(OBJS) | sort
 
findincl:
	@ echo "These includes were found in sources.  Paste into makefile."
	@ echo "Dont forget to remove \\ on last line and add"
	@ echo "  separator (line after last include starting with space)."
	@ echo "Also may need to add/change directory (relative to make) before file."
	@ echo "INCS=\\ "
	@ grep -i 'include ' $(SRCS) | awk '{print $$3}' | \
	  grep "'" | sed "s/'//;s/'/ \\\/" | sort -u
 
# showincl uses .DEFAULT rule.  There may be a cleaner way to do this?
showincl: $(INCS:.def=.de) $(INCS:.inc=.in)
	@ echo "The above dependencies on INCS (includes) were found in SRCS (.f).
	@ echo "Paste into makefile.  Adjust wrapped lines."
 
.DEFAULT:
	@ echo `grep -l $(@F) $(SRCS) | sed "s/\\.f/.o/"` : \
	  $@ | sed "s/\.in/\.inc/" | sed "s/\.de/\.def/"
 
# added by hgrie: clean up, delete .mod files in common-densities and twobody

clean:
#	rm *.mod *.o $(OBJS) $(COMMAND) $(MODS) $(MODDIR)/*.mod
	rm *.mod *.o $(OBJS) $(COMMAND)

