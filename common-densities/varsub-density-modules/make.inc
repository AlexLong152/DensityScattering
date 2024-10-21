# define positions of submodules required 

UNAME := $(shell uname -s)


#include $(BASE)/compile-options/dir.inc

# define the compilers and standard libraries for the current computer 

LIBSTD1 = -llapack
LIBSTD2 = -L /usr/lib/x86_64-linux-gnu/hdf5/openmpi/ 
LIBSTD3 = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
LIBSTD4 = -llapack -lgomp -ldl -lz -lsz

#ifeq ($(USER),hgrie)
#     LIBSTD2 = -L /usr/lib/x86_64-linux-gnu/hdf5/openmpi/ 
#else ifeq ($(USER),alex) # your username
#     LIBSTD2 =	 #whatever your path is
#     #     LIBSTD2 = -L /usr/lib/x86_64-linux-gnu/hdf5/openmpi/ 
#     LIBSTD3 = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
#     LIBSTD4 = -llapack -lgomp -ldl -lz -lsz
#else ifeq ($(USER),along) # your username
#     LIBSTD2 = -L /usr/lib/x86_64-linux-gnu/hdf5/openmpi/ 
#     LIBSTD3 = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
#     LIBSTD4 = -llapack -lgomp -ldl -lz -lsz
#endif	
#     LIBSTD3 = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
#     LIBSTD4 = -llapack -lgomp -ldl -lz -lsz

LIBSTD5 = 
LIBSTD6 = 
LIBSTD7 = 
LIBSTD8 =

LIBSTD = $(LIBSTD1) $(LIBSTD2) $(LIBSTD3) $(LIBSTD4) $(LIBSTD5) $(LIBSTD6) $(LIBSTD7) $(LIBSTD8)

CC = mpicc
COPT = -g -O3 -fopenmp
FC = mpif90 
F90 = mpif90
COMPILE_OPT1= -fopenmp -cpp -fcheck=all  -g -O5 
# -fno-automatic
COMPILE_OPT2= -ffixed-line-length-0

COMPILE_OPT1= -fPIC  -cpp  -g -O3 -ffpe-summary=none  -finit-integer=-1000 -finit-real=inf -fcheck=all    
COMPILE_OPT2= -ffixed-line-length-0 -fopenmp

COMPILE_OPT3= -I /usr/include/hdf5/openmpi/	

COMPILE_OPT4= 
COMPILE_OPT4= 
COMPILE_OPT5= 
COMPILE_OPT6= 
COMPILE_OPT7= 
COMPILE_OPT8= 

FOPT = $(COMPILE_OPT1) $(COMPILE_OPT2) $(COMPILE_OPT3) $(COMPILE_OPT4) \
       $(COMPILE_OPT5) $(COMPILE_OPT6) $(COMPILE_OPT7) $(COMPILE_OPT8) \
      -DVERDATE="'$(VERDATE)'" -DVERREV="'$(VERREV)'" \
      -DMAKEF90="'$(F90)'" -DMAKEF77="'$(FC)'" \
      -DMAKEFOPT1="'$(COMPILE_OPT1)'" \
      -DMAKEFOPT2="'$(COMPILE_OPT2)'" \
      -DMAKEFOPT3="'$(COMPILE_OPT3)'" \
      -DMAKEFOPT4="'$(COMPILE_OPT4)'" \
      -DMAKEFOPT5="'$(COMPILE_OPT5)'" \
      -DMAKEFOPT6="'$(COMPILE_OPT6)'" \
      -DMAKEFOPT7="'$(COMPILE_OPT7)'" \
      -DMAKEFOPT8="'$(COMPILE_OPT8)'" \
      -DMAKECOPT="'$(COPT)'" -DMAKECC="'$(CC)'" \
      -DMAKELIBSTD1="'$(LIBSTD1)'" \
      -DMAKELIBSTD2="'$(LIBSTD2)'" \
      -DMAKELIBSTD3="'$(LIBSTD3)'" \
      -DMAKELIBSTD4="'$(LIBSTD4)'" \
      -DMAKELIBSTD5="'$(LIBSTD5)'" \
      -DMAKELIBSTD6="'$(LIBSTD6)'" \
      -DMAKELIBSTD7="'$(LIBSTD7)'" \
      -DMAKELIBSTD8="'$(LIBSTD8)'" \
      -DCOMPDATE="'$(COMPDATE)'"
FOPT+=-fallow-argument-mismatch
AR     = ar
ARFLAGS= rv 
RANLIB   = touch
#  -ffpe-trap=invalid,zero,overflow  -fcheck=all -fno-automatic -finit-integer=-1000 -finit-real=inf -fbounds-check\
# -ffpe-trap=invalid,zero,overflow  -fcheck=all -fno-automatic -finit-integer=-1000 -finit-real=zero\
#-fdefault-real-8 -fno-automatic  -ffpe-trap=invalid,zero,overflow  -fbounds-check -fcheck=all -fcheck=bounds 
#   -ffpe-trap=invalid,zero,overflow \
#      -finit-integer=-1000 -finit-real=inf\
#      -fcheck=all \
