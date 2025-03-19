C-----------------------------------------------------------------------
C     Multipole data initialization for SAID pion photoproduction
C
C     This file defines and initializes the data structures used to store
C     multipole amplitudes for pion photoproduction. It implements a 
C     Fortran-77 version of the Python dictionary structure from the original
C     code, optimized for efficient access and memory usage.
C
C     Data Organization:
C     ------------------
C     The multipole data is organized by:
C
C     1. Multipole Type:
C        - E: Electric multipoles (E0+, E1+, E1-, etc.)
C        - M: Magnetic multipoles (M1+, M1-, M2+, etc.)
C
C     2. Sign (Total Angular Momentum relation):
C        - plus: J = L + 1/2 multipoles
C        - minus: J = L - 1/2 multipoles
C
C     3. Target:
C        - p12: proton target (isospin 1/2, charge +1)
C        - n12: neutron target (isospin 1/2, charge 0)  
C        - 32q: isospin 3/2 channel
C
C     4. Angular Momentum (ell):
C        - 0: S-waves (e.g., E0+, S11)
C        - 1: P-waves (e.g., M1+, P33)
C        - 2: D-waves
C        - 3: F-waves
C        - 4: G-waves
C
C     Array Dimensions:
C     ----------------
C     The arrays use a consistent indexing scheme:
C     - POINTS(0:MAX_ELL, 3, 2): Number of energy points for each combination
C       * First index: angular momentum (0-5)
C       * Second index: target (1=p12, 2=n12, 3=32q)
C       * Third index: sign (1=plus, 2=minus)
C
C     - SQRTS & AMPL arrays include an additional first dimension for energy points
C
C     This structure is significantly more memory-efficient than the previous
C     approach, requiring only 6 arrays instead of 36.
C-----------------------------------------------------------------------
      BLOCK DATA INIT_SAID_DATA
      IMPLICIT NONE
      
C     Include constants
      INCLUDE 'parseConstants.def'
      
C-----------------------------------------------------------------------
C     Data dimensions - number of energy points for each multipole combination
C     POINTS(0:MAX_ELL, 3, 2) array where:
C     - First index: angular momentum ell (0-5)
C     - Second index: target type (1=p12, 2=n12, 3=32q)
C     - Third index: sign (1=plus, 2=minus)
C-----------------------------------------------------------------------
      INTEGER EPOLE_POINTS(0:MAX_ELL, 3, 2)
      INTEGER MPOLE_POINTS(0:MAX_ELL, 3, 2)
      
C-----------------------------------------------------------------------
C     Data arrays for sqrt(s) values (center-of-mass energy in MeV)
C     SQRTS(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2) where:
C     - First index: energy point index (up to MAX_ENERGY_POINTS)
C     - Second index: angular momentum ell (0-5)
C     - Third index: target type (1=p12, 2=n12, 3=32q)
C     - Fourth index: sign (1=plus, 2=minus)
C-----------------------------------------------------------------------
      DOUBLE PRECISION EPOLE_SQRTS(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      DOUBLE PRECISION MPOLE_SQRTS(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      
C-----------------------------------------------------------------------
C     Data arrays for complex amplitude values
C     AMPL(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2) where dimensions match SQRTS
C     Stores the complex multipole amplitudes in units of MeV^-1
C-----------------------------------------------------------------------
      DOUBLE COMPLEX EPOLE_AMPL(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      DOUBLE COMPLEX MPOLE_AMPL(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      
C-----------------------------------------------------------------------
C     Common blocks for sharing data across routines
C     These common blocks allow the data to be accessed by other routines
C-----------------------------------------------------------------------
      COMMON /EPOLE_DATA/ EPOLE_POINTS, EPOLE_SQRTS, EPOLE_AMPL
      COMMON /MPOLE_DATA/ MPOLE_POINTS, MPOLE_SQRTS, MPOLE_AMPL
      
C-----------------------------------------------------------------------
C     Initialize all point counts to 0
C     All entries in POINTS arrays are set to 0 initially and will be 
C     populated by PARSE_SAID_FILE when reading the data
C-----------------------------------------------------------------------
      DATA EPOLE_POINTS /36*0/
      DATA MPOLE_POINTS /36*0/
      
      END