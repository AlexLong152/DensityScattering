C     No dummy program needed

C-----------------------------------------------------------------------
C     SAID data structure and access functions for pion photoproduction
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     Maximum dimensions for data arrays
C-----------------------------------------------------------------------
      INTEGER FUNCTION MAX_ENERGY_POINTS()
      MAX_ENERGY_POINTS = 200
      RETURN
      END

      INTEGER FUNCTION MAX_ELL()
      MAX_ELL = 5
      RETURN
      END

C-----------------------------------------------------------------------
C     Data structure for multipole data
C     This is a Fortran-77 version of Python dictionaries used in the 
C     original code. We use separate arrays for "plus" and "minus" signs,
C     and for each target ("p12", "n12", "32q").
C
C     The arrays store:
C     - POLE_POINTS:   Number of energy points for each (sign,target,ell) combo
C     - POLE_SQRTS:    The sqrt(s) values for each energy point
C     - POLE_AMPL:     The complex amplitude values
C-----------------------------------------------------------------------
      BLOCK DATA INIT_SAID_DATA
      IMPLICIT NONE
      
      INTEGER MAX_ENG
      PARAMETER (MAX_ENG = 200)
      
      INTEGER MAX_L
      PARAMETER (MAX_L = 5)
      
C     Data dimensions
      INTEGER EPOLE_POINTS_PLUS_P12(0:MAX_L)
      INTEGER EPOLE_POINTS_PLUS_N12(0:MAX_L)
      INTEGER EPOLE_POINTS_PLUS_32Q(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_P12(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_N12(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_32Q(0:MAX_L)
      
      INTEGER MPOLE_POINTS_PLUS_P12(0:MAX_L)
      INTEGER MPOLE_POINTS_PLUS_N12(0:MAX_L)
      INTEGER MPOLE_POINTS_PLUS_32Q(0:MAX_L)
      INTEGER MPOLE_POINTS_MINUS_P12(0:MAX_L)
      INTEGER MPOLE_POINTS_MINUS_N12(0:MAX_L)
      INTEGER MPOLE_POINTS_MINUS_32Q(0:MAX_L)
      
C     Data arrays for E-pole
      DOUBLE PRECISION EPOLE_SQRTS_PLUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION EPOLE_SQRTS_PLUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION EPOLE_SQRTS_PLUS_32Q(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION EPOLE_SQRTS_MINUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION EPOLE_SQRTS_MINUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION EPOLE_SQRTS_MINUS_32Q(MAX_ENG,0:MAX_L)
      
      DOUBLE COMPLEX EPOLE_AMPL_PLUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX EPOLE_AMPL_PLUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX EPOLE_AMPL_PLUS_32Q(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX EPOLE_AMPL_MINUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX EPOLE_AMPL_MINUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX EPOLE_AMPL_MINUS_32Q(MAX_ENG,0:MAX_L)
      
C     Data arrays for M-pole
      DOUBLE PRECISION MPOLE_SQRTS_PLUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION MPOLE_SQRTS_PLUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION MPOLE_SQRTS_PLUS_32Q(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION MPOLE_SQRTS_MINUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION MPOLE_SQRTS_MINUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION MPOLE_SQRTS_MINUS_32Q(MAX_ENG,0:MAX_L)
      
      DOUBLE COMPLEX MPOLE_AMPL_PLUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX MPOLE_AMPL_PLUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX MPOLE_AMPL_PLUS_32Q(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX MPOLE_AMPL_MINUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX MPOLE_AMPL_MINUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX MPOLE_AMPL_MINUS_32Q(MAX_ENG,0:MAX_L)
      
      COMMON /EPOLE_DATA/ 
     &    EPOLE_POINTS_PLUS_P12, EPOLE_POINTS_PLUS_N12, 
     &    EPOLE_POINTS_PLUS_32Q, EPOLE_POINTS_MINUS_P12,
     &    EPOLE_POINTS_MINUS_N12, EPOLE_POINTS_MINUS_32Q,
     &    EPOLE_SQRTS_PLUS_P12, EPOLE_SQRTS_PLUS_N12,
     &    EPOLE_SQRTS_PLUS_32Q, EPOLE_SQRTS_MINUS_P12,
     &    EPOLE_SQRTS_MINUS_N12, EPOLE_SQRTS_MINUS_32Q,
     &    EPOLE_AMPL_PLUS_P12, EPOLE_AMPL_PLUS_N12,
     &    EPOLE_AMPL_PLUS_32Q, EPOLE_AMPL_MINUS_P12,
     &    EPOLE_AMPL_MINUS_N12, EPOLE_AMPL_MINUS_32Q
     
      COMMON /MPOLE_DATA/
     &    MPOLE_POINTS_PLUS_P12, MPOLE_POINTS_PLUS_N12, 
     &    MPOLE_POINTS_PLUS_32Q, MPOLE_POINTS_MINUS_P12,
     &    MPOLE_POINTS_MINUS_N12, MPOLE_POINTS_MINUS_32Q,
     &    MPOLE_SQRTS_PLUS_P12, MPOLE_SQRTS_PLUS_N12,
     &    MPOLE_SQRTS_PLUS_32Q, MPOLE_SQRTS_MINUS_P12,
     &    MPOLE_SQRTS_MINUS_N12, MPOLE_SQRTS_MINUS_32Q,
     &    MPOLE_AMPL_PLUS_P12, MPOLE_AMPL_PLUS_N12,
     &    MPOLE_AMPL_PLUS_32Q, MPOLE_AMPL_MINUS_P12,
     &    MPOLE_AMPL_MINUS_N12, MPOLE_AMPL_MINUS_32Q
      
C     Initialize all point counts to 0
      DATA EPOLE_POINTS_PLUS_P12 /0,0,0,0,0,0/
      DATA EPOLE_POINTS_PLUS_N12 /0,0,0,0,0,0/
      DATA EPOLE_POINTS_PLUS_32Q /0,0,0,0,0,0/
      DATA EPOLE_POINTS_MINUS_P12 /0,0,0,0,0,0/
      DATA EPOLE_POINTS_MINUS_N12 /0,0,0,0,0,0/
      DATA EPOLE_POINTS_MINUS_32Q /0,0,0,0,0,0/
      
      DATA MPOLE_POINTS_PLUS_P12 /0,0,0,0,0,0/
      DATA MPOLE_POINTS_PLUS_N12 /0,0,0,0,0,0/
      DATA MPOLE_POINTS_PLUS_32Q /0,0,0,0,0,0/
      DATA MPOLE_POINTS_MINUS_P12 /0,0,0,0,0,0/
      DATA MPOLE_POINTS_MINUS_N12 /0,0,0,0,0,0/
      DATA MPOLE_POINTS_MINUS_32Q /0,0,0,0,0,0/
      
      END

C-----------------------------------------------------------------------
C     Helper function to get E-pole amplitude for a given sqrtS
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION GET_EPOLE_FUNC(SIGN, TARGET, ELL, SQRTS)
      IMPLICIT NONE
      CHARACTER*(*) SIGN
      CHARACTER*3 TARGET
      INTEGER ELL
      DOUBLE PRECISION SQRTS
      
      INTEGER MAX_ENG, MAX_L
      PARAMETER (MAX_ENG = 200, MAX_L = 5)
      
C     Common blocks for pole data
      INTEGER EPOLE_POINTS_PLUS_P12(0:MAX_L)
      INTEGER EPOLE_POINTS_PLUS_N12(0:MAX_L)
      INTEGER EPOLE_POINTS_PLUS_32Q(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_P12(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_N12(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_32Q(0:MAX_L)
      
      DOUBLE PRECISION EPOLE_SQRTS_PLUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION EPOLE_SQRTS_PLUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION EPOLE_SQRTS_PLUS_32Q(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION EPOLE_SQRTS_MINUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION EPOLE_SQRTS_MINUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION EPOLE_SQRTS_MINUS_32Q(MAX_ENG,0:MAX_L)
      
      DOUBLE COMPLEX EPOLE_AMPL_PLUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX EPOLE_AMPL_PLUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX EPOLE_AMPL_PLUS_32Q(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX EPOLE_AMPL_MINUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX EPOLE_AMPL_MINUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX EPOLE_AMPL_MINUS_32Q(MAX_ENG,0:MAX_L)
      
      COMMON /EPOLE_DATA/ 
     &    EPOLE_POINTS_PLUS_P12, EPOLE_POINTS_PLUS_N12, 
     &    EPOLE_POINTS_PLUS_32Q, EPOLE_POINTS_MINUS_P12,
     &    EPOLE_POINTS_MINUS_N12, EPOLE_POINTS_MINUS_32Q,
     &    EPOLE_SQRTS_PLUS_P12, EPOLE_SQRTS_PLUS_N12,
     &    EPOLE_SQRTS_PLUS_32Q, EPOLE_SQRTS_MINUS_P12,
     &    EPOLE_SQRTS_MINUS_N12, EPOLE_SQRTS_MINUS_32Q,
     &    EPOLE_AMPL_PLUS_P12, EPOLE_AMPL_PLUS_N12,
     &    EPOLE_AMPL_PLUS_32Q, EPOLE_AMPL_MINUS_P12,
     &    EPOLE_AMPL_MINUS_N12, EPOLE_AMPL_MINUS_32Q
      
      INTEGER POINTS, I, IDX_MIN
      DOUBLE PRECISION MIN_DIFF, DIFF
      
C     Check if ELL is in valid range
      IF (ELL .LT. 0 .OR. ELL .GT. MAX_L) THEN
        GET_EPOLE_FUNC = (0.0D0, 0.0D0)
        RETURN
      ENDIF
      
C     Get array dimensions and data based on SIGN and TARGET
      IF (SIGN .EQ. 'plus') THEN
        IF (TARGET .EQ. 'p12') THEN
          POINTS = EPOLE_POINTS_PLUS_P12(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_EPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(EPOLE_SQRTS_PLUS_P12(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(EPOLE_SQRTS_PLUS_P12(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_EPOLE_FUNC = EPOLE_AMPL_PLUS_P12(IDX_MIN,ELL)
          
        ELSE IF (TARGET .EQ. 'n12') THEN
          POINTS = EPOLE_POINTS_PLUS_N12(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_EPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(EPOLE_SQRTS_PLUS_N12(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(EPOLE_SQRTS_PLUS_N12(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_EPOLE_FUNC = EPOLE_AMPL_PLUS_N12(IDX_MIN,ELL)
          
        ELSE IF (TARGET .EQ. '32q') THEN
          POINTS = EPOLE_POINTS_PLUS_32Q(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_EPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(EPOLE_SQRTS_PLUS_32Q(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(EPOLE_SQRTS_PLUS_32Q(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_EPOLE_FUNC = EPOLE_AMPL_PLUS_32Q(IDX_MIN,ELL)
        ELSE
          GET_EPOLE_FUNC = (0.0D0, 0.0D0)
        ENDIF
      ELSE IF (SIGN .EQ. 'minus') THEN
        IF (TARGET .EQ. 'p12') THEN
          POINTS = EPOLE_POINTS_MINUS_P12(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_EPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(EPOLE_SQRTS_MINUS_P12(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(EPOLE_SQRTS_MINUS_P12(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_EPOLE_FUNC = EPOLE_AMPL_MINUS_P12(IDX_MIN,ELL)
          
        ELSE IF (TARGET .EQ. 'n12') THEN
          POINTS = EPOLE_POINTS_MINUS_N12(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_EPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(EPOLE_SQRTS_MINUS_N12(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(EPOLE_SQRTS_MINUS_N12(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_EPOLE_FUNC = EPOLE_AMPL_MINUS_N12(IDX_MIN,ELL)
          
        ELSE IF (TARGET .EQ. '32q') THEN
          POINTS = EPOLE_POINTS_MINUS_32Q(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_EPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(EPOLE_SQRTS_MINUS_32Q(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(EPOLE_SQRTS_MINUS_32Q(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_EPOLE_FUNC = EPOLE_AMPL_MINUS_32Q(IDX_MIN,ELL)
        ELSE
          GET_EPOLE_FUNC = (0.0D0, 0.0D0)
        ENDIF
      ELSE
        GET_EPOLE_FUNC = (0.0D0, 0.0D0)
      ENDIF
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Helper function to get M-pole amplitude for a given sqrtS
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION GET_MPOLE_FUNC(SIGN, TARGET, ELL, SQRTS)
      IMPLICIT NONE
      CHARACTER*(*) SIGN
      CHARACTER*3 TARGET
      INTEGER ELL
      DOUBLE PRECISION SQRTS
      
      INTEGER MAX_ENG, MAX_L
      PARAMETER (MAX_ENG = 200, MAX_L = 5)
      
C     Common blocks for pole data
      INTEGER MPOLE_POINTS_PLUS_P12(0:MAX_L)
      INTEGER MPOLE_POINTS_PLUS_N12(0:MAX_L)
      INTEGER MPOLE_POINTS_PLUS_32Q(0:MAX_L)
      INTEGER MPOLE_POINTS_MINUS_P12(0:MAX_L)
      INTEGER MPOLE_POINTS_MINUS_N12(0:MAX_L)
      INTEGER MPOLE_POINTS_MINUS_32Q(0:MAX_L)
      
      DOUBLE PRECISION MPOLE_SQRTS_PLUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION MPOLE_SQRTS_PLUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION MPOLE_SQRTS_PLUS_32Q(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION MPOLE_SQRTS_MINUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION MPOLE_SQRTS_MINUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE PRECISION MPOLE_SQRTS_MINUS_32Q(MAX_ENG,0:MAX_L)
      
      DOUBLE COMPLEX MPOLE_AMPL_PLUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX MPOLE_AMPL_PLUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX MPOLE_AMPL_PLUS_32Q(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX MPOLE_AMPL_MINUS_P12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX MPOLE_AMPL_MINUS_N12(MAX_ENG,0:MAX_L)
      DOUBLE COMPLEX MPOLE_AMPL_MINUS_32Q(MAX_ENG,0:MAX_L)
      
      COMMON /MPOLE_DATA/
     &    MPOLE_POINTS_PLUS_P12, MPOLE_POINTS_PLUS_N12, 
     &    MPOLE_POINTS_PLUS_32Q, MPOLE_POINTS_MINUS_P12,
     &    MPOLE_POINTS_MINUS_N12, MPOLE_POINTS_MINUS_32Q,
     &    MPOLE_SQRTS_PLUS_P12, MPOLE_SQRTS_PLUS_N12,
     &    MPOLE_SQRTS_PLUS_32Q, MPOLE_SQRTS_MINUS_P12,
     &    MPOLE_SQRTS_MINUS_N12, MPOLE_SQRTS_MINUS_32Q,
     &    MPOLE_AMPL_PLUS_P12, MPOLE_AMPL_PLUS_N12,
     &    MPOLE_AMPL_PLUS_32Q, MPOLE_AMPL_MINUS_P12,
     &    MPOLE_AMPL_MINUS_N12, MPOLE_AMPL_MINUS_32Q
      
      INTEGER POINTS, I, IDX_MIN
      DOUBLE PRECISION MIN_DIFF, DIFF
      
C     Check if ELL is in valid range
      IF (ELL .LT. 0 .OR. ELL .GT. MAX_L) THEN
        GET_MPOLE_FUNC = (0.0D0, 0.0D0)
        RETURN
      ENDIF
      
C     Get array dimensions and data based on SIGN and TARGET
      IF (SIGN .EQ. 'plus') THEN
        IF (TARGET .EQ. 'p12') THEN
          POINTS = MPOLE_POINTS_PLUS_P12(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_MPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(MPOLE_SQRTS_PLUS_P12(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(MPOLE_SQRTS_PLUS_P12(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_MPOLE_FUNC = MPOLE_AMPL_PLUS_P12(IDX_MIN,ELL)
          
        ELSE IF (TARGET .EQ. 'n12') THEN
          POINTS = MPOLE_POINTS_PLUS_N12(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_MPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(MPOLE_SQRTS_PLUS_N12(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(MPOLE_SQRTS_PLUS_N12(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_MPOLE_FUNC = MPOLE_AMPL_PLUS_N12(IDX_MIN,ELL)
          
        ELSE IF (TARGET .EQ. '32q') THEN
          POINTS = MPOLE_POINTS_PLUS_32Q(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_MPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(MPOLE_SQRTS_PLUS_32Q(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(MPOLE_SQRTS_PLUS_32Q(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_MPOLE_FUNC = MPOLE_AMPL_PLUS_32Q(IDX_MIN,ELL)
        ELSE
          GET_MPOLE_FUNC = (0.0D0, 0.0D0)
        ENDIF
      ELSE IF (SIGN .EQ. 'minus') THEN
        IF (TARGET .EQ. 'p12') THEN
          POINTS = MPOLE_POINTS_MINUS_P12(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_MPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(MPOLE_SQRTS_MINUS_P12(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(MPOLE_SQRTS_MINUS_P12(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_MPOLE_FUNC = MPOLE_AMPL_MINUS_P12(IDX_MIN,ELL)
          
        ELSE IF (TARGET .EQ. 'n12') THEN
          POINTS = MPOLE_POINTS_MINUS_N12(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_MPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(MPOLE_SQRTS_MINUS_N12(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(MPOLE_SQRTS_MINUS_N12(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_MPOLE_FUNC = MPOLE_AMPL_MINUS_N12(IDX_MIN,ELL)
          
        ELSE IF (TARGET .EQ. '32q') THEN
          POINTS = MPOLE_POINTS_MINUS_32Q(ELL)
          
          IF (POINTS .EQ. 0) THEN
            GET_MPOLE_FUNC = (0.0D0, 0.0D0)
            RETURN
          ENDIF
          
C         Find the closest sqrtS value
          MIN_DIFF = ABS(MPOLE_SQRTS_MINUS_32Q(1,ELL) - SQRTS)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(MPOLE_SQRTS_MINUS_32Q(I,ELL) - SQRTS)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
          GET_MPOLE_FUNC = MPOLE_AMPL_MINUS_32Q(IDX_MIN,ELL)
        ELSE
          GET_MPOLE_FUNC = (0.0D0, 0.0D0)
        ENDIF
      ELSE
        GET_MPOLE_FUNC = (0.0D0, 0.0D0)
      ENDIF
      
      RETURN
      END