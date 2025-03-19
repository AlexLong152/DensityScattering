C-----------------------------------------------------------------------
C     Test program for pion photoproduction
C-----------------------------------------------------------------------
      PROGRAM TEST_GETRAWM
      IMPLICIT NONE
      
C     External function declaration
      DOUBLE PRECISION GET_M_SQUARE
      EXTERNAL GET_M_SQUARE
      
C     Local variables
      DOUBLE PRECISION SQRTS, THETA_DEG, THETA_RAD
      CHARACTER*3 NUCS
      DOUBLE PRECISION MSQUARE1, MSQUARE2
      INTEGER BASIC_CHECK, I
      DOUBLE COMPLEX M(2,2)
      
C     Debug variables
      INTEGER MAX_ENG, MAX_L
      PARAMETER (MAX_ENG = 200, MAX_L = 5)
      
C     Common blocks for pole data
      INTEGER EPOLE_POINTS_PLUS_P12(0:MAX_L)
      INTEGER EPOLE_POINTS_PLUS_N12(0:MAX_L)
      INTEGER EPOLE_POINTS_PLUS_32Q(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_P12(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_N12(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_32Q(0:MAX_L)
      
      COMMON /EPOLE_DATA/ 
     &    EPOLE_POINTS_PLUS_P12, EPOLE_POINTS_PLUS_N12, 
     &    EPOLE_POINTS_PLUS_32Q, EPOLE_POINTS_MINUS_P12,
     &    EPOLE_POINTS_MINUS_N12, EPOLE_POINTS_MINUS_32Q
      
C     Set parameters (matching Python basicCheck function)
      SQRTS = 1100.0D0
      THETA_DEG = 40.0D0
      THETA_RAD = THETA_DEG * 3.14159265358979D0 / 180.0D0
      NUCS = 'pp0'
      
C     Initialize physics constants and SAID data
      CALL INITIALIZE_PION_PHOTO('said-SM22.txt')
      
C     Get the raw matrix element
      CALL GET_RAW_M(SQRTS, THETA_RAD, NUCS, M)
      
      
C     Calculate the first method using GET_M_SQUARE
      MSQUARE1 = GET_M_SQUARE(SQRTS, THETA_RAD, NUCS)
      
C     Calculate the second method (same as first, just to demonstrate both methods)
      CALL BASIC_CHECK_2(SQRTS, THETA_RAD, NUCS, MSQUARE2)
      
      
C     Print the results, matching the Python output format
      WRITE(*,*) 'Msquare1=', MSQUARE1
      WRITE(*,*) '--------------------------------------------------'
      WRITE(*,*) 'MSquare2=', MSQUARE2
      
      END
      
C-----------------------------------------------------------------------
C     Implementation of basicCheck function from the Python code
C-----------------------------------------------------------------------
      SUBROUTINE BASIC_CHECK_2(SQRTS, THETA, NUCS, MSQUARE)
      IMPLICIT NONE
      DOUBLE PRECISION SQRTS, THETA, MSQUARE
      CHARACTER*3 NUCS
      
C     Local variables
      DOUBLE COMPLEX M(2,2), M_CT(2,2), M_MS(2,2)
      DOUBLE PRECISION X
      
C     External functions
      DOUBLE PRECISION TRACE
      EXTERNAL TRACE
      
C     Calculate cosine of theta
      X = COS(THETA)
      
C     Always calculate the matrix fresh (don't reuse)
C     Get the raw matrix element
      CALL GET_RAW_M(SQRTS, THETA, NUCS, M)
      
C     Calculate the conjugate transpose
      CALL CONJUGATE_TRANSPOSE(M, M_CT)
      
C     Calculate M * M^\dagger
      CALL MAT_MUL(M, M_CT, M_MS)
      
C     Take the trace, which gives |M|^2
      MSQUARE = TRACE(M_MS)
      
      RETURN
      END