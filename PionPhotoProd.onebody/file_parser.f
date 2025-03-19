C     No dummy program needed

C-----------------------------------------------------------------------
C     File parsing module for SAID multipole data
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     Parse a spin string like 'S11pE' into its components
C-----------------------------------------------------------------------
      SUBROUTINE PARSE_SPIN_STRING(SPINSTRING, PLUSMINUS, ELL, ISOSPIN, 
     &                             SUBCHAN, SUCCESS)
      IMPLICIT NONE
      CHARACTER*(*) SPINSTRING
      CHARACTER*5 PLUSMINUS
      INTEGER ELL
      DOUBLE PRECISION ISOSPIN
      CHARACTER*2 SUBCHAN
      LOGICAL SUCCESS
      
C     Local variables
      CHARACTER*3 PW
      CHARACTER*1 LETTER
      INTEGER TWOI, TWOJ
      DOUBLE PRECISION DIFF, J
      LOGICAL IS_BAD_STRING
      
      
C     Check if string is one of the known bad strings
      IS_BAD_STRING = .FALSE.
      IF (SPINSTRING .EQ. "S11pM" .OR. SPINSTRING .EQ. "S11nM" .OR.
     &    SPINSTRING .EQ. "S31nM" .OR. SPINSTRING .EQ. "S31pM" .OR.
     &    SPINSTRING .EQ. "P11pE" .OR. SPINSTRING .EQ. "P11nE" .OR.
     &    SPINSTRING .EQ. "P31pE" .OR. SPINSTRING .EQ. "P31nE") THEN
        WRITE(*,*) "Warning: Unphysical spin string: ", SPINSTRING
        IS_BAD_STRING = .TRUE.
      ENDIF
      
C     Initialize output
      PLUSMINUS = 'none'
      ELL = -1
      ISOSPIN = -1.0D0
      SUBCHAN = '??'
      SUCCESS = .FALSE.
      
C     Check string length
      IF (LEN_TRIM(SPINSTRING) .LT. 5) THEN
        WRITE(*,*) 'Error: String too short: ', SPINSTRING
        RETURN
      ENDIF
      
C     Extract components: first 3 chars are wave (S11, P33, etc.)
C     and remaining 2 chars are subchannel (pE, nM, etc.)
      PW = SPINSTRING(1:3)
      SUBCHAN = SPINSTRING(4:5)
      
      
C     Extract first letter (S, P, D, F, G, ...)
      LETTER = PW(1:1)
      
C     Extract twoI and twoJ
      READ(PW(2:2), '(I1)', ERR=100) TWOI
      READ(PW(3:3), '(I1)', ERR=100) TWOJ
      
C     Map letter to angular momentum L
      IF (LETTER .EQ. 'S' .OR. LETTER .EQ. 's') THEN
        ELL = 0
      ELSE IF (LETTER .EQ. 'P' .OR. LETTER .EQ. 'p') THEN
        ELL = 1
      ELSE IF (LETTER .EQ. 'D' .OR. LETTER .EQ. 'd') THEN
        ELL = 2
      ELSE IF (LETTER .EQ. 'F' .OR. LETTER .EQ. 'f') THEN
        ELL = 3
      ELSE IF (LETTER .EQ. 'G' .OR. LETTER .EQ. 'g') THEN
        ELL = 4
      ELSE
        WRITE(*,*) 'Error: Invalid letter: ', LETTER
        GOTO 100
      ENDIF
      
C     Calculate isospin and total angular momentum
      ISOSPIN = 0.5D0 * DBLE(TWOI)
      J = 0.5D0 * DBLE(TWOJ)
      
C     Determine plus/minus from J-L
      DIFF = J - DBLE(ELL)
      IF (ABS(DIFF - 0.5D0) .LT. 1.0D-9) THEN
        PLUSMINUS = 'plus'
      ELSE IF (ABS(DIFF + 0.5D0) .LT. 1.0D-9) THEN
        PLUSMINUS = 'minus'
      ELSE
        WRITE(*,*) 'Error: Invalid J-L difference: ', DIFF
        GOTO 100
      ENDIF
      
      SUCCESS = .TRUE.
      RETURN
      
100   CONTINUE
C     Error return point
      WRITE(*,*) 'Error parsing spin string'
      SUCCESS = .FALSE.
      RETURN
      END
      
C-----------------------------------------------------------------------
C     Simpler manual parser for the SAID data file format
C-----------------------------------------------------------------------
      SUBROUTINE PARSE_SAID_FILE(FILENAME, UNITS_FACTOR)
      IMPLICIT NONE
      CHARACTER*(*) FILENAME
      DOUBLE PRECISION UNITS_FACTOR
      
C     External function declaration
      DOUBLE PRECISION LAB_E_SQRTS
      EXTERNAL LAB_E_SQRTS
      
C     Common blocks for physics constants
      DOUBLE PRECISION MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      COMMON /PHYS_CONST/ MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      
C     Local variables
      INTEGER MAX_LINE_LEN
      PARAMETER (MAX_LINE_LEN = 256)
      CHARACTER*(MAX_LINE_LEN) LINE
      INTEGER IOS, LINE_NUM, I, J, K
      
C     Current state variables
      CHARACTER*10 WAVE1, WAVE2, SUBCH
      CHARACTER*5 SPINSTRING, PLUSMINUS
      INTEGER ELL
      DOUBLE PRECISION ISOSPIN
      CHARACTER*3 TARGET
      CHARACTER*1 AMP_TYPE
      LOGICAL IS_DATA_SECTION
      LOGICAL PARSE_SUCCESS
      CHARACTER*2 SUBCHAN
      
C     Data arrays
      INTEGER MAX_ENG, MAX_L
      PARAMETER (MAX_ENG = 200, MAX_L = 5)
      
C     Data dimensions for E-poles
      INTEGER EPOLE_POINTS_PLUS_P12(0:MAX_L)
      INTEGER EPOLE_POINTS_PLUS_N12(0:MAX_L)
      INTEGER EPOLE_POINTS_PLUS_32Q(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_P12(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_N12(0:MAX_L)
      INTEGER EPOLE_POINTS_MINUS_32Q(0:MAX_L)
      
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
      
C     Data dimensions for M-poles
      INTEGER MPOLE_POINTS_PLUS_P12(0:MAX_L)
      INTEGER MPOLE_POINTS_PLUS_N12(0:MAX_L)
      INTEGER MPOLE_POINTS_PLUS_32Q(0:MAX_L)
      INTEGER MPOLE_POINTS_MINUS_P12(0:MAX_L)
      INTEGER MPOLE_POINTS_MINUS_N12(0:MAX_L)
      INTEGER MPOLE_POINTS_MINUS_32Q(0:MAX_L)
      
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
      
C     Data line parsing variables
      DOUBLE PRECISION LAB_E, RE_VAL, IM_VAL
      DOUBLE PRECISION SQRTS_VAL
      DOUBLE COMPLEX AMP_VAL
      
C     Open the file
      OPEN(UNIT=10, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
        WRITE(*,*) 'Error opening file: ', FILENAME
        RETURN
      ENDIF
      
C     Initialize state
      LINE_NUM = 0
      PLUSMINUS = 'none'
      ELL = -1
      TARGET = '   '
      AMP_TYPE = ' '
      IS_DATA_SECTION = .FALSE.
      
      
C     Read the file line by line
      DO WHILE (.TRUE.)
        READ(10, '(A)', IOSTAT=IOS) LINE
        IF (IOS .NE. 0) EXIT
        
        LINE_NUM = LINE_NUM + 1
        
        IF (INDEX(LINE, 'PI0P ') .GT. 0 .OR. 
     &      INDEX(LINE, 'PI0N ') .GT. 0) THEN
C         This is a header line - try to extract wave info
          
C         Read up to 3 space-separated values
          READ(LINE, *, IOSTAT=IOS) WAVE1, WAVE2, SUBCH
          
          IF (IOS .EQ. 0) THEN
            
C           Create spinstring
            SPINSTRING = TRIM(WAVE2) // TRIM(SUBCH)
            
C           Parse the spinstring
            CALL PARSE_SPIN_STRING(SPINSTRING, PLUSMINUS, ELL, ISOSPIN,
     &                           SUBCHAN, PARSE_SUCCESS)
            
            IF (PARSE_SUCCESS) THEN
              
C             Determine target from isospin
              IF (ABS(ISOSPIN - 1.5D0) .LT. 1.0D-9) THEN
                TARGET = '32q'
              ELSE
                IF (SUBCHAN(1:1) .EQ. 'p') THEN
                  TARGET = 'p12'
                ELSE
                  TARGET = 'n12'
                ENDIF
              ENDIF
              
C             Determine amplitude type
              AMP_TYPE = SUBCHAN(2:2)
              
C             Set flag for data section
              IS_DATA_SECTION = .TRUE.
            ELSE
              WRITE(*,*) 'Error: Failed to parse spinstring'
              PLUSMINUS = 'none'
              ELL = -1
              TARGET = '   '
              AMP_TYPE = ' '
              IS_DATA_SECTION = .FALSE.
            ENDIF
          ENDIF
        ELSE IF (IS_DATA_SECTION) THEN
C         Check if this is a header line to skip
          IF (INDEX(LINE, 'EG') .GT. 0 .OR. 
     &        INDEX(LINE, 'EMreal') .GT. 0 .OR.
     &        INDEX(LINE, 'GWU') .GT. 0) THEN
C           Skip column headers
          ELSE
C           Try to parse as data line
            READ(LINE, *, IOSTAT=IOS) LAB_E, RE_VAL, IM_VAL
            
            IF (IOS .EQ. 0) THEN
              
C             Convert lab energy to sqrt(s)
              SQRTS_VAL = LAB_E_SQRTS(LAB_E)
              
C             Convert amplitude to MeV^-1
              AMP_VAL = DCMPLX(RE_VAL, IM_VAL) * UNITS_FACTOR
              
C             Store in appropriate array
              IF (AMP_TYPE .EQ. 'E' .OR. AMP_TYPE .EQ. 'e') THEN
                IF (PLUSMINUS .EQ. 'plus') THEN
                  IF (TARGET .EQ. 'p12') THEN
                    EPOLE_POINTS_PLUS_P12(ELL) = 
     &                EPOLE_POINTS_PLUS_P12(ELL) + 1
                    K = EPOLE_POINTS_PLUS_P12(ELL)
                    EPOLE_SQRTS_PLUS_P12(K, ELL) = SQRTS_VAL
                    EPOLE_AMPL_PLUS_P12(K, ELL) = AMP_VAL
                  ELSE IF (TARGET .EQ. 'n12') THEN
                    EPOLE_POINTS_PLUS_N12(ELL) = 
     &                EPOLE_POINTS_PLUS_N12(ELL) + 1
                    K = EPOLE_POINTS_PLUS_N12(ELL)
                    EPOLE_SQRTS_PLUS_N12(K, ELL) = SQRTS_VAL
                    EPOLE_AMPL_PLUS_N12(K, ELL) = AMP_VAL
                  ELSE IF (TARGET .EQ. '32q') THEN
                    EPOLE_POINTS_PLUS_32Q(ELL) = 
     &                EPOLE_POINTS_PLUS_32Q(ELL) + 1
                    K = EPOLE_POINTS_PLUS_32Q(ELL)
                    EPOLE_SQRTS_PLUS_32Q(K, ELL) = SQRTS_VAL
                    EPOLE_AMPL_PLUS_32Q(K, ELL) = AMP_VAL
                  ENDIF
                ELSE IF (PLUSMINUS .EQ. 'minus') THEN
                  IF (TARGET .EQ. 'p12') THEN
                    EPOLE_POINTS_MINUS_P12(ELL) = 
     &                EPOLE_POINTS_MINUS_P12(ELL) + 1
                    K = EPOLE_POINTS_MINUS_P12(ELL)
                    EPOLE_SQRTS_MINUS_P12(K, ELL) = SQRTS_VAL
                    EPOLE_AMPL_MINUS_P12(K, ELL) = AMP_VAL
                  ELSE IF (TARGET .EQ. 'n12') THEN
                    EPOLE_POINTS_MINUS_N12(ELL) = 
     &                EPOLE_POINTS_MINUS_N12(ELL) + 1
                    K = EPOLE_POINTS_MINUS_N12(ELL)
                    EPOLE_SQRTS_MINUS_N12(K, ELL) = SQRTS_VAL
                    EPOLE_AMPL_MINUS_N12(K, ELL) = AMP_VAL
                  ELSE IF (TARGET .EQ. '32q') THEN
                    EPOLE_POINTS_MINUS_32Q(ELL) = 
     &                EPOLE_POINTS_MINUS_32Q(ELL) + 1
                    K = EPOLE_POINTS_MINUS_32Q(ELL)
                    EPOLE_SQRTS_MINUS_32Q(K, ELL) = SQRTS_VAL
                    EPOLE_AMPL_MINUS_32Q(K, ELL) = AMP_VAL
                  ENDIF
                ENDIF
              ELSE IF (AMP_TYPE .EQ. 'M' .OR. AMP_TYPE .EQ. 'm') THEN
                IF (PLUSMINUS .EQ. 'plus') THEN
                  IF (TARGET .EQ. 'p12') THEN
                    MPOLE_POINTS_PLUS_P12(ELL) = 
     &                MPOLE_POINTS_PLUS_P12(ELL) + 1
                    K = MPOLE_POINTS_PLUS_P12(ELL)
                    MPOLE_SQRTS_PLUS_P12(K, ELL) = SQRTS_VAL
                    MPOLE_AMPL_PLUS_P12(K, ELL) = AMP_VAL
                  ELSE IF (TARGET .EQ. 'n12') THEN
                    MPOLE_POINTS_PLUS_N12(ELL) = 
     &                MPOLE_POINTS_PLUS_N12(ELL) + 1
                    K = MPOLE_POINTS_PLUS_N12(ELL)
                    MPOLE_SQRTS_PLUS_N12(K, ELL) = SQRTS_VAL
                    MPOLE_AMPL_PLUS_N12(K, ELL) = AMP_VAL
                  ELSE IF (TARGET .EQ. '32q') THEN
                    MPOLE_POINTS_PLUS_32Q(ELL) = 
     &                MPOLE_POINTS_PLUS_32Q(ELL) + 1
                    K = MPOLE_POINTS_PLUS_32Q(ELL)
                    MPOLE_SQRTS_PLUS_32Q(K, ELL) = SQRTS_VAL
                    MPOLE_AMPL_PLUS_32Q(K, ELL) = AMP_VAL
                  ENDIF
                ELSE IF (PLUSMINUS .EQ. 'minus') THEN
                  IF (TARGET .EQ. 'p12') THEN
                    MPOLE_POINTS_MINUS_P12(ELL) = 
     &                MPOLE_POINTS_MINUS_P12(ELL) + 1
                    K = MPOLE_POINTS_MINUS_P12(ELL)
                    MPOLE_SQRTS_MINUS_P12(K, ELL) = SQRTS_VAL
                    MPOLE_AMPL_MINUS_P12(K, ELL) = AMP_VAL
                  ELSE IF (TARGET .EQ. 'n12') THEN
                    MPOLE_POINTS_MINUS_N12(ELL) = 
     &                MPOLE_POINTS_MINUS_N12(ELL) + 1
                    K = MPOLE_POINTS_MINUS_N12(ELL)
                    MPOLE_SQRTS_MINUS_N12(K, ELL) = SQRTS_VAL
                    MPOLE_AMPL_MINUS_N12(K, ELL) = AMP_VAL
                  ELSE IF (TARGET .EQ. '32q') THEN
                    MPOLE_POINTS_MINUS_32Q(ELL) = 
     &                MPOLE_POINTS_MINUS_32Q(ELL) + 1
                    K = MPOLE_POINTS_MINUS_32Q(ELL)
                    MPOLE_SQRTS_MINUS_32Q(K, ELL) = SQRTS_VAL
                    MPOLE_AMPL_MINUS_32Q(K, ELL) = AMP_VAL
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      
C     Close the file
      CLOSE(10)
      
      
      RETURN
      END