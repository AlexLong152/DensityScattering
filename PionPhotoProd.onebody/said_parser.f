C-----------------------------------------------------------------------
C     SAID Data Parser - Access functions for pion photoproduction multipoles
C
C     This file provides functions to parse and access SAID (Scattering Analysis 
C     Interactive Dial-in) multipole data for pion photoproduction. The data 
C     is organized by:
C
C     1. Multipole type (E or M for electric or magnetic)
C     2. Angular momentum (L) and total spin (J) configuration (plus or minus)
C     3. Target nucleon type (proton, neutron, or isospin-3/2 channel)
C
C     The main components are:
C     - PARSE_SPIN_STRING: Parse SAID-style multipole identifiers (e.g., "S11pE")
C     - GET_EPOLE_FUNC: Get electric multipole amplitudes for given kinematics
C     - GET_MPOLE_FUNC: Get magnetic multipole amplitudes for given kinematics
C     - PARSE_SAID_FILE: Parse and load the SAID data file into memory
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     Parse a SAID format spin string into its physical components
C
C     This subroutine parses SAID-format strings like 'S11pE' which encode:
C     - Angular momentum (S=0, P=1, D=2, F=3, G=4)
C     - Isospin and total angular momentum (first and second digit)
C     - Target nucleon type (p=proton, n=neutron)
C     - Multipole type (E=electric, M=magnetic)
C
C     Examples:
C     - "S11pE": L=0 (S-wave), I=1/2, J=1/2, proton target, Electric multipole
C     - "P33nM": L=1 (P-wave), I=3/2, J=3/2, neutron target, Magnetic multipole
C
C     Inputs:
C       SPINSTRING - The spin string to parse (e.g., "S11pE")
C
C     Outputs:
C       PLUSMINUS - Whether J=L+1/2 ("plus") or J=L-1/2 ("minus")
C       ELL - Orbital angular momentum value (0,1,2,3,4)
C       ISOSPIN - Isospin value (0.5 or 1.5)
C       SUBCHAN - Two-character string with target and multipole type (e.g., "pE")
C       SUCCESS - Logical flag indicating successful parsing
C-----------------------------------------------------------------------
      SUBROUTINE PARSE_SPIN_STRING(SPINSTRING, PLUSMINUS, ELL, ISOSPIN, 
     &                          SUBCHAN, SUCCESS)
      IMPLICIT NONE
C     Include constants definitions
      INCLUDE 'parseConstants.def'
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
C     Helper function to map TARGET string to index
C
C     Inputs:
C       TARGET - 'p12' (proton), 'n12' (neutron), or '32q' (isospin 3/2)
C
C     Returns:
C       1 for 'p12', 2 for 'n12', 3 for '32q', or 0 for invalid
C-----------------------------------------------------------------------
      INTEGER FUNCTION GET_TARGET_IDX(TARGET)
      IMPLICIT NONE
      CHARACTER*(*) TARGET
      
      IF (TARGET .EQ. 'p12') THEN
        GET_TARGET_IDX = 1
      ELSE IF (TARGET .EQ. 'n12') THEN
        GET_TARGET_IDX = 2
      ELSE IF (TARGET .EQ. '32q') THEN
        GET_TARGET_IDX = 3
      ELSE
        GET_TARGET_IDX = 0
      ENDIF
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Helper function to map SIGN string to index
C
C     Inputs:
C       SIGN - 'plus' (J=L+1/2) or 'minus' (J=L-1/2)
C
C     Returns:
C       1 for 'plus', 2 for 'minus', or 0 for invalid
C-----------------------------------------------------------------------
      INTEGER FUNCTION GET_SIGN_IDX(SIGN)
      IMPLICIT NONE
      CHARACTER*(*) SIGN
      
      IF (SIGN .EQ. 'plus') THEN
        GET_SIGN_IDX = 1
      ELSE IF (SIGN .EQ. 'minus') THEN
        GET_SIGN_IDX = 2
      ELSE
        GET_SIGN_IDX = 0
      ENDIF
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Generic function to retrieve multipole amplitude at specified energy
C
C     This function looks up either electric or magnetic multipole amplitudes
C     for a given partial wave at a specific center-of-mass energy. It 
C     automatically finds the closest available energy point in the data.
C
C     Inputs:
C       TYPE_CHAR - 'E' for electric or 'M' for magnetic multipole
C       SIGN_STR  - 'plus' (J=L+1/2) or 'minus' (J=L-1/2)
C       TARGET_STR - 'p12' (proton), 'n12' (neutron), or '32q' (isospin 3/2)
C       ELL       - Angular momentum (0=S, 1=P, 2=D, 3=F, 4=G)
C       SQRTS_VAL - Center-of-mass energy in MeV
C
C     Returns:
C       Complex amplitude value at the closest available energy.
C       If no data exists for the specified combination, returns 0.
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION GET_POLE_FUNC(TYPE_CHAR, SIGN_STR, 
     &                                      TARGET_STR, ELL, SQRTS_VAL)
      IMPLICIT NONE
      CHARACTER*1 TYPE_CHAR
      CHARACTER*(*) SIGN_STR
      CHARACTER*3 TARGET_STR 
      INTEGER ELL
      DOUBLE PRECISION SQRTS_VAL
      
C     Include constants definitions
      INCLUDE 'parseConstants.def'
      
C     External functions
      INTEGER GET_TARGET_IDX, GET_SIGN_IDX
      EXTERNAL GET_TARGET_IDX, GET_SIGN_IDX
      
C     Common blocks for pole data
      INTEGER EPOLE_POINTS(0:MAX_ELL, 3, 2)
      DOUBLE PRECISION EPOLE_SQRTS(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      DOUBLE COMPLEX EPOLE_AMPL(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      
      INTEGER MPOLE_POINTS(0:MAX_ELL, 3, 2)
      DOUBLE PRECISION MPOLE_SQRTS(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      DOUBLE COMPLEX MPOLE_AMPL(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      
      COMMON /EPOLE_DATA/ EPOLE_POINTS, EPOLE_SQRTS, EPOLE_AMPL
      COMMON /MPOLE_DATA/ MPOLE_POINTS, MPOLE_SQRTS, MPOLE_AMPL
      
C     Local variables
      INTEGER POINTS, I, IDX_MIN, TARGET_IDX, SIGN_IDX
      DOUBLE PRECISION MIN_DIFF, DIFF
      
C     Initialize default return value
      GET_POLE_FUNC = (0.0D0, 0.0D0)
      
C     Check if ELL is in valid range
      IF (ELL .LT. 0 .OR. ELL .GT. MAX_ELL) THEN
        RETURN
      ENDIF
      
C     Get target and sign indices
      TARGET_IDX = GET_TARGET_IDX(TARGET_STR)
      SIGN_IDX = GET_SIGN_IDX(SIGN_STR)
      
C     Check if indices are valid
      IF (TARGET_IDX .EQ. 0 .OR. SIGN_IDX .EQ. 0) THEN
        RETURN
      ENDIF
      
C     Get number of points based on TYPE
      IF (TYPE_CHAR .EQ. 'E' .OR. TYPE_CHAR .EQ. 'e') THEN
        POINTS = EPOLE_POINTS(ELL, TARGET_IDX, SIGN_IDX)
        
        IF (POINTS .GT. 0) THEN
C         Find the closest sqrtS value
          MIN_DIFF = ABS(EPOLE_SQRTS(1, ELL, TARGET_IDX, SIGN_IDX) 
     &                   - SQRTS_VAL)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(EPOLE_SQRTS(I, ELL, TARGET_IDX, SIGN_IDX) 
     &                - SQRTS_VAL)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
C         Return the amplitude at the closest energy point
          GET_POLE_FUNC = EPOLE_AMPL(IDX_MIN, ELL, TARGET_IDX, SIGN_IDX)
        ENDIF
      ELSE IF (TYPE_CHAR .EQ. 'M' .OR. TYPE_CHAR .EQ. 'm') THEN
        POINTS = MPOLE_POINTS(ELL, TARGET_IDX, SIGN_IDX)
        
        IF (POINTS .GT. 0) THEN
C         Find the closest sqrtS value
          MIN_DIFF = ABS(MPOLE_SQRTS(1, ELL, TARGET_IDX, SIGN_IDX) 
     &                   - SQRTS_VAL)
          IDX_MIN = 1
          
          DO I = 2, POINTS
            DIFF = ABS(MPOLE_SQRTS(I, ELL, TARGET_IDX, SIGN_IDX) 
     &                - SQRTS_VAL)
            IF (DIFF .LT. MIN_DIFF) THEN
              MIN_DIFF = DIFF
              IDX_MIN = I
            ENDIF
          ENDDO
          
C         Return the amplitude at the closest energy point
          GET_POLE_FUNC = MPOLE_AMPL(IDX_MIN, ELL, TARGET_IDX, SIGN_IDX)
        ENDIF
      ENDIF
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Function to retrieve Electric multipole amplitude at specified energy
C
C     This is a wrapper around the generic GET_POLE_FUNC for electric multipoles
C
C     Inputs:
C       SIGN_STR  - 'plus' (J=L+1/2) or 'minus' (J=L-1/2)
C       TARGET_STR - 'p12' (proton), 'n12' (neutron), or '32q' (isospin 3/2)
C       ELL       - Angular momentum (0=S, 1=P, 2=D, 3=F, 4=G)
C       SQRTS_VAL - Center-of-mass energy in MeV
C
C     Returns:
C       Complex amplitude value at the closest available energy.
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION GET_EPOLE_FUNC(SIGN_STR, TARGET_STR, 
     &                                      ELL, SQRTS_VAL)
      IMPLICIT NONE
      CHARACTER*(*) SIGN_STR
      CHARACTER*3 TARGET_STR
      INTEGER ELL
      DOUBLE PRECISION SQRTS_VAL
      
C     External function
      DOUBLE COMPLEX GET_POLE_FUNC
      EXTERNAL GET_POLE_FUNC
      
      GET_EPOLE_FUNC = GET_POLE_FUNC('E', SIGN_STR, TARGET_STR, 
     &                              ELL, SQRTS_VAL)
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Function to retrieve Magnetic multipole amplitude at specified energy
C
C     This is a wrapper around the generic GET_POLE_FUNC for magnetic multipoles
C
C     Inputs:
C       SIGN_STR  - 'plus' (J=L+1/2) or 'minus' (J=L-1/2)
C       TARGET_STR - 'p12' (proton), 'n12' (neutron), or '32q' (isospin 3/2)
C       ELL       - Angular momentum (0=S, 1=P, 2=D, 3=F, 4=G)
C       SQRTS_VAL - Center-of-mass energy in MeV
C
C     Returns:
C       Complex amplitude value at the closest available energy.
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION GET_MPOLE_FUNC(SIGN_STR, TARGET_STR, 
     &                                      ELL, SQRTS_VAL)
      IMPLICIT NONE
      CHARACTER*(*) SIGN_STR
      CHARACTER*3 TARGET_STR
      INTEGER ELL
      DOUBLE PRECISION SQRTS_VAL
      
C     External function
      DOUBLE COMPLEX GET_POLE_FUNC
      EXTERNAL GET_POLE_FUNC
      
      GET_MPOLE_FUNC = GET_POLE_FUNC('M', SIGN_STR, TARGET_STR, 
     &                              ELL, SQRTS_VAL)
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Parser for SAID multipole data files
C
C     This subroutine reads and parses a SAID data file containing pion 
C     photoproduction multipole amplitudes. It extracts the amplitude values
C     at various energy points for all multipoles and organizes them in the
C     data structures for later lookup.
C
C     The SAID data file format contains sections for each multipole, with
C     headers like "PI0P S11 pE" followed by tables of energy (EG), real part,
C     and imaginary part values.
C
C     The parsing process:
C     1. Identifies multipole headers (e.g., "PI0P S11 pE")
C     2. Parses the header to determine multipole characteristics
C     3. Reads the data lines that follow
C     4. Converts the data to proper units and stores in arrays
C
C     Inputs:
C       FILENAME     - Path to the SAID data file (e.g., "said-SM22.txt")
C       UNITS_FACTOR - Scale factor to convert amplitudes to appropriate units
C                      (typically needed to convert from SAID units to MeV^-1)
C
C     Output:
C       The data is stored in the EPOLE_* and MPOLE_* arrays via common blocks,
C       making it available to the GET_EPOLE_FUNC and GET_MPOLE_FUNC functions.
C-----------------------------------------------------------------------
      SUBROUTINE PARSE_SAID_FILE(FILENAME, UNITS_FACTOR)
      IMPLICIT NONE
      CHARACTER*(*) FILENAME
      DOUBLE PRECISION UNITS_FACTOR
      
C     Include constants definitions
      INCLUDE 'parseConstants.def'
      
C     External functions
      DOUBLE PRECISION LAB_E_SQRTS
      EXTERNAL LAB_E_SQRTS
      
C     Common blocks for physics constants
      DOUBLE PRECISION MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      COMMON /PHYS_CONST/ MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      
C     Local variables
      CHARACTER*(MAX_LINE_LEN) LINE
      INTEGER IOS, LINE_NUM, I, J, K
      
C     External function declarations
      INTEGER GET_TARGET_IDX, GET_SIGN_IDX
      EXTERNAL GET_TARGET_IDX, GET_SIGN_IDX

C     Common blocks for pole data with compact, indexed arrays
      INTEGER EPOLE_POINTS(0:MAX_ELL, 3, 2)
      DOUBLE PRECISION EPOLE_SQRTS(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      DOUBLE COMPLEX EPOLE_AMPL(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      
      INTEGER MPOLE_POINTS(0:MAX_ELL, 3, 2)
      DOUBLE PRECISION MPOLE_SQRTS(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      DOUBLE COMPLEX MPOLE_AMPL(MAX_ENERGY_POINTS, 0:MAX_ELL, 3, 2)
      
      COMMON /EPOLE_DATA/ EPOLE_POINTS, EPOLE_SQRTS, EPOLE_AMPL
      COMMON /MPOLE_DATA/ MPOLE_POINTS, MPOLE_SQRTS, MPOLE_AMPL
      
C     File parsing variables - consolidated to reduce variable count
      CHARACTER*10 WAVE1, WAVE2, SUBCH
      CHARACTER*5 SPINSTRING, PLUSMINUS
      CHARACTER*3 TARGET
      CHARACTER*2 SUBCHAN
      CHARACTER*1 AMP_TYPE
      DOUBLE PRECISION ISOSPIN, LAB_E, RE_VAL, IM_VAL, SQRTS_VAL
      DOUBLE COMPLEX AMP_VAL
      INTEGER ELL, TARGET_IDX, SIGN_IDX, POINTS_COUNT
      LOGICAL IS_DATA_SECTION, PARSE_SUCCESS
      
C     Open the file
      OPEN(UNIT=10, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
        WRITE(*,*) 'Error opening file: ', FILENAME
        RETURN
      ENDIF
      
C     Initialize state
      LINE_NUM = 0
      IS_DATA_SECTION = .FALSE.
      ELL = -1
      TARGET_IDX = 0
      SIGN_IDX = 0
      
C     Read the file line by line
      DO WHILE (.TRUE.)
        READ(10, '(A)', IOSTAT=IOS) LINE
        IF (IOS .NE. 0) EXIT
        
        LINE_NUM = LINE_NUM + 1
        
        IF (INDEX(LINE, 'PI0P ') .GT. 0 .OR. 
     &      INDEX(LINE, 'PI0N ') .GT. 0) THEN
C         This is a header line - extract wave info
          READ(LINE, *, IOSTAT=IOS) WAVE1, WAVE2, SUBCH
          
          IF (IOS .EQ. 0) THEN
C           Create and parse spinstring
            SPINSTRING = TRIM(WAVE2) // TRIM(SUBCH)
            CALL PARSE_SPIN_STRING(SPINSTRING, PLUSMINUS, ELL, ISOSPIN,
     &                           SUBCHAN, PARSE_SUCCESS)
            
            IF (PARSE_SUCCESS) THEN
C             Determine target from isospin and subchannel
              IF (ABS(ISOSPIN - 1.5D0) .LT. 1.0D-9) THEN
                TARGET = '32q'
              ELSE
                IF (SUBCHAN(1:1) .EQ. 'p') THEN
                  TARGET = 'p12'
                ELSE
                  TARGET = 'n12'
                ENDIF
              ENDIF
              
              TARGET_IDX = GET_TARGET_IDX(TARGET)
              SIGN_IDX = GET_SIGN_IDX(PLUSMINUS)
              AMP_TYPE = SUBCHAN(2:2)
              IS_DATA_SECTION = .TRUE.
            ELSE
C             Invalid spinstring, reset state
              ELL = -1
              TARGET_IDX = 0
              SIGN_IDX = 0
              IS_DATA_SECTION = .FALSE.
            ENDIF
          ENDIF
        ELSE IF (IS_DATA_SECTION .AND. TARGET_IDX .GT. 0 
     &          .AND. SIGN_IDX .GT. 0 .AND. ELL .GE. 0) THEN
C         Check if this is a header line to skip
          IF (INDEX(LINE, 'EG') .GT. 0 .OR. 
     &        INDEX(LINE, 'EMreal') .GT. 0 .OR.
     &        INDEX(LINE, 'GWU') .GT. 0) THEN
C           Skip column headers
            CONTINUE
          ELSE
C           Try to parse as data line
            READ(LINE, *, IOSTAT=IOS) LAB_E, RE_VAL, IM_VAL
            
            IF (IOS .EQ. 0) THEN
C             Process valid data line
              SQRTS_VAL = LAB_E_SQRTS(LAB_E)
              AMP_VAL = DCMPLX(RE_VAL, IM_VAL) * UNITS_FACTOR
              
C             Store in appropriate array based on amplitude type
              IF (AMP_TYPE .EQ. 'E' .OR. AMP_TYPE .EQ. 'e') THEN
                POINTS_COUNT = EPOLE_POINTS(ELL, TARGET_IDX, SIGN_IDX)+1
                EPOLE_POINTS(ELL, TARGET_IDX, SIGN_IDX) = POINTS_COUNT
                EPOLE_SQRTS(POINTS_COUNT, ELL, TARGET_IDX, SIGN_IDX) = 
     &              SQRTS_VAL
                EPOLE_AMPL(POINTS_COUNT, ELL, TARGET_IDX, SIGN_IDX) = 
     &              AMP_VAL
              ELSE IF (AMP_TYPE .EQ. 'M' .OR. AMP_TYPE .EQ. 'm') THEN
                POINTS_COUNT = MPOLE_POINTS(ELL, TARGET_IDX, SIGN_IDX)+1
                MPOLE_POINTS(ELL, TARGET_IDX, SIGN_IDX) = POINTS_COUNT
                MPOLE_SQRTS(POINTS_COUNT, ELL, TARGET_IDX, SIGN_IDX) = 
     &              SQRTS_VAL
                MPOLE_AMPL(POINTS_COUNT, ELL, TARGET_IDX, SIGN_IDX) = 
     &              AMP_VAL
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      
C     Close the file
      CLOSE(10)
      
      RETURN
      END