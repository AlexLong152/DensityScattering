C-----------------------------------------------------------------------
C     Physics functions for pion photoproduction
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     Convert lab-frame photon energy E_lab -> center-of-mass energy sqrt(s).
C     sqrt(s) = sqrt(mN^2 + 2*mN*E_lab),
C     assuming the target nucleon is at rest (proton or neutron ~ same mass).
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION LAB_E_SQRTS(E_LAB)
      IMPLICIT NONE
      DOUBLE PRECISION E_LAB
      
C     Common block for physical constants
      DOUBLE PRECISION MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      COMMON /PHYS_CONST/ MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      
      LAB_E_SQRTS = SQRT(MN*MN + 2.0D0 * MN * E_LAB)
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Calculate kinematics for pion photoproduction
C
C     This subroutine sets up the kinematics for pion photoproduction in
C     the center-of-mass (CM) frame. It calculates the photon and pion
C     momenta based on the CM energy and scattering angle.
C
C     The calculation assumes the photon is moving along the z-axis,
C     and the pion scattering occurs in the y-z plane.
C
C     Inputs:
C       SQRTS - Center-of-mass energy in MeV
C       X     - Cosine of the pion scattering angle in CM frame
C       NUCS  - Process code indicating the nucleon and pion types:
C               'pp0': γ + p → p + π⁰
C               'nn0': γ + n → n + π⁰
C               'pn+': γ + p → n + π⁺
C               'np-': γ + n → p + π⁻
C
C     Outputs:
C       S     - Mandelstam variable s = SQRTS² (center-of-mass energy squared)
C       KVEC  - 3D photon momentum vector
C       QVEC  - 3D pion momentum vector
C-----------------------------------------------------------------------
      SUBROUTINE GET_KINEMATICS(SQRTS, X, NUCS, S, KVEC, QVEC)
      IMPLICIT NONE
      DOUBLE PRECISION SQRTS, X, S
      CHARACTER*3 NUCS
      DOUBLE PRECISION KVEC(3), QVEC(3)
      
C     Common block for physical constants
      DOUBLE PRECISION MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      COMMON /PHYS_CONST/ MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      
C     Local variables
      DOUBLE PRECISION MNUCL, MPION   ! Masses for the specific process
      DOUBLE PRECISION OMEGA          ! Photon energy in CM frame
      DOUBLE PRECISION EPI            ! Pion energy in CM frame
      DOUBLE PRECISION ABSQ           ! Magnitude of pion momentum
      DOUBLE PRECISION SINTH          ! Sine of scattering angle
      
C     1. Select appropriate masses based on process code
      CALL SELECT_MASSES(NUCS, MNUCL, MPION)
      
C     2. Calculate center-of-mass quantities
      S = SQRTS * SQRTS
      
C     3. Calculate photon energy in CM frame (E_γ = (s - m_N²)/(2√s))
      OMEGA = (S - MNUCL*MNUCL) / (2.0D0 * SQRTS)
      
C     4. Set photon momentum vector (along z-axis)
      KVEC(1) = 0.0D0
      KVEC(2) = 0.0D0
      KVEC(3) = OMEGA
      
C     5. Calculate pion energy in CM frame (E_π = (s + m_π² - m_N²)/(2√s))
      EPI = (S + MPION*MPION - MNUCL*MNUCL) / (2.0D0 * SQRTS)
      
C     6. Calculate magnitude of pion momentum (|q| = √(E_π² - m_π²))
      ABSQ = SQRT(EPI*EPI - MPION*MPION)
      
C     7. Calculate sine of scattering angle (sin(theta) = sqrt(1-x^2))
      SINTH = SQRT(1.0D0 - X*X)
      
C     8. Set pion momentum vector
C        x-component = 0 (in y-z plane)
C        y-component = |q|*sin(theta)
C        z-component = |q|*cos(theta)
      QVEC(1) = 0.0D0
      QVEC(2) = SINTH * ABSQ
      QVEC(3) = X * ABSQ
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Helper function to select appropriate masses for the process
C-----------------------------------------------------------------------
      SUBROUTINE SELECT_MASSES(NUCS, MNUCL, MPION)
      IMPLICIT NONE
      CHARACTER*3 NUCS
      DOUBLE PRECISION MNUCL, MPION
      
C     Common block for physical constants
      DOUBLE PRECISION MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      COMMON /PHYS_CONST/ MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      
C     Select masses based on process code
      IF (NUCS .EQ. 'pp0') THEN
C       γ + p → p + π⁰
        MNUCL = MPROTON
        MPION = MPI
      ELSE IF (NUCS .EQ. 'nn0') THEN
C       γ + n → n + π⁰
        MNUCL = MNEUTRON
        MPION = MPI
      ELSE IF (NUCS .EQ. 'pn+') THEN
C       γ + p → n + π⁺
        MNUCL = MPROTON
        MPION = MPIPLUS
      ELSE IF (NUCS .EQ. 'np-') THEN
C       γ + n → p + π⁻
        MNUCL = MNEUTRON
        MPION = MPIPLUS
      ELSE
C       Default to average nucleon and neutral pion
        MNUCL = MN
        MPION = MPI
      ENDIF
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Calculate full form factor matrix for pion photoproduction
C
C     This subroutine calculates the 2×2 matrix form factor F needed for
C     pion photoproduction amplitudes. The form factor incorporates four
C     distinct terms (F1-F4) and combines them with kinematic factors.
C
C     The full formula is:
C     F = i*σ·ε*F1 + (σ·q)(σ·(k×ε))*F2/(|q|*|k|) + 
C         i*(σ·k)*(q·ε)*F3/(|q|*|k|) + i*(σ·q)*(q·ε)*F4/(|q|*|k|)
C
C     Inputs:
C       X      - cos(θ) where θ is the pion scattering angle in CM frame
C       SQRTS  - Center-of-mass energy in MeV
C       QVEC   - Pion momentum vector (3D)
C       KVEC   - Photon momentum vector (3D)
C       EPSVEC - Photon polarization vector (3D)
C       TARGET - Target nucleon specification: 'p12', 'n12', or '32q'
C
C     Output:
C       RESULT - 2×2 complex matrix representing the form factor
C-----------------------------------------------------------------------
      SUBROUTINE F_FUNC(X, SQRTS, QVEC, KVEC, EPSVEC, TARGET, RESULT)
      IMPLICIT NONE
      DOUBLE PRECISION X, SQRTS
      DOUBLE PRECISION QVEC(3), KVEC(3), EPSVEC(3)
      CHARACTER*3 TARGET
      DOUBLE COMPLEX RESULT(2,2)
      
C     External functions
      DOUBLE PRECISION VEC_ABS, DOT_PRODUCT
      EXTERNAL VEC_ABS, DOT_PRODUCT
      
C     Local variables
      DOUBLE COMPLEX F1, F2, F3, F4
      DOUBLE COMPLEX TERM1(2,2), TERM2(2,2), TERM3(2,2), TERM4(2,2)
      DOUBLE COMPLEX TEMP1(2,2), TEMP2(2,2)
      DOUBLE COMPLEX I_UNIT, FACTOR
      DOUBLE PRECISION KVEC_ABS, QVEC_ABS, DOT_PROD, DENOMINATOR
      DOUBLE PRECISION CROSS_PROD(3)
      INTEGER I, J
      
C     Initialize result and temporary matrices to zero
      DO I = 1, 2
        DO J = 1, 2
          RESULT(I,J) = (0.0D0, 0.0D0)
          TERM1(I,J) = (0.0D0, 0.0D0)
          TERM2(I,J) = (0.0D0, 0.0D0)
          TERM3(I,J) = (0.0D0, 0.0D0)
          TERM4(I,J) = (0.0D0, 0.0D0)
          TEMP1(I,J) = (0.0D0, 0.0D0)
          TEMP2(I,J) = (0.0D0, 0.0D0)
        ENDDO
      ENDDO
      
C     Initialize imaginary unit
      I_UNIT = (0.0D0, 1.0D0)
      
C     Calculate kinematic factors
      KVEC_ABS = VEC_ABS(KVEC)
      QVEC_ABS = VEC_ABS(QVEC)
      DOT_PROD = DOT_PRODUCT(QVEC, EPSVEC)  ! q·ε
      DENOMINATOR = QVEC_ABS * KVEC_ABS
      
C     Calculate cross product k×ε
      CALL CROSS_PRODUCT(KVEC, EPSVEC, CROSS_PROD)
      
C     Get individual form factors F1-F4
      CALL GET_FORM_FACTOR(X, SQRTS, 'plus', TARGET, 1, F1)
      CALL GET_FORM_FACTOR(X, SQRTS, 'plus', TARGET, 2, F2)
      CALL GET_FORM_FACTOR(X, SQRTS, 'plus', TARGET, 3, F3)
      CALL GET_FORM_FACTOR(X, SQRTS, 'plus', TARGET, 4, F4)
      
C     Calculate F1 term: i*σ·ε*F1
      CALL MAT_DOT_VEC(EPSVEC, TERM1)
      FACTOR = I_UNIT * F1
      DO I = 1, 2
        DO J = 1, 2
          TERM1(I,J) = TERM1(I,J) * FACTOR
        ENDDO
      ENDDO
      
C     Calculate F2 term: (σ·q)(σ·(k×ε))*F2 / (|q|*|k|)
      CALL MAT_DOT_VEC(QVEC, TEMP1)
      CALL MAT_DOT_VEC(CROSS_PROD, TEMP2)
      CALL MAT_MUL(TEMP1, TEMP2, TERM2)
      FACTOR = F2 / DENOMINATOR
      DO I = 1, 2
        DO J = 1, 2
          TERM2(I,J) = TERM2(I,J) * FACTOR
        ENDDO
      ENDDO
      
C     Calculate F3 term: i*(σ·k)*(q·ε)*F3 / (|q|*|k|)
      CALL MAT_DOT_VEC(KVEC, TERM3)
      FACTOR = I_UNIT * DOT_PROD * F3 / DENOMINATOR
      DO I = 1, 2
        DO J = 1, 2
          TERM3(I,J) = TERM3(I,J) * FACTOR
        ENDDO
      ENDDO
      
C     Calculate F4 term: i*(σ·q)*(q·ε)*F4 / (|q|*|k|)
      CALL MAT_DOT_VEC(QVEC, TERM4)
      FACTOR = I_UNIT * DOT_PROD * F4 / DENOMINATOR
      DO I = 1, 2
        DO J = 1, 2
          TERM4(I,J) = TERM4(I,J) * FACTOR
        ENDDO
      ENDDO

C     Sum all terms
      DO I = 1, 2
        DO J = 1, 2
          RESULT(I,J) = TERM1(I,J) + TERM2(I,J) + 
     &                   TERM3(I,J) + TERM4(I,J)
        ENDDO
      ENDDO
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Calculate form factor F(i) for pion photoproduction
C
C     This subroutine calculates the form factors F1-F4 used in pion 
C     photoproduction amplitude calculations, expressed in terms of 
C     electric and magnetic multipole amplitudes.
C
C     Inputs:
C       X      - cos(θ) where θ is the pion scattering angle in CM frame
C       SQRTS  - Center-of-mass energy in MeV
C       SIGN   - 'plus' or 'minus' (J=L+1/2 or J=L-1/2)
C       TARGET - 'p12', 'n12', or '32q' (proton, neutron, or isospin 3/2)
C       FI     - Form factor index (1-4)
C
C     Output:
C       RESULT - Complex value of the requested form factor
C-----------------------------------------------------------------------
      SUBROUTINE GET_FORM_FACTOR(X, SQRTS, SIGN, TARGET, FI, RESULT)
      IMPLICIT NONE
      DOUBLE PRECISION X, SQRTS
      CHARACTER*(*) SIGN
      CHARACTER*3 TARGET
      INTEGER FI
      DOUBLE COMPLEX RESULT
      
C     External functions
      DOUBLE PRECISION LEGP
      EXTERNAL LEGP
      
C     Local variables
      DOUBLE COMPLEX EPLUS, MPLUS, EMINUS, MMINUS
      DOUBLE COMPLEX TMP
      INTEGER ELL, DERIV, N1, N2
      DOUBLE PRECISION LP1, LP2
      
C     Initialize result
      RESULT = (0.0D0, 0.0D0)
      
C     Loop over partial waves (typically L=0-3 is sufficient)
      DO ELL = 0, 3
C       Get the multipole amplitudes for both + and - cases
        CALL GET_EPOLE(SIGN, TARGET, ELL, SQRTS, EPLUS)
        CALL GET_MPOLE(SIGN, TARGET, ELL, SQRTS, MPLUS)
        CALL GET_EPOLE('minus', TARGET, ELL, SQRTS, EMINUS)
        CALL GET_MPOLE('minus', TARGET, ELL, SQRTS, MMINUS)
        
C       Calculate form factor based on FI
        SELECT CASE (FI)
          CASE (1)
C           F1 term: First derivative Legendre polynomials
            DERIV = 1
            N1 = ELL + 1
            N2 = ELL - 1
            LP1 = LEGP(X, N1, DERIV)
            LP2 = LEGP(X, N2, DERIV)
            
            TMP = (ELL * MPLUS + EPLUS) * LP1 + 
     &            ((ELL + 1) * MMINUS + EMINUS) * LP2
            
          CASE (2)
C           F2 term: First derivative Legendre polynomial
            DERIV = 1
            N1 = ELL
            LP1 = LEGP(X, N1, DERIV)
            
            TMP = ((ELL + 1) * MPLUS + (ELL * MMINUS)) * LP1
            
          CASE (3)
C           F3 term: Second derivative Legendre polynomials
            DERIV = 2
            N1 = ELL + 1
            N2 = ELL - 1
            LP1 = LEGP(X, N1, DERIV)
            LP2 = LEGP(X, N2, DERIV)
            
            TMP = (EPLUS - MPLUS) * LP1 + (EMINUS + MMINUS) * LP2
            
          CASE (4)
C           F4 term: Second derivative Legendre polynomial
            DERIV = 2
            N1 = ELL
            LP1 = LEGP(X, N1, DERIV)
            
            TMP = (MPLUS - EPLUS - MMINUS - EMINUS) * LP1
            
          CASE DEFAULT
C           Invalid form factor index
            TMP = (0.0D0, 0.0D0)
        END SELECT
        
C       Add contribution from this partial wave
        RESULT = RESULT + TMP
      END DO
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Get multipole amplitudes (E and M) for given parameters
C
C     These wrapper subroutines simplify accessing the multipole amplitudes
C     from the data storage. They provide a consistent interface for both
C     electric (E) and magnetic (M) multipoles.
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     Generic function to retrieve multipole amplitude (E-pole or M-pole)
C
C     Inputs:
C       TYPE_CHAR - 'E' for electric or 'M' for magnetic multipole
C       SIGN_STR  - 'plus' (J=L+1/2) or 'minus' (J=L-1/2)
C       TARGET_STR - 'p12' (proton), 'n12' (neutron), or '32q' (isospin 3/2)
C       ELL       - Angular momentum (0=S, 1=P, 2=D, 3=F, 4=G)
C       SQRTS_VAL - Center-of-mass energy in MeV
C
C     Output:
C       RESULT    - Complex amplitude value
C-----------------------------------------------------------------------
      SUBROUTINE GET_POLE(TYPE_CHAR, SIGN_STR, TARGET_STR, ELL, 
     &                    SQRTS_VAL, RESULT)
      IMPLICIT NONE
      CHARACTER*1 TYPE_CHAR
      CHARACTER*(*) SIGN_STR
      CHARACTER*3 TARGET_STR
      INTEGER ELL
      DOUBLE PRECISION SQRTS_VAL
      DOUBLE COMPLEX RESULT
      
C     External function
      DOUBLE COMPLEX GET_POLE_FUNC
      EXTERNAL GET_POLE_FUNC
      
      RESULT = GET_POLE_FUNC(TYPE_CHAR, SIGN_STR, TARGET_STR, ELL, 
     &                      SQRTS_VAL)
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Get Electric multipole amplitude (E-pole)
C
C     This is a wrapper around the generic GET_POLE function for electric multipoles
C
C     Inputs:
C       SIGN   - 'plus' (J=L+1/2) or 'minus' (J=L-1/2)
C       TARGET - 'p12' (proton), 'n12' (neutron), or '32q' (isospin 3/2)
C       ELL    - Angular momentum (0=S, 1=P, 2=D, 3=F, 4=G)
C       SQRTS  - Center-of-mass energy in MeV
C
C     Output:
C       RESULT - Complex amplitude value
C-----------------------------------------------------------------------
      SUBROUTINE GET_EPOLE(SIGN, TARGET, ELL, SQRTS, RESULT)
      IMPLICIT NONE
      CHARACTER*(*) SIGN
      CHARACTER*3 TARGET
      INTEGER ELL
      DOUBLE PRECISION SQRTS
      DOUBLE COMPLEX RESULT
      
C     Call the generic function with type='E'
      CALL GET_POLE('E', SIGN, TARGET, ELL, SQRTS, RESULT)
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Get Magnetic multipole amplitude (M-pole)
C
C     This is a wrapper around the generic GET_POLE function for magnetic multipoles
C
C     Inputs:
C       SIGN   - 'plus' (J=L+1/2) or 'minus' (J=L-1/2)
C       TARGET - 'p12' (proton), 'n12' (neutron), or '32q' (isospin 3/2)
C       ELL    - Angular momentum (0=S, 1=P, 2=D, 3=F, 4=G)
C       SQRTS  - Center-of-mass energy in MeV
C
C     Output:
C       RESULT - Complex amplitude value
C-----------------------------------------------------------------------
      SUBROUTINE GET_MPOLE(SIGN, TARGET, ELL, SQRTS, RESULT)
      IMPLICIT NONE
      CHARACTER*(*) SIGN
      CHARACTER*3 TARGET
      INTEGER ELL
      DOUBLE PRECISION SQRTS
      DOUBLE COMPLEX RESULT
      
C     Call the generic function with type='M'
      CALL GET_POLE('M', SIGN, TARGET, ELL, SQRTS, RESULT)
      
      RETURN
      END