C     No dummy program needed

C-----------------------------------------------------------------------
C     Legendre polynomial calculations for pion photoproduction
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     Calculate Legendre polynomial P_n(x) and derivatives
C     DERIV = 0: regular Legendre polynomial
C     DERIV = 1: first derivative
C     DERIV = 2: second derivative
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION LEGP(X, N, DERIV)
      IMPLICIT NONE
      DOUBLE PRECISION X, X2, X3, X4
      INTEGER N, DERIV
      
C     Store common powers to avoid repeating calculations
      X2 = X * X
      X3 = X2 * X
      X4 = X2 * X2
      
      IF (DERIV .EQ. 0) THEN
C       Regular Legendre polynomial
        IF (N .EQ. -1 .OR. N .EQ. 0) THEN
          LEGP = 1.0D0
        ELSE IF (N .EQ. 1) THEN
          LEGP = X
        ELSE IF (N .EQ. 2) THEN
          LEGP = 1.5D0 * X2 - 0.5D0
        ELSE IF (N .EQ. 3) THEN
          LEGP = 2.5D0 * X3 - 1.5D0 * X
        ELSE IF (N .EQ. 4) THEN
          LEGP = (1.0D0 / 8.0D0) * (3.0D0 - 30.0D0 * X2 + 35.0D0 * X4)
        ELSE
          WRITE(*,*) 'Error: Legendre polynomial not implemented'
          WRITE(*,*) 'n=', N
          STOP 1
        ENDIF
      ELSE IF (DERIV .EQ. 1) THEN
C       First derivative
        IF (N .EQ. -1 .OR. N .EQ. 0) THEN
          LEGP = 0.0D0
        ELSE IF (N .EQ. 1) THEN
          LEGP = 1.0D0
        ELSE IF (N .EQ. 2) THEN
          LEGP = 3.0D0 * X
        ELSE IF (N .EQ. 3) THEN
          LEGP = 0.5D0 * (15.0D0 * X2 - 3.0D0)
        ELSE IF (N .EQ. 4) THEN
          LEGP = (-60.0D0 * X + 140.0D0 * X3) / 8.0D0
        ELSE
          WRITE(*,*) 'Error: Legendre derivative not implemented'
          WRITE(*,*) 'n=', N
          STOP 1
        ENDIF
      ELSE IF (DERIV .EQ. 2) THEN
C       Second derivative
        IF (N .EQ. -1 .OR. N .EQ. 0 .OR. N .EQ. 1) THEN
          LEGP = 0.0D0
        ELSE IF (N .EQ. 2) THEN
          LEGP = 3.0D0
        ELSE IF (N .EQ. 3) THEN
          LEGP = 15.0D0 * X
        ELSE IF (N .EQ. 4) THEN
          LEGP = -7.5D0 + 52.5D0 * X2
        ELSE
          WRITE(*,*) 'Error: Legendre 2nd derivative not implemented'
          WRITE(*,*) 'n=', N
          STOP 1
        ENDIF
      ELSE
        WRITE(*,*) 'Error: Invalid derivative order'
        WRITE(*,*) 'deriv=', DERIV
        STOP 1
      ENDIF
      
      RETURN
      END