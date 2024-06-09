cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for One/Twobody Contributions to Few-Nucleon Processes Calculated Via 1N/2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              gausleg :
c              trns    : populates fields xp and ap with ransformed Gauss-Legendre points and weights
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, identical to file of same name in common-densities/ of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gausleg(N,x1,x2,X,W)
      IMPLICIT NONE
      INTEGER N        
      REAL(8) x1,x2,X(N),W(N)
      REAL(8) z1,z,xm,xl,pp,p3,p2,p1,pi,tol
      INTEGER m,i,j
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      pi=acos(-1.0)
      tol=100.0*1.E-10
      m=(n+1)/2
      xm=0.5*(x2+x1)
      xl=0.5*(x2-x1)
      
      DO 10 i=1,m
         z=cos(pi*(i-0.25)/(N+0.5))
         
 20      CONTINUE
         p1=1.0E0
         p2=0.0E0
         DO 30 j=1,N
            p3=p2
            p2=p1
            p1=((2*j-1)*z*p2-(j-1)*p3)/j
 30      CONTINUE
         pp=N*(z*p1-p2)/(z*z-1.0E0)
         z1=z
         z=z1-p1/pp
         IF( abs(z1-z) .GT. tol) GOTO 20 ! Scheifenende
         
         X(i) = xm - xl*z
         X(n+1-i) = xm + xl*z
         W(i) = 2.E0*xl/((1.0-z*z)*pp*pp)
         W(n+1-i) = W(i)
 10   CONTINUE
      END SUBROUTINE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE TRNS(NP1,NP2,NP,P1,P2,P3,XP,AP)
      IMPLICIT REAL(8)(a-h,o-z)
      IMPLICIT INTEGER(i-n)
c     TRNS distributes the arrays XP and AP with transformed Gauß-Legendre points and weights.

c     NP1 points are mapped onto the interval [0;2π] via the hyperbolic transformation

c     X --> (1.+X) / (1./P1-(1./P1-2./P2)*X)

c     with

c     NP1/2 points in [0;π] and
c     NP1/2 points in [π;2π]

c     NP2 points are mapped onto the interval [p2;p3] via the linear transformation

c     X --> (P3+P2)/2. + (P3-P2)/2.*X

c     Finally,

c     NP = NP1 + NP2
C     ===============
C     
C     TRNS BELEGT DIE FELDER XP UND AP MIT TRANSFORMIERTEN
C     GAUSS-LEGENDRE-PUNKTEN UND GEWICHTEN
C     
C     NP1 PUNKTE WERDEN UEBER DIE HYPERBOLISCHE TRANSFORMATION
C     
C     X --> (1.+X) / (1./P1-(1./P1-2./P2)*X)
C
C     AUF DAS INTERVALL (0.;P2) ABGEBILDET, WOBEI
C     NP1/2 PUNKTE IN (0.;P1) UND
C     NP1/2 PUNKTE IN (P1;P2) LIEGEN
C
C     NP2 PUNKTE WERDEN UEBER DIE LINEARE TRANSFORMATION
C
C     X --> (P3+P2)/2. + (P3-P2)/2.*X
C
C     AUF DAS INTERVALL (P2;P3) ABGEBILDET
C
C     NP = NP1 + NP2
C     
      DIMENSION XP1(1000),AP1(1000),XP2(1000),AP2(1000)
      DIMENSION XP(NP),AP(NP)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      IF(NP1.GT.1000) STOP 'NP1 in TRNS'
      IF(NP2.GT.1000) STOP 'NP2 in TRNS'
      
      
      CALL gausleg(NP1,REAL(-1.0,8),REAL(1.0,8),XP1,AP1)
      
      DO 1 I=1,NP1
         X=XP1(I)
         A=AP1(I)
         XX=1./P1-(1./P1-2./P2)*X
         XP1(I)=(1.+X) / XX
 1       AP1(I)=(2./P1-2./P2)*A / XX**2
C     
         IF(NP2 .NE. 0) THEN
            
            CALL gausleg(NP2,REAL(-1.0,8),REAL(1.0,8),XP2,AP2)
            
            
            DO 2 I=1,NP2
               X=XP2(I)
               A=AP2(I)
               DELPH=(P3-P2)/2.
               XP2(I)=(P3+P2)/2. + DELPH*X
 2             AP2(I)=DELPH*A
            ENDIF
C     
            DO 3 I=1,NP1
               XP(I)=XP1(I)
 3             AP(I)=AP1(I)
C     
      IF(NP2 .NE. 0) THEN
         DO 4 I=1,NP2
            XP(I+NP1)=XP2(I)
 4          AP(I+NP1)=AP2(I)
         ENDIF
C     RETURN
         END SUBROUTINE

      
