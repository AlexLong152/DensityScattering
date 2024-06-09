C ***********************************************************************
C MODUL enthaelt Clebsch-Gordan-Koeffizienten
C  vor benutzen CALL initclebsch im Hauptprogramm
C
C ***********************************************************************
      MODULE clebsch
       USE precision
       PRIVATE
     
       PARAMETER (NFAK=150) ! vorsicht: origin and ibm haben kein overflow
       REAL(dpreal) FFAK(0:NFAK)
       REAL(dpreal) WFAK(0:NFAK)

       REAL(dpreal),PUBLIC ::  DCH(0:NFAK),SDCH(0:NFAK),SD2CH(0:NFAK),
     X                 SD2CHINV(0:NFAK),
     X                 GHI(0:NFAK),WFAKINV(0:NFAK),FFAKINV(0:NFAK)

       PUBLIC C6J,C6J2,C9J,C9J2,CG,CG000,RM1H,FFAK,WFAK,
     X        DACH,GAMMAHI,NFAK

       PUBLIC initclebsch,printclebsch
       
      CONTAINS

! subroutine that prints out NFAK and version number 
      SUBROUTINE printclebsch
       IMPLICIT NONE
       
       WRITE(*,*) 'clebsch version: ',VERREV
!       WRITE(*,*) 'Date: ',VERDATE
       WRITE(*,*) 'clebsch NFAK = ',NFAK
       
       
      END SUBROUTINE printclebsch
      
      FUNCTION C6J (I,J,K,L,M,N)
      implicit real(dpreal) (a-h,o-z)
C*******************************************************************************
C
C  6-J-SYMBOL   ( I/2  J/2  K/2 )
C               ( L/2  M/2  N/2 )
C
C*******************************************************************************
      INTEGER T,G
      LOGICAL FAC
      DIMENSION JF(2,3)
      COMMON /FAKLOG/ FL(0:100)
      SAVE /FAKLOG/
      DATA FAC /.FALSE./
      IF (.NOT. FAC) THEN
       FAC=.TRUE.
       FL(0)=0._dpreal
       DO 20 IF=1,100
   20   FL(IF)=LOG(REAL(IF,dpreal))+FL(IF-1)
      ENDIF
      C6J=0._dpreal
      CALL DELMIN(I,J,K,IJK,MIJK,DIJK)
      CALL DELMIN(I,M,N,IMN,MIMN,DIMN)
      CALL DELMIN(L,J,N,LJN,MLJN,DLJN)
      CALL DELMIN(L,M,K,LMK,MLMK,DLMK)
      NT=MIN(MIJK,MIMN,MLJN,MLMK)
      IF (NT .LT. 0) RETURN
      IF (MIN(I,J,K,L,M,N) .EQ. 0) THEN
       JF(1,1)=I
       JF(1,2)=J
       JF(1,3)=K
       JF(2,1)=L
       JF(2,2)=M
       JF(2,3)=N
       DO 10 JJ=1,3
        DO 10 II=1,2
         IF (JF(II,JJ) .EQ. 0) THEN
          I0=II
          J0=JJ
         ENDIF
   10    CONTINUE
       G=JF(3-I0,J0)
       IF (J0 .EQ. 3) THEN
        K0=1
       ELSE
        K0=3
       ENDIF
       JK=JF(1,K0)
       JG=JF(2,K0)
       C6J=RM1H((JK+JG+G)/2)/SQRT((JK+1._dpreal)*(JG+1._dpreal))
      ELSE
       NTMIN=MAX(IJK,IMN,LJN,LMK)
       W=0.5_dpreal*(DIJK+DIMN+DLJN+DLMK)
       IJLM=(I+J+L+M)/2
       JKMN=(J+K+M+N)/2
       KINL=(K+I+N+L)/2
       DO 30 T=NTMIN,NTMIN+NT
        SL=W+FL(T+1)-FL(T-IJK)-FL(T-IMN)-FL(T-LJN)-FL(T-LMK)-FL(IJLM-T)
     .     -FL(JKMN-T)-FL(KINL-T)
        S=RM1H(T)*EXP(SL)
   30   C6J=C6J+S
      ENDIF
      END FUNCTION 
      FUNCTION C6J2 (I,J,K,L,M)
      implicit real(dpreal) (a-h,o-z)
C*******************************************************************************
C
C  6-J-SYMBOL   ( I/2  J/2  1/2 )
C
C               ( K/2  L/2  M/2 )
C
C*******************************************************************************
      INTEGER G
      JK=MIN(I,J)
      JG=MIN(K,L)
      IF ((J-I)*(L-K) .GE. 0) THEN
       G=M-1
       F=(1+(G+JK-JG)/2)*(1+(G-JK+JG)/2)
      ELSE
       G=M
       F=(1+(JK+JG-G)/2)*(2+(JK+JG+G)/2)
      ENDIF
      IF (ABS(JK-JG) .LE. G .AND. G .LE. JK+JG) THEN
       C6J2=RM1H(1+(JK+JG+G)/2)*SQRT(F/((JK+1)*(JK+2)*(JG+1)
     .        *(JG+2)))
      ELSE
       C6J2=0._dpreal
      ENDIF
      END FUNCTION
 
      FUNCTION C9J (J1,J2,J5,J3,J4,J6,J7,J8,J9)
      implicit real(dpreal) (a-h,o-z)
C     ===============
CCCCCC
CCCCCC          ALL L9 MUST BE TWICE AS LARGE AS TRUE ARGUMENTS
CCCCCC
CCCC  TAMURAS NUMBERING CONVENTION IS USED HERE.
CCCC                       1  2  5
CCCC                       3  4  6
CCCC                       7  8  9
      DIMENSION LT(9)
      U9=0.000_dpreal
      LT(1)=J1
      LT(2)=J2
      LT(3)=J3
      LT(4)=J4
      LT(5)=J5
      LT(6)=J6
      LT(7)=J7
      LT(8)=J8
      LT(9)=J9
      LMIN=LT(1)
      IMIN=1
      DO20I=2,9
      IF(LT(I)-LMIN)15,20,20
   15 LMIN=LT(I)
      IMIN=I
   20 CONTINUE
      KEX=0
      GOTO(110,110,110,110,150,150,170,170,190),IMIN
  110 MM=(IMIN-1)/2+1
      M1=MM+MM-1
      M2=M1+1
      M3=MM+4
      L1=LT(7)
      LT(7)=LT(M1)
      LT(M1)=L1
      L1=LT(8)
      LT(8)=LT(M2)
      LT(M2)=L1
      L1=LT(9)
      LT(9)=LT(M3)
      LT(M3)=L1
      IMIN=IMIN+(7-M1)
      GOTO175
  150 KEX=1
      M1=7
      M2=8
      M3=IMIN+IMIN-9
      M4=M3+1
      GOTO 180
  170 KEX=1
  175 M1=5
      M2=6
      M3=IMIN-6
      M4=M3+2
  180 L1=LT(M1)
      L1=LT(M1)
      LT(M1)=LT(M3)
      LT(M3)=L1
      L1=LT(M2)
      LT(M2)=LT(M4)
      LT(M4)=L1
      L1=LT(9)
      LT(9)=LT(IMIN)
      LT(IMIN)=L1
  190 IF(LT(9))200,200,300
  200 IF(LT(5)-LT(6))1000,210,1000
  210 IF(LT(7)-LT(8))1000,220,1000
  220 RT=(LT(5)+1)*(LT(7)+1)
      K=(LT(5)+LT(7)+LT(2)+LT(3))/2
      RAC=RACAH (LT(1),LT(2),LT(5),LT(4),LT(3),LT(7))
      U9=(RAC/ SQRT(RT))*RM1H(K)
      GOTO370
  300 K1=IABS(LT(2)-LT(7))
      K2=IABS(LT(3)-LT(5))
      K3=IABS(LT(4)-LT(9))
      NMIN=MAX0(K1,K2,K3)
      K1=LT(2)+LT(7)
      K2=LT(3)+LT(5)
      K3=LT(4)+LT(9)
      NMAX=MIN0(K1,K2,K3)
      IF(NMIN-NMAX)320,320,1000
  320 DO350N=NMIN,NMAX,2
      W1=N+1
      RAC=RACAH (LT(2),LT(5),LT(1),LT(3),LT(7),N)
      IF(RAC)321,350,321
  321 W1=W1*RAC
      RAC=RACAH (LT(2),LT(4),LT(8),LT(9),LT(7),N)
      IF(RAC)322,350,322
  322 W1=W1*RAC
      RAC=RACAH (LT(3),LT(4),LT(6),LT(9),LT(5),N)
      IF(RAC)323,350,323
  323 U9=U9+W1*RAC
  350 CONTINUE
      U9=EIMIHO (LT(2)+LT(3)+LT(4)+LT(5)+LT(7)+LT(9),U9)
  370 IF(KEX)400,1000,400
  400 KP=0
      DO410I=1,9
  410 KP=KP+LT(I)
      U9=U9*RM1H(KP/2)
 1000 C9J=U9
      RETURN
      END FUNCTION
      FUNCTION C9J2 (I,J,K,L,M,N,O,P)
      implicit real(dpreal) (a-h,o-z)
C*******************************************************************************
C
C               ( I/2  J/2  K/2 )
C  9-J-SYMBOL   ( L/2  1/2  M/2 )
C               ( N/2  O/2  P/2 )
C
C*******************************************************************************
      INTEGER O,P,G,GMIN,GMAX
      C9J2=0._dpreal
      GMIN=MAX(ABS(L-P),ABS(K-1),ABS(O-I))
      GMAX=MIN(L+P,K+1,O+I)
      DO 10 G=GMIN,GMAX,2
   10  C9J2=C9J2+RM1H(G)*(G+1)*C6J2(L,M,K,G,P)*C6J2(K,G,O,J,I)
     .      *C6J(N,O,P,G,L,I)
      C9J2=RM1H((I+J+K+L+1+M+N+O+P)/2)*C9J2
      END FUNCTION
 
      FUNCTION CG (I,J,K,L,M,N)
      implicit real(dpreal) (a-h,o-z)
C*******************************************************************************
C
C  CLEBSCH-GORDON   < I/2 J/2 L/2 M/2 / K/2 N/2 >
C
C*******************************************************************************
      INTEGER T
      IF (L+M .NE. N .OR. ABS(I-J) .GT. K .OR. K .GT. I+J .OR.
     .    ABS(L) .GT. I .OR. ABS(M) .GT. J .OR. ABS(N) .GT. K) THEN
       CG=0._dpreal
      RETURN
      ENDIF
      I1=(I+L)/2
      I2=(I-L)/2
      I3=(J+M)/2
      I4=(J-M)/2
      I5=(K-N)/2
      I6=(K+N)/2
      I7=(I+J-K)/2
      I8=(J+K-I)/2
      I9=(K+I-J)/2
      I10=(I+J+K+2)/2
      XX=WFAK(I7)*WFAK(I8)*WFAK(I9)*WFAK(I1)*WFAK(I2)*
     X     WFAK(I3)*WFAK(I4)*WFAK(I5)*
     .   WFAK(I6)*WFAKINV(I10)
      J1=(K-J+L)/2
      J2=(K-I-M)/2
      NT=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9)+1
      IT=0
      SUM=0._dpreal
      T=-1
   10   T=T+1
        L1=J1+T
        L2=J2+T
        L3=I7-T
        L4=I2-T
        L5=I3-T
        IF (MIN(L1,L2,L3,L4,L5) .LT. 0) GOTO 10
       SUM=SUM+RM1H(T)*FFAKINV(T)*FFAKINV(L1)*FFAKINV(L2)
     X       *FFAKINV(L3)*FFAKINV(L4)*FFAKINV(L5)
       IT=IT+1
       IF (IT .LT. NT) GOTO 10
      CG=XX*SUM*SQRT(K+1._dpreal)
      END FUNCTION
      FUNCTION CG000 (I,J,K)
      implicit real(dpreal) (a-h,o-z)
C*******************************************************************************
C
C  CLEBSCH-GORDON   < I J 0 0 / K 0 >
C
C  DREIECKSUNGLEICHUNG UND PHASE VON (I+J+K) WERDEN NICHT GETESTET.
C
C*******************************************************************************
 
      IP=(I+J+K)/2
      CG000=RM1H(IP+I-J)*WFAK(I+J-K)*WFAK(J+K-I)*
     X     WFAK(K+I-J)*WFAKINV(I+J+K+1)
     .      *FFAK(IP)*FFAKINV(IP-I)*FFAKINV(IP-J)*FFAKINV(IP-K)*SDCH(K)
      END FUNCTION
      SUBROUTINE DELMIN (J1,J2,J3,MX,MNPER,DELTA)
      implicit real(dpreal) (a-h,o-z)
      COMMON /FAKLOG/ FL(0:100)
      MX=J1+J2+J3
      IF (MOD(MX,2) .NE. 0) THEN
       MNPER=-1
      ELSE
       MX=MX/2
       I=(J1+J2-J3)/2
       J=(J2+J3-J1)/2
       K=(J3+J1-J2)/2
       MNPER=MIN(I,J,K)
       IF (MNPER .GE. 0) DELTA=FL(I)+FL(J)+FL(K)-FL(MX+1)
      ENDIF
      END SUBROUTINE
 
      FUNCTION EIMIHO(K,FAC)
      implicit real(dpreal) (a-h,o-z)
C     ===============
      EIMIHO=FAC
      IF( (K/2) * 2 .NE. K) EIMIHO=-FAC
      RETURN
      END FUNCTION
 
      SUBROUTINE initclebsch
      implicit real(dpreal) (a-h,o-z)
C*******************************************************************************
C
C  I!  SQRT(I!)   0 <= I <= NFAK
C
C*******************************************************************************
      FFAK(0)=1.0_dpreal
      FFAKINV(0)=1.0_dpreal
      DO 10 I=1,NFAK
       FFAK(I)=I*FFAK(I-1)
 10    FFAKINV(I)=1.0/FFAK(I)
      DO 20 I=0,NFAK
       WFAK(I)=SQRT(FFAK(I))
 20    WFAKINV(I)=1.0/WFAK(I)
       CALL DACH
      END SUBROUTINE
 
      FUNCTION RACAH (IA,IB,IE,ID,IC,IF)
      implicit real(dpreal) (a-h,o-z)
C     ==============
CCC  ALL  L*  MUST BE TWICE AS LARGE AS TRUE ARGUMENTS
      DIMENSION LT(6),FACLOG(NFAK),FACLOGEXP(NFAK),FACLOGEXPINV(NFAK),
     X                WFACLOGEXP(NFAK),WFACLOGEXPINV(NFAK)
      DATA FN/-1.0_dpreal/
      REAL(dpreal) SUM
      INTEGER F1,F2
      SAVE FN,FACLOG,FACLOGEXP,FACLOGEXPINV,WFACLOGEXP,WFACLOGEXPINV
      IF(FN)1,1,11
    1 FACLOG(1)=0.000_dpreal
      FACLOG(2)=0.000_dpreal
      FACLOGEXP(1)=1.000_dpreal
      FACLOGEXP(2)=1.000_dpreal
      FACLOGEXPINV(1)=1.000_dpreal
      FACLOGEXPINV(2)=1.000_dpreal
      WFACLOGEXP(1)=1.000_dpreal
      WFACLOGEXP(2)=1.000_dpreal
      WFACLOGEXPINV(1)=1.000_dpreal
      WFACLOGEXPINV(2)=1.000_dpreal
      FN=1.000_dpreal
      DO 5 N=3,NFAK
      FN=FN+1.000_dpreal
      FACLOG(N)=FACLOG(N-1)+LOG(FN)
      FACLOGEXP(N)=EXP(FACLOG(N))
      FACLOGEXPINV(N)=1/FACLOGEXP(N)
      WFACLOGEXP(N)=SQRT(FACLOGEXP(N))
    5 WFACLOGEXPINV(N)=SQRT(FACLOGEXPINV(N))
   11 K1=IA+IB-IE
      K3=IC+ID-IE
      K5=IA+IC-IF
      K7=IB+ID-IF
      K2=IE-IABS(IA-IB)
      K4=IE-IABS(IC-ID)
      K6=IF-IABS(IA-IC)
      K8=IF-IABS(IB-ID)
      K9=MIN0(K1,K2,K3,K4,K5,K6,K7,K8)
      SUM=0.0_dpreal
      IF(K9)4000,20,20
   20 K2=K1-2*(K1/2)
      K4=K3-2*(K3/2)
      K6=K5-2*(K5/2)
      K8=K7-2*(K7/2)
      IF(MAX0(K2,K4,K6,K8))4000,25,4000
   25 LTMIN=MIN0(IA,IB,IC,ID,IE,IF)
      IF(LTMIN)4000,30,150
   30 LT(1)=IA
      LT(2)=IB
      LT(3)=IC
      LT(4)=ID
      LT(5)=IE
      LT(6)=IF
      LTMIN=LT(1)
      KMIN=1
      DO40N=2,6
      IF(LT(N)-LTMIN)35,40,40
   35 LTMIN=LT(N)
      KMIN=N
   40 CONTINUE
      S1=1.000_dpreal
      F1=IE
      F2=IF
      GOTO(55,55,55,55,45,50),KMIN
   45 F1=IA
      F2=IC
      S1=RM1H(K5/2)
      GOTO55
   50 F1=IA
      F2=IB
      S1=RM1H(K1/2)
   55 SUM=S1*SD2CHINV(F1)*SD2CHINV(F2)
      GOTO4000
  150 IABEP=(IA+IB+IE)/2+1
      ICDEP=(IC+ID+IE)/2+1
      IACFP=(IA+IC+IF)/2+1
      IBDFP=(IB+ID+IF)/2+1
      IABE=IABEP-IE
      IEAB=IABEP-IB
      IBEA=IABEP-IA
      ICDE=ICDEP-IE
      IECD=ICDEP-ID
      IDEC=ICDEP-IC
      IACF=IACFP-IF
      IFAC=IACFP-IC
      ICFA=IACFP-IA
      IBDF=IBDFP-IF
      IFBD=IBDFP-ID
      IDFB=IBDFP-IB
      NZMAX=MIN0(IABE,ICDE,IACF,IBDF)
      IABCD1=(IA+IB+IC+ID+4)/2
      IEFMAD=(IE+IF-IA-ID)/2
      IEFMBC=(IE+IF-IB-IC)/2
      NZMI1=-IEFMAD
      NZMI2=-IEFMBC
      NZMIN=MAX0(0,NZMI1,NZMI2)+1
      SQLOG=WFACLOGEXP(IABE)*WFACLOGEXP(IEAB)
     X     *WFACLOGEXP(IBEA)*WFACLOGEXP(ICDE)
     1     *WFACLOGEXP(IECD)*WFACLOGEXP(IACF)
     X     *WFACLOGEXP(IFAC)*WFACLOGEXP(IDEC)
     2     *WFACLOGEXP(ICFA)*WFACLOGEXP(IBDF)
     X     *WFACLOGEXP(IFBD)*WFACLOGEXP(IDFB)
     3     *WFACLOGEXPINV(IABEP+1)*WFACLOGEXPINV(ICDEP+1)
     X     *WFACLOGEXPINV(IACFP+1)*WFACLOGEXPINV(IBDFP+1)

      NZM1=NZMIN-1
      K1=IABCD1-NZM1
      K2=IABE-NZM1
      K3=ICDE-NZM1
      K4=IACF-NZM1
      K5=IBDF-NZM1
      K6=NZMIN
      K7=IEFMAD+NZMIN
      K8=IEFMBC+NZMIN
      SSTERM=-(-1)**NZMIN*SQLOG
      DO 200 NZ=NZMIN,NZMAX
      SUM=SUM+SSTERM
     X     *FACLOGEXP(K1)*FACLOGEXPINV(K2)
     X     *FACLOGEXPINV(K3)*FACLOGEXPINV(K4)*FACLOGEXPINV(K5)
     X     *FACLOGEXPINV(K6)*FACLOGEXPINV(K7)*FACLOGEXPINV(K8)
      SSTERM=-SSTERM
      K1=K1-1
      K2=K2-1
      K3=K3-1
      K4=K4-1
      K5=K5-1
      K6=K6+1
      K7=K7+1
      K8=K8+1
  200 CONTINUE
4000  RACAH=SUM*RM1H((IA+IB+IC+ID)/2)
      RETURN
      END FUNCTION
      FUNCTION RM1H (I)
      implicit real(dpreal) (a-h,o-z)
C*******************************************************************************
C
C  RM1H = (-1) ** I
C
C*******************************************************************************
      RM1H=(-1.0_dpreal)**I    ! 1-2*MOD(ABS(I),2)
      END FUNCTION

 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  new add for 3nf CCCCCCCCCCCCCCCCCCCCCCCC


      SUBROUTINE DACH
      implicit real(dpreal) (a-h,o-z)
C**********************************************************************
C
C  D(I) = 2*I-1   0 <= I <= NFAK
C
C**********************************************************************
      DO 10 I=0,NFAK
       DCH(I)=2*I+1
       SDCH(I)=SQRT(DCH(I))
       SD2CH(I)=SQRT(real(I+1,dpreal))
   10  SD2CHINV(I)=1.0/SQRT(real(I+1,dpreal))
      END SUBROUTINE


      SUBROUTINE GAMMAHI
      implicit real(dpreal) (a-h,o-z)
C**********************************************************************
C
C  Gamma(I+1/2)   0 <= I <= NFAK  (AS 6.1.12)
C
C  (Gamma(I) = (I-1)!)             (AS 6.1.6)
C
C**********************************************************************
      FAK=SQRT(ACOS(-1.0_dpreal))
      GHI(0)=1.0_dpreal
      GHI(1)=.5_dpreal
      DO 10 I=2,NFAK
   10  GHI(I)=GHI(I-1)*DCH(I-1)*.5_dpreal
      DO 20 I=0,NFAK
   20  GHI(I)=GHI(I)*FAK
      END SUBROUTINE 

CCCCCCCCCCCCCCCCCCC   new add for 3nf end CCCCCCCCCCCCCCCCCCCCCCCCCCCC 
     
      END MODULE clebsch

