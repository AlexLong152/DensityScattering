C *** spherical besselfunctions up to l=16
C *** based on analytical formulars and low x expansion

      MODULE besselfu
       USE precision
      PRIVATE

      PUBLIC SBF1


      CONTAINS

      SUBROUTINE SBF1(LM,P,NR,BL,RP)
C**** FOR SUB. BESSEL (940825)
      REAL(dpreal) BL(NR),RP(NR)
      REAL(dpreal)  P,X
      INTEGER LM,NR
      INTEGER MRM,IR,III

      REAL(dpreal) DF(19)
      DATA DF/1.,3.,15.,105.,945.,10395.,135135. ,2027025.,
     X        34459425.,654729075.,13749310575.,316234143225.,
     X        7905853580625.,213458046676875.,6190283353629375.,
     X        191898783962510625.,6332659870762850625.,
     X        221643095476699771875.,8200794532637891559375./
      SAVE DF
C
      IF(NR.EQ.0) RETURN  ! nothing to be done on this processor 
        
      III=1
      IF(P*RP(1).EQ.0. .AND. LM.EQ.0 ) THEN
        III=2
        BL(1)=1.0
      END IF
      IF(P*RP(1).EQ.0. .AND. LM.NE.0 ) THEN
        III=2
        BL(1)=0.0
      END IF
C
      MRM=NR+1
      DO IR=III,NR
       IF(P*RP(IR).GT.(LM*0.3+0.05)) THEN
         MRM=IR
         GOTO 1000
       END IF
      END DO
 1000 CONTINUE

      DO IR=III,MRM-1
       X=P*RP(IR)
       BL(IR)
     &      =(X**LM/DF(LM+1))*(1.-0.5*X**2/(2.*LM+3.)+
     1      (0.5*X**2)**2/(2.*(2.*LM+3.)*(2.*LM+5.))-
     2      (0.5*X**2)**3/(6.*(2.*LM+3.)*(2.*LM+5.)*(2.*LM+7.))+
     3      (0.5*X**2)**4/(24.*(2.*LM+3.)*(2.*LM+5.)*(2.*LM+7.)*(2.*
     4      LM+9.)))
      END DO
 
 1001 GOTO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)LM+1
 
      STOP 'LM not implemented'

 1    DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=SIN(X)/X
      END DO
      GOTO 2000

 2    DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=SIN(X)/X**2-COS(X)/X
      END DO
      GOTO 2000

 3    DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=(3./X**3-1./X)*SIN(X)-3.*COS(X)/X**2
      END DO
      GOTO 2000

 4    DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=(15./X**4-6./X**2)*SIN(X)-(15./X**3-1./X)*COS(X)
      END DO
      GOTO 2000

 5    DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=(105./X**5-45./X**3+1./X)*SIN(X)+(10./X**2-
     1         105./X**4)*COS(X)
      END DO
      GOTO 2000

 6    DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=(15./X**2-420./X**4+945./X**6)*SIN(X)+(-1./X+
     1        105./X**3-945./X**5)*COS(X)
      END DO
      GOTO 2000

 7    DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     -      ((-10395./x**5 + 1260./x**3 - 21./x)*COS(X) +(- 1. + 
     -      10395./x**6 - 4725./x**4 + 210./x**2)
     -      *Sin(x) )/x
      END DO
      GOTO 2000

 8    DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     -      ((1. - 135135./x**6 + 17325./x**4 - 378./x**2)*Cos(x) + 
     -      (  135135./x**7 - 62370./x**5 + 3150./x**3 - 
     -      28./x)*Sin(x))/x
      END DO
      GOTO 2000

 9    DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     X      (9*x*(-225225.0+30030.0*x**2
     X              -770.0*x**4+4.0*x**6)*cos(x)
     X       +(2027025.0-945945.0*x**2+51975.0*x**4
     X         -630.0*x**6+x**8)*sin(x))/x**9
      END DO
      GOTO 2000

 10   DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     X      (-x*(34459425.0-4729725.0*x**2
     X           +135135.0*x**4-990.0*x**6+x**8)*cos(x)
     X       +45.0*(765765.0-360360.0*x**2+21021.0*x**4
     X         -308.0*x**6+x**8)*sin(x))/x**10
      END DO
      GOTO 2000

 11   DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     X      -(55.0*x*(11904165.0-1670760.0*x**2
     X           +51597.0*x**4-468.0*x**6+x**8)*cos(x)
     X       +(-654729075.0+310134825.0*x**2-18918900.0*x**4
     X         +315315.0*x**6-1485.0*x**8+x**10)*sin(x))/x**11
      END DO
      GOTO 2000

 12   DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     X      (x*(-13749310575.0+1964187225.0*x**2
     X           -64324260.0*x**4+675675.0*x**6-2145.0*x**8
     X          +x**10)*cos(x)
     X       -33.0*(-416645775.0+198402750.0*x**2-12530700.0*x**4
     X         +229320.0*x**6-1365.0*x**8+2.0*x**10)*sin(x))/x**12
      END DO
      GOTO 2000

 13   DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     X      (39.0*x*(-8108567775.0+1175154750.0*x**2
     X           -40291020.0*x**4+471240.0*x**6-1925.0*x**8
     X          +2.0*x**10)*cos(x)
     X       +(316234143225.0-151242416325.0*x**2+9820936125.0*x**4
     X         -192972780.0*x**6+1351350.0*x**8-3003.0*x**10
     X         +x**12)*sin(x))/x**13
      END DO
      GOTO 2000

 14   DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     X      (-x*(7905853580625.0-1159525191825.0*x**2
     X           +41247931725.0*x**4-523783260.0*x**6
     X           +2552550.0*x**8
     X          -4095.0*x**10+x**12)*cos(x)
     X       +91.0*(86877511875.0-41701205700.0*x**2
     X         +2770007625.0*x**4
     X         -57558600.0*x**6+454410.0*x**8-1320.0*x**10
     X         +x**12)*sin(x))/x**14
      END DO
      GOTO 2000

 15   DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     X      -(105.0*x*(2032933777875.0-301175374500.0*x**2
     X           +11043097065.0*x**4-149652360.0*x**6
     X           +831402.0*x**8
     X          -1768.0*x**10+x**12)*cos(x)
     X       +(-213458046676875.0+102776096548125.0*x**2
     X         -6957151150950.0*x**4
     X         +151242416325.0*x**6-1309458150.0*x**8
     X         +4594590.0*x**10
     X         -5460.0*x**12+x**14)*sin(x))/x**15
      END DO
      GOTO 2000

 16   DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     X      (x*(-6190283353629375.0+924984868933125.0*x**2
     X           -34785755754750.0*x**4+496939367925.0*x**6
     X           -3055402350.0*x**8
     X          +7936110.0*x**10-7140.0*x**12+x**14)*cos(x)
     X       -15.0*(-412685556908625.0+199227510231750.0*x**2
     X         -13703479539750.0*x**4
     X         +309206717820.0*x**6-2880807930.0*x**8
     X         +11639628.0*x**10
     X         -18564.0*x**12+8.0*x**14)*sin(x))/x**16
      END DO
      GOTO 2000

 17   DO IR=MRM,NR
       X=P*RP(IR)
       BL(IR)=
     X      (17.0*x*(-11288163762500625.0+1699293469623750.0*x**2
     X           -65293049571750.0*x**4+974390917500.0*x**6
     X           -6495939450.0*x**8
     X          +19606860.0*x**10-23940.0*x**12+8.0*x**14)*cos(x)
     X       +(191898783962510625.0-92854250304440625.0*x**2
     X         +6474894082531875.0*x**4
     X         -150738274937250.0*x**6+1490818103775.0*x**8
     X         -6721885170.0*x**10
     X         +13226850.0*x**12-9180.0*x**14
     X         +x**16)*sin(x))/x**17
      END DO
      GOTO 2000

 2000 CONTINUE
      RETURN
      END SUBROUTINE sbf1

      END MODULE besselfu
