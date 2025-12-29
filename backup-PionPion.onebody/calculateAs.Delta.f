c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
      subroutine CalculateAsDelta(A1,A2,A3,A4,A5,A6,wx,xq,Nx,t,omega,verbosity)
c
c     Calculate O(e^2 delta^3) piece of RCS amplitude on a nucleon, i.e., 
c     leading Delta-tree and Delta-loop effects. Returns A1-A6, given input kinematics.
c
c     Written by Harald Griesshammer Feb. 2013, see README for 
c                     version history to that point, previous authors, etc.
c
c     wrapper main-programme to test routine is commented-out at end of this file
c      
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c     hgrie May 2022: set nucleon masses to Mnucleon in Dfac, delbeta,
c     to comply with deuteron implementation.
c      
c     Version 7 hgrie Feb 2013 -- compatibility to previous versions broken
c   
c    CLARIFICATION by hgrie Nov 2014:
c         calculates the Compton amplitudes A1 to A6
c         coming from Δ(1232) and Δπ loops:
c         1-loop results of [HHKK] in "strict" HBChiPT
c         (nonrelativistic Δ propagator without recoil; no E2 coupling),
c             no shift for pion-production threshold.
c         These amplitudes are "structure" only; there is also no "poles" part.
c         
c         They result in the contributions to the *static* polarisabilities
c         which are coded at the end of this file [HildeDiplom (5.10/11)].
c   
c         Only betaM1 and gammaM1M1 get contributions from Δ(1232) tree:
c     
c         betaM1Δ     =  alphaem*2*b1**2/(9*delta*Mnucleon**2) = (2*delta)*gammaM1M1Δ 
c     
c         alphaE1Δ    =  gammaE1E1Δ  = gammaE1M2Δ  = gammaM1E2Δ  = 0
c     
c         
c         All polarisabilities get contributions from Δπ loop diagrams:
c         
c         alphaE1Δπ   = -alphaem*delalpha ! (defined below)
c         
c         betaM1Δπ    = -alphaem*delbeta ! (defined below)
c         
c         gammaE1E1Δπ = [longish --defined at end of file]
c   
c         gammaM1M1Δπ = [longish --defined at end of file]
c   
c         gammaE1M2Δπ = gammaM1E2Δπ = -gammaM1M1Δπ
c   
c         Since the static values for alphaE1 and betaM1 are huge,
c         subtract these static contributions from alpha and beta --
c         but keep the static contributions to the spin-polarisabilities.
c         This is done at the end of this subroutine. 
c
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
c      
c**********************************************************************
c
      IMPLICIT NONE
      INCLUDE "../common-densities/constants.def"
c
c
c**********************************************************************
c
c  Output variables:
c
c     A1,A2,A3,A4,A5,A6-six invariant functions for gamma-N
c       scattering, calculated as per BKM review, p.66.
c
      REAL*8 A1,A2,A3,A4,A5,A6
      
c      offset ADDED to scalar polarisabilities so that their contribution from
c      Delta and Deltapi diagrams is ZERO -- included by hgrie
c
      real*8 delalpha,delbeta
c
c**********************************************************************
c
c  Input variables:
c
c     t-Mandelstam variable
c     omega-photon energy
c     xq,wx,Nx-quadratures for angular integration
c
      REAL*8 t,omega,xq(7000),wx(7000)
      INTEGER Nx
      integer verbosity
c
c**********************************************************************
c
c  Local variables:
c
C     BB1, BB2 - convenient functions
c     UDASINH - double precision inverse arcsinh function 
c                (was not an intrinsic function with the
c                 f77 compiler used when writing this code)
c     FTN - a commonly occuring function - t (1 - x) (1 - z) z
c     RFTN1/2 - a commonly occuring function - [-1+((rdelpi-+omega x (1-z))^2/
c                                               1-FTN(x,z) ]
c     RR - convenient function
c     awx - convenient function - SQRT((rdelpi-omegas*x)**2 -1.0d0))
c     z,x-arguments of BB1, BB2
c     AZ1..AZ6-integrands
c     sumAZ1..sumAZ6-results of integrals
c     Dpole(:) - contribution from pole
c     dum - a dummy variable
c     i,j-loop variables
c
c     omegas and ts are scsaled versions of omega and t
c     omegas=omega/mpi		ts=t/mpi**2

c      REAL*8 zf
c      REAL*8 sqt,theta

      REAL*8 BB1, BB2
      REAL*8 UDASINH
      REAL*8 FTN
      REAL*8 RFTN1,RFTN2
      REAL*8 RR
      REAL*8 awx
      REAL*8 z,x
      REAL*8 AZ11,AZ12,AZ2,AZ31,AZ32,AZ4,AZ5,AZ6
      REAL*8 sumAZ11,sumAZ12,sumAZ2
      REAL*8 sumAZ31,sumAZ32,sumAZ4,sumAZ5,sumAZ6
      REAL*8 A(6)
      REAL*8 Dpole(6)
      REAL*8 dum
      INTEGER i,j
c
c
      REAL*8 pi2
      REAL*8 mn
      REAL*8 rpin
      REAL*8 rdelpi
      REAL*8 omegas
      REAL*8 ts
      REAL*8 Afac
      REAL*8 Dfac
      real*8 Rparam
c      
c----------------------------------------------------------------------
c
      PARAMETER (pi2= pi*pi, mn=Mnucleon)
      PARAMETER (Afac   = gpind**2*mpi/fpi**2)
      PARAMETER (Dfac   = b1**2/(mn**2*18.0d0))
      PARAMETER (Rparam = delta/mpi + dsqrt((delta/mpi)**2 - 1))
c
c----------------------------------------------------------------------
c
c   Auxilliary functions
c
c----------------------------------------------------------------------
c
c
      FTN(x,z)=ts*(1.0d0-x)*(1.0d0-z)*z
      RFTN1(x,z)= -1.0d0 + (rdelpi-omegas*x*(1.0d0-z))**2 /
     &                            (1.0d0-FTN(x,z))
      RFTN2(x,z)= -1.0d0 + (rdelpi+omegas*x*(1.0d0-z))**2 /
     &                            (1.0d0-FTN(x,z))
      UDASINH(x)=DLOG(x+DSQRT(1.0d0+x**2))
      RR(dum)=dum+DSQRT(dum**2-1.0d0)
      awx(omegas,x)=DSQRT( (rdelpi-omegas*x)**2 -1.0d0)
      
      BB1(x,z)=DSQRT((rdelpi-omegas*x*(1.0d0-z)) * 
     &               (rdelpi-omegas*x*(1.0d0-z)) 
     &                  / (1.0d0-FTN(x,z)) - 1.0d0)
     
      BB2(x,z)=DSQRT((rdelpi+omegas*x*(1.0d0-z)) *
     &               (rdelpi+omegas*x*(1.0d0-z)) 
     &                  / (1-FTN(x,z)) - 1.0d0)
     
c
c----------------------------------------------------------------------
c
c   First set functions
c
c----------------------------------------------------------------------
c
c
 
      AZ11(x) =  awx( omegas,1.0d0)*UDASINH(awx( omegas,1.0d0)) 
     &          + awx(-omegas,1.0d0)*UDASINH(awx(-omegas,1.0d0))
     &   -2.0d0*( awx( omegas,x)*UDASINH(awx( omegas,x))
     &          + awx(-omegas,x)*UDASINH(awx(-omegas,x)) )
 
      AZ12(x,z)=
     & (1.0d0-z)/DSQRT(1.0d0 - FTN(x,z) ) * (
     &   ((0.5d0*ts*x*(1.0d0-z)+(omegas*x*(1.0d0-z))**2 +FTN(x,z)+
     &    5.0d0 * (
     &      (rdelpi-omegas*x*(1.0d0-z))**2 + FTN(x,z) -1.0d0
     &    )) *
     &    UDASINH( DSQRT(RFTN1(x,z)) ) / DSQRT(RFTN1(x,z))
     &    ) +
     &   ((0.5d0*ts*x*(1.0d0-z)+(omegas*x*(1.0d0-z))**2 +FTN(x,z)+
     &    5.0d0 * (
     &      (rdelpi+omegas*x*(1.0d0-z))**2 + FTN(x,z) -1.0d0
     &    )) *
     &    UDASINH( DSQRT(RFTN2(x,z)) ) / DSQRT(RFTN2(x,z))
     &    ) ) +
     & (3.0d0-3.0d0*rdelpi**2-4.0d0*ts*(1.0d0-x)*x) *
     &  DLOG( (rdelpi+DSQRT(-1.0d0+rdelpi**2+ts*(1.0d0-x)*x))
     &                  /DSQRT( 1.0d0-ts*(1.0d0-x)*x) ) /
     &  DSQRT(-1.0d0+rdelpi**2+ts*(1.0d0-x)*x) +
     &  rdelpi*(1.0d0-z)*(
     &    2.0d0*DLOG(1.0d0-FTN(x,z))+3.0d0*DLOG((1.0d0-FTN(x,z))/
     &  (1.0d0-ts*(1.0d0-x)*x) )
     &  )

      AZ2(x,z)= (1.0d0-z)* (
     &  ( (UDASINH( DSQRT(RFTN1(x,z)) ) / DSQRT(RFTN1(x,z))) +
     &         (UDASINH( DSQRT(RFTN2(x,z)) ) / DSQRT(RFTN2(x,z))) )*
     &          ((1.0d0-x)*(1.0d0-7.0d0*z)*(1.0d0-z)+z) /
     &            DSQRT(1.0d0-FTN(x,z))
     &       +
     &     (((-DSQRT((RFTN1(x,z)+1.0d0)) + 
     &       UDASINH( DSQRT(RFTN1(x,z))) / DSQRT(RFTN1(x,z)) )/ 
     &         RFTN1(x,z) ) +
     &      ((-DSQRT((RFTN2(x,z)+1.0d0)) + 
     &       UDASINH( DSQRT(RFTN2(x,z))) / DSQRT(RFTN2(x,z)) )/ 
     &         RFTN2(x,z) )) *
     &      (1.0d0-x)*(1.0d0-z)**2*z* 
     &        (0.5d0*ts*x+(omegas*x)**2*(1.0d0-z)+ts*(1.0d0-x)*z) /
     &       (1.0d0-FTN(x,z))**(1.5d0)
     &       )


      AZ31(x,z)=
     & -omegas**4*(1.0d0-(1.0d0+ts/(2.0d0*omegas**2))**2) *
     &             (1.0d0-x)*x*(1.0d0-z)**3*z *
     &  (( DSQRT(RFTN1(x,z)+1.0d0) - UDASINH(DSQRT(RFTN1(x,z))) / 
     &                                     DSQRT(RFTN1(x,z)) ) /
     &    RFTN1(x,z) - 
     &  ( DSQRT(RFTN2(x,z)+1.0d0) - UDASINH(DSQRT(RFTN2(x,z))) / 
     &                                     DSQRT(RFTN2(x,z)) ) /
     &    RFTN2(x,z) ) / (1.0d0-FTN(x,z))**(1.5d0)

      AZ32(x)=2.0d0 * (
     &   DSQRT( (rdelpi-omegas*x)**2-1.0d0 ) *
     &         UDASINH(DSQRT( (rdelpi-omegas*x)**2-1.0d0 )) -
     &   DSQRT( (rdelpi+omegas*x)**2-1.0d0 ) *
     &         UDASINH(DSQRT( (rdelpi+omegas*x)**2-1.0d0 ))  )

      AZ4(x,z)=x*(1-z)*(1-z)/DSQRT(1-ts*z*(1-z)*(1-x)) * 
     &      (UDASINH(BB1(x,z))/BB1(x,z) - UDASINH(BB2(x,z))/BB2(x,z)) 
     
     
      AZ5(x,z)= (1.0d0-z)*z* (
     &          ((UDASINH(DSQRT(RFTN1(x,z))) / DSQRT(RFTN1(x,z)) -
     &            UDASINH(DSQRT(RFTN2(x,z))) / DSQRT(RFTN2(x,z)) ) /
     &               DSQRT(1.0d0-FTN(x,z)) ) -
     &     (((DSQRT(RFTN1(x,z)+1.0d0) -
     &             UDASINH(DSQRT(RFTN1(x,z))) / DSQRT(RFTN1(x,z)) ) /
     &                RFTN1(x,z) )-
     &       (DSQRT(RFTN2(x,z)+1.0d0) -
     &             UDASINH(DSQRT(RFTN2(x,z))) / DSQRT(RFTN2(x,z))) /
     &                 RFTN2(x,z) ) *
     &     omegas**2*(1.0d0+ts/(2.0d0*omegas**2))*
     &                      (1.0d0-x)*x*(1.0d0-z)**2/ 
     &    (1.0d0-FTN(x,z))**(1.5d0)
     &     )
     
      AZ6(x,z)= (1.0d0-z)*z* (
     &          ((UDASINH(DSQRT(RFTN2(x,z))) / DSQRT(RFTN2(x,z)) -
     &            UDASINH(DSQRT(RFTN1(x,z))) / DSQRT(RFTN1(x,z)) ) /
     &               DSQRT(1.0d0-FTN(x,z)) ) + 
     &     (((DSQRT(RFTN1(x,z)+1.0d0) -
     &             UDASINH(DSQRT(RFTN1(x,z))) / DSQRT(RFTN1(x,z)) ) /
     &                RFTN1(x,z) ) -
     &       (DSQRT(RFTN2(x,z)+1.0d0) -
     &             UDASINH(DSQRT(RFTN2(x,z))) / DSQRT(RFTN2(x,z))) /
     &                 RFTN2(x,z) ) *
     &     omegas**2*(1.0d0-x)*x*(1.0d0-z)**2/ 
     &    (1.0d0-FTN(x,z))**(1.5d0)
     &     )
c----------------------------------------------------------------------
c
c   Some useful quantities
c
c----------------------------------------------------------------------

c      zf=1.0d0+t/(2.0d0*omega**2)
c      sqt     =  dsqrt((-1.0d0*t)) 
c      theta   =  dacos(1.0d0 + t/(2.0d0*omega**2))

      omegas =  omega/mpi
      ts     =  t/mpi**2
      rpin   =  mpi/mn
      rdelpi =  delta/mpi
c----------------------------------------------------------------------
c
c   Now calculate integrals
c
c----------------------------------------------------------------------
c
c   A1 to A6 are calculated using scaled omega and t, 
c   omegas=omega/mpi and ts=t/mpi**2
c
      
      sumAZ11 = 0.0d0
      sumAZ12 = 0.0d0
      sumAZ2   = 0.0d0
      sumAZ31 = 0.0d0
      sumAZ32 = 0.0d0
      sumAZ4   = 0.0d0
      sumAZ5   = 0.0d0
      sumAZ6   = 0.0d0
     
      do i=1,Nx
        sumAZ11 = sumAZ11 + wx(i)*AZ11(xq(i))
        sumAZ32 = sumAZ32 + wx(i)*AZ32(xq(i))	
        do j=1,Nx
	  sumAZ12 = sumAZ12 + wx(i)*wx(j)*AZ12(xq(i),xq(j))
	  sumAZ2   = sumAZ2   + wx(i)*wx(j)*AZ2(xq(i),xq(j))
	  sumAZ31 = sumAZ31 + wx(i)*wx(j)*AZ31(xq(i),xq(j))	  
	  sumAZ4   = sumAZ4   + wx(i)*wx(j)*AZ4(xq(i),xq(j))
	  sumAZ5   = sumAZ5   + wx(i)*wx(j)*AZ5(xq(i),xq(j))
	  sumAZ6   = sumAZ6   + wx(i)*wx(j)*AZ6(xq(i),xq(j))
	end do
      end do


      A(1)=(sumAZ11+sumAZ12)*2.0d0/(9.0d0*pi2)
      A(2)=-sumAZ2*2.0d0*omegas**2/(9.0d0*pi2)
      A(3)=(sumAZ31+sumAZ32 - DSQRT( (rdelpi-omegas)**2 -1.0d0) *
     &                        DLOG(RR(rdelpi-omegas))
     &                      + DSQRT( (rdelpi+omegas)**2 -1.0d0) *
     &                        DLOG(RR(rdelpi+omegas)) ) /
     &                 (9.0d0*pi2)
      A(4)=-sumAZ4*omegas**2/(9.0d0*pi2)
      A(5)= sumAZ5*omegas**2/(9.0d0*pi2)
      A(6)= sumAZ6*omegas**2/(9.0d0*pi2)     
      
c
c   Dpole is calculated in terms of unscaled omega and t
c
      Dpole(1)=-2.0d0*delta*(2.0d0*omega**2+t) / 
     &                 (omega**2-delta**2)
      Dpole(2)=4.0d0*delta*omega**2/(omega**2-delta**2)
      Dpole(3)=(2.0d0*omega**2+t)*omega/(omega**2-delta**2)
      Dpole(4)=2.0d0*omega**3/(omega**2-delta**2)
      Dpole(5)=-Dpole(4)
      Dpole(6)=0.d0 
      
      A1=(Afac*A(1) + Dfac*Dpole(1))
      A2=(Afac*A(2) + Dfac*Dpole(2))/omega**2
      A3=(Afac*A(3) + Dfac*Dpole(3))
      A4=(Afac*A(4) + Dfac*Dpole(4))/omega**2
      A5=(Afac*A(5) + Dfac*Dpole(5))/omega**2
      A6=(Afac*A(6) + Dfac*Dpole(6))/omega**2

c    hgrie Nov 2014: following used to be in constructamps*.f, but now here:
      
c    Subtract static values of Delta+DeltaPi diagrams from A1 and A2,
c               as given in Eq.(3.20) of Hildebrandt Thesis.
c    This (usually) sets the static values of the _scalar dipole polarisabilities_
c         to the order Q^3-values, as derived by BKM.
c         Spin-polarisabilities remain as predicted in epsilon^3.
c         
c    REASON: epsilon^3 delta-contributions give huge positive shift to alpha and beta,
c         while exp numbers much closer to Q^3 values.
c    commented out: implementation when static scalar pols set to some specific value   
c   
         delalpha=-(gpind**2/(54*fpi**2*pi**2)*(9*delta/(delta**2-mpi**2)+
     &             (delta**2-10*mpi**2)/((delta**2-mpi**2)*dsqrt(delta**2-mpi**2))*
     &             dlog((delta+dsqrt(delta**2-mpi**2))/mpi)))
         delbeta=-(2*b1**2/(9*delta*Mnucleon**2) + 
     &          gpind**2/(54*fpi**2*pi**2)*1/(dsqrt(delta**2-mpi**2))*
     &                    dlog((delta+dsqrt(delta**2-mpi**2))/mpi))
c         
c    NOW add these offsets to pols
      
         A1=A1+(delalpha+delbeta*(1.0d0 + t/(2.0d0*omega**2)))*omega**2
         A2=A2-delbeta
      
      if (verbosity.eq.2) then
           write (*,*) "Single-Nucleon Amplitude of Delta subroutine (iso-scalar):"
           write (*,*) "   [pre-boost, basis as in review (2.1), no factor e^2]"
           write (*,'(A,E15.7,A)') "   A1 = ",A1," MeV^-1"
           write (*,'(A,E15.7,A)') "   A2 = ",A2*omega**2," MeV^-1"
           write (*,'(A,E15.7,A)') "   A3 = ",A3," MeV^-1"
           write (*,'(A,E15.7,A)') "   A4 = ",A4*omega**2," MeV^-1"
           write (*,'(A,E15.7,A)') "   A5 = ",A5*omega**2," MeV^-1"
           write (*,'(A,E15.7,A)') "   A6 = ",A6*omega**2," MeV^-1"
           write (*,*) "Static value of Delta and DeltaPi contributions to alphaE1 and betaM1:"
           write (*,'(A,F5.2,A)') "   alphaE1DeltaDeltaPi = ",-delalpha*HC**3*10000.0*alphaem," x 10^-4 fm^3"
           write (*,'(A,F5.2,A)') "   betaM1DeltaDeltaPi  = ",-delbeta*HC**3*10000.0*alphaem," x 10^-4 fm^3"
           write (*,*) "**Subtracted static value of Delta and DeltaPi contributions to alphaE1 and betaM1.**"
           write (*,*) "Contribution to static polarisabilities from Delta subroutine (iso-scalar):"
           write (*,'(A,F5.2,A)') "   alphaE1DeltaDeltapi    = ",0.," x 10^-4 fm^3 [DeltaPi part subtracted]"
           write (*,'(A,F5.2,A)') "   betaM1DeltaDeltapi     = ",0.," x 10^-4 fm^3 [Delta+DeltaPi parts subtracted]"
           write (*,'(A,F5.2,A)') "   gammaE1E1DeltaDeltapi  = ",
     &          (alphaem*gpiND**2*((delta**2 + 5*mpi**2)/(delta**2 - mpi**2)**2 + 
     &          (delta*(delta**2 - 7*mpi**2)*Log(Rparam))/(delta**2 - mpi**2)**2.5))/(108.*fpi**2*Pi**2)*10000*(HC)**4,
     &          " x 10^-4 fm^4"
           write (*,'(A,F5.2,A)') "   gammaM1M1DeltaDeltapi  = ",
     &          ((alphaem*b1**2)/(9*delta**2*Mnucleon**2) + ! Delta pole only
     &          (alphaem*gpiND**2*(-(1/(delta**2 - mpi**2)) + (delta*dlog(Rparam))/(delta**2 - mpi**2)**1.5))/
     &          (108.*fpi**2*Pi**2))*10000*(HC)**4," x 10^-4 fm^4"
           write (*,'(A,F5.2,A)') "   gammaE1M2DeltaDeltapi  = ",
     &          -((alphaem*gpiND**2*(-(1/(delta**2 - mpi**2)) + (delta*dlog(Rparam))/(delta**2 - mpi**2)**1.5))/
     &          (108.*fpi**2*Pi**2))*10000*(HC)**4," x 10^-4 fm^4"
           write (*,'(A,F5.2,A)') "   gammaM1E2DeltaDeltapi  = ",
     &          -((alphaem*gpiND**2*(-(1/(delta**2 - mpi**2)) + (delta*dlog(Rparam))/(delta**2 - mpi**2)**1.5))/
     &          (108.*fpi**2*Pi**2))*10000*(HC)**4," x 10^-4 fm^4"
      end if   
      RETURN
      END
c end subroutine CalculateAsDelta
c      
c     DP, 1/19: This next portion is a stub that could be pulled
c     out and used to form a main program with which the routine could be tested
c      
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
c    Beginning of stub       
c         PROGRAM MAIN
c   c         
c         IMPLICIT NONE
c         
c         REAL*8 A1,A2,A3,A4,A5,A6
c         REAL*8 t,omega,theta
c         REAL*8 pi,mpi
c         REAL*8 wx(7000), xq(7000)
c         
c         INTEGER Nx
c         INTEGER i, j
c         
c         pi      = 4.0d0*datan(1.0d0)
c         mpi     = 139.6d0
c         theta   = 140.0d0*pi/180.
c         omega   = 0.8d0*mpi
c         t       = 2.0d0*omega**2*(dcos(theta)-1.0d0)
c   
c   c  Define number of steps used to perform numerical integration
c   c  with Nx=200, results appear to converge to 4dp (for omega<mpi)
c   c
c         Nx      = 1000
c         
c         write(6,*) "omega=",omega," t=",t," pi=",pi
c   
c         write(*,*)"Nx",Nx
c   
c         do i=1,Nx
c           wx(i)=1.0d0/(Nx-1)
c           xq(i)=(i-1.0d0)/(Nx-1)
c         end do
c         wx(1)=wx(1)/2.0d0
c         wx(Nx)=wx(Nx)/2.0d0
c   
c   c         
c         call CalculateAsDelta(A1,A2,A3,A4,A5,A6,wx,xq,Nx,t,omega)
c         
c         write(*,*)
c         write(*,*)"================================"
c         write(*,10)"A1 is ",A1
c         write(*,10)"A2 is ",A2
c         write(*,10)"A3 is ",A3
c         write(*,10)"A4 is ",A4
c         write(*,10)"A5 is ",A5
c         write(*,10)"A6 is ",A6
c         write(*,*)"================================"
c   
c   
c    10   FORMAT(A7,E15.8)
c         stop
c         END
c      End of stub
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
c   
c        checked Feb 2012 hgrie
c   

