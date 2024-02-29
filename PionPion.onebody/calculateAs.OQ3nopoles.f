c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
      subroutine CalculateAsOQ3nopoles(A1,A2,A3,A4,A5,A6,
     &     xq,wx,Nx,t,omega,thetacm,kappa,q,
     &     Dalpha,Dbeta,DgE1E1,DgM1M1,DgE1M2,DgM1E2,verbosity)
c
c Initial routine header:
c   REAL COMPTON SCATTERING (NOVEMBER 1997, MANUEL MALHEIRO)
c CALCULATION OF THE INVARIANT A1,...,A6 FUNCTIONS OF THE T- MATRIX
c     IN CHIRAL PERTURBATION THEORY AT O(q³), FOR GIVEN PHOTON
c     ENERGY AND SCATTERING ANGLE AND GIVEN KAPPA AND CHARGE.
c
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c      
c     hgrie May 2022: set nucleon masses to Mnucleon in A[2-6]
c      
c     Version 7 hgrie Feb 2013 -- compatibility to previous versions broken
c         consult README for explanation of variables, kinematics,
c                            version history, original authors etc.       
c   
c     Split by hgrie Feb 2012: Used to be part of calculateAs.OQ4.delta.nopoles.f
c         
c
c Checked 7/98-D.P. checked Feb 2012 hgrie
c         
c    CLARIFICATION by hgrie Nov 2014:
c         calculates the "nopoles" part of Compton amplitudes A1 to A6
c         at order e²δ² = Q³:
c         1-loop results of [BKM] in "strict" HBChiPT,
c             no shift for pion-production threshold.
c         As described in [HildeDiplom (5.10/11)] and [review], this is 
c         NOT ONLY a polarisability contribution but also contains
c         non-structure parts.
c         The polarisabilities parts are thus the 1-loop results of [BKM].
c         
c         They result in the static polarisabilities 
c         which are coded at the end of this file [BKM, review (4.13)].
c         
c         alphaE1 = 10 betaM1 =  10 alphaem ga**2/(192*Pi*mpi*fpi**2)
c         
c         gammaE1E1 = 5 gammaM1M1 = -5 gammaE1M2 = -5 gammaM1E2 =
c                   = -5*alphaem*ga**2/(96*Pi**2*mpi**2*fpi**2)
c         
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c         
c**********************************************************************
c
      implicit none
c
c**********************************************************************
c
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
c      include '../common-densities/calctype.def'
c
c**********************************************************************
c
c  Output variables:
c
c     A1,A2,A3,A4,A5,A6-six invariant functions for gamma-N
c       scattering, calculated as per BKM review, p.66.
c
      real*8 A1,A2,A3,A4,A5,A6
c
c**********************************************************************
c
c  Input variables:
c
c     t,omega-Mandelstam t and c.m. energy of photon scattering
c     thetacm-c.m. photon scattering angle
c     xq,wx,Nx-define quadratures for Feynman-parameter integration
c     kappa,q-anomalous magnetic moment and charge
c     Dalpha,Dbeta,DgE1E1,DgM1M1,DgE1M2,DgM1E2-Shifts in static
c     dipole polarizabilities, see above
c     verbosity-if set to 2 prints out results as we go along
c      
      real*8 t,omega,xq(Nxmax),wx(Nxmax)
      integer Nx
      real*8 kappa,q,thetacm
      real*8 Dalpha,Dbeta,DgE1E1,DgM1M1,DgE1M2,DgM1E2
      integer verbosity
c
c**********************************************************************
c
c  Local variables:
c
c     R,W-convenient functions
c     z,x-arguments of R,W
c     theta-scattering angle
c     sqt=sqrt(-t)
c     AZ1..AZ6-integrands
c     sumAZ1..sumAZ6-results of integrals
c     i,j-loop variables
c
      real*8 sqt,theta
      real*8 R,W,z,x
      real*8 AZ1,AZ2,AZ3,AZ4,AZ5,AZ6
      real*8 sumAZ1,sumAZ2,sumAZ3,sumAZ4,sumAZ5,sumAZ6
      integer i,j
c----------------------------------------------------------------------
c
c   First set functions
c
c----------------------------------------------------------------------

      W(x,z) = dsqrt(mpi2 - (omega*z)**2. + t*(1.0d0 - z)**2.*x*(x - 1.0d0))
      R(x,z) = dsqrt(mpi2 + t*(1.0d0-z)**2.*x*(x - 1.0d0))
      AZ1(z) = datan((1.0d0 - z)*sqt/
     &	                (2.0d0*dsqrt(mpi2 - (omega*z)**2.))) 
      AZ2(z) = 2.0d0*(1-z)*dsqrt(t*((omega*z)**2.-mpi2))/
     &      (4.0d0*(mpi2 - (omega*z)**2.) - t*(1.d0 - z)**2.) 
      AZ3(x,z) = x*(1.0d0-x)*z*(1-z)**3./(W(x,z)**3.)*(dasin(omega*z/R(x,z))
     &    + omega*z*W(x,z)/(R(x,z)**2.))
      AZ4(x,z) = z*(1-z)/W(x,z)*dasin(omega*z/R(x,z))
      AZ5(x,z) = (1.0d0-z)**2/W(x,z)*dasin(omega*z/R(x,z))
c----------------------------------------------------------------------
      sqt = dsqrt((-1.0d0*t)) 
      theta=thetacm !dacos(1.0 + t/(2.0*omega**2))
c----------------------------------------------------------------------
c
c   Now calculate integrals
c
c----------------------------------------------------------------------
      sumAZ1 = 0.d0
      sumAZ2 = 0.d0
      sumAZ3 = 0.d0
      sumAZ4 = 0.d0
      sumAZ5 = 0.d0
      do i=1,Nx         
         sumAZ1 = sumAZ1 + wx(i)*AZ1(xq(i))
         sumAZ2 = sumAZ2 + wx(i)*AZ2(xq(i))
         do j=1,Nx
            sumAZ3 = sumAZ3 + wx(i)*wx(j)*AZ3(xq(i),xq(j))
            sumAZ4 = sumAZ4 + wx(i)*wx(j)*AZ4(xq(i),xq(j))
            sumAZ5 = sumAZ5 + wx(i)*wx(j)*AZ5(xq(i),xq(j))
         end do
      enddo
c     add those piece which need no integration
c A1 INVARIANT FUNCTION -- -q**2/M  already in OQ2 amplitude!
      A1 = + gafac*(mpi - dsqrt(mpi2 - omega**2)
     &   +(2.0d0*mpi2 - t)/sqt * (datan(sqt/(2.d0*mpi))/2.0d0 - sumAZ1))
c A2 INVARIANT FUNCTION 
      A2 = gafac*(t - 2.0d0*mpi2)/(sqt**3)*(sumAZ1 - sumAZ2) 
c     &     +q**2/(Mnucl**2*omega)   ! Thomson boost, now included in ??????????????????????????????
c     Note (DP 1/19): This, and the last line of each subsequent A, are now taken
c     care of in calculateAs.OQ3poles.f
c A3 INVARIANT FUNCTION
      A3 = omega/(2.0d0*Mnucleon**2.)*(q + 2.d0*kappa)*q ! spin-orbit seagull
c     &	 -  (q + kappa)**2*dcos(theta)) ! double mag moment
     &   + (2.d0*q - 1.d0)*ga*t*omega/(8.0d0*fpi2*pi**2.*(mpi02 - t)) ! pi0 pole
     &   + gafac/pi*(mpi2/omega*(dasin(omega/mpi)**2.) - omega)
     &   + 2.0d0*gafac/pi*(omega**4.)*(dsin(theta)**2.)*sumAZ3
c A4 INVARIANT FUNCTION
      A4 = 2.0d0*gafac/pi*sumAZ4
c     &   - (q+kappa)**2/(2.d0*omega*Mnucleon**2) !  double mag moment
c A5 INVARIANT FUNCTION
      A5 = - (2.d0*q -1.d0)*ga*omega/(8.0d0*fpi2*pi**2.*(mpi02 - t)) ! pi0 pole
     &   + gafac/pi*(-sumAZ5 + 2.d0*omega**2*dcos(theta)*sumAZ3)
c     &   + (q + kappa)**2/(2.d0*omega*Mnucleon**2.) !  double mag moment
c A6 INVARIANT FUNCTION 
      A6 =  (2.d0*q - 1.d0)*ga*omega/(8.0d0*fpi2*pi**2.*(mpi02 - t)) ! pi0 pole
     &   + gafac/pi*(sumAZ5 - 2.d0*omega**2*sumAZ3)
c     &   - q*(q + kappa)/(2.0d0*omega*Mnucleon**2.) ! 1 min subst + 1 mag moment in ant-diagram
      
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c following MODIFIED by Deepshikha 05/08 to include shifts in polarizabilities
c         re-modified by hgrie to correct signs of DgE1E1 in A3.
c         consistent with Review 2012 and ragusa definition and Babusci and HGHP and...
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      A1=A1+(Dalpha+Dbeta*dcos(thetacm))*omega**2
      A2=A2-Dbeta
      A3=A3-(DgE1E1+DgE1M2+(DgM1E2+DgM1M1)*dcos(thetacm))*omega**3
      A4=A4+(DgM1E2-DgM1M1)*omega
      A5=A5+DgM1M1*omega
      A6=A6+DgE1M2*omega
c
      if (verbosity.eq.2) then
           write (*,'(A,I1,A)') " Single-Nucleon Amplitudes of OQ3nopoles subroutine for nucleon of charge ",nint(q),":"
           write (*,*) "   [pre-boost, basis as in review (2.1), no factor e²]"
           write (*,'(A,E15.7,A)') "   A1 = ",A1," MeV¯¹"
           write (*,'(A,E15.7,A)') "   A2 = ",A2*omega**2," MeV¯¹"
           write (*,'(A,E15.7,A)') "   A3 = ",A3," MeV¯¹"
           write (*,'(A,E15.7,A)') "   A4 = ",A4*omega**2," MeV¯¹"
           write (*,'(A,E15.7,A)') "   A5 = ",A5*omega**2," MeV¯¹"
           write (*,'(A,E15.7,A)') "   A6 = ",A6*omega**2," MeV¯¹"
           write (*,'(A,I1)') "Contribution to static polarisabilities of nucleon with charge ",nint(q)
           write (*,'(A)')    "             from OQ3nopoles subroutine (iso-scalar):"
           write (*,'(A,F5.2,A)') "   alphaE1Npi    = ",10*alphaem*ga**2/(192* pi* fpi**2* mpi)*10000*(HC)**3," x 10¯⁴ fm³"
           write (*,'(A,F5.2,A)') "   betaM1Npi     = ",alphaem*ga**2/(192* pi* fpi**2* mpi)*10000*(HC)**3," x 10¯⁴ fm³"
           write (*,'(A,F5.2,A)') "   gammaE1E1Npi  = ",-20*alphaem*ga**2/(384 *pi**2 *fpi**2 *mpi**2 )*10000*(HC)**4,
     &          " x 10¯⁴ fm⁴"
           write (*,'(A,F5.2,A)') "   gammaM1M1Npi  = ",-4*alphaem*ga**2/(384 *pi**2 *fpi**2 *mpi**2 )*10000*(HC)**4,
     &          " x 10¯⁴ fm⁴"
           write (*,'(A,F5.2,A)') "   gammaE1M2Npi  = ",4*alphaem*ga**2/(384 *pi**2 *fpi**2 *mpi**2 )*10000*(HC)**4,
     &          " x 10¯⁴ fm⁴"
           write (*,'(A,F5.2,A)') "   gammaM1E2Npi  = ",4*alphaem*ga**2/(384 *pi**2 *fpi**2 *mpi**2 )*10000*(HC)**4,
     &          " x 10¯⁴ fm⁴"
      end if   
      if (AZ6.eq.1000) continue ! unused variable kept for future use
      if (sumAZ6.eq.1000) continue ! unused variable kept for future use
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c END Modification 05/08
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
      end
c end subroutine CalculateAsOQ3nopoles
