c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
      subroutine CalculateAsOQ3poles(A1,A2,A3,A4,A5,A6,
     &     omega,thetacm,kappa,q,massnumber,verbosity)
c
c Initial routine header:
c     REAL COMPTON SCATTERING (NOVEMBER 1997, MANUEL MALHEIRO)
c     CALCULATION OF THE INVARIANT A1,...,A6 FUNCTIONS OF THE T- MATRIX
c     IN CHIRAL PERTURBATION THEORY AT O(q^3), FOR GIVEN PHOTON
c     ENERGY AND SCATTERING ANGLE AND GIVEN KAPPA AND CHARGE.
c
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$     
c     hgrie May 2022: set nucleon masses to Mnucleon in A[2-6]
c      
c     Version 7 hgrie Feb 2013 -- compatibility to previous versions broken
c         consult README for explanation of variables, kinematics,
c                            version history, original authors etc.       
c   
c     Split by hgrie Feb 2012: Used to be part of calculateAs.OQ4.delta.f
c     AMPLITUDES which are pole in Born and must NOT be included when rescattering
c            is implemented for correct Thomson limit.
c         
c Checked 7/98-D.P. checked Feb 2012 hgrie
c    following expressions sum with the expression of calculateAs.OQ3nopoles.f
c    indeed to the total BKM amplitudes --
c         see also Hildebrandt PhD App. B. with explanation of contribs
c         on page following (B.1): "true" pole contributions are only these
c         associated with s- & u-channel diagrams Fig. B.1(a+b). 
c    hgrie Feb 2012, supplemented hgrie Oct 2014 for clarity
c    hgrie Nov 2014: corrected boost correction to A2 and added variable
c         "massnumber" = target nucleus mass in units of nucleon mass.
c         boost correction for A2 depends on target mass
c           see notes Daniel boost-corrections.phillips20141028.pdf
c           and [Beane2005 eq. (D.6)]
c         following email discussions, this is implemented as pole-contribution
c         see email 29 Oct 2014 from hgrie to Judith, Daniel, Bruno. 
c
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$     
c
c**********************************************************************
c
      implicit none
c
c**********************************************************************
c
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
c
c**********************************************************************
c
c  Output variables:
c
c     A1,A2,A3,A4,A5,A6-six invariant functions for gamma-N
c
      real*8 A1,A2,A3,A4,A5,A6
c
c**********************************************************************
c
c  Input variables:
c
      real*8 omega
      real*8 kappa,q,thetacm
      real*8 massnumber ! target nucleus mass in units of nucleon mass
      integer verbosity
c
c**********************************************************************
c
c  Local variables:
c
c     theta-scattering angle
c
      real*8 theta
c
      theta=thetacm !dacos(1.0 + t/(2.0*omega**2))

      A1 = 0
c    hgrie Nov 2014: added target mass in A2 boost correction, described above
      A2 = q**2/(massnumber*Mnucleon**2.*omega)
c     rest implemented following [Beane2005] and [HildePhD App. B]:
c     not "boost" contributions, so not factor "massnumber".
      A3 = -omega/(2.0d0*Mnucleon**2.)*(q + kappa)**2*dcos(theta) ! 2 mag moment
      A4 = -(q+kappa)**2/(2.d0*omega*Mnucleon**2.) ! 2 mag moment
      A5 = +(q+kappa)**2/(2.d0*omega*Mnucleon**2.) ! 2 mag moment
      A6 = -q*(q +kappa)/(2.d0*omega*Mnucleon**2.) ! 1 min subst, 1 mag moment
   
      if (verbosity.eq.2) then
           write (*,*) "Single-Nucleon Amplitudes of OQ3poles subroutine for nucleon with charge ",q,":"
           write (*,*) "   [pre-boost, basis as in review (2.1), no factor e^2]"
           write (*,'(A,E15.7,A)') "   A1 = ",A1," MeV^-1"
           write (*,'(A,E15.7,A)') "   A2 = ",A2*omega**2," MeV^-1"
           write (*,'(A,E15.7,A)') "   A3 = ",A3," MeV^-1"
           write (*,'(A,E15.7,A)') "   A4 = ",A4*omega**2," MeV^-1"
           write (*,'(A,E15.7,A)') "   A5 = ",A5*omega**2," MeV^-1"
           write (*,'(A,E15.7,A)') "   A6 = ",A6*omega**2," MeV^-1"
           write (*,*) " NO contribution to static polarisabilities from OQ3poles subroutine."
      end if   
      return
      end
c end subroutine CalculateAsOQ3poles
