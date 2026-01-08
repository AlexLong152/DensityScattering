cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of KERNEL code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              CalculateQs      : compute vectors needed by kernel/diagrams
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, loosely based on subroutine which was part of common-densities/calcmomenta.f of Compton density code v2.0 hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     specific to process considered
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine calculateqsmass(kVec,kpVec,thetacm,mPion,mNucl,verbosity)
c     Kinematics for pion photoproduction
c     Derivation for the kinematics can be found in pionpionAngle.pdf
c     OneDrive/thesis/Kinematics/pionpionAngle.pdf
c     The conversion from lenkewitz to our variables is found in
c     DensityScattering/documentation/pionpionangle.pdf

      implicit none
c
c**********************************************************************
c
c  Input variables:
c
      real*8, intent(in) :: kVec(3),mPion,mNucl, thetacm
      integer, intent(in) ::verbosity 
c  Output variables
      real*8, intent(out) :: kpVec(3)
c**********************************************************************
c
c     temporary variables
    
      real*8 Epion, mandalS, ENuc, kpsq, kpAbs, omegaThreshold,k
c**********************************************************************      
c
c     
c     s =(p+k)^2=[(E_nuc, 0,0,-omega) + (omega, 0,0,omega)]^2
c     s = (E_nuc+omega)^2
c     mpi and mpi0 are pion mass in MeV defined in ../common-densities/constants.def
c     Internal Variables first     

c     VARIABLE DESCRIPTIONS
c     -----------------------------------------------
c     q: propgator for diagram A, note q=q_2, the second propogator for diagram B
c     q1:First propogatior for diagram B, q1=q-k

c     A full derivation of the kinematics can be found in the particle physics bookley (pdg)
c     https://pdg.lbl.gov/2023/download/db2022.pdf pg 259, see figure 49.6
      k=sqrt(DOT_PRODUCT(kVec,kVec))
      omegaThreshold=(mPion*(mPion+2*mNucl))/(2*(mPion+mNucl))
      if (abs(k-omegaThreshold).le.2) then
c             assume that if within 2MeV we are trying to run at threshold
c             write(*,*) "reset k"
          k=omegaThreshold
      else
        write(*,*) "kp^2<0 -> given masses/energy are incompatable"
        stop
      end if
      ENuc=sqrt((mNucl**2) + (k**2))
      mandalS=(ENuc + k)**2 !lab frame
      Epion=(mandalS+(mPion**2)-(mNucl**2))/(2*sqrt(mandalS))
      kpsq=(((mandalS+mPion**2-mNucl**2)**2)/(4*mandalS))-mPion**2
      kpAbs=sqrt(kpsq)
      kpVec=(/0.d0,kpAbs*sin(thetacm), kpAbs*cos(thetacm)/)!without loss of generality - phi undetermined

      if (abs(Epion-sqrt(mPion**2+kpsq)).ge.1) then
        write(*,*) "something wrong in kinematics"
        stop
      end if

      if (kpsq.lt.-10) then
          write(*,*) "kp^2<0 -> given masses/energy are incompatable"
          write(*,*) "This error should have been caught already.... something went wrong"
          stop
      end if


c     q = (p-pp)+((kVec+kp)/2)
c     q1 = (p-pp)+((kp-kVec)/2)
c     q1 = q-k

c     write(*,*) "#################################################################################"
c     write(*,*) "In calcmomenta.f: Epion=",Epion 
c     write(*,*) "In calcmomenta.f: mandalS=",mandalS 
c     write(*,*) "In calcmomenta.f: kpsq=",kpsq 
c     write(*,*) "kp=",kp 
c     write(*,*) "#################################################################################"
c     write(*,*) "In calcmomenta kpAbs=", kpAbs
c     write(*,*) "In calcmomenta q=", q
c     write(*,*) ""
c     write(*,*) "In common-densities/calcmomenta.f" 
c     write(*,*) "Check equality of the next few"
c     write(*,*) "Check from density: k?=omega:  k=", k
c     write(*,*) ""
c     write(*,*) "omega check with mandal: omega = k= s-M^2/(2sqrt(s))"
c     write(*,*) "k?=(mandalS-M*M)/(2*sqrt(s))",k,"?=",(mandalS-mNucl*mNucl)/(2*sqrt(mandalS))
c     write(*,*) ""
c     write(*,*) "E_pion check with mandalstam"
c     write(*,*) "E_pi= sqrt(m^2+k'^2)=(s+m^2-M^2)/(2sqrt(s)) -- Next line"
c     write(*,*) sqrt(mPion**2+kpsq),"?=",(mandalS+mPion**2-mNucl**2)/(2*sqrt(mandalS))
c     write(*,*) "#################################################################################"
c     write(*,*) ""
      if (verbosity.eq.1000) continue
      return
      end 
