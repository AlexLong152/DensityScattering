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
      subroutine CalculateQs(qx,qy,qz,q12x,q12y,q12z,qpx,qpy,qpz,
     &     qp12x,qp12y,qp12z,qppx,qppy,qppz,qpp12x,qpp12y,qpp12z,
     &     qpppx,qpppy,qpppz,qppp12x,qppp12y,qppp12z,
     &     qsq,qpsq,qppsq,qpppsq,q12sq,qp12sq,qpp12sq,qppp12sq,px,py,pz,
     &     ppx,ppy,ppz,
     &     k,thetacm,verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
c
      real*8,intent(in)  :: px,py,pz,ppx,ppy,ppz,k,thetacm
      integer,intent(in) :: verbosity
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES:
c
      real*8,intent(out) :: qx,qy,qz,qpx,qpy,qpz,qppx,qppy,qppz,qsq,qpsq,qppsq
      real*8,intent(out) :: q12x,q12y,q12z,qp12x,qp12y,qp12z,qpp12x,qpp12y,qpp12z
      real*8,intent(out) :: q12sq,qp12sq,qpp12sq
      real*8,intent(out) :: qpppx,qpppy,qpppz, qpppsq
      real*8,intent(out) :: qppp12x,qppp12y,qppp12z, qppp12sq
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      qx=px - ppx + k*dsin(thetacm)/2.d0
      qy=py - ppy
      qz=pz - ppz + k*(1.d0 + dcos(thetacm))/2.d0
      qpppx=px + ppx - k*dsin(thetacm)/2.d0
      qpppy=py + ppy
      qpppz=pz + ppz - k*(1.d0 + 3.d0*dcos(thetacm))/6.d0
      qppp12x=-px - ppx - k*dsin(thetacm)/2.d0
      qppp12y=-py - ppy
      qppp12z=-pz - ppz - k*(1.d0 + 3.d0*dcos(thetacm))/6.d0
      q12x=-px + ppx + k*dsin(thetacm)/2.d0
      q12y=-py + ppy
      q12z=-pz + ppz + k*(1.d0 + dcos(thetacm))/2.d0
      qpx=px - ppx + k*dsin(thetacm)/2.d0
      qpy=py - ppy
      qpz=pz - ppz + k*(dcos(thetacm)- 1.d0)/2.d0
      qppx=px - ppx - k*dsin(thetacm)/2.d0
      qppy=py - ppy 
      qppz=pz - ppz + k*(1.d0 - dcos(thetacm))/2.d0
      qp12x=-px + ppx + k*dsin(thetacm)/2.d0
      qp12y=-py + ppy
      qp12z=-pz + ppz + k*(dcos(thetacm) - 1.d0)/2.d0
      qpp12x=-px + ppx - k*dsin(thetacm)/2.d0
      qpp12y=-py + ppy
      qpp12z=-pz + ppz + k*(1.d0 - dcos(thetacm))/2.d0
      qsq=qx**2 + qy**2 + qz**2
      q12sq=q12x**2 + q12y**2 + q12z**2
      qpsq=qpx**2 + qpy**2 + qpz**2
      qppsq=qppx**2 + qppy**2 + qppz**2
      qp12sq=qp12x**2 + qp12y**2 + qp12z**2
      qpp12sq=qpp12x**2 + qpp12y**2 + qpp12z**2
      qpppsq=qpppx**2 + qpppy**2 + qpppz**2
      qppp12sq=qppp12x**2 + qppp12y**2 + qppp12z**2


      
      if (verbosity.eq.1000) continue
      return
      end


      subroutine calculateqsmass(p,pp,q,k,q1,kVec,kp,thetacm,mPion,mNucl,verbosity)
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
      real*8 p(3), pp(3), k,thetacm, mPion, mNucl
      integer verbosity 
c
c**********************************************************************
c
c     temporary variables
    
      real*8 Epion, mandalS, ENuc, kpsq, kpAbs, omegaThreshold
      real*8 kp(3), q(3)
      real*8 q1(3), kVec(3)
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
      omegaThreshold=(mPion*(mPion+2*mNucl))/(2*(mPion+mNucl))
      ENuc=sqrt((mNucl**2) + (k**2))
      mandalS=(ENuc + k)**2 !lab frame
      Epion=(mandalS+(mPion**2)-(mNucl**2))/(2*sqrt(mandalS))
      kpsq=(((mandalS+mPion**2-mNucl**2)**2)/(4*mandalS))-mPion**2
      kpAbs=sqrt(kpsq)

      kVec=(/0.d0,0.d0,k/)
      kp=(/0.d0,kpAbs*sin(thetacm), kpAbs*cos(thetacm)/)
      write(*,*) "pion photoproduction"
      write(*,*) "In calculateQs.PionPion.f: kp=",kp 
      q = (p-pp)+((kVec+kp)/2)
      q1 = q-k
      if (verbosity.eq.1000) continue
      return
      end 

      subroutine triangle(val,a,b,c)
      implicit none
      real*8 a,b,c, val
      val=a*a+b*b+c*c-2*(a*b+a*c+b*c)
      end subroutine triangle

      subroutine calculateqs2Mass(p,pp,k,kp,m1,m2,m3,m4,thetacm,verbosity)
c     m1, is nucleous, m2 is incoming particle
c     m3 is nucleous, m4 is outgoing particle 
c     current implimentation is for elastic scattering based on the below
c     https://edu.itp.phys.ethz.ch/hs10/ppp1/PPP1_2.pdf equations 2.9 and 2.10
c     But its in many places including the PDG
c     Kinematics are for arbitrary inelastic scattering
c     p_1=(E_1,\vec{p}) p_2=(E_2,-\vec{p}) p_3=(E_3,\vec{p}') p_4=(E_4,-\vec{p}')
c     \vec{p}_1+\vec{p}_2=\vec{p}_3+\vec{p}_4=0
c     Derivation for the kinematics can be found in pionpionAngle.pdf
c     OneDrive/thesis/Kinematics/pionpionAngle.pdf
c     The conversion from lenkewitz to our variables is found in
c     DensityScattering/documentation/pionpionangle.pdf
      implicit none
c**********************************************************************
c     Input variables:
      real*8 p(3), pp(3),k(3),m1,m2,m3,m4,thetacm
      integer verbosity 

c**********************************************************************
c     Output variables       
      real*8 kp(3)

c**********************************************************************
c     Internal variables

c     real*8 val
      real*8 mandalS, sqrtS
      real*8 kpAbs,kpSquare!, kpNuclAbs
      real*8 E1,E2,E3,E4
      real*8 E1check,E2check,E4check
      E1=sqrt(m1*m1+DOT_PRODUCT(k,k))
      E2=sqrt(m2*m2+DOT_PRODUCT(k,k))!don't need the negative sign on k bc it cancels out anyways
      sqrtS=E1+E2
      mandalS=sqrtS*sqrtS

      E1check=(1/(2*sqrt(mandalS)))*(mandalS+m1*m1-m2*m2)
      E2check=(1/(2*sqrt(mandalS)))*(mandalS+m2*m2-m1*m1)

      E3=(1/(2*sqrt(mandalS)))*(mandalS+m3*m3-m4*m4)
      E4=E1+E2-E3
      E4check=(1/(2*sqrt(mandalS)))*(mandalS+m4*m4-m3*m3)
      kpSquare=E4*E4-m4*m4
      kpAbs=sqrt(kpSquare)
      
c     write(*,*) "E1?=E1check",E1,E1check
c     write(*,*) "E2?=E2check",E2,E2check
c     write(*,*) "E2=",E2
c     write(*,*) "E4?=E4check",E4,E4check
c     write(*,*) "Total check", E1+E2,E3+E4
c     write(*,*) "kpSquare=",kpSquare
      kp=(/0.d0,kpAbs*sin(thetacm), kpAbs*cos(thetacm)/)
      if (verbosity.eq.1000) continue
      return
      end 
