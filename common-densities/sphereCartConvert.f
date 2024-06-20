cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              CalculatePVector : compute Cartesian components given spherical polar co-ordinates 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, was part of common-densities/calcmomenta.f of Compton density code v2.0 hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     specific to process considered
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sphere2cart(qx,qy,qz,q,thq,phiq,verbosity)
c
c  Calculate Cartesian components given spherical polar co-ordinates
c     
c**********************************************************************
c
      implicit none
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
c
      real*8,intent(in)  :: q,thq,phiq
      integer,intent(in) :: verbosity
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES:
      
      real*8,intent(out)  :: qx,qy,qz
c
c**********************************************************************
c
      qx=q*dsin(thq)*dcos(phiq)
      qy=q*dsin(thq)*dsin(phiq)
      qz=q*dcos(thq)
      if (verbosity.eq.1000) continue
      return
      end

      subroutine cart2sphere(qx,qy,qz,q,thq,phiq,verbosity)
c
c  Calculate spherical components given cartesian 
c     
c**********************************************************************
c
      implicit none
      include '../common-densities/constants.def'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
c
      real*8,intent(in)  :: qx,qy,qz
      integer,intent(in) :: verbosity
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES:
      
      real*8,intent(out)  :: q,thq,phiq
c
c**********************************************************************
c
c     write(*,*) "In sphereCartConvert: PI=",PI 
c     write(*,*) "In sphereCartConvert: dcos(PI/2.d0)=",dcos(PI/2.d0)
c     write(*,*) "In sphereCartConvert: acos(0)=",acos(0.d0) 
c     write(*,*) "In sphereCartConvert: acos(1)=",acos(1.d0) 
c     write(*,*) "In sphereCartConvert: sign(1.d0,-0.2d0)=",sign(1.d0,-0.2d0) 
c     write(*,*) "In sphereCartConvert: sign(1.d0,0.2d0)=",sign(1.d0,0.2d0) 
c     write(*,*) "In sphereCartConvert: sign(1.d0,0.d0)=",sign(1.d0,0.d0) 
c     stop
      q=sqrt(qx*qx+qy*qy+qz*qz)
      thq=acos(qz/q)
      phiq=sign(1.d0,qy)*acos(qx/sqrt(qx*qx+qy*qy))

      if(thq.lt.0) then !same as taking mod 2pi
          thq=thq+2.d0*PI
      end if

      if(phiq.lt.0) then !same as taking mod 2pi
          phiq=phiq+2.d0*PI
      end if

      if ((0.d0.eq.qx).and.(0.d0.eq.qy)) then
          phiq=0.d0
      end if

      if (abs(thq-PI).le.1e-6) then
          phiq=0.d0!ambiguous edge case, assigning to zero makes this match with Lebedev Laikov integration routine
      end if

      if (verbosity.eq.1000) continue
      return
      end

