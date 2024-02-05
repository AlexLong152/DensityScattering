cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of KERNEL code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              CalcKernel2BAsym : set up (1↔2) symmetric piece of diagram A
c              CalcKernel2BBsym : set up (1↔2) symmetric piece of diagram B
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     hgrie 17 Nov 2023: split the following subroutines into new file spinstructures.f and renamed two for more intuitive names:

c         singlesigma => singlesigmasym
c         Calchold    => doublesigmasym
c      
c     This way, spintricks*f only contains individual diagram
c     contributions and not these routines which are generally relevant for spin structures.
c    
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      subroutine CalcKernel2BAsym(Kernel2B,
     &     factor,
     &     Sp,S,extQnumlimit,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram A
c     
c********************************************************************
c     
      implicit none
      include '../common-densities/constants.def'
c     
c********************************************************************
c     INPUT/OUTPUT VARIABLES:
c     
      complex*16,intent(inout) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c      complex*16 Kernel2Bpx(0:1,-1:1,0:1,-1:1),Kernel2Bpy(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      complex*16,intent(in) :: factor
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     LOCAL VARIABLES:
c      
c     complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms, extQnum
      real*8 tmp,tmp2, ddelta
     
c     hold = identity in this case since theres no explicit spin depedence, only isospin
      do extQnum=1,3
      do Msp=-Sp,Sp
      do Ms=-S,S
            tmp=ddelta(Sp,S)
            tmp2=ddelta(Msp,Ms)
            Kernel2B(extQnum,Sp,Msp,S,Ms) = Kernel2B(extQnum,Sp,Msp,S,Ms) + factor*tmp*tmp2
      end do
      end do
      end do 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (verbosity.eq.1000) continue
      return
      end

      subroutine CalcKernel2BBsym(Kernel2B,q,
     &     factor,
     &     Sp,S,extQnumlimit,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram A
c     
c********************************************************************
c     
      implicit none
      include '../common-densities/constants.def'
c     
c********************************************************************
c     INPUT/OUTPUT VARIABLES:
c     
      complex*16,intent(inout) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c      complex*16 Kernel2Bpx(0:1,-1:1,0:1,-1:1),Kernel2Bpy(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8, intent(in) :: q(3)
      complex*16,intent(in) :: factor
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     LOCAL VARIABLES:
c      
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms, extQnum
c     real*8 tmp,tmp2, ddelta
     
c     hold = identity in this case since theres no explicit spin depedence, only isospin
      call doublesigmasym(hold,q(1),q(2),q(3),q(1),q(2),q(3),Sp,S,verbosity)
      do extQnum=1,3
      do Msp=-Sp,Sp
      do Ms=-S,S
            Kernel2B(extQnum,Sp,Msp,S,Ms) = Kernel2B(extQnum,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
      end do
      end do
      end do 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (verbosity.eq.1000) continue
      return
      end
