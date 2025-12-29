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
     &           factor,
     &           s12p,s12,t12,mt12,extQnumlimit,verbosity)
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
      integer,intent(in) :: s12p,s12,t12,mt12
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     LOCAL VARIABLES:
c      
c     complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms, extQnum
      real*8 tmp,tmp2, ddelta
      real*8 isospin

      isospin=((-1)**(t12))*ddelta(mt12,0)
      do extQnum=1,extQnumlimit!no dependence on this either
      do Msp=-s12p,s12p
      do Ms=-s12,s12
            tmp=ddelta(s12p,s12)
            tmp2=ddelta(Msp,Ms)
            Kernel2B(extQnum,s12p,Msp,s12,Ms) = Kernel2B(extQnum,s12p,Msp,s12,Ms) +
     &              factor*isospin*tmp*tmp2

      end do
      end do
      end do 
c     Note the factor 2*((-1)**(2*t12+1)) only appears in the case where t12=t12p and mt12=mt12p but this is taken care
c     of in 2Bkernel.PionPion.f
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (verbosity.eq.1000) continue
      return
      end

      subroutine CalcKernel2BBsym(Kernel2B,qVec,
     &           factor,
     &           s12p,s12,t12,extQnumlimit,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram B
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
      real*8, intent(in) :: qVec(3)
      complex*16,intent(in) :: factor
      integer,intent(in) :: s12p,s12,t12
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     LOCAL VARIABLES:
c      
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms, extQnum
      real*8 isospin
     
      isospin=(2*t12*(t12+1))-3
      call doublesigmasym(hold,qVec(1),qVec(2),qVec(3),qVec(1),qVec(2),qVec(3),s12p,s12,verbosity)
      do extQnum=1,extQnumlimit
      do Msp=-s12p,s12
      do Ms=-s12,s12
            Kernel2B(extQnum,s12p,Msp,s12,Ms) = Kernel2B(extQnum,s12p,Msp,s12,Ms) + factor*hold(s12p,Msp,s12,Ms)*isospin
      end do
      end do
      end do 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (verbosity.eq.1000) continue
      return
      end

      subroutine CalcKernel2BCsym(Kernel2B,qVec,
     &           factor,
     &           s12p,s12,t12,extQnumlimit,verbosity)
c********************************************************************
c     
c     Calculates diagram C
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
      real*8, intent(in) :: qVec(3)
      complex*16,intent(in) :: factor
      integer,intent(in) :: s12p,s12,t12
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     LOCAL VARIABLES:
c      
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms, extQnum
      real*8 isospin

      isospin= (2*t12*(t12+1))-3
      do extQnum=1,extQnumlimit
      do Msp=-s12p,s12
      do Ms=-s12,s12
            call doublesigmasym(hold,qVec(1),qVec(2),qVec(3),qVec(1),qVec(2),qVec(3),s12p,s12,verbosity)
            Kernel2B(extQnum,s12p,Msp,s12,Ms) = Kernel2B(extQnum,s12p,Msp,s12,Ms) + factor*hold(s12p,Msp,s12,Ms)*isospin
      end do
      end do
      end do 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (verbosity.eq.1000) continue
      return
      end


      subroutine CalcKernel2BDsym(Kernel2B,qVec,
     &           factor,
     &           s12p,s12,mt12,extQnumlimit,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram D
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
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8, intent(in) :: qVec(3)
      complex*16,intent(in) :: factor
      integer,intent(in) :: s12p,s12,mt12
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     LOCAL VARIABLES:
c      
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms, extQnum
      real*8 isospin

      isospin=(-1)**(mt12)
      call doublesigmasym(hold,qVec(1),qVec(2),qVec(3),qVec(1),qVec(2),qVec(3),s12p,s12,verbosity)

      do extQnum=1,extQnumlimit
      do Msp=-s12p,s12p
      do Ms=-s12,s12
            Kernel2B(extQnum,s12p,Msp,s12,Ms) = Kernel2B(extQnum,s12p,Msp,s12,Ms) + factor*hold(s12p,Msp,s12,Ms)*isospin
      end do
      end do
      end do 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (verbosity.eq.1000) continue
      return
      end
