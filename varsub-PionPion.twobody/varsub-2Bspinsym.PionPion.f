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
     &           s12p,s12,t12,mt12,t12p,mt12p,extQnum,verbosity)
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
      complex*16,intent(inout) :: Kernel2B(0:1,-1:1,0:1,-1:1)
c      complex*16 Kernel2Bpx(0:1,-1:1,0:1,-1:1),Kernel2Bpy(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      complex*16,intent(in) :: factor
      integer,intent(in) :: s12p,s12,t12,mt12,t12p,mt12p
      integer,intent(in) :: extQnum
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     LOCAL VARIABLES:
c      
c     complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms
      real*8 tmp,tmp2, ddelta
      complex*16 isospin

      call PionPionA(t12,mt12,t12p,mt12p,extQnum,isospin)
      do Msp=-s12p,s12p
      do Ms=-s12,s12
            tmp=ddelta(s12p,s12)
            tmp2=ddelta(Msp,Ms)
            Kernel2B(s12p,Msp,s12,Ms) = Kernel2B(s12p,Msp,s12,Ms) +
     &              factor*tmp*tmp2*isospin

      end do
      end do 
c     Note the factor 2*((-1)**(2*t12+1)) only appears in the case where t12=t12p and mt12=mt12p but this is taken care
c     of in 2Bkernel.PionPion.f
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (verbosity.eq.1000) continue
      return
      end


      subroutine CalcKernel2BBCsym(Kernel2B,qVec,
     &           factor,
     &           s12p,s12,t12,mt12, t12p, mt12p, extQnum,verbosity)
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
      complex*16,intent(inout) :: Kernel2B(0:1,-1:1,0:1,-1:1)
c      complex*16 Kernel2Bpx(0:1,-1:1,0:1,-1:1),Kernel2Bpy(0:1,-1:1,0:1,-1:1)
c
c********************************************************************
c     INPUT VARIABLES:
c
      real*8, intent(in) :: qVec(3)
      complex*16,intent(in) :: factor
      integer,intent(in) :: s12p,s12,t12,mt12, t12p, mt12p,extQnum
      integer,intent(in) :: verbosity
      
c     
c********************************************************************
c     LOCAL VARIABLES:
c      
      complex*16 hold(0:1,-1:1,0:1,-1:1), isospin
      integer Msp,Ms

      call doublesigmasym(hold,qVec(1),qVec(2),qVec(3),qVec(1),qVec(2),qVec(3),s12p,s12,verbosity)
      call PionPionBC(t12, mt12, t12p, mt12p, extQnum, isospin)
      do Msp=-s12p,s12p
      do Ms=-s12,s12
            Kernel2B(s12p,Msp,s12,Ms) = Kernel2B(s12p,Msp,s12,Ms) + factor*hold(s12p,Msp,s12,Ms)*isospin
      end do
      end do 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (verbosity.eq.1000) continue
      return
      end


      subroutine CalcKernel2BDsym(kernel2b,
     &           factor,
     &           s12p,s12,t12,mt12, t12p, mt12p, extqnum,verbosity)
c********************************************************************
c
c     calculates diagram D
c
c********************************************************************
c
      implicit none
      include '../common-densities/constants.def'
c
c********************************************************************
c     input/output variables:
c
      complex*16,intent(inout) :: kernel2b(0:1,-1:1,0:1,-1:1)
c      complex*16 kernel2bpx(0:1,-1:1,0:1,-1:1),kernel2bpy(0:1,-1:1,0:1,-1:1)
c
c********************************************************************
c     input variables:
c
      complex*16,intent(in) :: factor
      integer,intent(in) :: s12p,s12,t12,mt12, t12p, mt12p,extqnum
      integer,intent(in) :: verbosity
      
c     
c********************************************************************
c     local variables:
c      
      complex*16 hold(0:1,-1:1,0:1,-1:1), isospin
      integer msp,ms
      real*8 tmp,tmp2, ddelta

      call PionPionD(t12,mt12,t12p,mt12p,extqnum,isospin)
      do msp=-s12p,s12p
      do ms=-s12,s12
            tmp=ddelta(s12p,s12)
            tmp2=ddelta(msp,ms)
            kernel2b(s12p,msp,s12,ms) = kernel2b(s12p,msp,s12,ms) +
     &              factor*tmp*tmp2*isospin

      end do
      end do 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (verbosity.eq.1000) continue
      return
      end
