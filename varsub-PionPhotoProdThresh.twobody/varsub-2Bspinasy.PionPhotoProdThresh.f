cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of KERNEL code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              CalcKernel2BAasy : set up (1↔2) asymmetric part of the kernel
c              CalcKernel2BBasy : set up (1↔2) asymmetric part of the kernel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     hgrie 17 Nov 2023: split the following subroutines into new file spinstructures.f and renamed two for more intuitive names:

c         singlesigmaasy => singlesigmaasy
c         Calcholdasy    => doublesigmaasy
c      
c     This way, spintricks*f only contains individual diagram
c     contributions and not these routines which are generally relevant for spin structures.
c      
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      

      subroutine StaticKernelEasym(Kernel,
     &     factor,qVec,qpVec,pVec,
     &     Sp,S,extQnumlimit,verbosity)
      !B.31 from Lenkewitz thesis
      implicit none
      include '../common-densities/constants.def'

      complex*16,intent(inout) :: Kernel(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
      real*8,intent(in)  :: factor
      real*8,intent(in)  :: qVec(3),qpVec(3),pVec(3)
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity

      real*8 eps(3,3), epsVec(3)
      integer ieps
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      complex*16 Ihold
      external Ihold
c     real*8 tmpVec(3),kCrossEps(3)
      integer Msp,Ms

      eps = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      hold=c0

      do ieps=1,3
        epsVec=eps(ieps,:)
        call singlesigmaasy(hold,qpVec(1),qpVec(2),qpVec(3),Sp,S,verbosity)
        do Msp=-Sp,Sp
        do Ms=-S,S
             Kernel(ieps,Sp,Msp,S,Ms) = Kernel(ieps,Sp,Msp,S,Ms) + factor*(
     &        hold(Sp,Msp,S,Ms)*dot_product(epsVec,qVec+qpVec))
        end do
        end do  

      end do
      end subroutine StaticKernelEasym 

      subroutine StaticKernelDasym(Kernel,
     &     factor,qVec,ppVec,kVec,
     &     Sp,S,extQnumlimit,verbosity)
      !B.32 from Lenkewitz thesis
      implicit none
      include '../common-densities/constants.def'

      complex*16,intent(inout) :: Kernel(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
      real*8,intent(in)  :: factor
      real*8,intent(in)  :: qVec(3),ppVec(3),kVec(3)
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity

      real*8 eps(3,3), epsVec(3)
      integer ieps
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      complex*16 Ihold
      external Ihold
      real*8 tmpVec(3),kCrossEps(3)
      integer Msp,Ms

      eps = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      hold=c0
      tmpVec=qVec+2.d0*ppVec-kVec

      do ieps=1,3
        epsVec=eps(ieps,:)
        call cross(qVec-kVec,epsVec,kCrossEps)
        do Msp=-Sp,Sp
        do Ms=-S,S
             !TODO: check this, just added verbosity and S, Sp arguments to this
             call doublesigmaasy(hold,kCrossEps(1),kCrossEps(2),kCrossEps(3),qVec(1),qVec(2),qVec(3),Sp,S,verbosity)
             Kernel(ieps,Sp,Msp,S,Ms) = Kernel(ieps,Sp,Msp,S,Ms) + factor*(ci*hold(Sp,Msp,S,Ms)
     &      -Ihold(Sp,Msp,S,Ms)*2.d0*dot_product(epsVec,tmpVec))
        end do
        end do  

      end do
      end subroutine StaticKernelDasym 

      subroutine StaticKernelCasym(Kernel,
     &     factor,qVec,ppVec,kVec,
     &     Sp,S,extQnumlimit,verbosity)
      !B.30 from Lenkewitz thesis
      implicit none
      include '../common-densities/constants.def'

      complex*16,intent(inout) :: Kernel(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
      real*8,intent(in)  :: factor
      real*8,intent(in)  :: qVec(3),ppVec(3),kVec(3)
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity

      real*8 eps(3,3), epsVec(3)
      integer ieps
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      complex*16 Ihold
      external Ihold
      real*8 tmpVec(3),kCrossEps(3)
      integer Msp,Ms

      eps = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      hold=c0
      tmpVec=qVec+2.d0*ppVec-kVec

      do ieps=1,3
        epsVec=eps(ieps,:)
        call cross(kVec,epsVec,kCrossEps)
        kCrossEps=kCrossEps*(1+kappanu)
        do Msp=-Sp,Sp
        do Ms=-S,S
          !TODO: check this, just added Sp, S, verbossity to doublesigma
             call doublesigmasym(hold,kCrossEps(1),kCrossEps(2),kCrossEps(3),qVec(1),qVec(2),qVec(3),Sp,S,verbosity)
             Kernel(ieps,Sp,Msp,S,Ms) = Kernel(ieps,Sp,Msp,S,Ms) + factor*(ci*hold(Sp,Msp,S,Ms)
     &      +Ihold(Sp,Msp,S,Ms)*2.d0*dot_product(epsVec,tmpVec))
        end do
        end do

      end do

      if (verbosity.eq.1000) continue ! unused variable, kept for future use
      end subroutine StaticKernelCasym 

      subroutine StaticKernelBasym(Kernel,
     &     factor,qVec,ppVec,kVec,
     &     Sp,S,extQnumlimit,verbosity)

      implicit none
      include '../common-densities/constants.def'

      complex*16,intent(inout) :: Kernel(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
      real*8,intent(in)  :: factor
      real*8,intent(in)  :: qVec(3),ppVec(3),kVec(3)
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity

      real*8 eps(3,3), epsVec(3)
      integer ieps
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      complex*16 Ihold
      external Ihold
      real*8 tmpVec(3),tmpVec2(3)
      integer Msp,Ms

      eps = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      hold=c0
      tmpVec=qVec+2.d0*ppVec-kVec
    
      call singlesigmaasy(hold,tmpVec(1),tmpVec(2),tmpVec(3),Sp,S,verbosity)
      call cross(qVec,kvec,tmpVec2)

      tmpVec2=tmpVec2*(1+kappanu)
      do ieps=1,3!extQnumlimit=3
        epsVec=eps(ieps,:)

        do Msp=-Sp,Sp
        do Ms=-S,S
             Kernel(ieps,Sp,Msp,S,Ms) = Kernel(ieps,Sp,Msp,S,Ms) + factor*(
     &       hold(Sp,Msp,S,Ms)*dot_product(epsVec,tmpVec) +
     &      ci*dot_product(tmpVec2,epsVec)*Ihold(Sp,Msp,S,Ms)
     &      )
        end do
        end do  

      end do
      end subroutine StaticKernelBasym 


      subroutine StaticKernelAasym(Kernel2B,
     &     factor,
     &     Sp,S,extQnumlimit,verbosity)

      implicit none
      include '../common-densities/constants.def'


      complex*16,intent(inout) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
      real*8,intent(in)  :: factor
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity

      complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms
c     εx:
      call singlesigmasym(hold,1.d0,0.d0,0.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(1,Sp,Msp,S,Ms) = Kernel2B(1,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εy:
      call singlesigmasym(hold,0.d0,1.d0,0.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(2,Sp,Msp,S,Ms) = Kernel2B(2,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εz:
      call singlesigmasym(hold,0.d0,0.d0,1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(3,Sp,Msp,S,Ms) = Kernel2B(3,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      
      if (verbosity.eq.1000) continue
      return
      end subroutine StaticKernelAasym
      subroutine CalcKernel2BAasy(Kernel2B,
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
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     INPUT/OUTPUT VARIABLE:
c     
      complex*16,intent(inout) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c      complex*16 Kernel2Bpx(0:1,-1:1,0:1,-1:1),Kernel2Bpy(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8,intent(in)  :: factor
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     LOCAL VARIABLES:
c      
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     εx:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      call singlesigmaasy(hold,1.d0,0.d0,0.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(1,Sp,Msp,S,Ms) = Kernel2B(1,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εy: 
      call singlesigmaasy(hold,0.d0,1.d0,0.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(2,Sp,Msp,S,Ms) = Kernel2B(2,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εz:  
      call singlesigmaasy(hold,0.d0,0.d0,1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(3,Sp,Msp,S,Ms) = Kernel2B(3,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
      if (verbosity.eq.1000) continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine CalcKernel2BBasyVec(Kernel2B,
     &     factor,
     &     A,B, ! A.σ, B.ε
     &     Sp,S,extQnumlimit,verbosity)    
c     
c********************************************************************
c     
c     Calculates diagram B
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     INPUT/OUTPUT VARIABLES:
c     
      complex*16,intent(inout) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c      complex*16 Kernel2Bpx(0:1,-1:1,0:1,-1:1),Kernel2Bpy(0:1,-1:1,0:1,-1:1)
      
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8,intent(in)  :: factor
      real*8,intent(in)  :: A(3), B(3)
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     LOCAL VARIABLES:
c      
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call singlesigmaasy(hold,A(1),A(2),A(3),Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
c     εx:
            Kernel2B(1,Sp,Msp,S,Ms) = Kernel2B(1,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)*B(1)
c     εy:
            Kernel2B(2,Sp,Msp,S,Ms) = Kernel2B(2,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)*B(2)
c     εz:
            Kernel2B(3,Sp,Msp,S,Ms) = Kernel2B(3,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)*B(3)
         end do
      end do  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      
      if (verbosity.eq.1000) continue
      return
      end
