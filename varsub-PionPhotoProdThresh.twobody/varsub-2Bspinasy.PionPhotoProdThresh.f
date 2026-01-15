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
     &     factor,qVec,qpVec,
     &     Sp,S,extQnumlimit,verbosity)
      !B.31 from Lenkewitz thesis
      implicit none
      include '../common-densities/constants.def'

      complex*16,intent(inout) :: Kernel(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
      real*8,intent(in)  :: factor
      real*8,intent(in)  :: qVec(3),qpVec(3)
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
      complex*16 holdPart1(0:1,-1:1,0:1,-1:1), holdPart2(0:1,-1:1,0:1,-1:1)  
      complex*16 part1(0:1,-1:1,0:1,-1:1), part2(0:1,-1:1,0:1,-1:1) 
      complex*16 Ihold
      external Ihold
      real*8 tmpVec(3),CrossVec(3), qMinusK(3)
      integer Msp,Ms

      eps = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      holdPart2=c0
      tmpVec=qVec+2.d0*ppVec-kVec
      qMinusK=qVec-kVec
c     This has a similar structure to diagram C
c     [-ε·(q+2p'-k)](σ₂·q) + i[ε·(σ₁×(q-k))](σ₂·q)
      do ieps=1,3
        epsVec=eps(ieps,:)
c       Want to compute epsVec dot ( sigmaVec x kVec)=sigmaVec dot (kVec x epsVec)
        
        call singlesigmaasy(holdPart1,qVec(1),qVec(2),qVec(3),Sp,S,verbosity)!holdPart1=(σ₂·q)
        part1= holdPart1*dot_product(epsVec,-1.d0*tmpVec)!part1 = [ε·(q+2p'-k)](σ₂·q)

        call cross(qMinusK,epsVec,CrossVec)
c       σ_1·[(q-k)x ε](σ₂·q)i
        call doublesigmaasy(holdPart2,CrossVec(1),CrossVec(2),CrossVec(3),
     &    qVec(1),qVec(2),qVec(3),Sp,S,verbosity)
        part2 = ci*holdPart2
        do Msp=-Sp,Sp
        do Ms=-S,S
          Kernel(ieps,Sp,Msp,S,Ms) = Kernel(ieps,Sp,Msp,S,Ms) + factor*(
     &        part1(Sp,Msp,S,Ms)+part2(Sp,Msp,S,Ms))
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
      complex*16 holdPart1(0:1,-1:1,0:1,-1:1), holdPart2(0:1,-1:1,0:1,-1:1)  
      complex*16 part1(0:1,-1:1,0:1,-1:1), part2(0:1,-1:1,0:1,-1:1) 
      complex*16 Ihold
      external Ihold
      real*8 tmpVec(3),kCrossEps(3)
      integer Msp,Ms

      eps = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      holdPart2=c0
      tmpVec=qVec+2.d0*ppVec-kVec
c     This expands to
c     [ε·(q+2p'-k)](σ₂·q) + i[ε·(σ₁×k)](1+κᵥ)(σ₂·q)
c    =[ε·(q+2p'-k)](σ₂·q) + [σ₁·(k×ε)](σ₂·q) i(1+κᵥ)
      do ieps=1,3
        epsVec=eps(ieps,:)
c       Want to compute epsVec dot ( sigmaVec x kVec)=sigmaVec dot (kVec x epsVec)
        
        call singlesigmaasy(holdPart1,qVec(1),qVec(2),qVec(3),Sp,S,verbosity)!holdPart1=(σ₂·q)
        part1= holdPart1*dot_product(epsVec,tmpVec)!part1 = [ε·(q+2p'-k)](σ₂·q)

        call cross(kVec,epsVec,kCrossEps)
        kCrossEps=kCrossEps
c       [σ₁·(k×ε)](σ₂·q)(1+κᵥ)
        call doublesigmaasy(holdPart2,kCrossEps(1),kCrossEps(2),kCrossEps(3),
     &    qVec(1),qVec(2),qVec(3),Sp,S,verbosity)
        part2 = (1+kappanu)*ci*holdPart2
        do Msp=-Sp,Sp
        do Ms=-S,S
          Kernel(ieps,Sp,Msp,S,Ms) = Kernel(ieps,Sp,Msp,S,Ms) + factor*(
     &        part1(Sp,Msp,S,Ms)+part2(Sp,Msp,S,Ms))
        end do
        end do  

      end do
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
      complex*16 hold2(0:1,-1:1,0:1,-1:1)
      complex*16 Ihold
      external Ihold
      real*8 tmpVec(3),tmpVec2(3)
      integer Msp,Ms
      eps = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      hold=c0
      tmpVec=qVec+2.d0*ppVec-kVec
    
      call singlesigmaasy(hold,qVec(1),qVec(2),qVec(3),Sp,S,verbosity)
      do ieps=1,3!extQnumlimit=3
        epsVec=eps(ieps,:)
        call cross(epsVec,kvec,tmpVec2)
        tmpVec2=tmpVec2*(1+kappanu)
        !use a dot( b x c)= b dot (c x a)- > eps dot (sigma x k)= sigma dot( k x eps) = sigma dot tmpVec2
        call singlesigmaasy(hold2,tmpVec2(1),tmpVec2(2),tmpVec2(3),Sp,S,verbosity)
        do Msp=-Sp,Sp
        do Ms=-S,S
             Kernel(ieps,Sp,Msp,S,Ms) = Kernel(ieps,Sp,Msp,S,Ms) + factor*(
     &       hold(Sp,Msp,S,Ms)!\sigma dot qVec
     &       *dot_product(epsVec,tmpVec) +!eps dot qVec+2*ppVec-kVec
     &       ci*hold2(Sp,Msp,S,Ms) ! eps dot tmpVec2
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
