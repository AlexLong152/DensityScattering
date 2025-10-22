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

      subroutine StaticKernelDsym(Kernel,
     &     factor,qVec,ppVec,kVec,
     &     Sp,S,extQnumlimit,verbosity)
      !B.31 from Lenkewitz thesis
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
        kCrossEps=kCrossEps*(1+kappanu)
        call doublesigmasym(hold,tmpVec(1),tmpVec(2),tmpVec(3),qVec(1),qVec(2),qVec(3))
        do Msp=-Sp,Sp
        do Ms=-S,S
             Kernel(ieps,Sp,Msp,S,Ms) = Kernel(ieps,Sp,Msp,S,Ms) + factor*(ci*hold(Sp,Msp,S,Ms)
     &      -Ihold(Sp,Msp,S,Ms)*2.d0*dot_product(epsVec,tmpVec))
        end do
        end do  

      end do
      end subroutine StaticKernelDsym 


      subroutine StaticKernelCsym(Kernel,
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
        call doublesigmasym(hold,tmpVec(1),tmpVec(2),tmpVec(3),qVec(1),qVec(2),qVec(3))
        do Msp=-Sp,Sp
        do Ms=-S,S
             Kernel(ieps,Sp,Msp,S,Ms) = Kernel(ieps,Sp,Msp,S,Ms) + factor*(ci*hold(Sp,Msp,S,Ms)
     &      +Ihold(Sp,Msp,S,Ms)*2.d0*dot_product(epsVec,tmpVec))
        end do
        end do  

      end do
      end subroutine StaticKernelCsym 

      subroutine StaticKernelBsym(Kernel,
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
      logical holdnan
      eps = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      hold=c0
      tmpVec=qVec+2.d0*ppVec-kVec
    
      call singlesigmasym(hold,tmpVec(1),tmpVec(2),tmpVec(3),Sp,S,verbosity)
      call cross(qVec,kvec,tmpVec2)

      if(HoldNan(hold)) then
        write(*,*) "hold=", hold 
        write(*,*) ""
        write(*,*) "tmpVec=", tmpVec 
        write(*,*) "Sp=", Sp 
        write(*,*) "S=", S 
      end if
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
      end subroutine StaticKernelBsym 


      subroutine StaticKernelAsym(Kernel2B,
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
      logical HoldNan !checks if Hold is NaN
c     εx:
      call singlesigmasym(hold,1.d0,0.d0,0.d0,Sp,S,verbosity)

      if(HoldNan(hold)) then

        write(*,*) "In diag A static"
        write(*,*) "hold=", hold 
        write(*,*) "eps=1,0,0" 
        write(*,*) "Sp=", Sp 
        write(*,*) "S=", S 
        write(*,*) "NaN on eps=1,0,0"
        stop
      end if
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
      end subroutine StaticKernelAsym

      logical function HoldNan(hold)
      ! returns false if hold has no NaN values
      ! returns true if hold has NaN values
      implicit none
      complex*16,intent(in) :: hold(0:1,-1:1,0:1,-1:1)
      complex*16 tmp
      integer S,Sp,Ms,Msp

      HoldNan=.false.
      do S=0, 1
      do Ms = -1, 1
      do Sp=0,1
      do Msp=-1,1
        if(hold(S,Ms,Sp,Msp).ne.hold(S,Ms,Sp,Msp)) then
        write(*,*) "FOUND IT"
        write(*,*) "S,Ms,Sp,Msp=", S,Ms,Sp,Msp 
        write(*,*) "hold(S,Ms,Sp,Msp)=", hold(S,Ms,Sp,Msp) 
        HoldNan=.true.
        end if
      end do
      end do
      end do
      end do
      end function HoldNan

      subroutine CalcKernel2BAsym(Kernel2B,
     &     factor,!qVec,pVec
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
c     INPUT/OUTPUT VARIABLES:
c     
      complex*16,intent(inout) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8,intent(in)  :: factor
      ! real*8,intent(in)  :: qVec(3), pVec(3)
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
      call singlesigmasym(hold,1.d0,0.d0,0.d0,Sp,S,verbosity)
      ! write(*,*) "In diag A O(q^3)"
      ! write(*,*) "hold=", hold 
      ! write(*,*) "eps=1,0,0" 
      ! write(*,*) "Sp=", Sp 
      ! write(*,*) "S=", S 
      ! write(*,*) ""
      ! write(*,*) ""
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
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine CalcKernel2BBsymVec(Kernel2B,
     &     factor, A,! A.σ
     &     B, !B.ε
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
c     INPUT/OUTPUT VARIABLE:
c     
      complex*16, intent(inout) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c      complex*16 Kernel2Bpx(0:1,-1:1,0:1,-1:1),Kernel2Bpy(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8,intent(in)  :: factor
c     real*8,intent(in)  :: Ax,Ay,Az,Bx,By,Bz
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
      call singlesigmasym(hold,A(1),A(2),A(3),Sp,S,verbosity)
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
      end

      subroutine cross(a, b,output)
        real*8, intent(in) :: a(3), b(3)
        real*8, intent(out):: output(3)
        output(1) = a(2) * b(3) - a(3) * b(2)
        output(2) = a(3) * b(1) - a(1) * b(3)
        output(3) = a(1) * b(2) - a(2) * b(1)
      end subroutine cross
      
      complex*16 function Ihold(Sp,Msp,S,Ms)
        integer Sp,Msp, S,Ms  
        integer delta
        Ihold = dcmplx(delta(Sp,S)*delta(Msp,Ms))
      end function

      integer function delta(a,b)
          integer, intent(in) :: a,b
          delta = merge(1, 0, a == b)
      end function delta
