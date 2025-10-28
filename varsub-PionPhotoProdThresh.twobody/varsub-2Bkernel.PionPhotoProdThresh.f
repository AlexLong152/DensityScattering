cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of KERNEL code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              KernelGreeting         : output to stdout with message which process computed and its version number
c              KernelFarewell         : output to stdout with message which described computed quantity,
c                                       symmetry/-ies used and mapping of extQnum
c              Calc2Bspinisospintrans : compute total kernel by calling all diagrams up to given order
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, loosely based on 2Bspinisospintrans.f of Compton density code v2.0 hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c      
c     Organisation of orders:
c     First Odelta0 computation -- terminates after that if not more needed.
c     Else moves on to Odelta2 computation -- terminates after that if not more needed.
c     Else moves on to Odelta3 computation -- etc.
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie Nov 2023: show kernel process and version
c     included here since will change when kernel changes
      subroutine KernelGreeting(Egamma,EProbe,Mnucl,verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8, intent(in) :: Egamma, Mnucl
      real*8, intent(out) :: Eprobe
      integer,intent(in) :: verbosity         ! verbosity index for stdout
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) "--------------------------------------------------------------------------------"
      write(*,*) "Kernel: Twobody Pion Photoproduction at Threshold"
      write(*,*) "--------------------------------------------------------------------------------"
      write(*,*) "   Kernel Code Version 1.0"
      write(*,*) "      Alexander Long/hgrie starting November 2023   "
      write(*,*)
      Eprobe=Egamma !in pion photoproduction, incoming photon energy = probe energy, not the case in pion scattering
      if (verbosity.eq.1000) continue
      end
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie Nov 2023: show kernel process and version
c     included here since will change when kernel changes
      subroutine KernelFarewell(extQnumlimit,symmetry,verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer,intent(in) :: extQnumlimit,symmetry
      integer,intent(in) :: verbosity         ! verbosity index for stdout
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     show what is outputted, and its units
      write(*,*) "Output: F_{L/T} ε.S^Mprime_M = < Mprime | kernel | M > Lenkewitz Eur. Phys. J. A (2013) 49:20 eq (10)"
      write(*,*) "        with ε: incoming-photon polarisation, S: nucleus spin; [F_{L/T}]=[fm]¯¹"
      write(*,*) "        Mapping of extQnum: 1 = ε_x, 2 = ε_y (both transversal); 3 = ε_z (longitudinal)"
c     characterise symmetry/-ies, if any used.
      If (symmetry.eq.0) then
         write(*,*) "        No symmetries used."
c     following is a STUMP/SKETCH what symmetry outoput could be
c     else if (symmetry.eq.1) then
c         write(*,*) "        Symmetry imposed: ME(extQnum=1) =  ME(extQnum=1) up to sign." ! better specify signs!!
      end if
      write(*,*)
      
      if (extQnumlimit.eq.1000) continue
      if (verbosity.eq.1000) continue
      end
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Calc2Bspinisospintrans(Kernel2B,ppVecs,Mnucl,
     &     extQnumlimit,ml12,ml12p,
     &     t12,mt12,t12p,mt12p,l12,s12,
     &     l12p,s12p,thetacm,Eprobe,pVec,uVec,calctype,numDiagrams,verbosity)
c     !Alex Long 2024:
c     !pVec, is the  physical momenta, but uVec is the generic integration variable which may be transformed

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc    NOTE ON UNITS: hgrie Nov 2023
cc    Overall units here determine the units of the total twobody output ME Result().
cc            If kernel given in MeV^-n, then output ME will be given in MeV^(3-n). 
cc            Multiplying by powers of HC translates into output of final MEs "Result" into powers of fm.
cc            EXAMPLE: Compton has twobody kernel in MeV^-4 (n=4) ==> result ME in MeV^-1. To convert to fm, multiply by HC.
cc                     In Compton, that multiplication by HC isnot done in the fortran code, but later in the mathematica processing files.
cc            EXAMPLE: Pion Photoproduction kernel has units MeV^-2 if the output should be the twobody functions f_TL.
cc                     n=-2 => Results() output in MeV^1. But F_TL output in fm^-1, so divide here in kernel by HC to get fm^-1 units in Results().
cc
cc    (2π)³ is a factor of the twobody integration. TO INCLUDE IT OR NOT DEPENDS ON DEFINITIONS OF TWOBODY KERNELS!
cc            In Compton, we insert it so that onebody and twobody Result() immediately hve same size and can be added: ottal=onebody+twobody. 
cc            In pion photoproduction, it is part of the prefactor K2n of a diagram.
cc            ==> If you want the twobody Result() output to be F_TL, you must un-compensate it here by *(2*Pi)**3.
cc                But if you want the twobody Result() output to be E_+ etc, so that you can simply add total=onebody+twobody,
cc                then the prefactor K2n shouldNOT contain the 1/(2π)³, i.e. multiply NOT with *(2*Pi)**3/HC, but with
cc            K2n = sqrt(4*Pi*alpaEM)*gA*mpi**2/(16*Pi*fpi**3)*10**3 to get result() in the canonical units of 10^-3/mπplus.
cc     
cc    ==> Set your kernel up here so that your Result() has the desired units and factors of (2π)³. Do NOT make unit changes outside this file!
cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc    !CALLS:
cc    calculateqs: to calculate momenta
cc    CalcKernel2B...: to calculate kernel amplitudes. Much of the work is done in those routines via CalcKernel2B...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     !OUTPUT VARIABLES:
      
      complex*16,intent(out) :: Kernel2B(1:numDiagrams,1:extQnumlimit,0:1,-1:1,0:1,-1:1) ! was Comp2Bxx/xy/yx/yy
      real*8, intent(out) :: ppVecs(1:numDiagrams,1:3)

c     Note that Kernel2B.. computes the amplitude for extQnums
c     Indices: 1st: extQnum
c              2nd: NN spin of final state: S=0 or S=1 (s12p)
c              3rd: NN spin projection of final state ms12p      
c              4th: NN spin of initial state: S=0 or S=1 (s12)
c              5th: NN spin projection of initial state ms12      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
      
      integer,intent(in) :: calctype
      real*8,intent(in)  :: thetacm,Eprobe
      integer,intent(in) :: extQnumlimit, numDiagrams
      integer,intent(in) :: t12,mt12,t12p,mt12p,l12,l12p,s12,s12p, ml12,ml12p
      real*8, intent(in) :: pVec(3), uVec(3)
c!     real*8,intent(in)  :: px,py,pz,ppx,ppy,ppz
               
      integer,intent(in) :: verbosity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c!     LOCAL VARIABLES:
      integer diagNumber
c!     real*8 qVec(3)!, qpVec(3)
      real*8 ppVec(3)
c!     real*8 qVec(3)
      real*8 dl12by2
      complex*16 factorAsym!,factorBsym,factorB2,factorA2
      complex*16 factorAasy!,factorBasy
      real*8 kVec(3),kpVec(3)
      real*8 mPion, Mnucl
      real*8 r,theta,phi
      complex*16 Yl12(-5:5)
      complex*16 KernelA(1:extQnumlimit,0:1,-1:1,0:1,-1:1) !contribution just from diagram A
      complex*16 KernelStatic(1:extQnumlimit,0:1,-1:1,0:1,-1:1) !contribution just from diagram A
      integer diagNum
c!     
c     First a little initialization:
c     
      ppVecs=0.d0
      kVec=(/0.d0,0.d0,real(Eprobe,8)/)!pion photoproduction, Eprobe=omega=momentum since photon has no mass
      diagNumber=0
      Kernel2B=c0
      kernelA=c0
      dl12by2=(l12-l12p)/2.d0   !to check if l12-l12p is  even or odd
c     
c     
      
c      pVec=(/px,py,pz/)
       ppVec=uVec!just for the calculateqs call
       ! mPion=134.976d0

c      !subroutine calculateqsmass is available for kpVec calculation
c      Calculate momenta q,q',q':
       call calculateqsmass(pVec,ppVec,kVec,kpVec,thetacm,mpi0,mNucl,verbosity)
c      kpVec=(/0.d0,0.d0,0.d0/) !need to read in more precison in input file energy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Odelta0 2N contributions: NONE
c     !<if they were nonzero, enter diagrams here>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Odelta2 2N contributions
cccccccccccccccccccccccTrue
      diagNumber=1
      call getDiagAB(KernelA,pVec,uVec,ppVecs(diagNumber,:),kVec,kpVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
      !getDiagABfinite is not completed at the moment, 
c     call getDiagABfinite(KernelA,pVec,uVec,ppVecs(diagNumber,:),kVec,kpVec,ppVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)

      ! if (calctype.eq.OQ3) then
      !    return
      ! end if
      Kernel2B(diagNumber,:,:,:,:,:)=KernelA
      if (calctype.eq.Odelta2) return

      !StaticDiags O(q^4) uses some variable substitution as diagAB
      !would need to reassign ppVecs(diagNumber,:), and diagNumber if this wasn't the case
      call getStaticDiags(KernelStatic,pVec,uVec,kVec,kpVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
      Kernel2B(diagNumber,:,:,:,:,:)=Kernel2B(diagNumber,:,:,:,:,:)+KernelStatic*(-1/(3.d0*mNucl*mpi))

      if (calctype.eq.Odelta4) return 
      end


      subroutine getStaticDiags(Kerneltmp,pVec,uVec,kVec,kpVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
      implicit none
      !in language of lenkewitz thesis, what he calls q' we call q, and what he calls q'' we call q'
      include '../common-densities/constants.def'
c     Parameters-------------------------------------
      complex*16,intent(out) :: Kerneltmp(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c     integer, intent(out) :: diagNumber
      ! real*8, intent(out) :: ppVecA(3)
      integer,intent(in) :: extQnumlimit
      real*8, intent(in) :: pVec(3), kVec(3), uVec(3),kpVec(3)
      integer, intent(in) :: t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12, verbosity

c     Internal variables
      real*8 qVec(3), ppVec(3),ell(3), qpVec(3)!, r, theta, phi,
      complex*16 factorAsym, factorAasy
      complex*16 factorBsym, factorBasy
      complex*16 factorCsym, factorCasy
      complex*16 factorDsym, factorDasy
      complex*16 factorEsym, factorEasy
      logical useTransform
      real*8 Jacobian, unitsFact, factor
      real*8 vecsquare
      external vecsquare
      logical Nancheck

      Kerneltmp=c0
      useTransform=.true.
      if (useTransform) then
c       uVec=pVec-ppVec+kVec/2!-> ppVec= pVec-uVec+kVec/2 -> jacobian on the integration gives a factor of -1
        ppVec=pVec-uVec+kVec/2
        Jacobian=-1.d0!TODO: check if this is needed
      else
          ppVec=uVec
          Jacobian=1.d0
      end if

      qVec=pVec-ppVec+(kVec/2)!qVec=uVec with the substitution
      qpVec=pVec-ppVec-(kVec/2)


      if (DOT_PRODUCT(qVec,qVec).le.0.001) then
          write(*,*) "In 2Bkernel"
          write(*,*) "DOT_PRODUCT(qVec,qVec).le.0.001 evaluated true, stopping"
          stop
      end if

      unitsFact=(2*Pi)**3/HC
      factor=-1*(-1)**(t12)*unitsFact

      factorAsym=factor*(1-2.d0*ga*ga)*(1.d0+2.d0*(dot_product(qVec,ppVec)/vecsquare(qVec)))
      factorBsym=2.d0*factor* (1/vecsquare(qVec))
      factorCsym=2.d0*factor*(1/(vecsquare(qVec)+mpi*mpi))

      factorDsym = factor * (2.d0*ga*ga) * (1.d0 / (vecsquare(qVec) + mpi*mpi))

      factorEsym=factor*(-1.d0+2*ga*ga)*(1.d0+2.d0*(dot_product(qVec,pVec)/vecsquare(qVec)))
     & *(1/(vecsquare(qpVec)+mpi**2))

      factorAasy=factorAsym
      factorBasy=factorBsym
      factorCasy=factorCsym
      factorDasy=factorDsym
      factorEasy=factorEsym
      ! write(*,*) "In static"
c     write(*,*) "1-2.d0*ga*ga=", 1-2.d0*ga*ga 
c     write(*,*) "qVec=", qVec 
c     write(*,*) "ppVec=", ppVec 
c     write(*,*) "dot_product(qVec,ppVec)=", dot_product(qVec,ppVec) 
c     write(*,*) "vecsquare(qVec)=", vecsquare(qVec) 
c     write(*,*) "dot_product(qVec,ppVec)/vecsquare(qVec)=", dot_product(qVec,ppVec)/vecsquare(qVec) 
c     write(*,*) "(1-2.d0*ga*ga)*(1.d0+2.d0*(dot_product(qVec,ppVec)/vecsquare(qVec)))",
c    &  (1-2.d0*ga*ga)*(1.d0+2.d0*(dot_product(qVec,ppVec)/vecsquare(qVec)))
c     write(*,*) "factorAsym=", factorAsym !Issue with factorAsym maybe
c     write(*,*) "factorBsym=", factorBsym 
c     write(*,*) "factorCsym=", factorCsym 
c     write(*,*) "factorDsym=", factorDsym 
c     write(*,*) ""
c     write(*,*) ""

      if ((t12 .eq. t12p) .and. (mt12 .eq. 0) .and.(mt12p .eq. 0)) then
         if (s12p .eq. s12) then ! s12-s12p=0 => l12-l12p is even; spin symmetric part only

            call StaticKernelAsym(Kerneltmp,
     &           factorAsym,!qVec,pVec,
     &           s12p,s12,extQnumlimit,verbosity)

            if(Nancheck(Kerneltmp)) then
              write(*,*) "NaN on Diag Asym"
              stop
            end if

            call StaticKernelBsym(Kerneltmp,
     &           factorBsym, qVec,
     &           pVec,kVec,
     &           s12p,s12,extQnumlimit,verbosity)

            if(Nancheck(Kerneltmp)) then
              write(*,*) "factorBsym=", factorBsym 
              write(*,*) "Kerneltmp=", Kerneltmp 
              write(*,*) "NaN on Diag Bsym"
              stop
            end if
            call StaticKernelCsym(Kerneltmp,
     &           factorCsym, qVec,
     &           ppVec,kVec,
     &           s12p,s12,extQnumlimit,verbosity)

            if(Nancheck(Kerneltmp)) then
              write(*,*) "NaN on Diag Csym"
              stop
            end if
            call StaticKernelDsym(Kerneltmp,
     &           factorDsym, qVec,
     &           ppVec,kVec,
     &           s12p,s12,extQnumlimit,verbosity)

            if(Nancheck(Kerneltmp)) then
              write(*,*) "NaN on Diag Dsym"
              stop
            end if

            call StaticKernelEsym(Kerneltmp,
     &           factorEsym, qVec,
     &           qpVec,pVec,
     &           s12p,s12,extQnumlimit,verbosity)

            if(Nancheck(Kerneltmp)) then
              write(*,*) "NaN on Diag Dsym"
              stop
            end if

         else                   ! s12 question: s12-s12p=±1 => l12-l12p is odd; spin anti-symmetric part only

              call StaticKernelAasym(Kerneltmp,
     &           factorAsym,!qVec,pVec,
     &           s12p,s12,extQnumlimit,verbosity)
      
              if(Nancheck(Kerneltmp)) then
                write(*,*) "NaN on Diag A asym"
                stop
              end if
              call StaticKernelBasym(Kerneltmp,
     &           factorBsym, qVec,
     &           pVec,kVec,
     &           s12p,s12,extQnumlimit,verbosity)
      
              if(Nancheck(Kerneltmp)) then
                write(*,*) "NaN on Diag B asym"
                stop
              end if
              call StaticKernelCasym(Kerneltmp,
     &           factorCasy, qVec,
     &           ppVec,kVec,
     &           s12p,s12,extQnumlimit,verbosity)

              if(Nancheck(Kerneltmp)) then
                write(*,*) "NaN on Diag C asym"
                stop
              end if
              call StaticKernelDasym(Kerneltmp,
     &           factorDasy, qVec,
     &           ppVec,kVec,
     &           s12p,s12,extQnumlimit,verbosity)
      
              if(Nancheck(Kerneltmp)) then
                write(*,*) "NaN on Diag D asym"
                stop
              end if

            call StaticKernelEasym(Kerneltmp,
     &           factorEsym, qVec,
     &           qpVec,pVec,
     &           s12p,s12,extQnumlimit,verbosity)

              if(Nancheck(Kerneltmp)) then
                write(*,*) "NaN on Diag E asym"
                stop
              end if
         end if                 ! s12 question
      else                      ! t12!=t12p
         continue
c     diagrams (A/B) have no components with t12!=t12p. 
      end if                    !t12 question
      end subroutine getStaticDiags



      logical function NanCheck(Kerneltmp)

      implicit none
      complex*16,intent(in) :: Kerneltmp(1:3,0:1,-1:1,0:1,-1:1)
      integer i,S,Sp,Ms,Msp

      NanCheck=.false.
      do i = 1, 3
      do Sp=0,1
      do Msp=-1,1
      do S=0, 1
      do Ms = -1, 1
        if(Kerneltmp(i,Sp,Msp,S,Ms).ne.Kerneltmp(i,Sp,Msp,S,Ms)) then
        NanCheck=.true.
        end if
      end do
      end do
      end do
      end do
      end do
      end function NanCheck

      subroutine getDiagAB(Kerneltmp,pVec,uVec,ppVecA,kVec,kpVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
c     Threshold pion photoproduction - first two diagrams
c     Diagram A and Diagram B actually have the same integration variable, so combine them
c     and use only "one" diagram in the input file numDiagrams=1

      implicit none

      include '../common-densities/constants.def'
c     Parameters-------------------------------------
      complex*16,intent(out) :: Kerneltmp(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c     integer, intent(out) :: diagNumber
      real*8, intent(out) :: ppVecA(3)
      integer,intent(in) :: extQnumlimit
      real*8, intent(in) :: pVec(3), kVec(3), uVec(3),kpVec(3)
      integer, intent(in) :: t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12, verbosity

c     Internal variables
      real*8 qVec(3), ppVec(3),ell(3), qpVec(3)!, r, theta, phi,
      complex*16 factorAsym, factorAasy
      complex*16 factorBsym, factorBasy
      logical useTransform
      real*8 Jacobian, prefactor

      if (.not.(all(ppVecA.eq.0))) then
          write(*,*) "ppVec assigned elsewhere, stopping"
          write(*,*) "In 2Bkernel.PionPhotoProdThresh.f: ppVecA=",ppVecA 
          stop
      end if

      useTransform=.true.
      if (useTransform) then
c       uVec=pVec-ppVec+kVec/2!-> ppVec= pVec-uVec+kVec/2 -> jacobian on the integration gives a factor of -1
        ppVec=pVec-uVec+kVec/2
        Jacobian=-1.d0!TODO: check if this is needed
      else
          ppVec=uVec
          Jacobian=1.d0
      end if
      qVec=pVec-ppVec+(kVec/2)!qVec=uVec with the substitution
      qpVec=pVec-ppVec-(kVec/2)


      if (DOT_PRODUCT(qVec,qVec).le.0.001) then
          write(*,*) "In 2Bkernel"
          write(*,*) "DOT_PRODUCT(qVec,qVec).le.0.001 evaluated true, stopping"
          stop
      end if

c     Sign difference comes from convention. In Lenkewitz paper the result is 
c     F^a-F^b (Table 2 in Lenkewitz 2011) whereas we define the result to be 
c     the addition of quantities, so for us its F^a+F^b

c     For imlicit cancelation
      factorAsym=-1*(-1)**(t12)*(1.d0/(DOT_PRODUCT(qVec,qVec)))*(2*Pi)**3/HC
c     factorBsym=+2*(-1)**(t12)*(1.d0/(
c    &        DOT_PRODUCT(qVec,qVec)))*
c    &        (1.d0/(DOT_PRODUCT(qpVec,qpVec)+mpi2))
c    &     *(2*Pi)**3/HC

      factorBsym=2*(-1)**(t12)*(1.d0/(
     &        DOT_PRODUCT(qVec,qVec)))*
     &        (1.d0/(DOT_PRODUCT(qpVec,qpVec)+mpi0**2))
     &     *(2*Pi)**3/HC

      factorAasy=factorAsym
      factorBasy=factorBsym

      if ((t12 .eq. t12p) .and. (mt12 .eq. 0) .and.(mt12p .eq. 0)) then
         if (s12p .eq. s12) then ! s12-s12p=0 => l12-l12p is even; spin symmetric part only

            call CalcKernel2BAsym(Kerneltmp,
     &           factorAsym,
     &           s12p,s12,extQnumlimit,verbosity)

            call CalcKernel2BBsymVec(Kerneltmp,
     &           factorBsym,
     &           pVec-ppVec-(kVec/2), ! preceding is vector dotted with σ
     &           pVec-ppVec, ! preceding is vector dotted with ε
     &           s12p,s12,extQnumlimit,verbosity)
         else                   ! s12 question: s12-s12p=±1 => l12-l12p is odd; spin anti-symmetric part only

            call CalcKernel2BAasy(Kerneltmp,
     &           factorAasy,
     &           s12p,s12,extQnumlimit,verbosity)

            call CalcKernel2BBasyVec(Kerneltmp,
     &           factorBasy,
     &           pVec-ppVec-(kVec/2), ! preceding is vector dotted with σ
     &           pVec-ppVec, ! preceding is vector dotted with ε
     &           s12p,s12,extQnumlimit,verbosity)


         end if                 ! s12 question
      else                      ! t12!=t12p
         continue
c     diagrams (A/B) have no components with t12!=t12p. 
      end if                    !t12 question
      ppVecA=ppVec
c     write(*,*) "In varsub-2Bkernel: ppVecA=",ppVecA 
c     write(*,*) ""
      Kerneltmp=Kerneltmp*Jacobian
      end subroutine getDiagAB



      function vecsquare(vec)
          implicit none
          real*8 vec(3)
          real*8 vecsquare

          vecsquare=DOT_PRODUCT(vec,vec)
      end function vecsquare


      subroutine getDiagABfinite(Kerneltmp,pVec,uVec,ppVecA,kVec,kpVec,ppVec,
     &    t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
c     !Diagram A and Diagram B actually have the same integration variable, so combine them
c     !and use only "one" diagram in the input file numDiagrams=1

      implicit none

      include '../common-densities/constants.def'
c     Parameters-------------------------------------
      complex*16,intent(out) :: Kerneltmp(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c     integer, intent(out) :: diagNumber
      real*8, intent(out) :: ppVecA(3)
      integer,intent(in) :: extQnumlimit
      real*8, intent(in) :: pVec(3), kVec(3), uVec(3),kpVec(3)
      integer, intent(in) :: t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12, verbosity

c     Internal variables
      real*8 qVec(3), ppVec(3),ell(3), qpVec(3)!, r, theta, phi,
      complex*16 factorAsym, factorAasy
      complex*16 factorBsym, factorBasy
      logical useTransform
      real*8 Jacobian, prefactor
      real*8 q0,kp0, vecsquare,qp0,omega

c     if (.not.(all(ppVecA.eq.0))) then
c         write(*,*) "ppVec assigned elsewhere, stopping"
c         write(*,*) "In 2Bkernel.PionPhotoProdThresh.f: ppVecA=",ppVecA 
c         stop
c     end if

      useTransform=.true.
      if (useTransform) then
c       uVec=pVec-ppVec+kVec/2!-> ppVec= pVec-uVec+kVec/2 -> jacobian on the integration gives a factor of -1
        ppVec=pVec-uVec+kVec/2 !TODO: I think this should have +kpVec/2 in it as well for above threshold, test this
        Jacobian=-1.d0
      else 
          ppVec=uVec
          Jacobian=1.d0
      end if

      qVec=pVec-ppVec+(kVec/2)!qVec=uVec with the substitution
      qpVec=pVec-ppVec-(kVec/2)

      if (DOT_PRODUCT(qVec,qVec).le.0.001) then
          write(*,*) "In 2Bkernel"
          write(*,*) "DOT_PRODUCT(qVec,qVec).le.0.001 evaluated true, stopping"
          stop
      end if


      omega=sqrt(vecsquare(kVec))! omega=k0

      prefactor=(-1)**(t12)*(2*Pi)**3/(HC)
c   In the threshold result q0+kp0 cancels with 2*mpi0, but we need it here above threshold
      prefactor=prefactor/(2*mpi0)

c     The lines below assigning q0=mpi0, qp0=0, and kp0=mpi0 recover the threshold case
c     exactly as it appears in lenkewitz 2011
      q0=mpi0
      qp0=0
      kp0=mpi0

c     The lines below assigning q0, qp0 and kp0 recover higher order correctionsto lenkewitz 
      q0=omega+1/(2*Mnucleon) *(
     &  vecsquare(pVec-kVec/2)-vecsquare(ppVec-kpVec/2))
      qp0=1/(2*Mnucleon)*(vecsquare(pVec-kVec/2)-vecsquare(ppVec-kpVec/2))
      qp0=q0-omega
      kp0=sqrt(mpi0**2+vecsquare(kpVec))
      
c     write(*,*) "qp0=", qp0 , "close to 0?"

      factorAsym=(q0+kp0)/(q0**2-vecsquare(qVec)-mpi0**2)
      factorBsym=factorAsym*(1/(qp0**2 -vecsquare(qpVec)-mpi0**2))

      factorAsym=factorAsym*prefactor
      factorBsym=factorBsym*prefactor
      factorAasy=factorAsym
      factorBasy=factorBsym
      
      if ((t12 .eq. t12p) .and. (mt12 .eq. 0) .and.(mt12p .eq. 0)) then
         if (s12p .eq. s12) then ! s12-s12p=0 => l12-l12p is even; spin symmetric part only

            call CalcKernel2BAsym(Kerneltmp,
     &           factorAsym,
     &           s12p,s12,extQnumlimit,verbosity)

            call CalcKernel2BBsymVec(Kerneltmp,
     &           factorBsym,
     &           qpVec, ! preceding is vector dotted with σ
     &           qVec+qpVec, ! preceding is vector dotted with ε
     &           s12p,s12,extQnumlimit,verbosity)
         else                   ! s12 question: s12-s12p=±1 => l12-l12p is odd; spin anti-symmetric part only
c     
            call CalcKernel2BAasy(Kerneltmp,
     &           factorAasy,
     &           s12p,s12,extQnumlimit,verbosity)

            call CalcKernel2BBasyVec(Kerneltmp,
     &           factorBasy,
     &           qpVec, ! preceding is vector dotted with σ
     &           qVec+qpVec, ! preceding is vector dotted with ε
     &           s12p,s12,extQnumlimit,verbosity)


         end if                 ! s12 question
      else                      ! t12!=t12p
         continue
c     diagrams (A/B) have no components with t12!=t12p. 
      end if                    !t12 question
      ppVecA=ppVec
      Kerneltmp=Kerneltmp*Jacobian
      end subroutine getDiagABfinite
