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
      subroutine KernelGreeting(verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer,intent(in) :: verbosity         ! verbosity index for stdout
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) "--------------------------------------------------------------------------------"
      write(*,*) "Kernel: Twobody Pion Photoproduction at Threshold"
      write(*,*) "--------------------------------------------------------------------------------"
      write(*,*) "   Kernel Code Version 1.0"
      write(*,*) "      Alexander Long/hgrie starting November 2023   "
      write(*,*)
      
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
c     &     Comp2Bx,Comp2By,Comp2Bpx,Comp2Bpy, ! for STUMP, see below
     &     extQnumlimit,ml12,ml12p,
     &     t12,mt12,t12p,mt12p,l12,s12,
     &     l12p,s12p,thetacm,k,pVec,uVec,calctype,numDiagrams,verbosity)
c     Alex Long 2024:
c     pVec, is the  physical momenta, but uVec is the generic integration variable which may be transformed

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NOTE ON UNITS: hgrie Nov 2023
c     Overall units here determine the units of the total twobody output ME Result().
c             If kernel given in MeV^-n, then output ME will be given in MeV^(3-n). 
c             Multiplying by powers of HC translates into output of final MEs "Result" into powers of fm.
c             EXAMPLE: Compton has twobody kernel in MeV^-4 (n=4) ==> result ME in MeV^-1. To convert to fm, multiply by HC.
c                      In Compton, that multiplication by HC isnot done in the fortran code, but later in the mathematica processing files.
c             EXAMPLE: Pion Photoproduction kernel has units MeV^-2 if the output should be the twobody functions f_TL.
c                      n=-2 => Results() output in MeV^1. But F_TL output in fm^-1, so divide here in kernel by HC to get fm^-1 units in Results().
c
c     (2π)³ is a factor of the twobody integration. TO INCLUDE IT OR NOT DEPENDS ON DEFINITIONS OF TWOBODY KERNELS!
c             In Compton, we insert it so that onebody and twobody Result() immediately hve same size and can be added: ottal=onebody+twobody. 
c             In pion photoproduction, it is part of the prefactor K2n of a diagram.
c             ==> If you want the twobody Result() output to be F_TL, you must un-compensate it here by *(2*Pi)**3.
c                 But if you want the twobody Result() output to be E_+ etc, so that you can simply add total=onebody+twobody,
c                 then the prefactor K2n shouldNOT contain the 1/(2π)³, i.e. multiply NOT with *(2*Pi)**3/HC, but with
c             K2n = sqrt(4*Pi*alpaEM)*gA*mpi**2/(16*Pi*fpi**3)*10**3 to get result() in the canonical units of 10^-3/mπplus.
c      
c     ==> Set your kernel up here so that your Result() has the desired units and factors of (2π)³. Do NOT make unit changes outside this file!
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CALLS:
c     calculateqs: to calculate momenta
c     CalcKernel2B...: to calculate kernel amplitudes. Much of the work is done in those routines via CalcKernel2B...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES:
      
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
      real*8,intent(in)  :: thetacm,k
      integer,intent(in) :: extQnumlimit, numDiagrams
      integer,intent(in) :: t12,mt12,t12p,mt12p,l12,l12p,s12,s12p, ml12,ml12p
      real*8, intent(in) :: pVec(3), uVec(3)
c     real*8,intent(in)  :: px,py,pz,ppx,ppy,ppz
               
      integer,intent(in) :: verbosity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
      integer diagNumber
      real*8 tmpVec(3)!, tmpVec2(3)
      real*8 ppVec(3)
      real*8 qVec(3)
      real*8 dl12by2
      complex*16 factorAsym!,factorBsym,factorB2,factorA2
      complex*16 factorAasy!,factorBasy
      real*8 kVec(3),kpVec(3)
      real*8 mPion, Mnucl
      real*8 r,theta,phi
      complex*16 Yl12(-5:5)
      complex*16 KernelA(1:extQnumlimit,0:1,-1:1,0:1,-1:1) !contribution just from diagram A
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Definitions of momenta repeated here for convenience
c     (All quantities in this comment to be read as vectors)
c
c     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Factors:
c      
c     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*************************************************************************************
c     
c     First a little initialization:
c     
      ppVecs=0.d0
      kVec=(/0.d0,0.d0,real(k,8)/)
      diagNumber=0
      Kernel2B=c0
      kernelA=c0
      dl12by2=(l12-l12p)/2.d0   !to check if l12-l12p is  even or odd
c     
c     Calculate momenta q,q',q':
c     
      
c      pVec=(/px,py,pz/)
       ppVec=uVec
       mPion=134.976d0

c      subroutine calculateqsmass is available for kpVec calculation
c      call calculateqsmass(pVec,ppVec,kVec,kpVec,thetacm,mPion,mNucl,verbosity)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Odelta0 2N contributions: NONE
c     <if they were nonzero, enter diagrams here>
      if (calctype.eq.Odelta0) return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Odelta2 2N contributions
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      diagNumber=1
c     call getDiagmpi(KernelA,pVec,uVec,ppVecs(diagNumber,:),kVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
      call getmpiEZsub(KernelA,pVec,uVec,ppVecs(diagNumber,:),kVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
c     call getDiagAB(KernelA,pVec,uVec,ppVecs(diagNumber,:),kVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
      Kernel2B(diagNumber,:,:,:,:,:)=KernelA

c     diagNumber=2
c     Kernel2B(diagNumber,:,:,:,:,:)=0.d0
c     ppVecs(diagNumber,:)=ppVecs(1,:)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end Odelta2 2N contributions
      if (calctype.eq.Odelta2) return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Odelta3 2N contributions: NONE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     <if they were nonzero, enter diagrams here>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end Odelta3 2N contributions
      if (calctype.eq.Odelta3) return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Odelta4 2N contributions: NONE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     <if they were nonzero, enter diagrams here>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end Odelta2 2N contributions
      if (calctype.eq.Odelta4) return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      if (verbosity.eq.1000) continue
      end


      subroutine getDiagAB(Kerneltmp,pVec,uVec,ppVecA,kVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
c     Diagram A and Diagram B actually have the same integration variable, so combine them
c     and use only "one" diagram in the input file numDiagrams=1

      implicit none

      include '../common-densities/constants.def'
c     Parameters-------------------------------------
      complex*16,intent(out) :: Kerneltmp(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c     integer, intent(out) :: diagNumber
      real*8, intent(out) :: ppVecA(3)
      integer,intent(in) :: extQnumlimit
      real*8, intent(in) :: pVec(3), kVec(3), uVec(3)
      integer, intent(in) :: t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12, verbosity

c     Internal variables
      real*8 tmpVec(3), ppVec(3),ell(3), tmpVec2(3)!, r, theta, phi,
      complex*16 factorAsym, factorAasy
      complex*16 factorBsym, factorBasy
      logical useTransform

      if (.not.(all(ppVecA.eq.0))) then
          write(*,*) "ppVec assigned elsewhere, stopping"
          write(*,*) "In 2Bkernel.PionPhotoProdThresh.f: ppVecA=",ppVecA 
          stop
      end if

      useTransform=.true.
      if (useTransform) then
c       uVec=pVec-ppVec+kVec/2!-> ppVec= pVec-uVec+kVec/2 -> jacobian on the integration gives a factor of -1
        ppVec=pVec-uVec+kVec/2
      end if

      if (.not.useTransform) then
          ppVec=uVec
      end if

      tmpVec=pVec-ppVec+(kVec/2)


      if (DOT_PRODUCT(tmpVec-uVec,tmpVec-uVec).ge.1e-5) then
        write(*,*) "tmpVec!=uVec"
        write(*,*) "In 2Bkernel.PionPhotoProdThresh.f: uVec=",uVec 
        write(*,*) "In 2Bkernel.PionPhotoProdThresh.f: tmpVec=",tmpVec 
        stop
      end if

      if (DOT_PRODUCT(tmpVec,tmpVec).le.0.001) then
          write(*,*) "In 2Bkernel"
          write(*,*) "DOT_PRODUCT(tmpVec,tmpVec).le.0.01 evaluated true, stopping"
          stop
      end if

      factorAsym=-(-1)**(t12)*(1.d0/(DOT_PRODUCT(tmpVec,tmpVec)))*(2*Pi)**3/HC

      factorBsym=+2*(-1)**(t12)*(1.d0/(
     &        DOT_PRODUCT(tmpVec,tmpVec)))*
     &        (1.d0/(DOT_PRODUCT(tmpVec2,tmpVec2)+mpi2))
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
c     
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

      end subroutine getDiagAB



      subroutine getDiagmpi(KernelA,pVec,uVec,ppVecA,kVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
      implicit none

      include '../common-densities/constants.def'
c     Parameters-------------------------------------
      complex*16,intent(out) :: KernelA(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c     integer, intent(out) :: diagNumber
      real*8, intent(out) :: ppVecA(3)
      integer,intent(in) :: extQnumlimit
      real*8, intent(in) :: pVec(3), kVec(3), uVec(3)
      integer, intent(in) :: t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12, verbosity

c     Internal variables
      real*8 tmpVec(3), ppVec(3)!, r, theta, phi
      complex*16 factorAsym, factorAasy
      logical useTransform
      real*8 Jacobian

      useTransform=.true.!For easy changes by the 2Bkernel "user", may want to change this in the final version for speed

      if (.not.(all(ppVecA.eq.0))) then
          write(*,*) "ppVec assigned elsewhere, stopping"
          write(*,*) "In 2Bkernel.PionPhotoProdThresh.f: ppVecA=",ppVecA 
          stop
      end if

c     Use this block if you actually want to use the transformation
      if (useTransform) then
c         u=pVec-ppVec+kVec/2 -> ppVec= pVec-uVec+kVec/2 -> jacobian on the integration gives a factor of -1?
c         -> tmpVec=uVec, denominator is just uVec
          ppVec=pVec-uVec+(kVec/2)
          Jacobian=1.d0
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This line below ppVec=uVec if you don't want to use the transformation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (.not.useTransform) then
          ppVec=uVec
          Jacobian=1.d0
      end if


      tmpVec=pVec-ppVec+(kVec/2)!=uVec in case of transform

      if (DOT_PRODUCT(tmpVec-uVec,tmpVec-uVec).ge.1e-6) then
          write(*,*) "DOT_PRODUCT(tmpVec-uVec,tmpVec-uVec).ge.1e-6 evaluated true, stopping"
          stop
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (DOT_PRODUCT(tmpVec,tmpVec).le.0.001) then
          write(*,*) "In 2Bkernel"
          write(*,*) "DOT_PRODUCT(tmpVec,tmpVec).le.0.001 evaluated true, stopping"
          write(*,*) "In 2Bkernel.PionPhotoProdThresh.f: DOT_PRODUCT(tmpVec,tmpVec)=", DOT_PRODUCT(tmpVec,tmpVec)
          stop
      end if

      factorAsym=-(-1)**(t12)*(1.d0/(mpi2+DOT_PRODUCT(uVec,uVec)))*(2*Pi)**3/HC
      factorAasy=factorAsym
      factorAsym=factorAsym*Jacobian
      factorAasy=factorAasy*Jacobian
      
      if ((t12 .eq. t12p) .and. (mt12 .eq. 0) .and.(mt12p .eq. 0)) then
         if (s12p .eq. s12) then ! s12-s12p=0 => l12-l12p is even; spin symmetric part only

            call CalcKernel2BAsym(KernelA,
     &           factorAsym,
     &           s12p,s12,extQnumlimit,verbosity)
         else                   ! s12 question: s12-s12p=±1 => l12-l12p is odd; spin anti-symmetric part only
            call CalcKernel2BAasy(KernelA,
     &           factorAasy,
     &           s12p,s12,extQnumlimit,verbosity)
         end if                 ! s12 question
      else                      ! t12!=t12p
         continue
c     diagrams (A/B) have no components with t12!=t12p. 
      end if                    !t12 question
      ppVecA=ppVec

      end subroutine getDiagmpi

      subroutine getmpiEZsub(KernelA,pVec,uVec,ppVecA,kVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
      implicit none
c     Uses the non-physical propogator mpi[ m_\pi^2 + (ppVec+kGamma/2)^2 + pVec^2]^(-2)
      include '../common-densities/constants.def'
c     Parameters-------------------------------------
      complex*16,intent(out) :: KernelA(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c     integer, intent(out) :: diagNumber
      real*8, intent(out) :: ppVecA(3)
      integer,intent(in) :: extQnumlimit
      real*8, intent(in) :: pVec(3), kVec(3), uVec(3)
      integer, intent(in) :: t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12, verbosity

c     Internal variables
      real*8 tmpVec(3), ppVec(3)!, r, theta, phi
      complex*16 factorAsym, factorAasy
      logical useTransform
      real*8 Jacobian
c     For easy changes by the 2Bkernel "user", may want to change this in the final version for speed, since "if"
c     statements are slow
      useTransform=.true.

      if (.not.(all(ppVecA.eq.0))) then
          write(*,*) "ppVec assigned elsewhere, stopping"
          write(*,*) "In 2Bkernel.PionPhotoProdThresh.f: ppVecA=",ppVecA 
          stop
      end if

c     Use this block if you actually want to use the transformation
      if (useTransform) then
c         uVec=ppVec+kVec/2 -> ppVec= uVec-kVec/2 
          ppVec=uVec-kVec/2
          Jacobian=1.d0
c         So later when we set tmpVec=ppVec+kVec/2  then
c         tmpVec=uVec-kVec/2+kVec/2=uVec
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This line below ppVec=uVec if you don't want to use the transformation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (.not.useTransform) then
          ppVec=uVec
          Jacobian=1.d0
      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      tmpVec=ppVec+(kVec/2.d0)
      factorAsym=-(-1)**(t12)*mpi*
     &    (mpi2+DOT_PRODUCT(tmpVec,tmpVec)+DOT_PRODUCT(pVec,pVec))**(-2.d0)
     &      *(2*Pi)**3/HC
      factorAasy=factorAsym
      factorAsym=factorAsym*Jacobian
      factorAasy=factorAasy*Jacobian
      
      if ((t12 .eq. t12p) .and. (mt12 .eq. 0) .and.(mt12p .eq. 0)) then
         if (s12p .eq. s12) then ! s12-s12p=0 => l12-l12p is even; spin symmetric part only

            call CalcKernel2BAsym(KernelA,
     &           factorAsym,
     &           s12p,s12,extQnumlimit,verbosity)
         else                   ! s12 question: s12-s12p=±1 => l12-l12p is odd; spin anti-symmetric part only
            call CalcKernel2BAasy(KernelA,
     &           factorAasy,
     &           s12p,s12,extQnumlimit,verbosity)
         end if                 ! s12 question
      else                      ! t12!=t12p
         continue
c     diagrams (A/B) have no components with t12!=t12p. 
      end if                    !t12 question
      ppVecA=ppVec

      end subroutine getmpiEZsub
