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
      real*8 kSquare, mPion
      real*8 term1, term2, term3
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) "--------------------------------------------------------------------------------"
      write(*,*) "Kernel: Twobody Pion Pion at Threshold"
      write(*,*) "--------------------------------------------------------------------------------"
      write(*,*) "   Kernel Code Version 1.0"
      write(*,*) "      Alexander Long starting September 2025   "
      write(*,*)
      mPion=134.976d0
c     write(*,*) "DEBUG KernelGreeting: Egamma=", Egamma
c     write(*,*) "DEBUG KernelGreeting: Mnucl=", Mnucl
c     write(*,*) "DEBUG KernelGreeting: mPion=", mPion
      term1=Egamma*Egamma

      term2=-1*mPion**2*(Mnucl**2+Egamma**2-Egamma*sqrt(Mnucl**2+Egamma**2))
      term2=term2/Mnucl**2

      term3=(mPion**4)*(Mnucl**2+2*Egamma**2-2*Egamma*sqrt(Mnucl**2+Egamma**2))
      term3=term3/(4.d0*Mnucl**4)
      kSquare=term1+term2+term3
c     write(*,*) "kSquare=", kSquare 
      if ((kSquare.lt.0).and.(kSquare.gt.-10.d0)) then
        kSquare=0.d0
      end if
c     write(*,*) "DEBUG KernelGreeting: term1=", term1
c     write(*,*) "DEBUG KernelGreeting: term2=", term2
c     write(*,*) "DEBUG KernelGreeting: term3=", term3
c     write(*,*) "DEBUG KernelGreeting: kSquare=", kSquare
      Eprobe=sqrt(kSquare + mPion**2)
c     write(*,*) "DEBUG KernelGreeting: Eprobe=", Eprobe
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
      write(*,*) "Based on BKM review equation 5.30"
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

c     subroutine Calc2Bspinisospintrans(Kernel2B,ppVecs,Mnucl,
c    &     extQnumlimit,ml12,ml12p,
c    &     t12,mt12,t12p,mt12p,l12,s12,
c    &     l12p,s12p,thetacm,Eprobe,pVec,uVec,numDiagrams,verbosity)

      subroutine Calc2Bspinisospintrans(Kernel2B,ppVecs,Mnucl,
     &      extQnumlimit,ml12,ml12p,
     &      t12,mt12,t12p,mt12p,l12,
     &      s12,l12p,s12p,thetacm,Eprobe,pVec,
     &      uVec,numDiagrams,calctype,verbosity)
c     !Alex Long 2024:
c     pVec, is the  physical momenta, but uVec is the generic integration variable which may be transformed
c     taking uVec=ppVec works when there isn't a singularity.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc    Momentum vectors passed to this function are in MeV
cc    NOTE ON UNITS: hgrie Nov 2023
cc    Overall units here determine the units of the total twobody output ME Result().
cc            If kernel given in MeV^-n, then output ME will be given in MeV^(3-n). 
cc            Multiplying by powers of HC translates into output of final MEs "Result" into powers of fm.
cc            EXAMPLE: Compton has twobody kernel in MeV^-4 (n=4) ==> result ME in MeV^-1. To convert to fm, multiply by HC.
cc                     In Compton, that multiplication by HC isnot done in the fortran code, but later in the mathematica processing files.
cc            EXAMPLE: Pion Photoproduction kernel has units MeV^-2 if the output should be the twobody functions f_TL.
cc                     n=2 => Results() output in MeV^1. But F_TL output in fm^-1, so divide in kernel by HC to get fm^-1 units in Results().
cc            EXAMPLE: Pion scattering has units of MeV^-2 (recall f_pi) has units of MeV.
cc
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

      real*8,intent(in)  :: thetacm,Eprobe
      integer,intent(in) :: extQnumlimit, numDiagrams
      integer,intent(in) :: t12,mt12,t12p,mt12p,l12,l12p,s12,s12p, ml12,ml12p
      real*8, intent(in) :: pVec(3), uVec(3)
c!     real*8,intent(in)  :: px,py,pz,ppx,ppy,ppz

      integer,intent(in) :: calctype, verbosity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c!     LOCAL VARIABLES:
      integer diagNumber
c!     real*8 qVec(3)!, qpVec(3)
      real*8 ppVec(3)
c!     real*8 qVec(3)
      real*8 dl12by2
      complex*16 factorAsym!,factorBsym,factorB2,factorA2
      complex*16 factorAasy!,factorBasy
      real*8 kVec(3),kpVec(3),k
      real*8 mPion, Mnucl
      real*8 r,theta,phi
      complex*16 Yl12(-5:5)
      complex*16 KernelA(1:extQnumlimit,0:1,-1:1,0:1,-1:1) !contribution just from diagram A
      real*8 m1,m2,m3,m4
      real*8 sqrtS,prefactor
      real*8 ppVecA(3)
c!
      mPion=mpi0
      ppVecs=0.d0
      k=real(sqrt(Eprobe*Eprobe - mPion*mPion),8)
      kVec=(/0.d0,0.d0,real(k,8)/)
      diagNumber=0
      Kernel2B=c0
      kernelA=c0
      dl12by2=(l12-l12p)/2.d0   !to check if l12-l12p is  even or odd
c
      ppVec=uVec !just for the calculateqs call

c      !subroutine calculateqsmass is available for kpVec calculation
      m1=mNucl
      m2=mPion
      m3=mNucl
      m4=mPion

      call calculateqs2Mass(pVec,ppVec,kVec,kpVec,m1,m2,m3,m4,thetacm,verbosity)
      sqrtS=sqrt(mNucl*mNucl+k*k)+sqrt(mPion*mPion+k*k)
c     assign prefactor=8*Pi*sqrtS to convert from scattering length to matrix element
c     setting prefactor=1.d0=HC returns scattering length in fm
      prefactor=(1/(1+mPion/mNucl))*(1/(4*Pi))!matrix elements are unitless
      prefactor=prefactor*1000*mpi0
      diagNumber=1
      ppVecA=0.d0
      call getDiagABC(KernelA,pVec,uVec,ppVecA,kVec,kpVec,t12,t12p,
     &      mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,mNucl,verbosity)
      ppVecs(diagNumber,:)=ppVecA
c     KernelA computed in getDiagABC has units [MeV^-4]
c     Multiply by 1000*mpi0 -> units MeV^-3 -> units of `Result` "unitless" but actually are in (1000*mpi0)^-1 units
c     prefactor = prefactor*8*pi*sqrtS*HC
      Kernel2B(diagNumber,:,:,:,:,:)=prefactor*KernelA

      end




      subroutine getDiagABC(Kerneltmp,pVec,uVec,ppVecA,kVec,kpVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,mNucl,verbosity)
      implicit none
c     diagrams B and C combine to one diagram https://arxiv.org/abs/1003.3826v1
      include '../common-densities/constants.def'
c     Parameters-------------------------------------
      complex*16,intent(out) :: Kerneltmp(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c     integer, intent(out) :: diagNumber
      real*8, intent(out) :: ppVecA(3)
      integer,intent(in) :: extQnumlimit
      real*8, intent(in) :: pVec(3), kVec(3), uVec(3),kpVec(3), mNucl
      integer, intent(in) :: t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12, verbosity

c     Internal variables
      real*8 qVec(3), ppVec(3),ell(3), qpVec(3)!, r, theta, phi,
      complex*16 factorAsym, factorAasy
      complex*16 factorBCsym, factorBCasy

      logical useTransform
      real*8 Jacobian, prefactor
      real*8 vecsquare
c     write(*,*) "shape(Kerneltmp)=", shape(Kerneltmp)
c     write(*,*) "getDiagAB: extQnumlimit=", extQnumlimit, "s12p=", s12p, "s12=", s12
      ! write(*,*) "varsub-2Bkernel.PionPion.f:216 shape(Kerneltmp)=", shape(Kerneltmp) 
      if (.not.(all(ppVecA.eq.0))) then
          write(*,*) "ppVec assigned elsewhere, stopping"
          write(*,*) "In 2Bkernel.PionPhotoProdThresh.f: ppVecA=",ppVecA 
          stop
      end if

      useTransform=.true.
      if (useTransform) then
c       uVec=qVec=pVec-ppVec+kVec/2+kpVec/2
        ppVec=pVec-uVec+(kVec+kpVec)/2
        Jacobian=-1.d0
      else
          ppVec=uVec
          Jacobian=1.d0
      end if

      qVec=pVec-ppVec+kVec/2+kpVec/2 !=uVec
      qpVec=qVec-kVec

c     fpi=92.42 MeV defined in constants.def
c     from https://arxiv.org/abs/1003.3826v1 the factor (1/4pi) (1/(1+mpi/mNucl)) is put on the outside
c     
c
c     DIMENSIONAL ANALYSIS: prefactor = mpi0^2/(4*fpi^4) = [MeV^2]/[MeV^4] = [MeV^-2]
      prefactor = mpi0*mpi0/(4*(fpi**4.d0))

      if (DOT_PRODUCT(qVec,qVec).lt.1.d-10) then
         write(*,*) "DOT_PRODUCT(qVec,qVec)=", DOT_PRODUCT(qVec,qVec)
         stop
      end if

      factorAsym=prefactor*(1/(DOT_PRODUCT(qVec,qVec)))! units MeV^-4
c     factorBCsym=-1*gA*gA*prefactor*(1.d0/((DOT_PRODUCT(qVec,qVec)+mpi0**2)**2))
c      Fact
      factorBCsym=-1*gA*gA*prefactor*(1.d0/((DOT_PRODUCT(qVec,qVec)+mpi0**2)**2))!units MeV^-6, but CalcKernel2BBsym/asy adds factor MeV^2, gives result MeV^-4

      factorAasy=factorAsym
      factorBCasy=factorBCsym
      if (factorAsym.ne.factorAsym) then
        write(*,*) "factorAsym=", factorAsym 
        stop
      end if

      if (factorBCsym.ne.factorBCsym) then
        write(*,*) "factorBCsym=", factorBCsym 
        stop
      end if
c     if ((t12 .eq. t12p) .and. (mt12 .eq. mt12p)) then
         if (s12p .eq. s12) then ! spin symmetric part only; s12-s12p=0 => l12-l12p is even
            call CalcKernel2BAsym(Kerneltmp,
     &           factorAsym,
     &           s12p,s12,t12,mt12,t12p,mt12p,extQnumlimit,verbosity)

            if (any(Kerneltmp.ne.Kerneltmp)) then
              write(*,*) "Kernel was NaN after sym A stopping"
              error stop
            end if
            call CalcKernel2BBCsym(Kerneltmp,qVec,
     &           factorBCsym,
     &           s12p,s12,t12,mt12, t12p, mt12p, extQnumlimit,verbosity)

            if (any(Kerneltmp.ne.Kerneltmp)) then
              write(*,*) "Kernel was NaN after sym B stopping"
              error stop
            end if
         else                   !  spin anti-symmetric part only; s12 question: s12-s12p=±1 => l12-l12p is odd


            call CalcKernel2BAasy(Kerneltmp,
     &           factorAasy,
     &           s12p,s12,t12,mt12,t12p,mt12p,extQnumlimit,verbosity)

            if (any(Kerneltmp.ne.Kerneltmp)) then
              write(*,*) "Kernel was NaN after asym A stopping"
              error stop
            end if
            call CalcKernel2BBCasy(Kerneltmp,qVec,
     &           factorBCasy,
     &           s12p,s12,t12,mt12, t12p, mt12p, extQnumlimit,verbosity)

            if (any(Kerneltmp.ne.Kerneltmp)) then
              write(*,*) "Kernel was NaN after asym B stopping"
              error stop
            end if
         end if                 ! s12 question
c     else                      ! t12!=t12p
c       continue
c     end if                    !t12 question

      ppVecA=ppVec
      Kerneltmp=Kerneltmp*Jacobian
      if (any(Kerneltmp.ne.Kerneltmp)) then
        write(*,*) "Kernel was Nan, stopping"
        error stop
      end if
      end subroutine getDiagABC



      function vecsquare(vec)
          implicit none
          real*8 vec(3)
          real*8 vecsquare

          vecsquare=DOT_PRODUCT(vec,vec)
      end function vecsquare
