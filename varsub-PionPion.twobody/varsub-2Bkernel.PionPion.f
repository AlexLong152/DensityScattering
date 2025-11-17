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
      real*8 kVec(3),kpVec(3),k
      real*8 mPion, Mnucl
      real*8 r,theta,phi
      complex*16 Yl12(-5:5)
      complex*16 KernelA(1:extQnumlimit,0:1,-1:1,0:1,-1:1) !contribution just from diagram A
      real*8 m1,m2,m3,m4
      integer diagNum
      real*8 sqrtS,prefactor
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
      sqrtS=sqrt(mNucl*mNucl+k*k)+sqrt(mpion*mpion+k*k)
      prefactor=8*Pi*sqrtS
      ! write(*,*) "varsub-2Bkernel.PionPion.f:184 extQnumlimit=", extQnumlimit 
      diagNumber=1
      call getDiagAB(KernelA,pVec,uVec,ppVecs(diagNumber,:),kVec,kpVec,t12,t12p,
     &      mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
c     KernelA is in units MeV^-4, prefactor is units MeV -> Kernel2B in units MeV^-3 -> Result in units MeV^0
      Kernel2B(diagNumber,:,:,:,:,:)=prefactor*KernelA

      end




      subroutine getDiagAB(Kerneltmp,pVec,uVec,ppVecA,kVec,kpVec,t12,t12p,mt12,mt12p,l12p,ml12p,s12p,s12,extQnumlimit,verbosity)
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
      complex*16 factorCsym, factorCasy

      logical useTransform
      real*8 Jacobian, prefactor,reducedMass
      real*8 vecsquare
c     write(*,*) "shape(Kerneltmp)=", shape(Kerneltmp)
c     write(*,*) "getDiagAB: extQnumlimit=", extQnumlimit, "s12p=", s12p, "s12=", s12
      ! write(*,*) "varsub-2Bkernel.PionPion.f:216 shape(Kerneltmp)=", shape(Kerneltmp) 
      if (.not.(all(ppVecA.eq.0))) then
          write(*,*) "ppVec assigned elsewhere, stopping"
          write(*,*) "In 2Bkernel.PionPhotoProdThresh.f: ppVecA=",ppVecA 
          stop
      end if

c     write(*,*) "DEBUG INPUT: pVec=", pVec
c     write(*,*) "DEBUG INPUT: uVec=", uVec
c     write(*,*) "DEBUG INPUT: kVec=", kVec
c     write(*,*) "DEBUG INPUT: kpVec=", kpVec
      useTransform=.true.
      if (useTransform) then
c       uVec=qVec=pVec-ppVec+kVec/2+kpVec/2
        ppVec=pVec-uVec+(kVec+kpVec)/2
c       write(*,*) "DEBUG: ppVec=", ppVec
        Jacobian=1.d0!check this
      else
          ppVec=uVec
          Jacobian=1.d0
      end if

      qVec=pVec-ppVec+kVec/2+kpVec/2 !=uVec
      qpVec=qVec-kVec
c     if (vecsquare(qVec-uVec).ge.1.0) then
c         write(*,*) "Error in getDiagAB: qVec-uVec too large:",vecsquare(qVec-uVec)
c         write(*,*) "qVec=", qVec
c         write(*,*) "uVec=", uVec
c         write(*,*) "kinemetics inconsistent,  stopping"
c         stop
c     end if

c     write(*,*) "DEBUG: qVec=", qVec
c     write(*,*) "DEBUG: DOT_PRODUCT(qVec,qVec)=", DOT_PRODUCT(qVec,qVec)
c     write(*,*) "DEBUG: |qVec|=", sqrt(DOT_PRODUCT(qVec,qVec))

      ! q0=(mPion**2 + DOT_PRODUCT(qVec,qVec))
      ! Epi=(mPion**2 + DOT_PRODUCT(kpVec,kpVec))

c     fpi=92.42 defined in constants.def
c     Define from BKM review 
      reducedMass=mpi0/mNucleon
c     prefactor=(1/(32*(1+reducedMass)*(Pi*fpi)**4))*((2*Pi)**3/HC)
      prefactor=(1/(32*(1+reducedMass)*(Pi*fpi)**4))
      
      factorAsym=-4*mpi0*mpi0*prefactor*(1/(DOT_PRODUCT(qVec,qVec)))*((-1)*(-1)**(t12))
      factorBsym=-1*gA*gA*prefactor/((DOT_PRODUCT(qVec,qVec)+mpi0**2))*((2*t12*(t12+1))-3)
      factorCsym=gA*gA*prefactor*(1/((DOT_PRODUCT(qVec,qVec)+mpi0**2)**2))*((2*t12*(t12+1))-3)
      ! factorDsym=-2*gA*gA*prefactor*mPion*mPion*(1/((DOT_PRODUCT(qVec,qVec)+mPion**2)**2))


      factorAasy=factorAsym
      factorBasy=factorBsym
      factorCasy=factorCsym
      ! factorDasy=factorDsym
      ! write(*,*) "kVec=", kVec 
      ! write(*,*) "kpVec=", kpVec 
      ! write(*,*) "qVec=", qVec 
      ! write(*,*) "prefactor=", prefactor 
      ! write(*,*) "mpi0=", mpi0 
      ! write(*,*) "factorAsym=", factorAsym 
      ! write(*,*) "factorBsym=", factorBsym 
      ! write(*,*) "factorCsym=", factorCsym 
      if ((t12 .eq. t12p) .and. (mt12 .eq. mt12p)) then

         if (s12p .eq. s12) then ! spin symmetric part only; s12-s12p=0 => l12-l12p is even
            call CalcKernel2BAsym(Kerneltmp,
     &           factorAsym,
     &           s12p,s12,t12,mt12,extQnumlimit,verbosity)

            if (ANY(Kerneltmp.Ne.Kerneltmp)) then
              write(*,*) "Sym A Kernel contains NaN values"
              write(*,*) "factorAsym=", factorAsym 
              stop
            end if

            call CalcKernel2BBsym(Kerneltmp,qVec,
     &           factorBsym,
     &           s12p,s12,t12,extQnumlimit,verbosity)

            if (ANY(Kerneltmp.Ne.Kerneltmp)) then
              write(*,*) "Sym B Kernel contains NaN values"
              stop
            end if
            call CalcKernel2BCsym(Kerneltmp,qVec,
     &           factorCsym,
     &           s12p,s12,t12,extQnumlimit,verbosity)

            if (ANY(Kerneltmp.Ne.Kerneltmp)) then
              write(*,*) "Sym c Kernel contains NaN values"
              stop
            end if
c           call CalcKernel2BDsym(Kernel2B,qVec,
c    &           factorDsym,
c    &           s12p,s12,mt12,extQnumlimit,verbosity)

         else                   !  spin anti-symmetric part only; s12 question: s12-s12p=±1 => l12-l12p is odd

            call CalcKernel2BAasy(Kerneltmp,
     &           factorAsym,
     &           s12p,s12,t12,extQnumlimit,verbosity)

            if (ANY(Kerneltmp.Ne.Kerneltmp)) then
              write(*,*) "Asym A Kernel contains NaN values"
              stop
            end if
            call CalcKernel2BBasy(Kerneltmp,qVec,
     &           factorBsym,
     &           s12p,s12,t12,extQnumlimit,verbosity)

            if (ANY(Kerneltmp.Ne.Kerneltmp)) then
              write(*,*) "Asym B Kernel contains NaN values"
              stop
            end if
            call CalcKernel2BCasy(Kerneltmp,qVec,
     &           factorCsym,
     &           s12p,s12,t12,extQnumlimit,verbosity)

            if (ANY(Kerneltmp.Ne.Kerneltmp)) then
              write(*,*) "Asym C Kernel contains NaN values"
              stop
            end if
c           call CalcKernel2BDasy(Kernel2B,qVec,
c    &           factorDsym,
c    &           s12p,s12,mt12,extQnumlimit,verbosity)
         end if                 ! s12 question
      else                      ! t12!=t12p
        continue
      end if                    !t12 question

      ppVecA=ppVec
      Kerneltmp=Kerneltmp*Jacobian
      end subroutine getDiagAB



      function vecsquare(vec)
          implicit none
          real*8 vec(3)
          real*8 vecsquare

          vecsquare=DOT_PRODUCT(vec,vec)
      end function vecsquare
