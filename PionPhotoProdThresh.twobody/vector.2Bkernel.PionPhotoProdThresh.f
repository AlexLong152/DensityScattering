cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of KERNEL code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              KernelGreeting         : output to stdout with message which process computed and its version number
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
      subroutine Calc2Bspinisospintrans(Kernel2B,Mnucl,
c     &     Comp2Bx,Comp2By,Comp2Bpx,Comp2Bpy, ! for STUMP, see below
     &     extQnumlimit,
     &     t12,mt12,t12p,mt12p,l12,s12,
     &     l12p,s12p,thetacm,k,px,py,pz,ppx,ppy,ppz,calctype,verbosity)
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
      
      complex*16,intent(out) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1) ! was Comp2Bxx/xy/yx/yy
c     
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
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: t12,mt12,t12p,mt12p,l12,l12p,s12,s12p
      real*8,intent(in)  :: px,py,pz,ppx,ppy,ppz
               
      integer,intent(in) :: verbosity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
      real*8 tmpVec(3), tmpVec2(3)
      real*8 pVec(3), ppVec(3)
      real*8 qpx,qpy,qpz, qpVec(3)
      real*8 qppx,qppy,qppz, qppVec(3)
      real*8 qx,qy,qz, qVec(3)
      real*8 q12x,q12y,q12z, q12Vec(3)
      real*8 qp12x,qp12y,qp12z, qp12(3)
      real*8 qpp12x,qpp12y,qpp12z, qpp12(3)
      real*8 qsq,qpsq,qppsq,q12sq,qp12sq,qpp12sq
      real*8 qpppx,qpppy,qpppz, qppp(3)
      real*8 qppp12x,qppp12y,qppp12z, qppp12(3)
      real*8 qpppsq,qppp12sq
      real*8 dl12by2
      real*8 factorA,factorB
      real*8 factorAvec,factorBvec
      real*8 factorAasy,factorBasy
      real*8 kVec(3),q1Vec(3)
      real*8 mPion, Mnucl
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
      Kernel2B=c0
      dl12by2=(l12-l12p)/2.d0   !to check if l12-l12p is  even or odd
c     
c     Calculate momenta q,q',q':
c     

       pVec=(/px,py,pz/)
       ppVec=(/ppx,ppy,ppz/)
       kVec=(/0.d0,0.d0,real(k,8)/)
       mPion=134.976
       call calculateqsmass(pVec,ppVec,qVec,k,q1Vec,kVec,thetacm,mPion,Mnucl,verbosity)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Odelta0 2N contributions: NONE
c     <if they were nonzero, enter diagrams here>
      if (calctype.eq.Odelta0) return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Odelta2 2N contributions
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      tmpVec=pVec-ppVec+(kVec/2)
      tmpVec2=pVec-ppVec-(kVec/2)

      factorA=-(-1)**(t12)*(1.d0/(DOT_PRODUCT(tmpVec,tmpVec)))*(2*Pi)**3/HC
      factorB=+2*(-1)**(t12)*(1.d0/
     &            DOT_PRODUCT(tmpVec,tmpVec))*
     &            (1.d0/(DOT_PRODUCT(tmpVec2,tmpVec2)+mpi2))
     &         *(2*Pi)**3/HC

c     antisymmetric part: turns out to be the same, only the vaue of t12 will be different
      factorAasy=factorA
      factorBasy=factorB
      
      if ((t12 .eq. t12p) .and. (mt12 .eq. 0) .and.(mt12p .eq. 0)) then
         if (s12p .eq. s12) then ! s12-s12p=0 => l12-l12p is even; spin symmetric part only

            call CalcKernel2BAsym(Kernel2B,
     &           factorA,
     &           s12p,s12,extQnumlimit,verbosity)
            call CalcKernel2BBsymVec(Kernel2B,
     &           factorB,
     &           pVec-ppVec-(kVec/2), ! preceding is vector dotted with σ
     &           pVec-ppVec, ! preceding is vector dotted with ε
     &           s12p,s12,extQnumlimit,verbosity)
         else                   ! s12 question: s12-s12p=±1 => l12-l12p is odd; spin anti-symmetric part only
c     
            call CalcKernel2BAasy(Kernel2B,
     &           factorAasy,
     &           s12p,s12,extQnumlimit,verbosity)
            call CalcKernel2BBasyVec(Kernel2B,
     &           factorBasy,
     &           pVec-ppVec-(kVec/2), ! preceding is vector dotted with σ
     &           pVec-ppVec, ! preceding is vector dotted with ε
     &           s12p,s12,extQnumlimit,verbosity)


         end if                 ! s12 question
      else                      ! t12!=t12p
         continue
c     diagrams (A/B) have no components with t12!=t12p. 
      end if                    !t12 question
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
