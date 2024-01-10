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
      write(*,*) "Kernel: Twobody Pion Pion Scattering"
      write(*,*) "--------------------------------------------------------------------------------"
      write(*,*) "   Kernel Code Version 1.0"
      write(*,*) "      Alexander Long/hgrie starting December 2023"
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
c     write(*,*) "Output: F_{L/T} ε.S^Mprime_M = < Mprime | kernel | M > Lenkewitz Eur. Phys. J. A (2013) 49:20 eq (10)"
c     write(*,*) "        with ε: incoming-photon polarisation, S: nucleus spin; [F_{L/T}]=[fm]¯¹"
c     write(*,*) "        Mapping of extQnum: 1 = ε_x, 2 = ε_y (both transversal); 3 = ε_z (longitudinal)"
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
      subroutine Calc2Bspinisospintrans(Kernel2B,Mnucl,
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
      real*8 identity(3,3), tmp(3)
      integer i
      real*8 isospinFactor(0:1,-1:1,0:1,-1:1)
      real*8 tmpVec(3), tmpVec2(3)
      real*8 pVec(3), ppVec(3)
      real*8 qVec(3)
      real*8 dl12by2
      real*8 factorAsym,factorBsym
      real*8 factorAasy,factorBasy
      real*8 kVec(3),q1Vec(3)
      real*8 kpVec(3)
      real*8 mPion, Mnucl,Epi,q0!, Fpi
      real*8 isospin(3)
      real*8 mu
      real*8 m1,m2,m3,m4
      real*8 q(3), q1(3)
      real*8 kp(3), kp2Mass
      real*8 kpPion(3)
      real*8 tmpk
c     
c     TODO: Fix this jank
      tmpk=k
      if (tmpk.le.131.89) then
          tmpk=131.89
      end if
      mPion=134.976
c     call calculateqs2Mass(pVec,ppVec,qVec,k,kVec,kpVec,thetacm,mPion,mNucl,verbosity)
      m1=mNucl
      m2=0.d0!photon
      m3=mNucl
      m4=mPion
      kVec=(/0.d0,0.d0,real(tmpk,8)/)! TODO: Fix tmpk here when fixed the jank above

      pVec=(/px,py,pz/)
      ppVec=(/ppx,ppy,ppz/)
      kp2Mass=0.d0
      kpPion=0.d0
      call calculateqs2Mass(pVec,ppVec,kVec,kp2Mass,m1,m2,m3,m4,thetacm,verbosity)
      call calculateqsmass(pVec,ppVec,q,tmpk,q1,kVec,kpPion,thetacm,mPion,mNucl,verbosity)!pion photoproduction specific
      write(*,*) ""
      write(*,*) ""
      if (verbosity.eq.1000) continue
      end
