c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018:
c     subroutine constructs the 1N Compton amplitudes for proton and neutron
c     modified from onebody/constructAmps.f
c
c     twoSmax/twoMz dependence: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c      
c     hgrie May 2022: set nucleon mass to Mnucleon in Thomson term
      
c     hgrie May 2018: more extensive description of changes in main.*.f
c                     rewritten such that magnetic quantum numbers are "2xMz" etc
c      
c     hgrie May 2018: output now provides Ai's for
c     proton (mt1N=+1/2) and neutron (mt1N=-1/2)
c     for convolution with 1Ndensity quantum numbers.  
c     Old output was for isoscalar and isovector. 
c      
c     hgrie June 2014:
c     modified to include piN and Delta routines of
c     deuteron code -- as much as possble followed structure of
c     deuteron file constructAmps.f
c     deleted all reference to potential shift of polarisabilities
c     split original file 1Bspinisospintrans.f
c     This file only contains 1N calculation
c     
c     BS: Added variable 'variedA' to the routine. If variedA != 0 the
c     full calculation is not performed, but rather the simple variations
c     for  amplitude indicated by variedA (1..6) are calculated, other
c     amplitudes forced to 0

c     hgrie 19 Oct 2014: modified varyA => varyAp and varyAn:
c     hgrie 19 Oct 2014:
c     if calctype = VaryAn, proton amplitudes A1p-A6p set to zero
c     and only one of the neutron amps (number variedA)
c     is nonzero.
c     vice versa for VaryAp.
c     hgrie Nov 2014: implemented corrected boost routine for A2 in
c     calculateAsOQ3poles:
c     depends on mass of target nucleus in units of nucleon mass 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine Construct1NAmps(A1p,A2p,A3p,A4p,A5p,A6p,
     &     A1n,A2n,A3n,A4n,A5n,A6n,
     &     thetacm,xq,wx,Nx,t,omega,calctype,print1Namps,
     &     variedA,verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      implicit none
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
      real*8 A1p,A2p,A3p,A4p,A5p,A6p! sum of all contribs, proton
      real*8 A1n,A2n,A3n,A4n,A5n,A6n! sum of all contribs, neutron
 
      real*8 A1pOQ3,A2pOQ3,A3pOQ3,A4pOQ3,A5pOQ3,A6pOQ3! non-poles contribs at Odelta2, proton
      real*8 A1nOQ3,A2nOQ3,A3nOQ3,A4nOQ3,A5nOQ3,A6nOQ3! non-poles contribs at Odelta2, neutron
      real*8 A1pOQ3p,A2pOQ3p,A3pOQ3p,A4pOQ3p,A5pOQ3p,A6pOQ3p! poles contribs at Odelta2, proton
      real*8 A1nOQ3p,A2nOQ3p,A3nOQ3p,A4nOQ3p,A5nOQ3p,A6nOQ3p! poles contribs at Odelta2, neutron
      real*8 ADelta1,ADelta2,ADelta3,ADelta4,ADelta5,ADelta6! Delta contribs at Odelta3; isoscalar

      integer verbosity
      
      logical print1Namps
      
      integer variedA

      real*8 thetacm
c     
      integer Nx,calctype
      real*8 xq(Nxmax),wx(Nxmax)
      real*8 omega,t

      character (len=30) units ! output in MeV^-1 for A1-6, but MeV^0 for VaryA
      
c     for mathematica-friendly output, define numbers as strings.
c     not elegant, but works
      character(len=29) string1,string2,string3,string4,string5,string6
c----------------------------------------------------------------------
c     define units for screen output
      if (( calctype.eq.VaryAp ).or.( calctype.eq.VaryAp )) then
         units = "MeV⁰: insertion of 1/αEM"
      else
         units = "MeV⁻¹"
      end if
c      
c----------------------------------------------------------------------
c     
c     Calculate the photon amplitudes for the particular gamma-N energy
c     
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c     BS: Add option for calctype = VaryA. Then none of the calc's
c     above are performed, instead appropriate CalculateDA is called
c     hgrie 19 Oct 2014: extended to allow variation of n or p amps, see above.
c----------------------------------------------------------------------
      if ( calctype.eq.VaryAn ) then

         A1p=0.d0
         A2p=0.d0
         A3p=0.d0
         A4p=0.d0
         A5p=0.d0
         A6p=0.d0

         if (variedA .eq. 1 ) then
            call CalculateDA1(A1n,A2n,A3n,A4n,A5n,A6n,verbosity)
         else if (variedA .eq. 2) then
            call CalculateDA2(A1n,A2n,A3n,A4n,A5n,A6n,verbosity)
         else if (variedA .eq. 3) then
            call CalculateDA3(A1n,A2n,A3n,A4n,A5n,A6n,verbosity)
         else if (variedA .eq. 4) then
            call CalculateDA4(A1n,A2n,A3n,A4n,A5n,A6n,verbosity)
         else if (variedA .eq. 5) then
            call CalculateDA5(A1n,A2n,A3n,A4n,A5n,A6n,verbosity)
         else if (variedA .eq. 6) then
            call CalculateDA6(A1n,A2n,A3n,A4n,A5n,A6n,verbosity)
         endif
         
      else if ( calctype.eq.VaryAp ) then

         A1n=0.d0
         A2n=0.d0
         A3n=0.d0
         A4n=0.d0
         A5n=0.d0
         A6n=0.d0

         if (variedA .eq. 1 ) then
            call CalculateDA1(A1p,A2p,A3p,A4p,A5p,A6p,verbosity)
         else if (variedA .eq. 2) then
            call CalculateDA2(A1p,A2p,A3p,A4p,A5p,A6p,verbosity)
         else if (variedA .eq. 3) then
            call CalculateDA3(A1p,A2p,A3p,A4p,A5p,A6p,verbosity)
         else if (variedA .eq. 4) then
            call CalculateDA4(A1p,A2p,A3p,A4p,A5p,A6p,verbosity)
         else if (variedA .eq. 5) then
            call CalculateDA5(A1p,A2p,A3p,A4p,A5p,A6p,verbosity)
         else if (variedA .eq. 6) then
            call CalculateDA6(A1p,A2p,A3p,A4p,A5p,A6p,verbosity)
         endif

      else ! If it's NOT a VaryA calculation, it MUST be at least OQ2:
c     define OQ2 polarisabilities 
         A1p = -ep**2/Mnucleon 
         A2p = 0.d0
         A2p = 0.d0
         A3p = 0.d0
         A4p = 0.d0
         A5p = 0.d0
         A6p = 0.d0
      
         A1n = 0.d0 
         A2n = 0.d0
         A2n = 0.d0
         A3n = 0.d0
         A4n = 0.d0
         A5n = 0.d0
         A6n = 0.d0
      endif
c------------------------------------------------------------------
c     for all that is not OQ2 calculation: include OQ3 Npi diagrams
c     hgrie 19 Oct 2014: VaryA => VaryAn or VaryAp as described above.
c------------------------------------------------------------------------------
      
      if ( calctype.gt.OQ2 ) then ! setting ".gt." includes ALL Oe²δ^(3,4,..), but not VaryA or OQ2=Oe²δ⁰
c     "no-poles" contributions (i.e. piN diagrams and pi-pole)
c     hardwire that pols are NOT shifted: 6 zeroes at nex-to last position
c     last: verbosity irrelevant for subroutine: arb. value
         call CalculateAsOQ3nopoles(A1pOQ3,A2pOQ3,A3pOQ3,A4pOQ3,A5pOQ3,A6pOQ3,xq,wx,Nx,
     &        t,omega,thetacm,kappap,ep,0,0,0,0,0,0,verbosity)
         call CalculateAsOQ3nopoles(A1nOQ3,A2nOQ3,A3nOQ3,A4nOQ3,A5nOQ3,A6nOQ3,xq,wx,Nx,
     &        t,omega,thetacm,kappan,en,0,0,0,0,0,0,verbosity)
c     "poles" contributions (i.e. Nucleon-pole diagrams)
c     last: verbosity irrelevant for subroutine: arb. value
         call CalculateAsOQ3poles(A1pOQ3p,A2pOQ3p,A3pOQ3p,A4pOQ3p,A5pOQ3p,A6pOQ3p,
     &        omega,thetacm,kappap,ep,3.d0 ! target nucleus mass in units of nucleon mass: for A2 boost correction
     &        ,verbosity)
         call CalculateAsOQ3poles(A1nOQ3p,A2nOQ3p,A3nOQ3p,A4nOQ3p,A5nOQ3p,A6nOQ3p,
     &        omega,thetacm,kappan,en,3.d0 ! target nucleus mass in units of nucleon mass: for A2 boost correction
     &        ,verbosity)
         
c     Now add p and n contributions (they differ by pi0-piece only)
c     Note that calculateAsOQ3XXX gives amplitudes without omega^2 factors,
c     while 3He code uses definitions including omega^2 -- so adjust here.
         A1p=A1p+A1pOQ3+A1pOQ3p 
         A2p=A2p+(A2pOQ3+A2pOQ3p)*omega**2
         A3p=A3p+A3pOQ3+A3pOQ3p
         A4p=A4p+(A4pOQ3+A4pOQ3p)*omega**2
         A5p=A5p+(A5pOQ3+A5pOQ3p)*omega**2
         A6p=A6p+(A6pOQ3+A6pOQ3p)*omega**2
         
         A1n=A1n+A1nOQ3+A1nOQ3p
         A2n=A2n+(A2nOQ3+A2nOQ3p)*omega**2
         A3n=A3n+A3nOQ3+A3nOQ3p 
         A4n=A4n+(A4nOQ3+A4nOQ3p)*omega**2
         A5n=A5n+(A5nOQ3+A5nOQ3p)*omega**2
         A6n=A6n+(A6nOQ3+A6nOQ3p)*omega**2
      end if
c     end of calculation of amplitudes at Q^3  
c-------------------------------------------------------------------------------
c     add Delta from O(epsilon^3) when appropriate --
c     all structure, no N-pole contributions!

      if (calctype.eq.Oepsilon3) then
         call CalculateAsDelta(ADelta1,ADelta2,ADelta3,ADelta4,ADelta5,ADelta6,
     &        wx,xq,Nx,t,omega,verbosity)
         
c     Now add p and n contributions -- proton and neutron same
c     Note that calculateADelta gives amplitudes without omega^2 factors,
c     while 3He code uses definitions including omega^2 -- so adjust here.
         A1p=A1p+ADelta1 
         A2p=A2p+ADelta2*omega**2
         A3p=A3p+ADelta3
         A4p=A4p+ADelta4*omega**2
         A5p=A5p+ADelta5*omega**2
         A6p=A6p+ADelta6*omega**2
         
         A1n=A1n+ADelta1 
         A2n=A2n+ADelta2*omega**2
         A3n=A3n+ADelta3
         A4n=A4n+ADelta4*omega**2
         A5n=A5n+ADelta5*omega**2
         A6n=A6n+ADelta6*omega**2

      endif         
c     end of calculation of amplitudes with Delta 
c     
c     end hgrie modification 2014
c------------------------------------------------------------------------------
c     option: 
c     output proton and neutron single-N amplitudes used (before boost)
c
      if (print1Namps) then
         write(*,*) "   Single-nucleon amplitude used for convolution:"
         write(*,*) "   [pre-boost, basis as in review (2.1), no factor e^2]"
         write(*,*) "   Proton:"
         write(*,*) "       A1p = ",A1p," ",units
         write(*,*) "       A2p = ",A2p," ",units
         write(*,*) "       A3p = ",A3p," ",units
         write(*,*) "       A4p = ",A4p," ",units
         write(*,*) "       A5p = ",A5p," ",units
         write(*,*) "       A6p = ",A6p," ",units
         write(*,*) "   Neutron:"
         write(*,*) "       A1n = ",A1n," ",units
         write(*,*) "       A2n = ",A2n," ",units
         write(*,*) "       A3n = ",A3n," ",units
         write(*,*) "       A4n = ",A4n," ",units
         write(*,*) "       A5n = ",A5n," ",units
         write(*,*) "       A6n = ",A6n," ",units
c     
         write(*,*) "   Mathematica-friendly: {A1,A2,A3,A4,A5,A6} [MeV^-1]"
         write(string1,*) A1p
         write(string2,*) A2p
         write(string3,*) A3p
         write(string4,*) A4p
         write(string5,*) A5p
         write(string6,*) A6p
         write(*,*) "   Proton:"
         write(*,*) '   {',string1(1:21),'*10^(0',string1(23:26),') , ',
     &        string2(1:21),'*10^(0',string2(23:26),') , ',
     &        string3(1:21),'*10^(0',string3(23:26),') , ',
     &        string4(1:21),'*10^(0',string4(23:26),') , ',
     &        string5(1:21),'*10^(0',string5(23:26),') , ',
     &        string6(1:21),'*10^(0',string6(23:26),') }'
         write(string1,*) A1n
         write(string2,*) A2n
         write(string3,*) A3n
         write(string4,*) A4n
         write(string5,*) A5n
         write(string6,*) A6n
         write(*,*) "   Neutron:"
         write(*,*) '   {',string1(1:21),'*10^(0',string1(23:26),') , ',
     &        string2(1:21),'*10^(0',string2(23:26),') , ',
     &        string3(1:21),'*10^(0',string3(23:26),') , ',
     &        string4(1:21),'*10^(0',string4(23:26),') , ',
     &        string5(1:21),'*10^(0',string5(23:26),') , ',
     &        string6(1:21),'*10^(0',string6(23:26),') }'
      endif   

      return
      end
