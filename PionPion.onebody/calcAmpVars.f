c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c***************************************************************
c     Bruno Strandberg: copied from Deepshikha's CalculateDeltaAs.f
c     
c     These subroutines calculate simple variations to amplitudes
c     A1 from A6. These variations are added to the full calculated
c     amplitudes in a different subroutine to estimate how
c     polarisability shifts affect cross-sections
c     hgrie 19 Oct 2014: original file took A2 to be A2 = omega/Mnucl**2
c     This is a residual of some test runs. Corrected to read
c     A2 = + omega**2*(1.d0/(alphaem*HC**3.d0)), identical to A1 variation.
c     [It may be that D. used same with "-" sign, see commented bit!]
c     hgrie 27 Oct 2014: CHANGED variation to "unity",
c     i.e. changed the one varied (nonzero) Ai:
c     Ai = ω^(2/3)/(αem hbarc^(3,4)) (OLD) ==> Ai = 1/αem (NEW)
c     While old variation led directly to polarisabilities, the NEW ONE
c     will alleviate future variations, cross checks, and adding boosts.
c     NOTE: Amplitudes still divided by αem=1/137 because OVERALL amplitudes
c     are normalised that way.
c     => In order to calculate e.g. magnetic proton dipole polarisability
c     change by 2 canonical units,
c     calculate -4 Pi omega**2 [betaM1=2]*10**(-4)/hbarc**3*VaryA2p
c     etc; see [GMPF review (2.11)].
c     
c     ***************************************************************
c     
      subroutine CalculateDA1(A1,A2,A3,A4,A5,A6,verbosity)
c     
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
c     
c**********************************************************************
c     
c     Output variables:
c     
c     A1,A2,A3,A4,A5,A6-six invariant functions for gamma-N
c     A1 gets a value corresponding to unit variation, rest forced to 0
c     
      real*8 A1,A2,A3,A4,A5,A6
      integer verbosity
c     
c**********************************************************************
      A1 = 1.d0/alphaem         ! was: omega**2*(1.d0/(alphaem*HC**3.d0))
      A2 = 0.d0
      A3 = 0.d0
      A4 = 0.d0
      A5 = 0.d0
      A6 = 0.d0
      if (verbosity.eq.1000) continue
      return
      end
      
c     
c     ***************************************************************
c     
      subroutine CalculateDA2(A1,A2,A3,A4,A5,A6,verbosity)
c     
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
c     
c**********************************************************************
c     
c     Output variables:
c     
c     A1,A2,A3,A4,A5,A6-six invariant functions for gamma-N
c     A2 gets a value corresponding to unit variation, rest forced to 0
c     
      real*8 A1,A2,A3,A4,A5,A6
      integer verbosity
c     
c**********************************************************************
c     
      A1 = 0.d0
      A2 = 1.d0/alphaem         ! was: omega**2*(1.d0/(alphaem*HC**3.d0))!- omega**2*(1.d0/(alphaem*HC**3.d0))
      A3 = 0.d0
      A4 = 0.d0
      A5 = 0.d0
      A6 = 0.d0
      if (verbosity.eq.1000) continue
      return
      end
      
c     
c     ***************************************************************
c     
      subroutine CalculateDA3(A1,A2,A3,A4,A5,A6,verbosity)
c     
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
c     
c**********************************************************************
c     
c     Output variables:
c     
c     A1,A2,A3,A4,A5,A6-six invariant functions for gamma-N
c     A3 gets a value corresponding to unit variation, rest are forced to 0
c     
      real*8 A1,A2,A3,A4,A5,A6
      integer verbosity
c     
c**********************************************************************
      A1 = 0.d0
      A2 = 0.d0
      A3 = 1.d0/alphaem         ! was: omega**3*(1.d0/(alphaem*HC**4.d0))
      A4 = 0.d0
      A5 = 0.d0
      A6 = 0.d0
      if (verbosity.eq.1000) continue
      return
      end
      
c
c ***************************************************************
c     
      subroutine CalculateDA4(A1,A2,A3,A4,A5,A6,verbosity)
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
c     
c**********************************************************************
c     
c     Output variables:
c     
c     A1,A2,A3,A4,A5,A6-six invariant functions for gamma-N
c     A4 gets a value corresponding to unit variation, rest are forced to 0
c     
      real*8 A1,A2,A3,A4,A5,A6
      integer verbosity
c     
c**********************************************************************
      A1 = 0.d0
      A2 = 0.d0
      A3 = 0.d0
      A4 = 1.d0/alphaem         ! was: omega**3*(1.d0/(alphaem*HC**4.d0))
      A5 = 0.d0
      A6 = 0.d0
      if (verbosity.eq.1000) continue
      return
      end
      
c     
c ***************************************************************
c
      subroutine CalculateDA5(A1,A2,A3,A4,A5,A6,verbosity)
c     
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
c     
c**********************************************************************
c     
c     Output variables:
c     
c     A1,A2,A3,A4,A5,A6-six invariant functions for gamma-N
c     A5 gets a value corresponding to unit variation, rest are forced to 0
c     
c     
      real*8 A1,A2,A3,A4,A5,A6
      integer verbosity
c     
c**********************************************************************
      A1 = 0.d0
      A2 = 0.d0
      A3 = 0.d0
      A4 = 0.d0
      A5 = 1.d0/alphaem         ! was: omega**3*(1.d0/(alphaem*HC**4.d0))
      A6 = 0.d0
      if (verbosity.eq.1000) continue
      return
      end
      
c     ***************************************************************
c     
      subroutine CalculateDA6(A1,A2,A3,A4,A5,A6,verbosity)
c     
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
c     
c**********************************************************************
c     
c     Output variables:
c     
c     A1,A2,A3,A4,A5,A6-six invariant functions for gamma-N
c     A6 gets a value corresponding to unit variation, rest are forced to 0
c     
      real*8 A1,A2,A3,A4,A5,A6
      integer verbosity
c     
c**********************************************************************
      A1 = 0.d0
      A2 = 0.d0
      A3 = 0.d0
      A4 = 0.d0
      A5 = 0.d0
      A6 = 1.d0/alphaem         ! was: omega**3*(1.d0/(alphaem*HC**4.d0))
      if (verbosity.eq.1000) continue
      return
      end
      
