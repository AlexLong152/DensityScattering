c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018, based on common/readinput.f
c
c     subroutines to read input files for density calculations:
c
c     ReadinputCommon: input for parameters which are identical for onebody and twobody
c     ReadinputOnebody: input for parameters for onebody
c     ReadinputTwobody: input for parameters for twobody
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     hgrie Sep 2020: in readinputTwobody():
c        -- Removed NP12p = 3rd variable in integration-grid
c           line of input file. It is actually never used for anything!
c           That also reduced number of arguments of readinputTwobody()!
c        -- Changed read of angle parameters such that input file does not need to
c           contain 5 numbers when only 2 are needed for LL
c     hgrie Aug 2020: added readinputCommonComments()
c     hgrie Aug 2020: corected couting of independent MEs:
c      now ...+2*mod(twoSnucl+1,2) for correct number -- was 2*(mod(twoSnucl,2)-1)
c      
c     hgrie June 2018: renamed "parity" to "symmetry
c     -- see notes in usesymmetry+writeoutput.densities.f
c      
c     hgrie May 2018: used to be part of 3HeCompt/common
c     now part of common-densities, backward compatibility deliberately broken
c     no changes yet
c      
c     hgrie May 2018: more extensive description of changes in main.*.f
c                     rewritten such that magnetic quantum numbers are "2xMz" etc
c      
c     -------------------------------------------------------------------
c     LEGACY: Relevant notes on original common/redinput.f:
c     this combined the 3 subroutines (with different ordering in input.dat) into one
c      
c     Written by D.P.-11/97
c     
c     Bruno Strandberg
c     ********************************************************************
c     Rev: modified 4 June 2014 from Deepshikha's 3He code
c     1. Got rid of densitytype, deuttype
c     2. Got rid of palpha,pbeta,nalpha,nbeta
c     3. Got rid of dg1p,dg2p,dg3p,dg4p,dg1n,dg2n,dg3n,dg4n
c     4. Got rid of firsttime - remove from other parts of code!
c     5. Got rid of ampfile,amp1Bfile,amp2Bfile
c     6. Got rid of ampUnitno,amp1BUnitno,amp2BUnitno
c     7. Created thetaLow, thetaHigh
c     8. Changed interval-->thetaInterval
c     9. Added Oepsilon3 option to calcstring parsing
c     10. Added whichbody option
c     11. Added variables to control quadrature settings
c     ********************************************************************
c     hgrie 20 June 2014: modified for use of LebedevLaikov or Gaussian
c     integration for theta & phi separately,
c     for solid angle integral in (12) system
c     
c     This routine reads in the values of the external photon momentum &
c     scattering angle in the lab frame. It takes them from the I/O
c     unit denoted by inUnitno, where they are stored in MeV and degrees.
c     
c     *******************************************************************
c     hgrie Feb 2016: added OQ4 and Odelta4 options,
c     but activated only for twobody!     
c     ********************************************************************
c
c     hgrie May 2017: implemented that j12max (max total ang mom in (12) subsystem)
c                     can be set in input file.
c                     Defaults are j12max=2 for onebody and j12max=1 for twobody.
c                     That is enough for convergence on the <1% level in amplitudes.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018
c     Now the routine which reads input specific for onebody

      subroutine ReadinputOnebody(inUnitno,calctype,variedA,descriptors,
c---- Variable to control Feynman quadrature settings------------------------
     &     Nx,
     &     verbosity)
c     
      implicit none
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
c     
c     
c     VARIABLES PASSED INOUT:
      character*200,intent(inout) :: descriptors ! additional descriptors of calculation for outputfilename
      
c     VARIABLES PASSED OUT:
c     
      integer,intent(out) :: calctype
      integer,intent(out) :: variedA
      integer,intent(out) :: Nx                ! number of quadratures for loop integration in one-body amp.
      
c     ---------------------------------------------------------------
c     BS: Ultimately variable def's in Readme?
c     ---------------------------------------------------------------
c     
c     calctype    - which calculation to do
c     hgrie 19 Oct 2014: modified to VaryAp or VaryAn
c     to indicate if p or n amp varied
c     hgrie 19 Oct 2014: changed IA assignation since unused. VaryA was 5;
c     now VaryAp is 0, VaryAn is 1. Leaves numbers >=5 open for future.
c     -1: uninitialised value.
c     0=VaryAp, vary a proton amplitude defined by variedA, all other amps 0
c     1=VaryAp, vary a neutron amplitude defined by variedA, all other amps 0
c     2=OQ2=Odelta0, the O(q²) calculation Chi PT calculation;
c     3=OQ3=Odelta2, the full O(q³) Chi PT calculation;
c     4=Oepsilon3=Odelta3, calculation with delta diagrams.
c     
c     BS:----
c     variedA to indicate which A is varied by calctype=VaryA
c     *******************************************************************
c     
c     VARIABLES PASSED IN:
c     
      integer,intent(in) :: inUnitno ! I/O unit containing all the input information
      integer,intent(in) :: verbosity
c     
c     *******************************************************************
c     
c     INTERNAL VARIABLES:
c     
      character*500 calcstring ! used to hold info. on type of calculation being done
c     
c     *******************************************************************
c     
      if (verbosity.eq.1000) continue ! keep for future use
      if (descriptors.ne.'a') continue ! keep for future use
      
c     Read in the input file line-by-line, write output info to terminal
c     
      read (inUnitno,'(A80)') calcstring
c     Read in Feynman quadrature settings------------------------------------
      read (inUnitno,*) Nx
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
c     Determine the order of the calculation to be done

      if ((index(calcstring,'OQ2').ne. 0).or.(index(calcstring,'Odelta0').ne. 0)) then 
         calctype=OQ2
         write (*,*) 'O(e²delta⁰)=O(Q²) calculation.'
      else if ((index(calcstring,'OQ3').ne. 0).or.(index(calcstring,'Odelta2').ne. 0) ) then
         calctype=OQ3
         write (*,*) 'O(e²delta²)=O(Q³) calculation.'
      else if ((index(calcstring,'Oepsilon3').ne. 0).or.(index(calcstring,'Odelta3').ne. 0) ) then
         calctype=Oepsilon3
         write (*,*) 'O(e²delta³)=O(epsilon³) calculation.'
         write (*,*) '  Delta and DeltaPi parts of static scalar polarisabilities are subtracted.'
         write (*,*) "  Delta parameters:"
         write (*,*) "    Delta mass                  Mdelta    = ",Mdelta," MeV"
         write (*,*) "      => Delta-N mass splitting Delta     = ",delta," MeV"
         write (*,*) "    PionNDelta coupling         gPiNDelta = ",gpind
         write (*,*) "    DeltaGamma coupl. (nonrel.) b1        = ",b1
      else if ((index(calcstring,'OQ4').ne. 0).or.(index(calcstring,'Odelta4').ne. 0) ) then
         
c     hgrie note Feb 2017: when implemented, need to use different LECs c_i, for onebody δ⁴ vs Q⁴

         if (index(calcstring,'OQ4').ne. 0) then
            calctype=OQ4
            write (*,*) 'O(Q⁴) calculation.'
            write (*,*) '   -- NOT YET IMLPEMENTED FOR ONEBODY. -- Exiting.'
            stop
         else
            calctype=Odelta4
            write (*,*) 'O(e²delta⁴) calculation.'
            write (*,*) '   -- NOT YET IMLPEMENTED FOR ONEBODY. -- Exiting.'
            stop
         end if
c     variation of proton amplitudes         
      else if (index(calcstring,'VaryA1p').ne. 0) then
         calctype=VaryAp
         variedA=1
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA2p').ne. 0) then
         calctype=VaryAp
         variedA=2
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA3p').ne. 0) then
         calctype=VaryAp
         variedA=3
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA4p').ne. 0) then
         calctype=VaryAp
         variedA=4
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA5p').ne. 0) then
         calctype=VaryAp
         variedA=5
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA6p').ne. 0) then
         calctype=VaryAp
         variedA=6
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
c     variation of neutron amplitudes         
      else if (index(calcstring,'VaryA1n').ne. 0) then
         calctype=VaryAn
         variedA=1
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA2n').ne. 0) then
         calctype=VaryAn
         variedA=2
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA3n').ne. 0) then
         calctype=VaryAn
         variedA=3
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA4n').ne. 0) then
         calctype=VaryAn
         variedA=4
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA5n').ne. 0) then
         calctype=VaryAn
         variedA=5
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA6n').ne. 0) then
         calctype=VaryAn
         variedA=6
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else
c     hgrie 19 Oct 2014: if calctype not yet specified, things went wrong. 
         write (*,*) '*** ERROR: Calculation type unknown. -- Exiting.'
         stop
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     write quadrature variables to terminal
c     
      write (*,*) "No integrals in onebody."
c     onebody only: number of points of Feynman parameter integration
      write (*,180) Nx

c     end hgrie mod         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     formats
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 180  format (1X,"Number of quadratures in Feynman parameter integral of single-N amplitude = ",I2)
      
      return
      end
      
