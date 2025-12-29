cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of KERNEL code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              ReadinputTwobody : read parts of input specific to twobody kernel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, based on common/readinput-densities.f of Compton density code v2.0 hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     ReadinputTwobody: read input file for parameters for twobody
c      
c     j12max (max total ang mom in (12) subsystem) can be set in input file.
c     Defaults are j12max=2 for onebody and j12max=1 for twobody.
c     That is enough for convergence on the <1% level in amplitudes.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ReadinputTwobody(inUnitno,calctype,descriptors,
c---- Variables to control radial quadrature settings------------------------
     &     NP12A,NP12B,P12A,P12B,P12C,
c---- Variables to control angular quadrature settings------------------------
     &     AngularType12,Nanggrid12,
     &     Nordth12,Nordphi12,
     &     NthBins12,NphiBins12,
     &     j12max,verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read input file for parameters for twobody
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      implicit none

      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
c     
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INOUT/OUTPUT VARIABLES:
c      
      character*200,intent(inout) :: descriptors ! additional descriptors of calculation for outputfilename
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES:
c     
      integer,intent(out) :: calctype ! which order/type of calculation to do
      integer,intent(out) :: AngularType12,Nordth12,Nordphi12,Nanggrid12 ! define quadrature distribution for theta and phi integral
      integer,intent(out) :: NthBins12,NphiBins12
      integer,intent(out) :: j12max            !hgrie May 2017: maximum total ang mom in (12) system 
      
      integer,intent(out) :: NP12A,NP12B
      real*8,intent(out)  :: P12A,P12B,P12C
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
c     
      integer,intent(in) :: inUnitno ! I/O unit containing all the input information
      integer,intent(in) :: verbosity
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
c     
      character*500 calcstring ! holds info on type of calculation being done
      character*1  dummy ! a variable only for converting the integer j12max to a character string
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      if (verbosity.eq.1000) continue ! keep for future use
      
c     Read in the input file line-by-line, write output info to terminal
c     
      read (inUnitno,'(A80)') calcstring
      read (inUnitno,*) NP12A,NP12B
      read (inUnitno,*) P12A,P12B,P12C
c     Read in angular quadrature settings------------------------------------
c     hgrie Sep 2020: first determine angular type, then read only those parameters needed for that type
      read (inUnitno,*) AngularType12
      backspace inUnitno                ! reset to start of line
      if (AngularType12.eq.1) then      ! Gauss
         read (inUnitno,*) AngularType12,Nordth12,Nthbins12,Nordphi12,Nphibins12
      else if (AngularType12.eq.2) then ! LebedevLaikov
         read (inUnitno,*) AngularType12,Nordth12
      else                              ! none of the two: continue -- error message will be produced below
         continue
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     Determine the order of the calculation to be done
      
      if ((index(calcstring,'OQ2').ne. 0).or.(index(calcstring,'Odelta0').ne. 0)) then 
         write (*,*) 'Twobody has no O(e delta⁰)=O(Q²) MEC contribution.'
         write (*,*) '   -- Use at least O(e delta²)=O(Q³). -- Exiting.'
         stop
      else if ((index(calcstring,'OQ3').ne. 0).or.(index(calcstring,'Odelta2').ne. 0) ) then
         calctype=OQ3
         write (*,*) 'O(e delta²)=O(Q³) calculation. (12) subsystem is (pn) since only charged MECs.'
      else if ((index(calcstring,'Oepsilon3').ne. 0).or.(index(calcstring,'Odelta3').ne. 0) ) then
         write (*,*) 'Twobody has no O(e delta³)=O(epsilon³) MEC contribution.'
         write (*,*) '   -- Use O(e delta²)=O(Q³) or O(e delta⁴)=O(Q⁴) instead. -- Exiting.'
         stop
      else if ((index(calcstring,'OQ4').ne. 0).or.(index(calcstring,'Odelta4').ne. 0) ) then
c     hgrie note Feb 2017: once implemented, need to use different LECs c_i, for onebody δ⁴ vs Q⁴
         write (*,*) 'O(e delta⁴)=O(Q⁴) MECs identical [no Δ(1232)]. (12) subsystem is (pn) since only charged MECs.'
         calctype=Odelta4           ! again: MECs for Δ-less and Δ-ful are identical, so need no switch.
      else
c     hgrie 19 Oct 2014: if calctype not yet specified, things went wrong. 
         write (*,*) '*** ERROR: Calculation type unknown. -- Exiting.'
         stop
      end if
c
c     Now start on parameters of remaining integration
      write (*,*) "Integration parameters of (12) subsystem:"

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Determine maximum total angular momentum in (12) subsystem from input file
c     
      if (index(calcstring,'j12max=').eq.0) then
         j12max = 1             ! twobody converges faster
         write(*,'(A,I4)') " Total angular momentum in (12) subsystem not specified -- using default j12max =",j12max
      else if (index(calcstring,'j12max=').ne.0) then
         if (index(calcstring,'j12max=1').ne.0) then
            j12max = 1
         else  if (index(calcstring,'j12max=0').ne.0) then
            j12max = 0
         else  if (index(calcstring,'j12max=2').ne.0) then
            j12max = 2
         else  if (index(calcstring,'j12max=3').ne.0) then
            j12max = 3
         else  if (index(calcstring,'j12max=4').ne.0) then
            j12max = 4
         else  if (index(calcstring,'j12max=5').ne.0) then
            j12max = 5
         else
            write(*,*) "*** ERROR: Input attempted to set j12max to value which is not 0,1,2,3,4 or 5. -- Exiting."
            stop
         end if    
         write(*,'(A,I4)') " Total angular momentum in (12) subsystem set to j12max = ",j12max
         write(dummy,'(I1)') j12max
         descriptors = trim(descriptors) // "j12max=" // dummy // "-"
      end if
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     write quadrature variables to terminal
c        radial integration routines & parameters used    

      write (*,*) "Radial Integration in (12) subsystem:"
      write (*,'(A,I4,A)') "   total of NP12 = NP12A+NP12B = ",NP12A+NP12B," points"
      write (*,126) NP12A
      write (*,127) NP12A/2,P12A
      write (*,128) NP12A/2,P12A,P12B
      write (*,129) NP12B,P12B,P12C
      write (*,'(A,I4,A)') "   NP12p == ( NP12A+NP12B = ",NP12A+NP12B,") cannot be dialled."
c        angular integration routines & parameters used 
      if (AngularType12.eq.1) then
         write (*,*) "Angular Integration in (12) subsystem:"
         write (*,*) "   by Gauss-Legendre for theta & phi separately"
         write (*,132) Nordth12, Nthbins12, Nordth12*Nthbins12
         write (*,140) Nordphi12, Nphibins12, Nordphi12* Nphibins12
         write (*,142) (Nordth12*NthBins12)*(Nordphi12*NphiBins12)
         write (*,*) "           Nanggrid12 not used."
      else if (AngularType12.eq.2) then
         write (*,*) "Angular Integration in (12) subsystem"
         write (*,*) "   by Lebedev-Laikov for theta & phi combined:"
         Nanggrid12 = Nordth12
         if (Nanggrid12.gt.Nangmax) then
            write (*,*) '*** ERROR: Nanggrid12 (set to Nordth12) >',Nangmax,' too large -- Exiting.'
            stop
         end if
         write (*,144) Nanggrid12
         write (*,*) "   Nthbins12, Nordphi12, Nphibins12 have no meaning."
      else
         write (*,*) "*** ERROR: Illegal AngularType12 -- Exiting."
         stop
      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     formats
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 126  format (' ','   1 bin with hyperbolic map: (NP12A = ', I2,') points:')
 127  format (' ','      (NP12A/2 = ', I2,') points in p12 interval [0; (P12A = ',F10.4,')] fm^-1')
 128  format (' ','      (NP12A/2 = ', I2,') points in p12 interval [(P12A = ',F10.4,');(P12B = ',F10.4,')] fm^-1')
 129  format (' ','   1 bin with linear map: (NP12B=', I2,') points in p12 interval [(P12B = ',F10.4,');(P12C = ',F10.4,')] fm^-1')

 132  format (' ','   theta: (Nordth12 = ',I4,') points per (NthBins12 = ', I2,') bins = ',I16,' points.')
 140  format (' ','   phi:   (Nordphi12  = ',I4,') points per (NphiBins12  = ', I2,') bins = ',I16,' points.')
 142  format (' ','   ==> size of solid angle grid = (Nordth12*NthBins12)*(Nordphi12*NphiBins12) = ',I16,' points.')
 144  format (' ','   ==> preliminary size of solid angle grid Nanggrid12(set to Nordth12) = 'I4,' points.')

      return
      end
