cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for One/Twobody Contributions to Few-Nucleon Processes Calculated Via 1N/2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              ReadinputCommon         : process input for parameters which are identical for onebody and twobody
c              ReadinputCommonComments : process input comments (last line of input file)
c     CONTAINS FUNCTIONS:
c              stringtolower           : function to convert string to all-lowercase
c              halfinteger             : function to divide integer by 2 to string representing half-integer
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, based on file of same name in common-densities/ of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c      
c     subroutines to read input files for density calculations:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ReadinputCommon(Elow,Ehigh,Einterval,
     &     thetaLow,thetaHigh,thetaInterval,
     &     outfile,descriptors,densityFileName,inUnitno,
     &     nucleus,Anucl,twoSnucl,Mnucl,extQnumlimit,
     &     symmetry,verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This routine reads in the values of the external photon momentum &
c     scattering angle in the lab frame. It takes them from the I/O
c     unit denoted by inUnitno, where they are stored in MeV and degrees.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
c     
c     
c     
c     VARIABLES PASSED OUT:
c     
      real*8,intent(out)        :: Elow,Ehigh,Einterval
      real*8,intent(out)        :: thetaLow,thetaHigh,thetaInterval
      character*500,intent(out) :: outfile
      character*200,intent(out) :: descriptors ! additional descriptors of calculation for outputfilename
      character*500,intent(out) :: densityFileName

      character*3, intent(out)  :: nucleus ! name of nucleus to be considered, for output file name
      integer,intent(out) :: Anucl               ! mass number of target nucleus
      integer,intent(out) :: twoSnucl            ! spin of target nucleus
      real*8,intent(out)  :: Mnucl               ! mass of target nucleus
      
      integer,intent(out) :: extQnumlimit ! number of combined external quantum numbers of in and out state
      integer,intent(out) :: symmetry     ! whether or not a specific package Usesymmetry() should be used for process -- NOT YET IMPLEMENTED!!!
      
      integer,intent(out) :: verbosity
      
      character*2  dummied ! a variable only for converting the integer extQnumlimit to a character string
c     
c     ---------------------------------------------------------------
c     BS: Ultimately variable def's in Readme?
c     ---------------------------------------------------------------
c     
c     Elow, Ehigh - beam-in low, beam-in high
c     Einterval   - step to move from Elow to Ehigh
c     
c     thetaLow, thetaHigh - low and high limits of theta angle
c     thetaInterval       - step to move from thetaLow-->thetaHigh       
c     
c     outfile - name of file to which output is to be written
c     *******************************************************************
c     
c     VARIABLES PASSED IN:
c     
      integer,intent(in) :: inUnitno ! I/O unit containing all the input information
c     
c     *******************************************************************
c     
c     INTERNAL VARIABLES:
c     
      character*500 calcstring ! holds info on type of calculation being done hgrie Aug 2020: increased from 80
      character*500 stringtolower ! function to convert string to all-lowercase, for comparisons -- defined at file end
      character*5 halfinteger     ! function to divide integer by 2 to string representing half-integer-- defined at end
c
c     *******************************************************************
c     
c     set logical variables to FALSE, other parameters to "0": defines default values!
      descriptors = "-"
      
c     Read in the input file line-by-line, write output info to terminal
c     
      read (inUnitno,*) Elow,Ehigh,Einterval
      write (*,*) '   Lowest beam energy           = ', Elow, 'MeV'
      write (*,*) '   Highest beam energy          = ', Ehigh, 'MeV'
      write (*,*) '   Energy step from low to high = ', Einterval, 'MeV'
      read (inUnitno,*) thetaLow,thetaHigh,thetaInterval
      write (*,*) '   Lowest theta                 = ', thetaLow, 'deg'
      write (*,*) '   Highest theta                = ', thetaHigh, 'deg'
      write (*,*) '   Theta step from low to high  = ', thetaInterval, 'deg'

      read (inUnitno,'(A500)') outfile
      outfile = TRIM(outfile)
      write (*,*) 'First take on output filename (placeholders to be replaced below): '
      write (*,*) '      ',TRIM(outfile)
      
c     density filename and determination of target nucleus
      read (inUnitno, *) densityFileName
      write (*,*) 'First take on density filename (placeholders to be replaced below): '
      write (*,*) '      ',TRIM(densityFileName)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     determine target nucleus from density filename, covering 3HE, 3he, 3He,...
c     Mnucl = 0. ! initialise Mnucl so that we can later ask if it was adapted by target determination
      if (index(stringtolower(densityFileName),'3he').ne.0) then
         nucleus = "3He"
         Anucl = 3
         twoSnucl = 1           ! 2*target spin
         Mnucl = M3He
         write (*,'(A,I2,A,A,A,G13.7,A)') '   Target Nucleus: 3He, Anucl = ',Anucl,
     &        ', spin Snucl = ', halfinteger(twoSnucl),', mass Mnucl = ',Mnucl,'MeV'
      else if (index(stringtolower(densityFileName),'3h').ne.0) then
         nucleus = "3H"
         Anucl = 3
         twoSnucl = 1           ! 2*target spin
c         Mnucl = M3H
         write (*,'(A,I2,A,A,A,G13.7,A)') '   Target Nucleus: 3H, Anucl = ',Anucl,
     &        ', spin Snucl = ', halfinteger(twoSnucl),', mass Mnucl = ',Mnucl,'MeV'
      else if ((index(stringtolower(densityFileName),'deuteron').ne.0).or.(index(stringtolower(densityFileName),'2h').ne.0)) then
         nucleus = "2H"
         Anucl = 2
         twoSnucl = 2           ! 2*target spin
         Mnucl = Md
         write (*,'(A,I2,A,A,A,G13.7,A)') '   Target Nucleus: Deuteron, Anucl = ',Anucl,
     &        ', spin Snucl = ', halfinteger(twoSnucl),', mass Mnucl = ',Mnucl,'MeV'
      else if (index(stringtolower(densityFileName),'4he').ne.0) then
         nucleus = "4He"
         Anucl = 4
         twoSnucl = 0           ! 2*target spin
         Mnucl = M4He
         write (*,'(A,I2,A,A,A,G13.7,A)') '   Target Nucleus: 4He, Anucl = ',Anucl,
     &        ', spin Snucl = ', halfinteger(twoSnucl),', mass Mnucl = ',Mnucl,'MeV'
      else if (index(stringtolower(densityFileName),'6li').ne.0) then
         nucleus = "6Li"
         Anucl = 6
         twoSnucl = 2           ! 2*target spin
         Mnucl = M6Li
         write (*,'(A,I2,A,A,A,G13.7,A)') '   Target Nucleus: 6Li, Anucl = ',Anucl,
     &        ', spin Snucl = ', halfinteger(twoSnucl),', mass Mnucl = ',Mnucl,'MeV'
c     I could also implement a routine "search string for A=7 and Z=3"...
      else
         stop "*** ERROR: Could not determine target nucleus from density filename: Abort."
      end if

      if (Mnucl.eq.0.) stop "*** ERROR: Target Nucleus Mass not defined: Abort."
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     finally, read string of common parameters
      read (inUnitno,'(A80)') calcstring

c      set up different levels of verbosity for debugging
      if (index(calcstring,'verbose').eq.0) then
         verbosity = 0
      else if (index(calcstring,'verbose4').ne.0) then 
         verbosity = 4
         write (*,*) '********** Verbose Mode 4 **********'
      else if (index(calcstring,'verbose3').ne.0) then
         verbosity = 3
         write (*,*) '********** Verbose Mode 3 **********'
      else if (index(calcstring,'verbose2').ne.0) then
         verbosity = 2
         write (*,*) '********** Verbose Mode 2 **********'
      else 
         verbosity = 1
c     print1Namps=.true.
         write (*,*) '********** Verbose Mode 1 **********'
      end if

      write (*,*) 'Input and output in centre-of-mass frame.'
      
c     Determine number of external quantum numbers
      
      if (index(calcstring,'extQnumlimit=').eq.0) then
         extQnumlimit = 4             ! default to Compton case: 2 circpols in, 2 circpols out
         write(*,'(A,I4)') " Number of external quantum numbers not specified -- using default extQnumlimit =",extQnumlimit
      else if (index(calcstring,'extQnumlimit=').ne.0) then
c     first need to work over double-digist, then single-digits -- or else "extQnumlimit=16" is set to 1, not 16.
         if (index(calcstring,'extQnumlimit=11').ne.0) then
            extQnumlimit = 11
         else if (index(calcstring,'extQnumlimit=12').ne.0) then
            extQnumlimit = 12
         else if (index(calcstring,'extQnumlimit=13').ne.0) then
            extQnumlimit = 13
         else if (index(calcstring,'extQnumlimit=14').ne.0) then
            extQnumlimit = 14
         else if (index(calcstring,'extQnumlimit=15').ne.0) then
            extQnumlimit = 15
         else if (index(calcstring,'extQnumlimit=16').ne.0) then
            extQnumlimit = 16
         else if (index(calcstring,'extQnumlimit=1').ne.0) then
            extQnumlimit = 1
         else if (index(calcstring,'extQnumlimit=2').ne.0) then
            extQnumlimit = 2
         else if (index(calcstring,'extQnumlimit=3').ne.0) then
            extQnumlimit = 3
         else if (index(calcstring,'extQnumlimit=4').ne.0) then
            extQnumlimit = 4
         else if (index(calcstring,'extQnumlimit=5').ne.0) then
            extQnumlimit = 5
         else  if (index(calcstring,'extQnumlimit=6').ne.0) then
            extQnumlimit = 6
         else if (index(calcstring,'extQnumlimit=7').ne.0) then
            extQnumlimit = 7
         else if (index(calcstring,'extQnumlimit=8').ne.0) then
            extQnumlimit = 8
         else if (index(calcstring,'extQnumlimit=9').ne.0) then
            extQnumlimit = 9
         else
            write(*,*) "*** ERROR: Input attempted to set extQnumlimit to number not natural and between 1 and 16. -- Exiting."
            stop
         end if    
         write(*,'(A,I4)') " Number of external quantum numbers set to extQnumlimit = ",extQnumlimit
         write(dummied,'(I2.2)') extQnumlimit
         descriptors = trim(descriptors) // "extQnumlimit=" // dummied // "-"
      end if
c      
c     Determine if/which symmetry subroutine should be used (must be specific to process! -- set up inside kernel directory)
c
      if (index(calcstring,'usesymmetry1').ne.0) then
         symmetry=1
      else if (index(calcstring,'usesymmetry2').ne.0) then
         symmetry=2
      else if (index(calcstring,'usesymmetry3').ne.0) then
         symmetry=3
      else if (index(calcstring,'usesymmetry4').ne.0) then
         symmetry=4
      else if (index(calcstring,'usesymmetry5').ne.0) then
         symmetry=5
      else if (index(calcstring,'usesymmetry6').ne.0) then
         symmetry=6
      else
         symmetry=0
         write(*,'(A,I4,A)') " Calculate all",extQnumlimit*(twoSnucl+1)**2,
     &        " amplitudes independently, not using symmetry() of process."
      end if
      
      if (symmetry.ne.0) then
         write(*,'(A,I3)') "Calculate amplitudes using symmetry() of process, variant number ",symmetry
         write(*,*) " === Usesymmetry() NOT YET IMPLEMENTED. ==="
      else
c     add here call to a subroutine for specific message what symmetry implemented -- needs to be part of kernel directory!
c        call symmetrymessage(symmetry)
         continue
      end if
      
c add symmetry status to descriptor      
      write(dummied,'(I1)') symmetry
      descriptors = trim(descriptors) // "symmetry=" // trim(dummied) // "-"
      
      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ReadinputCommonComments(descriptors,inUnit,verbosity)     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie Aug 2020
c     Now the routine which reads last line as comment and appends to descriptors
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     VARIABLES PASSED INOUT:
      character*200,intent(inout) :: descriptors ! additional descriptors of calculation for outputfilename
      
c     VARIABLES PASSED OUT:
      
c     VARIABLES PASSED IN:
      integer,intent(in)     :: inUnit
      integer,intent(in)     :: verbosity
c     
c     INTERNAL VARIABLES:
      character*200 comments
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      if (verbosity.eq.1000) continue ! keep for future use
      
      read (inUnit,*) comments
      write (*,*) "ADDITIONAL ",trim(comments)
c     add comments to end of descriptors
      descriptors = trim(descriptors) // "." // comments(11:) ! strips out "COMMENTS:_"

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     function to convert string to all-lowercase, hgrie May 2018
c     from https://groups.google.com/forum/#!msg/comp.lang.fortran/CKx1L2Ahkxg/HH_kMoHAffcJ
      
      function stringtolower( string ) result (new)
      character(len=500)           :: string

      character(len=len(string)) :: new
      
      integer                    :: i
      integer                    :: k

      length = len(string)
      new    = string
      do i = 1,len(string)
         k = iachar(string(i:i))
         if ( k >= iachar('A') .and. k <= iachar('Z') ) then
            k = k + iachar('a') - iachar('A')
            new(i:i) = achar(k)
         endif
      enddo
      end function stringtolower 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     function to make half-integer -- self-made hgrie May 2018

      function halfinteger(number) result (string)
      integer :: number
      character*5 string

      if (mod(number,2).eq.0) then
         write (string,'(I2)') number/2
      else
         write (string,'(I2,A)') number,"/2"
      end if
      end function halfinteger
