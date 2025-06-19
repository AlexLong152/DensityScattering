c     Global variables to cache the SAID data
      module said_data_cache
        implicit none
        
c       Flag to track if the SAID data has been loaded
        logical, save :: said_data_loaded = .false.
        
c       Maximum dimensions for the cache arrays
        integer, parameter :: MAX_WAVES = 50     ! Maximum number of waves to store
        integer, parameter :: MAX_ENERGIES = 200 ! Maximum energies per wave
        
c       Wave information
        character(len=3), save :: wave_labels(MAX_WAVES)
        character(len=2), save :: target_parts(MAX_WAVES)
        integer, save :: wave_count = 0
        
c       Energy and amplitude data for each wave
        double precision, save :: energy_data(MAX_WAVES, MAX_ENERGIES)
        double complex, save :: amplitude_data(MAX_WAVES, MAX_ENERGIES)
        integer, save :: energy_counts(MAX_WAVES)
        
c       Constants
        include '../common-densities/constants.def'
        
c       Nucleon mass parameter (needed for CM energy conversion)
        double precision, parameter :: mN = Mnucleon  ! Average nucleon mass from constants.def
        
c       Units conversion factor
        double precision, parameter :: UNITS_FACTOR = 1.0d0 / (HC * 1000.0d0)
      end module said_data_cache

c     Implementation of the getPoles function that uses cached SAID data
      subroutine getPoles(target_str, ell, sqrtS,
     &      Eplus, Mplus, Eminus, Mminus)
c-----------------------------------------------------------------------
c     Retrieves the E and M amplitudes for a given target, angular momentum ell,
c     and center-of-mass energy sqrtS by using cached data from the said-SM22.txt file.
c     
c     Based on the Python implementation in readData.py, particularly the functions:
c     - parseSpinString
c     - buildspinstring
c     - getPoles
c
c     Parameters:
c         target_str - String, one of 'p12', 'n12', '32q'
c         ell - Integer, the angular momentum (0 for S-wave, 1 for P-wave, etc.)
c         sqrtS - Double precision, the center-of-mass energy in MeV
c
c     Returns (through arguments):
c         Eplus, Mplus, Eminus, Mminus - Complex values representing amplitudes
c
      use said_data_cache
      implicit none
      
      character*10 target_str
      integer ell
      double precision sqrtS
      double complex Eplus, Mplus, Eminus, Mminus
      
c     Constants
      integer UNITNO
      parameter (UNITNO = 15)
      
c     Local variables for pattern matching
      character*20 spinStrings(4) ! Will store the spin strings for each amplitude type
      character*20 wavePart
      character*10 target_prefix
      real isospin
      
c     Variables for finding closest energy
      double precision minDiff, currDiff
      integer minIndex, waveIndex
      
c     Loop indices
      integer i, j, k
      
c     Initialize output values to zero
      Eplus = dcmplx(0.0d0, 0.0d0)
      Mplus = dcmplx(0.0d0, 0.0d0)
      Eminus = dcmplx(0.0d0, 0.0d0)
      Mminus = dcmplx(0.0d0, 0.0d0)
c     Units conversion factor is now defined in the module
      
c     Check if we need to load the SAID data
      if (.not. said_data_loaded) then
        call load_said_data()
      endif
      
c     Build spin strings using Python-like logic
c     First, determine isospin and target_prefix from target_str
      if (target_str .eq. 'p12') then
        isospin = 0.5
        target_prefix = 'p'
      else if (target_str .eq. 'n12') then
        isospin = 0.5
        target_prefix = 'n'
      else if (target_str .eq. '32q') then
        isospin = 1.5
        target_prefix = 'p'  ! Default to proton for 3/2
      else
        write(*,*) 'ERROR: Invalid target:', target_str
        write(*,*) 'Valid targets are: p12, n12, 32q'
        stop
      endif
      
c     Build the spin strings for each amplitude type using buildspinstring logic
c     For ell=0 (S-wave), there are no "minus" amplitudes
c     
c     1. Eplus
      call buildSpinString('plus', ell, isospin, 
     &                    target_prefix // 'E', spinStrings(1))
      
c     2. Mplus
      call buildSpinString('plus', ell, isospin, 
     &                    target_prefix // 'M', spinStrings(2))
      
c     3. Eminus - Only valid for ell > 0
      if (ell .gt. 0) then
        call buildSpinString('minus', ell, isospin, 
     &                      target_prefix // 'E', spinStrings(3))
      else
        spinStrings(3) = '???'
      endif
      
c     4. Mminus - Only valid for ell > 0
      if (ell .gt. 0) then
        call buildSpinString('minus', ell, isospin, 
     &                      target_prefix // 'M', spinStrings(4))
      else
        spinStrings(4) = '???'
      endif
      
c     Extract just the wave part (first 3 chars) from each spin string
c     This is what we'll match against the file headers
c     IMPORTANT: Convert to uppercase for matching with the file
      do i = 1, 4
         if (len_trim(spinStrings(i)) .ge. 3) then
            wavePart = spinStrings(i)(1:3)
            
c           Convert to uppercase 
            call strToUpper(wavePart)
            
            spinStrings(i) = wavePart
         endif
      enddo
      
c     For each amplitude type, search the cached data to find the closest energy
      do i = 1, 4
        if (len_trim(spinStrings(i)) .ge. 3) then
          do j = 1, wave_count
c           Match the wave part (e.g., S11, P33)
            if (spinStrings(i) .eq. wave_labels(j)) then
c             Match the target part (e.g., pE, nM)
              if (target_str .eq. '32q' .or. 
     &            target_parts(j)(1:1) .eq. target_prefix) then
                if (i .eq. 1 .or. i .eq. 3) then
c                 For E amplitudes, check for 'E' in target part
                  if (target_parts(j)(2:2) .eq. 'E') then
                    waveIndex = j
                    goto 10  ! Found a match, proceed to energy search
                  endif
                else if (i .eq. 2 .or. i .eq. 4) then
c                 For M amplitudes, check for 'M' in target part
                  if (target_parts(j)(2:2) .eq. 'M') then
                    waveIndex = j
                    goto 10  ! Found a match, proceed to energy search
                  endif
                endif
              endif
            endif
          enddo
          
c         No match found for this amplitude type
          cycle
          
c         Found a match, now find the closest energy
10        if (energy_counts(waveIndex) .gt. 0) then
            minDiff = abs(energy_data(waveIndex, 1) - sqrtS)
            minIndex = 1
            
            do k = 2, energy_counts(waveIndex)
              currDiff = abs(energy_data(waveIndex, k) - sqrtS)
              if (currDiff .lt. minDiff) then
                minDiff = currDiff
                minIndex = k
              endif
            enddo
            
c           Assign to the appropriate output variable
            if (i .eq. 1) then
              Eplus = amplitude_data(waveIndex, minIndex)
            else if (i .eq. 2) then
              Mplus = amplitude_data(waveIndex, minIndex)
            else if (i .eq. 3) then
              Eminus = amplitude_data(waveIndex, minIndex)
            else if (i .eq. 4) then
              Mminus = amplitude_data(waveIndex, minIndex)
            endif
          endif
        endif
      enddo
      
c     If no data was found, print error message but don't stop the program
c     This allows for cases like ell=0 with minus amplitudes, which are physically impossible
      if (Eplus .eq. dcmplx(0.0d0, 0.0d0) .and.
     &    Mplus .eq. dcmplx(0.0d0, 0.0d0) .and.
     &    Eminus .eq. dcmplx(0.0d0, 0.0d0) .and.
     &    Mminus .eq. dcmplx(0.0d0, 0.0d0)) then
        write(*,*) 'WARNING: No amplitude data found for:'
        write(*,*) '  target = ', target_str
        write(*,*) '  ell = ', ell
        write(*,*) '  sqrtS = ', sqrtS
        write(*,*) 'This may be normal for certain combinations'
        write(*,*) 'Returning zero amplitudes'
      endif
      
c     Final values are returned through the arguments
      return
      end

c-----------------------------------------------------------------------
      subroutine load_said_data()
c     Load the SAID data from the text file into the global cache
      use said_data_cache
      implicit none
      
c     Constants
      integer UNITNO
      parameter (UNITNO = 15)
      
c     Local variables for file reading
      character*500 filename, line
      integer iostat, lineCount, i
      logical foundHeader, inDataBlock, inGWUBlock
      character*20 wavePart, targetPart
      integer currWaveIndex
      
c     Variables for amplitude data
      double precision labE, emr_val, emi_val
      double precision sqrtS_point
      
c     Initialize
      wave_count = 0
      energy_counts = 0
      wave_labels = ''
      target_parts = ''
      
c     Open the SAID data file
      filename = 'said-SM22.txt'
      open(unit=UNITNO, file=filename, status='old', iostat=iostat)
      if (iostat .ne. 0) then
        write(*,*) 'ERROR: Cannot open SAID data file:'
        write(*,*) '  ', filename
        write(*,*) 'Please check that the file exists and is readable.'
        stop
      endif
      
c     Read through the file and collect all data
      lineCount = 0
      foundHeader = .false.
      inDataBlock = .false.
      inGWUBlock = .false.
      currWaveIndex = 0
      
      do while (.true.)
        read(UNITNO, '(A)', iostat=iostat) line
        if (iostat .ne. 0) exit  ! End of file or error
        lineCount = lineCount + 1
        
c       Look for partial wave headers like "PI0P S11 pE" or "PI0N S11 nE"
        if (index(line, 'PI0P') .gt. 0 .or. 
     &      index(line, 'PI0N') .gt. 0) then
          foundHeader = .false.
          inDataBlock = .false.
          inGWUBlock = .false.
          
c         Extract the wave part (S11) and target part (pE/nE)
          call extractWaveAndTarget(line, wavePart, targetPart)
          
c         Check if we already have this wave in our cache
          currWaveIndex = 0
          do i = 1, wave_count
            if (wave_labels(i) .eq. wavePart .and.
     &          target_parts(i) .eq. targetPart) then
              currWaveIndex = i
              exit
            endif
          enddo
          
c         If not found, add a new wave
          if (currWaveIndex .eq. 0) then
            if (wave_count .lt. MAX_WAVES) then
              wave_count = wave_count + 1
              currWaveIndex = wave_count
              wave_labels(currWaveIndex) = wavePart(1:3)  ! Explicitly use first 3 chars
              target_parts(currWaveIndex) = targetPart(1:2)  ! Explicitly use first 2 chars
              energy_counts(currWaveIndex) = 0
            else
              write(*,*) 'WARNING: Maximum number of waves exceeded'
              write(*,*) 'Ignoring wave:', wavePart, targetPart
              cycle
            endif
          endif
          
          foundHeader = .true.
          inDataBlock = .true.
          
        else if (inDataBlock) then
c         Check for GWU block
          if (index(line, 'GWU') .gt. 0 .and.
     &        index(line, 'Single energy values') .gt. 0) then
            inGWUBlock = .true.
          
c         Process data lines
          else if (len_trim(line) .gt. 0) then
c           Try to extract energy and amplitude values
            read(line, *, iostat=iostat) labE, emr_val, emi_val
            
            if (iostat .eq. 0) then
c             Convert lab energy to center-of-mass energy
              sqrtS_point = sqrt(mN**2 + 2.0*mN*labE)
              
c             Store the energy and amplitude values
              if (foundHeader .and. currWaveIndex .gt. 0) then
                if (energy_counts(currWaveIndex) .lt. MAX_ENERGIES) then
                  energy_counts(currWaveIndex) = 
     &                energy_counts(currWaveIndex) + 1
                  i = energy_counts(currWaveIndex)
                  energy_data(currWaveIndex, i) = sqrtS_point
                  amplitude_data(currWaveIndex, i) = 
     &                dcmplx(emr_val, emi_val) * UNITS_FACTOR
                else
                  write(*,*) 'WARNING: Maximum number of energies',
     &                       ' exceeded for wave:', 
     &                       wave_labels(currWaveIndex),
     &                       target_parts(currWaveIndex)
                endif
              endif
            endif
          endif
        endif
      enddo
      
c     Close the file
      close(UNITNO)
      
c     Mark the data as loaded
      said_data_loaded = .true.
      
      write(*,*) 'SAID data loaded with:', wave_count, 'waves'
      
      return
      end

c-----------------------------------------------------------------------
      subroutine buildSpinString(plusminus, ell, i, subchan, spinString)
c     This subroutine implements the Python buildspinstring function
c     
c     Multipole to spinstring, builds strings like 'p33pm', 's11pe', etc.
c     based on input parameters. This is the inverse of parseSpinString.
c
c     Parameters:
c         plusminus - String, either 'plus' or 'minus'
c         ell - Integer, the partial wave (0 for S-wave, 1 for P-wave, etc.)
c         i - Real, the isospin (0.5 or 1.5)
c         subchan - String, the subchannel (e.g., 'pE', 'nM')
c         spinString - Output string, the built spin string
c
c     Example mapping:
c         ell=1, i=1.5, plusminus='plus', subchan='pE' => 'p33pE'

      character*(*) plusminus, subchan
      integer ell
      real i
      character*(*) spinString
      
c     Local variables
      character*1 letter
      integer twoi, twoj
      
c     1. Map ell to letter (s,p,d,f,g) - lowercase for Python compatibility
      if (ell .eq. 0) then
        letter = 's'
      else if (ell .eq. 1) then
        letter = 'p'
      else if (ell .eq. 2) then
        letter = 'd'
      else if (ell .eq. 3) then
        letter = 'f'
      else if (ell .eq. 4) then
        letter = 'g'
      else
        write(*,*) 'ERROR: ell value not supported: ', ell
        letter = '?'
      endif
      
c     2. Calculate 2*i (twoi)
      twoi = nint(2.0 * i)  ! Nearest integer to 2*i
      
c     3. Calculate 2*j (twoj) based on ell and plusminus
      if (plusminus .eq. 'plus') then
        twoj = 2 * ell + 1  ! j = ell + 0.5
      else if (plusminus .eq. 'minus') then
        twoj = 2 * ell - 1  ! j = ell - 0.5
        
c       Check if twoj is valid (must be >= 1)
        if (twoj .lt. 1) then
          write(*,*) 'WARNING: Invalid j value for ell=', ell, 
     &               ' and sign=minus'
          write(*,*) 'This would give j < 0, which is unphysical'
          write(*,*) 'This is expected for ell=0 (S-wave)'
          spinString = '???'
          return
        endif
      else
        write(*,*) 'ERROR: plusminus must be "plus" or "minus", got:',
     &             plusminus
        spinString = '????'
        return
      endif
      
c     4. Construct partial wave label (e.g., 'p33')
      write(spinString, '(a1,i1,i1)') letter, twoi, twoj
      
c     5. Append subchannel
      spinString = trim(spinString) // trim(subchan)
      
c     Check for unphysical combinations (as in Python version)
      if (spinString .eq. 's11pm' .or.
     &    spinString .eq. 's11nm' .or.
     &    spinString .eq. 's31nm' .or.
     &    spinString .eq. 's31pm' .or.
     &    spinString .eq. 'p11pe' .or.
     &    spinString .eq. 'p11ne' .or.
     &    spinString .eq. 'p31pe' .or.
     &    spinString .eq. 'p31ne') then
        write(*,*) 'WARNING: Unphysical pole:', spinString
      endif
      
      return
      end

c-----------------------------------------------------------------------
      subroutine parseSpinString(spinString, plusMinus, ell, isospin, 
     &                          subChan, success)
c     This subroutine implements the Python parseSpinString function
c
c     Parse something like 'S11pE' or 'P33nM' into components
c     (plusMinus, ell, isospin, subChan)
c
c     Parameters:
c         spinString - Input string to parse (e.g., 'S11pE')
c         plusMinus - Output, 'plus' or 'minus'
c         ell - Output, integer (0 for S-wave, 1 for P-wave, etc.)
c         isospin - Output, real (0.5 or 1.5)
c         subChan - Output, the subchannel part (e.g., 'pE')
c         success - Output, logical indicating if parsing succeeded

      character*(*) spinString, plusMinus, subChan
      integer ell
      real isospin
      logical success
      
c     Local variables
      character*3 pw
      character*1 letter
      integer twoI, twoJ, pwlen, spinlen
      real J, diff
      
c     Initialize outputs
      success = .false.
      plusMinus = ''
      ell = -1
      isospin = 0.0
      subChan = ''
      
c     Validate input string length
      spinlen = len_trim(spinString)
      if (spinlen .lt. 4) then
        write(*,*) 'ERROR: spinString too short:', spinString
        return
      endif
      
c     Extract partial wave label (first 3 chars) and subChan
      pw = spinString(1:3)
      pwlen = len_trim(pw)
      if (pwlen .ne. 3) then
        write(*,*) 'ERROR: Invalid partial wave format:', pw
        return
      endif
      
      subChan = spinString(4:spinlen)
      
c     Extract letter, twoI, and twoJ from pw
      letter = pw(1:1)
      
c     Convert twoI and twoJ from characters to integers
      read(pw(2:2), '(i1)', err=100) twoI
      read(pw(3:3), '(i1)', err=100) twoJ
      
c     Map letter to ell
      if (letter .eq. 'S' .or. letter .eq. 's') then
        ell = 0
      else if (letter .eq. 'P' .or. letter .eq. 'p') then
        ell = 1
      else if (letter .eq. 'D' .or. letter .eq. 'd') then
        ell = 2
      else if (letter .eq. 'F' .or. letter .eq. 'f') then
        ell = 3
      else if (letter .eq. 'G' .or. letter .eq. 'g') then
        ell = 4
      else
        goto 100  ! Invalid letter
      endif
      
c     Calculate isospin and J
      isospin = twoI / 2.0
      J = twoJ / 2.0
      
c     Check if plusMinus is valid (J - ell = ±0.5)
      diff = J - ell
      if (abs(diff - 0.5) .lt. 0.01) then
        plusMinus = 'plus'
      else if (abs(diff + 0.5) .lt. 0.01) then
        plusMinus = 'minus'
      else
        goto 100  ! Invalid J-ell relationship
      endif
      
c     Success!
      success = .true.
      return
      
c     Error handling
  100 continue
      write(*,*) 'ERROR: Failed to parse spin string:', spinString
      success = .false.
      return
      end

c-----------------------------------------------------------------------
      subroutine extractWaveAndTarget(line, wavePart, targetPart)
c     Extracts the wave part (e.g., S11) and target part (e.g., pE)
c     from a header line like "        PI0P S11  pE        1/21/25"
      
      character*(*) line
      character*(*) wavePart, targetPart
      
c     Local variables
      integer i, pi0p_pos, length
      character*500 temp_line
      
c     Initialize
      wavePart = '     '
      targetPart = '  '
      
c     Find position of "PI0P" or "PI0N"
      pi0p_pos = index(line, 'PI0P')
      if (pi0p_pos .le. 0) then
        pi0p_pos = index(line, 'PI0N')
        if (pi0p_pos .le. 0) then
          write(*,*) 'ERROR: Header format error in SAID data file.'
          write(*,*) 'Expected line with PI0P or PI0N, but found: '
          write(*,*) line
          stop
        endif
      endif
      
c     Create a temporary string with the part after PI0P
      temp_line = line(pi0p_pos+4:)
      
c     Separate the parts using string parsing
      length = len_trim(temp_line)
      
c     Extract wave and target parts - using a simple string format
c     The format is "  S11  pE  " with varying number of spaces
      i = 1
      do while (i .le. length .and. temp_line(i:i) .eq. ' ')
        i = i + 1
      enddo
      
c     Now i points to the first non-space char, extract the next 3 chars
      if (i+2 .le. length) then
        wavePart = temp_line(i:i+2)
      endif
      
c     Skip to the next word
      i = i + 3
      do while (i .le. length .and. temp_line(i:i) .eq. ' ')
        i = i + 1
      enddo
      
c     Extract the next 2 chars for target part
      if (i+1 .le. length) then
        targetPart = temp_line(i:i+1)
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine strToUpper(str)
c     Converts a string to uppercase
      
      character*(*) str
      integer i, length
      
      length = len_trim(str)
      
      do i = 1, length
        if (str(i:i) .ge. "a" .and. str(i:i) .le. "z") then
          str(i:i) = char(ichar(str(i:i)) - ichar("a") + ichar("A"))
        endif
      enddo
      
      return
      end
