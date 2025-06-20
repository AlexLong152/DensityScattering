c     parseFile.f - Fortran translation of pyVersion/parseFile.py
c     Author: Translated from Python version by alexl

      module parseFileData
      implicit none
      
c     Include common physics constants and parameters
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
      include '../common-densities/constants.def'
      
      integer, parameter :: MAXDATA = 50000
      
      integer :: totalEntries
      integer :: letterIndex(MAXDATA)
      integer :: twoI(MAXDATA)
      integer :: twoL(MAXDATA)  
      real*8 :: wcmValues(MAXDATA)
      real*8 :: delValues(MAXDATA)
      real*8 :: srValues(MAXDATA)
      
      logical :: dataInitialized = .false.
      
      contains
      
      integer function letterToIndex(letter)
      implicit none
      character*1 letter
      
      select case (letter)
         case ('S', 's')
            letterToIndex = 0
         case ('P', 'p')
            letterToIndex = 1
         case ('D', 'd')
            letterToIndex = 2
         case ('F', 'f')
            letterToIndex = 3
         case ('G', 'g')
            letterToIndex = 4
         case default
            letterToIndex = -1
      end select
      end function letterToIndex
      
      subroutine initializeFileData(filename, verbosity)
      implicit none
      character*(*) filename
      integer verbosity
      
      integer :: unit, ios
      character*200 :: line
      character*1 :: currentLetter
      integer :: currentLetterIdx, current2I, current2L
      real*8 :: wcm, del, sr, dummy
      
      if (dataInitialized) return
      
      totalEntries = 0
      currentLetterIdx = -1
      current2I = -1
      current2L = -1
      
      unit = 20
      open(unit=unit, file=filename, status='old', iostat=ios)
      if (ios .ne. 0) then
         write(*,*) 'ERROR: Could not open file: ', filename
         stop
      endif
      
      do
         read(unit, '(A)', iostat=ios) line
         if (ios .ne. 0) exit
         
         if (len_trim(line) .le. 5) cycle
         
c        Check for header line containing "PI-N"
         if (index(line, 'PI-N') .gt. 0) then
            call parseHeader(line, currentLetter, current2I, current2L)
            currentLetterIdx = letterToIndex(currentLetter)
            if (verbosity .ge. 2 .and. currentLetterIdx .ge. 0) then
               write(*,*) 'Header: ', currentLetter, current2I, 
     &                   current2L
            endif
            cycle
         endif
         
c        Skip column header lines
         if (index(line, 'WCM') .gt. 0) cycle
         if (index(line, 'Del') .gt. 0) cycle
         
c        Try to parse data line
         if (currentLetterIdx .ge. 0 .and. current2I .ge. 0 .and. 
     &       current2L .ge. 0) then
            read(line, *, iostat=ios) wcm, del, dummy, sr
            if (ios .eq. 0) then
               if (verbosity .ge. 3) then
                  write(*,*) 'Data: ', currentLetter, wcm, del, sr
               endif
               
               if (totalEntries .lt. MAXDATA) then
                  totalEntries = totalEntries + 1
                  letterIndex(totalEntries) = currentLetterIdx
                  twoI(totalEntries) = current2I
                  twoL(totalEntries) = current2L
                  wcmValues(totalEntries) = wcm
                  delValues(totalEntries) = del * Pi / 180.0d0
                  srValues(totalEntries) = sr
               endif
            endif
         endif
      enddo
      
      close(unit)
      dataInitialized = .true.
      
      ! if (verbosity .ge. 1) then
      !    write(*,*) 'Total entries: ', totalEntries
      ! endif
      
      end subroutine initializeFileData
      
      subroutine parseHeader(line, letter, n2I, n2L)
      implicit none
      character*(*) line
      character*1 letter
      integer n2I, n2L
      
      integer i, pos
      character*3 qnStr
      
      letter = '?'
      n2I = -1
      n2L = -1
      
      pos = index(line, 'PI-N')
      if (pos .gt. 0) then
         do i = pos + 4, len_trim(line) - 2
            if (line(i:i) .eq. ' ') then
               letter = line(i+1:i+1)
               if (letterToIndex(letter) .ge. 0) then
                  qnStr = line(i+1:i+3)
                  read(qnStr(2:2), '(I1)', iostat=pos) n2I
                  read(qnStr(3:3), '(I1)', iostat=pos) n2L
                  return
               endif
            endif
         enddo
      endif
      
      end subroutine parseHeader
      
      subroutine getScatteringData(targetLetterIndex, target2I, 
     &     target2L, wcm_target, del_result, sr_result, found)
      implicit none
      integer targetLetterIndex, target2I, target2L
      real*8 wcm_target, del_result, sr_result
      logical found
      
      integer :: i, best_match
      real*8 :: wcm_diff, min_diff
      
      found = .false.
      del_result = 0.0d0
      sr_result = 0.0d0
      best_match = 0
      min_diff = 1.0d30
      
      if (.not. dataInitialized) return
      
      do i = 1, totalEntries
         if (letterIndex(i) .eq. targetLetterIndex .and.
     &       twoI(i) .eq. target2I .and.
     &       twoL(i) .eq. target2L) then
            
            wcm_diff = abs(wcmValues(i) - wcm_target)
            if (wcm_diff .lt. min_diff) then
               min_diff = wcm_diff
               best_match = i
            endif
            
            if (wcm_diff .lt. 1.0d-6) exit
         endif
      enddo
      
      if (best_match .gt. 0 .and. min_diff .lt. 50.0d0) then
         del_result = delValues(best_match)
         sr_result = srValues(best_match)
         found = .true.
      endif
      
      end subroutine getScatteringData
      
      subroutine getScatteringDataByLetter(letter, target2I, target2L, 
     &     wcm_target, del_result, sr_result, found)
      implicit none
      character*1 letter
      integer target2I, target2L
      real*8 wcm_target, del_result, sr_result
      logical found
      character*32 :: arg1, arg2, arg3, arg4
      character*200 :: cmd
      character*50 :: outname
      
c     Initialize variables
      found = .false.
      del_result = 0.0d0
      sr_result = 0.0d0
      
      call getScatteringData(letterToIndex(letter), target2I, target2L, 
     &     wcm_target, del_result, sr_result, found)

      write(arg1,'(I0)') letterToIndex(letter)
      write(arg2,'(I0)') target2I
      write(arg3,'(I0)') target2L
      write(arg4,'(F0.2)') wcm_target
      
      write(*,'(A1,I1,I1,A,F8.2, A)') letter, target2I, target2L, 
     &     ' at ', wcm_target, ' MeV'

      cmd = 'python3 pyVersion/parseFile.py ' // trim(arg1) // ' ' // 
     &      trim(arg2) // ' ' // trim(arg3) // ' ' // trim(arg4) 
      call EXECUTE_COMMAND_LINE(cmd)
      if(found) then
        write(*,'(A,ES12.5,A,ES12.5)') 'Fortran result: del=',
     &   del_result, ' sr=', sr_result
        write(*,'(A)') '---------------------------------------------'
      else
        write(*,'(A)') 'Fortran result - Entry not found'
        write(*,'(A)') '---------------------------------------------'
      endif
      end subroutine getScatteringDataByLetter
      
      end module parseFileData
      
