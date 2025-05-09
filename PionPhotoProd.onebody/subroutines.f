      subroutine strip_string(string)
c     Strips leading and trailing whitespace from a string
      character*(*) string
      
c     Local variables
      integer i, j, len_str
      
c     Strip leading whitespace
      len_str = len(string)
      i = 1
      do while (i .le. len_str .and. 
     &        (string(i:i) .eq. ' ' .or. string(i:i) .eq. char(9)))
         i = i + 1
      enddo
      
c     If i > 1, shift string to the left
      if (i .gt. 1) then
         do j = 1, len_str - i + 1
            string(j:j) = string(j+i-1:j+i-1)
         enddo
         do j = len_str - i + 2, len_str
            string(j:j) = ' '
         enddo
      endif
      
c     Strip trailing whitespace
      len_str = len_trim(string)
      string(len_str+1:) = ''
      
      return
      end

      subroutine add_unique(array, count, value)
c     Adds a value to an array if it's unique
      character*10 array(*)
      integer count
      character*10 value
      
c     Local variables
      integer i
      logical found
      
c     Check if value is already in array
      found = .false.
      do i = 1, count
         if (array(i) .eq. value) then
            found = .true.
            goto 10
         endif
      enddo
      
   10 continue
c     Add if not found
      if (.not. found) then
         count = count + 1
         array(count) = value
      endif
      
      return
      end

      subroutine parse_gwu_line(line, eg_val, rc_val, re_val, 
     &                          ic_val, ie_val, status)
c     Parses a line in GWU format like "182.00   14.94(  0.57)    1.43(  0.06)"
      character*500 line
      double precision eg_val, rc_val, re_val, ic_val, ie_val
      integer status
      
c     Local variables
      integer pos, open_pos, close_pos, len_line
      character*20 value_str, err_str
      
c     Initialize
      status = 1  ! Error
      eg_val = 0.0
      rc_val = 0.0
      re_val = 0.0
      ic_val = 0.0
      ie_val = 0.0
      
      len_line = len_trim(line)
      if (len_line .lt. 10) return  ! Not enough characters
      
c     First, extract the energy value
      pos = 1
      do while (pos .le. len_line .and. 
     &        (line(pos:pos) .eq. ' ' .or. line(pos:pos) .eq. char(9)))
         pos = pos + 1
      enddo
      
      if (pos .gt. len_line) return  ! End of line
      
c     Extract energy value
      read(line(pos:), *, iostat=status) eg_val
      if (status .ne. 0) return
      
c     Advance to first number after energy
      do while (pos .le. len_line .and. 
     &        (line(pos:pos) .ne. ' ' .and. line(pos:pos) .ne. char(9)))
         pos = pos + 1
      enddo
      
      do while (pos .le. len_line .and. 
     &        (line(pos:pos) .eq. ' ' .or. line(pos:pos) .eq. char(9)))
         pos = pos + 1
      enddo
      
      if (pos .gt. len_line) return  ! End of line
      
c     Find the opening parenthesis for the real part error
      open_pos = index(line(pos:), '(')
      if (open_pos .eq. 0) return  ! No open parenthesis
      open_pos = open_pos + pos - 1
      
c     Extract the real center value
      value_str = line(pos:open_pos-1)
      read(value_str, *, iostat=status) rc_val
      if (status .ne. 0) return
      
c     Find the closing parenthesis
      close_pos = index(line(open_pos:), ')')
      if (close_pos .eq. 0) return  ! No close parenthesis
      close_pos = close_pos + open_pos - 1
      
c     Extract the real error value
      err_str = line(open_pos+1:close_pos-1)
      read(err_str, *, iostat=status) re_val
      if (status .ne. 0) return
      
c     Advance to the imaginary part
      pos = close_pos + 1
      do while (pos .le. len_line .and. 
     &        (line(pos:pos) .eq. ' ' .or. line(pos:pos) .eq. char(9)))
         pos = pos + 1
      enddo
      
      if (pos .gt. len_line) return  ! End of line
      
c     Find the opening parenthesis for the imaginary part error
      open_pos = index(line(pos:), '(')
      if (open_pos .eq. 0) return  ! No open parenthesis
      open_pos = open_pos + pos - 1
      
c     Extract the imaginary center value
      value_str = line(pos:open_pos-1)
      read(value_str, *, iostat=status) ic_val
      if (status .ne. 0) return
      
c     Find the closing parenthesis
      close_pos = index(line(open_pos:), ')')
      if (close_pos .eq. 0) return  ! No close parenthesis
      close_pos = close_pos + open_pos - 1
      
c     Extract the imaginary error value
      err_str = line(open_pos+1:close_pos-1)
      read(err_str, *, iostat=status) ie_val
      if (status .ne. 0) return
      
c     Success
      status = 0
      
      return
      end