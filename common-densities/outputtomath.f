cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for One/Twobody Contributions to Few-Nucleon Processes Calculated Via 1N/2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              outputtomath         : output to stdout in mathematica-frendly format
c              ConvertRealToMath    : convert a string which used to be a Fortran Real into a mathematica-readable string
c              ConvertComplexToMath : convert a string which used to be a Fortran Complex into a mathematica-readable string
c              StripSpaces          : strip spaces from string
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, identical to file of same name in common-densities/ of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility
c     Order of output changed from Conpton code: run first over extQnum=1 starting at Mzp=max and Mz=mad, then Mz=max-1 etc,
c     then Mzp=max-1 etc, until extQnum=1, Mzp=-max, Mz=-max;
c     -- then continue with extQnum=2 by same fashion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outputtomath(Result,twoSnucl,extQnumlimit,verbosity)
c**********************************************************************
      IMPLICIT NONE
c**********************************************************************
c     input variables

      integer,intent(in)    :: extQnumlimit
      integer,intent(in)    :: twoSnucl
      complex*16,intent(in) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      
      integer,intent(in)    :: verbosity
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     intrinsic variables

      integer :: extQnum,twoMzp,twoMz
c     for mathematica-friendly output, define numbers as strings. not elegant, but works
      character(len=64) string
      character(len=4*(twoSnucl+1)**2*68) longstring
      longstring = "" ! initialise auxiliary
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      if (verbosity.eq.1000) continue
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c output ALL amplitudes when symmetry NOT invoked        
      write(*,'(A,I2,A,I2,A,I4,A)') "Mathematica-friendly output: all ",extQnumlimit,
     &     "*",(twoSnucl+1)**2," = ", extQnumlimit*(twoSnucl+1)**2," amplitudes (some related by symmetries):"
      write(*,*)    "   [Sequence: extQnum∈[1;extQnumlimit] (outmost), Mzp∈[S;-S], Mz∈[S;-S] (innermost)]"
      string = ""               ! initialise
      do extQnum=1,extQnumlimit
         do  twoMzp=twoSnucl,-twoSnucl,-2
            do twoMz=twoSnucl,-twoSnucl,-2
               write(string,'(SP,"(",E24.18,",",E24.18,")")') Result(extQnum,twoMzp,twoMz)
               call ConvertComplexToMath(string)
               longstring = trim(adjustl(longstring)) // string // ","
               call StripSpaces(longstring)
            end do              ! twoMz
         end do                 ! twoMzp
      end do                    ! extQnum

      longstring = '{' // trim(adjustl(longstring)) // '}'    
      longstring = longstring(:index(longstring,",}")-1) // "}"
      write(*,*) '        ',trim(adjustl(longstring))    
      
      return
      
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ConvertRealToMath(string)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     convert a string which used to be a Fortran Real into a mathematica-readable string,
c     e.g. 1.2345E003 => 1.2345*10*(003)
      character(len=*),intent(inout) :: string
      if ( index(string,"E").ne.0 ) then
         string = trim(string(1:index(string,"E")-1) // '*10^(' // string(index(string,"E")+1:)) // ')'
      else
         string = trim(string)
      end if
      string = trim(adjustl(string))
      call StripSpaces(string)
      end ! subroutine ConvertRealToMath

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ConvertComplexToMath(string)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     convert a string which used to be a Frortran Compex into a mathematica-readable string,
c     e.g. (1.2E3,3.4E3) => 1.2*10^3+3.4*10^3*I
      character(len=*),intent(inout) :: string
      character(len=30) :: restring,imstring
      write(restring,*) trim(adjustl(string(index(string,"(")+1:index(string,",")-1)))
      write(imstring,*) trim(adjustl(string(index(string,",")+1:index(string,")")-1)))
      call ConvertRealToMath(restring)
      call ConvertRealToMath(imstring)
      string = trim(adjustl(restring // imstring // "*I"))
      call StripSpaces(string)
      end ! subroutine ConvertComplexToMath
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine StripSpaces(string)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      character(len=*),intent(inout) :: string
      integer :: stringLen 
      integer :: last, actual
      
      stringLen = len (string)
      last = 1
      actual = 1
      
      do while (actual < stringLen)
         if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
         else
            last = last + 1
            if (actual.lt.last) actual = last
         endif
      end do
      
      end ! subroutine StripSpaces

      
