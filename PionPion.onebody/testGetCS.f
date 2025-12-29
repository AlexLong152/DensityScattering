c     ================================================================
c     Test program to call getCS and output results
c     ================================================================
      program testGetCS
      use pionScatLib
      implicit none
      double precision sqrtS, sqrtSReal, x, theta, mNucl
      integer isospin, piCharge
      double precision CrossSec

c     Initialize the SAID data file
      call initializeFileData('said-pi.txt', 0)

c     Read input parameters from stdin
      read(*,*) sqrtS, theta, isospin, piCharge

c     Convert theta to cos(theta)
      x = cos(theta * 3.14159265358979323846d0 / 180.0d0)

c     Set nucleon mass based on isospin
      if (isospin .eq. -1) then
         mNucl = Mneutron
      else if (isospin .eq. 1) then
         mNucl = Mproton
      endif

c     Call getCS to compute cross section
      call getCS(sqrtS, x, isospin, piCharge, CrossSec, mNucl)

c     Output result in a format easy to parse from Python
      write(*,'(A,E23.15)') 'CrossSec=', CrossSec/2.d0

      end program testGetCS
