c     Test program for getRawM subroutine
c     Reads input from stdin, calls getRawM, outputs result to stdout
c
c     Input format (one line):
c       sqrtS x nucs epsVec(1) epsVec(2) epsVec(3)
c     where nucs is one of: pp0, nn0, pn+, np-
c     and epsVec components are real (imaginary parts assumed 0)
c
c     Output format:
c       Mout(-1,-1)= real imag
c       Mout(-1, 1)= real imag
c       Mout( 1,-1)= real imag
c       Mout( 1, 1)= real imag

      program testRawM
      implicit none
      include '../common-densities/constants.def'

c     Variables
      double precision sqrtS, sqrtSReal, x, mNucl
      character*3 nucs
      integer MaxEll
      complex*16 epsVec(3)
      complex*16 Mout(-1:1,-1:1)
      double precision eps1, eps2, eps3
      integer i, j

c     Read input
      read(*,*) sqrtS, x, nucs, eps1, eps2, eps3

c     Set epsVec (real values only for simplicity)
      epsVec(1) = dcmplx(eps1, 0.0d0)
      epsVec(2) = dcmplx(eps2, 0.0d0)
      epsVec(3) = dcmplx(eps3, 0.0d0)

c     Set nucleon mass based on reaction type
      if (nucs .eq. 'pp0' .or. nucs .eq. 'pn+') then
        mNucl = Mproton
      else
        mNucl = Mneutron
      endif

c     sqrtSReal = sqrtS for this test
      sqrtSReal = sqrtS

c     MaxEll=5 to match Python's range(6)
      MaxEll = 5

c     Initialize output
      Mout = c0

c     Call getRawM
      call getRawM(sqrtS, x, nucs, Mout, mNucl, sqrtSReal, MaxEll,
     &             epsVec)

c     Output results
      do i = -1, 1, 2
        do j = -1, 1, 2
          write(*,'(A,I2,A,I2,A,E25.16,1X,E25.16)')
     &      'Mout(', i, ',', j, ')=', real(Mout(i,j)), aimag(Mout(i,j))
        end do
      end do

      end program
