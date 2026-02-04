c     Test program for getRawM subroutine
      program testRawM
      implicit none
      include '../common-densities/constants.def'

c     Variables
      double precision sqrtS, sqrtSReal, x, mNucl
      character*3 nucs
      integer MaxEll
      complex*16 epsVec(3)
      complex*16 Mout(-1:1,-1:1)
      integer i, j

c     Set MaxEll (matching Python's getF which uses range(6) -> 0..5)
      MaxEll = 5


c     Test case 1
      sqrtS = 1078.2423100000d0
      sqrtSReal = sqrtS
      x = 1.0000000000d0
      nucs = 'pp0'
      mNucl = 938.2723100000d0
      epsVec(1) = dcmplx(1.0000000000d0, 0.0000000000d0)
      epsVec(2) = dcmplx(0.0000000000d0, 0.0000000000d0)
      epsVec(3) = dcmplx(0.0000000000d0, 0.0000000000d0)

      Mout = c0
      call getRawM(sqrtS, x, nucs, Mout, mNucl, sqrtSReal, MaxEll, epsVec)

      write(*,'(A,I3)') 'TEST_CASE:', 1
      write(*,'(A,F15.10)') 'sqrtS:', sqrtS
      write(*,'(A,F15.10)') 'x:', x
      write(*,'(A,A3)') 'nucs:', nucs
      do i = -1, 1, 2
        do j = -1, 1, 2
          write(*,'(A,I2,A,I2,A,E25.16,A,E25.16)')
     &      'Mout(', i, ',', j, '):', real(Mout(i,j)), ' ', aimag(Mout(i,j))
        end do
      end do
      write(*,'(A)') 'END_TEST_CASE'


c     Test case 2
      sqrtS = 1100.0000000000d0
      sqrtSReal = sqrtS
      x = 1.0000000000d0
      nucs = 'pp0'
      mNucl = 938.2723100000d0
      epsVec(1) = dcmplx(1.0000000000d0, 0.0000000000d0)
      epsVec(2) = dcmplx(0.0000000000d0, 0.0000000000d0)
      epsVec(3) = dcmplx(0.0000000000d0, 0.0000000000d0)

      Mout = c0
      call getRawM(sqrtS, x, nucs, Mout, mNucl, sqrtSReal, MaxEll, epsVec)

      write(*,'(A,I3)') 'TEST_CASE:', 2
      write(*,'(A,F15.10)') 'sqrtS:', sqrtS
      write(*,'(A,F15.10)') 'x:', x
      write(*,'(A,A3)') 'nucs:', nucs
      do i = -1, 1, 2
        do j = -1, 1, 2
          write(*,'(A,I2,A,I2,A,E25.16,A,E25.16)')
     &      'Mout(', i, ',', j, '):', real(Mout(i,j)), ' ', aimag(Mout(i,j))
        end do
      end do
      write(*,'(A)') 'END_TEST_CASE'


c     Test case 3
      sqrtS = 1150.0000000000d0
      sqrtSReal = sqrtS
      x = 1.0000000000d0
      nucs = 'pp0'
      mNucl = 938.2723100000d0
      epsVec(1) = dcmplx(1.0000000000d0, 0.0000000000d0)
      epsVec(2) = dcmplx(0.0000000000d0, 0.0000000000d0)
      epsVec(3) = dcmplx(0.0000000000d0, 0.0000000000d0)

      Mout = c0
      call getRawM(sqrtS, x, nucs, Mout, mNucl, sqrtSReal, MaxEll, epsVec)

      write(*,'(A,I3)') 'TEST_CASE:', 3
      write(*,'(A,F15.10)') 'sqrtS:', sqrtS
      write(*,'(A,F15.10)') 'x:', x
      write(*,'(A,A3)') 'nucs:', nucs
      do i = -1, 1, 2
        do j = -1, 1, 2
          write(*,'(A,I2,A,I2,A,E25.16,A,E25.16)')
     &      'Mout(', i, ',', j, '):', real(Mout(i,j)), ' ', aimag(Mout(i,j))
        end do
      end do
      write(*,'(A)') 'END_TEST_CASE'


c     Test case 4
      sqrtS = 1150.0000000000d0
      sqrtSReal = sqrtS
      x = 0.0000000000d0
      nucs = 'pp0'
      mNucl = 938.2723100000d0
      epsVec(1) = dcmplx(1.0000000000d0, 0.0000000000d0)
      epsVec(2) = dcmplx(0.0000000000d0, 0.0000000000d0)
      epsVec(3) = dcmplx(0.0000000000d0, 0.0000000000d0)

      Mout = c0
      call getRawM(sqrtS, x, nucs, Mout, mNucl, sqrtSReal, MaxEll, epsVec)

      write(*,'(A,I3)') 'TEST_CASE:', 4
      write(*,'(A,F15.10)') 'sqrtS:', sqrtS
      write(*,'(A,F15.10)') 'x:', x
      write(*,'(A,A3)') 'nucs:', nucs
      do i = -1, 1, 2
        do j = -1, 1, 2
          write(*,'(A,I2,A,I2,A,E25.16,A,E25.16)')
     &      'Mout(', i, ',', j, '):', real(Mout(i,j)), ' ', aimag(Mout(i,j))
        end do
      end do
      write(*,'(A)') 'END_TEST_CASE'


c     Test case 5
      sqrtS = 1250.0000000000d0
      sqrtSReal = sqrtS
      x = 1.0000000000d0
      nucs = 'pp0'
      mNucl = 938.2723100000d0
      epsVec(1) = dcmplx(1.0000000000d0, 0.0000000000d0)
      epsVec(2) = dcmplx(0.0000000000d0, 0.0000000000d0)
      epsVec(3) = dcmplx(0.0000000000d0, 0.0000000000d0)

      Mout = c0
      call getRawM(sqrtS, x, nucs, Mout, mNucl, sqrtSReal, MaxEll, epsVec)

      write(*,'(A,I3)') 'TEST_CASE:', 5
      write(*,'(A,F15.10)') 'sqrtS:', sqrtS
      write(*,'(A,F15.10)') 'x:', x
      write(*,'(A,A3)') 'nucs:', nucs
      do i = -1, 1, 2
        do j = -1, 1, 2
          write(*,'(A,I2,A,I2,A,E25.16,A,E25.16)')
     &      'Mout(', i, ',', j, '):', real(Mout(i,j)), ' ', aimag(Mout(i,j))
        end do
      end do
      write(*,'(A)') 'END_TEST_CASE'


c     Test case 6
      sqrtS = 1079.5356300000d0
      sqrtSReal = sqrtS
      x = 1.0000000000d0
      nucs = 'nn0'
      mNucl = 939.5656300000d0
      epsVec(1) = dcmplx(1.0000000000d0, 0.0000000000d0)
      epsVec(2) = dcmplx(0.0000000000d0, 0.0000000000d0)
      epsVec(3) = dcmplx(0.0000000000d0, 0.0000000000d0)

      Mout = c0
      call getRawM(sqrtS, x, nucs, Mout, mNucl, sqrtSReal, MaxEll, epsVec)

      write(*,'(A,I3)') 'TEST_CASE:', 6
      write(*,'(A,F15.10)') 'sqrtS:', sqrtS
      write(*,'(A,F15.10)') 'x:', x
      write(*,'(A,A3)') 'nucs:', nucs
      do i = -1, 1, 2
        do j = -1, 1, 2
          write(*,'(A,I2,A,I2,A,E25.16,A,E25.16)')
     &      'Mout(', i, ',', j, '):', real(Mout(i,j)), ' ', aimag(Mout(i,j))
        end do
      end do
      write(*,'(A)') 'END_TEST_CASE'


c     Test case 7
      sqrtS = 1100.0000000000d0
      sqrtSReal = sqrtS
      x = 1.0000000000d0
      nucs = 'nn0'
      mNucl = 939.5656300000d0
      epsVec(1) = dcmplx(1.0000000000d0, 0.0000000000d0)
      epsVec(2) = dcmplx(0.0000000000d0, 0.0000000000d0)
      epsVec(3) = dcmplx(0.0000000000d0, 0.0000000000d0)

      Mout = c0
      call getRawM(sqrtS, x, nucs, Mout, mNucl, sqrtSReal, MaxEll, epsVec)

      write(*,'(A,I3)') 'TEST_CASE:', 7
      write(*,'(A,F15.10)') 'sqrtS:', sqrtS
      write(*,'(A,F15.10)') 'x:', x
      write(*,'(A,A3)') 'nucs:', nucs
      do i = -1, 1, 2
        do j = -1, 1, 2
          write(*,'(A,I2,A,I2,A,E25.16,A,E25.16)')
     &      'Mout(', i, ',', j, '):', real(Mout(i,j)), ' ', aimag(Mout(i,j))
        end do
      end do
      write(*,'(A)') 'END_TEST_CASE'


c     Test case 8
      sqrtS = 1150.0000000000d0
      sqrtSReal = sqrtS
      x = 0.0000000000d0
      nucs = 'nn0'
      mNucl = 939.5656300000d0
      epsVec(1) = dcmplx(1.0000000000d0, 0.0000000000d0)
      epsVec(2) = dcmplx(0.0000000000d0, 0.0000000000d0)
      epsVec(3) = dcmplx(0.0000000000d0, 0.0000000000d0)

      Mout = c0
      call getRawM(sqrtS, x, nucs, Mout, mNucl, sqrtSReal, MaxEll, epsVec)

      write(*,'(A,I3)') 'TEST_CASE:', 8
      write(*,'(A,F15.10)') 'sqrtS:', sqrtS
      write(*,'(A,F15.10)') 'x:', x
      write(*,'(A,A3)') 'nucs:', nucs
      do i = -1, 1, 2
        do j = -1, 1, 2
          write(*,'(A,I2,A,I2,A,E25.16,A,E25.16)')
     &      'Mout(', i, ',', j, '):', real(Mout(i,j)), ' ', aimag(Mout(i,j))
        end do
      end do
      write(*,'(A)') 'END_TEST_CASE'


      end program
