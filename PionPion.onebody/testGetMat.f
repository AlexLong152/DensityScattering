c     ================================================================
c     Test program to call getMat and output results
c     ================================================================
      program testGetMat
      use pionScatLib
      implicit none
      double precision sqrtS, sqrtSReal, x, theta
      integer isospin, piCharge
      double complex resultmat(-1:1,-1:1)

c     Initialize the SAID data file
      call initializeFileData('said-pi.txt', 0)

c     Read input parameters from command line arguments would be complex
c     in Fortran, so we'll read from stdin
      read(*,*) sqrtS, theta, isospin, piCharge

c     Convert theta to cos(theta)
      x = cos(theta * 3.14159265358979323846d0 / 180.0d0)

c     Set sqrtSReal = sqrtS as requested
      sqrtSReal = sqrtS

c     Call getMat
      call getMat(sqrtS, x, isospin, piCharge, resultmat, sqrtSReal)

c     Output results in a format easy to parse from Python
c     Format: real(mat[-1,-1]) imag(mat[-1,-1]) real(mat[-1,1]) ...
      write(*,'(A,E23.15,1X,E23.15)') 'mat(-1,-1)=', resultmat(-1,-1)
      write(*,'(A,E23.15,1X,E23.15)') 'mat(-1, 1)=', resultmat(-1, 1)
      write(*,'(A,E23.15,1X,E23.15)') 'mat( 1,-1)=', resultmat( 1,-1)
      write(*,'(A,E23.15,1X,E23.15)') 'mat( 1, 1)=', resultmat( 1, 1)

      end program testGetMat
