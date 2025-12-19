c     ================================================================
c     Main program for pion scattering calculations  
c     ================================================================
      program pionScattering
      use parseFileData
      use pionScatLib
      implicit none

      integer :: verbosity = 1
      double precision :: sqrtS, sqrtSReal, x, DSG, MatDSG, ratio, mNucl
      integer :: isospin, piCharge, theta

      call initializeFileData('said-pi.txt', verbosity)

      sqrtS = 1162.0d0
      isospin = 1
      piCharge = 1

c     Set nucleon mass based on isospin
      if (isospin .eq. -1) then
         mNucl = Mneutron
      else if (isospin .eq. 1) then
         mNucl = Mproton
      endif

c     For pion-nucleon scattering, sqrtSReal = sqrtS
      sqrtSReal = sqrtS

      write(*,*) '  theta      CC       Mat CS'
      do theta = 0, 180, 30
         x = dcos(theta * 3.141592653589793d0 / 180.0d0)
         call getcsGH(sqrtS, x, isospin, piCharge, DSG, mNucl,
     &                sqrtSReal)
         call getCS(sqrtS, x, isospin, piCharge, MatDSG, mNucl)
         ratio = DSG / MatDSG
         write(*,'(F7.2,2F11.6,F11.6)') dble(theta), DSG, MatDSG, ratio
      enddo
      
      end program pionScattering
