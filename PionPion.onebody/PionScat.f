      module pionScatLib
      use parseFileData
      implicit none
      
c     Physical constants (explicit values to ensure they're not zero)
      double precision, parameter :: mpiPlus = 139.5675d0
      double precision, parameter :: mN = 938.919d0
      double precision, parameter :: MeVtofm = 197.327d0 !197.3 fm= MeV^-1
      double precision, parameter :: MYPI = 3.141592653589794d0
      
c     Pauli matrices (as 2x2 complex arrays)
      double complex, parameter :: sigx(2,2) = reshape((/
     &    (0.0d0,0.0d0), (1.0d0,0.0d0),
     &    (1.0d0,0.0d0), (0.0d0,0.0d0) /), (/2,2/))
      double complex, parameter :: sigy(2,2) = reshape((/
     &    (0.0d0,0.0d0), (0.0d0,-1.0d0),
     &    (0.0d0,1.0d0), (0.0d0,0.0d0) /), (/2,2/))
      double complex, parameter :: sigz(2,2) = reshape((/
     &    (1.0d0,0.0d0), (0.0d0,0.0d0),
     &    (0.0d0,0.0d0), (-1.0d0,0.0d0) /), (/2,2/))
      double complex, parameter :: iden(2,2) = reshape((/
     &    (1.0d0,0.0d0), (0.0d0,0.0d0),
     &    (0.0d0,0.0d0), (1.0d0,0.0d0) /), (/2,2/))

c     Isospin projection matrices for I=1
      double complex, parameter :: P1plus(2,2) = reshape((/
     &    (0.0d0,0.0d0), (0.0d0,0.0d0),
     &    (0.0d0,0.0d0), (0.6666666666666666d0,0.0d0) /), (/2,2/))
      double complex, parameter :: P1minus(2,2) = reshape((/
     &    (0.6666666666666666d0,0.0d0), (0.0d0,0.0d0),
     &    (0.0d0,0.0d0), (0.0d0,0.0d0) /), (/2,2/))
      double complex, parameter :: P1neutral(2,2) = reshape((/
     &    (0.3333333333333333d0,0.0d0), (0.0d0,0.0d0),
     &    (0.0d0,0.0d0), (0.3333333333333333d0,0.0d0) /), (/2,2/))

c     Isospin projection matrices for I=3
      double complex, parameter :: P3plus(2,2) = reshape((/
     &    (1.0d0,0.0d0), (0.0d0,0.0d0),
     &    (0.0d0,0.0d0), (0.3333333333333333d0,0.0d0) /), (/2,2/))
      double complex, parameter :: P3minus(2,2) = reshape((/
     &    (0.3333333333333333d0,0.0d0), (0.0d0,0.0d0),
     &    (0.0d0,0.0d0), (1.0d0,0.0d0) /), (/2,2/))
      double complex, parameter :: P3neutral(2,2) = reshape((/
     &    (0.6666666666666666d0,0.0d0), (0.0d0,0.0d0),
     &    (0.0d0,0.0d0), (0.6666666666666666d0,0.0d0) /), (/2,2/))
      
      contains

c     ================================================================
c     ================================================================
      subroutine getCS(sqrtS, x, isospin, piCharge, CrossSec,
     &                        mNucl)
c     Calculate cross section from scattering matrix
c     ================================================================
      implicit none
      double precision sqrtS, x, CrossSec,mNucl
      integer isospin, piCharge

      double complex mat(-1:1,-1:1), matDag(-1:1,-1:1)
      double precision fudgeFactor, normFactor
      double complex trace
      integer i, j

      fudgeFactor = 2.0d0

      call getMat(sqrtS, x, isospin, piCharge, mat, sqrtS,mNucl)

c     Divide out the 8*pi*sqrtSReal factor that getMat includes
      normFactor = 8.0d0 * MYPI * sqrtS
      mat = mat / normFactor

c     Preform conjugate transpose
      do i = -1, 1, 2
        do j = -1, 1, 2
          matDag(i, j) = dconjg(mat(j, i))
        end do
      end do

c     Calculate trace of mat * matDag
      trace = (0.0d0, 0.0d0)
      do i = -1, 1, 2
         do j = -1, 1, 2
            trace = trace + mat(i,j)*matDag(j,i)
         enddo
      enddo

      CrossSec = 10.0d0 * dreal(trace) / 4.0d0

      CrossSec = CrossSec * fudgeFactor*MeVtofm*MeVtofm
      write(*,*) "CrossSec=", CrossSec 
      return
      end subroutine getCS

c     ================================================================
      subroutine getMat(sqrtS, x, isospin, piCharge, resultmat,sqrtSReal,mNucl)
c     Calculate scattering matrix
c     ================================================================
      implicit none
      double precision sqrtS, x,sqrtSReal
      integer isospin, piCharge
      double complex mat(2,2)
      double complex resultmat(-1:1,-1:1)
      
      double complex g, h
      double precision m1, m2
      double precision S, qVec(3), qpVec(3)
      double precision qHat(3), qpHat(3), cross(3)
      double complex matFactor(2,2)
      double complex sigVec(3,2,2)
      double precision piMassArr(-1:1), mNucl

      integer i, j, k
      
      sigVec(1,:,:)=sigx
      sigVec(2,:,:)=sigy
      sigVec(3,:,:)=sigz
      
      call getGH(sqrtS, x, isospin, piCharge, g, h,sqrtSReal,mNucl) 
c     Get pion and nucleon masses
      piMassArr(-1)=mpi
      piMassArr(0)=mpi0
      piMassArr(1)=mpi
      m1=piMassArr(piCharge)
      m2=mNucl
      call getKinematics(sqrtSReal, x, m1, m2, m1, m2, S, qVec, qpVec)

c     Debug output
c     write(*,*) "DEBUG getMat: isospin=", isospin,
c    &           " piCharge=", piCharge
c     write(*,*) "  sqrtS=", sqrtS, " sqrtSReal=", sqrtSReal
c     write(*,*) "  qVec=", qVec(1), qVec(2), qVec(3)
c     write(*,*) "  qpVec=", qpVec(1), qpVec(2), qpVec(3)
c     write(*,*) "  |qVec|=", vecAbs(qVec), " |qpVec|=", vecAbs(qpVec)

c     Normalize vectors
      call normalize(qVec, qHat)
      call normalize(qpVec, qpHat)
c     Calculate cross product qpHat x qHat
      cross(1) = qpHat(2)*qHat(3) - qpHat(3)*qHat(2)
      cross(2) = qpHat(3)*qHat(1) - qpHat(1)*qHat(3)
      cross(3) = qpHat(1)*qHat(2) - qpHat(2)*qHat(1)
c     Calculate matFactor = -i * (sig · (qpHat x qHat))
      matFactor=0.d0
      matFactor=sigVec(1,:,:)*cross(1)+sigVec(2,:,:)*cross(2)+sigVec(3,:,:)*cross(3)
      matFactor=matFactor*cmplx(0,-1,kind=8)
      
c     Calculate final matrix: mat = iden * g + matFactor * h
      mat=iden*g+matFactor*h
      resultmat(-1,-1)=mat(1,1)
      resultmat(-1, 1)=mat(1,2)
      resultmat( 1,-1)=mat(2,1)
      resultmat( 1, 1)=mat(2,2)! resultMat has units of MeV^-1 here
      !factor 1/sqrt(2) could come from spin averaging?
      resultmat= resultmat*8*pi*sqrtSReal!/sqrt(2.d0)!unitless,sqrtS real has units of MeV
      return
      end subroutine getMat

c     ================================================================
      subroutine getcsGH(sqrtS, x, isospin, piCharge, DSG, mNucl,sqrtSReal)
c     Calculate cross section using g,h amplitudes
c     ================================================================
      implicit none
      double precision sqrtS, x, DSG, mNucl, sqrtSReal
      integer isospin, piCharge

      double complex g, h
      double precision sintheta

      call getGH(sqrtS, x, isospin, piCharge, g, h, sqrtSReal, mNucl)
      sintheta = dsqrt(1.0d0 - x*x)
      
      DSG = cdabs(g)**2 + cdabs(h * sintheta)**2
      DSG = DSG * 10.0d0*(8*pi*sqrtSReal)**2
      
      return
      end subroutine getcsGH

c     ================================================================
      subroutine getGH(sqrtS, x, isospin, piCharge, g, h,sqrtSReal,mNucl)
c     Calculate g and h scattering amplitudes in fm
c     ================================================================
      implicit none
      double precision sqrtS, x, sqrtSReal, mNucl
      integer isospin, piCharge
      double complex g, h
      
      integer twoIs(2), ells(5)
      integer twoI, ell, chargeIdx
      double precision m1, m2, S, qVec(3), qpVec(3)
      double complex gTerm, hTerm, gTermTmp, hTermTmp
      double complex fPlus, fMinus
      double precision poly0, poly1, weight
      double complex projII(2,2)
      double precision isoVec(2)
      integer i, j
      
      data twoIs /1, 3/! `data` keyword sets these values at compile time
      data ells /0, 1, 2, 3, 4/
      
c     Set isospin vector
      if (isospin .eq. 1) then
         isoVec(1) = 1.0d0
         isoVec(2) = 0.0d0
      else if (isospin .eq. -1) then
         isoVec(1) = 0.0d0
         isoVec(2) = 1.0d0
      endif
      
c     Get masses
      if (piCharge .eq. -1) then
         m1 = mpiPlus
      else if (piCharge .eq. 0) then
         m1 = mpi0
      else if (piCharge .eq. 1) then
         m1 = mpiPlus
      endif

c     if (isospin .eq. -1) then
c        m2 = Mneutron
c     else if (isospin .eq. 1) then
c        m2 = Mproton
c     endif
      m2=mNucl
      call getKinematics(sqrtSReal, x, m1, m2, m1, m2, S, qVec, qpVec)

      gTerm = (0.0d0, 0.0d0)
      hTerm = (0.0d0, 0.0d0)
      
c     Loop over isospin values I=1/2, 3/2
      do i = 1, 2
         twoI = twoIs(i)
         gTermTmp = (0.0d0, 0.0d0)
         hTermTmp = (0.0d0, 0.0d0)
         
c        Get charge index (0→π⁺, 1→π⁻, 2→π⁰)
         if (piCharge .eq. 1) then
            chargeIdx = 0
         else if (piCharge .eq. -1) then
            chargeIdx = 1
         else if (piCharge .eq. 0) then
            chargeIdx = 2
         endif
         
c        Get isospin projection matrix
         call getIsospinProjection(twoI, chargeIdx, projII)
         
c        Calculate weight = isoVec · projII · isoVec
         weight = 0.0d0
         do j = 1, 2
            weight = weight + isoVec(j) * dreal(projII(j,1)) * isoVec(1)
            weight = weight + isoVec(j) * dreal(projII(j,2)) * isoVec(2)
         enddo
         
c        Loop over partial waves ℓ = 0…4
         do j = 1, 5
            ell = ells(j)
            
            call getF(qVec, twoI, ell, 1, sqrtS, fPlus)
            call getF(qVec, twoI, ell, -1, sqrtS, fMinus)
            
            call legendre_poly(x, ell, 0, poly0)
            call legendre_poly(x, ell, 1, poly1)
            
            gTermTmp = gTermTmp + ((ell + 1) * fPlus + ell * fMinus) 
     &               * poly0
            hTermTmp = hTermTmp + (fPlus - fMinus) * poly1
         enddo
         
         gTerm = gTerm + gTermTmp * weight
         hTerm = hTerm + hTermTmp * weight
      enddo
      
      g = gTerm
      h = hTerm
      
      return
      end subroutine getGH

c     ================================================================
      subroutine getF(qVec, twoI, ell, sign, sqrtS, fOut)
c     Calculate partial wave amplitude with linear interpolation
c     ================================================================
      implicit none
      double precision qVec(3), sqrtS
      integer twoI, ell, sign
      double complex fOut

      character*1 letter
      integer x1, x2, target2L
      double precision diff, slope
      double complex y1, y2

c     Map ell to letter
      if (ell .eq. 0) then
         letter = 'S'
      else if (ell .eq. 1) then
         letter = 'P'
      else if (ell .eq. 2) then
         letter = 'D'
      else if (ell .eq. 3) then
         letter = 'F'
      else if (ell .eq. 4) then
         letter = 'G'
      else
         fOut = (0.0d0, 0.0d0)
         return
      endif

c     Linear interpolation between even sqrtS values
c     Calculate diff = sqrtS mod 2 (remainder when divided by 2)
      diff = dmod(sqrtS, 2.0d0)

c     x1 is the lower even value, x2 is the upper even value
      x1 = int(sqrtS - diff)
      x2 = x1 + 2

c     Calculate target2L = 2*ell + sign
      target2L = 2 * ell + sign

c     Get values at x1 and x2
      call getFAtValue(qVec, twoI, ell, sign, target2L, letter,
     &                 dble(x1), y1)
      call getFAtValue(qVec, twoI, ell, sign, target2L, letter,
     &                 dble(x2), y2)

c     Linear interpolation: fOut = y1 + slope * diff
c     where slope = (y2 - y1) / (x2 - x1) = (y2 - y1) / 2
      slope = (y1-y2)/(x1-x2)
      fOut = slope*diff+y1

      return
      end subroutine getF

c     ================================================================
      subroutine getFAtValue(qVec, twoI, ell, sign, target2L, letter,
     &                       sqrtS, fOut)
c     Calculate partial wave amplitude at a specific sqrtS value, returns fOut in outs of MeV^-1
c     ================================================================
      implicit none
      double precision qVec(3), sqrtS
      integer twoI, ell, sign, target2L
      character*1 letter
      double complex fOut

      double precision del_result, sr_result, eta, qAbs
      logical found

c     Get scattering data at the specified sqrtS
      call getScatteringData(getLetterIndex(letter), twoI, target2L,
     &     sqrtS, del_result, sr_result, found)

      if (.not. found) then
         del_result = 0.0d0
         sr_result = 0.0d0
      endif

      eta = dsqrt(1.0d0 - sr_result)
      qAbs = vecAbs(qVec)

c     Calculate fOut = (eta * exp(2i*deltaRe) - 1) / (2i*qAbs) * MeVtofm
      fOut = eta * cdexp((0.0d0, 2.0d0) * del_result) - (1.0d0, 0.0d0)
      fOut = fOut / ((0.0d0, 2.0d0) * qAbs)
c     fOut = fOut * MeVtofm

      return
      end subroutine getFAtValue

c     ================================================================
      subroutine getKinematics(sqrtS, x, m1, m2, m3, m4, S, qVec, qpVec)
c     Calculate kinematic variables
c     ================================================================
      implicit none
      double precision sqrtS, x, m1, m2, m3, m4, S
      double precision qVec(3), qpVec(3)
      
      double precision E1cm, E2cm, E3cm, E4cm, absp, absQp
      
      S = sqrtS * sqrtS
      
      E1cm = (S + m1*m1 - m2*m2) / (2.0d0 * sqrtS)
      E2cm = (S - m1*m1 + m2*m2) / (2.0d0 * sqrtS)
      E3cm = (S + m3*m3 - m4*m4) / (2.0d0 * sqrtS)
      E4cm = (S - m3*m3 + m4*m4) / (2.0d0 * sqrtS)

c     Handle near-threshold case: if E is within 5 MeV of m, set to m (momentum = 0)
c     If E < m - 5 MeV, we're significantly below threshold - this is an error
      if (E1cm .lt. m1) then
        if (abs(E1cm - m1) .le. 5.0d0) then
          absp = 0.0d0
        else
          write(*,*) 'ERROR in getKinematics: E1cm significantly below m1'
          write(*,*) 'E1cm =', E1cm, ' MeV'
          write(*,*) 'm1 =', m1, ' MeV'
          write(*,*) 'E1cm - m1 =', E1cm - m1, ' MeV'
          write(*,*) 'sqrtS =', sqrtS, ' MeV'
          stop
        end if
      else
        absp = dsqrt(E1cm*E1cm - m1*m1)
      end if
      if (absp.le.1) then
        absp=1.d0
      end if
      qVec(1) = 0.0d0
      qVec(2) = 0.0d0
      qVec(3) = absp

c     Handle near-threshold case for final state
      if (E3cm .lt. m3) then
        if (abs(E3cm - m3) .le. 5.0d0) then
          absQp = 0.0d0
        else
          write(*,*) 'ERROR in getKinematics: E3cm significantly below m3'
          write(*,*) 'E3cm =', E3cm, ' MeV'
          write(*,*) 'm3 =', m3, ' MeV'
          write(*,*) 'E3cm - m3 =', E3cm - m3, ' MeV'
          write(*,*) 'sqrtS =', sqrtS, ' MeV'
          stop
        end if
      else
        absQp = dsqrt(E3cm*E3cm - m3*m3)
      end if

      if (absQp.le.1) then
        absQp=1.d0
      end if
      qpVec(1) = 0.0d0
      qpVec(2) = dsqrt(1.0d0 - x*x) * absQp
      qpVec(3) = x * absQp
      
      return
      end subroutine getKinematics

c     ================================================================
      subroutine normalize(vec, vecNorm)
c     Normalize a 3D vector
c     ================================================================
      implicit none
      double precision vec(3), vecNorm(3)
      double precision norm
      
      norm = dsqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
      vecNorm(1) = vec(1) / norm
      vecNorm(2) = vec(2) / norm
      vecNorm(3) = vec(3) / norm
      
      return
      end subroutine normalize

c     ================================================================
      subroutine crossProduct(a, b, result)
c     Calculate cross product of two 3D vectors
c     ================================================================
      implicit none
      double precision a(3), b(3), result(3)
      
      result(1) = a(2)*b(3) - a(3)*b(2)
      result(2) = a(3)*b(1) - a(1)*b(3)
      result(3) = a(1)*b(2) - a(2)*b(1)
      
      return
      end subroutine crossProduct

c     ================================================================
      subroutine getIsospinProjection(twoI, chargeIdx, projII)
c     Get isospin projection matrix
c     ================================================================
      implicit none
      integer twoI, chargeIdx
      double complex projII(2,2)
      
      if (twoI .eq. 1) then
         if (chargeIdx .eq. 0) then
            projII = P1plus
         else if (chargeIdx .eq. 1) then
            projII = P1minus
         else if (chargeIdx .eq. 2) then
            projII = P1neutral
         endif
      else if (twoI .eq. 3) then
         if (chargeIdx .eq. 0) then
            projII = P3plus
         else if (chargeIdx .eq. 1) then
            projII = P3minus
         else if (chargeIdx .eq. 2) then
            projII = P3neutral
         endif
      endif
      
      return
      end subroutine getIsospinProjection

c     ================================================================
      function vecAbs(v)
      implicit none
      double precision vecAbs, v(3)
      vecAbs = dsqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
      return
      end

c     ================================================================
      subroutine dag(m, mdag, lower, upper)
c     Calculate conjugate transpose (dagger) of a matrix
c     ================================================================
      implicit none
      double complex m(2, 2), mdag(2, 2)
      integer i, j
      integer lower, upper
      
      do i = lower, upper
        do j = lower, upper
          mdag(i, j) = dconjg(m(j, i))
        end do
      end do
      
      return
      end subroutine dag

c     ================================================================
      function getLetterIndex(letter)
      implicit none
      character*1 letter
      integer getLetterIndex
      
      select case (letter)
         case ('S', 's')
            getLetterIndex = 0
         case ('P', 'p')
            getLetterIndex = 1
         case ('D', 'd')
            getLetterIndex = 2
         case ('F', 'f')
            getLetterIndex = 3
         case ('G', 'g')
            getLetterIndex = 4
         case default
            getLetterIndex = -1
      end select
      return
      end function getLetterIndex

c     ================================================================
      subroutine legendre_poly(x, n, deriv, result)
c     Calculate Legendre polynomial P_n(x) or its derivatives
c     ================================================================
      implicit none
      double precision x, result
      integer n, deriv
      double precision tmp, x2, x3, x4
      
      if (deriv .eq. 0) then
c       Calculate P_n(x)
        if (n .eq. -1 .or. n .eq. 0) then
          result = 1.0d0
        else if (n .eq. 1) then
          result = x
        else if (n .eq. 2) then
          result = 1.5d0 * (x*x) - 0.5d0
        else if (n .eq. 3) then
          result = 2.5d0 * (x*x*x) - 1.5d0 * x
        else if (n .eq. 4) then
          x2 = x*x
          x4 = x2*x2
          tmp = 3.0d0 - 30.0d0 * x2 + 35.0d0 * x4
          result = (1.0d0 / 8.0d0) * tmp
        else if (n .eq. 5) then
          result=(1.d0 / 8.d0) * (15.d0 * x - 70.d0 * x**(3.d0)
     &    + 63.d0 * x**5.d0)
        else
          write(*,*) 'Error: legendreP not implemented for n=', n
          stop
        end if
        
      else if (deriv .eq. 1) then
c       Calculate first derivative dP_n(x)/dx
        if (n .eq. -1 .or. n .eq. 0) then
          result = 0.0d0
        else if (n .eq. 1) then
          result = 1.0d0
        else if (n .eq. 2) then
          result = 3.0d0 * x
        else if (n .eq. 3) then
          x2 = x*x
          tmp = 15.0d0 * x2 - 3.0d0
          result = 0.5d0 * tmp
        else if (n .eq. 4) then
          x2 = x*x
          x3 = x2*x
          tmp = -60.0d0 * x + 140.0d0 * x3
          result = tmp / 8.0d0
        else if (n .eq. 5) then
          result= (15.d0 - 210.d0 * x**2.d0 + 315.d0 * x**4.d0) / 8.d0
        else
          write(*,*) 'Error: legendreP not implemented for n=', n
          stop
        end if
      else
        write(*,*) 'Error: Invalid derivative order in legP'
        stop
      end if
      
      return
      end subroutine legendre_poly

      end module pionScatLib
