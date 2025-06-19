      module pionScatLib
      use parseFileData
      implicit none
      
c     Physical constants (explicit values to ensure they're not zero)
      double precision, parameter :: mpiPlus = 139.5675d0
      double precision, parameter :: mN = 938.919d0
      double precision, parameter :: MeVtofm = 197.327d0
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
      subroutine getCS(sqrtS, x, isospin, piCharge, CrossSec)
c     Calculate cross section from scattering matrix
c     ================================================================
      implicit none
      double precision sqrtS, x, CrossSec
      integer isospin, piCharge
      
      double complex mat(2,2), matDag(2,2)
      double precision fudgeFactor
      double complex trace
      integer i
      
      call getMat(sqrtS, x, isospin, piCharge, mat, sqrtS)
      call dag(mat, matDag, 1, 2)
      
c     Calculate trace of mat * matDag
      trace = (0.0d0, 0.0d0)
      do i = 1, 2
         trace = trace + mat(i,1)*matDag(1,i) + mat(i,2)*matDag(2,i)
      enddo
      
      CrossSec = 10.0d0 * dreal(trace) / 4.0d0
      
c     TODO: find where the missing factor of 2 comes from
      fudgeFactor = 2.0d0
      CrossSec = CrossSec * fudgeFactor
      
      return
      end subroutine getCS

c     ================================================================
      subroutine getMat(sqrtS, x, isospin, piCharge, mat,sqrtSReal)
c     Calculate scattering matrix
c     ================================================================
      implicit none
      double precision sqrtS, x,sqrtSReal
      integer isospin, piCharge
      double complex mat(2,2)
      
      double complex g, h
      double precision m1, m2
      double precision S, qVec(3), qpVec(3)
      double precision qHat(3), qpHat(3), cross(3)
      double complex matFactor(2,2)
      double complex sigVec(3,2,2)
      double precision piMassArr(-1:1)
      double precision nucMassArr(2)
      integer i, j, k
      
c     Initialize sigma vector
      do i = 1, 2
         do j = 1, 2
            sigVec(1,i,j) = sigx(i,j)
            sigVec(2,i,j) = sigy(i,j)
            sigVec(3,i,j) = sigz(i,j)
         enddo
      enddo
      
      call getGH(sqrtS, x, isospin, piCharge, g, h)
      
c     write(*,*) "g=", g 
c     write(*,*) "h=", h 
c     Get pion and nucleon masses
      piMassArr(-1)=mpi
      piMassArr(0)=mpi0
      piMassArr(1)=mpi
      m1=piMassArr(piCharge)
c     nucMassArr(1)=Mneutron
c     nucMassArr(2)=Mproton
c     TODO: impliment this nucMassArray
c     if (piCharge .eq. -1) then
c        m1 = mpiPlus
c     else if (piCharge .eq. 0) then
c        m1 = mpi0
c     else if (piCharge .eq. 1) then
c        m1 = mpiPlus
c     endif
      
      if (isospin .eq. -1) then
         m2 = 939.56563d0
      else if (isospin .eq. 1) then
         m2 = 938.27231d0
      endif
      
      call getKinematics(sqrtSReal, x, m1, m2, m1, m2, S, qVec, qpVec)
      
c     Normalize vectors
      call normalize(qVec, qHat)
      call normalize(qpVec, qpHat)
      
c     Calculate cross product qpHat x qHat
      cross(1) = qpHat(2)*qHat(3) - qpHat(3)*qHat(2)
      cross(2) = qpHat(3)*qHat(1) - qpHat(1)*qHat(3)
      cross(3) = qpHat(1)*qHat(2) - qpHat(2)*qHat(1)
      
c     Calculate matFactor = -i * (sig · (qpHat x qHat))
      do i = 1, 2
        do j = 1, 2
          matFactor(i,j) = (0.0d0, 0.0d0)
          do k = 1, 3
            matFactor(i,j) = matFactor(i,j) + sigVec(k,i,j) * cross(k)
          enddo
          matFactor(i,j) = matFactor(i,j) * (0.0d0, -1.0d0)
        enddo
      enddo
      
c     Calculate final matrix: mat = iden * g + matFactor * h
      do i = 1, 2
         do j = 1, 2
            mat(i,j) = iden(i,j) * g + matFactor(i,j) * h
         enddo
      enddo
      
      return
      end subroutine getMat

c     ================================================================
      subroutine getcsGH(sqrtS, x, isospin, piCharge, DSG)
c     Calculate cross section using g,h amplitudes
c     ================================================================
      implicit none
      double precision sqrtS, x, DSG
      integer isospin, piCharge
      
      double complex g, h
      double precision sintheta
      
      call getGH(sqrtS, x, isospin, piCharge, g, h)
      sintheta = dsqrt(1.0d0 - x*x)
      
      DSG = cdabs(g)**2 + cdabs(h * sintheta)**2
      DSG = DSG * 10.0d0
      
      return
      end subroutine getcsGH

c     ================================================================
      subroutine getGH(sqrtS, x, isospin, piCharge, g, h)
c     Calculate g and h scattering amplitudes
c     ================================================================
      implicit none
      double precision sqrtS, x
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
      
      data twoIs /1, 3/
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
         m1 = mpi
      else if (piCharge .eq. 1) then
         m1 = mpiPlus
      endif
      
      if (isospin .eq. -1) then
         m2 = Mneutron
      else if (isospin .eq. 1) then
         m2 = Mproton
      endif
      
      call getKinematics(sqrtS, x, m1, m2, m1, m2, S, qVec, qpVec)
      
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
c     Calculate partial wave amplitude
c     ================================================================
      implicit none
      double precision qVec(3), sqrtS
      integer twoI, ell, sign
      double complex fOut
      
      character*1 letter
      integer sqrtStmp, target2L
      double precision del_result, sr_result, eta, qAbs
      logical found
      
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
      
c     Round sqrtS to nearest even integer
      sqrtStmp = int(sqrtS / 2.0d0 + 0.5d0)
      sqrtStmp = 2 * sqrtStmp
      
c     Calculate target2L = 2*ell + sign
      target2L = 2 * ell + sign
      
c     Get scattering data  
      call getScatteringData(getLetterIndex(letter), twoI, target2L,
     &     dble(sqrtStmp), del_result, sr_result, found)
      
      if (.not. found) then
         del_result = 0.0d0
         sr_result = 0.0d0
      endif
      
      eta = dsqrt(1.0d0 - sr_result)
      qAbs = vecAbs(qVec)
      
c     Calculate fOut = (eta * exp(2i*deltaRe) - 1) / (2i*qAbs) * MeVtofm
      fOut = eta * cdexp((0.0d0, 2.0d0) * del_result) - (1.0d0, 0.0d0)
      fOut = fOut / ((0.0d0, 2.0d0) * qAbs)
      fOut = fOut * MeVtofm
      
      return
      end subroutine getF

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
      
      absp = dsqrt(E1cm*E1cm - m1*m1)
      qVec(1) = 0.0d0
      qVec(2) = 0.0d0
      qVec(3) = absp
      
      absQp = dsqrt(E3cm*E3cm - m3*m3)
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
