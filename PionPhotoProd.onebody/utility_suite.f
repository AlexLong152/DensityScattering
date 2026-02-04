c ================================
c Utility Functions for Physics Calculations
c ================================
      
c     ==================================================
      double precision function vecAbs(v)
c     Calculate the absolute value of a vector
c     ==================================================
      implicit none
      double precision v(3)
      vecAbs = dsqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
      return
      end
      
c     ==================================================
      double precision function mandT(v1, v2)
c     t=(v1-v2)^2 = v1^2 +v2^2 -2 v1 dot v2
c     ==================================================
      implicit none
      double precision v1(4), v2(4)
      double precision t1, t2, t3
      double precision fourVecDot
      external fourVecDot
      
      t1 = fourVecDot(v1, v1)
      t2 = fourVecDot(v2, v2)
      t3 = fourVecDot(v1, v2)
      
      mandT = t1 + t2 - 2.0d0 * t3
      return
      end
      
c     ==================================================
      double precision function fourVecSquare(v1)
c     Calculate square of a 4-vector (v1)^2
c     ==================================================
      implicit none
      double precision v1(4)
      double precision fourVecDot
      external fourVecDot
      
      fourVecSquare = fourVecDot(v1, v1)
      return
      end
      
c     ==================================================
      double precision function fourVecDot(v1, v2)
c     Calculate dot product of two 4-vectors with Minkowski metric
c     ==================================================
      implicit none
      double precision v1(4), v2(4)
      
      fourVecDot = v1(1)*v2(1) - v1(2)*v2(2) - v1(3)*v2(3) - v1(4)*v2(4)
      return
      end
      
c     ==================================================
      subroutine matDotVec(matVec, vec, result)
c     Multiply a matrix vector by a vector
c     matVec is a 3x2x2 array (3 2x2 matrices)
c     vec is a 3-component vector
c     result is a 2x2 matrix
c     ==================================================
      implicit none
      double complex matVec(3, 2, 2), result(2, 2)
      double precision vec(3)
      integer i, j, k
      
      do i = 1, 2
        do j = 1, 2
          result(i, j) = (0.0d0, 0.0d0)
          do k = 1, 3
            result(i, j) = result(i, j) + matVec(k, i, j) * vec(k)
          end do
        end do
      end do
      return
      end
      
c     ==================================================
      subroutine dag(m, mdag, lower,upper)
c     Calculate conjugate transpose (dagger) of a matrix
c     ==================================================
      implicit none
      double complex m(2, 2), mdag(2, 2)
      integer i, j
      integer lower,upper
      
      do i = lower, upper
        do j = upper, lower
          mdag(i, j) = dconjg(m(j, i))
        end do
      end do
      
      return
      end
      
c     ==================================================
      subroutine crossProduct(a, b, result)
c     Calculate cross product of two 3D vectors A×B
c     For complex vector B, result is also complex
c     ==================================================
      implicit none
c     Input parameters
      double precision a(3)
      double complex b(3)
      double complex result(3)
      
c     Calculate cross product components
      result(1) = a(2)*b(3) - a(3)*b(2)
      result(2) = a(3)*b(1) - a(1)*b(3)
      result(3) = a(1)*b(2) - a(2)*b(1)
      
      return
      end
      
c     ==================================================
      double precision function legP(x, n, deriv)
c     Calculate Legendre polynomial P_n(x) or its derivatives
c     ==================================================
      implicit none
      double precision x
      integer n, deriv
      double precision tmp, x2, x3, x4
      
      if (deriv .eq. 0) then
c       Calculate P_n(x)
        if (n .eq. -1 .or. n .eq. 0) then
          legP = 1.0d0
        else if (n .eq. 1) then
          legP = x
        else if (n .eq. 2) then
          legP = 1.5d0 * (x*x) - 0.5d0
        else if (n .eq. 3) then
          legP = 2.5d0 * (x*x*x) - 1.5d0 * x
        else if (n .eq. 4) then
          x2 = x*x
          x4 = x2*x2
          tmp = 3.0d0 - 30.0d0 * x2 + 35.0d0 * x4
          legP = (1.0d0 / 8.0d0) * tmp
        else if (n .eq. 5) then
          legP=(1.d0 / 8.d0) * (15.d0 * x - 70.d0 * x**(3.d0)
     &    + 63.d0 * x**5.d0)
        else if (n .eq. 6) then
          legP = (1.d0/16.d0) * (231.d0*x**6 - 315.d0*x**4
     &         + 105.d0*x**2 - 5.d0)
        else
          write(*,*) 'Error: legendreP not implemented for n=', n
          stop
        end if

      else if (deriv .eq. 1) then
c       Calculate first derivative dP_n(x)/dx
        if (n .eq. -1 .or. n .eq. 0) then
          legP = 0.0d0
        else if (n .eq. 1) then
          legP = 1.0d0
        else if (n .eq. 2) then
          legP = 3.0d0 * x
        else if (n .eq. 3) then
          x2 = x*x
          tmp = 15.0d0 * x2 - 3.0d0
          legP = 0.5d0 * tmp
        else if (n .eq. 4) then
          x2 = x*x
          x3 = x2*x
          tmp = -60.0d0 * x + 140.0d0 * x3
          legP = tmp / 8.0d0
        else if (n .eq. 5) then
          legP= (15.d0 - 210.d0 * x**2.d0 + 315.d0 * x**4.d0) / 8.d0
        else if (n .eq. 6) then
          legP = (693.d0*x**5)/8.d0 - (315.d0*x**3)/4.d0
     &         + (105.d0*x)/8.d0
        else
          write(*,*) 'Error: legendreP not implemented for n=', n
          stop
        end if

      else if (deriv .eq. 2) then
c       Calculate second derivative d^2P_n(x)/dx^2
        if (n .eq. -1 .or. n .eq. 0 .or. n .eq. 1) then
          legP = 0.0d0
        else if (n .eq. 2) then
          legP = 3.0d0
        else if (n .eq. 3) then
          legP = 15.0d0 * x
        else if (n .eq. 4) then
          x2 = x*x
          tmp = -7.5d0 + 52.5d0 * x2
          legP = tmp
        else if (n .eq. 5) then
          legP= (-420.d0 * x + 1260.d0 * x**3.d0) / 8.d0
        else if (n .eq. 6) then
          legP = (3465.d0*x**4)/8.d0 - (945.d0*x**2)/4.d0
     &         + 105.d0/8.d0
        else
          write(*,*) 'Error: legendreP not implemented for n=', n
          stop
        end if

      else
        write(*,*) 'Error: Invalid derivative order in legP'
        stop
      end if
      
      return
      end
