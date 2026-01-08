cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Isospin Mapping Subroutines for Pion-Pion Scattering
c     Based on IsospinMap.py
c     Alexander Long 2025
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES/FUNCTIONS:
c              PionPionBC            : compute matrix element of isotensor operator
c              combineOpers          : combine two 2x2 operators into 4x4 coupled basis
c              combinePaulis         : convenience function to combine Pauli matrices
c              getCoupledIndex       : map quantum numbers to matrix indices
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Pauli matrices (τ^a) stored as module data
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getPauliMatrix(a, tau)
      implicit none
      integer, intent(in) :: a              ! Which Pauli matrix (1, 2, or 3)
      complex*16, intent(out) :: tau(2,2)   ! Output 2x2 matrix

      tau = dcmplx(0.d0, 0.d0)

      if (a .eq. 1) then
         ! τ^1 = σx
         tau(1,2) = dcmplx(1.d0, 0.d0)
         tau(2,1) = dcmplx(1.d0, 0.d0)
      else if (a .eq. 2) then
         ! τ^2 = σy
         tau(1,2) = dcmplx(0.d0, -1.d0)
         tau(2,1) = dcmplx(0.d0, 1.d0)
      else if (a .eq. 3) then
         ! τ^3 = σz
         tau(1,1) = dcmplx(1.d0, 0.d0)
         tau(2,2) = dcmplx(-1.d0, 0.d0)
      else
         write(*,*) "Error: a must be 1, 2, or 3"
         stop
      end if

      end subroutine getPauliMatrix


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Get the unitary transformation matrix U
c     Transforms between uncoupled (|++>, |+->, |-+>, |-->)
c     and coupled ((1,+1), (1,0), (1,-1), (0,0)) bases
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getUnitaryMatrix(U)
      implicit none
      complex*16, intent(out) :: U(4,4)
      real*8 inv_sqrt2

      inv_sqrt2 = 1.d0 / dsqrt(2.d0)

      U = dcmplx(0.d0, 0.d0)

      ! Row 1: |++>
      U(1,1) = dcmplx(1.d0, 0.d0)

      ! Row 2: |+->
      U(2,2) = dcmplx(inv_sqrt2, 0.d0)
      U(2,4) = dcmplx(0.d0, inv_sqrt2)

      ! Row 3: |-+>
      U(3,2) = dcmplx(inv_sqrt2, 0.d0)
      U(3,4) = dcmplx(0.d0, -inv_sqrt2)

      ! Row 4: |-->
      U(4,3) = dcmplx(1.d0, 0.d0)
      write(*,*) "Tried to get unitary matrix"
      error stop
      end subroutine getUnitaryMatrix


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Map coupled quantum numbers (t, mt), (t', mt') to matrix indices
c     Coupled basis ordering: (1,+1), (1,0), (1,-1), (0,0)
c     Returns indices as 1-based for Fortran
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getCoupledIndex(t, mt, tp, mtp, irow, icol)
      implicit none
      integer, intent(in) :: t, mt, tp, mtp
      integer, intent(out) :: irow, icol

      ! Map (tp, mtp) to row index (bra)
      if (tp .eq. 1 .and. mtp .eq. 1) then
         irow = 1
      else if (tp .eq. 1 .and. mtp .eq. 0) then
         irow = 2
      else if (tp .eq. 1 .and. mtp .eq. -1) then
         irow = 3
      else if (tp .eq. 0 .and. mtp .eq. 0) then
         irow = 4
      else
         write(*,*) "Error: Invalid (tp, mtp) = (", tp, ",", mtp, ")"
         stop
      end if

      ! Map (t, mt) to column index (ket)
      if (t .eq. 1 .and. mt .eq. 1) then
         icol = 1
      else if (t .eq. 1 .and. mt .eq. 0) then
         icol = 2
      else if (t .eq. 1 .and. mt .eq. -1) then
         icol = 3
      else if (t .eq. 0 .and. mt .eq. 0) then
         icol = 4
      else
         write(*,*) "Error: Invalid (t, mt) = (", t, ",", mt, ")"
         stop
      end if

      end subroutine getCoupledIndex


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Kronecker product of two 2x2 matrices resulting in 4x4 matrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine kroneckerProduct(A, B, result)
      implicit none
      complex*16, intent(in) :: A(2,2), B(2,2)
      complex*16, intent(out) :: result(4,4)
      integer i, j, k, l, row, col

      do i = 1, 2
         do j = 1, 2
            do k = 1, 2
               do l = 1, 2
                  row = (i-1)*2 + k
                  col = (j-1)*2 + l
                  result(row, col) = A(i,j) * B(k,l)
               end do
            end do
         end do
      end do

      end subroutine kroneckerProduct


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Combine two single-particle 2x2 operators into coupled basis 4x4 operator
c     Steps: 1) Kronecker product oper1 ⊗ oper2
c            2) Rotate to coupled basis: U† (oper1 ⊗ oper2) U
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine combineOpers(oper1, oper2, total_cpl)
      implicit none
      complex*16, intent(in) :: oper1(2,2), oper2(2,2)
      complex*16, intent(out) :: total_cpl(4,4)
      complex*16 total_unc(4,4), U(4,4), UH(4,4), temp(4,4)
      integer i, j, k

      ! Get Kronecker product
      call kroneckerProduct(oper1, oper2, total_unc)

      ! Get unitary matrix U
c     call getUnitaryMatrix(U)

      ! Compute U† (Hermitian conjugate)
      do i = 1, 4
         do j = 1, 4
            UH(i,j) = conjg(U(j,i))
         end do
      end do

      ! Compute temp = U† * total_unc
      temp = dcmplx(0.d0, 0.d0)
      do i = 1, 4
         do j = 1, 4
            do k = 1, 4
               temp(i,j) = temp(i,j) + UH(i,k) * total_unc(k,j)
            end do
         end do
      end do

      ! Compute total_cpl = temp * U
      total_cpl = temp
c     do i = 1, 4
c        do j = 1, 4
c           do k = 1, 4
c              total_cpl(i,j) = total_cpl(i,j) + temp(i,k) * U(k,j)
c           end do
c        end do
c     end do

      end subroutine combineOpers

      subroutine physicalPionOper(extQnum,oper)
      implicit none
      include '../common-densities/constants.def'
      integer, intent(in) :: extQnum
      complex*16, intent(out) :: oper(2,2)
      integer :: charge
      real*8 invSqrt
      complex*16 :: tau1(2,2), tau2(2,2)


      charge = extQnum-2
      invSqrt=1/sqrt(2.d0)
      if (charge.eq.-1) then
        call getPauliMatrix(1, tau1)
        call getPauliMatrix(2, tau2)
        oper=invSqrt*tau1-ci*invSqrt*tau2

      else if (charge.eq.0) then
        call getPauliMatrix(3, tau1)
        oper=tau1

      else if (charge.eq.1) then
        call getPauliMatrix(1, tau1)
        call getPauliMatrix(2, tau2)
        oper=invSqrt*tau1+ci*invSqrt*tau2

      else
        write(*,*) "Charge passed to subroutine `physicalPions` doesn't work"
        write(*,*) "charge=", charge 
        write(*,*) "extQnum=", extQnum 
        error stop
      end if

      end subroutine physicalPionOper
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Convenience function to combine two Pauli matrices in coupled basis
c     Returns: U† (τ^iso ⊗ τ^isoPrime) U
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function combinePaulis(iso, isoPrime) result(contrib)
      implicit none
      integer, intent(in) :: iso, isoPrime
      complex*16 :: contrib(4,4)
      complex*16 tau(2,2), taup(2,2)

      call physicalPionOper(iso, tau)
      call physicalPionOper(isoPrime, taup)

c     call getPauliMatrix(iso, tau)
c     call getPauliMatrix(isoPrime, taup)
      call combineOpers(tau, taup, contrib)

      end function combinePaulis


      function combinePaulisMathBasis(iso, isoPrime) result(contrib)
      implicit none
      integer, intent(in) :: iso, isoPrime
      complex*16 :: contrib(4,4)
      complex*16 tau(2,2), taup(2,2)


      call getPauliMatrix(iso, tau)
      call getPauliMatrix(isoPrime, taup)
      call combineOpers(tau, taup, contrib)
      end function combinePaulisMathBasis


      subroutine PionPionA(t, mt, tp, mtp, extQnum, result)
c     Computes the PionPion 2 body diagram A isospin contribution
c     Assumes no charge exchange, calculates
c
c     2 (\vec{\tau}_1 \cdot \vec{\tau}_2 - \tau_1^a\tau_2^a)
c
      implicit none
      integer, intent(in) :: t, mt, tp, mtp,extQnum
      complex*16, intent(out) :: result
      complex*16 total_cpl(4,4), contrib(4,4)
      integer i, irow, icol

      interface
         function combinePaulis(iso, isoPrime) result(contrib)
            integer, intent(in) :: iso, isoPrime
            complex*16 :: contrib(4,4)
         end function combinePaulis
      end interface


      ! Initialize total operator
      total_cpl = dcmplx(0.d0, 0.d0)

      ! Add Σ_{i=1}^3 τ_1^i ⊗ τ_2^i
      do i = 1, 3
         total_cpl = total_cpl + 2.d0*combinePaulis(i, i)
      end do

      ! Get matrix element
      call getCoupledIndex(t, mt, tp, mtp, irow, icol)
      result = total_cpl(irow, icol)

      end subroutine PionPionA

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Compute matrix element in coupled basis of isotensor-like structure:
c        Σ_{i=1}^3 (τ_1^i τ_2^i) - 2 (τ_1^a τ_2^a)
c     where a = extQnum
c
c     Returns: <t' mt' | O | t mt>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine PionPionBC(t, mt, tp, mtp, extQnum, result)! ,useMathBasis)
      implicit none
      integer, intent(in) :: t, mt, tp, mtp,extQnum
      complex*16, intent(out) :: result
      complex*16 total_cpl(4,4), contrib(4,4)
      integer i, irow, icol
c     logical, optional :: useMathBasis

      interface
         function combinePaulis(iso, isoPrime) result(contrib)
            integer, intent(in) :: iso, isoPrime
            complex*16 :: contrib(4,4)
         end function combinePaulis

         function combinePaulisMathBasis(iso, isoPrime) result(contrib)
            integer, intent(in) :: iso, isoPrime
            complex*16 :: contrib(4,4)
         end function combinePaulisMathBasis
      end interface

      ! Initialize total operator
      total_cpl = dcmplx(0.d0, 0.d0)
      ! Add Σ_{i=1}^3 τ_1^i ⊗ τ_2^i
      do i = 1, 3
        total_cpl = total_cpl + combinePaulisMathBasis(i, i)
      end do

       contrib = combinePaulis(extQnum, extQnum)
       total_cpl = total_cpl - 2.d0 * contrib
      ! Add -2 τ_1^a ⊗ τ_2^a

      ! Get matrix element
      call getCoupledIndex(t, mt, tp, mtp, irow, icol)
      result = total_cpl(irow, icol)

      end subroutine PionPionBC


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Compute matrix element in coupled basis of:
c        (τ_1 · τ_2 - τ_1^z τ_2^z) = τ_1^1 τ_2^1 + τ_1^2 τ_2^2
c
c     Returns: <t' mt' | O | t mt>
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine piPhotoOper(t, mt, tp, mtp, result)!,useMathBasis)
      implicit none
      integer, intent(in) :: t, mt, tp, mtp
      complex*16, intent(out) :: result
      complex*16 total_cpl(4,4), sub

      integer i, irow, icol
!     logical, optional :: useMathBasis
      logical :: useMath

      interface
         function combinePaulis(iso, isoPrime) result(contrib)
            integer, intent(in) :: iso, isoPrime
            complex*16 :: contrib(4,4)
         end function combinePaulis

         function combinePaulisMathBasis(iso, isoPrime) result(contrib)
            integer, intent(in) :: iso, isoPrime
            complex*16 :: contrib(4,4)
         end function combinePaulisMathBasis
      end interface

      ! Initialize total operator
      total_cpl = dcmplx(0.d0, 0.d0)

      ! Add τ_1^1 ⊗ τ_2^1 + τ_1^2 ⊗ τ_2^2 (x and y components only)
c     This is only for neutral pion scattering.
c     The operator is \delta^{a,b} \vec{\tau_1}\cdot \vec{\tau}_2 -(tau_1^a\tau_2^b + \tau_1^b \tau_2^a)
c     if (present(useMathBasis)) then
c       useMath = useMathBasis
c     else
c       useMath = .False.
c     end if
      useMath=.False.

      total_cpl = combinePaulisMathBasis(1,1)+combinePaulisMathBasis(2,2)+combinePaulisMathBasis(3,3)
      total_cpl=total_cpl-combinePaulis(3,3)
c     total_cpl = combinePaulis(1, 1)+combinePaulis(2, 2)
      ! Get matrix element
      call getCoupledIndex(t, mt, tp, mtp, irow, icol)
      result = total_cpl(irow, icol)

      end subroutine piPhotoOper
