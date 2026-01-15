c     =================================================================
c     Isospin mapping subroutines for Pion-Pion scattering
c     Based on tests/IsospinMap.py
c     This could in theory be used for combination of two spin 1/2 matricies in general
c     =================================================================

c     ----------------------------------------------------------------
c     FUNCTION: mapping
c     Maps coupled quantum numbers (t, mt), (t', mt') to matrix indices
c     Coupled basis order: (1,+1), (1,0), (1,-1), (0,0)
c     Returns indices (1-based for Fortran)
c     ----------------------------------------------------------------
      subroutine mapping(t, mt, tp, mtp, irow, icol)
      implicit none
      integer t, mt, tp, mtp, irow, icol

c     Map (tp, mtp) to row index
      if (tp .eq. 1) then
         if (mtp .eq. 1) then
            irow = 1
         else if (mtp .eq. 0) then
            irow = 2
         else if (mtp .eq. -1) then
            irow = 3
         endif
      else if (tp .eq. 0 .and. mtp .eq. 0) then
         irow = 4
      endif

c     Map (t, mt) to column index
      if (t .eq. 1) then
         if (mt .eq. 1) then
            icol = 1
         else if (mt .eq. 0) then
            icol = 2
         else if (mt .eq. -1) then
            icol = 3
         endif
      else if (t .eq. 0 .and. mt .eq. 0) then
         icol = 4
      endif

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: getPauliMatrix
c     Returns Pauli matrix tau^i (i=1,2,3) as 2x2 complex matrix
c     ----------------------------------------------------------------
      subroutine getPauliMatrix(i, tau)
      implicit none
      integer i
      complex*16 tau(2,2)
      complex*16 im
      parameter (im = (0.d0, 1.d0))

      if (i .eq. 1) then
c        tau^1 = [[0, 1], [1, 0]]
c        sigmax
         tau(1,1) = (0.d0, 0.d0)
         tau(1,2) = (1.d0, 0.d0)
         tau(2,1) = (1.d0, 0.d0)
         tau(2,2) = (0.d0, 0.d0)
      else if (i .eq. 2) then
c        tau^2 = [[0, -i], [i, 0]]
c        sigmay
         tau(1,1) = (0.d0, 0.d0)
         tau(1,2) = -im 
         tau(2,1) = im
         tau(2,2) = (0.d0, 0.d0)
      else if (i .eq. 3) then
c        tau^3 = [[1, 0], [0, -1]]
c        sigmaz
         tau(1,1) = (1.d0, 0.d0)
         tau(1,2) = (0.d0, 0.d0)
         tau(2,1) = (0.d0, 0.d0)
         tau(2,2) = (-1.d0, 0.d0)
      endif

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: kronecker
c     Compute Kronecker product of two 2x2 matrices -> 4x4 matrix
c     ----------------------------------------------------------------
      subroutine kronecker(A, B, C)
      implicit none
      complex*16 A(2,2), B(2,2), C(4,4)
      integer i, j, k, l, row, col

      do i = 1, 2
         do j = 1, 2
            do k = 1, 2
               do l = 1, 2
                  row = (i-1)*2 + k
                  col = (j-1)*2 + l
                  C(row, col) = A(i,j) * B(k,l)
               enddo
            enddo
         enddo
      enddo

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: combineOpers
c     Combine two spin 1/2 operators into coupled basis
c     Returns 4x4 matrix in coupled basis: |1,1>, |1,0>, |1,-1>, |0,0>
c     ----------------------------------------------------------------
      subroutine combineOpers(oper1, oper2, result)
      implicit none
      complex*16 oper1(2,2), oper2(2,2), result(4,4)
      complex*16 O_prod(4,4), U(4,4), Udag(4,4)
      complex*16 temp(4,4)
      real*8 invSqrt
      integer i, j

      invSqrt = 1.d0 / dsqrt(2.d0)

c     Define unitary transformation matrix U
c     Columns: |1,1>, |1,0>, |1,-1>, |0,0>
      U = 0.d0

c     Column 1: |1,1> = |++> = [1,0,0,0]
      U(1,1) = (1.d0, 0.d0)

c     Column 2: |1,0> = (1/√2)(|+-> + |-+>) = invSqrt*[0,1,1,0]
      U(2,2) = invSqrt
      U(3,2) = invSqrt

c     Column 3: |1,-1> = |--> = [0,0,0,1]
      U(4,3) = (1.d0, 0.d0)

c     Column 4: |0,0> = (1/√2)(|+-> - |-+>) = invSqrt*[0,-1,1,0]
      U(2,4) = -invSqrt
      U(3,4) = invSqrt

c     Compute Kronecker product: O_prod = oper1 ⊗ oper2
      call kronecker(oper1, oper2, O_prod)

c     Compute U† (Hermitian conjugate)
      do i = 1, 4
         do j = 1, 4
            Udag(i,j) = dconjg(U(j,i))
         enddo
      enddo

c     Compute result = U† * O_prod * U using matmul
      temp = matmul(Udag, O_prod)
      result = matmul(temp, U)

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: getPhysicalOper
c     Returns physical isospin operators (raising/lowering/z-component)
c     Based on charge = extQnum - 2:
c        charge = -1: (1/√2)(τ¹ - iτ²) -- lowering operator
c        charge =  0: τ³               -- z-component
c        charge = +1: (1/√2)(τ¹ + iτ²) -- raising operator
c
c     Inputs:
c        extQnum - external quantum number (1, 2, or 3)
c     Output:
c        physOper - 2x2 complex physical operator
c     ----------------------------------------------------------------
      subroutine getPhysicalOper(extQnum, physOper)
      implicit none
      integer extQnum
      complex*16 physOper(2,2)
      complex*16 tau1(2,2), tau2(2,2)
      real*8 invSqrt
      integer charge
      complex*16 im
      parameter (im = (0.d0, 1.d0))

      invSqrt = 1.d0 / dsqrt(2.d0)
      charge = extQnum - 2

      if (charge .eq. -1) then
c        Lowering operator: (1/√2)(τ¹ - iτ²)
         call getPauliMatrix(1, tau1)
         call getPauliMatrix(2, tau2)
         physOper = invSqrt * (tau1 - im * tau2)
      else if (charge .eq. 0) then
c        Z-component: τ³
         call getPauliMatrix(3, physOper)
      else if (charge .eq. 1) then
c        Raising operator: (1/√2)(τ¹ + iτ²)
         call getPauliMatrix(1, tau1)
         call getPauliMatrix(2, tau2)
         physOper = invSqrt * (tau1 + im * tau2)
      endif

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: combinePhysical
c     Helper function that combines physical operators
c     Given extQnum (1, 2, or 3), returns τ_1^phys ⊗ τ_2^phys as 4x4
c     where the physical operator depends on charge = extQnum - 2
c
c     Inputs:
c        extQnum - external quantum number (1, 2, or 3)
c     Output:
c        result  - 4x4 complex matrix in coupled basis
c     ----------------------------------------------------------------
      subroutine combinePhysical(extQnum, result)
      implicit none
      integer extQnum
      complex*16 result(4,4)
      complex*16 physOper(2,2)

c     Get the physical operator for this charge
      call getPhysicalOper(extQnum, physOper)

c     Combine the operator with itself: τ_1^phys ⊗ τ_2^phys
      call combineOpers(physOper, physOper, result)

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: PionPionBC
c     Matrix element in coupled basis of:
c        Σ_{i=1}^3 (τ_1^i τ_2^i) - 2 (τ_1^a τ_2^a)
c     where a = charge=extQnum-2
c
c     Inputs:
c        t, mt   - ket quantum numbers (total isospin, z-projection)
c        tp, mtp - bra quantum numbers
c        extQnum - external quantum number (1, 2, or 3)
c     Output:
c        val     - complex matrix element <t',mt'|O|t,mt>
c     ----------------------------------------------------------------
      subroutine PionPionBC(t, mt, tp, mtp, extQnum, val)
      implicit none
      integer t, mt, tp, mtp, extQnum
      complex*16 val
      integer i, irow, icol
      complex*16 tau1(2,2), tau2(2,2)
      complex*16 combined(4,4), total(4,4)

c     Initialize total matrix to zero
      total = 0.d0

c     Sum over i = 1, 2, 3: add τ^i ⊗ τ^i
      do i = 1, 3
         call getPauliMatrix(i, tau1)
         call getPauliMatrix(i, tau2)
         call combineOpers(tau1, tau2, combined)
         total = total + combined
      enddo

c     Subtract 2 * τ^extQnum ⊗ τ^extQnum
c     call getPauliMatrix(extQnum, tau1)
c     call getPauliMatrix(extQnum, tau2)
      call combinePhysical(extQnum, combined)
      total = total - 2.d0 * combined

c     Get matrix indices from quantum numbers
      call mapping(t, mt, tp, mtp, irow, icol)

c     Return matrix element
      val = total(irow, icol)

      return
      end
c     ----------------------------------------------------------------
c     SUBROUTINE: PionPionA
c     Matrix element in coupled basis of:
c        Σ_{i=1}^3 2*(τ_1^i τ_2^i) - 2 (τ_1^a τ_2^a)
c     where a = extQnum
c
c     Note: This differs from PionPionBC by the factor of 2 in the sum
c
c     Inputs:
c        t, mt   - ket quantum numbers (total isospin, z-projection)
c        tp, mtp - bra quantum numbers
c        extQnum - external quantum number (1, 2, or 3)
c     Output:
c        val     - complex matrix element <t',mt'|O|t,mt>
c     ----------------------------------------------------------------
      subroutine PionPionA(t, mt, tp, mtp, extQnum, val)
      implicit none
      integer t, mt, tp, mtp, extQnum
      complex*16 val
      integer i, irow, icol
      complex*16 tau1(2,2), tau2(2,2)
      complex*16 combined(4,4), total(4,4)

c     Initialize total matrix to zero
      total = 0.d0

c     Sum over i = 1, 2, 3: add 2*τ^i ⊗ τ^i
      do i = 1, 3
         call getPauliMatrix(i, tau1)
         call getPauliMatrix(i, tau2)
         call combineOpers(tau1, tau2, combined)
         total = total + combined
      enddo

c     Subtract 2 * τ^extQnum ⊗ τ^extQnum
c     call getPauliMatrix(extQnum, tau1)
c     call getPauliMatrix(extQnum, tau2)
      call combinePhysical(extQnum, combined)
      total = total -  combined
      total= 2.d0 * total
c     Get matrix indices from quantum numbers
      call mapping(t, mt, tp, mtp, irow, icol)

c     Return matrix element
      val = total(irow, icol)

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: piPhotoOper
c     Matrix element in coupled basis of:
c        (τ_1 · τ_2 - τ_1^z τ_2^z)
c      = τ_1^x τ_2^x + τ_1^y τ_2^y
c      = Σ_{i=1}^3 τ_1^i τ_2^i - τ_1^extQnum τ_2^extQnum
c     This isn't actually used in this folder, but its good for comparison
c     Inputs:
c        t, mt   - ket quantum numbers (total isospin, z-projection)
c        tp, mtp - bra quantum numbers
c        extQnum - external quantum number (typically 3 for z-direction)
c     Output:
c        val     - complex matrix element <t',mt'|O|t,mt>
c     ----------------------------------------------------------------
      subroutine piPhotoOper(t, mt, tp, mtp, extQnum, val)
      implicit none
      integer t, mt, tp, mtp, extQnum
      complex*16 val
      integer i, irow, icol
      complex*16 tau1(2,2), tau2(2,2)
      complex*16 combined(4,4), total(4,4)

c     Initialize total matrix to zero
      total = 0.d0

c     Sum over i = 1, 2, 3: add τ^i ⊗ τ^i
      do i = 1, 3
         call getPauliMatrix(i, tau1)
         call getPauliMatrix(i, tau2)
         call combineOpers(tau1, tau2, combined)
         total = total + combined
      enddo

c     Subtract τ^extQnum ⊗ τ^extQnum
      call getPauliMatrix(extQnum, tau1)
      call getPauliMatrix(extQnum, tau2)
      call combineOpers(tau1, tau2, combined)
      total = total - combined

c     Get matrix indices from quantum numbers
      call mapping(t, mt, tp, mtp, irow, icol)

c     Return matrix element
      val = total(irow, icol)

      return
      end

c     ----------------------------------------------------------------
c     FUNCTION: delta
c     Kronecker delta function: returns 1 if a==b, 0 otherwise
c     ----------------------------------------------------------------
      integer function delta(a, b)
      implicit none
      integer a, b

      if (a .eq. b) then
         delta = 1
      else
         delta = 0
      endif

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: ClosedPiPhoto
c     Closed form expression for pion photoproduction operator
c     val = 2 * (-1)^(t+1) * δ(t,tp) * δ(mt,mtp) * δ(mtp,0)
c
c     Note: This assumes extQnum == 3
c
c     Inputs:
c        t, mt   - ket quantum numbers (total isospin, z-projection)
c        tp, mtp - bra quantum numbers
c        extQnum - external quantum number (must be 3)
c     Output:
c        val     - matrix element <t',mt'|O|t,mt>
c     ----------------------------------------------------------------
      subroutine ClosedPiPhoto(t, mt, tp, mtp, extQnum, val)
      implicit none
      integer t, mt, tp, mtp, extQnum
      real*8 val
      integer delta
      real*8 sign_factor
      if (extQnum.ne.2) then
        write(*,*) "Closed Pion photo requires NEUTRAL pion photoproduction with extQnum=2"
        write(*,*) "passed extQnum=", extQnum 
        error stop
      end if 

      sign_factor=(-1)**(t+1)

c     val = 2 * (-1)^(t+1) * δ(t,tp) * δ(mt,mtp) * δ(mtp,0)
      val = 2.d0 * sign_factor * dble(delta(t, tp))
     &    * dble(delta(mt, mtp)) * dble(delta(mtp, 0))

      return
      end
