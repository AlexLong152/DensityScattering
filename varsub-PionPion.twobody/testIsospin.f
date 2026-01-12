c     =================================================================
c     Test program for Isospin mapping subroutines
c     Based on tests/IsospinMap.py
c     =================================================================

      program testIsospin
      implicit none

      call test1()
      call test2()
      call testGetInd()
      call testPiPhoto()
      call testPiPiA()

      end program testIsospin

c     ----------------------------------------------------------------
c     FUNCTION: getInd
c     Maps twom values to array indices
c     twom=1 -> index=1, twom=-1 -> index=2
c     ----------------------------------------------------------------
      integer function getInd(twom)
      implicit none
      integer twom

      getInd = (1 - twom) / 2 + 1

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: pipiDecoupled
c     Calculates matrix elements in decoupled basis
c     Implements: Σ_i τ_1^i τ_2^i - 2*τ_1^a τ_2^a
c     ----------------------------------------------------------------
      subroutine pipiDecoupled(twom1, twom2, twom1p, twom2p,
     &                         extQnum, val)
      implicit none
      integer twom1, twom2, twom1p, twom2p, extQnum
      complex*16 val
      integer i, ind1, ind2, ind1p, ind2p, getInd
      complex*16 tau1(2,2), tau2(2,2), tau1a(2,2), tau2a(2,2)

      ind1 = getInd(twom1)
      ind2 = getInd(twom2)
      ind1p = getInd(twom1p)
      ind2p = getInd(twom2p)

      val = 0.d0

c     Sum over i = 1, 2, 3
      do i = 1, 3
         call getPauliMatrix(i, tau1)
         call getPauliMatrix(i, tau2)
         val = val + tau1(ind1p, ind1) * tau2(ind2p, ind2)
      enddo

c     Subtract 2*τ^extQnum
      call getPauliMatrix(extQnum, tau1a)
      call getPauliMatrix(extQnum, tau2a)
      val = val - 2.d0 * tau1a(ind1p, ind1) * tau2a(ind2p, ind2)

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: coupledToUncoupled
c     Converts coupled state |T,MT> to uncoupled basis coefficients
c     Basis order: |++>, |+->, |-+>, |-->
c     ----------------------------------------------------------------
      subroutine coupledToUncoupled(T, MT, coefs)
      implicit none
      integer T, MT
      real*8 coefs(4)
      real*8 invSqrt

      invSqrt = 1.d0 / dsqrt(2.d0)
      coefs = 0.d0

      if (T .eq. 1) then
         if (MT .eq. 1) then
c           |1,1> = |++>
            coefs(1) = 1.d0
         else if (MT .eq. 0) then
c           |1,0> = (1/√2)(|+-> + |-+>)
            coefs(2) = invSqrt
            coefs(3) = invSqrt
         else if (MT .eq. -1) then
c           |1,-1> = |-->
            coefs(4) = 1.d0
         endif
      else if (T .eq. 0 .and. MT .eq. 0) then
c        |0,0> = (1/√2)(|+-> - |-+>)
         coefs(2) = -invSqrt
         coefs(3) = invSqrt
      endif

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: DecoupledFromCoupled
c     Calculates matrix element by expanding in decoupled basis
c     ----------------------------------------------------------------
      subroutine DecoupledFromCoupled(t, mt, tp, mtp, extQnum, result)
      implicit none
      integer t, mt, tp, mtp, extQnum
      complex*16 result
      real*8 bra_coefs(4), ket_coefs(4)
      integer twom_states(4,2)
      integer i, j, twom1, twom2, twom1p, twom2p
      real*8 ci, cj
      complex*16 me_unc

c     Define quantum number mapping
      data twom_states / 1, 1, -1, -1,
     &                   1, -1, 1, -1 /

c     Get expansion coefficients
      call coupledToUncoupled(tp, mtp, bra_coefs)
      call coupledToUncoupled(t, mt, ket_coefs)

      result = 0.d0

c     Sum over all combinations
      do i = 1, 4
         ci = bra_coefs(i)
         if (abs(ci) .lt. 1.d-15) cycle

         twom1p = twom_states(i, 1)
         twom2p = twom_states(i, 2)

         do j = 1, 4
            cj = ket_coefs(j)
            if (abs(cj) .lt. 1.d-15) cycle

            twom1 = twom_states(j, 1)
            twom2 = twom_states(j, 2)

            call pipiDecoupled(twom1, twom2, twom1p, twom2p,
     &                         extQnum, me_unc)
            result = result + ci * cj * me_unc
         enddo
      enddo

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: test1
c     Tests DecoupledFromCoupled vs PionPionBC
c     ----------------------------------------------------------------
      subroutine test1()
      implicit none
      integer test_states(4,2)
      integer extQnum, t, mt, tp, mtp, i, j
      complex*16 coupled_result, decoupled_result, diff
      real*8 tolerance
      logical all_passed

      data test_states / 1, 1, 1, 0,
     &                   1, 0, -1, 0 /

      print *, 'Testing DecoupledFromCoupled vs PionPionBC...'
      print *, '============================================',
     &         '===================='

      tolerance = 1.d-5
      all_passed = .true.

      do extQnum = 1, 3
         do i = 1, 4
            t = test_states(i, 1)
            mt = test_states(i, 2)
            do j = 1, 4
               tp = test_states(j, 1)
               mtp = test_states(j, 2)

c              Calculate using coupled basis method
               call PionPionBC(t, mt, tp, mtp, extQnum, coupled_result)

c              Calculate using decoupled basis method
               call DecoupledFromCoupled(t, mt, tp, mtp, extQnum,
     &                                   decoupled_result)

c              Compare results
               diff = coupled_result - decoupled_result
               if (abs(diff) .ge. tolerance) then
                  all_passed = .false.
                  print *, 'FAIL <', tp, ',', mtp, '|O|', t, ',', mt,
     &                     '>: Coupled=', coupled_result,
     &                     ' != ', decoupled_result
               endif
            enddo
         enddo
      enddo

      if (all_passed) then
         print *, 'Test 1 passed: check explicit decoupled vs',
     &            ' coupled result'
      else
         print *, 'In Test 1 Some tests failed'
      endif
      print *, ''

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: test2
c     Tests explicit index and Kronecker product
c     ----------------------------------------------------------------
      subroutine test2()
      implicit none
      integer extQnum, bra, ket
      integer twom_states(4,2)
      integer m1, m2, m1p, m2p, getInd
      complex*16 val, tau2(2,2), tau2_me
      complex*16 Op(4,4), kron_prod(4,4)
      complex*16 tau(2,2)
      integer i
      logical all_passed

      data twom_states / 1, 1, -1, -1,
     &                   1, -1, 1, -1 /

c     Test Kronecker products
      all_passed = .true.
      do extQnum = 1, 3
c        Build Op = Σ_i τ^i ⊗ τ^i - 2*τ^a ⊗ τ^a
         Op = 0.d0
         do i = 1, 3
            call getPauliMatrix(i, tau)
            call kronecker(tau, tau, kron_prod)
            Op = Op + kron_prod
         enddo
         call getPauliMatrix(extQnum, tau)
         call kronecker(tau, tau, kron_prod)
         Op = Op - 2.d0 * kron_prod

c        Check all matrix elements
         do bra = 1, 4
            m1p = twom_states(bra, 1)
            m2p = twom_states(bra, 2)
            do ket = 1, 4
               m1 = twom_states(ket, 1)
               m2 = twom_states(ket, 2)

               call pipiDecoupled(m1, m2, m1p, m2p, extQnum, val)

               if (abs(val - Op(bra, ket)) .ge. 1.d-5) then
                  if (all_passed) then
                     print *, 'Test 2 Failed: explicit index checking',
     &                        ' and Kronecker product'
                  endif
                  all_passed = .false.
                  print *, 'Mismatch at (', m1p, ',', m2p, '|O|',
     &                     m1, ',', m2, ')'
               endif
            enddo
         enddo
      enddo

c     Test tau2 matrix elements
      call getPauliMatrix(2, tau2)
      tau2_me = tau2(getInd(1), getInd(-1))
      if (abs(tau2_me - (-1.d0 * (0.d0, 1.d0))) .ge. 1.d-5) then
         if (all_passed) then
            print *, 'Test 2 Failed: explicit index checking',
     &               ' and Kronecker product'
         endif
         all_passed = .false.
         print *, 'tau2_me(1,-1)=', tau2_me, ' != -i'
      endif

      tau2_me = tau2(getInd(-1), getInd(1))
      if (abs(tau2_me - (0.d0, 1.d0)) .ge. 1.d-5) then
         if (all_passed) then
            print *, 'Test 2 Failed: explicit index checking',
     &               ' and Kronecker product'
         endif
         all_passed = .false.
         print *, 'tau2_me(-1,1)=', tau2_me, ' != i'
      endif

      if (all_passed) then
         print *, 'Test 2 Passed: explicit index checking and',
     &            ' Kronecker product'
      endif
      print *, ''

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: testGetInd
c     Tests for index consistency
c     ----------------------------------------------------------------
      subroutine testGetInd()
      implicit none
      integer extQnum, bra, ket
      integer twom_states(4,2)
      integer m1, m2, m1p, m2p
      complex*16 val, Op(4,4), kron_prod(4,4), tau(2,2)
      integer i
      logical all_pass

      data twom_states / 1, 1, -1, -1,
     &                   1, -1, 1, -1 /

      all_pass = .true.

      do extQnum = 1, 3
c        Build Op matrix
         Op = 0.d0
         do i = 1, 3
            call getPauliMatrix(i, tau)
            call kronecker(tau, tau, kron_prod)
            Op = Op + kron_prod
         enddo
         call getPauliMatrix(extQnum, tau)
         call kronecker(tau, tau, kron_prod)
         Op = Op - 2.d0 * kron_prod

c        Check all matrix elements
         do bra = 1, 4
            m1p = twom_states(bra, 1)
            m2p = twom_states(bra, 2)
            do ket = 1, 4
               m1 = twom_states(ket, 1)
               m2 = twom_states(ket, 2)

               call pipiDecoupled(m1, m2, m1p, m2p, extQnum, val)

               if (abs(val - Op(bra, ket)) .ge. 1.d-5) then
                  if (all_pass) then
                     print *, 'In test_getInd'
                  endif
                  print *, 'getInd broken at <', m1p, ',', m2p,
     &                     '|O|', m1, ',', m2, '>'
                  all_pass = .false.
               endif
            enddo
         enddo
      enddo

      if (all_pass) then
         print *, 'Test 3 Passed: getInd, for index check'
      else
         print *, 'Test Failed: getInd -- for index check'
      endif
      print *, ''

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: testPiPhoto
c     Tests piPhotoOper vs ClosedPiPhoto
c     ----------------------------------------------------------------
      subroutine testPiPhoto()
      implicit none
      integer test_states(4,2)
      integer extQnum, t, mt, tp, mtp, i, j
      complex*16 coupled_result
      real*8 closed_result, diff, tolerance
      logical all_passed

      data test_states / 1, 1, 1, 0,
     &                   1, 0, -1, 0 /

      all_passed = .true.
      tolerance = 1.d-5
      extQnum = 3

      do i = 1, 4
         t = test_states(i, 1)
         mt = test_states(i, 2)
         do j = 1, 4
            tp = test_states(j, 1)
            mtp = test_states(j, 2)

c           Calculate using coupled basis method
            call piPhotoOper(t, mt, tp, mtp, extQnum, coupled_result)

c           Calculate using closed form
            call ClosedPiPhoto(t, mt, tp, mtp, extQnum, closed_result)

c           Compare results
            diff = abs(dble(coupled_result) - closed_result)
            if (diff .ge. tolerance) then
               if (all_passed) then
                  print *, 'Test 4 Failed: closed form vs',
     &                     ' coupled pion photo'
               endif
               all_passed = .false.
               print *, 'FAIL <', tp, ',', mtp, '|O|', t, ',', mt,
     &                  '>: Coupled=', dble(coupled_result),
     &                  ' != ', closed_result
            endif
         enddo
      enddo

      if (all_passed) then
         print *, 'Test 4 passed: closed form vs coupled pion photo'
      endif
      print *, ''

      return
      end

c     ----------------------------------------------------------------
c     SUBROUTINE: testPiPiA
c     Tests PionPionA vs 2*ClosedPiPhoto
c     ----------------------------------------------------------------
      subroutine testPiPiA()
      implicit none
      integer test_states(4,2)
      integer extQnum, t, mt, tp, mtp, i, j
      complex*16 coupled_result
      real*8 closed_result, diff, tolerance
      logical all_passed

      data test_states / 1, 1, 1, 0,
     &                   1, 0, -1, 0 /

      all_passed = .true.
      tolerance = 1.d-5
      extQnum = 3

      do i = 1, 4
         t = test_states(i, 1)
         mt = test_states(i, 2)
         do j = 1, 4
            tp = test_states(j, 1)
            mtp = test_states(j, 2)

c           Calculate using coupled basis method
            call PionPionA(t, mt, tp, mtp, extQnum, coupled_result)

c           Calculate using closed form (factor of 2)
            call ClosedPiPhoto(t, mt, tp, mtp, extQnum, closed_result)
            closed_result = 2.d0 * closed_result

c           Compare results
            diff = abs(dble(coupled_result) - closed_result)
            if (diff .ge. tolerance) then
               if (all_passed) then
                  print *, 'Test 5 Failed: closed form vs',
     &                     ' coupled pi pi A'
               endif
               all_passed = .false.
               print *, 'FAIL <', tp, ',', mtp, '|O|', t, ',', mt,
     &                  '>: Coupled=', dble(coupled_result),
     &                  ' != ', closed_result
            endif
         enddo
      enddo

      if (all_passed) then
         print *, 'Test 5 passed: closed form vs coupled pion photo'
      endif
      print *, ''

      return
      end
