c     Helper subroutines needed for testing poles.f
c     These are extracted from main.onebodyvia1Ndensity.f

      subroutine commutator(C, A, B, twoSnucl)
c     Calculates the commutator [A,B] = A.B - B.A
      implicit none
      integer, intent(in) :: twoSnucl
      complex*16, intent(in) :: A(-twoSnucl:twoSnucl,
     &                             -twoSnucl:twoSnucl)
      complex*16, intent(in) :: B(-twoSnucl:twoSnucl,
     &                             -twoSnucl:twoSnucl)
      complex*16, intent(out) :: C(-twoSnucl:twoSnucl,
     &                              -twoSnucl:twoSnucl)

      C = matmul(A,B) - matmul(B,A)

      end subroutine commutator


      subroutine checkCommutes(SpinVec, twoSnucl)
c     Checks that the spin matrices in SpinVec satisfy the
c     SU(2) commutation relations
      implicit none
      include '../common-densities/constants.def'
      integer, intent(in) :: twoSnucl
      complex*16, intent(in) :: SpinVec(3,-twoSnucl:twoSnucl,
     &                                    -twoSnucl:twoSnucl)

      complex*16 :: Spinx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 :: Spiny(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 :: Spinz(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 :: t1(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 :: t2(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 :: t3(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      logical passes
      passes=.True.
      Spinx = SpinVec(1,:,:)
      Spiny = SpinVec(2,:,:)
      Spinz = SpinVec(3,:,:)

      call commutator(t1, Spinx, Spiny, twoSnucl)
      t1 = t1 - ci*Spinz

      call commutator(t2, Spiny, Spinz, twoSnucl)
      t2 = t2 - ci*Spinx

      call commutator(t3, Spinz, Spinx, twoSnucl)
      t3 = t3 - ci*Spiny

      if (maxval(abs(t1)) .gt. 1D-10) then
        write(*,*) "ERROR: [Sx,Sy] - i*Sz != 0"
        passes=.False.
      end if
      if (maxval(abs(t2)) .gt. 1D-10) then
        write(*,*) "ERROR: [Sy,Sz] - i*Sx != 0"
        passes=.False.
      end if
      if (maxval(abs(t3)) .gt. 1D-10) then
        write(*,*) "ERROR: [Sz,Sx] - i*Sy != 0"
        passes=.False.
      end if
      if(.not.passes) then
        stop
      end if

      end subroutine checkCommutes
