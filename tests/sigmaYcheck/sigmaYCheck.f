c     =================================================================
c     Test program to compare sigma_y calculations between:
c     1. singlesigmasym (from varsub-spinstructures.f)
c     2. combineOpers with sigma_y matrices (from varsub-IsospinMap.PionPion.f)
c     =================================================================
      program sigmaYCheck
      implicit none
      include '../../common-densities/constants.def'

c     Variables for singlesigmasym
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      real*8 Ax, Ay, Az
      integer Sp, S, verbosity

c     Variables for IsospinMap approach
      complex*16 sigy(2,2), identity(2,2)
      complex*16 sig1y_tens_I(4,4), I_tens_sig2y(4,4)
      complex*16 combined(4,4)

c     Variables for comparison
      integer t, mt, tp, mtp, irow, icol
      complex*16 value1, value2
      real*8 diff
      integer nDiff, nCompared, nTotal

c     Set parameters for singlesigmasym: only Ay component
      Ax = 0.d0
      Ay = 1.d0
      Az = 0.d0
      verbosity = 0
      write(*,*) ""
      write(*,*) "==================================================="
      write(*,*) "Sigma_y Comparison Test"
      write(*,*) "==================================================="
      write(*,*)

c     ====================================================================
c     Method 1: Use IsospinMap routines (calculate once, independent of Sp,S)
c     ====================================================================
c     Get sigma_y matrix using getPauliMatrix (i=2 for sigma_y)

      write(*,*) "Method 1: Kroncker product of pauli matrices"
      write(*,*) "           using IsospinMap"
      call getPauliMatrix(2, sigy)

c     Create identity matrix
      identity = c0
      identity(1,1) = (1.d0, 0.d0)
      identity(2,2) = (1.d0, 0.d0)

c     Compute sigma_1^y ⊗ I using combineOpers
      call combineOpers(sigy, identity, sig1y_tens_I)

c     Compute I ⊗ sigma_2^y using combineOpers
      call combineOpers(identity, sigy, I_tens_sig2y)

c     Add them to get (sigma_1 + sigma_2)_y
      combined = sig1y_tens_I + I_tens_sig2y

c     ====================================================================
c     Method 2: Loop over Sp, S and compare using mapping function
c     ====================================================================
      write(*,*)
      write(*,*) "Method 2: Direct comparison using mapping function"
      write(*,*) "==================================================="

      nTotal = 0
      nDiff = 0

c     Loop over Sp and S values
      do Sp = 0, 1
         do S = 0, 1
            write(*,*)
            write(*,'(A,I1,A,I1)')
     &         "Testing Sp=", Sp, ", S=", S

c           ====================================================================
c           Method 1: Call singlesigmasym with Ay=1
c           This calculates sigma_1^y + sigma_2^y
c           ====================================================================
            call singlesigmasym(hold, Ax, Ay, Az, Sp, S, verbosity)

            nCompared = 0

c           Loop over all possible quantum number combinations
c           t, tp = 0 or 1 (total spin/isospin)
c           mt ranges from -t to +t
c           mtp ranges from -tp to +tp
            do tp = 0, 1
               do t = 0, 1
                  do mtp = -tp, tp
                     do mt = -t, t
c                       Only compare if (tp,t) = (Sp,S)
                        if (tp .eq. Sp .and. t .eq. S) then
                           nCompared = nCompared + 1
                           nTotal = nTotal + 1

c                          Get mapping from (t,mt,tp,mtp) to (irow,icol)
                           call mapping(t, mt, tp, mtp, irow, icol)

c                          Get values from both methods
                           value1 = hold(tp, mtp, t, mt)
                           value2 = combined(irow, icol)

c                          Calculate difference
                           diff = abs(value1 - value2)

c                          Print if different (threshold 1e-8)
                           if (diff .gt. 1.d-8) then
                              nDiff = nDiff + 1
                              write(*,'(A,I1,A,I2,A,I1,A,I2,A)')
     &                           "  DIFF: t=", t, " mt=", mt,
     &                           " tp=", tp, " mtp=", mtp, ":"
                              write(*,'(A,F12.8,A,F12.8,A)')
     &                           "    hold     = (", dreal(value1),
     &                           ",", dimag(value1), ")"
                              write(*,'(A,F12.8,A,F12.8,A)')
     &                           "    combined = (", dreal(value2),
     &                           ",", dimag(value2), ")"
                              write(*,'(A,E12.4)')
     &                           "    |diff|   = ", diff
                              write(*,*)
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo

            write(*,'(A,I4,A)')
     &         "  Compared ", nCompared, " matrix elements for this (Sp,S)"
         enddo
      enddo

      write(*,*)
      write(*,*) "==================================================="
      write(*,'(A,I4,A)')
     &   " Total matrix elements compared: ", nTotal
      write(*,'(A,I4,A)')
     &   " Total differences found:        ", nDiff
      write(*,*)

      if (nDiff .eq. 0) then
         write(*,*) "PASS: All matrix elements match"
         write(*,*) "         for all (Sp,S) combinations"
      else
         write(*,*) "FAIL: Some matrix elements differ"
      endif
      write(*,*)

      end program sigmaYCheck
