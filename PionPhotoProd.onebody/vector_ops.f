C-----------------------------------------------------------------------
C     Vector operations for pion photoproduction
C
C     This file contains a collection of vector and matrix operations
C     commonly used in the calculation of pion photoproduction amplitudes.
C     These include:
C     - Basic vector operations (magnitude, dot product, cross product)
C     - Matrix operations for 2x2 matrices (multiplication, transpose)
C     - Specialized operations with Pauli matrices
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     Calculate the absolute value (magnitude) of a vector
C
C     Input:
C       VEC(3) - 3D vector [x,y,z]
C
C     Returns:
C       |VEC| = sqrt(x^2 + y^2 + z^2)
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VEC_ABS(VEC)
      IMPLICIT NONE
      DOUBLE PRECISION VEC(3)
      
      VEC_ABS = SQRT(VEC(1)*VEC(1) + VEC(2)*VEC(2) + VEC(3)*VEC(3))
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Calculate the dot product of two vectors
C
C     Input:
C       VEC1(3) - First 3D vector [x1,y1,z1]
C       VEC2(3) - Second 3D vector [x2,y2,z2]
C
C     Returns:
C       VEC1·VEC2 = x1*x2 + y1*y2 + z1*z2
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DOT_PRODUCT(VEC1, VEC2)
      IMPLICIT NONE
      DOUBLE PRECISION VEC1(3), VEC2(3)
      
      DOT_PRODUCT = VEC1(1)*VEC2(1) + VEC1(2)*VEC2(2) + VEC1(3)*VEC2(3)
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Calculate the cross product of two vectors
C
C     Input:
C       VEC1(3) - First 3D vector [x1,y1,z1]
C       VEC2(3) - Second 3D vector [x2,y2,z2]
C
C     Output:
C       RESULT(3) - VEC1 × VEC2 = [y1*z2-z1*y2, z1*x2-x1*z2, x1*y2-y1*x2]
C-----------------------------------------------------------------------
      SUBROUTINE CROSS_PRODUCT(VEC1, VEC2, RESULT)
      IMPLICIT NONE
      DOUBLE PRECISION VEC1(3), VEC2(3), RESULT(3)
      
      RESULT(1) = VEC1(2)*VEC2(3) - VEC1(3)*VEC2(2)
      RESULT(2) = VEC1(3)*VEC2(1) - VEC1(1)*VEC2(3)
      RESULT(3) = VEC1(1)*VEC2(2) - VEC1(2)*VEC2(1)
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Calculate the matrix-vector dot product with Pauli matrices
C
C     This function calculates the dot product of the Pauli matrices
C     with a 3D vector: σ·v = σx*vx + σy*vy + σz*vz
C     
C     This is an important operation in relativistic quantum mechanics
C     for spin-1/2 particles, used in various scattering calculations.
C     
C     Input:
C       VEC(3) - 3D vector [x,y,z]
C
C     Output:
C       RESULT(2,2) - 2x2 complex matrix representing σ·v
C
C     The Pauli matrices are:
C       σx = ( 0 1 )    σy = ( 0 -i )    σz = ( 1  0 )
C            ( 1 0 )         ( i  0 )         ( 0 -1 )
C-----------------------------------------------------------------------
      SUBROUTINE MAT_DOT_VEC(VEC, RESULT)
      IMPLICIT NONE
      DOUBLE PRECISION VEC(3)
      DOUBLE COMPLEX RESULT(2,2)
      INTEGER I, J
      
C     Common block for Pauli matrices
      DOUBLE COMPLEX SIGX(2,2), SIGY(2,2), SIGZ(2,2)
      COMMON /PAULI_MAT/ SIGX, SIGY, SIGZ
      
C     Initialize the result matrix to zero
      DO I = 1, 2
        DO J = 1, 2
          RESULT(I,J) = (0.0D0, 0.0D0)
        ENDDO
      ENDDO
      
C     Calculate σ·v = σx*vx + σy*vy + σz*vz using loops
      DO I = 1, 2
        DO J = 1, 2
          RESULT(I,J) = SIGX(I,J)*VEC(1) + SIGY(I,J)*VEC(2) 
     &                 + SIGZ(I,J)*VEC(3)
        ENDDO
      ENDDO
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Calculate matrix multiplication for 2x2 matrices
C
C     This function performs matrix multiplication for 2x2 complex matrices:
C     RESULT = MAT1 * MAT2
C
C     Input:
C       MAT1(2,2) - First 2x2 complex matrix
C       MAT2(2,2) - Second 2x2 complex matrix
C
C     Output:
C       RESULT(2,2) - The product matrix MAT1 * MAT2
C-----------------------------------------------------------------------
      SUBROUTINE MAT_MUL(MAT1, MAT2, RESULT)
      IMPLICIT NONE
      DOUBLE COMPLEX MAT1(2,2), MAT2(2,2), RESULT(2,2)
      INTEGER I, J, K
      
C     Initialize the result matrix to zero
      DO I = 1, 2
        DO J = 1, 2
          RESULT(I,J) = (0.0D0, 0.0D0)
        ENDDO
      ENDDO
      
C     Perform matrix multiplication using the formula:
C     RESULT(i,j) = SUM_k MAT1(i,k) * MAT2(k,j)
      DO I = 1, 2
        DO J = 1, 2
          DO K = 1, 2
            RESULT(I,J) = RESULT(I,J) + MAT1(I,K) * MAT2(K,J)
          ENDDO
        ENDDO
      ENDDO
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Calculate the conjugate transpose (Hermitian conjugate) of a matrix
C
C     For a complex matrix A, the conjugate transpose A† is defined as:
C     A†(i,j) = conjugate(A(j,i))
C
C     Input:
C       MAT(2,2) - Complex 2x2 matrix
C
C     Output:
C       RESULT(2,2) - Conjugate transpose of the input matrix
C-----------------------------------------------------------------------
      SUBROUTINE CONJUGATE_TRANSPOSE(MAT, RESULT)
      IMPLICIT NONE
      DOUBLE COMPLEX MAT(2,2), RESULT(2,2)
      INTEGER I, J
      
C     Calculate conjugate transpose by taking the complex conjugate
C     and swapping indices
      DO I = 1, 2
        DO J = 1, 2
          RESULT(I,J) = CONJG(MAT(J,I))
        ENDDO
      ENDDO
      
      RETURN
      END

C-----------------------------------------------------------------------
C     Calculate the trace of a matrix
C
C     The trace of a matrix is the sum of its diagonal elements.
C     For a complex matrix, the trace function returns the real part
C     of the sum of diagonal elements.
C
C     Input:
C       MAT(2,2) - Complex 2x2 matrix
C
C     Returns:
C       Real part of (MAT(1,1) + MAT(2,2))
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION TRACE(MAT)
      IMPLICIT NONE
      DOUBLE COMPLEX MAT(2,2)
      DOUBLE COMPLEX TEMP
      INTEGER I
      
C     Sum the diagonal elements
      TEMP = (0.0D0, 0.0D0)
      DO I = 1, 2
        TEMP = TEMP + MAT(I,I)
      ENDDO
      
C     Return the real part of the trace
      TRACE = DBLE(TEMP)
      
      RETURN
      END