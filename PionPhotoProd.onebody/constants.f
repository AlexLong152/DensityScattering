C-----------------------------------------------------------------------
C     Physical constants for pion photoproduction
C
C     This file initializes essential physical constants and mathematical
C     objects needed for pion photoproduction calculations, including:
C     
C     1. Physical constants (masses, conversion factors)
C     2. Pauli spin matrices for spin-1/2 particle calculations
C
C     Constants are stored in common blocks for access across the codebase.
C-----------------------------------------------------------------------
      SUBROUTINE INIT_CONSTANTS
      IMPLICIT NONE
      INTEGER I, J

C     Common block for physical constants
      DOUBLE PRECISION MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN
      COMMON /PHYS_CONST/ MEVTOFM, MPI, MPIPLUS, MPROTON, MNEUTRON, MN

C     Common block for Pauli matrices
      DOUBLE COMPLEX SIGX(2,2), SIGY(2,2), SIGZ(2,2)
      COMMON /PAULI_MAT/ SIGX, SIGY, SIGZ

C-----------------------------------------------------------------------
C     Initialize physical constants (all in MeV unless noted)
C-----------------------------------------------------------------------
      MEVTOFM = 197.3D0       ! Conversion factor: hc ≈ 197.3 MeV·fm
      MPI = 134.97D0          ! Neutral pion (π⁰) mass 
      MPIPLUS = 139.57D0      ! Charged pion (π⁺) mass
      MPROTON = 938.272D0     ! Proton mass
      MNEUTRON = 939.565D0    ! Neutron mass
      MN = (MPROTON + MNEUTRON) / 2.0D0  ! Average nucleon mass

C-----------------------------------------------------------------------
C     Initialize Pauli spin matrices
C     
C     The Pauli matrices are:
C       σx = ( 0 1 )    σy = ( 0 -i )    σz = ( 1  0 )
C            ( 1 0 )         ( i  0 )         ( 0 -1 )
C
C     These matrices represent spin operators for spin-1/2 particles
C     and satisfy the commutation relations:
C     [σi,σj] = 2i·εijk·σk, where εijk is the Levi-Civita symbol
C-----------------------------------------------------------------------

C     Initialize all matrices to zero
      DO I = 1, 2
        DO J = 1, 2
          SIGX(I,J) = (0.0D0, 0.0D0)
          SIGY(I,J) = (0.0D0, 0.0D0)
          SIGZ(I,J) = (0.0D0, 0.0D0)
        ENDDO
      ENDDO

C     Set non-zero elements for σx
      SIGX(1,2) = (1.0D0, 0.0D0)
      SIGX(2,1) = (1.0D0, 0.0D0)

C     Set non-zero elements for σy
      SIGY(1,2) = (0.0D0, -1.0D0)
      SIGY(2,1) = (0.0D0, 1.0D0)

C     Set non-zero elements for σz
      SIGZ(1,1) = (1.0D0, 0.0D0)
      SIGZ(2,2) = (-1.0D0, 0.0D0)

      RETURN
      END