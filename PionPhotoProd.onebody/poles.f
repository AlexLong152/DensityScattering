      subroutine TraceFromMats(Mats,trace)

      implicit none
      include '../common-densities/constants.def'
      complex*16 Mats(3,-1:1,-1:1), Mat(-1:1,-1:1),MatDag(-1:1,-1:1)
      complex*16 matSquare(-1:1,-1:1)
      real*8, intent(out):: trace
      integer i,j,k,ii
      trace=c0

      do i=1,3
        Mat=Mats(i,:,:)

        do j = -1, 1
        do k = -1, 1
          MatDag(j,k)=dconjg(Mat(k,j))
        end do
        end do

c       Compute M*M† for this polarization
        matSquare=matmul(mat,MatDag)

c       Add trace of this polarization's contribution
        do ii=-1,1
          trace=trace+real(matSquare(ii,ii),8)
        end do!ii

      end do!i
      end subroutine


      subroutine getRawM(sqrtS, x, nucs, Mout, mNucl,sqrtSReal,
     &                   MaxEll,epsVec,threshold)
c     Get the raw matrix element.
c
c     Parameters
c     ----------
c     sqrtS: double precision [MeV]
c         The equivalent sqrtS for the kinematics we have
c         if the reaction was single nucleon scattering
c     x: double precision [dimensionless]
c         x=cos(theta)
c     nucs: character*3
c         Specifies the reaction ("pp0", "nn0", "pn+", "np-")
c     Mout: double complex, dimension(-1:1,-1:1) (output) [dimensionless]
c         2x2 complex matrix representing the scattering amplitude
c     mNucl: double precision [MeV]
c         Mass of the target nucleus
c     sqrtSReal: double precision [MeV]
c         The actual square root of the Mandelstam variable S
c     MaxEll: integer
c         Maximum ell value for pole summation
c     epsVec: complex*16, dimension(3) [dimensionless]
c         Photon polarization vector
c
c     UNITS FLOW:
c       fResult from f(): [MeV^-1]
c       prefactor = 8*pi*sqrtS: [MeV]
c       Mout = prefactor * sum(coefs * fResult): [MeV * MeV^-1] = [dimensionless]
c     ==================================================
      implicit none
      include '../common-densities/constants.def'

c     Input parameters
      double precision sqrtS, x, mNucl,sqrtSReal
      character*3 nucs
      integer MaxEll
      logical threshold

c     Output parameter
      complex*16, intent(out) :: Mout(-1:1,-1:1)

c     Local variables

      double complex Mmat(2,2)

      double precision S
      double precision kVec(3), qVec(3)
c     double complex epsVecs(2,3)
      complex*16, intent(in) :: epsVec(3)

c     double complex eps(3)
      double precision coefs(3)
      double precision prefactor
      character*10 targets(3)
      character*10 targ

      double complex polContribution(2,2)
      double complex fResult(2,2)
c     double complex scaled(2,2)
      double complex E0pSum
      double precision ampThresh
      double precision E_prot, E_neut, L_prot, L_neut

      integer  j, k
      double precision sqrtTwo
      
c     Get kinematics
      

      targets(1) = 'p12'
      targets(2) = 'n12'
      targets(3) = '32q'
      sqrtTwo=2**(0.5d0)
      
c     Set coefficients based on reaction type
      if (nucs .eq. 'pp0') then
          coefs(1) = 1.0d0
          coefs(2) = 0.0d0
          coefs(3) = 2.0d0 / 3.0d0
      else if (nucs .eq. 'nn0') then
          coefs(1) = 0.0d0
          coefs(2) = -1.0d0
          coefs(3) = 2.0d0 / 3.0d0
      else if (nucs .eq. 'pn+') then
          coefs(1) = 1.0d0 * sqrtTwo
          coefs(2) = 0.0d0
          coefs(3) = -1.0d0 / 3.0d0 * sqrtTwo
      else if (nucs .eq. 'np-') then
          coefs(1) = 0.0d0
          coefs(2) = 1.0d0 * sqrtTwo
          coefs(3) = 1.0d0 / 3.0d0 * sqrtTwo
      else
          write(*,*) 'Error: nucs value not supported:', nucs
          stop
          Mmat=c0
          return
      endif
c     Kinematics defined in terms of sqrtSReal, but the poles defined in
c     Terms of the nucleON - photon system

      call getKinematics(sqrtSReal, x, nucs, S, kVec, qVec, mNucl)

c     use sqrtS here bc thats what is built into the single nucleon multipoles
c     the single nucleon multipoles don't "know" about the true nucleus mass
c     UNITS: prefactor has units [MeV], which converts fResult [MeV^-1] to dimensionless Mmat
      prefactor = 8.0d0 * Pi * sqrtS
c     Initialize Mmat
      Mmat = dcmplx(0.0d0, 0.0d0)

c     Initialize this polarization's contribution
      polContribution = dcmplx(0.0d0, 0.0d0)
      
c     At threshold, use hardcoded ChPT amplitudes (Lenkewitz Eqs. 5.23, 5.26)
c     instead of SAID. Above threshold, use SAID as usual.
      if (threshold) then
        if (.not.(Maxell.eq.0)) then!This only works for S-wave scattering
           write(*,*) "Assertion failed: (maxEll.eq.0) evaluated False stopping"
           stop
        end if
c       Hardcoded values in MeV^{-1} (Lenkewitz Eqs. 5.23, 5.26)
c       Original: value × 10^{-3}/M_{pi+}, converted to MeV^{-1}
        E_prot = -1.16d-3 / mpi
        E_neut = 2.13d-3 / mpi
        L_prot = -1.35d-3 / mpi
        L_neut = -2.41d-3 / mpi
c       Select amplitude based on channel and polarization direction
c       Photon travels along z: eps with z-component = longitudinal (L_{0+})
c       eps perpendicular to z = transverse (E_{0+})
        if (nucs .eq. 'pp0') then

          if (abs(epsVec(3)) .eq. 1.d0) then
            ampThresh = L_prot
          else
            ampThresh = E_prot
          end if

        else if (nucs .eq. 'nn0') then

          if (epsVec(3) .eq. 1.d0) then
            ampThresh = L_neut
          else
            ampThresh = E_neut
          end if

        else
          write(*,*) "Error: nucs value not supported for threshold amplitudes: ", nucs
          stop
          ampThresh = 0.0d0
        end if
c       CGLN at threshold: F = i*(eps.sigma)*A_{0+}
c       Scattering matrix: M = 8*pi*sqrtS * F
c       ampThresh is in MeV^{-1}, prefactor in MeV, so Mmat is dimensionless

        call complexVecDotSigma(epsVec, fResult)
        Mmat = fResult * ci * ampThresh * prefactor
      else
        if (epsVec(3).ne.0.d0) then
          write(*,*) "z component of polarization vector component detected above threshold. Check epsVec input."
          write(*,*) "CGLN amplitudes must be extended to F1-F6 with L_+- for longitudinal component, and SAID multipoles only include E_0+, M_1-, etc. for transverse components. Stopping to avoid incorrect results."
          stop
        end if
c       Above threshold: use SAID multipoles via isospin decomposition
        polContribution = dcmplx(0.0d0, 0.0d0)
        do j = 1, 3

          targ = targets(j)
          call f(x, sqrtS, qVec, kVec, epsVec, targ, fResult, MaxEll,
     &           threshold)
          polContribution = polContribution + fResult*coefs(j)

        enddo
        Mmat = Mmat + polContribution * prefactor
      end if
c     enddo
      
c     Cast the "regular" matricies generated from f to the spin indicies 
c     we assigned in main

      call map(Mmat, Mout)
      return
      end subroutine




c     ==================================================
      subroutine f(x, sqrtS, qVec, kVec, epsVec, target, result, MaxEll,
     &             threshold)
c     Calculate F term for given inputs
c     Valid targets are: "p12", "n12", "32q"
c
c     Parameters (with units)
c     -----------------------
c     x: [dimensionless] cos(theta)
c     sqrtS: [MeV] center-of-mass energy
c     qVec: [MeV] pion 3-momentum
c     kVec: [MeV] photon 3-momentum
c     epsVec: [dimensionless] photon polarization vector
c     target: isospin channel identifier
c     result: [MeV^-1] output 2x2 amplitude matrix
c     MaxEll: maximum partial wave
c
c     UNITS FLOW:
c       f1,f2,f3,f4 from getF(): [MeV^-1]
c       epsVec·σ: [dimensionless]
c       qVec·σ, kVec·σ: [MeV]
c       qAbs, kAbs: [MeV]
c       F1Term = (epsVec·σ)*f1*i: [MeV^-1]
c       F2Term = (qVec·σ)@(kVec×epsVec)·σ*f2/(qAbs*kAbs): [MeV*MeV*MeV^-1/MeV^2] = [MeV^-1]
c       F3Term = (kVec·σ)*(qVec·epsVec)*f3/(qAbs*kAbs): [MeV*MeV*MeV^-1/MeV^2] = [MeV^-1]
c       F4Term = (qVec·σ)*(qVec·epsVec)*f4/(qAbs^2): [MeV*MeV*MeV^-1/MeV^2] = [MeV^-1]
c       result: [MeV^-1]
c     ==================================================
      implicit none
c     Input parameters
      double precision x, sqrtS
      double precision qVec(3), kVec(3)

      double complex epsVec(3), qVecComplex(3)
      character*10 target
      integer MaxEll
      logical threshold
      double complex result(2, 2)
      
c     Local variables
      double complex f1, f2, f3, f4
      double complex f1Term(2, 2)
      double complex f2Term(2, 2)
      double complex f3Term(2, 2)
      double complex f4Term(2, 2)
      double complex crossTmp(3)
      double complex matRes1(2, 2)
      double complex matRes2(2, 2)
      double complex epsDotP
      double precision qAbs, kAbs
c     double complex imag

c     External functions
      double precision vecAbs
      external vecAbs, crossProduct, realVecDotSigma, complexvecDotSigma
      
c     External declarations for getF (with correct type declarations)
      external getF
      
      include '../common-densities/constants.def'

c     Initialize imaginary unit
      
c     Calculate absolute values of vectors
      qAbs = vecAbs(qVec)
      kAbs = vecAbs(kVec)
      f1Term=0.d0
      f2Term=0.d0
      f3Term=0.d0
      f4Term=0.d0

      f1=cmplx(0.d0,KIND=8)
      f2=cmplx(0.d0,KIND=8)
      f3=cmplx(0.d0,KIND=8)
      f4=cmplx(0.d0,KIND=8)
      call getF(x, sqrtS, 1, target, f1, MaxEll, threshold)
      call getF(x, sqrtS, 2, target, f2, MaxEll, threshold)
      call getF(x, sqrtS, 3, target, f3, MaxEll, threshold)
      call getF(x, sqrtS, 4, target, f4, MaxEll, threshold)
c     f1=cmplx(real(f1),0.d0)
c     f2=cmplx(real(f2),0.d0)
c     f3=cmplx(real(f3),0.d0)
c     f4=cmplx(real(f4),0.d0)



c     Calculate cross product of kVec and epsVec
      call complexVecDotSigma(epsVec, matRes1)

c     Calculate F1 term: vec · σ (epsVec dot product with Pauli matrices) * F1 * 1j
      f1Term = matRes1 * f1 * ci

      if (qAbs.ne.0.d0) then!avoid divide by zero errors
        write(*,*) "Above threshold"
        stop
        call crossProduct(kVec, epsVec, crossTmp)
        
        qVecComplex=dcmplx(qVec,0.d0)

        epsDotP=dot_product(qVec,epsVec)
c       write(*,*) "epsVec=", real(epsVec)
c       Calculate F2 term: (qVec · σ) @ [(kVec × epsVec) · σ] * F2
        call realVecDotSigma(qVec, matRes1)
        call complexVecDotSigma(crossTmp, matRes2)


        f2Term=matmul(matRes1,matRes2)
        f2Term = f2Term * f2 / (qAbs * kAbs)

c       Calculate F3 term: 1j * (kVec · σ) * dot(qVec, epsVec) * F3
        call realVecDotSigma(kVec, matRes1)
        f3Term = matRes1 * ci * epsdotP * f3 / (qAbs * kAbs)

c     Calculate F4 term: 1j * (qVec · σ) * dot(qVec, epsVec) * F4
        call realVecDotSigma(qVec, matRes1)

        f4Term = matRes1 * ci * epsdotP * f4 / (qAbs * qAbs)

        result = f1Term + f2Term + f3Term + f4Term
c       Poles already converted to MeV^-1 by UNITS_FACTOR in said_subs.f
      else
        result=f1Term
c       Poles already converted to MeV^-1 by UNITS_FACTOR in said_subs.f
      end if      
      
      return
      end

c     ==================================================
      subroutine getF(x, sqrtS, fi, target, result, MaxEll, threshold)
c     Calculate F value for given index by summing over partial waves
c
c     Parameters (with units)
c     -----------------------
c     x: [dimensionless] cos(theta)
c     sqrtS: [MeV] center-of-mass energy
c     fi: which F value (1,2,3,4)
c     target: isospin channel ('p12', 'n12', '32q')
c     result: [MeV^-1] output amplitude
c     MaxEll: maximum ell value for summation
c
c     UNITS FLOW:
c       ePlus, mPlus, eMinus, mMinus from getPoles(): [MeV^-1]
c       Legendre polynomials legP(): [dimensionless]
c       result = sum over ell of (poles * Legendre): [MeV^-1]
c     ==================================================
      implicit none
c     Input parameters
      double precision x, sqrtS
      integer fi, MaxEll
      character*10 target
      logical threshold
      double complex result
      
c     Local variables
      integer ell
c     double complex ePlus, mPlus, eMinus, mMinus  ! Must match type in getPoles
      double complex tmpF, ePlus, mPlus, eMinus, mMinus
      double precision legPVal
      
c     External functions
      double precision legP  !defined in utility_suite.f
      external legP, getPoles !getPoles in said_subs.f
      
c     Initialize result
      result = dcmplx(0.0d0, 0.0d0)

c     Calculate sum over ell
      do ell = 0, MaxEll
        call getPoles(target, ell, sqrtS,
     &               ePlus, mPlus, eMinus, mMinus)

        if (threshold) then
          ePlus  = dcmplx(real(ePlus),  0.d0)
          mPlus  = dcmplx(real(mPlus),  0.d0)
          eMinus = dcmplx(real(eMinus), 0.d0)
          mMinus = dcmplx(real(mMinus), 0.d0)
        end if
        if (fi .eq. 1) then
c         Case 1: (ell * mPlus + ePlus) * legP(x, ell + 1, deriv=1) +
c                 ((ell + 1) * mMinus + eMinus) * legP(x, ell - 1, deriv=1)
          legPVal = legP(x, ell + 1, 1)
          tmpF = (ell * mPlus + ePlus) * legPVal
          
          legPVal = legP(x, ell - 1, 1)
          tmpF = tmpF + ((ell + 1) * mMinus + eMinus) * legPVal
          
          result = result + tmpF
          
        else if (fi .eq. 2) then
c         Case 2: ((ell + 1) * mPlus + (ell * mMinus)) * legP(x, ell, deriv=1)
          legPVal = legP(x, ell, 1)
          tmpF = ((ell + 1) * mPlus + (ell * mMinus)) * legPVal
          
          result = result + tmpF
          
        else if (fi .eq. 3) then
c         Case 3: (ePlus - mPlus) * legP(x, ell + 1, deriv=2) +
c                 (eMinus + mMinus) * legP(x, ell - 1, deriv=2)
          legPVal = legP(x, ell + 1, 2)
          tmpF = (ePlus - mPlus) * legPVal
          
          legPVal = legP(x, ell - 1, 2)
          tmpF = tmpF + (eMinus + mMinus) * legPVal
          
          result = result + tmpF
          
        else if (fi .eq. 4) then
c         Case 4: (mPlus - ePlus - mMinus - eMinus) * legP(x, ell, deriv=2)
          legPVal = legP(x, ell, 2)
          tmpF = (mPlus - ePlus - mMinus - eMinus) * legPVal
          
          result = result + tmpF
          
        else
c         Invalid fi value
          write(*,*) 'Error: Invalid fi value in getF:', fi
          stop
        end if
      end do
      
      return
      end

c     ==================================================
      subroutine realVecDotSigma(vec, result)
c     Multiply a real vector by the Pauli matrices (sigma)
c     Calculates vec · σ where σ = (σx, σy, σz) are the Pauli matrices
c     ==================================================
      implicit none
c     Input parameters
      double precision vec(3)
      double complex result(2, 2)

c     Local variables
      double complex sigmaX(2, 2), sigmaY(2, 2), sigmaZ(2, 2)
      double complex sigmaXmapped(-1:1, -1:1)
      double complex sigmaYmapped(-1:1, -1:1)
      double complex sigmaZmapped(-1:1, -1:1)
      double complex sigVec(3, -1:1, -1:1)
      external map

      include '../common-densities/constants.def'

c     Define Pauli matrices (same convention as complexVecDotSigma)
      sigmaX = c0
      sigmaY = c0
      sigmaZ = c0

      sigmaX(1, 2) = dcmplx(1.0d0, 0.0d0)
      sigmaX(2, 1) = dcmplx(1.0d0, 0.0d0)

      sigmaY(1, 2) = dcmplx(0.0d0, -1.0d0)
      sigmaY(2, 1) = dcmplx(0.0d0, 1.0d0)

      sigmaZ(1, 1) = dcmplx(1.0d0, 0.0d0)
      sigmaZ(2, 2) = dcmplx(-1.0d0, 0.0d0)

c     Map to (-1:1,-1:1) indexing using map subroutine
      call map(sigmaX, sigmaXmapped)
      call map(sigmaY, sigmaYmapped)
      call map(sigmaZ, sigmaZmapped)
    
c     Populate sigVec array with mapped matrices
      sigVec(1, :, :) = sigmaXmapped
      sigVec(2, :, :) = sigmaYmapped
      sigVec(3, :, :) = sigmaZmapped
      sigVec=sigVec/2.d0
c     Verify commutation relations (twoSnucl=1 for spin-1/2)
      call checkCommutes(sigVec, 1)

c     Compute result: vec · σ = vec(1)*σx + vec(2)*σy + vec(3)*σz
      result = vec(1)*sigmaX + vec(2)*sigmaY + vec(3)*sigmaZ

      return
      end

c     ==================================================
      subroutine map(matPoles, matMain)
c     Map a 2x2 matrix from (1:2,1:2) indexing to (-1:1,-1:1) indexing
c     Mapping: index 1 -> +1 (spin up), index 2 -> -1 (spin down)
c     ==================================================
      implicit none
      double complex, intent(in) :: matPoles(2, 2)
      double complex, intent(out) :: matMain(-1:1, -1:1)

      matMain = dcmplx(0.0d0, 0.0d0)
c     (1,1) -> (1,1), (1,2) -> (1,-1), (2,1) -> (-1,1), (2,2) -> (-1,-1)
      matMain(1, 1) = matPoles(1, 1)
      matMain(1, -1) = matPoles(1, 2)
      matMain(-1, 1) = matPoles(2, 1)
      matMain(-1, -1) = matPoles(2, 2)

      return
      end subroutine

c     ==================================================
      subroutine complexVecDotSigma(vec, result)
c     Multiply a complex vector by the Pauli matrices (sigma)
c     Calculates vec · σ where σ = (σx, σy, σz) are the Pauli matrices
c     ==================================================
      implicit none
c     Input parameters
      double complex vec(3)
      double complex result(2, 2)

c     Local variables
      double complex sigmaX(2, 2), sigmaY(2, 2), sigmaZ(2, 2)
      double complex sigmaXmapped(-1:1, -1:1)
      double complex sigmaYmapped(-1:1, -1:1)
      double complex sigmaZmapped(-1:1, -1:1)
      double complex sigVec(3, -1:1, -1:1)
c     integer i
      external map

      include '../common-densities/constants.def'
c     Define Pauli matrices in (1:2,1:2) indexing
      sigmaX = c0
      sigmaY = c0
      sigmaZ = c0

      sigmaX(1, 2) = dcmplx(1.0d0, 0.0d0)
      sigmaX(2, 1) = dcmplx(1.0d0, 0.0d0)

      sigmaY(1, 2) = dcmplx(0.0d0, -1.0d0)
      sigmaY(2, 1) = dcmplx(0.0d0, 1.0d0)

      sigmaZ(1, 1) = dcmplx(1.0d0, 0.0d0)
      sigmaZ(2, 2) = dcmplx(-1.0d0, 0.0d0)

c     Map to (-1:1,-1:1) indexing using map subroutine
      call map(sigmaX, sigmaXmapped)
      call map(sigmaY, sigmaYmapped)
      call map(sigmaZ, sigmaZmapped)
      if (sigmaZmapped(-1,-1).ne.dcmplx(-1,0)) then
        write(*,*) "sigmaZmapped(-1,-1)=", sigmaZmapped(-1,-1) ,"!=-1"
        write(*,*) "oops"
        stop
      end if

c     Populate sigVec array with mapped matrices
      sigVec(1, :, :) = sigmaXmapped
      sigVec(2, :, :) = sigmaYmapped
      sigVec(3, :, :) = sigmaZmapped

c     In case we're feeling extra paranoid...
c     call outputroutinelabeled(10,1,3,
c    &           sigVec,1,"Matrix")

c     Verify commutation relations (twoSnucl=1 for spin-1/2)
      sigVec=sigVec/2.d0 !commutation relations only hold for spin operator, which is sigmaVec/2
      call checkCommutes(sigVec, 1)

c     Compute result using original (1:2,1:2) matrices
      result = vec(1)*sigmaX + vec(2)*sigmaY + vec(3)*sigmaZ

      return
      end
c     ==================================================
      subroutine getKinematics(sqrtS, x, nucs, S, kVec, qVec, mNucl)
c     Calculate kinematic variables for a given reaction.
c
c     Parameters (with units)
c     -----------------------
c     sqrtS: [MeV] The square root of the Mandelstam variable S
c     x: [dimensionless] cos(theta)
c     nucs: Specifies the reaction ("pp0", "nn0", "pn+", "np-")
c     S: [MeV^2] (output) Mandelstam S
c     kVec: [MeV] (output) photon 3-momentum vector
c     qVec: [MeV] (output) pion 3-momentum vector
c     mNucl: [MeV] Mass of the target nucleus
c
c     Internal variables:
c       omega: [MeV] photon energy in CM frame
c       Epi: [MeV] pion energy in CM frame
c       mPion: [MeV] pion mass (mpi0 or mpi from constants.def)
c       absQ: [MeV] magnitude of pion 3-momentum
c     ==================================================
      implicit none
c     Input parameters
      double precision sqrtS, x, mNucl
      character*3 nucs
      
c     Output parameters
      double precision S, kVec(3), qVec(3), mPion
      
c     Local variables
c     double precision mPion
      double precision omega, Epi, absQ
      double precision sqrtTerm
      
      include '../common-densities/constants.def'

c     Constants from PionPhotoLib.py
c     Using mpi (charged) and mpi0 (neutral) from constants.def
c     Using mpi, Mproton, Mneutron from constants.def
      
c     Determine pion mass based on reaction type
      if (nucs .eq. 'pp0') then
          mPion = mpi0  
      else if (nucs .eq. 'nn0') then
          mPion = mpi0  
      else if (nucs .eq. 'pn+') then
          mPion = mpi   
      else if (nucs .eq. 'np-') then
          mPion = mpi   
      else
          write(*,*) 'Error: nucs value not supported: ', nucs
          stop
      endif
      
c     Calculate Mandelstam S
      S = sqrtS**2.d0
      
c     Calculate photon energy
      omega = (S - mNucl*mNucl) / (2.0d0 * sqrtS)
c     write(*,*) "In kinematics poles.f"
c     write(*,*) "S=", S 
c     write(*,*) "mNucl=", mNucl 
c     write(*,*) "sqrtS=", sqrtS 
c     write(*,*) "omega=", omega 
c     write(*,*) "x=",x
c     Set photon momentum vector
      kVec(1) = 0.0d0
      kVec(2) = 0.0d0
      kVec(3) = omega
c     Calculate pion energy
      Epi = (S + mPion**2 - mNucl**2) / (2.0d0 * sqrtS)

      if (abs(Epi-mPion).le.2.5d0) then
        Epi=mPion!Just assume the user meant to input threshold energy
      end if

c     Calculate absolute pion momentum
c     write(*,*) "poles.f:627 Epi**2-mPion**2=", Epi**2-mPion**2
      if (Epi.lt.mPion) then
        write(*,*) "mPion=", mPion
        write(*,*) "Epi=",Epi
        write(*,*) "Epi<mPion, stopping"
        stop
      end if
      absQ = sqrt(Epi**2 - mPion**2)
      
c     Calculate x and y components based on scattering angle (x = cos(theta))
      sqrtTerm = sqrt(1.0d0 - x**2)
      
c     Set pion momentum vector components
      qVec(1) = 0.0d0
      qVec(2) = sqrtTerm * absQ
      qVec(3) = x * absQ
c     write(*,*) "In getKinematics kVec=",kVec      
c     write(*,*) "In getKinematics qVec=",qVec      


c     write(*,*) "In getKinematics"
c     write(*,*) "sqrtS=", sqrtS 
c     write(*,*) "x=", x 
c     write(*,*) "nucs=", nucs
c     write(*,*) "Epi=",Epi
c
c     write(*, '(A,3(F10.4))') "qVec=",qVec
c     write(*, '(A,3(F10.4))') "kVec=",kVec
c     write(*,*) "End getKinematics prints"
c     stop

      return
      end

c     ==================================================
      subroutine printmat(mat, name)
      implicit none
      double complex mat(-1:1,-1:1)
      character(*) name
      integer i, j

      do i = -1, 1,2
        do j = -1, 1,2
c         write(*,'(A,"(",I1,",",I1,")=",F15.13,1x,F15.13,"j")')
c    &       trim(name), i, j, real(mat(i,j)), aimag(mat(i,j))
          write(*,'(A,"(",I2,",",I2,")=",E18.10,1x,E18.10,"j")')
     &       trim(name), i, j, real(mat(i,j)), aimag(mat(i,j))

        end do 
      end do
      write(*,*) "" 
      end subroutine printmat
