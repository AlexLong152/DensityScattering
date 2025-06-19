c ================================
c Poles Functions for SAID Data Analysis
c ================================
      subroutine getCrossSec(sqrtS, x, nuc, crossSec, mNucl,twoSnucl)
      implicit none
c     Input parameters
      double precision sqrtS, x, mNucl
      character*10 nuc
      integer twoSnucl
c     Return value
      double precision crossSec

c     External function
      double precision vecAbs
c     Internal variables
      double precision S
      double precision qVec(3), kVec(3)
      integer i,j

      double complex Mmat(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      double complex MmatDag(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      double complex MmatSquare(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
c     integer*8 theta


      include '../common-densities/constants.def'
c     write(*,*) "Warning, this cross section function only works for spin 1/2 right now"
      S=sqrtS*sqrtS
      
      call getKinematics(sqrtS, x, nuc, S, kVec, qVec, mNucl)
      call getRawM(sqrtS, x, nuc, Mmat, mNucl, twoSnucl,sqrtS)
      ! call printmat(Mmat,"Mmat")

c     write(*,'(A,E30.25)') "For x=",x

c     call dag(Mmat, MmatDag,-twoSnucl,twoSnucl)
      do i = -twoSnucl, twoSnucl
      do j = -twoSnucl, twoSnucl
        Mmatdag(i,j)=dconjg(Mmat(j,i))
      end do
      end do

      MmatSquare=matmul(Mmat,MmatDag)


      crossSec=0.d0
      do i=-twoSnucl,twoSnucl! take the trace
        crossSec=crossSec+real(MmatSquare(i,i),8)
      end do

c     crossSec=real(MmatSquare(1,1)+MmatSquare(2,2),8)!Trace
      crossSec=crossSec*vecAbs(qVec)/vecAbs(kVec)
      crossSec=crossSec*0.25*(HC*HC/(64*S*pi*pi))
c     write(*,*) "crossSec=", crossSec 
      crossSec=crossSec/100
      crossSec=crossSec/(10.d0**(-6.d0))

      end subroutine


      subroutine getRawM(sqrtS, x, nucs, Mout, mNucl,twoSnucl,sqrtSReal)
c     Get the raw matrix element.
c     
c     Parameters
c     ----------
c     sqrtS: double precision
c         The equivalent sqrtS for the kinetmatics we have 
c         if the reaction was single nucleon scattering
c     x: double precision
c         x=cos(theta)
c     nucs: character*3
c         Specifies the reaction ("pp0", "nn0", "pn+", "np-")
c     Mmat: double complex, dimension(2,2) (output)
c         2x2 complex matrix representing the amplitude
c     mNucl: double precision
c         Mass of the nucleon
c     sqrtS: double precision
c         The square root of the Mandelstam variable S
c     ==================================================
      implicit none
      include '../common-densities/constants.def'

c     Input parameters
      double precision sqrtS, x, mNucl,sqrtSReal
      character*3 nucs
      integer twoSnucl
      
c     Output parameter
      complex*16, intent(out) :: Mout(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      
c     Local variables

      double complex Mmat(2,2)

      double precision S
      double precision kVec(3), qVec(3)
      double complex epsVecs(2,3)
      double complex eps(3)
      double precision coefs(3)
      double precision prefactor
      character*10 targets(3)
      character*10 targ

      double complex polContribution(2,2)
      double complex fResult(2,2)
      double complex scaled(2,2)
      integer i, j, k
      double precision sqrtTwo
      
c     Get kinematics
      

c     Define polarization vectors (normalized by sqrt(2))
      sqrtTwo = sqrt(2.0d0)
c     There are two epsilon vectors, left and right circularly polarized
c     eps=(-1,-i,0), eps=(1,-i,0)

      epsVecs(1,1) = dcmplx(-1.0d0, 0.0d0) / sqrtTwo
      epsVecs(1,2) = dcmplx(0.0d0, -1.0d0) / sqrtTwo
      epsVecs(1,3) = dcmplx(0.0d0, 0.0d0) / sqrtTwo
      
      epsVecs(2,1) = dcmplx(1.0d0, 0.0d0) / sqrtTwo
      epsVecs(2,2) = dcmplx(0.0d0, -1.0d0) / sqrtTwo
      epsVecs(2,3) = dcmplx(0.0d0, 0.0d0) / sqrtTwo
      
c     Define target types
      targets(1) = 'p12'
      targets(2) = 'n12'
      targets(3) = '32q'
      
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
          !write(*,*) 'Error: nucs value not supported:', nucs
          Mmat=c0
          return
      endif
      
      call getKinematics(sqrtSReal, x, nucs, S, kVec, qVec, mNucl)
      ! write(*,*) "In poles.f"
      ! write(*,*) "sqrtS=", sqrtS 
      ! write(*,*) "x=", x 
      ! write(*,*) "nucs=", nucs 
      ! write(*,*) "kVec=", kVec 
      ! write(*,*) "qVec=", qVec 


c     Calculate prefactor
      prefactor = 8.0d0 * Pi * sqrtSReal
      
c     Initialize Mmat
      Mmat = dcmplx(0.0d0, 0.0d0)
      
c     Loop over polarizations
      do i = 1, 2
c         Initialize this polarization's contribution
          polContribution = dcmplx(0.0d0, 0.0d0)
          
c         Loop over targets
          do j = 1, 3
c             Calculate F for this target/polarization
              targ=targets(j)
              eps=epsVecs(i,:)
              call f(x, sqrtS, qVec, kVec, eps, targ, fResult)
c             write(*,*) "poles.f:156 fResult=", fResult 
c             write(*,*) "in getRawM"
c             write(*,*) "x,sqrtS=",x, sqrtS 
c             write(*,*) "targ=", targ 
c             write(*,*) "qVec=", qVec 
c             write(*,*) "kVec=", kVec 
c             write(*, '(A,3("(", F8.6, ",", F8.6, ")"))')
c    &        "epsVecs", epsVecs(i,:)
c             write(*,*) "target=", targets(j)
c             do ii = 1, 2
c             do jj = 1, 2
c               write(*,'(A,I1,A,I1,A,F15.10,sp,F15.10,"j")')
c    &          "Mat(",ii,",",jj,")=",  
c    &          real(fResult(ii,jj)), aimag(fResult(ii,jj))
c             end do
c             end do
c             write(*,*) "End getRawM prints"
c             stop

              
c             Scale by coefficient
              do k = 1, 2
                  scaled(k,:) = fResult(k,:) * coefs(j)
              enddo
              
c             Add to this polarization's contribution
              polContribution = polContribution + scaled
          enddo
          
c         Apply prefactor to this polarization's contribution
          do k = 1, 2
              polContribution(k,:) = polContribution(k,:) * prefactor
          enddo
          
c         Add to total matrix
          Mmat = Mmat + polContribution
      enddo
      
c     Cast the "regular" matricies generated from f to the spin indicies 
c     we assigned in main
      Mout(-1,-1)=Mmat(1,1)
      Mout(-1,1)=Mmat(1,2)
      Mout(1,-1)=Mmat(2,1)
      Mout(1,1)=Mmat(2,2)
      return
      end subroutine

      subroutine printmat(mat, name)
      implicit none
      double complex mat(-1:1,-1:1)
      character(*) name
      integer i, j

      do i = -1, 1,2
        do j = -1, 1,2
c         write(*,'(A,"(",I1,",",I1,")=",F15.13,1x,F15.13,"j")')
c    &       trim(name), i, j, real(mat(i,j)), aimag(mat(i,j))
          write(*,'(A,"(",I1,",",I1,")=",E18.10,1x,E18.10,"j")')
     &       trim(name), i, j, real(mat(i,j)), aimag(mat(i,j))

        end do 
      end do
      write(*,*) "" 
      end subroutine printmat


c     ==================================================
      subroutine f(x, sqrtS, qVec, kVec, epsVec, target, result)
c     Calculate F term for given inputs
c     Valid targets are: "p12", "n12", "32q"
c     ==================================================
      implicit none
c     Input parameters
      double precision x, sqrtS
      double precision qVec(3), kVec(3)

      double complex epsVec(3), qVecComplex(3)
      character*10 target
      double complex result(2, 2)
      
c     Local variables
      double complex f1, f2, f3, f4
      double complex f1Term(2, 2), f2Term(2,2)
      double complex f3Term(2, 2)
      double complex f4Term(2, 2)
      double complex crossTmp(3)
      double complex matRes1(2, 2)
      double complex matRes2(2, 2)
      double complex epsDotP
      double precision qAbs, kAbs
      double complex imag

c     External functions
      double precision vecAbs
      external vecAbs, crossProduct, realVecDotSigma, complexvecDotSigma
      
c     External declarations for getF (with correct type declarations)
      external getF
      
c     Initialize imaginary unit
      imag = dcmplx(0.0d0, 1.0d0)
c     Calculate absolute values of vectors
      qAbs = vecAbs(qVec)
      kAbs = vecAbs(kVec)
c     write(*,*) "kVec=", kVec 
c     write(*,*) "qVec=", qVec 
c     write(*,*) "qAbs=", qAbs  
c     write(*,*) "kAbs=", kAbs 
c     Calculate F terms
      if (qAbs.eq.0) then
        !Trick here, just change qAbs, but not q, so that the other terms vanish
        !And we don't get a divide by zero issue
        qAbs=1.d0 
      end if
      call getF(x, sqrtS, 1, target, f1)
      call getF(x, sqrtS, 2, target, f2)
      call getF(x, sqrtS, 3, target, f3)
      call getF(x, sqrtS, 4, target, f4)
c     write(*,*) "f1=", f1 
c     write(*,*) "f2=", f2 
c     write(*,*) "f3=", f3 
c     write(*,*) "f4=", f4 
c     Calculate cross product of kVec and epsVec
      call crossProduct(kVec, epsVec, crossTmp)
      
c     Calculate dot product of qVec and epsVec
c     dotP = qVec(1)*real(epsVec(1)) + qVec(2)*real(epsVec(2)) 
c    &     + qVec(3)*real(epsVec(3))
c     do i = 1, 3
c     qVecComplex=complex(qVec,0.d0)
c     end do
      qVecComplex=dcmplx(qVec,0.d0)

      epsDotP=dot_product(qVec,epsVec)
c     Calculate F1 term: vec · σ (epsVec dot product with Pauli matrices) * F1 * 1j
      call complexVecDotSigma(epsVec, matRes1)
      
      f1Term = matRes1 * f1 * imag
      
c     Calculate F2 term: (qVec · σ) @ [(kVec × epsVec) · σ] * F2
      call realVecDotSigma(qVec, matRes1)
      call complexVecDotSigma(crossTmp, matRes2)
c     call realVecDotSigma(crossTmp, matRes2)

      
      f2Term=matmul(matRes1,matRes2)
      f2Term = f2Term * f2 / (qAbs * kAbs)
      
c     Calculate F3 term: 1j * (kVec · σ) * dot(qVec, epsVec) * F3
      call realVecDotSigma(kVec, matRes1)
      f3Term = matRes1 * imag * epsdotP * f3 / (qAbs * kAbs)
      
c     Calculate F4 term: 1j * (qVec · σ) * dot(qVec, epsVec) * F4
      call realVecDotSigma(qVec, matRes1)
      
      f4Term = matRes1 * imag * epsdotP * f4 / (qAbs * kAbs)
      
      result = f1Term + f2Term + f3Term + f4Term
      
c     write(*,*) "In fortran f subroutine"
c     write(*,*) "target=", target 
c     write(*,*) "x=", x 
c     write(*,*) "sqrtS=", sqrtS 
c     write(*, '(A,3(F10.4))') "qVec=",qVec
c     write(*, '(A,3(F10.4))') "kVec=",kVec
c     write(*, '(A,3("(", F8.6, ",", F8.6, ")"))')
c    &        "epsVecs", epsVec
c     f1-f4 match python
c     write(*,*) "f1=", f1 
c     write(*,*) "f2=", f2 
c     write(*,*) "f3=", f3 
c     write(*,*) "f4=", f4 
c     write(*,*) "f1Term=", f1Term 
c     write(*,*) "f2Term=", f2Term 
c     write(*,*) "f3Term=", f3Term 
c     write(*,*) "f4Term=", f4Term 
c
c     write(*,*) "end fortran f subroutine writes"
c     stop

      return
      end

c     ==================================================
      subroutine getF(x, sqrtS, fi, target, result)
c     Calculate F value for given index
c     x=cos(theta)
c     Fi is which F value is being use, i=1,2,3,4
c     ==================================================
      implicit none
c     Input parameters
      double precision x, sqrtS
      integer fi
      character*10 target
      double complex result
      
c     Local variables
      integer ell
      double complex ePlus, mPlus, eMinus, mMinus  ! Must match type in getPoles
      double complex tmpF, dePlus, dmPlus, deMinus, dmMinus
      double precision legPVal
      
c     External functions
      double precision legP  !defined in utility_suite.f
      external legP, getPoles
      
c     Initialize result
      result = dcmplx(0.0d0, 0.0d0)
      
c     Calculate sum over ell
      do ell = 0, 4
        call getPoles(target, ell, sqrtS, 
     &               ePlus, mPlus, eMinus, mMinus)
        
c       Convert from complex to double complex for consistency
        dePlus = dcmplx(real(ePlus), aimag(ePlus))
        dmPlus = dcmplx(real(mPlus), aimag(mPlus))
        deMinus = dcmplx(real(eMinus), aimag(eMinus))
        dmMinus = dcmplx(real(mMinus), aimag(mMinus))
        
        if (fi .eq. 1) then
c         Case 1: (ell * mPlus + ePlus) * legP(x, ell + 1, deriv=1) +
c                 ((ell + 1) * mMinus + eMinus) * legP(x, ell - 1, deriv=1)
          legPVal = legP(x, ell + 1, 1)
          tmpF = (ell * dmPlus + dePlus) * legPVal
          
          legPVal = legP(x, ell - 1, 1)
          tmpF = tmpF + ((ell + 1) * dmMinus + deMinus) * legPVal
          
          result = result + tmpF
          
        else if (fi .eq. 2) then
c         Case 2: ((ell + 1) * mPlus + (ell * mMinus)) * legP(x, ell, deriv=1)
          legPVal = legP(x, ell, 1)
          tmpF = ((ell + 1) * dmPlus + (ell * dmMinus)) * legPVal
          
          result = result + tmpF
          
        else if (fi .eq. 3) then
c         Case 3: (ePlus - mPlus) * legP(x, ell + 1, deriv=2) +
c                 (eMinus + mMinus) * legP(x, ell - 1, deriv=2)
          legPVal = legP(x, ell + 1, 2)
          tmpF = (dePlus - dmPlus) * legPVal
          
          legPVal = legP(x, ell - 1, 2)
          tmpF = tmpF + (deMinus + dmMinus) * legPVal
          
          result = result + tmpF
          
        else if (fi .eq. 4) then
c         Case 4: (mPlus - ePlus - mMinus - eMinus) * legP(x, ell, deriv=2)
          legPVal = legP(x, ell, 2)
          tmpF = (dmPlus - dePlus - dmMinus - deMinus) * legPVal
          
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
      double complex sigVec(3, 2, 2)
      external matDotVec
      
c     Define Pauli matrices
      sigmaX(1, 1) = dcmplx(0.0d0, 0.0d0)
      sigmaX(1, 2) = dcmplx(1.0d0, 0.0d0)
      sigmaX(2, 1) = dcmplx(1.0d0, 0.0d0)
      sigmaX(2, 2) = dcmplx(0.0d0, 0.0d0)
      
      sigmaY(1, 1) = dcmplx(0.0d0, 0.0d0)
      sigmaY(1, 2) = dcmplx(0.0d0, -1.0d0)
      sigmaY(2, 1) = dcmplx(0.0d0, 1.0d0)
      sigmaY(2, 2) = dcmplx(0.0d0, 0.0d0)
      
      sigmaZ(1, 1) = dcmplx(1.0d0, 0.0d0)
      sigmaZ(1, 2) = dcmplx(0.0d0, 0.0d0)
      sigmaZ(2, 1) = dcmplx(0.0d0, 0.0d0)
      sigmaZ(2, 2) = dcmplx(-1.0d0, 0.0d0)
      
c     Populate Pauli vector
      sigVec(1, 1:2, 1:2) = sigmaX(1:2, 1:2)
      sigVec(2, 1:2, 1:2) = sigmaY(1:2, 1:2)
      sigVec(3, 1:2, 1:2) = sigmaZ(1:2, 1:2)
      
c     Perform matrix-vector multiplication
      call matDotVec(sigVec, vec, result)
      
      return
      end

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
      double complex sigVec(3, 2, 2), tmpVec(3)
      double precision realVec(3)
      integer i
      external matDotVec
      
c     Define Pauli matrices
      sigmaX(1, 1) = dcmplx(0.0d0, 0.0d0)
      sigmaX(1, 2) = dcmplx(1.0d0, 0.0d0)
      sigmaX(2, 1) = dcmplx(1.0d0, 0.0d0)
      sigmaX(2, 2) = dcmplx(0.0d0, 0.0d0)
      
      sigmaY(1, 1) = dcmplx(0.0d0, 0.0d0)
      sigmaY(1, 2) = dcmplx(0.0d0, -1.0d0)
      sigmaY(2, 1) = dcmplx(0.0d0, 1.0d0)
      sigmaY(2, 2) = dcmplx(0.0d0, 0.0d0)
      
      sigmaZ(1, 1) = dcmplx(1.0d0, 0.0d0)
      sigmaZ(1, 2) = dcmplx(0.0d0, 0.0d0)
      sigmaZ(2, 1) = dcmplx(0.0d0, 0.0d0)
      sigmaZ(2, 2) = dcmplx(-1.0d0, 0.0d0)
      
c     Populate Pauli vector
      sigVec(1, 1:2, 1:2) = sigmaX(1:2, 1:2)
      sigVec(2, 1:2, 1:2) = sigmaY(1:2, 1:2)
      sigVec(3, 1:2, 1:2) = sigmaZ(1:2, 1:2)
      
c     Handle complex vector by using real part with matDotVec
c     First handle real part
      do i = 1, 3
        realVec(i) = dreal(vec(i))
      end do
      
c     Initialize result
      result = dcmplx(0.0d0, 0.0d0)
      
c     Compute real part contribution
      call matDotVec(sigVec, realVec, result)
      
c     Now add imaginary part manually (since matDotVec only handles real vectors)
      do i = 1, 3
        tmpVec(i) = dcmplx(0.0d0, 1.0d0) * dimag(vec(i))
      end do
      
c     Manually add the imaginary part contributions
      result(1, 2) = result(1, 2) + tmpVec(1)
      result(2, 1) = result(2, 1) + tmpVec(1)
      
      result(1, 2) = result(1, 2) + tmpVec(2) * dcmplx(0.0d0, -1.0d0)
      result(2, 1) = result(2, 1) + tmpVec(2) * dcmplx(0.0d0, 1.0d0)
      
      result(1, 1) = result(1, 1) + tmpVec(3)
      result(2, 2) = result(2, 2) - tmpVec(3)
      
      return
      end

c     ==================================================
      subroutine getKinematics(sqrtS, x, nucs, S, kVec, qVec, mNucl)
c     Calculate kinematic variables for a given reaction.
c     
c     Parameters
c     ----------
c     sqrtS: double precision
c         The square root of the Mandelstam variable S
c     x: double precision
c         cos(theta)
c     nucs: character*3
c         Specifies the reaction ("pp0", "nn0", "pn+", "np-")
c     S: double precision (output)
c         Mandelstam S
c     kVec: double precision, dimension(3) (output)
c         photon momentum vector
c     qVec: double precision, dimension(3) (output)
c         pion momentum vector
c     mNucl: double precision
c         Mass of the nucleon
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

c     write(*,*) "Epi-mPion=", Epi-mPion 
      if ((Epi.ge.(mPion-2.d0)).and.(Epi.le.mPion+2.d0)) then
        Epi=mPion!Just assuume the user meant to input threshold energy
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
