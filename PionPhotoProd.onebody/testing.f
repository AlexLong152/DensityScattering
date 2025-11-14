c ================================
c Poles Functions for SAID Data Analysis
c ================================
      program CheckCrossSec
c     Compiles with `check` command
     &
      implicit none
c     Input parameters
      double precision sqrtS, x, mNucl
      character*10 nuc
      integer twoSnucl, MaxEll
c     Return value
      double precision crossSec

c     External function
      double precision vecAbs
c     Internal variables
      double precision S
      double precision qVec(3), kVec(3)
      integer t
      double complex epsVec(3)
      double complex Mmats(3,-1:1,-1:1)
      double complex Mmat(-1:1,-1:1)
c     integer*8 theta
      integer ieps,targetCalc,thetaIdx
      real*8 theta,omegalab
      real*8 trace
      complex*16 eps(3,3)
c     Experimental values from CompareToFortran.py
      real*8 ccs_pp0(19), ccs_nn0(19)
      real*8 ccExp, diff, pct_err

      include '../common-densities/constants.def'

c     Experimental cross section values for pp0 (proton + gamma -> proton + pi0)
      data ccs_pp0 /
     &  0.1007d0, 0.1027d0, 0.1084d0, 0.1175d0, 0.1293d0,
     &  0.1432d0, 0.1581d0, 0.1729d0, 0.1866d0, 0.1982d0,
     &  0.2070d0, 0.2126d0, 0.2151d0, 0.2151d0, 0.2133d0,
     &  0.2106d0, 0.2080d0, 0.2061d0, 0.2054d0 /

c     Experimental cross section values for nn0 (neutron + gamma -> neutron + pi0)
      data ccs_nn0 /
     &  0.2351d0, 0.2362d0, 0.2391d0, 0.2430d0, 0.2464d0,
     &  0.2477d0, 0.2456d0, 0.2387d0, 0.2265d0, 0.2091d0,
     &  0.1871d0, 0.1618d0, 0.1348d0, 0.1078d0, 0.8286d-1,
     &  0.6152d-1, 0.4522d-1, 0.3501d-1, 0.3154d-1 /

      omegalab=160.d0
      nuc="pp0"
      twoSnucl=1
      mNucl=938.272
      write(*,*) ""
      write(*,'(A)') "Results for:"
      write(*,'(A,F12.5)') "omegalab=",omegalab
      write(*,'(A)') "########################################"
      trace=0.d0

      crossSec=0.d0
      do targetCalc=0,1
        if (targetCalc.eq.1) then
        nuc="nn0"
        mNucl=939.565
        end if

      sqrtS = (2 * omegaLab * mNucl + mNucl * mNucl)**(0.5d0)
      write(*,*) ""
      write(*,'(A,F12.5)') "sqrtS=", sqrtS
      write(*,'(A,A)') "nuc= ", nuc
      thetaIdx = 1
      do t=0,180,10
        theta=real(t,KIND=8)
        x=cos(theta*3.14592d0/180.d0)

        call getKinematics(sqrtS, x, nuc, S, kVec, qVec, mNucl)
        maxEll=4
        eps = RESHAPE((/1,0,0,0,1,0,0,0,0/),(/3,3/))
        do ieps=1,3
          Mmat=c0
          epsVec=eps(:,ieps)
          call getRawM(sqrtS, x, nuc, Mmat, mNucl, twoSnucl,sqrtS,MaxEll,epsVec)
          Mmats(ieps,:,:)=Mmat

        end do!extQnumlimit ieps loop

        call TraceFromMats(Mmats,trace)
c       write(*,*) "poles.f:71 trace=", trace
c       stop
        crossSec=trace*vecAbs(qVec)/vecAbs(kVec)
        crossSec=crossSec*0.25*(HC*HC)/(64*S*pi*pi)
c       write(*,*) "crossSec=", crossSec
        crossSec=crossSec/100
        crossSec=crossSec/(10.d0**(-6.d0))

c       Get experimental value and compute percent difference
        if (targetCalc.eq.0) then
          ccExp = ccs_pp0(thetaIdx)
        else
          ccExp = ccs_nn0(thetaIdx)
        end if
        diff = abs(ccExp - crossSec)
        pct_err = 100.d0 * diff / ccExp

        if (t.eq.0) then
          write(*,'(A)') "theta   Cross Sec   Exp Val     " //
     &                   "Diff      %Error"
        end if
        write(*,'(I3,4X,F8.6,4X,F8.6,4X,F8.6,4X,F6.3,A)')
     &         int(theta), crossSec, ccExp, diff, pct_err, "%"

        thetaIdx = thetaIdx + 1
      end do
      write(*,'(A)') "########################################"
      end do
      end program
