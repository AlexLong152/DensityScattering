      PROGRAM onebodytest
      USE CompDens              ! needs module CompDens.mod

      IMPLICIT NONE


      include '../common-densities/constants.def'
c**********************************************************************
      real*8 mNucl
      real*8 masses(2)
      real*8 x,sqrtS,crossSec, omegalab,ccVal
      real*8 ccs(19)
      integer theta,i,j
      character(len=3) :: nucs(2)
      character*3 nuc
      masses=(/mProton,mNeutron/)
      nucs(1)="pp0"
      nucs(2)="nn0"

      do i=1,2
      write(*,*) ""
      if (i.eq.1) then
        write(*,*) "proton, neutral pion"
        ccs=(/
     &  0.1007E+00,
     &  0.1027E+00,
     &  0.1084E+00,
     &  0.1175E+00,
     &  0.1293E+00,
     &  0.1432E+00,
     &  0.1581E+00,
     &  0.1729E+00,
     &  0.1866E+00,
     &  0.1982E+00,
     &  0.2070E+00,
     &  0.2126E+00,
     &  0.2151E+00,
     &  0.2151E+00,
     &  0.2133E+00,
     &  0.2106E+00,
     &  0.2080E+00,
     &  0.2061E+00,
     &  0.2054E+00 /)
      end if

      if (i.eq.2) then
        write(*,*) "neutron, neutral pion"

        ccs=(/
     &  0.2351E+00,  
     &  0.2362E+00, 
     &  0.2391E+00,
     &  0.2430E+00,
     &  0.2464E+00,
     &  0.2477E+00,
     &  0.2456E+00,
     &  0.2387E+00,
     &  0.2265E+00,
     &  0.2091E+00,
     &  0.1871E+00,
     &  0.1618E+00,
     &  0.1348E+00,
     &  0.1078E+00,
     &  0.8286E-01,
     &  0.6152E-01,
     &  0.4522E-01,
     &  0.3501E-01,
     &  0.3154E-01/)
      end if
      nuc=nucs(i)
      mNucl=masses(i)

      omegalab=160.d0
      sqrtS=sqrt(2.d0*omegalab*mNucl+mNucl*mNucl)
      write(*,*) "sqrtS=", sqrtS,"omegaLab=",omegaLab
      write(*,*) ""
c     sqrtS=omega+sqrt(omega*omega+mNucl*mNucl)
      do theta=0,180,10
        j=1+theta/10
        x=cos(real(theta,8)*Pi/180.d0)
        call getCrossSec(sqrtS,x,nuc,crossSec,mNucl,1)
        ccVal=ccs(j)
        write(*,'(A,I4,A,F8.6,A,F8.6,A,F8.4,A,F8.3,A)')  'theta=',theta, ', cross sec=',crossSec,' Exp cc=',ccVal, ' ----- Expermental diff=',ccVal-crossSec,
     &     '-->',100*(ccVal-crossSec)/ccVal,'% error'
      end do!theta
      end do!i
      end PROGRAM
