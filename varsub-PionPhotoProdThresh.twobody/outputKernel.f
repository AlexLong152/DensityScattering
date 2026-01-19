c     subroutine outputKernel(Result,extQnumlimit,twoSnucl)
      subroutine outputKernel(outUnitno,twoSnucl,extQnumlimit,Result,Mnucl,verbosity)

      implicit none
      include '../common-densities/constants.def'
      complex*16,intent(in) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      real*8,intent(in) :: Mnucl
      complex*16 :: ResultMat(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      integer extQnumlimit, twoSnucl, outUnitno, verbosity
      real*8 K2N, units,numer,denom,reducedMass

      write(*,'(A)') "Now printing form factor matrix, in units fm^-1"
      call outputroutinelabeled(outUnitno,twoSnucl,extQnumlimit,
     &           Result,verbosity,"FormFactors")
c     In theory you can do math here to generate the matrix elements, but that isn't relevant to my work at the moment
c     Formulation of matrix element can be found in lenkewitz 2011:
c     E_{0+}^{2N}= K2N(F_T^a - F_T^b)
      units=(10.d0**(-3.d0))/mpi0
      reducedMass=Mnucl*mpi0/(Mnucl+mpi0)
      numer=echarge*ga
      denom=(fpi**3)*(16*PI)*((2*PI)**3.d0)
      K2N=(numer*HC/denom)*reducedMass/units  !K2N has 10^-3/mpi units, K2N=mpie g_A 
      resultMat=K2N*Result

      write(*,'(A)') ""
      write(*,'(A)') ""
      write(*,'(A)') "Now printing 2 Body values of E_{0+} and L_{0+}"
      write(*,'(A)') "Differs from matrix elements by 2i, in the Lenkewitz convention if you set:"
      write(*,'(A)') "\mathcal{M} = 2i E_{0+} (\vec{ε}_T· S) + 2i L_{0+} (\vec{ε}_L· S)"
      write(*,'(A)') "See Lenkewitz 2011 equation 2, arXiv:1103.3400 [nucl-th]"
      write(*,*) ""
      write(*,'(A)') "Using units 10^-3/mπ"
      write(*,'(A)') "E_{0+} in extQnum=1,2"
      write(*,'(A)') "L_{0+} in extQnum=3"
      call outputroutinelabeled(outUnitno,twoSnucl,extQnumlimit,
     &           ResultMat,verbosity,"Matrix")
      end subroutine outputKernel
