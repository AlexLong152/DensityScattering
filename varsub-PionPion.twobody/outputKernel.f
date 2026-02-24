c     subroutine outputKernel(Result,extQnumlimit,twoSnucl)
      subroutine outputKernel(outUnitno,twoSnucl,extQnumlimit,Result,Mnucl,verbosity)

      implicit none
      complex*16,intent(in) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 ScatMat(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      integer extQnumlimit, twoSnucl, outUnitno, verbosity
      real*8 Mnucl
      write(*,'(A)') "Pion charge = extQnum-2 : extQnum=1 -> π-, extQnum=2 -> π0, extQnum=3 -> π+ "
      write(*,'(A)') "Scattering length units 10^-3/Mpi"
      call outputroutinelabeled(outUnitno,twoSnucl,extQnumlimit,
     &           Result,verbosity,"Scattering Length")
      end subroutine outputKernel
