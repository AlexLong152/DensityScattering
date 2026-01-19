c     subroutine outputKernel(Result,extQnumlimit,twoSnucl)
      subroutine outputKernel(outUnitno,twoSnucl,extQnumlimit,Result,Mnucl,verbosity)

      implicit none
      complex*16,intent(in) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 ScatMat(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      integer extQnumlimit, twoSnucl, outUnitno, verbosity

      call outputroutinelabeled(outUnitno,twoSnucl,extQnumlimit,
     &           Result,verbosity,"Matrix")
      end subroutine outputKernel
