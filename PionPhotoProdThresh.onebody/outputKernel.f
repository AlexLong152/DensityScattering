c     subroutine outputKernel(Result,extQnumlimit,twoSnucl)
      subroutine outputKernel(outUnitno,twoSnucl,extQnumlimit,Result,verbosity,label)

      implicit none
      complex*16,intent(in) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 ScatMat(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      integer extQnumlimit, twoSnucl, outUnitno, verbosity

      character(len=*),intent(in) :: label

      call outputroutinelabeled(outUnitno,twoSnucl,extQnumlimit,
     &           Result,verbosity,label)
c     In theory you can do math here to generate the matrix elements, but that isn't relevant to my work at the moment
      end subroutine outputKernel
