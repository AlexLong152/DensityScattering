ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for One/Twobody Contributions to Few-Nucleon Processes Calculated Via 1N/2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              outputroutine        : output result to stdout and output file
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, near-identical to file of same name in common-densities/ of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c
c     twoSmax/twoMz dependence: via array sizes & do-loops

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outputroutine(outUnitno,twoSnucl,extQnumlimit,
     &     Result,verbosity)
c      
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
c     amplitudes with photon helicities 
c     
      complex*16,intent(in) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: twoSnucl
      integer,intent(in) :: outUnitno
      
      integer,intent(in) :: verbosity
      
      integer :: extQnum,twoMzp,twoMz
c     integer i
c**********************************************************************
c     
      if (verbosity.eq.1000) continue
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     output as human readable number to stdout            
c     if (verbosity.ge.0) then
c     write(*,'(A)') REPEAT("%",100) 
c     do extQnum=1,extQnumlimit
c        do twoMzp=twoSnucl,-twoSnucl,-2
c           do twoMz=twoSnucl,-twoSnucl,-2
c                 write (*,'(A,I4,A,I4,A,I4,A,F24.19,SP,F24.19," i")') !E30.19 for exponential form 0.123.....E-56
c    &                 "Result(extQnum=",extQnum,",twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Result(extQnum,twoMzp,twoMz)
c
c                 write (*,'(F24.19,SP,F24.19," i")') !E30.19 for exponential form 0.123.....E-56
c    &                 Result(extQnum,twoMzp,twoMz)
c              end do           ! extQnum
c           end do              ! twoMz
c        end do                 ! twoMzp
c     end if                    ! verbosity
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     output to file !This should be taken care of by run.system file

      write(*,'(A)') REPEAT("%",100) 
      write(outUnitno,'(A)') REPEAT("%",100) 
      do extQnum=1,extQnumlimit
         do twoMzp=twoSnucl,-twoSnucl,-2
            ! do twoMz=twoSnucl,-twoSnucl,-2
            !    write (outUnitno,*) Result(extQnum,twoMzp,twoMz)
            ! end do              ! twoMz
            do twoMz=twoSnucl,-twoSnucl,-2
                  write(outUnitno,'("( ",E24.17,SP,E24.17,"i)")') Result(extQnum,twoMzp,twoMz)
                  write (*,'(A,I4,A,I4,A,I4,A,F24.19,SP,F24.19," i")') !outputs to top of file
     &                 "Result(extQnum=",extQnum,",twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Result(extQnum,twoMzp,twoMz)
               end do           ! extQnum
         end do                 ! twoMzp
      end do                    ! extQnum

      write(*,'(A)') REPEAT("%",100) 
      write(outUnitno,'(A)') REPEAT("%",100) 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     hgrie Aug 2020: if so wanted, output first independent MEs also to screen in a form that can directly be pasted into mathematica
      if (verbosity.ge.0) call outputtomath(Result,twoSnucl,extQnumlimit,verbosity)

      return
      end                       ! outputroutine
