cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of KERNEL code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              <none>
c     CONTAINS FUNCTIONS:
c              invokesymmetrytwoMzp     : symmetry in Mzp
c              invokesymmetrytwoMz      : symmetry in Mz, may depend on Mzp
c              invokesymmetryextQnum    : symmetry in extQnum, may depend on Mzp and Mz
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New -- this is a STUMP which is not used thus far
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     To contain symmetry considerations and corresponding output messages
c     Can use variable "symmetry" to impose several constraints
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function invokesymmetrytwoMzp(symmetry,twoSnucl,twoMzp,verbosity)

      implicit none
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:

      integer,intent(in) :: symmetry          ! whether and which specific package Usesymmetry() should be used for process -- NOT YET IMPLEMENTED!!!
      integer,intent(in) :: twoSnucl,twoMzp
      
      integer,intent(in) :: verbosity

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set default: no symmetry
      invokesymmetrytwoMzp=.False.

      if ((symmetry.ge.1).and.(twoMzp.le.0.)) invokesymmetrytwoMzp = .True.

      if (twoSnucl.eq.1000) continue
      if (verbosity.eq.1000) continue
      end function invokesymmetrytwoMzp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function invokesymmetrytwoMz(symmetry,twoSnucl,twoMzp,twoMz,verbosity)

      implicit none
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:

      integer,intent(in) :: symmetry          ! whether and which specific package Usesymmetry() should be used for process -- NOT YET IMPLEMENTED!!!
      integer,intent(in) :: twoSnucl,twoMzp,twoMz
      
      integer,intent(in) :: verbosity

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set default: no symmetry
      invokesymmetrytwoMz=.False.

      if (symmetry.eq.1000) continue
      if (twoSnucl.eq.1000) continue
      if (twoMzp.eq.1000) continue
      if (twoMz.eq.1000) continue
      if (verbosity.eq.1000) continue
      end function invokesymmetrytwoMz
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function invokesymmetryextQnum(symmetry,extQnumlimit,extQnum,twoSnucl,twoMzp,twoMz,verbosity)

      implicit none
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
      
      integer,intent(in) :: symmetry          ! whether and which specific package Usesymmetry() should be used for process -- NOT YET IMPLEMENTED!!!
      integer,intent(in) :: extQnumlimit,extQnum,twoSnucl,twoMzp,twoMz
      
      integer,intent(in) :: verbosity

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set default: no symmetry
      invokesymmetryextQnum=.False.

      if ((symmetry.eq.2).and.(extQnum.eq.2)) invokesymmetryextQnum = .True.

      if (extQnumlimit.eq.1000) continue
      if (extQnum.eq.1000) continue
      if (twoSnucl.eq.1000) continue
      if (twoMzp.eq.1000) continue
      if (twoMz.eq.1000) continue
      if (verbosity.eq.1000) continue
      end function invokesymmetryextQnum
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
