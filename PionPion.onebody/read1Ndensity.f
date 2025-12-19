c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018,
c     based on Andreas Nogga's "template" files common-densities/density-modules/testcompdens.F90 in May 2018.
c                                           and common-densities/density-modules/CompDens.F90 in May 2018.
c      
c     Purpose: read in a set of 1N Compton densities.
c     Result is the available as "rho1b(rindx)"
      
c     rho1b(rindx) has units on [MeV]^0/dimensionless.
c     This can be seen from the normalization which is only a sum over quantum numbers.
c     [Andreas email 24 Moay 2018]

c     subroutine get1Nqnnum(rindx,ms3,mt1N,bm,ms3p,mt1Np,bmp,bk,k)
c     from CompDens-complex translates parameters from Andreas' testcompdens.F90:
c      
c     ANDREAS =>  MAIN CODE
c     ms3     = twom1N        : twice spin of in-nucleon
c     mt3     = twomt1N       : twice iso-spin of in-nucleon
c     bm      = twoMz         : twice spin of in-TARGET nucleus
c     ms3p    = twom1Np       : twice spin of out-nucleon
c     mt3p    = twomt1Np      : twice iso-spin of out-nucleon
c     bmp     = twoMzp        : twice spin of out-TARGET/RECOIL nucleus
c     bk      = L1N           : angular momentum of 1N operator (integer starts at 0, up to L1Nmax; "K" in Andreas' notes)
c     k       = ML1N          : mag. quantum of L1Nop (-L1N to L1N: 2L1N+1 integers; "κ" in Andreas' notes)
c
c     Andreas denotes as "3" the nucleon that interacts with photons, while we usually reserve "3" for spectator in 3He
c     
c     in Compton scattering, twomt1N = twomt1Np: isospin of struck nucleon does NOT change.

c     How to translate a set of quantum numbers to rindx (email Andreas Nogga 7 May 2018):

c     in HIS notation: rindx = (ms3+1)/2 + 2*(mt1N+1)/2 + 4*(bm+1)/2 + 8*(ms3p+1)/2 + 16*(mt1Np+1)/2 + 32*(bmp+1)/2 + 64*bk + 128*(k+1) + 1
c
c     HOWEVER I find the following actually matches the output of channels produced by file-read:
c
c     in HIS notation: rindx = 1+ (ms3+1)/2    + 2*(ms3p+1)/2    + 4*(mt3+1)/2     + 8*(Mz+1)/2 +
c                                    16*(Mzp+1)/2    + 32*L1N + 64*(ML1N+L1N)
c     in MY notation:  rindx = 1+ (twom1N+1)/2 + 2*(twom1Np+1)/2 + 4*(twomt1N+1)/2 + 8*(twoMz+1)/2 +
c                                    16*(twoMzp+1)/2 + 32*L1N + 64*(ML1N+L1N)
c
c     So I found NO mt3p, and I find a (k+bk) in the last entry, not a (k+1)!
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
c     TO do:
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c      
c     hgrie June 2018: renamed "parity" to "symmetry -- see notes in usesymmetry+*.f
c       
c     hgrie May 2018: more extensive description of changes in main.*.f
c                     rewritten such that magnetic quantum numbers are "2xMz" etc
c
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read1Ndensity(densityFileName,Anucl,twoSnucl,omega,theta,verbosity)
      
      USE CompDens              ! this needs module CompDens.mod

      IMPLICIT NONE
c**********************************************************************
      include '../common-densities/constants.def'
c**********************************************************************
      
      character*500,intent(in) :: densityFileName  ! abuse of language: for Compton with density, this is name of density file

      real*8,intent(in) :: omega, theta       ! energy (in MeV) and angle (in rad) of input file. Is input to subroutine
                                ! to check that they match what's available in 1Ndensity file. 
      
      integer,intent(in) :: Anucl,twoSnucl
      
      integer,intent(in) :: verbosity              ! verbosity of stdoutout
      
      integer test              ! a generic integer for testing i/o

      character*500      :: mathoutputfilename
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     following are local variables for possible cross checks below
      
      real*8,parameter :: eps=1.d-1 ! a "small parameter" to compare omega and theta to be numerically identical
      
      integer,parameter :: L1Nmax=1  ! hard-wired!
      integer L1N, ML1N              ! orb ang mom & projection of 1N operator
      integer twomt1N,twom1N,twom1Np ! projections of single-nucleon's isospin, in-spin, out-spin
      integer twomt1Np               ! only for cross check: in Compton, twomt1Np = twomt1N
      integer twoMz,twoMzp
      integer rindx                  ! Andreas' index which combines all quantum numbers
      real*8 :: ff,ffn,ffp
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) "*********************** 1N DENSITY MATRIX PARAMETERS ***************************"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      CALL openread_rho1bfile(densityFileName)

c now assume that there is a consistent data set for the one body case (same q etc))
c     read this set
      CALL read_rho1bset(Anucl) ! Anucl is passed into routine to return FF result to STdoUT

      write(*,*) "        q [fm^-1]      F_n            F_p            F_N"
      write(*,*) "Computed above are Form Factors of Nucleon, normalised to F(q²=0)=1 for A=3,4; =2S+1 else."
      
c     output form factors * number of nucleons: cross check not using the densities grid but the integration grid and j12max
      ff    = 0.d0
      ffn   = 0.d0
      ffp   = 0.d0

      write(*,*) "Computing Form Factors of N, normalised to F_pairs(q²=0)=1."
      write(*,'(AI3A)') " Total # nucleons:",Anucl,"    n              p              nucleons"
      do rindx=1,maxrho1bindex
         CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
         if ((L1N.eq.0).and.(twom1N.eq.twom1Np).and.(twomt1N.eq.twomt1Np).and.(twoMz.eq.twoMzp)) then
            if (twomt1N.eq.1) then
               ffp = ffp + rho1b(rindx)/(twoSnucl+1)
            else if (twomt1N.eq.-1) then
               ffn = ffn + rho1b(rindx)/(twoSnucl+1)
            end if
            ff     = ff + rho1b(rindx)/(twoSnucl+1)
         end if
      end do ! rindx
      write(*,'(A(4F15.6))') "  Form Factors:   ",ffn,ffp,ff
      write(*,'(A(4F15.6))') "  Eff. # nucleons:",Anucl*ffn,Anucl*ffp,Anucl*ff

      
      if (verbosity.ge.1) then
         write(*,*) "  1N density matrix: (only mt1Np=mt1N for Compton: no 1N-isospin change)"
         write(*,*) "   rindx:   twoMzp,    twoMz, twomt1Np,  twomt1N,  twom1Np,   twom1N,      L1N,     ML1N:     ρ1(rindx)"
         do rindx=1,maxrho1bindex
c following do-loops are alternatives, where we only show the MEs necessary by symmetry!            
c         do twoMzp=twoSnucl,twoMzplimit,-2        ! only the first few: symmetry!
c         for Mzp=0, run only over Mz>=0 -- that's still 2 more than necessary, but good enough
c            if ((twoMzp.eq.0).and.(twoMzplimit.eq.0)) then
c                twoMzlimit = 0
c                else
c                twoMzlimit = -twoSnucl
c            end if   
c            do twoMz=twoSnucl,twoMzlimit,-2
c               do twomt1Np=1,-1,-2
c                  do twomt1N=1,-1,-2
c                     do twom1Np=1,-1,-2
c                        do twom1N=1,-1,-2
c                           do L1N=0,L1Nmax
c                              do ML1N=-L1N,L1N
c                                 write(*,30) "  ",rindx,": ",twoMzp,twoMz,twomt1Np,twomt1N,twom1Np,twom1N,L1N,ML1N,rho1b(rindx)
            CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
c     don't forget to convert Andreas' quantum numbers to ours
            write(*,30) "  ",rindx,":",twoMzp,twoMz,twomt1Np,twomt1N,twom1Np,twom1N,L1N,ML1N,":",rho1b(rindx)
c                                 write(*,30) "  ",rindx,": ",twoMzp,twoMz,twomt1Np,twomt1N,twom1Np,twom1N,L1N,ML1N,": ",rho1b(rindx)
c                              end do !ML1N
c                           end do    !L1N
c                        end do       !twom1N
c                     end do          !twom1Np
c                  end do             !twomt1N
c               end do                !twomt1Np
c            end do                   !twoMz
c         end do                      !twoMzp
        end do !rindx
      end if
c      
c     for large verbosity number, create file "1Ndensity-formath.XXX.m" of the density which can be read in by mathematica
c
      if (verbosity.ge.4) then
         mathoutputfilename = TRIM('1Ndensity-for-mathematica.'//
     &        TRIM(densityFileName(INDEX(densityFileName,'/',back=.True.)+1:INDEX(densityFileName,'.',back=.True.)-1))//'.m')
         open(unit=15, file=mathoutputfilename,iostat=test)
         if (test .ne. 0) stop "*** ERROR: Could not open output file!!! Aborting."
c     following makes sure we get nice output that can be pasted into mathematica
         write(15,*) 'ToExpression[StringReplace["{'
         write(15,'(9(A,","),A)')  " {index","twoMp","twoM","twomtp","twomt","twomp"," twom","    K","    κ","    ρ1N}"
         do rindx=1,maxrho1bindex
c     don't forget to convert Andreas' quantum numbers to ours
            CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
            write(15,26) "  ",rindx,",",twoMzp,twoMz,twomt1Np,twomt1N,twom1Np,twom1N,L1N,ML1N," ",rho1b(rindx)
         end do
         write(15,*) '}","E" -> "*10^"]]'
         close(unit=15,iostat=test)
         if (test .ne. 0) stop "*** ERROR: Could not close output file of 1N densities!!! Aborting."
         write(*,*) '   ******* Wrote 1Ndensities to mathematica-readable file named: '
         write(*,*) '       ',mathoutputfilename
      end if
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     perform checks that density file parameters match requested parameters -- exit if not.

c         if (Anuclhere.ne.Anucl) then
c            write(*,*) "*** ERROR: Anucl of 1Ndensity file does not match:"
c            write(*,*)             "    Anucl = ",Anuclhere,"   <===>   Anucl1N = ",Anucl
c            write(*,*) "-- Exiting."
c            stop
c         end if

c         if (nsets.ne.1) then
c            write(*,*) "*** ERROR: 1Ndensity file contains more than 1 energy/angle combination --- "
c            write(*,*)             " --- THAT IS NOT YET IMPLEMENTED ---"
c            write(*,*) "-- Exiting."
c            stop
c         end if   
      
         if (abs(omega-omval*hc).ge.eps) then
            write(*,*) "*** ERROR: Energy of 1Ndensity file does not match:"
            write(*,*)             "    omega = ",omega," MeV   <===>   omega1N = ",omval*hc," MeV"
            write(*,'(A,F15.6,A)') "    omega - omega(1Ndensity) = ",omega-omval*hc," MeV"
            write(*,*) "-- Exiting."
            stop
         else if (abs(theta-thetaval).ge.eps) then
            write(*,*) "*** ERROR: Angle of 1Ndensity file does not match:"
            write(*,*)             "    theta = ",theta*180.0d0/Pi," deg   <===>   theta1N = ",thetaval*180.0d0/Pi," deg"
            write(*,'(A,F15.6,A)') "    theta - theta(1Ndensity) = ",(theta-thetaval)*180.0d0/Pi," deg"
            write(*,*) "-- Exiting."
            stop
         else
            write(*,*)             "   Energy and angle of 1Ndensity file match request." 
         end if

      if (verbosity.eq.1000) continue
      if (twoSnucl.eq.1000) continue
      
      return
          
 26   format(',{',A,I3,A,8(I5,","),A,E24.15,"}")  
 30   format(' ',A,I4,A,8I10,A,F20.13,SP,F21.13," I")
      
      end
