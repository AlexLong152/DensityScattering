c     hgrie Oct 2022: v2.0 fewbody-Compton
c     hgrie Aug 2020: v1.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
! hgrie May 2020, adapted from testcompdens.F90
! produce a mathematica-readable .m file with all quantum numbers for specific
! combination of (p12p,p12) and also (p12,p12p), as given by indices
! produce executable via "make produce2Ndensityformath"


! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   testcompdens.F90
! Author: andreasnogga
!
! Created on 20. März 2017, 14:18
!
#include "fdefs.h"

!> program reads in a set of Compton 2Ndensities and produces mathematica-readable file
PROGRAM produce2Ndensityformath
 USE CompDens
 USE precision
 USE meshpoints 
 USE gauss 
 USE parallel
 USE hdf_tool
 USE constants
 USE amplitudes
 USE pwbook
 
 IMPLICIT NONE
  
 INTEGER,PARAMETER :: NP1=48,NP2=1,NP=NP1+NP2
 REAL(dpreal),PARAMETER :: P1=1.0,P2=5.0,P3=15.0
 REAL(dpreal),ALLOCATABLE :: PP(:),PW(:)

 INTEGER,PARAMETER ::  Anucl=3
 INTEGER nsets,iset,np12,ierr
 REAL(dpreal),ALLOCATABLE :: omsets(:),thetasets(:)
 REAL(dpreal) ::  munucleon(-1:1),magmom(-1:1)
 COMPLEX(dpreal) :: op(-1:1,-1:1,-1:1)   ! op(ms3p,ms3,mt3)


 character*200 :: wfFileName  ! abuse of language: for Compton with density, this is name of density file
 
 integer rindx
 integer ip12,l12,s12,j12,t12,mt12,m12 ! automatic integers, so not mutiplied with 2, unlike twoMz 
 integer ip12p,l12p,s12p,j12p,t12p,mt12p,m12p ! automatic integers, so not mutiplied with 2, unlike twoMz
 integer twoMz,twoMzp  ! magnetic quantum numbers of in/out target nucleus, times 2. 
 integer alpha2N,alpha2Np    
!     for outputting 2Ndensity at particular values of momenta (NOT ALL!!) in mathematica friendly format
!     at two momentum-pairs with indices (idxp12,idxp12p) and (idxp12p,idxp12) 
 character*220    :: mathoutputfilename
 integer          :: test        ! a generic integer for testing i/o
 logical          :: existence   ! logical varable whether file exists
 integer          :: idxp12      ! index of first momentum ip12
 integer          :: idxp12p     ! index of second momentum ip12p
 integer          :: jump        ! dummy
 character*22     :: mompair     ! dummy to add momentum pair to output filename
!
 
LOGICAL,PARAMETER :: master=.true.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! prepare grid points to be used 
  ALLOCATE(PP(NP),PW(NP))
  CALL TRNS(NP1,NP2,NP,P1,P2,P3,PP,PW)
  ! initialization of the module precision
  ! prepares machine precisions
  call init_eps

  ! prepare MPI and distributions of processors
  ! parameters are in parallel.dat 
  CALL initparallel
  
  
  ! first printout version running
  IF(master) THEN
   WRITE(*,*) 'ComptonDensity version: ',VERREV
   WRITE(*,*) 'Date                  : ',VERDATE
  END IF

  ! print the machine precision  on master 
  ! and prints compile options etc.   
  IF(master) CALL print_eps
    
  ! printout statistics on the processors used 
  ! info on parallel module   
  CALL printparallel

  ! initialize hdf   
  CALL init_hdf 
  
  ! call preparation routine for bookkeeping
  ! this prepares the partial wave channels for the 
  ! NN,3N and 4N system (even in 3N runs)  
  ! constraints on channels are read in from book-para.dat   
  CALL preppwbook
  
  ! and print the tables with channel numbers 
  IF(master) CALL printpwbook
  

  ! now generate an input file for meshinit that matches 
  ! the generated points 
  
  CALL writemeshpara(PP,PW,NP)
  
  ! call the meshpoint initialization routine
  ! prepares the arrays with grid points and 
  ! integration weights based on the definitions in meshpoints.dat
  ! 4N is prepared too. Should be set to one point per q4 to switch to 
  ! 3N mode   
  CALL meshinit

  ! call the meshpoint output routine to printout grids 
  IF(master) CALL meshprint
    
  ! call the initialization of the amplitudes modules 
  ! this call mainly prepares the distribution of the 
  ! grid points over processors   
  CALL initamp
 
  write(*,*)  "Specify 2Ndensity file name now (in single-inverted-commas!):"
  read(*,*) wfFileName
!  wfFileName = &
!       '../../2Ndensities/tests/compton-dens-3he-idaho-500-empot-3nf-b-om=1.05E+02-th=&
!       5.70E+01-nx=20-nphi=40-np12=48+12-np3=48+12-jmax=2-rho.h5'
 
  INQUIRE(FILE=wfFileName, EXIST=existence)
  if (.not.existence) stop "*** ERROR: Could not open file!!! Aborting."
  
  CALL openread_rho_hdffile(wfFileName)
  np12=SIZE(P12P_density)
  CALL readsets(omsets,thetasets,nsets)
  IF(master) WRITE(*,*) 'Sets in file:' 
  DO iset=1,nsets
   IF(master) WRITE(*,'(A,I4,2E15.6)') 'SETS:',iset,omsets(iset)*hbarc,thetasets(iset)*180.0/pi
END DO
  
  do iset=1,nsets
     CALL readhdf_rhoset(Anucl,omsets(iset),thetasets(iset)) 
  end do
!  np12=NP

  ! now list al momenta in file with grid index
  write(*,*) "  momentum index, momentum [fm^-1] (?)" 
  do jump=1,SIZE(P12P_density)
     write(*,*) jump,P12P_density(jump)
  end do
  write(*,*) 'Number of points in density grid = allowed range of idpx12 and idx12p: ',SIZE(P12P_density)
  
! read data sets until end of file is reached 

!     ****SET NOW WHAT YOU WANT***: momentum-pairs with indices (idxp12,idxp12p) and (idxp12p,idxp12)
  write(*,*) "Specify now momentum-pairs with indices (idxp12,idxp12p) and (idxp12p,idxp12):"
  write(*,*) "Specify first (unprimed) momentum index:  idxp12:"
  read(*,*) idxp12
  write(*,*) "Specify second (primed)  momentum index: idxp12p:"
  read(*,*) idxp12p
!  idxp12  = 2  ! ****SET HERE WHAT YOU WANT*** 
!  jumpidx = 1  ! ****SET HERE WHAT YOU WANT***
!     check if these grid points are actually inside the rage of allowed grid points
  if ((idxp12.gt.SIZE(P12P_density)).or.((idxp12p).gt.SIZE(P12P_density)).or.(idxp12.lt.1).or.((idxp12p).lt.1)) then
     write(*,*) " **** Indices of pair of requested momenta not inside input file -- ABORT. ****"
     stop
  end if
  write(mompair,'(".mompoints-",I0.5,"-",I0.5)') idxp12p,idxp12
!  write(*,*) mompair
  mathoutputfilename = TRIM('2Ndensity-for-mathematica.'//&
       TRIM(wfFileName(INDEX(wfFileName,'/',back=.True.)+1:INDEX(wfFileName,'.',back=.True.)-1))//mompair//'.m')
  write(*,*) '******* Write 2Ndensities to mathematica-readable file: '
  open(unit=15, file=mathoutputfilename,iostat=test)
  if (test .ne. 0) stop "*** ERROR: Could not open output file!!! Aborting."
!     following makes sure we get nice output that can be pasted into mathematica
  
  write(15,*) 'ToExpression[StringReplace["{'
  write(15,'(A)') &
       " {{ ip,  i},{        pp,         p},  index,{αp,  α},{lp,sp,jp,tp,mtp,mp,twoMp},{ l, s, j, t,mt,m,twoM},  ρ2N}"
!     uncomment following if you want range of momenta -- would be VERY sizable file!!!
!         do ip12=7,7     ! 1,P12N   ! unprimed momentum grid ****SET HERE WHAT YOU WANT***
!           do ip12p=4,4        ! 1,P12N   ! primed momentum grid   ****SET HERE WHAT YOU WANT***
!     now produce output
  do jump=0,1
     ip12  = idxp12 + jump* ( idxp12p - idxp12 )
     ip12p = idxp12 + ( 1 - jump) * ( idxp12p - idxp12 )
     write(*,*) "        indices of momentum pair:   (              ",ip12p,",              ",ip12,")"
     write(*,*) "        translate to momentum pair: (",P12P_density(ip12p),",",P12P_density(ip12),")"
     DO alpha2Np=1,num2Nchan
        DO alpha2N=1,num2Nchan
        ! rhoindex determines index for a pair of alpha2N channels
           rindx=rhoindx(alpha2N,alpha2Np)
           IF(rindx.NE.0) THEN  ! if rindx = 0 no coupling of channels, or matrix element not available
              !     produce quantum numbers for channel 
              ! get a two body channel including third components of j12 and Jtot and t12   
              CALL getalpha2N(alpha2Np,l12p,s12p,j12p,t12p,mt12p,m12p,twoMzp)
              ! get a two body channel including third components of j12 and Jtot and t12   
              CALL getalpha2N(alpha2N,l12,s12,j12,t12,mt12,m12,twoMz)
              if (rho(ip12,ip12p,rindx).ne.0) then
                 write(15,25) &
                      ',{{',ip12p,ip12,'},{',P12P_density(ip12p),P12P_density(ip12),'},',& ! {ip12p,ip12},{p12p,p12},
                      rindx,',{',alpha2Np,alpha2N,'},{',& ! rindx,{alpha2np,alpha2n},
                      l12p,s12p,j12p,t12p,mt12p,m12p,twoMzp,'},{',& ! primed quantum #s from alpha12np
                      l12 ,s12 ,j12 ,t12 ,mt12 ,m12 ,twoMz, '},',& ! unprimed quantum #s from alpha12n
                      rho(ip12,ip12p,rindx),'}' ! value of density
              end if
           END IF ! rindx
        END DO ! alpha 2N   
     END DO ! alpha 2Np  
  end do              ! index-pairing
  
  write(15,*) '}","E" -> "*10^"]]'
  close(unit=15,iostat=test)
  if (test .ne. 0) stop "*** ERROR: Could not close output file of 2N densities!!! Aborting."
  write(*,*) '        file name: ',mathoutputfilename
  write(*,*) '        with format:'
  write(*,'(A)') &
       " {{ ip,  i},{     p(ip),      p(i)},  index,{αp,  α},{lp,sp,jp,tp,mtp,mp,twoMp},{ l, s, j, t,mt,m,twoM},  ρ2N}"
  write(*,*) '******* Successfully written 2Ndensities to mathematica-readable file.'
  
! The format is tweaked such that it synchronises with the "header",    || <= here and || <= here
  
25 format(A,I3,",",I3,A,F10.6,",",F10.6,A,I6,A,I3,",",I3,A,6(I2,","),I6,A,6(I2,","),I3,A,E24.15,A)
  
END PROGRAM produce2Ndensityformath
