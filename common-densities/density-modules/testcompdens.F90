!
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

!> program testcompdens reads in a set of Compton densities 
PROGRAM testcompdens
 USE CompDens
 USE precision
 USE meshpoints 
 USE gauss 
#ifdef RHODISTR
 USE parallel
 USE mpi_const
 USE hdf_tool
 USE constants
 USE amplitudes
 USE pwbook
#endif
 
 IMPLICIT NONE
  
 INTEGER,PARAMETER :: NP1=20,NP2=10,NP=NP1+NP2
 REAL(dpreal),PARAMETER :: P1=1.0,P2=5.0,P3=15.0
 REAL(dpreal),ALLOCATABLE :: PP(:),PW(:)

 INTEGER,PARAMETER ::  Anucl=3
 INTEGER ip12,alpha2N,l12,s12,j12,t12,mt12,m12,mjtot,rindx,ip12loc
 REAL(dpreal) sum,sumnp,sumpp,sumnn,sumgl,sumglnp,sumglpp,sumglnn
 REAL(dpreal) sumn,sump
 INTEGER nsets,iset,np12,shiftp12,ierr
 INTEGER ms3,mt3,bm,ms3p,mt3p,bmp,bk,k
 REAL(dpreal),ALLOCATABLE :: omsets(:),thetasets(:)
 REAL(dpreal) ::  munucleon(-1:1),magmom(-1:1)
 COMPLEX(dpreal) :: op(-1:1,-1:1,-1:1)   ! op(ms3p,ms3,mt3)
#ifndef RHODISTR
LOGICAL,PARAMETER :: master=.true.
#endif
  
  ! prepare grid points to be used 
  ALLOCATE(PP(NP),PW(NP))
  CALL TRNS(NP1,NP2,NP,P1,P2,P3,PP,PW)

#ifdef RHODISTR  
  ! initialization of the module precision
  ! prepares machine precisions
  call init_eps

  ! prepare MPI and distributions of processors
  ! parameters are in parallel.dat 
  CALL initparallel
  
  
  ! first printout version running
  IF(master) THEN
   WRITE(*,*) 'ComptonDensity version: ',VERREV
   ! WRITE(*,*) 'Date                  : ',VERDATE
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
 
  CALL openread_rho_hdffile('compton-rho.h5')
  CALL readsets(omsets,thetasets,nsets)
  IF(master) WRITE(*,*) 'Sets in file:' 
  DO iset=1,nsets
   IF(master) WRITE(*,'(A,I4,2E15.6)') 'SETS:',iset,omsets(iset)*hbarc,thetasets(iset)*180.0/pi
  END DO
#else
  CALL openread_rhofile('compton-rho.dat',PP,NP)
#endif

  CALL openread_rho1bfile('compton-rho1b.dat')
  
#ifndef RHODISTR   
 ! read data sets until end of file is reached 
 DO WHILE(.NOT. rho_eof)
  CALL read_rhoset(Anucl)
  IF(rho_eof) EXIT 
  np12=NP
  shiftp12=0
#else  
 
 DO iset=1,nsets
  CALL readhdf_rhoset(Anucl,omsets(iset),thetasets(iset)) 
  np12=mynp12
  shiftp12=myp12-1
#endif  

  ! check normalization
  sum=0.0_dpreal   ! total 
  sumnp=0.0_dpreal  ! mt12=0
  sumpp=0.0_dpreal  ! mt12=1
  sumnn=0.0_dpreal ! mt12=-1 

  ! this loops only access the diagonal elements of rho 
  DO alpha2N=1,num2Nchan
    ! get a two body channel including third components of j12 and Jtot and t12   
    CALL getalpha2N(alpha2N,l12,s12,j12,t12,mt12,m12,mjtot)
    ! rhoindex determines index for a pair of alpha2N channels
    !  (here only diagonal)
    rindx=rhoindx(alpha2N,alpha2N)
    IF(rindx.NE.0) THEN  ! if rindx = 0 no coupling of channels, or matrix element not available 
     IF(mt12.EQ.-1) THEN ! nn should not exist
      DO ip12loc=1,np12
       ip12=ip12loc+shiftp12
       ! access rho(p12',p12,rindx) here only diagonal     
       sumnn=sumnn+rho(ip12,ip12loc,rindx)*PP(ip12)**2*PW(ip12)
       sum=sum+rho(ip12,ip12loc,rindx)*PP(ip12)**2*PW(ip12) 
      END DO
     ELSE IF(mt12.EQ.0) THEN ! np -> spectator is p 
      DO ip12loc=1,np12
       ip12=ip12loc+shiftp12
       sumnp=sumnp+rho(ip12,ip12loc,rindx)*PP(ip12)**2*PW(ip12) 
       sum=sum+rho(ip12,ip12loc,rindx)*PP(ip12)**2*PW(ip12) 
      END DO
     ELSE IF(mt12.EQ.1) THEN ! pp -> spectator is n 
      DO ip12loc=1,np12
       ip12=ip12loc+shiftp12
       sumpp=sumpp+rho(ip12,ip12loc,rindx)*PP(ip12)**2*PW(ip12) 
       sum=sum+rho(ip12,ip12loc,rindx)*PP(ip12)**2*PW(ip12) 
      END DO          
     END IF  ! mt12
    END IF ! rindx
  END DO ! alpha 2N 
  
#ifdef RHODISTR
  CALL MPI_ALLREDUCE(sum,sumgl,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sumnp,sumglnp,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sumnn,sumglnn,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sumpp,sumglpp,1,MPI_REAL8,MPI_SUM,commp12,ierr)     
#else
  sumgl=sum
  sumglnp=sumnp
  sumglpp=sumpp
  sumglnn=sumnn
#endif  
  
  
  IF(master) THEN   
   IF(Anucl.EQ.3) THEN  ! A=3 includes sum over M=-1/2, +1/2 
    WRITE(*,'(A,5E15.6)')  'FF: ',qval/2.0,sumglnn/2.0_dpreal, &
    &   sumglnp/2.0_dpreal,sumglpp/2.0_dpreal,sumgl/2.0_dpreal
   ELSE
    WRITE(*,'(A,5E15.6)')  'FF: ',qval,sumglnn,sumglnn,sumglnp,sumglpp
   END IF ! Anucl 
  END IF
  
  ! now assume that there is a consistent data set for the one body case (same q etc))
  ! read this set   
  CALL read_rho1bset(Anucl)
!   IF(.NOT. rho1b_eof) EXIT  ! just in case something is wrong 
  
  IF(master) THEN ! rho1b is not distribute, calculate norm on master 
  sump = 0._DPREAL
  sumn = 0._DPREAL
  sum  = 0._DPREAL
  
    ! test the data by calculating norm/ff
   DO rindx=1,maxrho1bindex
    CALL get1Nqnnum(rindx,ms3,mt3,bm,ms3p,mt3p,bmp,bk,k)

    IF(bk.eq.0 .and. k.eq.0 & 
   & .and. ms3.eq.ms3p .and. mt3.eq.mt3p &
   & .and. bm.eq.bmp ) THEN
     IF(mt3.eq.1) THEN ! proton 
      sump=sump+rho1b(rindx)
     END IF 
     IF(mt3.eq.-1) THEN ! neutron 
      sumn=sumn+rho1b(rindx)
     END IF 
     sum=sum+rho1b(rindx)
    END IF
   END DO ! rindx
  
 
   IF(Anucl.EQ.3) THEN  ! A=3 includes sum over M=-1/2, +1/2 
    WRITE(*,'(A,5E15.6)')  'FF1B (qval,neutron,proton,all): ',qval,0.5*sumn,0.5*sump,0.5*sum
   ELSE
    WRITE(*,'(A,5E15.6)')  'FF1B (qval,neutron,proton,all): ',qval,sumn,sump,sum
   END IF
     
  END IF ! master 
 
  IF(master) THEN ! rho1b is not distributed, calculate magnetic moment on master   
!  first I need to the define the matrix elements of the one body operator 
! proton and neutron magnetic moment
   munucleon(1)=2.793     ! proton
   munucleon(-1)=-1.913   ! neutron
! the matrix element of O is nonzero only for K=0 and k=0 since 
! there is no dependence on the nucleon momentum of this operator
!  <ms3p,mt3p | munucleon*sigma_z | ms3,mt3 > = delta_ms3ms3p delta_mt3mt3p *munucleon(mt3)*(-1)^((1-mt3)/2)  
   op=0.0
   op(1,1,1)=munucleon(1)
   op(-1,-1,1)=-munucleon(1)
   op(1,1,-1)=munucleon(-1)
   op(-1,-1,-1)=-munucleon(-1)
 
   magmom = 0._DPREAL
  
    ! test the data by calculating norm/ff
   DO rindx=1,maxrho1bindex
    CALL get1Nqnnum(rindx,ms3,mt3,bm,ms3p,mt3p,bmp,bk,k)

    IF(bk.eq.0 .and. k.eq.0 & 
   & .and. mt3.eq.mt3p &
   & .and. bm.eq.bmp ) THEN
     magmom(bm)=magmom(bm)+rho1b(rindx)*op(ms3p,ms3,mt3)
    END IF
   END DO ! rindx
  
   WRITE(*,'(2(A,E15.6),E15.6)')  'MAGNETIC MOMENT from inserting σz at qval = ',qval,": ",3*magmom(1),3*magmom(-1)
   
! now the same for sigma_x
! the matrix element of O is nonzero only for K=0 and k=0 since 
! there is no dependence on the nucleon momentum of this operator
!  <ms3p,mt3p | munucleon*sigma_x | ms3,mt3 > =  delta_mt3mt3p *munucleon(mt3) for ms3.ne.ms3p
   op=0.0
!  op(ms3p,ms3,mt3)
   op(1,-1,1)=munucleon(1)
   op(-1,1,1)=munucleon(1)

   op(1,-1,-1)=munucleon(-1)
   op(-1,1,-1)=munucleon(-1)
! i.e. indeed op(ms3p,ms3=ms3p,X) = 0, as required for σx.
 
   magmom = 0._DPREAL
  
    ! test the data by calculating norm/ff
   DO rindx=1,maxrho1bindex
    CALL get1Nqnnum(rindx,ms3,mt3,bm,ms3p,mt3p,bmp,bk,k)

    IF(bk.eq.0 .and. k.eq.0 & 
   & .and. mt3.eq.mt3p ) THEN   ! sum over all bm,bmp combinations 
     magmom(1)=magmom(1)+rho1b(rindx)*op(ms3p,ms3,mt3)
    END IF
   END DO ! rindx
  
   WRITE(*,'(2(A,E15.6))')  'MAGNETIC MOMENT from inserting σx at qval = ',qval,": ",3.0/2.0*magmom(1)  
  
! now the same for sigma_y   
! the matrix element of O is nonzero only for K=0 and k=0 since 
! there is no dependence on the nucleon momentum of this operator
!  <ms3p,mt3p | munucleon*sigma_y | ms3,mt3 > =  delta_mt3mt3p *munucleon(mt3) for ms3.ne.ms3p
   op=0.0
!  op(ms3p,ms3,mt3)
   op(1,-1,1)=-(0.0,1.0)*munucleon(1)
   op(-1,1,1)=(0.0,1.0)*munucleon(1)

   op(1,-1,-1)=-(0.0,1.0)*munucleon(-1)
   op(-1,1,-1)=(0.0,1.0)*munucleon(-1)
! i.e. indeed op(ms3p,ms3=ms3p,X) = 0, as required for σx.
 
   magmom = 0._DPREAL
  
    ! test the data by calculating norm/ff
   DO rindx=1,maxrho1bindex
    CALL get1Nqnnum(rindx,ms3,mt3,bm,ms3p,mt3p,bmp,bk,k)

    IF(bk.eq.0 .and. k.eq.0 .and. mt3.eq.mt3p ) THEN   ! sum over all bm,bmp combinations 
      IF(bmp.eq.1 .and.bm.eq.1) THEN
       magmom(1)=magmom(1)+rho1b(rindx)*op(ms3p,ms3,mt3)
      ELSE IF(bmp.eq.-1 .and.bm.eq.1) THEN
       magmom(1)=magmom(1)+(0.0,1.0)*rho1b(rindx)*op(ms3p,ms3,mt3)
      ELSE IF(bmp.eq.1 .and.bm.eq.-1) THEN
       magmom(1)=magmom(1)-(0.0,1.0)*rho1b(rindx)*op(ms3p,ms3,mt3)
      ELSE IF(bmp.eq.-1 .and.bm.eq.-1) THEN
       magmom(1)=magmom(1)+rho1b(rindx)*op(ms3p,ms3,mt3)
      END IF
    END IF
   END DO ! rindx
  
   WRITE(*,'(2(A,E15.6))')  'MAGNETIC MOMENT from inserting σy at qval = ',qval,": ",-3.0/2.0*magmom(1)  
   
  END IF ! master 
  
 END DO  ! while not eof, or do loop

END PROGRAM testcompdens
