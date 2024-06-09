!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!
#include "fdefs.h"
!     
! File:   CompDens.F90
! Author: andreasnogga
!
! Created on 20. März 2017, 13s:07
!

!> module CompDens provides routines to read and write the density 
!! and the corresponding bookkeeping. 
MODULE CompDens
  USE precision
  USE pwbook
  USE meshpoints
  USE mpi_const
  USE constants 
  USE spline
#ifdef RHODISTR 
  USE HDF5
  USE parallel
  USE amplitudes
  USE hdf_tool
#endif
  
  IMPLICIT NONE
  PRIVATE
  !> number of 2N channels including mt12 and m12 dependence.
  INTEGER,PUBLIC ::  num2Nchan  
 
  !> rho defines array for the density:
  !! rho(p12p,p12,chanindex=1,...,#maxrhoindex)
  REAL(spreal),ALLOCATABLE,PUBLIC :: rho(:,:,:)
  REAL(spreal),ALLOCATABLE,PUBLIC :: rho_tmp(:,:,:)
  REAL(spreal),ALLOCATABLE:: rho_p12p(:,:,:)
  REAL(spreal),ALLOCATABLE,PUBLIC :: rhoDensity(:,:,:)
! REAL(spreal),ALLOCATABLE,PUBLIC :: p12Read(:)

  !>  array to keep bookkeeping for 2N channels (including Jz of nucleus). <br>
  !! qn2Nchan(1:7,alpha)  contains l12,s12,j12,t12,mt12,m12,M depending on alpha where <br>
  !! alpha = 1, ..., #num2Nchan 
  INTEGER,ALLOCATABLE :: qn2Nchan(:,:) 

  
  !> keeps the number of coupled 1N channels 
  INTEGER,PUBLIC ::  maxrho1bindex
  !> rho1bindx(1:8,rindx) contains the set of quantum numbers 2m_s,2m_t,2M,2m_s',2m_t',2M',K,k 
  !! for a specific rindx combination. The maximal value 
  !! is #maxrho1bindex.
  INTEGER,ALLOCATABLE,PUBLIC :: rho1bindx(:,:)

  !> rho defines the density:
  !! rho(rindx=1,...,#maxrho1bindex)
  REAL(spreal),ALLOCATABLE,PUBLIC :: rho1b(:) 
  
  !> number of P12 grid points for rho
  INTEGER P12N_density
  !> P12 grid points for rho: P12P_density(1:#p12n_density)
  REAL(dpreal),ALLOCATABLE,PUBLIC :: P12P_density(:)
  !> P12 weights for rho: P12W_density(1:#p12n_density)
  REAL(dpreal),ALLOCATABLE,PUBLIC :: P12W_density(:)
  
  !> keeps the number of coupled 2N channels 
  INTEGER,PUBLIC ::  maxrhoindex
  !> rhoindx(alphap,alpha) contains the index for a specific alphap,alpha combination. The maximal value
  !! is #maxrhoindex.
  INTEGER,ALLOCATABLE,PUBLIC :: rhoindx(:,:)

  

  !> number of omega (incoming photon energy) points 
  INTEGER,PUBLIC :: numomega  
  !> number of scattering angles  
  INTEGER,PUBLIC :: numtheta
  
  !> omega(1:#numomega) keep photon energies.
  !! we assume the photon-nucleus CM
  REAL(dpreal),ALLOCATABLE :: omega(:)
  !> theta(1:#numtheta) keep photon scattering angles.
  !! we assume the photon-nucleus CM.  
  REAL(dpreal),ALLOCATABLE :: theta(:)
  !> actual photon energy.   
  REAL(dpreal),PUBLIC :: omval
  !> actual photon scattering angle 
  REAL(dpreal),PUBLIC :: thetaval
  !> actual magnitude of momentum transfer
  REAL(dpreal),PUBLIC :: qval
  !> actual angle of momentum transfer
  REAL(dpreal),PUBLIC :: thetaqval  ! omega,theta and Q and thetaQ for specific data set
   
  ! define namelists for output of data 
  ! first dimensions 
  NAMELIST /book2Ndim/ num2Nchan,maxrhoindex
  ! then data set for bookkeeping
  NAMELIST /book2Ndata/ qn2Nchan,rhoindx
  ! then data for rho 
  NAMELIST /rhodim/ P12N_density
  NAMELIST /p12mesh/ p12p_density,p12w_density
  NAMELIST /rhodata/ omval,thetaval,qval,thetaqval,rho 
 
  ! define namelists for output of data 
  ! first dimensions 
  NAMELIST /book1Ndim/ maxrho1bindex
  ! then data set for bookkeeping
  NAMELIST /book1Ndata/ rho1bindx
  NAMELIST /rho1Nnml/  rho1b
  
  ! then data for rho 
   
  !> maximal j12 use in rho
  INTEGER,PUBLIC ::  j12max_rho
  !> twice minimal Mtot 
  INTEGER,PUBLIC :: Mmin
  !> twice maximal Mtot 
  INTEGER,PUBLIC :: Mmax
  
  !> array keeping the channel numbers corresponding to specific quantum 
  !! numbers. <br>
  !! alpha2Nind(j12_m12=1..j12max_rho**2+2*j12max_rho+1,l12-j12=-1..1,s12=0..1,mt12=-1..1,Mtot=0..(Mmax-Mmin)/2) <br>
  !! the combined j12_m12 index is given by j12_m12=j12**2+j12+1+m12
  INTEGER,ALLOCATABLE :: alpha2Nind(:,:,:,:,:)
  
  INTERFACE openread_rhofile
   MODULE PROCEDURE openread_rhofile_meshdef,openread_rhofile_womeshdef
  END INTERFACE openread_rhofile

  PUBLIC generaterhochan,printrhochan,getalpha2N
  PUBLIC openwrite_rhofile,openread_rhofile
  PUBLIC write_rhoset,read_rhoset
  PUBLIC get2Nchannum

  PUBLIC openwrite_rho1bfile,openread_rho1bfile
  PUBLIC write_rho1bset,read_rho1bset 
  PUBLIC get1Nqnnum
  PUBLIC writemeshpara
! PUBLIC rhoRead
  PUBLIC printrhoread
     
#ifdef RHODISTR  
  PUBLIC openwrite_rho_hdffile,openread_rho_hdffile
  PUBLIC writehdf_rhoset,readhdf_rhoset,readsets 
  
  INTERFACE openread_rho_hdffile
   MODULE PROCEDURE openread_rho_hdffile_meshdef,openread_rho_hdffile_womeshdef
  END INTERFACE openread_rho_hdffile    
  
  !> id of the HDF5 file for rho, used externally for closeing the file. 
  INTEGER(HID_T),PUBLIC :: rhofileid
#endif  
   
  LOGICAL,PUBLIC ::  rho_eof=.false.
  LOGICAL,PUBLIC ::  rho1b_eof=.false.
CONTAINS

!> subroutine opens a file for writing and write grid and 
!! bookkeeping information to file 
!! @param filename filename of rho data file
!! @pre 
!! #generaterhochan needs to be called before
!! @post
!! file is prepared for the writing to rho data sets <br>
!! unit 10 is used for rho file 
SUBROUTINE openwrite_rhofile(filename)
 IMPLICIT NONE
 CHARACTER(LEN=*) filename 
  ! open output file for writing 
  ! and write channel bookkeeping and p12 mesh 

  IF(master) THEN
   OPEN(10,FILE=trim(filename),STATUS='NEW')
  
   ! then replace density grid points by grid points 
   ! used for the calculations which should be in 
   ! P12P and P12W
   P12N_density=P12N
   IF(allocated(P12P_density)) DEALLOCATE(P12P_density)
   IF(allocated(P12W_density)) DEALLOCATE(P12W_density)
   ALLOCATE(P12P_density(P12N_density))
   ALLOCATE(P12W_density(P12N_density))
   P12P_density=P12P
   P12W_density=P12W
   
   ! write the namelists for bookkeeping and grid points to file 
   WRITE(10,nml=book2Ndim)
   ! write bookkeeping and rho without namelist for now 
   !WRITE(10,nml=book2Ndata)
   WRITE(10,'(7I5)') qn2Nchan
   WRITE(10,'(8I6)') rhoindx
   WRITE(10,nml=rhodim)
   WRITE(10,nml=p12mesh)
  END IF
  
END SUBROUTINE openwrite_rhofile

!> subroutine opens a file for writing and writes  
!! bookkeeping information to the file 
!! @param filename filename of rho data file
!! @pre 
!! #generaterhochan needs to be called before
!! @post
!! file is prepared for  writing of rho data sets <br>
!! unit 20 is used for the rho1b file 
SUBROUTINE openwrite_rho1bfile(filename)
 IMPLICIT NONE
 CHARACTER(LEN=*) filename 
  ! open output file for writing 
  ! and write channel bookkeeping and p12 mesh 
  IF(master) THEN
   OPEN(20,FILE=trim(filename),STATUS='NEW')
     
   ! write the namelists for bookkeeping and grid points to file 
   WRITE(20,nml=book1Ndim)
   WRITE(20,nml=book1Ndata)
  END IF
  
END SUBROUTINE openwrite_rho1bfile

#ifdef RHODISTR
!> subroutine opens a file for writing and write grid and 
!! bookkeeping information to an HDF file 
!! @param filename filename of rho data file
!! @pre 
!! #generaterhochan needs to be called before
!! @post
!! file is prepared for the writing to rho data sets 
SUBROUTINE openwrite_rho_hdffile(filename)
 IMPLICIT NONE
 CHARACTER(LEN=*) filename 
  ! open output file for writing 
  ! and write channel bookkeeping and p12 mesh 
 
  CALL open_write_h5file(filename,commall,rhofileid)
  ! write mesh points directly to root of file  
  CALL writep12p(rhofileid)
  ! write number of channels and rho index directly to file 
  CALL write_scalar_int(rhofileid,'num2Nchan',num2Nchan,master)
  CALL write_scalar_int(rhofileid,'maxrhoindex',maxrhoindex,master)
  ! finally write bookkeeping arrays to file 
  CALL write_2D_int(rhofileid,'qn2Nchan',7,num2Nchan,qn2Nchan,master)
  CALL write_2D_int(rhofileid,'rhoindx',num2Nchan,num2Nchan,rhoindx,master)
  
END SUBROUTINE openwrite_rho_hdffile
#endif

!> subroutine opens a file for reading and reads grid and 
!! bookkeeping information from file 
!! @param[in] filename filename of rho data file
!! @param[in] np12 number of p12 grid points to be used 
!! @param[in] p12grid p12grid(1:#np12) given grid points  
!! @post
!! file is opened and 
!! mesh points and bookkeeping for rho data is read in  <br>
!! unit 10 is used for rho file 
!! the P12 grids points are replace by the given ones 

SUBROUTINE openread_rhofile_meshdef(filename,p12grid,np12)
 IMPLICIT NONE
 REAL(dpreal) p12grid(np12)
 INTEGER np12
 CHARACTER(LEN=*) filename
 
 ! remove old grid points 
 IF(allocated(P12P)) DEALLOCATE(P12P)
 IF(allocated(P12W)) DEALLOCATE(P12W)
 ! use dimensions and new grid points 
 P12N=NP12
 ALLOCATE(P12P(P12N),P12W(P12N))
 P12P=p12grid
 P12W=0.0_dpreal  ! do not define grid points anymore 
 
 ! now read file as usual
 CALL openread_rhofile_womeshdef(trim(filename))
 
END SUBROUTINE openread_rhofile_meshdef

#ifdef RHODISTR
!> subroutine opens a file for reading and reads grid and 
!! bookkeeping information from an HDF file 
!! @param[in] filename filename of rho data file
!! @param[in] np12 number of p12 grid points to be used 
!! @param[in] p12grid p12grid(1:#np12) given grid points  
!! @post
!! file is opened and 
!! mesh points and bookkeeping for rho data is read in  <br>
!! unit 10 is used for rho file 
!! the P12 grids points are replace by the given ones 

SUBROUTINE openread_rho_hdffile_meshdef(filename,p12grid,np12)
 IMPLICIT NONE
 REAL(dpreal) p12grid(np12)
 INTEGER np12
 CHARACTER(LEN=*) filename
 
 ! remove old grid points 
 IF(allocated(P12P)) DEALLOCATE(P12P)
 IF(allocated(P12W)) DEALLOCATE(P12W)
 ! use dimensions and new grid points 
 P12N=NP12
 ALLOCATE(P12P(P12N),P12W(P12N))
 write(*,*) "In CompDens.F90: p12grid=",p12grid 
 P12P=p12grid
 P12W=0.0_dpreal  ! do not define grid points anymore 
 
 ! now read file as usual
 CALL openread_rho_hdffile_womeshdef(trim(filename))
 
END SUBROUTINE openread_rho_hdffile_meshdef
#endif

!> subroutine opens a file for reading and reads grid and 
!! bookkeeping information from file 
!! @param[in] filename filename of rho data file
!! @post
!! file is opened and 
!! mesh points and bookkeeping for rho data is read in  <br>
!! unit 10 is used for rho file
SUBROUTINE openread_rhofile_womeshdef(filename)
 IMPLICIT NONE
 INTEGER alpha,l12,s12,j12,t12,mt12,m12,mjtot,j12_m12
 CHARACTER(LEN=*) filename 
  ! open output file for writing 
  ! and write channel bookkeeping and p12 mesh 
  write(*,*) "in CompDens.F90 in openread_rhofile_womeshdef"
  OPEN(10,FILE=trim(filename),STATUS='OLD')
  
  READ(10,nml=book2Ndim)
  IF(allocated(qn2Nchan)) DEALLOCATE(qn2Nchan)
  IF(allocated(rhoindx)) DEALLOCATE(rhoindx)
  ALLOCATE(qn2Nchan(7,num2Nchan))
  ALLOCATE(rhoindx(num2Nchan,num2Nchan))
  ! read bookkeeping data without namelist for now
  ! READ(10,nml=book2Ndata)
  READ(10,'(7I5)') qn2Nchan
  READ(10,'(8I6)') rhoindx
  
  READ(10,nml=rhodim)
  IF(allocated(P12P_density)) DEALLOCATE(P12P_density)
  IF(allocated(P12W_density)) DEALLOCATE(P12W_density)
  ALLOCATE(P12P_density(P12N_density),P12W_density(P12N_density))
  
  READ(10,nml=p12mesh)
  
 ! extract limits of quantum numbers
  Mmin=+320000  ! use really large limits 
  Mmax=-320000
  j12max_rho=0
  DO alpha=1,num2Nchan
    CALL getalpha2N(alpha,l12,s12,j12,t12,mt12,m12,mjtot)
    j12max_rho=max(j12max_rho,j12)
    Mmax=max(Mmax,mjtot)
    Mmin=min(Mmin,mjtot) 
  END DO ! alpha2N
  
  
  ! generate array that allows to get alpha channels for a given set 
  ! of l12p,s12p,j12p,mt12p,m12p,mjtotp quantum numbers
  IF(allocated(alpha2Nind)) DEALLOCATE(alpha2Nind)
  ALLOCATE(alpha2Nind(j12max_rho**2+2*j12max_rho+1,-1:1,0:1,-1:1,0:(Mmax-Mmin)/2))
  alpha2Nind=0
  DO alpha=1,num2Nchan
    CALL getalpha2N(alpha,l12,s12,j12,t12,mt12,m12,mjtot)
    j12_m12=j12**2+j12+1+m12
    alpha2Nind(j12_m12,l12-j12,s12,mt12,(mjtot-Mmin)/2)=alpha
  END DO ! alpha2N  
  
  IF(master) THEN
   CALL printrhochan
  END IF  
  
  rho_eof=.false.
END SUBROUTINE openread_rhofile_womeshdef

!> subroutine opens a 1body rho file for reading and reads 
!! bookkeeping information from file 
!! @param[in] filename filename of rho data file
!! @post
!! file is opened and 
!! bookkeeping for rho1b data is read in  <br>
!! unit 20 is used for rho file 

SUBROUTINE openread_rho1bfile(filename)
 IMPLICIT NONE
 CHARACTER(LEN=*) filename
  ! open output file for writing 
  ! and write channel bookkeeping and p12 mesh 

  OPEN(20,FILE=trim(filename),STATUS='OLD')
  
  READ(20,nml=book1Ndim)
  IF(allocated(rho1bindx)) DEALLOCATE(rho1bindx)
  ALLOCATE(rho1bindx(8,maxrho1bindex+1))
  ! read bookkeeping data without namelist for now
  ! READ(10,nml=book2Ndata)
  READ(20,nml=book1Ndata)
      
  IF(master) THEN
   CALL printrhochan
  END IF  
  
  rho1b_eof=.false. 
END SUBROUTINE openread_rho1bfile

#ifdef RHODISTR 
!> subroutine opens a file for reading and reads grid and 
!! bookkeeping information from an HDF file 
!! @param[in] filename filename of rho data file
!! @post
!! file is opened and 
!! mesh points and bookkeeping for rho data is read in  <br>
!! unit 10 is used for rho file
SUBROUTINE openread_rho_hdffile_womeshdef(filename)
 IMPLICIT NONE
 INTEGER alpha,l12,s12,j12,t12,mt12,m12,mjtot,j12_m12
 CHARACTER(LEN=*) filename 
 
  ! open input file for reading 
  ! and read channel bookkeeping and p12 mesh 
 
  CALL open_read_h5file(filename,commall,rhofileid)
  ! read mesh points directly from root of file
  IF(allocated(P12P_density)) DEALLOCATE(P12P_density)
  IF(allocated(P12W_density)) DEALLOCATE(P12W_density)
    
  CALL readp12p(rhofileid,P12P_density,P12W_density,P12N_density)

  ! read number of channels and rho index directly from file 
  CALL read_scalar_int(rhofileid,'num2Nchan',num2Nchan,commall)
  CALL read_scalar_int(rhofileid,'maxrhoindex',maxrhoindex,commall)
  ! bookkeeping arrays can be allocated 
  IF(allocated(qn2Nchan)) DEALLOCATE(qn2Nchan)
  ALLOCATE(qn2Nchan(7,num2Nchan))
  IF(allocated(rhoindx)) DEALLOCATE(rhoindx)
  ALLOCATE(rhoindx(num2Nchan,num2Nchan))
  ! and read in   
  CALL read_2D_int(rhofileid,'qn2Nchan',7,num2Nchan,qn2Nchan,commall)
  CALL read_2D_int(rhofileid,'rhoindx',num2Nchan,num2Nchan,rhoindx,commall)
  
 ! extract limits of quantum numbers
  Mmin=+320000  ! use really large limits 
  Mmax=-320000
  j12max_rho=0
  DO alpha=1,num2Nchan
    CALL getalpha2N(alpha,l12,s12,j12,t12,mt12,m12,mjtot)
    j12max_rho=max(j12max_rho,j12)
    Mmax=max(Mmax,mjtot)
    Mmin=min(Mmin,mjtot) 
  END DO ! alpha2N
  
  
  ! generate array that allows to get alpha channels for a given set 
  ! of l12p,s12p,j12p,mt12p,m12p,mjtotp quantum numbers
  IF(allocated(alpha2Nind)) DEALLOCATE(alpha2Nind)
  ALLOCATE(alpha2Nind(j12max_rho**2+2*j12max_rho+1,-1:1,0:1,-1:1,0:(Mmax-Mmin)/2))
  alpha2Nind=0
  DO alpha=1,num2Nchan
    CALL getalpha2N(alpha,l12,s12,j12,t12,mt12,m12,mjtot)
    j12_m12=j12**2+j12+1+m12
    alpha2Nind(j12_m12,l12-j12,s12,mt12,(mjtot-Mmin)/2)=alpha
  END DO ! alpha2N  
  
  IF(master) THEN
   CALL printrhochan
  END IF  
  
  rho_eof=.false.
END SUBROUTINE openread_rho_hdffile_womeshdef
#endif

!> function gets the alpha 2N channel number 
!! for a given set of quantum numbers. It is assumed that the quantum 
!! numbers are within the allowed range of values. 
!! @param l12  orbital angular momentum
!! @param s12  spin 
!! @param j12  total angular momentum
!! @param mt12 third component of total angular momentum 
!! @param mjtot third component of total angular momentum of the nucleus
!! @return channel number or zero when channel does not exist in 
!!   rho data

 FUNCTION get2Nchannum(l12,s12,j12,mt12,m12,mjtot)
  IMPLICIT NONE
  INTEGER get2Nchannum,l12,s12,j12,m12,mt12,mjtot,j12_m12
  
  IF(j12.LE.j12max_rho .AND. mjtot.GE.Mmin .AND. mjtot.LE.Mmax) THEN
   j12_m12=j12**2+j12+1+m12
   get2Nchannum=alpha2Nind(j12_m12,l12-j12,s12,mt12,(mjtot-Mmin)/2)
  ELSE
   get2Nchannum=0
  END IF
 END FUNCTION get2Nchannum

!> subroutine writes a set of rho data to the file 
!! open before using #openwrite_rhofile
!! @pre 
!! file has to be opened before
!! @param Anucl mass number for nucleus: A=3 or 4
SUBROUTINE write_rhoset(Anucl) 
  IMPLICIT NONE
  INTEGER Anucl
  INTEGER ip12,alpha2N,l12,s12,j12,t12,mt12,m12,mjtot,rindx,ip12loc
  REAL(dpreal) sum,sumn,sump,sumnn,sumgl,sumgln,sumglp,sumglnn
  REAL(spreal), ALLOCATABLE :: rhocollect(:,:,:) 
  INTEGER nump12,shiftp12
  INTEGER ierr
  
#ifdef RHODISTR  
  ! first collect the rho, if the code runs in distributed mode 
  ! finally collect rho from result
  ! collect all p12 points   
  ALLOCATE(rhocollect(P12N,P12N,maxrhoindex))
  CALL collect_part_spreal(rho,rhocollect,P12N,P12N,maxrhoindex,commp12)

#endif
  
  ! write the rho data as a new data set to the output file 
  IF(master) THEN
   ! write rhodata without namelist for now 
   !WRITE(10,nml=rhodata)
   WRITE(10,'(4E20.10)') omval,thetaval,qval,thetaqval

#ifdef RHODISTR  
  ! first collect the rho, if the code runs in distributed mode 
  ! finally collect rho from result
  ! collect all p12 points   
  WRITE(10,'(6E20.10)') rhocollect 
  DEALLOCATE(rhocollect)
#else
  WRITE(10,'(6E20.10)') rho  
#endif  
  ENDIF
  
#ifdef RHODISTR  
  nump12=mynp12
  shiftp12=myp12-1
#else
  nump12=P12N
  shiftp12=0
#endif  
  
  ! check normalization
  sum=0.0_dpreal
  sump=0.0_dpreal
  sumn=0.0_dpreal
  sumnn=0.0_dpreal
  
  DO alpha2N=1,num2Nchan
   CALL getalpha2N(alpha2N,l12,s12,j12,t12,mt12,m12,mjtot)
   rindx=rhoindx(alpha2N,alpha2N)
   IF(rindx.NE.0) THEN
    IF(mt12.EQ.-1) THEN ! nn should not exist
     DO ip12loc=1,nump12
      ip12=ip12loc+shiftp12
      sumnn=sumnn+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
      sum=sum+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
     END DO
    ELSE IF(mt12.EQ.0) THEN ! np -> spectator is p 
     DO ip12loc=1,nump12 
      ip12=ip12loc+shiftp12
      sump=sump+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
      sum=sum+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
     END DO
    ELSE IF(mt12.EQ.1) THEN ! pp -> spectator is n 
     DO ip12loc=1,nump12 
      ip12=ip12loc+shiftp12
      sumn=sumn+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
      sum=sum+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
     END DO          
    END IF  ! mt12
   END IF ! rindx
  END DO ! alpha 2N 
  
#ifdef RHODISTR
  CALL MPI_ALLREDUCE(sum,sumgl,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sumn,sumgln,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sumnn,sumglnn,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sump,sumglp,1,MPI_REAL8,MPI_SUM,commp12,ierr)     
#else
  sumgl=sum
  sumgln=sumn
  sumglp=sump
  sumglnn=sumnn
#endif  
 
  IF(master) THEN   
   IF(Anucl.EQ.3) THEN  ! A=3 includes sum over M=-1/2, +1/2 
    WRITE(*,'(A,5E15.6)')  'FF: ',qval,sumglnn/2.0_dpreal, &
    &   sumglp/2.0_dpreal,sumgln/2.0_dpreal,sumgl/2.0_dpreal
   ELSE
    WRITE(*,'(A,5E15.6)')  'FF: ',qval,sumglnn,sumglp,sumgln,sumgl
   END IF
  END IF
    
END SUBROUTINE write_rhoset

!> subroutine writes a set of rho1b data to the file 
!! open before using #openwrite_rho1bfile
!! @pre 
!! file has to be opened before
!! @param Anucl mass number for nucleus: A=3 or 4
SUBROUTINE write_rho1bset(Anucl) 
  IMPLICIT NONE
  INTEGER Anucl
  INTEGER rindx,ms3,ms3p,bm,bmp,mt3,mt3p,k,bk
  REAL(dpreal) sum,sumn,sump
    
  ! write the rho data as a new data set to the output file 
  IF(master) THEN
   WRITE(20,'(4E20.10)') omval,thetaval,qval,thetaqval
   WRITE(20,'(6E20.10)') rho1b  
  ENDIF
  
  
  ! check normalization
  sum=0.0_dpreal
  sump=0.0_dpreal
  sumn=0.0_dpreal
  
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
  
 
  IF(master) THEN   
   IF(Anucl.EQ.3) THEN  ! A=3 includes sum over M=-1/2, +1/2 
    WRITE(*,'(A,5E15.6)')  'FF1B: ',qval,0.5*sumn,0.5*sump,0.5*sum
   ELSE
    WRITE(*,'(A,5E15.6)')  'FF1B: ',qval,sumn,sump,sum
   END IF
  END IF
    
END SUBROUTINE write_rho1bset

#ifdef RHODISTR 
!> subroutine writes a set of rho data to an HDF file 
!! open before using #openwrite_rho_hdffile
!! @pre 
!! file has to be opened before
!! @param Anucl mass number for nucleus: A=3 or 4
SUBROUTINE writehdf_rhoset(Anucl) 
  IMPLICIT NONE
  INTEGER Anucl
  INTEGER ip12,alpha2N,l12,s12,j12,t12,mt12,m12,mjtot,rindx,ip12loc
  REAL(dpreal) sum,sumn,sump,sumnn,sumgl,sumgln,sumglp,sumglnn
  REAL(spreal), ALLOCATABLE :: rhocollect(:,:,:) 
  CHARACTER(LEN=50) :: blockname
  INTEGER(HID_T) group_rho_id
  INTEGER ierr
  INTEGER(HID_T) dset_rho_id,dsp_rho_id,dsp_hyper_id
  INTEGER(HID_T) msp_rho_id,msp_hyper_id
  INTEGER(HSIZE_T) dims(3)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
  INTEGER(HSIZE_T) start(3)     ! storage for offset in output arrays
  INTEGER(HSIZE_T) block(3)     ! storage for blocksize in output arrays
  INTEGER(HSIZE_T) stride(3)    ! storage for stride in output arrays
  INTEGER(HSIZE_T) count(3)     ! storage for count in output arrays
  INTEGER myrhoindx,myerhoindx,mynrhoindx     ! distribution of indx for writing  
  INTEGER myp12p_rho,myep12p_rho,mynp12p_rho  ! distribution of p12p for writing
  INTEGER(lint) blocksize
  
  ! first create new group for given omval,thetaval
  ! generate group name out of omval and thetaval 
  
  CALL set_rho_name(omval,thetaval,blockname)
  CALL h5gcreate_f(rhofileid,trim(blockname),group_rho_id, ierr)
  
  ! then write parameters of rho
  
  CALL write_scalar_dpreal(group_rho_id,'omval',omval,master)
  CALL write_scalar_dpreal(group_rho_id,'thetaval',thetaval,master)
  CALL write_scalar_dpreal(group_rho_id,'qval',qval,master)
  CALL write_scalar_dpreal(group_rho_id,'thetaqval',thetaqval,master)

  ! prepare writing of rho set 
  
  CALL distr_block(P12N,npe_p3,myp3id,myp12p_rho,myep12p_rho,mynp12p_rho) 
  CALL distr_block(maxrhoindex,npe_alpha,myalphaid,myrhoindx,myerhoindx,mynrhoindx) 

  ! now I can finally write the block 
  hdftime=hdftime-MPI_WTIME() 
      ! this needs to be written in parallel since the data 
      ! is distributed and because of performance

      !     need to create a data set in the file 
      !     with complete number amp-elements 
      !     outline shape of global set of data 
  dims(1)=P12N       
  dims(2)=P12N
  dims(3)=maxrhoindex          !  write needs dimensions of arrays 
     !     create corresponding data space            
  CALL h5screate_simple_f(3,dims,dsp_rho_id,ierr);

  !     describe the memory layout of the data
  dims(1)=P12N       
  dims(2)=mynp12
  dims(3)=maxrhoindex          !  write needs dimensions of arrays 
     !     create corresponding mem space    
  CALL h5screate_simple_f(3,dims,msp_rho_id,ierr);

     ! select the hyperslap in the file  that corresponds to the local 
     ! 3N data 
  start(1)=myp12p_rho-1
  start(2)=myp12-1
  start(3)=myrhoindx-1
  block(1)=mynp12p_rho
  block(2)=mynp12
  block(3)=mynrhoindx
  stride(1)=1
  stride(2)=1
  stride(3)=1
  count(1)=1
  count(2)=1
  count(3)=1
  
#ifdef HDF_MPIO
      ! check for consistency with MPI
  blocksize=4_lint*mynrhoindx*mynp12*mynp12p_rho   
  IF(blocksize.GE.blocklimit) THEN
    WRITE(*,*) 'block too large to be handled with MPI-IO'
    WRITE(*,*) 'blocksize = ',blocksize
    WRITE(*,*) 'blocklimit = ',blocklimit
    CALL abort
  END IF
#endif


  CALL h5scopy_f(dsp_rho_id,dsp_hyper_id,ierr)
  CALL h5sselect_hyperslab_f(dsp_hyper_id,H5S_SELECT_SET_F,&
    &                      start, count,ierr,stride,block)

     ! select the hyperslap in memorythat corresponds to the local 
     ! 3N data 
  start(1)=myp12p_rho-1
  start(2)=0
  start(3)=myrhoindx-1
  block(1)=mynp12p_rho
  block(2)=mynp12
  block(3)=mynrhoindx
  stride(1)=1
  stride(2)=1
  stride(3)=1
  count(1)=1
  count(2)=1
  count(3)=1
  
  CALL h5scopy_f(msp_rho_id,msp_hyper_id,ierr)
  CALL h5sselect_hyperslab_f(msp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)


   !    create the dataset in file 
  CALL h5dcreate_f(group_rho_id,"RHO",H5T_NATIVE_REAL,&
     &             dsp_rho_id,dset_rho_id,ierr)

   ! write data of amplitude (use here collective communication)       
  CALL h5dwrite_f(dset_rho_id, H5T_NATIVE_REAL,& 
        &      rho,dims,ierr,&
        &      mem_space_id = msp_hyper_id,&
        &      file_space_id = dsp_hyper_id,& 
        &      xfer_prp = pcollectwrite_id)
        
        
        
      ! close data set   
  CALL h5dclose_f(dset_rho_id,ierr)
      ! do not need mem_space anymore      
  CALL h5sclose_f(msp_hyper_id,ierr)
  CALL h5sclose_f(msp_rho_id,ierr)
      ! and not data space anymore    
  CALL h5sclose_f(dsp_hyper_id,ierr)
  CALL h5sclose_f(dsp_rho_id,ierr)   
  CALL h5gclose_f(group_rho_id,ierr) 
    
  hdftime=hdftime+MPI_WTIME()
      
      
    ! check normalization
  sum=0.0_dpreal
  sump=0.0_dpreal
  sumn=0.0_dpreal
  sumnn=0.0_dpreal
  
  DO alpha2N=1,num2Nchan
   CALL getalpha2N(alpha2N,l12,s12,j12,t12,mt12,m12,mjtot)
   rindx=rhoindx(alpha2N,alpha2N)
   IF(rindx.NE.0) THEN
    IF(mt12.EQ.-1) THEN ! nn should not exist
     DO ip12loc=1,mynp12
      ip12=ip12loc+myp12-1
      sumnn=sumnn+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
      sum=sum+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
     END DO
    ELSE IF(mt12.EQ.0) THEN ! np -> spectator is p 
     DO ip12loc=1,mynp12 
      ip12=ip12loc+myp12-1
      sump=sump+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
      sum=sum+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
     END DO
    ELSE IF(mt12.EQ.1) THEN ! pp -> spectator is n 
     DO ip12loc=1,mynp12 
      ip12=ip12loc+myp12-1
      sumn=sumn+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
      sum=sum+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
     END DO          
    END IF  ! mt12
   END IF ! rindx
  END DO ! alpha 2N 
  

  CALL MPI_ALLREDUCE(sum,sumgl,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sumn,sumgln,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sumnn,sumglnn,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sump,sumglp,1,MPI_REAL8,MPI_SUM,commp12,ierr)     
  
 
  IF(master) THEN   
   IF(Anucl.EQ.3) THEN  ! A=3 includes sum over M=-1/2, +1/2 
    WRITE(*,'(A,5E15.6)')  'FF: ',qval,sumglnn/2.0_dpreal, &
    &   sumglp/2.0_dpreal,sumgln/2.0_dpreal,sumgl/2.0_dpreal
   ELSE
    WRITE(*,'(A,5E15.6)')  'FF: ',qval,sumglnn,sumglp,sumgln,sumgl
   END IF
  END IF
    
END SUBROUTINE writehdf_rhoset


 ! subroutine produces strings from omval and thetaval  
 ! to get blockname and name of the data set for hdf 
SUBROUTINE set_rho_name(omval,thetaval,blockname)
 IMPLICIT NONE
 REAL(dpreal) omval,thetaval
 CHARACTER(LEN=*) blockname
   
 WRITE(blockname,'(A,SP,E13.6,A,E13.6)') & 
      &        'RHO_om=',omval,'_th=',thetaval
       
END SUBROUTINE set_rho_name
 

#endif  


!> subroutine reads in a rho data set from a file 
!! opened before with #openread_rhofile.
!! This routine replaces rho by a not distributed rho independent of the mode 
!! the code was compiled with. 
!! @param Anucl mass number for nucleus: A=3 or 4

SUBROUTINE read_rhoset(Anucl)
  IMPLICIT NONE
  INTEGER Anucl
  INTEGER ip12,alpha2N,l12,s12,j12,t12,mt12,m12,mjtot,rindx
  REAL(dpreal) sum,sumn,sump,sumnn
  REAL(spreal), ALLOCATABLE :: rho_tmp(:,:,:)
  REAL(dpreal), ALLOCATABLE :: splp12(:,:)
  INTEGER, ALLOCATABLE      :: indxp12(:,:)
  INTEGER ind1,ind2,ind3,ind4,ind1p,ind2p,ind3p,ind4p
  INTEGER jp12
  REAL(spreal) spl11,spl21,spl31,spl41,spl12,spl22,spl32,spl42,&
            &  spl13,spl23,spl33,spl43,spl14,spl24,spl34,spl44
  
  ! read the rho data as a new data set from the output file 
  IF(allocated(rho)) DEALLOCATE(rho)
  ALLOCATE(rho(P12N_density,P12N_density,maxrhoindex))
  
  ! read rhodata without namelist for now
  !READ(10,END=100,nml=rhodata)
  READ(10,'(4E20.10)',end=100) omval,thetaval,qval,thetaqval
  READ(10,'(6E20.10)',end=100) rho
   
  ! now generate intermediate array to store rho at original grid points 
  ALLOCATE(rho_tmp(P12N_density,P12N_density,maxrhoindex))
  rho_tmp=rho 
  
    
  ! the allocate rho for new grid points
  DEALLOCATE(rho)
  ALLOCATE(rho(P12N,P12N,maxrhoindex))
  
  ! get splines for interpolating from original to 
  ! new grid points 
  ALLOCATE(splp12(4,P12N),indxp12(4,P12N))
  CALL cubfast_dp(P12P_density,P12N_density,p12p,P12N,splp12,indxp12)
  
  rho=0.0_spreal 
   
  DO ip12=1,P12N
   DO jp12=1,P12N
    spl11=splp12(1,jp12)*splp12(1,ip12)
    spl12=splp12(1,jp12)*splp12(2,ip12)
    spl13=splp12(1,jp12)*splp12(3,ip12)
    spl14=splp12(1,jp12)*splp12(4,ip12)

    spl21=splp12(2,jp12)*splp12(1,ip12)
    spl22=splp12(2,jp12)*splp12(2,ip12)
    spl23=splp12(2,jp12)*splp12(3,ip12)
    spl24=splp12(2,jp12)*splp12(4,ip12)

    spl31=splp12(3,jp12)*splp12(1,ip12)
    spl32=splp12(3,jp12)*splp12(2,ip12)
    spl33=splp12(3,jp12)*splp12(3,ip12)
    spl34=splp12(3,jp12)*splp12(4,ip12)

    spl41=splp12(4,jp12)*splp12(1,ip12)
    spl42=splp12(4,jp12)*splp12(2,ip12)
    spl43=splp12(4,jp12)*splp12(3,ip12)
    spl44=splp12(4,jp12)*splp12(4,ip12)

    ind1=indxp12(1,ip12)
    ind2=indxp12(2,ip12)
    ind3=indxp12(3,ip12)
    ind4=indxp12(4,ip12)
    
    ind1p=indxp12(1,jp12)
    ind2p=indxp12(2,jp12)
    ind3p=indxp12(3,jp12)
    ind4p=indxp12(4,jp12)
    
    
    DO rindx=1,maxrhoindex
     rho(jp12,ip12,rindx) &
       & = spl11*rho_tmp(ind1p,ind1,rindx) &
       &  +spl21*rho_tmp(ind2p,ind1,rindx) &
       &  +spl31*rho_tmp(ind3p,ind1,rindx) &
       &  +spl41*rho_tmp(ind4p,ind1,rindx) &
       &  +spl12*rho_tmp(ind1p,ind2,rindx) &
       &  +spl22*rho_tmp(ind2p,ind2,rindx) &
       &  +spl32*rho_tmp(ind3p,ind2,rindx) &
       &  +spl42*rho_tmp(ind4p,ind2,rindx) &
       &  +spl13*rho_tmp(ind1p,ind3,rindx) &
       &  +spl23*rho_tmp(ind2p,ind3,rindx) &
       &  +spl33*rho_tmp(ind3p,ind3,rindx) &
       &  +spl43*rho_tmp(ind4p,ind3,rindx) &
       &  +spl14*rho_tmp(ind1p,ind4,rindx) &
       &  +spl24*rho_tmp(ind2p,ind4,rindx) &
       &  +spl34*rho_tmp(ind3p,ind4,rindx) &
       &  +spl44*rho_tmp(ind4p,ind4,rindx) 
    END DO ! rindx
   END DO ! jp12
  END DO ! ip12
  
! DEALLOCATE(rho_tmp,splp12,indxp12)
  DEALLOCATE(splp12,indxp12)
  ! test the data by calculating norm/ff
  IF(master) THEN 
   ! check normalization
   sum=0.0_dpreal
   sump=0.0_dpreal
   sumn=0.0_dpreal
   sumnn=0.0_dpreal
   
   DO alpha2N=1,num2Nchan
    CALL getalpha2N(alpha2N,l12,s12,j12,t12,mt12,m12,mjtot)
    rindx=rhoindx(alpha2N,alpha2N)
    IF(rindx.NE.0) THEN
     IF(mt12.EQ.-1) THEN ! nn should not exist
      DO ip12=1,P12N 
       sumnn=sumnn+rho(ip12,ip12,rindx)*P12P(ip12)**2*P12W(ip12) 
       sum=sum+rho(ip12,ip12,rindx)*P12P(ip12)**2*P12W(ip12) 
      END DO
     ELSE IF(mt12.EQ.0) THEN ! np -> spectator is p 
      DO ip12=1,P12N 
       sump=sump+rho(ip12,ip12,rindx)*P12P(ip12)**2*P12W(ip12) 
       sum=sum+rho(ip12,ip12,rindx)*P12P(ip12)**2*P12W(ip12) 
      END DO
     ELSE IF(mt12.EQ.1) THEN ! pp -> spectator is n 
      DO ip12=1,P12N 
       sumn=sumn+rho(ip12,ip12,rindx)*P12P(ip12)**2*P12W(ip12) 
       sum=sum+rho(ip12,ip12,rindx)*P12P(ip12)**2*P12W(ip12) 
      END DO          
     END IF  ! mt12
    END IF ! rindx
   END DO ! alpha 2N 
    
   IF(Anucl.EQ.3) THEN  ! A=3 includes sum over M=-1/2, +1/2 
    WRITE(*,'(A,5E15.6)')  'FF: ',qval,sumnn/2.0_dpreal, &
    &   sump/2.0_dpreal,sumn/2.0_dpreal,sum/2.0_dpreal
   ELSE
    WRITE(*,'(A,5E15.6)')  'FF: ',qval,sumnn,sump,sumn,sum
   END IF ! Anucl 
   
   ! also write kinetatics 
      
   WRITE(*,'(A)') 'Calculation of rho for ...'
   WRITE(*,'(A,E15.6)') '  omega     = ',omval*hbarc 
   WRITE(*,'(A,E15.6)') '  qval      = ',qval*hbarc
   WRITE(*,'(A,E15.6)') '  thetaval  = ',thetaval*180.0_dpreal/pi
   WRITE(*,'(A,E15.6)') '  thetaqval = ',thetaqval*180.0_dpreal/pi
  END IF ! master 
  RETURN
100 CONTINUE
  WRITE(*,*) 'end reached'
  rho_eof=.true.
END SUBROUTINE read_rhoset   

!> subroutine reads in a rho data set from a file 
!! opened before with #openread_rho1bfile.
!! This routine replaces rho1b by a not distributed rho independent of the mode 
!! the code was compiled with. 
!! @param Anucl mass number for nucleus: A=3 or 4
SUBROUTINE read_rho1bset(Anucl)
  IMPLICIT NONE
  INTEGER Anucl
  INTEGER rindx,ms3,ms3p,bm,bmp,mt3,mt3p,k,bk
  REAL(dpreal) sum,sumn,sump
  
  ! read the rho data as a new data set from the output file 
  IF(allocated(rho1b)) DEALLOCATE(rho1b)
  ALLOCATE(rho1b(maxrho1bindex))
  
  ! read rhodata without namelist for now
  !READ(10,END=100,nml=rhodata)
  READ(20,'(4E20.10)',end=100) omval,thetaval,qval,thetaqval
  READ(20,'(6E20.10)',end=100) rho1b
   
  IF(master) THEN   

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
    WRITE(*,'(A,5E15.6)')  'FF1B: ',qval,0.5*sumn,0.5*sump,0.5*sum
   ELSE
    WRITE(*,'(A,5E15.6)')  'FF1B: ',qval,sumn,sump,sum
   END IF
   ! also write kinetatics 
      
   WRITE(*,'(A)') 'Calculation of rho for ...'
   WRITE(*,'(A,E15.6)') '  omega     = ',omval*hbarc 
   WRITE(*,'(A,E15.6)') '  qval      = ',qval*hbarc
   WRITE(*,'(A,E15.6)') '  thetaval  = ',thetaval*180.0_dpreal/pi
   WRITE(*,'(A,E15.6)') '  thetaqval = ',thetaqval*180.0_dpreal/pi
  END IF ! master 
  RETURN
100 CONTINUE
  WRITE(*,*) 'end reached'
  rho1b_eof=.true.
END SUBROUTINE read_rho1bset   

#ifdef RHODISTR

!> subroutine reads array om and theta of all datasets in hdf file 
!! opened before with #openread_rho_hdffile.
!! @param omsets omsets(nsets) omega values used 
!! @param thetasets thetasets(nsets) omega values used 
!! @param nsets number of data sets in file 
SUBROUTINE readsets(omsets,thetasets,nsets)
 INTEGER nsets
 REAL(dpreal),ALLOCATABLE :: omsets(:),thetasets(:)
 INTEGER nmembers,ntype,ierr
 INTEGER i,pos
 REAL(dpreal) om,theta
 CHARACTER(LEN=100) :: name_buffer
 CHARACTER(LEN=20) :: dummy1,dummy2
 ! get the number of elements in root group
 CALL h5gn_members_f(rhofileid, "/", nmembers, ierr)
 
 nsets=0
 DO i=0,nmembers-1
  CALL h5gget_obj_info_idx_f(rhofileid, "/",i,name_buffer,ntype,ierr)
  
  pos=index(trim(name_buffer),"RHO")
  IF(pos.NE.0) THEN
    nsets=nsets+1
    READ(name_buffer,'(A7,SP,E13.6,A4,E13.6)') dummy1,om,dummy2,theta
  END IF
 END DO
 
 IF(ALLOCATED(omsets)) DEALLOCATE(omsets)
 IF(ALLOCATED(thetasets)) DEALLOCATE(thetasets)
 ALLOCATE(omsets(nsets),thetasets(nsets))
 nsets=0
 DO i=0,nmembers-1
  CALL h5gget_obj_info_idx_f(rhofileid, "/",i,name_buffer,ntype,ierr)
  
  pos=index(trim(name_buffer),"RHO")
  IF(pos.NE.0) THEN
    nsets=nsets+1
    READ(name_buffer,'(A7,SP,E13.6,A4,E13.6)') dummy1,om,dummy2,theta
    omsets(nsets)=om
    thetasets(nsets)=theta
  END IF
 END DO
 
 
END SUBROUTINE readsets

!> subroutine reads in a rho data set from an hdf file 
!! opened before with #openread_rho_hdffile.
!! and stores is in rho. The parameter set is identified 
!! by omset and thetaset. The actual omval and thetaval is read in 
!! from the file, too. 
!! @param Anucl mass number for nucleus: A=3 or 4
!! @param omset omega value of data set 
!! @param thetaset thetaval of data set 
SUBROUTINE readhdf_rhoset(Anucl,omset,thetaset)!rhoRead,p12Read)
  IMPLICIT NONE
  INTEGER Anucl
  REAL(dpreal) omset,thetaset
  INTEGER ip12,alpha2N,l12,s12,j12,t12,mt12,m12,mjtot,rindx
  REAL(dpreal) sumall,sumn,sump,sumnn,sumgl,sumgln,sumglp,sumglnn
  REAL(spreal), ALLOCATABLE :: rho_tmp(:,:,:)!,rho_p12p(:,:,:)!,rhoRead(:,:,:), p12Read(:)
  REAL(dpreal), ALLOCATABLE :: splp12(:,:)
  INTEGER, ALLOCATABLE      :: indxp12(:,:)
  INTEGER ind1,ind2,ind3,ind4
  INTEGER jp12,ip12loc
  REAL(spreal) spl11,spl12,spl13,spl14
  CHARACTER(LEN=100) blockname
  INTEGER(HID_T) group_rho_id
  INTEGER ierr
  INTEGER(HID_T) dset_rho_id,dsp_rho_id,dsp_hyper_id
  INTEGER(HID_T) msp_rho_id,msp_hyper_id
  INTEGER(HSIZE_T) dims(3)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
  INTEGER(HSIZE_T) start(3)     ! storage for offset in output arrays
  INTEGER(HSIZE_T) block(3)     ! storage for blocksize in output arrays
  INTEGER(HSIZE_T) stride(3)    ! storage for stride in output arrays
  INTEGER(HSIZE_T) count(3)     ! storage for count in output arrays
  INTEGER myrhoindx,myerhoindx,mynrhoindx     ! distribution of indx for reading  
  INTEGER myp12p_rho,myep12p_rho,mynp12p_rho  ! distribution of p12p for reading
  INTEGER myp12_rho,myep12_rho,mynp12_rho  ! distribution of p12p for reading
  INTEGER(lint) blocksize
  REAL(spreal) res          
            
  ! first open new group for given omval,thetaval
  ! generate group name out of omval and thetaval 
  
  CALL set_rho_name(omset,thetaset,blockname)
  CALL h5gopen_f(rhofileid,trim(blockname),group_rho_id, ierr)
  
  ! then read parameters of rho 
  
  CALL read_scalar_dpreal(group_rho_id,'omval',omval,commall)
  CALL read_scalar_dpreal(group_rho_id,'thetaval',thetaval,commall)
  CALL read_scalar_dpreal(group_rho_id,'qval',qval,commall)
  CALL read_scalar_dpreal(group_rho_id,'thetaqval',thetaqval,commall)
             
  ! prepare reading of rho set (distribution of input data)
  CALL distr_block(P12N_density,npe_p3,myp3id,myp12p_rho,myep12p_rho,mynp12p_rho) 
  CALL distr_block(P12N_density,npe_p12,myp12id,myp12_rho,myep12_rho,mynp12_rho) 
  CALL distr_block(maxrhoindex,npe_alpha,myalphaid,myrhoindx,myerhoindx,mynrhoindx) 

  ! first prepare local tmp. array 
  ALLOCATE(rho_tmp(mynp12p_rho,mynp12_rho,mynrhoindx))
! write(*,*) "(mynp12p_rho,mynp12_rho,mynrhoindx)=",mynp12p_rho,mynp12_rho,mynrhoindx 
! ALLOCATE(rhoRead(mynp12p_rho,mynp12_rho,mynrhoindx))
! ALLOCATE(p12Read(P12N_density))
  
  ! now I can finally write the block 
  hdftime=hdftime-MPI_WTIME() 
      ! this needs to be written in parallel since the data 
      ! is distributed and because of performance

  CALL h5dopen_f(group_rho_id,"RHO",dset_rho_id,ierr)
  
      ! get the data space descriptor 
  CALL h5dget_space_f(dset_rho_id, dsp_rho_id, ierr)
  
  
      !     need to create a data set in the file 
      !     with complete number amp-elements 
      !     outline shape of global set of data 

  !     describe the memory layout of the data
  dims(1)=mynp12p_rho       
  dims(2)=mynp12_rho
  dims(3)=mynrhoindx          !  write needs dimensions of arrays 
     !     create corresponding mem space    
  CALL h5screate_simple_f(3,dims,msp_rho_id,ierr);

     ! select the hyperslap in the file  that corresponds to the local 
     ! 3N data 
  start(1)=myp12p_rho-1
  start(2)=myp12_rho-1
  start(3)=myrhoindx-1
  block(1)=mynp12p_rho
  block(2)=mynp12_rho
  block(3)=mynrhoindx
  stride(1)=1
  stride(2)=1
  stride(3)=1
  count(1)=1
  count(2)=1
  count(3)=1
  
#ifdef HDF_MPIO
      ! check for consistency with MPI
  blocksize=4_lint*mynrhoindx*mynp12_rho*mynp12p_rho   
  IF(blocksize.GE.blocklimit) THEN
    WRITE(*,*) 'block too large to be handled with MPI-IO'
    WRITE(*,*) 'blocksize = ',blocksize
    WRITE(*,*) 'blocklimit = ',blocklimit
    CALL abort
  END IF
#endif

  
  CALL h5scopy_f(dsp_rho_id,dsp_hyper_id,ierr)
  CALL h5sselect_hyperslab_f(dsp_hyper_id,H5S_SELECT_SET_F,&
    &                      start, count,ierr,stride,block)

     ! select the hyperslap in memorythat corresponds to the local 
     ! 3N data 
  start(1)=0
  start(2)=0
  start(3)=0
  block(1)=mynp12p_rho
  block(2)=mynp12_rho
  block(3)=mynrhoindx
  stride(1)=1
  stride(2)=1
  stride(3)=1
  count(1)=1
  count(2)=1
  count(3)=1
  
  CALL h5scopy_f(msp_rho_id,msp_hyper_id,ierr)
  CALL h5sselect_hyperslab_f(msp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)

  CALL h5dread_f(dset_rho_id, H5T_NATIVE_REAL,& 
     &      rho_tmp,dims,ierr,&
     &      mem_space_id = msp_hyper_id,&
     &      file_space_id = dsp_hyper_id,& 
     &      xfer_prp = pcollectread_id) 
       
      ! close data set   
  CALL h5dclose_f(dset_rho_id,ierr)
      ! do not need mem_space anymore      
  CALL h5sclose_f(msp_hyper_id,ierr)
  CALL h5sclose_f(msp_rho_id,ierr)
      ! and not data space anymore    
  CALL h5sclose_f(dsp_hyper_id,ierr)
  CALL h5sclose_f(dsp_rho_id,ierr)   
  CALL h5gclose_f(group_rho_id,ierr) 
    
  hdftime=hdftime+MPI_WTIME()
  !mark 
  ! write(*,*) ""
  ! write(*,*) "in readhdf_rhoset"         
  ! write(*,*) "In CompDens.F90: shape(rhoRead)=",shape(rhoRead) 
  ! data set read now prepare interpolation 
  
  ! get splines for interpolating from original to 
  ! new grid points 

  ALLOCATE(splp12(4,P12N),indxp12(4,P12N))
  CALL cubfast_dp(P12P_density,P12N_density,p12p,P12N,splp12,indxp12)
  
  ! collect first p12p points 
  ALLOCATE(rho_p12p(mynp12p_rho,P12N_density,mynrhoindx))
  CALL collect_part_spreal(rho_tmp,rho_p12p,mynp12p_rho,P12N_density,mynrhoindx,commp12)
  DEALLOCATE(rho_tmp)
  ALLOCATE(rho_tmp(mynp12p_rho,mynp12,mynrhoindx))


  
  DO ip12loc=1,mynp12
   ip12=ip12loc+myp12-1
   spl11=splp12(1,ip12)
   spl12=splp12(2,ip12)
   spl13=splp12(3,ip12)
   spl14=splp12(4,ip12)
   ind1=indxp12(1,ip12)
   ind2=indxp12(2,ip12)
   ind3=indxp12(3,ip12)
   ind4=indxp12(4,ip12)
   
   DO rindx=1,mynrhoindx
    DO jp12=1,mynp12p_rho
     rho_tmp(jp12,ip12loc,rindx) &
       & = spl11*rho_p12p(jp12,ind1,rindx) &
       &  +spl12*rho_p12p(jp12,ind2,rindx) &
       &  +spl13*rho_p12p(jp12,ind3,rindx) &
       &  +spl14*rho_p12p(jp12,ind4,rindx)
    END DO
   END DO
  END DO ! i12loc
! write(*,*) "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
! write(*,*) "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
! write(*,*) "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
! write(*,*) "Before interpolation"
! write(*,*) "In CompDens.F90: shape(rho_p12p)=",shape(rho_p12p) 
! write(*,*) "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
! write(*,*) "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
  rhoDensity=rho_p12p
  DEALLOCATE(rho_p12p)
    ! data set read now prepare interpolation 

  
  
  ! collect first p12p points 
  ALLOCATE(rho_p12p(P12N_density,mynp12,mynrhoindx))
  CALL collect_part_spreal(rho_tmp,rho_p12p,1,P12N_density,mynp12*mynrhoindx,commp3)
  DEALLOCATE(rho_tmp)
  ALLOCATE(rho_tmp(P12N,mynp12,mynrhoindx))

  
  DO jp12=1,P12N
   spl11=splp12(1,jp12)
   spl12=splp12(2,jp12)
   spl13=splp12(3,jp12)
   spl14=splp12(4,jp12)
   ind1=indxp12(1,jp12)
   ind2=indxp12(2,jp12)
   ind3=indxp12(3,jp12)
   ind4=indxp12(4,jp12)
   
   DO rindx=1,mynrhoindx
    DO ip12loc=1,mynp12
     rho_tmp(jp12,ip12loc,rindx) &
       & = spl11*rho_p12p(ind1,ip12loc,rindx) &
       &  +spl12*rho_p12p(ind2,ip12loc,rindx) &
       &  +spl13*rho_p12p(ind3,ip12loc,rindx) &
       &  +spl14*rho_p12p(ind4,ip12loc,rindx)
    END DO
   END DO
  END DO ! i12loc
  DEALLOCATE(rho_p12p)
    ! data set read now prepare interpolation 
    
  ! finally collect all channels 

   
  ! read the rho data as a new data set from the output file 
  IF(allocated(rho)) DEALLOCATE(rho)
  ALLOCATE(rho(P12N,mynp12,maxrhoindex))
  
  CALL collect_part_spreal(rho_tmp,rho,P12N*mynp12,maxrhoindex,1,commalpha)   
  DEALLOCATE(rho_tmp,splp12,indxp12)
  

  
  ! check normalization
  sumall=0.0_dpreal
  sump=0.0_dpreal
  sumn=0.0_dpreal
  sumnn=0.0_dpreal
  DO alpha2N=1,num2Nchan
   CALL getalpha2N(alpha2N,l12,s12,j12,t12,mt12,m12,mjtot)
   rindx=rhoindx(alpha2N,alpha2N)
   IF(rindx.NE.0) THEN
    IF(mt12.EQ.-1) THEN ! nn should not exist
     DO ip12loc=1,mynp12
      ip12=ip12loc+myp12-1
      sumnn=sumnn+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
      sumall=sumall+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
     END DO
    ELSE IF(mt12.EQ.0) THEN ! np -> spectator is p 
     DO ip12loc=1,mynp12 
      ip12=ip12loc+myp12-1
      sump=sump+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
      sumall=sumall+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
     END DO
    ELSE IF(mt12.EQ.1) THEN ! pp -> spectator is n 
     DO ip12loc=1,mynp12 
      ip12=ip12loc+myp12-1
      sumn=sumn+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
      sumall=sumall+rho(ip12,ip12loc,rindx)*P12P(ip12)**2*P12W(ip12) 
     END DO          
    END IF  ! mt12
   END IF ! rindx
  END DO ! alpha 2N 
  
  CALL MPI_ALLREDUCE(sumall,sumgl,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sumn,sumgln,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sumnn,sumglnn,1,MPI_REAL8,MPI_SUM,commp12,ierr) 
  CALL MPI_ALLREDUCE(sump,sumglp,1,MPI_REAL8,MPI_SUM,commp12,ierr)       
   
  IF(master) THEN  
   IF(Anucl.EQ.3) THEN  ! A=3 includes sum over M=-1/2, +1/2 
    WRITE(*,'(A,5E15.6)')  'FF: ',qval,sumglnn/2.0_dpreal, &
    &   sumglp/2.0_dpreal,sumgln/2.0_dpreal,sumgl/2.0_dpreal
   ELSE
    WRITE(*,'(A,5E15.6)')  'FF: ',qval,sumglnn,sumglp,sumgln,sumgl
   END IF ! Anucl 

   
   ! also write kinetatics 
      
   WRITE(*,'(A)') 'Calculation of rho for ...'
   WRITE(*,'(A,E15.6)') '  omega     = ',omval*hbarc 
   WRITE(*,'(A,E15.6)') '  qval      = ',qval*hbarc
   WRITE(*,'(A,E15.6)') '  thetaval  = ',thetaval*180.0_dpreal/pi
   WRITE(*,'(A,E15.6)') '  thetaqval = ',thetaqval*180.0_dpreal/pi
  END IF ! master 
END SUBROUTINE readhdf_rhoset   

#endif

!> subroutine generates the set of 2N channels from pwbook 
!! adding m12 and M (or dependening on nucleus and j3max or j4max) 
!! also coupled channels bookkeeping is generated. Finally, a set of 
!! spectator channels is generated. 
!! @post 
!! 2N bookkeeping is stored in \ref qn2Nchan, and \ref num2Nchan <br>
!! coupled channel bookkeeping is stored in \ref rhoindx, and \ref maxrhoindx <br> 
!! 1body coupled channel bookkeeping is stored in \ref rho1bindx, and \ref maxrho1bindx <br> 
!! @param Anucl   Number of nucleons A=3 or 4 
!! @param bmtau twice third component of isospin
!! @param bjfix twice total angular momentum 
SUBROUTINE generaterhochan(Anucl,bmtau,bjfix)
 IMPLICIT NONE
 INTEGER Anucl,bmtau,bjfix
 INTEGER alpha,l12,s12,j12,t12,mt12
 INTEGER alpha2N,m12,MJtot,mcase,index
 INTEGER alphap,l12p,s12p,j12p,t12p,mt12p,m12p,mjtotp 
 INTEGER l3,I3,j3,tau3,mtau3,alpha2b,pari
 INTEGER l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34
 INTEGER j12_m12
 INTEGER k,bk,ms3,ms3p,mt3,mt3p
 IF(Anucl.EQ.3) THEN ! generate 3N case  
! first generate alpha2N bookkkeeping including m12,mt12,M
! start with counting    
 
   alpha2N=0
   DO alpha=1,alphaNNcdepmax 
    call getNNqn(alpha,l12,s12,j12,t12,mt12)
       ! all quantum numbers are actual values  
    DO m12=-j12,j12 
     DO MJtot=-bjfix,bjfix,2   ! bjfix is half integer 
       alpha2N=alpha2N+1 
     END DO ! MJtot
    END DO  ! m12 
   END DO   ! alpha 
   num2Nchan=alpha2N
! then set quantum numbers 

   ALLOCATE(qn2Nchan(7,num2Nchan))
   alpha2N=0
   j12max_rho=0
   Mmax=-bjfix
   Mmin=bjfix
   DO alpha=1,alphaNNcdepmax 
    call getNNqn(alpha,l12,s12,j12,t12,mt12)
       ! all quantum numbers are actual values  
    DO m12=-j12,j12            ! m12 is integer not multiplied by 2
     DO MJtot=-bjfix,bjfix,2   ! MJtot,j3max is half integer multiplied by 2 
       alpha2N=alpha2N+1 
       qn2Nchan(1,alpha2N)=l12
       qn2Nchan(2,alpha2N)=s12
       qn2Nchan(3,alpha2N)=j12
       qn2Nchan(4,alpha2N)=t12
       qn2Nchan(5,alpha2N)=mt12
       qn2Nchan(6,alpha2N)=m12
       qn2Nchan(7,alpha2N)=MJtot
       j12max_rho=max(j12max_rho,j12)
       Mmax=max(Mmax,mjtot)
       Mmin=min(Mmin,mjtot) 
     END DO ! MJtot
    END DO  ! m12 
   END DO   ! alpha    
   IF(alpha2N.NE.num2Nchan) STOP 'counting wrong'
 ELSE   ! assuming 4N system  
  
   
! first generate alpha2N bookkkeeping including m12,mt12,M
! start with counting    
 
   alpha2N=0
   DO alpha=1,alphaNNcdepmax 
    call getNNqn(alpha,l12,s12,j12,t12,mt12)
       ! all quantum numbers are actual values  
    DO m12=-j12,j12 
     DO MJtot=-bjfix,bjfix,2   ! j4max is integer nevertheless MJtot is twice its value
       alpha2N=alpha2N+1 
     END DO ! MJtot
    END DO  ! m12 
   END DO   ! alpha 
   num2Nchan=alpha2N
! then set quantum numbers 

   ALLOCATE(qn2Nchan(7,num2Nchan))
   alpha2N=0
   j12max_rho=0
   Mmax=-bjfix
   Mmin=bjfix
   DO alpha=1,alphaNNcdepmax 
    call getNNqn(alpha,l12,s12,j12,t12,mt12)
       ! all quantum numbers are actual values  
    DO m12=-j12,j12            ! m12 is integer not multiplied by 2
     DO MJtot=-bjfix,bjfix,2   ! j4max is integer, MJtot is nevertheless multiplied by 2 
       alpha2N=alpha2N+1 
       qn2Nchan(1,alpha2N)=l12
       qn2Nchan(2,alpha2N)=s12
       qn2Nchan(3,alpha2N)=j12
       qn2Nchan(4,alpha2N)=t12
       qn2Nchan(5,alpha2N)=mt12
       qn2Nchan(6,alpha2N)=m12
       qn2Nchan(7,alpha2N)=MJtot
       j12max_rho=max(j12max_rho,j12)
       Mmax=max(Mmax,mjtot)
       Mmin=min(Mmin,mjtot)
     END DO ! MJtot
    END DO  ! m12 
   END DO   ! alpha    
   IF(alpha2N.NE.num2Nchan) STOP 'counting wrong'
 END IF   
    
 ! generate index for coupled 2N channels 
 ! later this index can be used to related different mjtot to each other 
 ! this is the same for 3N and 4N 
 ! and counting and storing can be done at the same time 
 ALLOCATE(rhoindx(num2Nchan,num2Nchan))
 rhoindx=0  ! indicates not coupled 
 
 index=0
 DO alpha=1,num2Nchan
  CALL getalpha2N(alpha,l12,s12,j12,t12,mt12,m12,mjtot)
  DO alphap=1,num2Nchan
   CALL getalpha2N(alphap,l12p,s12p,j12p,t12p,mt12p,m12p,mjtotp)
   IF(mt12p.EQ.mt12) THEN ! charge conservation 
   index=index+1
   rhoindx(alphap,alpha)=index 
   END IF
  END DO ! alpha2N 
 END DO ! alpha2N
 maxrhoindex=index  

 ! generate index for coupled 1N channels 
 
 ! first count chanels
 index=0
 DO bk=0,1
  DO k=-bk,bk
   DO mjtotp=-bjfix,bjfix,2
    DO mjtot=-bjfix,bjfix,2
     DO mt3=-1,1,2 
      DO ms3p=-1,1,2
       DO ms3=-1,1,2
        index=index+1
       END DO ! ms3
      END DO ! ms3p
     END DO ! mt3=mt3p
    END DO ! mjtot
   END DO  ! mjtotp
  END DO ! k
 END DO  ! bk 
 
 maxrho1bindex=index
 
 ALLOCATE(rho1bindx(8,maxrho1bindex))
 index=0
 DO bk=0,1
  DO k=-bk,bk
   DO mjtotp=-bjfix,bjfix,2
    DO mjtot=-bjfix,bjfix,2
     DO mt3=-1,1,2 
      DO ms3p=-1,1,2
       DO ms3=-1,1,2
        index=index+1
        rho1bindx(1,index)=ms3
        rho1bindx(2,index)=mt3
        rho1bindx(3,index)=mjtot
        rho1bindx(4,index)=ms3p
        rho1bindx(5,index)=mt3  ! conservation of charge 
        rho1bindx(6,index)=mjtotp
        rho1bindx(7,index)=k
        rho1bindx(8,index)=bk
        
       END DO ! ms3
      END DO ! ms3p
     END DO ! mt3=mt3p
    END DO ! mjtot
   END DO  ! mjtotp
  END DO ! k
 END DO  ! bk 

 IF(maxrho1bindex.NE.index) STOP 'problem in generaterhochan'

  ! generate array that allows to get alpha channels for a given set 
  ! of l12p,s12p,j12p,mt12p,m12p,mjtotp quantum numbers
 
  IF(allocated(alpha2Nind)) DEALLOCATE(alpha2Nind)
  ALLOCATE(alpha2Nind(j12max_rho**2+2*j12max_rho+1,-1:1,0:1,-1:1,0:(Mmax-Mmin)/2))
  alpha2Nind=0
  DO alpha=1,num2Nchan
    CALL getalpha2N(alpha,l12,s12,j12,t12,mt12,m12,mjtot)
    j12_m12=j12**2+j12+1+m12
    alpha2Nind(j12_m12,l12-j12,s12,mt12,(mjtot-Mmin)/2)=alpha
  END DO ! alpha2N  
 
END SUBROUTINE generaterhochan
 
!> subroutine prints the set of 2N channels for rho. 
SUBROUTINE printrhochan
 IMPLICIT NONE
 INTEGER Anucl,mtau3min,mtau3max
 INTEGER alpha,l12,s12,j12,t12,mt12
 INTEGER alpha2N,m12,MJtot,mcase,index
 INTEGER alphap,l12p,s12p,j12p,t12p,mt12p,m12p,mjtotp 
 INTEGER l3,I3,j3,tau3,mtau3,alpha2b,pari
 INTEGER l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34
 INTEGER rindx,ms3,mt3,bm,ms3p,mt3p,bmp,bk,k
 
 IF(master.AND.allocated(qn2Nchan)) THEN ! now print the 2N bookkeeping for 3N or 4N 
  WRITE(*,'(A,I4)') 'j12max_rho = ',j12max_rho
  WRITE(*,'(A,I4,A,I4)') 'mjtot      = ',Mmin,"...",Mmax
 
! WRITE(*,'(A,2A5,2X,4A5,2X,3A6)') '        ','alpha','al2','l12','s12','j12','t12','mt12','m12','mjtot' 
! DO alpha=1,num2Nchan
!   CALL getalpha2N(alpha,l12,s12,j12,t12,mt12,m12,mjtot)
!   alpha2N=get2Nchannum(l12,s12,j12,mt12,m12,mjtot)
!   WRITE(*,'(A,2I5,2X,4I5,2X,3I6)') 'ALPHA2N: ',alpha,alpha2N,l12,s12,j12,t12,mt12,m12,mjtot
! END DO ! alpha2N
 END IF ! master 
 
 IF(master.AND.allocated(rho1bindx)) THEN ! now print the 1N bookkeeping for 3N or 4N 
  WRITE(*,'(A,I4)') 'j12max_rho = ',j12max_rho
 
  WRITE(*,'(A12,A6,2X,3A5,2X,3A5,2X,2A4)') '  ','rindx','ms3','mt3','bm','ms3p','mt3p','bmp','bk','k'
  DO rindx=1,maxrho1bindex
    CALL get1Nqnnum(rindx,ms3,mt3,bm,ms3p,mt3p,bmp,bk,k)

    WRITE(*,'(A12,I6,2X,3I5,2X,3I5,2X,2I4)') 'RHO1BINDX: ',rindx,ms3,mt3,bm,ms3p,mt3p,bmp,bk,k
  END DO ! rindx
 END IF ! master 
END SUBROUTINE printrhochan

!> Subroutine retrieves NN quantum from
!! the prepared table. 
!! @param[in] alpha index for m quantum number alpha = 1,...,num2Nchan
!! @param[out] l12 NN orbital angular momentum (actual value)
!! @param[out] s12 NN spin   (actual value)
!! @param[out] j12 NN total angular momentum (actual value)
!! @param[out] t12 NN isospin (actual value) 
!! @param[out] mt12 third component of NN isospin (actual value -1,..,1)
!! @param[out] m12 third component of NN angular momentum (actual value)
!! @param[out] mjtot third component of nuclear total angular momentum (twice actual value)
!! @pre 
!! subroutine #generate2nchan needs to be called before
SUBROUTINE getalpha2N(alpha,l12,s12,j12,t12,mt12,m12,mjtot)
  IMPLICIT NONE 
  INTEGER alpha,l12,s12,j12,t12,mt12,m12,mjtot
  l12=qn2Nchan(1,alpha)
  s12=qn2Nchan(2,alpha)
  j12=qn2Nchan(3,alpha)
  t12=qn2Nchan(4,alpha)
  mt12=qn2Nchan(5,alpha)
  m12=qn2Nchan(6,alpha)
  mjtot=qn2Nchan(7,alpha)
END SUBROUTINE getalpha2N 

!> Subroutine retrieves 1N quantum from
!! the prepared table. 
!! @param[in] rindx index for m quantum number rindx = 1,...,maxrho1bindex
!! @param[out] ms3 twice the third component of spin of the spectator (incoming)
!! @param[out] mt3 twice the third component of isospin of the spectator (incoming)
!! @param[out] bm  twice the third component of spin of the nucleus (incoming)
!! @param[out] ms3p twice the third component of spin of the spectator (outgoing)
!! @param[out] mt3p twice the third component of isospin of the spectator (outgoing)
!! @param[out] bmp  twice the third component of spin of the nucleus (outgoing)
!! @param[out] bk power of momentum of spectator  (0 or 1)
!! @param[out] k  third component of spherical components of momentum of spectator
!! @pre 
!! subroutine #generate2nchan needs to be called before
SUBROUTINE get1Nqnnum(rindx,ms3,mt3,bm,ms3p,mt3p,bmp,bk,k)
  IMPLICIT NONE
  INTEGER rindx,ms3,mt3,bm,ms3p,mt3p,bmp,bk,k
  ms3=rho1bindx(1,rindx)
  mt3=rho1bindx(2,rindx)
  bm=rho1bindx(3,rindx)
  ms3p=rho1bindx(4,rindx)
  mt3p=rho1bindx(5,rindx)
  bmp=rho1bindx(6,rindx)
  k=rho1bindx(7,rindx)
  bk=rho1bindx(8,rindx)
END SUBROUTINE get1Nqnnum 

!> Subroutine writes the file meshpoints.dat to the current directory 
!! based on the mesh points and weights given as parameter. 
!! @param[in] PP array with meshpoints in fm-1 
!! @param[in] PW array with corresponding weights
!! @param[in] NP number of meshpoints 
!! @post routine write a file "meshpoints.dat" to disk 

 SUBROUTINE writemeshpara(PP,PW,NP)
  IMPLICIT NONE
  INTEGER NP
  REAL(dpreal) PP(NP),PW(NP)
  INTEGER ip12 
  
  IF(master) THEN
   OPEN(45,FILE="meshpoints.dat",FORM="FORMATTED",STATUS="UNKNOWN")  
   WRITE(45,'(A)') 'FLEX'
   WRITE(45,'(3I4)') 1,1,1  ! write unnecessary angular grtid sizes 
   WRITE(45,'(A)') 'READ'   ! flexible P12 points 
   WRITE(45,'(I5)') NP     
   DO ip12=1,NP
    WRITE(45,'(2E20.10)') PP(ip12),PW(ip12)
   END DO
   WRITE(45,'(A)') 'TRNS'   ! set other points just to one small point P3
   WRITE(45,'(3E20.10,2I5)') 1.0E-10,5.0E-10,20.0E-10,1,0 
   WRITE(45,'(A)') 'TRNS'   ! set other points just to one small point Q4
   WRITE(45,'(3E20.10,2I5)') 1.0E-10,5.0E-10,20.0E-10,1,0 
   WRITE(45,'(A)') 'TRNS'   ! set other points just to one small point P34
   WRITE(45,'(3E20.10,2I5)') 1.0E-10,5.0E-10,20.0E-10,1,0 
   WRITE(45,'(A)') 'TRNS'   ! set other points just to one small point Q
   WRITE(45,'(3E20.10,2I5)') 1.0E-10,5.0E-10,20.0E-10,1,0  
   CLOSE(45)
   write(*,*) "just closed meshpoints.dat"
  END IF
  
  
 END SUBROUTINE writemeshpara

 subroutine printrhoread()
  implicit none
! write(*,*) "In CompDens.F90: rhoGrid=",rhoGrid 
 end subroutine
 
END MODULE CompDens
