#include "fdefs.h" 
MODULE mpi_const
 USE precision
#ifdef MPIMOD 
 USE mpi
#else
 INCLUDE 'mpif.h'
#endif
 
 INTERFACE collect_part
  MODULE PROCEDURE collect_part_spreal,collect_part_dpreal,&
              &     collect_part_spcmplx,collect_part_dpcmplx,&
              &     collect_part_int
 END INTERFACE
 
 INTERFACE collect_piece
  MODULE PROCEDURE collect_piece_spreal,collect_piece_dpreal,&
       &     collect_piece_spcmplx,collect_piece_dpcmplx
 END INTERFACE
 
 REAL(dpreal) :: collecttime=0.0_dpreal

 LOGICAL,PUBLIC :: master=.true.
 
 CONTAINS 
  
 !> subroutine prints out the measured memory usage  
 !! - maximum,minimum and average is printed for HEAP and STACK
 !! - used the flags #BGQ and #LINUX
 !! - on other computers it just prints zeros 
 subroutine print_mem_stat
  implicit none 
  real(dpreal) heap_size,stack_size,heap_avail,stack_avail,heap_max
  real(dpreal) mem(5),maxmem(5),avemem(5),minmem(5)
  INTEGER ierr,myid,npe 
  
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierr)

#if defined(BGQ) || defined(LINUX) || defined(MAC) 
  call get_memory_usage(heap_size,stack_size,heap_avail,&
     &                        stack_avail,heap_max)       

  mem(1)=heap_size
  mem(2)=stack_size
  mem(3)=heap_avail
  mem(4)=stack_avail
  mem(5)=heap_max  
#else
  mem(1)=0.0
  mem(2)=0.0
  mem(3)=0.0
  mem(4)=0.0
  mem(5)=0.0
#endif
  
  CALL MPI_ALLREDUCE(mem,avemem,5,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mem,maxmem,5,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(mem,minmem,5,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)


  avemem=avemem/REAL(npe,dpreal)


  IF(myid.EQ.0) THEN
   WRITE(*,*)
   WRITE(*,'(A)') 'MEMORY STATS: '
   WRITE(*,'(A,3E15.6)') 'HEAP SIZE:   ',minmem(1),avemem(1),maxmem(1)       
   WRITE(*,'(A,3E15.6)') 'STACK SIZE:  ',minmem(2),avemem(2),maxmem(2)       
   WRITE(*,'(A,3E15.6)') 'HEAP AVAIL:  ',minmem(3),avemem(3),maxmem(3)       
   WRITE(*,'(A,3E15.6)') 'STACK AVAIL: ',minmem(4),avemem(4),maxmem(4)       
   WRITE(*,'(A,3E15.6)') 'HEAP MAX:    ',minmem(5),avemem(5),maxmem(5)       
   WRITE(*,*)
  END IF
 END subroutine print_mem_stat


 
 !> subroutine prints out the measured runtime with tag 
 !! runtime is total local runtime measure,e.g., with MPI_WTIME 
 !! maximum and average over all CPU's is printed out on master
 SUBROUTINE printtime(tag,runtime)
  IMPLICIT NONE 
  REAL(dpreal) runtime,maxtime,avetime
  CHARACTER(LEN=*) tag
  INTEGER ierr,myid,npe
  
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierr)
  
  CALL MPI_ALLREDUCE(runtime,avetime,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(runtime,maxtime,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

  avetime=avetime/REAL(npe,dpreal)

  IF(myid.EQ.0) THEN
   WRITE(*,'(A,2E20.10)') tag,avetime,maxtime        
  END IF
 END SUBROUTINE printtime

 !> distribution of mesh/alpha in blocks over communicator with nsize PE 
 !! @param[in]  NP number of meshpoints/channels
 !! @param[in]  nsize size of communicator 
 !! @param[in]  id id of the task to be used  
 !! @param[out] locmyp   first point on PE id 
 !! @param[out] loclastp last point on PE id 
 !! @param[out] locmynp  number of local mesh points on PE id 
 SUBROUTINE distr_block(NP,nsize,id,locmyp,loclastp,locmynp)
  IMPLICIT NONE 
  INTEGER NP,nsize,id,locmyp,locmynp,loclastp
  INTEGER divnp,modnp

#ifdef DEBUG
  IF(id.GE.nsize) STOP 'inconsistent size' 
#endif

  divnp=NP/nsize
  modnp=mod(NP,nsize)

  IF(id.LT.modnp) THEN
   ! first modnp PE's get divnp+1 elements 
   locmynp=divnp+1
   IF(locmynp.NE.0) THEN
     locmyp=(divnp+1)*id+1
   ELSE
     locmyp=1
   END IF
   loclastp=locmyp-1+locmynp
  ELSE
   ! last nsize-modnp PE's get divnp elements 
   locmynp=divnp
   IF(locmynp.NE.0) THEN
     locmyp=divnp*id+modnp+1
   ELSE
     locmyp=1
   END IF
   loclastp=locmyp-1+locmynp
  END IF
  ! now modnp*(divnp+1) + (nsize-modnp)*dvinp = nsize*divnp+modnp=np 
  ! points are distributed 
 END SUBROUTINE distr_block
 
 !> Position of block distributed data in terms of  id and local grid point
 !! @param[in] ip  grid point
 !! @param[in] NP  total number of grid points
 !! @param[in] npe number of tasks in communicator
 !! @param[out] id  id number  
 !! @param[out] locip local grid point 
 SUBROUTINE pos_block(ip,np,npe,id,locip)
  IMPLICIT NONE
  INTEGER ip,np,npe,id,locip
  INTEGER locnp,modnp,ipbound

  locnp=NP/npe
  modnp = mod(NP,npe) 
  ipbound=(locnp+1)*modnp
  
  IF(ip.LE.ipbound) THEN 
   id=(ip-1)/(locnp+1)
   locip=ip-id*(locnp+1)
  ELSE
   id=(ip-1-modnp)/locnp
   locip=ip-modnp-id*locnp
  END IF
 END SUBROUTINE pos_block

 !> distribution of mesh/alpha piecewise over communicator 
 !! this distribution gives <br>
 !!                        element 1,nsize+1,2*nsize+1, ... to PE 0 <br>
 !!                        element 2,nsize+2,2*nsize+2, ... to PE 1 <br>
 !!     ..... <br>
 !!                        element nsize,2*nsize,3*nsize, ... to PE nsize-1 <br>
 !! in order to have altogether np elements, the first tasks 
 !! gets  1 element more than the last ones (if NP mod nsize .NE. 0) 
 !! @param[in]  NP number of meshpoints/channels
 !! @param[in]  nsize size of communicator 
 !! @param[in]  id id of the task to be used  
 !! @param[out] locmynp  number of grid points on id 
 SUBROUTINE distr_piece(NP,nsize,id,locmynp)
  IMPLICIT NONE 
  INTEGER NP,nsize,id,locmynp
  INTEGER divnp,modnp

#ifdef DEBUG
  IF(id.GE.nsize) STOP 'inconsistent size' 
#endif

  divnp=NP/nsize
  modnp=mod(NP,nsize)

  IF(id.LT.modnp) THEN
   ! first modnp PE's get divnp+1 elements 
   locmynp=divnp+1
  ELSE
   ! last nsize-modnp PE's get divnp elements 
   locmynp=divnp
  END IF
  ! now modnp*(divnp+1) + (nsize-modnp)*dvinp = nsize*divnp+modnp=np 
  ! points are distributed 
 END SUBROUTINE distr_piece
 
 !> Position of piecewise distributed data in terms of  id and local grid point
 !! @param[in] ip  grid point
 !! @param[in] np  total number of grid points
 !! @param[in] npe number of tasks in communicator
 !! @param[out] id  id number  
 !! @param[out] locip local grid point 
 SUBROUTINE pos_piece(ip,np,npe,id,locip)
  IMPLICIT NONE
  INTEGER ip,np,npe,id,locip
  INTEGER locnp,modnp,ipbound
  id=mod(ip-1,npe)
  locip=(ip-1-id)/npe+1
 END SUBROUTINE pos_piece
 

!> subroutine collects a distributed array Aloc into an array Aglob
!! it assumed that the local array has the shape 
!!     Aloc(dim1,dim2loc,dim3) 
!! where dim2loc may be different on the  processors 
!! in communicator comm according to the distr_block
!! dim1 and dim2 should be the same within communicator comm
!! the data of Aloc is collected in array
!!     Aglob(dim1,dim2tot,dim3) 
!! where dim2tot is the sum of all dim2loc in the communicator 
!! When debugging is enable, the routine checks the overall 
!! size of Aloc and Aglob. The exact shape is not tested, since 
!! it can be different to rank=3, if dim1 and dim2 are combinations of dimensions
!! implementation is for REAL(spreal) 
!! @param[in] Aloc  local data to be collected 
!! @param[out] Aglob on exit global data 
!! @param[in] dim1 leading dimension
!! @param[in] dim2tot dimension to be collected (global number)
!! @param[in] dim3  trailing dimension
!! @param[in] comm  MPI communicator  
 SUBROUTINE collect_part_spreal(Aloc,Aglob,dim1,dim2tot,dim3,comm) 
  IMPLICIT NONE
  REAL(spreal) :: Aloc(*),Aglob(*) 
  INTEGER :: dim1,dim2tot,dim3,comm
  INTEGER myid,dim2loc,dim2first,dim2last
  INTEGER nsize,ierr,id,iglob,ibuff,i3,i2
  REAL(spreal),ALLOCATABLE :: Abuff(:) 
  INTEGER,ALLOCATABLE :: nloc(:),ploc(:)
  INTEGER ncopy,nstride

  collecttime=collecttime-MPI_WTIME()
  ! first determine from communicator the local dimension of Aloc 
  CALL MPI_COMM_SIZE(comm,nsize,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(comm,myid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL distr_block(dim2tot,nsize,myid,dim2first,dim2last,dim2loc)

#ifdef DEBUG
  ! check whether sizes are consistent with dimensions 
!  IF(SIZE(Aloc).NE.dim1*dim2loc*dim3) STOP 'Aloc dimensions'
!  IF(SIZE(Aglob).NE.dim1*dim2tot*dim3) STOP 'Aglob dimensions'
#endif  

  ! the basic operation is an ALLGATHERV 
  ! to get the most efficient communication, we 
  ! want to have only one call of allgather
  ! for dim3.NE.1, one local part of Aglob, will 
  ! be spreaded over several parts in the global
  ! array. This situation cannot be handled properly 
  ! by standard MPI_ALLGATHERV calls. 
  ! therefore, for dim3.NE.1, the local data will 
  ! be gathered to a global temporary array Abuff 
  ! and then properly distributed

  ! first set up the vector of local dimensions and displacements  
  ALLOCATE(nloc(0:nsize-1),ploc(0:nsize-1))

  DO id=0,nsize-1
   CALL distr_block(dim2tot,nsize,id,ploc(id),dim2last,nloc(id))
   ploc(id)=(ploc(id)-1)*dim1*dim3
   nloc(id)=nloc(id)*dim1*dim3
  END DO
#ifdef DEBUG
  IF(master) THEN
   DO id=0,nsize-1
    WRITE(*,'(A,4I15)') 'DISTR: ',id,ploc(id),nloc(id),ploc(id)+nloc(id) 
   END DO
  END IF
#endif

  ! now do the gathering 
  ! either directly into outgoing Aglob 
  ! or first in Abuff and then copy in Aglob
  
  IF(dim3.EQ.1) THEN  
   ! simple allgatherv to Aglob 
   ! complete array is sequential ordering of parts of the PE's 
   CALL MPI_ALLGATHERV(Aloc,dim1*dim2loc*dim3,MPI_REAL4,& 
                       Aglob,nloc,ploc,MPI_REAL4,&
                       comm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
  ELSE
   ALLOCATE(Abuff(dim1*dim2tot*dim3))
   ! gather to Abuff
   CALL MPI_ALLGATHERV(Aloc,dim1*dim2loc*dim3,MPI_REAL4,& 
                       Abuff,nloc,ploc,MPI_REAL4,&
                       comm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'   
   ! Abuff is ordered as Abuff(dim1,dim2loc(id),dim3,id)
   ! Aglob is ordered as Aglob(dim1,dim2tot,dim3)
   ! reorder so that Aglob has the expected ordering 

   ibuff=1
   DO id=0,nsize-1
    ncopy=nloc(id)/dim3
    nstride=dim1*dim2tot-ncopy
    iglob=ploc(id)/dim3+1
    DO i3=1,dim3
     DO i2=1,ncopy
      Aglob(iglob)=Abuff(ibuff)
      ibuff=ibuff+1
      iglob=iglob+1
     END DO
     iglob=iglob+nstride   
    END DO
   END DO

   DEALLOCATE(Abuff)
  END IF

  DEALLOCATE(nloc,ploc)

  collecttime=collecttime+MPI_WTIME()
 END SUBROUTINE collect_part_spreal

 

!> subroutine collects a distributed array Aloc into an array Aglob
!! it assumed that the local array has the shape 
!!     Aloc(dim1,dim2loc,dim3) 
!! where dim2loc may be different on the  processors 
!! in communicator comm according to the distr_block
!! dim1 and dim2 should be the same within communicator comm
!! the data of Aloc is collected in array
!!     Aglob(dim1,dim2tot,dim3) 
!! where dim2tot is the sum of all dim2loc in the communicator 
!! When debugging is enable, the routine checks the overall 
!! size of Aloc and Aglob. The exact shape is not tested, since 
!! it can be different to rank=3, if dim1 and dim2 are combinations of dimensions
!! implementation is for INTEGER 
!! @param[in] Aloc  local data to be collected 
!! @param[out] Aglob on exit global data 
!! @param[in] dim1 leading dimension
!! @param[in] dim2tot dimension to be collected (global number)
!! @param[in] dim3  trailing dimension
!! @param[in] comm  MPI communicator
 SUBROUTINE collect_part_int(Aloc,Aglob,dim1,dim2tot,dim3,comm) 
  IMPLICIT NONE
  INTEGER :: Aloc(*),Aglob(*) 
  INTEGER :: dim1,dim2tot,dim3,comm
  INTEGER myid,dim2loc,dim2first,dim2last
  INTEGER nsize,ierr,id,iglob,ibuff,i3,i2
  INTEGER,ALLOCATABLE :: Abuff(:) 
  INTEGER,ALLOCATABLE :: nloc(:),ploc(:)

  collecttime=collecttime-MPI_WTIME()
  
  ! first determine from communicator the local dimension of Aloc 
  CALL MPI_COMM_SIZE(comm,nsize,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(comm,myid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL distr_block(dim2tot,nsize,myid,dim2first,dim2last,dim2loc)
  
#ifdef DEBUG
  ! check whether sizes are consistent with dimensions 
!  IF(SIZE(Aloc).NE.dim1*dim2loc*dim3) STOP 'Aloc dimensions'
!  IF(SIZE(Aglob).NE.dim1*dim2tot*dim3) STOP 'Aglob dimensions'
#endif  

  ! the basic operation is an ALLGATHERV 
  ! to get the most efficient communication, we 
  ! want to have only one call of allgather
  ! for dim3.NE.1, one local part of Aglob, will 
  ! be spreaded over several parts in the global
  ! array. This situation cannot be handled properly 
  ! by standard MPI_ALLGATHERV calls. 
  ! therefore, for dim3.NE.1, the local data will 
  ! be gathered to a global temporary array Abuff 
  ! and then properly distributed

  ! first set up the vector of local dimensions and displacements  
  ALLOCATE(nloc(0:nsize-1),ploc(0:nsize-1))

  DO id=0,nsize-1
   CALL distr_block(dim2tot,nsize,id,ploc(id),dim2last,nloc(id))
   ploc(id)=(ploc(id)-1)*dim1*dim3
   nloc(id)=nloc(id)*dim1*dim3
  END DO
  
#ifdef DEBUG
  IF(master) THEN
   DO id=0,nsize-1
    WRITE(*,'(A,4I15)') 'DISTR: ',id,ploc(id),nloc(id),ploc(id)+nloc(id) 
   END DO
  END IF
#endif

  ! now do the gathering 
  ! either directly into outgoing Aglob 
  ! or first in Abuff and then copy in Aglob
  
  IF(dim3.EQ.1) THEN  
   ! simple allgatherv to Aglob 
   ! complete array is sequential ordering of parts of the PE's 
   CALL MPI_ALLGATHERV(Aloc,dim1*dim2loc*dim3,MPI_INTEGER,& 
                       Aglob,nloc,ploc,MPI_INTEGER,&
                       comm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
  ELSE
   ALLOCATE(Abuff(dim1*dim2tot*dim3))
   ! gather to Abuff
   CALL MPI_ALLGATHERV(Aloc,dim1*dim2loc*dim3,MPI_INTEGER,& 
                       Abuff,nloc,ploc,MPI_INTEGER,&
                       comm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   
   ! Abuff is ordered as Abuff(dim1,dim2loc(id),dim3,id)
   ! Aglob is ordered as Aglob(dim1,dim2tot,dim3)
   ! reorder so that Aglob has the expected ordering 
   ibuff=0
   DO id=0,nsize-1
    DO i3=1,dim3
     iglob=(i3-1)*dim1*dim2tot+ploc(id)/dim3
     DO i2=1,nloc(id)/dim3
      ibuff=ibuff+1
      iglob=iglob+1
      Aglob(iglob)=Abuff(ibuff)
     END DO
    END DO
   END DO

   DEALLOCATE(Abuff)
  END IF

  DEALLOCATE(nloc,ploc)

  collecttime=collecttime+MPI_WTIME()
 END SUBROUTINE collect_part_INT

!> subroutine collects a distributed array Aloc into an array Aglob
!! it assumed that the local array has the shape 
!!     Aloc(dim1,dim2loc,dim3) 
!! where dim2loc may be different on the  processors 
!! in communicator comm according to the distr_block
!! dim1 and dim2 should be the same within communicator comm
!! the data of Aloc is collected in array
!!     Aglob(dim1,dim2tot,dim3) 
!! where dim2tot is the sum of all dim2loc in the communicator 
!! When debugging is enable, the routine checks the overall 
!! size of Aloc and Aglob. The exact shape is not tested, since 
!! it can be different to rank=3, if dim1 and dim2 are combinations of dimensions
!! implementation is for REAL(dpreal) 
!! @param[in] Aloc  local data to be collected 
!! @param[out] Aglob on exit global data 
!! @param[in] dim1 leading dimension
!! @param[in] dim2tot dimension to be collected (global number)
!! @param[in] dim3  trailing dimension
!! @param[in] comm  MPI communicator
 SUBROUTINE collect_part_dpreal(Aloc,Aglob,dim1,dim2tot,dim3,comm) 
  IMPLICIT NONE
  REAL(dpreal) :: Aloc(*),Aglob(*) 
  INTEGER :: dim1,dim2tot,dim3,comm
  INTEGER myid,dim2loc,dim2first,dim2last
  INTEGER nsize,ierr,id,iglob,ibuff,i3,i2
  REAL(dpreal),ALLOCATABLE :: Abuff(:) 
  INTEGER,ALLOCATABLE :: nloc(:),ploc(:)

  collecttime=collecttime-MPI_WTIME()
  
  ! first determine from communicator the local dimension of Aloc 
  CALL MPI_COMM_SIZE(comm,nsize,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(comm,myid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL distr_block(dim2tot,nsize,myid,dim2first,dim2last,dim2loc)
  
#ifdef DEBUG
  ! check whether sizes are consistent with dimensions 
!  IF(SIZE(Aloc).NE.dim1*dim2loc*dim3) STOP 'Aloc dimensions'
!  IF(SIZE(Aglob).NE.dim1*dim2tot*dim3) STOP 'Aglob dimensions'
#endif  

  ! the basic operation is an ALLGATHERV 
  ! to get the most efficient communication, we 
  ! want to have only one call of allgather
  ! for dim3.NE.1, one local part of Aglob, will 
  ! be spreaded over several parts in the global
  ! array. This situation cannot be handled properly 
  ! by standard MPI_ALLGATHERV calls. 
  ! therefore, for dim3.NE.1, the local data will 
  ! be gathered to a global temporary array Abuff 
  ! and then properly distributed

  ! first set up the vector of local dimensions and displacements  
  ALLOCATE(nloc(0:nsize-1),ploc(0:nsize-1))

  DO id=0,nsize-1
   CALL distr_block(dim2tot,nsize,id,ploc(id),dim2last,nloc(id))
   ploc(id)=(ploc(id)-1)*dim1*dim3
   nloc(id)=nloc(id)*dim1*dim3
  END DO
  
#ifdef DEBUG
  IF(master) THEN
   DO id=0,nsize-1
    WRITE(*,'(A,4I15)') 'DISTR: ',id,ploc(id),nloc(id),ploc(id)+nloc(id) 
   END DO
  END IF
#endif

  ! now do the gathering 
  ! either directly into outgoing Aglob 
  ! or first in Abuff and then copy in Aglob
  
  IF(dim3.EQ.1) THEN  
   ! simple allgatherv to Aglob 
   ! complete array is sequential ordering of parts of the PE's 
   CALL MPI_ALLGATHERV(Aloc,dim1*dim2loc*dim3,MPI_REAL8,& 
                       Aglob,nloc,ploc,MPI_REAL8,&
                       comm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
  ELSE
   ALLOCATE(Abuff(dim1*dim2tot*dim3))
   ! gather to Abuff
   CALL MPI_ALLGATHERV(Aloc,dim1*dim2loc*dim3,MPI_REAL8,& 
                       Abuff,nloc,ploc,MPI_REAL8,&
                       comm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   
   ! Abuff is ordered as Abuff(dim1,dim2loc(id),dim3,id)
   ! Aglob is ordered as Aglob(dim1,dim2tot,dim3)
   ! reorder so that Aglob has the expected ordering 
   ibuff=0
   DO id=0,nsize-1
    DO i3=1,dim3
     iglob=(i3-1)*dim1*dim2tot+ploc(id)/dim3
     DO i2=1,nloc(id)/dim3
      ibuff=ibuff+1
      iglob=iglob+1
      Aglob(iglob)=Abuff(ibuff)
     END DO
    END DO
   END DO

   DEALLOCATE(Abuff)
  END IF

  DEALLOCATE(nloc,ploc)

  collecttime=collecttime+MPI_WTIME()
 END SUBROUTINE collect_part_dpreal

!> subroutine collects a distributed array Aloc into an array Aglob
!! it assumed that the local array has the shape 
!!     Aloc(dim1,dim2loc,dim3) 
!! where dim2loc may be different on the  processors 
!! in communicator comm according to the distr_block
!! dim1 and dim2 should be the same within communicator comm
!! the data of Aloc is collected in array
!!     Aglob(dim1,dim2tot,dim3) 
!! where dim2tot is the sum of all dim2loc in the communicator 
!! When debugging is enable, the routine checks the overall 
!! size of Aloc and Aglob. The exact shape is not tested, since 
!! it can be different to rank=3, if dim1 and dim2 are combinations of dimensions
!! implementation is for COMPLEX(spreal) 
!! @param[in] Aloc  local data to be collected 
!! @param[out] Aglob on exit global data 
!! @param[in] dim1 leading dimension
!! @param[in] dim2tot dimension to be collected (global number)
!! @param[in] dim3  trailing dimension
!! @param[in] comm  MPI communicator
 SUBROUTINE collect_part_spcmplx(Aloc,Aglob,dim1,dim2tot,dim3,comm) 
  IMPLICIT NONE
  COMPLEX(spreal) :: Aloc(*),Aglob(*) 
  INTEGER :: dim1,dim2tot,dim3,comm
  INTEGER myid,dim2loc,dim2first,dim2last
  INTEGER nsize,ierr,id,iglob,ibuff,i3,i2
  COMPLEX(spreal),ALLOCATABLE :: Abuff(:) 
  INTEGER,ALLOCATABLE :: nloc(:),ploc(:)

  collecttime=collecttime-MPI_WTIME()
  
  ! first determine from communicator the local dimension of Aloc 
  CALL MPI_COMM_SIZE(comm,nsize,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(comm,myid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL distr_block(dim2tot,nsize,myid,dim2first,dim2last,dim2loc)
  
#ifdef DEBUG
  ! check whether sizes are consistent with dimensions 
  IF(SIZE(Aloc).NE.dim1*dim2loc*dim3) STOP 'Aloc dimensions'
  IF(SIZE(Aglob).NE.dim1*dim2tot*dim3) STOP 'Aglob dimensions'
#endif  

  ! the basic operation is an ALLGATHERV 
  ! to get the most efficient communication, we 
  ! want to have only one call of allgather
  ! for dim3.NE.1, one local part of Aglob, will 
  ! be spreaded over several parts in the global
  ! array. This situation cannot be handled properly 
  ! by standard MPI_ALLGATHERV calls. 
  ! therefore, for dim3.NE.1, the local data will 
  ! be gathered to a global temporary array Abuff 
  ! and then properly distributed

  ! first set up the vector of local dimensions and displacements  
  ALLOCATE(nloc(0:nsize-1),ploc(0:nsize-1))

  DO id=0,nsize-1
   CALL distr_block(dim2tot,nsize,id,ploc(id),dim2last,nloc(id))
   ploc(id)=(ploc(id)-1)*dim1*dim3
   nloc(id)=nloc(id)*dim1*dim3
  END DO

#ifdef DEBUG
  IF(master) THEN
   DO id=0,nsize-1
    WRITE(*,'(A,4I15)') 'DISTR: ',id,ploc(id),nloc(id),ploc(id)+nloc(id) 
   END DO
  END IF
#endif
  
  ! now do the gathering 
  ! either directly into outgoing Aglob 
  ! or first in Abuff and then copy in Aglob
  
  IF(dim3.EQ.1) THEN  
   ! simple allgatherv to Aglob 
   ! complete array is sequential ordering of parts of the PE's 
   CALL MPI_ALLGATHERV(Aloc,dim1*dim2loc*dim3,MPI_COMPLEX8,& 
                       Aglob,nloc,ploc,MPI_COMPLEX8,&
                       comm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
  ELSE
   ALLOCATE(Abuff(dim1*dim2tot*dim3))
   ! gather to Abuff
   CALL MPI_ALLGATHERV(Aloc,dim1*dim2loc*dim3,MPI_COMPLEX8,& 
                       Abuff,nloc,ploc,MPI_COMPLEX8,&
                       comm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   
   ! Abuff is ordered as Abuff(dim1,dim2loc(id),dim3,id)
   ! Aglob is ordered as Aglob(dim1,dim2tot,dim3)
   ! reorder so that Aglob has the expected ordering 
   ibuff=0
   DO id=0,nsize-1
    DO i3=1,dim3
     iglob=(i3-1)*dim1*dim2tot+ploc(id)/dim3
     DO i2=1,nloc(id)/dim3
      ibuff=ibuff+1
      iglob=iglob+1
      Aglob(iglob)=Abuff(ibuff)
     END DO
    END DO
   END DO

   DEALLOCATE(Abuff)
  END IF

  DEALLOCATE(nloc,ploc)

  collecttime=collecttime+MPI_WTIME()
 END SUBROUTINE collect_part_spcmplx

!> subroutine collects a distributed array Aloc into an array Aglob
!! it assumed that the local array has the shape 
!!     Aloc(dim1,dim2loc,dim3) 
!! where dim2loc may be different on the  processors 
!! in communicator comm according to the distr_block
!! dim1 and dim2 should be the same within communicator comm
!! the data of Aloc is collected in array
!!     Aglob(dim1,dim2tot,dim3) 
!! where dim2tot is the sum of all dim2loc in the communicator 
!! When debugging is enable, the routine checks the overall 
!! size of Aloc and Aglob. The exact shape is not tested, since 
!! it can be different to rank=3, if dim1 and dim2 are combinations of dimensions
!! implementation is for COMPLEX(dpreal) 
!! @param[in] Aloc  local data to be collected 
!! @param[out] Aglob on exit global data 
!! @param[in] dim1 leading dimension
!! @param[in] dim2tot dimension to be collected (global number)
!! @param[in] dim3  trailing dimension
!! @param[in] comm  MPI communicator
 SUBROUTINE collect_part_dpcmplx(Aloc,Aglob,dim1,dim2tot,dim3,comm) 
  IMPLICIT NONE
  COMPLEX(dpreal) :: Aloc(*),Aglob(*) 
  INTEGER :: dim1,dim2tot,dim3,comm
  INTEGER myid,dim2loc,dim2first,dim2last
  INTEGER nsize,ierr,id,iglob,ibuff,i3,i2
  COMPLEX(dpreal),ALLOCATABLE :: Abuff(:) 
  INTEGER,ALLOCATABLE :: nloc(:),ploc(:)

  collecttime=collecttime-MPI_WTIME()
  
  ! first determine from communicator the local dimension of Aloc 
  CALL MPI_COMM_SIZE(comm,nsize,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(comm,myid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL distr_block(dim2tot,nsize,myid,dim2first,dim2last,dim2loc)
  
#ifdef DEBUG
  ! check whether sizes are consistent with dimensions 
!  IF(SIZE(Aloc).NE.dim1*dim2loc*dim3) STOP 'Aloc dimensions'
!  IF(SIZE(Aglob).NE.dim1*dim2tot*dim3) STOP 'Aglob dimensions'
#endif  

  ! the basic operation is an ALLGATHERV 
  ! to get the most efficient communication, we 
  ! want to have only one call of allgather
  ! for dim3.NE.1, one local part of Aglob, will 
  ! be spreaded over several parts in the global
  ! array. This situation cannot be handled properly 
  ! by standard MPI_ALLGATHERV calls. 
  ! therefore, for dim3.NE.1, the local data will 
  ! be gathered to a global temporary array Abuff 
  ! and then properly distributed

  ! first set up the vector of local dimensions and displacements  
  ALLOCATE(nloc(0:nsize-1),ploc(0:nsize-1))

  DO id=0,nsize-1
   CALL distr_block(dim2tot,nsize,id,ploc(id),dim2last,nloc(id))
   ploc(id)=(ploc(id)-1)*dim1*dim3
   nloc(id)=nloc(id)*dim1*dim3
  END DO

#ifdef DEBUG
  IF(master) THEN
   DO id=0,nsize-1
    WRITE(*,'(A,4I15)') 'DISTR: ',id,ploc(id),nloc(id),ploc(id)+nloc(id) 
   END DO
  END IF
#endif
  
  ! now do the gathering 
  ! either directly into outgoing Aglob 
  ! or first in Abuff and then copy in Aglob
  
  IF(dim3.EQ.1) THEN  
   ! simple allgatherv to Aglob 
   ! complete array is sequential ordering of parts of the PE's 
   CALL MPI_ALLGATHERV(Aloc,dim1*dim2loc*dim3,MPI_COMPLEX16,& 
                       Aglob,nloc,ploc,MPI_COMPLEX16,&
                       comm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
  ELSE
   ALLOCATE(Abuff(dim1*dim2tot*dim3))
   ! gather to Abuff
   CALL MPI_ALLGATHERV(Aloc,dim1*dim2loc*dim3,MPI_COMPLEX16,& 
                       Abuff,nloc,ploc,MPI_COMPLEX16,&
                       comm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   
   ! Abuff is ordered as Abuff(dim1,dim2loc(id),dim3,id)
   ! Aglob is ordered as Aglob(dim1,dim2tot,dim3)
   ! reorder so that Aglob has the expected ordering 
   ibuff=0
   DO id=0,nsize-1
    DO i3=1,dim3
     iglob=(i3-1)*dim1*dim2tot+ploc(id)/dim3
     DO i2=1,nloc(id)/dim3
      ibuff=ibuff+1
      iglob=iglob+1
      Aglob(iglob)=Abuff(ibuff)
     END DO
    END DO
   END DO

   DEALLOCATE(Abuff)
  END IF

  DEALLOCATE(nloc,ploc)

  collecttime=collecttime+MPI_WTIME()
 END SUBROUTINE collect_part_dpcmplx



!> subroutine collects a piecewise distributed array Aloc into an array Aglob
!! it assumed that the local array has the shape 
!!     Aloc(dim1,dim2loc,dim3) 
!! where dim2loc may be different on the  processors 
!! in communicator comm according to the #distr_piece
!! dim1 and dim2 should be the same within communicator comm
!! the data of Aloc is collected in array
!!     Aglob(dim1,dim2tot,dim3) 
!! where dim2tot is the sum of all dim2loc in the communicator 
!! When debugging is enable, the routine checks the overall 
!! size of Aloc and Aglob. The exact shape is not tested, since 
!! it can be different to rank=3, if dim1 and dim2 are combinations of dimensions
!! implementation is for REAL(spreal) 
!! @param[in] Aloc  local data to be collected 
!! @param[out] Aglob on exit global data 
!! @param[in] dim1 leading dimension
!! @param[in] dim2tot dimension to be collected (global number)
!! @param[in] dim3  trailing dimension
!! @param[in] comm  MPI communicator
 subroutine collect_piece_spreal(aloc,aglob,dim1,dim2tot,dim3,comm)
  implicit none
  real(spreal) :: Aloc(*),Aglob(*)
  integer :: dim1,dim2tot,dim3,comm,myid,dim2loc,nsize,ierr,id,iglob,ibuff,i1,i2,i3
  real(spreal),allocatable :: Abuff(:)
  integer,allocatable :: nloc(:),nloc_tot(:),ploc(:)

  collecttime=collecttime-mpi_wtime()
  
  ! first determine from communicator the local dimension of Aloc
  call mpi_comm_size(comm,nsize,ierr)
  if(ierr.ne.0) stop 'mpi failure'
  call mpi_comm_rank(comm,myid,ierr)
  if(ierr.ne.0) stop 'mpi failure'
  call distr_piece(dim2tot,nsize,myid,dim2loc)
  
#ifdef DEBUG
  ! check whether sizes are consistent with dimensions
  ! if(SIZE(Aloc).NE.dim1*dim2loc*dim3) STOP 'Aloc dimensions'
  ! if(SIZE(Aglob).NE.dim1*dim2tot*dim3) STOP 'Aglob dimensions'
#endif

  ! the basic operation is an ALLGATHERV
  ! since the variables are distributed piece-wise,
  ! we have the same problem as in collect_part_spreal
  ! for dim3.NE.1 therefore, the local data will
  ! be gathered to a global temporary array Abuff
  ! and then properly distributed

  ! first set up the vector of local dimensions and displacements
  allocate(nloc(0:nsize-1),nloc_tot(0:nsize-1),ploc(0:nsize-1))

  ploc(0)=0
  do id=0,nsize-1
   call distr_piece(dim2tot,nsize,id,nloc(id))
   if (id.ge.1) then
    ploc(id)=ploc(id-1)+dim1*nloc(id-1)*dim3
   end if
   nloc_tot(id)=nloc(id)*dim1*dim3
  end do

#ifdef DEBUG
  IF(master) THEN
   DO id=0,nsize-1
    WRITE(*,'(A,4I15)') 'DISTR: ',id,ploc(id),nloc(id),ploc(id)+nloc(id) 
   END DO
  END IF
#endif
  
  ! now do the gathering 
  ! first in Abuff and then copy in Aglob
  
  allocate(Abuff(dim1*dim2tot*dim3))
  ! gather to Abuff
  call mpi_allgatherv(Aloc,dim1*dim2loc*dim3,MPI_REAL4,& 
       Abuff,nloc_tot,ploc,MPI_REAL4,comm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  deallocate(ploc,nloc,nloc_tot)
  
  ! Abuff is ordered as Abuff(dim1,dim2loc(id),dim3,id)
  ! Aglob is ordered as Aglob(dim1,dim2tot,dim3)
  ! reorder so that Aglob has the expected ordering

  ibuff=0
  do id=0,nsize-1
   do i3=1,dim3
    do i2=id+1,dim2tot,nsize
     iglob=(i2-1)*dim1+(i3-1)*dim1*dim2tot
     do i1=1,dim1
      ibuff=ibuff+1
      iglob=iglob+1
      Aglob(iglob)=Abuff(ibuff)
     end do
    end do
   end do
  end do

  deallocate(Abuff)

  collecttime=collecttime+MPI_WTIME()
 end subroutine collect_piece_spreal


 !> subroutine collects a piecewise distributed array Aloc into an array Aglob
!! it assumed that the local array has the shape 
!!     Aloc(dim1,dim2loc,dim3) 
!! where dim2loc may be different on the  processors 
!! in communicator comm according to the #distr_piece
!! dim1 and dim2 should be the same within communicator comm
!! the data of Aloc is collected in array
!!     Aglob(dim1,dim2tot,dim3) 
!! where dim2tot is the sum of all dim2loc in the communicator 
!! When debugging is enable, the routine checks the overall 
!! size of Aloc and Aglob. The exact shape is not tested, since 
!! it can be different to rank=3, if dim1 and dim2 are combinations of dimensions
!! implementation is for REAL(dpreal) 
!! @param[in] Aloc  local data to be collected 
!! @param[out] Aglob on exit global data 
!! @param[in] dim1 leading dimension
!! @param[in] dim2tot dimension to be collected (global number)
!! @param[in] dim3  trailing dimension
!! @param[in] comm  MPI communicator
 subroutine collect_piece_dpreal(aloc,aglob,dim1,dim2tot,dim3,comm)
  implicit none
  real(dpreal) :: Aloc(*),Aglob(*)
  integer :: dim1,dim2tot,dim3,comm,myid,dim2loc,nsize,ierr,id,iglob,ibuff,i1,i2,i3
  real(dpreal),allocatable :: Abuff(:)
  integer,allocatable :: nloc(:),nloc_tot(:),ploc(:)

  collecttime=collecttime-mpi_wtime()
  
  ! first determine from communicator the local dimension of Aloc
  call mpi_comm_size(comm,nsize,ierr)
  if(ierr.ne.0) stop 'mpi failure'
  call mpi_comm_rank(comm,myid,ierr)
  if(ierr.ne.0) stop 'mpi failure'
  call distr_piece(dim2tot,nsize,myid,dim2loc)
  
#ifdef DEBUG
  ! check whether sizes are consistent with dimensions
  ! if(SIZE(Aloc).NE.dim1*dim2loc*dim3) STOP 'Aloc dimensions'
  ! if(SIZE(Aglob).NE.dim1*dim2tot*dim3) STOP 'Aglob dimensions'
#endif

  ! the basic operation is an ALLGATHERV
  ! since the variables are distributed piece-wise,
  ! we have the same problem as in collect_part_dpreal
  ! for dim3.NE.1 therefore, the local data will
  ! be gathered to a global temporary array Abuff
  ! and then properly distributed

  ! first set up the vector of local dimensions and displacements
  allocate(nloc(0:nsize-1),nloc_tot(0:nsize-1),ploc(0:nsize-1))

  ploc(0)=0
  do id=0,nsize-1
   call distr_piece(dim2tot,nsize,id,nloc(id))
   if (id.ge.1) then
    ploc(id)=ploc(id-1)+dim1*nloc(id-1)*dim3
   end if
   nloc_tot(id)=nloc(id)*dim1*dim3
  end do

#ifdef DEBUG
  IF(master) THEN
   DO id=0,nsize-1
    WRITE(*,'(A,4I15)') 'DISTR: ',id,ploc(id),nloc(id),ploc(id)+nloc(id) 
   END DO
  END IF
#endif
  
  ! now do the gathering 
  ! first in Abuff and then copy in Aglob
  
  allocate(Abuff(dim1*dim2tot*dim3))
  ! gather to Abuff
  call mpi_allgatherv(Aloc,dim1*dim2loc*dim3,MPI_REAL8,& 
       Abuff,nloc_tot,ploc,MPI_REAL8,comm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  deallocate(ploc,nloc,nloc_tot)
  
  ! Abuff is ordered as Abuff(dim1,dim2loc(id),dim3,id)
  ! Aglob is ordered as Aglob(dim1,dim2tot,dim3)
  ! reorder so that Aglob has the expected ordering

  ibuff=0
  do id=0,nsize-1
   do i3=1,dim3
    do i2=id+1,dim2tot,nsize
     iglob=(i2-1)*dim1+(i3-1)*dim1*dim2tot
     do i1=1,dim1
      ibuff=ibuff+1
      iglob=iglob+1
      Aglob(iglob)=Abuff(ibuff)
     end do
    end do
   end do
  end do

  deallocate(Abuff)

  collecttime=collecttime+MPI_WTIME()
 end subroutine collect_piece_dpreal


 !> subroutine collects a piecewise distributed array Aloc into an array Aglob
!! it assumed that the local array has the shape 
!!     Aloc(dim1,dim2loc,dim3) 
!! where dim2loc may be different on the  processors 
!! in communicator comm according to the #distr_piece
!! dim1 and dim2 should be the same within communicator comm
!! the data of Aloc is collected in array
!!     Aglob(dim1,dim2tot,dim3) 
!! where dim2tot is the sum of all dim2loc in the communicator 
!! When debugging is enable, the routine checks the overall 
!! size of Aloc and Aglob. The exact shape is not tested, since 
!! it can be different to rank=3, if dim1 and dim2 are combinations of dimensions
!! implementation is for COMPLEX(spreal) 
!! @param[in] Aloc  local data to be collected 
!! @param[out] Aglob on exit global data 
!! @param[in] dim1 leading dimension
!! @param[in] dim2tot dimension to be collected (global number)
!! @param[in] dim3  trailing dimension
!! @param[in] comm  MPI communicator
 subroutine collect_piece_spcmplx(aloc,aglob,dim1,dim2tot,dim3,comm) 
  implicit none
  complex(spreal) :: Aloc(*),Aglob(*)
  integer :: dim1,dim2tot,dim3,comm,myid,dim2loc,nsize,ierr,id,iglob,ibuff,i1,i2,i3
  complex(spreal),allocatable :: Abuff(:) 
  integer,allocatable :: nloc(:),nloc_tot(:),ploc(:)

  collecttime=collecttime-mpi_wtime()
  
  ! first determine from communicator the local dimension of Aloc 
  call mpi_comm_size(comm,nsize,ierr)
  if(ierr.ne.0) stop 'mpi failure'
  call mpi_comm_rank(comm,myid,ierr)
  if(ierr.ne.0) stop 'mpi failure'
  call distr_piece(dim2tot,nsize,myid,dim2loc)
  
#ifdef DEBUG
  ! check whether sizes are consistent with dimensions
  ! if(SIZE(Aloc).NE.dim1*dim2loc*dim3) STOP 'Aloc dimensions'
  ! if(SIZE(Aglob).NE.dim1*dim2tot*dim3) STOP 'Aglob dimensions'
#endif

  ! the basic operation is an ALLGATHERV
  ! since the variables are distributed piece-wise,
  ! we have the same problem as in collect_part_spcmplx
  ! for dim3.NE.1 therefore, the local data will
  ! be gathered to a global temporary array Abuff
  ! and then properly distributed

  ! first set up the vector of local dimensions and displacements
  allocate(nloc(0:nsize-1),nloc_tot(0:nsize-1),ploc(0:nsize-1))

  ploc(0)=0
  do id=0,nsize-1
   call distr_piece(dim2tot,nsize,id,nloc(id))
   if (id.ge.1) then
    ploc(id)=ploc(id-1)+dim1*nloc(id-1)*dim3
   end if
   nloc_tot(id)=nloc(id)*dim1*dim3
  end do

#ifdef DEBUG
  IF(master) THEN
   DO id=0,nsize-1
    WRITE(*,'(A,4I15)') 'DISTR: ',id,ploc(id),nloc(id),ploc(id)+nloc(id) 
   END DO
  END IF
#endif
  
  ! now do the gathering 
  ! first in Abuff and then copy in Aglob
  
  allocate(Abuff(dim1*dim2tot*dim3))
  ! gather to Abuff
  call mpi_allgatherv(Aloc,dim1*dim2loc*dim3,MPI_COMPLEX8,& 
       Abuff,nloc_tot,ploc,MPI_COMPLEX8,comm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  deallocate(ploc,nloc,nloc_tot)
  
  ! Abuff is ordered as Abuff(dim1,dim2loc(id),dim3,id)
  ! Aglob is ordered as Aglob(dim1,dim2tot,dim3)
  ! reorder so that Aglob has the expected ordering

  ibuff=0
  do id=0,nsize-1
   do i3=1,dim3
    do i2=id+1,dim2tot,nsize
     iglob=(i2-1)*dim1+(i3-1)*dim1*dim2tot
     do i1=1,dim1
      ibuff=ibuff+1
      iglob=iglob+1
      Aglob(iglob)=Abuff(ibuff)
     end do
    end do
   end do
  end do

  deallocate(Abuff)

  collecttime=collecttime+MPI_WTIME()
 end subroutine collect_piece_spcmplx

!> subroutine collects a piecewise distributed array Aloc into an array Aglob
!! it assumed that the local array has the shape 
!!     Aloc(dim1,dim2loc,dim3) 
!! where dim2loc may be different on the  processors 
!! in communicator comm according to the #distr_piece
!! dim1 and dim2 should be the same within communicator comm
!! the data of Aloc is collected in array
!!     Aglob(dim1,dim2tot,dim3) 
!! where dim2tot is the sum of all dim2loc in the communicator 
!! When debugging is enable, the routine checks the overall 
!! size of Aloc and Aglob. The exact shape is not tested, since 
!! it can be different to rank=3, if dim1 and dim2 are combinations of dimensions
!! implementation is for COMPLEX(dpreal) 
!! @param[in] Aloc  local data to be collected 
!! @param[out] Aglob on exit global data 
!! @param[in] dim1 leading dimension
!! @param[in] dim2tot dimension to be collected (global number)
!! @param[in] dim3  trailing dimension
!! @param[in] comm  MPI communicator
 subroutine collect_piece_dpcmplx(aloc,aglob,dim1,dim2tot,dim3,comm) 
  implicit none
  complex(dpreal) :: Aloc(*),Aglob(*)
  integer :: dim1,dim2tot,dim3,comm,myid,dim2loc,nsize,ierr,id,iglob,ibuff,i1,i2,i3
  complex(dpreal),allocatable :: Abuff(:) 
  integer,allocatable :: nloc(:),nloc_tot(:),ploc(:)

  collecttime=collecttime-mpi_wtime()
  
  ! first determine from communicator the local dimension of Aloc 
  call mpi_comm_size(comm,nsize,ierr)
  if(ierr.ne.0) stop 'mpi failure'
  call mpi_comm_rank(comm,myid,ierr)
  if(ierr.ne.0) stop 'mpi failure'
  call distr_piece(dim2tot,nsize,myid,dim2loc)
  
#ifdef DEBUG
  ! check whether sizes are consistent with dimensions
  ! if(SIZE(Aloc).NE.dim1*dim2loc*dim3) STOP 'Aloc dimensions'
  ! if(SIZE(Aglob).NE.dim1*dim2tot*dim3) STOP 'Aglob dimensions'
#endif

  ! the basic operation is an ALLGATHERV
  ! since the variables are distributed piece-wise,
  ! we have the same problem as in collect_part_dpcmplx
  ! for dim3.NE.1 therefore, the local data will
  ! be gathered to a global temporary array Abuff
  ! and then properly distributed

  ! first set up the vector of local dimensions and displacements
  allocate(nloc(0:nsize-1),nloc_tot(0:nsize-1),ploc(0:nsize-1))

  ploc(0)=0
  do id=0,nsize-1
   call distr_piece(dim2tot,nsize,id,nloc(id))
   if (id.ge.1) then
    ploc(id)=ploc(id-1)+dim1*nloc(id-1)*dim3
   end if
   nloc_tot(id)=nloc(id)*dim1*dim3
  end do

#ifdef DEBUG
  IF(master) THEN
   DO id=0,nsize-1
    WRITE(*,'(A,4I15)') 'DISTR: ',id,ploc(id),nloc(id),ploc(id)+nloc(id) 
   END DO
  END IF
#endif
  
  ! now do the gathering 
  ! first in Abuff and then copy in Aglob
  
  allocate(Abuff(dim1*dim2tot*dim3))
  ! gather to Abuff
  call mpi_allgatherv(Aloc,dim1*dim2loc*dim3,MPI_COMPLEX16,& 
       Abuff,nloc_tot,ploc,MPI_COMPLEX16,comm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  deallocate(ploc,nloc,nloc_tot)

  ! Abuff is ordered as Abuff(dim1,dim2loc(id),dim3,id)
  ! Aglob is ordered as Aglob(dim1,dim2tot,dim3)
  ! reorder so that Aglob has the expected ordering

  ibuff=0
  do id=0,nsize-1
   do i3=1,dim3
    do i2=id+1,dim2tot,nsize
     iglob=(i2-1)*dim1+(i3-1)*dim1*dim2tot
     do i1=1,dim1
      ibuff=ibuff+1
      iglob=iglob+1
      Aglob(iglob)=Abuff(ibuff)
     end do
    end do
   end do
  end do

  deallocate(Abuff)

  collecttime=collecttime+MPI_WTIME()
 end subroutine collect_piece_dpcmplx

 
END MODULE mpi_const
