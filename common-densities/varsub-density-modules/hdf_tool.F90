#include "fdefs.h"

MODULE hdf_tool
 USE precision
 USE HDF5
 USE mpi_const
 IMPLICIT NONE 
 PRIVATE

#ifdef HDF_MPIO
 
!> defines limit of 2GB for MPI-IO, if used. 
!! not required for other drivers 
 INTEGER(lint),PUBLIC,PARAMETER :: blocklimit=2_lint**31
#endif
 
 PUBLIC writedata_dpreal,readdata_dpreal,&
   &    open_write_newh5file
 PUBLIC write_simple_dpreal,write_simple_spreal,write_simple_int
 PUBLIC write_2D_dpreal,write_2D_spreal,write_2D_int
 PUBLIC write_3D_int
 PUBLIC write_scalar_dpreal,write_scalar_spreal,write_scalar_int
 PUBLIC read_simple_dpreal,read_simple_spreal,read_simple_int
 PUBLIC read_2D_dpreal,read_2D_spreal,read_2D_int
 PUBLIC read_3D_int
 PUBLIC read_scalar_dpreal,read_scalar_spreal,read_scalar_int
 PUBLIC read_character,write_character
 
 PUBLIC  init_hdf, finalize_hdf, open_read_h5file, open_write_h5file, close_h5file
 
 INTEGER(HID_T) :: dsp_scalar_id
 INTEGER(HID_T),PUBLIC :: pcollectread_id,pcollectwrite_id ! property list  identifier 
 INTEGER(HID_T) :: compl_dpreal_id,compl_spreal_id  ! type for complex numbers (not implemented yet) 
 
 REAL(dpreal),PUBLIC :: hdftime=0.0_dpreal
 CONTAINS

! init_hdf starts up the hdf system and prepares calls
! called once 
 
  SUBROUTINE init_hdf
   IMPLICIT NONE 
   INTEGER ierr
   INTEGER(SIZE_T) re_size,im_size,complex_t_size,offset ! 4-b
   INTEGER(HSIZE_T) size ! long integer

   IF(master) THEN
    WRITE(*,*) 'hdf tools version: ',VERREV
!    WRITE(*,*) 'Date             : ',VERDATE
   END IF
   
   hdftime=hdftime-MPI_WTIME()
   
   ! first startup hdf5 
   CALL h5open_f(ierr)
   
!     set up property for collective read transfer
   CALL h5pcreate_f (H5P_DATASET_XFER_F, pcollectread_id, ierr)
#ifdef HDF_MPIO
   IF(master) THEN
    WRITE(*,*) 'HDFinit: Using MPIO in HDF for collective read.'
   END IF 
!   CALL h5pset_dxpl_mpio_f (pcollectxfer_id,H5FD_MPIO_INDEPENDENT_F,ierr)
   CALL h5pset_dxpl_mpio_f (pcollectread_id,H5FD_MPIO_COLLECTIVE_F,ierr)
#else
   IF(master) THEN
    WRITE(*,*) 'HDFinit: Using standard driver in HDF for all reading.'
   END IF    
#endif
   size=5*1024*1024  ! buffer size (standard = 1 MB)
   CALL h5pset_buffer_f(pcollectread_id, size, ierr)
   ! size=64*1024      ! sieve buffer size (standard = 64 kB) 
   ! CALL h5pset_sieve_buf_size_f(pcollectread_id, size, ierr)

   !     set up property for collective write  transfer
   CALL h5pcreate_f (H5P_DATASET_XFER_F, pcollectwrite_id, ierr)
#ifdef HDF_MPIO
   IF(master) THEN
    WRITE(*,*) 'HDFinit: Using MPIO in HDF for collective write.'
   ENDIF 
   ! CALL h5pset_dxpl_mpio_f (pcollectwrite_id,H5FD_MPIO_INDEPENDENT_F,ierr)
   CALL h5pset_dxpl_mpio_f (pcollectwrite_id,H5FD_MPIO_COLLECTIVE_F,ierr)
#else
   IF(master) THEN
    WRITE(*,*) 'HDFinit: Using standard driver in HDF for all writing.'
   END IF    
#endif   
   size=5*1024*1024  ! buffer size (standard = 1 MB)
   CALL h5pset_buffer_f(pcollectwrite_id, size, ierr)
   ! size=512*1024      ! sieve buffer size (standard = 64 kB) 
   ! CALL h5pset_sieve_buf_size_f(pcollectwrite_id, size, ierr)
 
   !     generate dataspace and memspace for single integer and single float       
   CALL h5screate_f(H5S_SCALAR_F,dsp_scalar_id,ierr)
   
   ! create HDF type for complex double precision
   
   CALL h5tget_size_f(H5T_NATIVE_DOUBLE,re_size,ierr)
   CALL h5tget_size_f(H5T_NATIVE_DOUBLE,im_size,ierr)
   complex_t_size = re_size + im_size
   CALL h5tcreate_f(H5T_COMPOUND_F, complex_t_size, compl_dpreal_id,ierr)
   offset = 0 ! define the location of first byte of the real part
   CALL h5tinsert_f(compl_dpreal_id,'real', offset, H5T_NATIVE_DOUBLE,ierr)
   offset = offset + re_size
   CALL h5tinsert_f(compl_dpreal_id,'imaginary', offset, H5T_NATIVE_DOUBLE,ierr)
   
   ! create HDF type for complex single precision
   
   CALL h5tget_size_f(H5T_NATIVE_REAL,re_size,ierr)
   CALL h5tget_size_f(H5T_NATIVE_REAL,im_size,ierr)
   complex_t_size = re_size + im_size
   CALL h5tcreate_f(H5T_COMPOUND_F, complex_t_size, compl_spreal_id,ierr)
   offset = 0 
   CALL h5tinsert_f(compl_spreal_id,'real', offset, H5T_NATIVE_REAL,ierr)
   offset = offset + re_size
   CALL h5tinsert_f(compl_spreal_id,'imaginary', offset, H5T_NATIVE_REAL,ierr)
   
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE init_hdf
 
! finalize_hdf shuts down the hdf system
! called once 

  SUBROUTINE finalize_hdf
   IMPLICIT NONE 
   INTEGER ierr
   hdftime=hdftime-MPI_WTIME()
   ! close complex types 
   CALL h5tclose_f(compl_spreal_id,ierr)
   CALL h5tclose_f(compl_dpreal_id,ierr)
   ! close property for collective transfer
   CALL h5pclose_f (pcollectread_id, ierr)
   CALL h5pclose_f (pcollectwrite_id, ierr)
   ! close dataspace and memspace for single integer and single float 
   CALL h5sclose_f(dsp_scalar_id,ierr)
   ! finally stop hdf5 
   CALL h5close_f(ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE finalize_hdf
  
  
  ! open file "file" using communicator comm 
  ! gives back file id f_id
  ! it is expected that a file named fname exists 
  ! file is opened in RW mode
  
  SUBROUTINE open_read_h5file(fname,comm,f_id)
   IMPLICIT NONE
   CHARACTER(LEN=*) fname
   INTEGER comm
   INTEGER(HID_T) f_id
   INTEGER ierr,info

   INTEGER(HID_T) :: plist_id
   hdftime=hdftime-MPI_WTIME()
!     set up parallel I/O
!     use alpha communicator for parallelization
!     it is meant to access file only from mybuchid=0 in commbuch 
   CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
!  use mpi-io with hints to improve performance 
   CALL MPI_Info_create(info,ierr)
!  use for GPFS   
   CALL MPI_Info_set(info,'IBM_largeblock_io','true',ierr)

!   CALL h5pset_fapl_mpio_f(plist_id,comm,info,ierr)
!   CALL h5pset_fapl_mpiposix_f(plist_id,comm,.true.,ierr)
   
!  open file in mpi_mode  
   CALL h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, f_id, ierr,access_prp = plist_id)

!  close property list
   CALL h5pclose_f(plist_id, ierr)  
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE open_read_h5file
  
  ! open file "file" using communicator comm 
  ! gives back file id f_id
  ! it is expected that a file named fname exists 
  ! file is opened in RW mode
 

  SUBROUTINE open_write_newh5file(fname,comm,f_id)
   IMPLICIT NONE
   CHARACTER(LEN=*) fname
   INTEGER comm
   INTEGER(HID_T) f_id
   INTEGER ierr,info

   INTEGER(HID_T) :: plist_id
   hdftime=hdftime-MPI_WTIME()
!     set up parallel I/O
!     use alpha communicator for parallelization
!     it is meant to access file only from mybuchid=0 in commbuch 
   CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
!  use mpi-io with hints to improve performance 
   CALL MPI_Info_create(info,ierr)
!  use for GPFS   
   CALL MPI_Info_set(info,'IBM_largeblock_io','true',ierr)

   !set file-open property in mpi_file_open to plist_id, collective mode
!   CALL h5pset_fapl_mpio_f(plist_id,comm,info,ierr)
!   CALL h5pset_fapl_mpiposix_(plist_id,comm,.true.,ierr)
   
!  open file in mpi_mode  
   !CALL h5fopen_f(trim(fname), H5F_ACC_RDWR_F, f_id, ierr,access_prp = plist_id)
   call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,f_id,ierr,access_prp=plist_id)

!  close property list
   CALL h5pclose_f(plist_id, ierr)  
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE open_write_newh5file




  SUBROUTINE open_write_h5file(fname,comm,f_id)
   IMPLICIT NONE
   CHARACTER(LEN=*) fname
   INTEGER comm
   INTEGER(HID_T) f_id
   INTEGER ierr,info
   LOGICAL fileexists

   INTEGER(HID_T) :: plist_id
   INQUIRE(FILE=trim(fname),EXIST=fileexists)   
   IF(.NOT. fileexists) THEN ! if the file does not exit
    ! create a new file 
    CALL open_write_newh5file(fname,comm,f_id)
   ELSE
    hdftime=hdftime-MPI_WTIME()
    ! open the old file 
!     set up parallel I/O
!     use alpha communicator for parallelization
!     it is meant to access file only from mybuchid=0 in commbuch 
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
!  use mpi-io with hints to improve performance 
    CALL MPI_Info_create(info,ierr)
!  use for GPFS   
    CALL MPI_Info_set(info,'IBM_largeblock_io','true',ierr)

   !set file-open property in mpi_file_open to plist_id, collective mode
!    CALL h5pset_fapl_mpio_f(plist_id,comm,info,ierr)
!    CALL h5pset_fapl_mpiposix_f(plist_id,comm,.true.,ierr)
   
!  open file in mpi_mode  
    CALL h5fopen_f(trim(fname), H5F_ACC_RDWR_F, f_id, ierr,access_prp = plist_id)

!  close property list
    CALL h5pclose_f(plist_id, ierr)  
    hdftime=hdftime+MPI_WTIME()
   END IF
  END SUBROUTINE open_write_h5file
  
! close file with ID f_id  
  SUBROUTINE close_h5file(f_id)
   IMPLICIT NONE
   INTEGER ierr
   INTEGER(HID_T) f_id
   hdftime=hdftime-MPI_WTIME()
   CALL h5fclose_f(f_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE close_h5file
  
  
! writes a 2D dimensional data set in double precision 
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_2D_dpreal(group_id,tag,n,m,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(2)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   REAL(dpreal) x(n,m)
   INTEGER n,m
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id,dsp_simple_id
   hdftime=hdftime-MPI_WTIME()

   dims(1)=n               !  write needs dimensions of arrays 
   dims(2)=m               !  write needs dimensions of arrays 
   ! 
   CALL h5screate_simple_f(2,dims,dsp_simple_id,ierr);
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_DOUBLE,dsp_simple_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,x,dims,ierr)
   END IF
   CALL h5dclose_f(dset_id, ierr)
   CALL h5sclose_f(dsp_simple_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_2D_dpreal
  
! reads a 2D single precision array
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_2D_dpreal(group_id,dsetname,n,m,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname  ! name of the data want to copy
   INTEGER n,m
   REAL(dpreal) x(n,m)
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id 
   INTEGER(HSIZE_T) dims(2) ! one dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=n   ! set dimensions 
   dims(2)=m
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE,x,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(x,n*m, MPI_REAL8,0,comm,ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_2D_dpreal
    
! writes a 2D dimensional data set in single precision 
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_2D_spreal(group_id,tag,n,m,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(2)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   REAL(spreal) x(n,m)
   INTEGER n,m
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id,dsp_simple_id
   hdftime=hdftime-MPI_WTIME()
   dims(1)=n               !  write needs dimensions of arrays 
   dims(2)=m               !  write needs dimensions of arrays 
   ! 
   CALL h5screate_simple_f(2,dims,dsp_simple_id,ierr);
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_REAL,dsp_simple_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,x,dims,ierr)
   END IF
   CALL h5dclose_f(dset_id, ierr)
   CALL h5sclose_f(dsp_simple_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_2D_spreal
  
! reads a 2D single precision array
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_2D_spreal(group_id,dsetname,n,m,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   INTEGER n,m
   REAL(spreal) x(n,m)
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id 
   INTEGER(HSIZE_T) dims(2) ! one dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=n   ! set dimensions 
   dims(2)=m
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_REAL,x,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(x,n*m, MPI_REAL4,0,comm,ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_2D_spreal
   
 
! writes a 2D dimensional data set in integer(4 bytes)
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_2D_int(group_id,tag,n,m,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(2)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   INTEGER x(n,m)
   INTEGER n,m
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id,dsp_simple_id
   hdftime=hdftime-MPI_WTIME()
   dims(1)=n               !  write needs dimensions of arrays 
   dims(2)=m               !  write needs dimensions of arrays 
   ! 
   CALL h5screate_simple_f(2,dims,dsp_simple_id,ierr);
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_INTEGER,dsp_simple_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,x,dims,ierr)
   END IF
   CALL h5dclose_f(dset_id, ierr)
   CALL h5sclose_f(dsp_simple_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_2D_int
  
! reads a 2D  integer   array
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_2D_int(group_id,dsetname,n,m,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   INTEGER n,m
   INTEGER x(n,m)
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id 
   INTEGER(HSIZE_T) dims(2) ! one dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=n   ! set dimensions 
   dims(2)=m
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER,x,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(x,n*m, MPI_INTEGER,0,comm,ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_2D_int
    
 
! writes a 3D dimensional data set (integer)
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_3D_int(group_id,tag,n,m,l,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(3)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   INTEGER x(n,m,l)
   INTEGER n,m,l
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id,dsp_simple_id
   hdftime=hdftime-MPI_WTIME()
   dims(1)=n               !  write needs dimensions of arrays 
   dims(2)=m               !  write needs dimensions of arrays 
   dims(3)=l               !  write needs dimensions of arrays 
   ! 
   CALL h5screate_simple_f(3,dims,dsp_simple_id,ierr);
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_INTEGER,dsp_simple_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,x,dims,ierr)
   END IF
   CALL h5dclose_f(dset_id, ierr)
   CALL h5sclose_f(dsp_simple_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_3D_int
  
! reads a 3D integer array
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_3D_int(group_id,dsetname,n,m,l,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   INTEGER n,m,l
   INTEGER x(n,m,l)
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id 
   INTEGER(HSIZE_T) dims(3) ! one dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=n   ! set dimensions 
   dims(2)=m
   dims(3)=l
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER,x,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(x,n*m*l, MPI_INTEGER,0,comm,ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_3D_int
    
 
! writes a 3D data set in single precision 
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_3D_spreal(group_id,tag,n,m,l,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(3)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   REAL(spreal) x(n,m,l)
   INTEGER n,m,l
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id,dsp_simple_id
   hdftime=hdftime-MPI_WTIME()
   dims(1)=n               !  write needs dimensions of arrays 
   dims(2)=m               !  write needs dimensions of arrays 
   dims(3)=l               !  write needs dimensions of arrays 
   ! 
   CALL h5screate_simple_f(3,dims,dsp_simple_id,ierr);
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_REAL,dsp_simple_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,x,dims,ierr)
   END IF
   CALL h5dclose_f(dset_id, ierr)
   CALL h5sclose_f(dsp_simple_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_3D_spreal
  
! reads a 3D single precision array
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_3D_spreal(group_id,dsetname,n,m,l,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   INTEGER n,m,l
   REAL(spreal) x(n,m,l)
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id 
   INTEGER(HSIZE_T) dims(3) ! 3-dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=n   ! set dimensions 
   dims(2)=m
   dims(3)=l
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_REAL,x,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(x,n*m*l, MPI_REAL4,0,comm,ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_3D_spreal
   
 
! writes a 3D dimensional data set in double precision 
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_3D_dpreal(group_id,tag,n,m,l,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(3)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   REAL(dpreal) x(n,m,l)
   INTEGER n,m,l
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id,dsp_simple_id
   hdftime=hdftime-MPI_WTIME()

   dims(1)=n               !  write needs dimensions of arrays 
   dims(2)=m               !  write needs dimensions of arrays 
   dims(3)=l               !  write needs dimensions of arrays 
   ! 
   CALL h5screate_simple_f(3,dims,dsp_simple_id,ierr);
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_DOUBLE,dsp_simple_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,x,dims,ierr)
   END IF
   CALL h5dclose_f(dset_id, ierr)
   CALL h5sclose_f(dsp_simple_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_3D_dpreal
  
! reads a 3D double precision array
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_3D_dpreal(group_id,dsetname,n,m,l,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   INTEGER n,m,l
   REAL(dpreal) x(n,m,l)
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id 
   INTEGER(HSIZE_T) dims(3) ! 3-dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=n   ! set dimensions 
   dims(2)=m
   dims(2)=l
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE,x,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(x,n*m*l, MPI_REAL8,0,comm,ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_3D_dpreal

  
! writes a 1 dimensional data set in double precision 
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_simple_dpreal(group_id,tag,n,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(1)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   REAL(dpreal) x(n)
   INTEGER n
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id,dsp_simple_id
   
   hdftime=hdftime-MPI_WTIME()
   
   dims(1)=n               !  write needs dimensions of arrays 
   ! 
   CALL h5screate_simple_f(1,dims,dsp_simple_id,ierr);
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_DOUBLE,dsp_simple_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,x,dims,ierr)
   END IF
   CALL h5dclose_f(dset_id, ierr)
   CALL h5sclose_f(dsp_simple_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_simple_dpreal

! reads a simple single precision array
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_simple_dpreal(group_id,dsetname,n,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   INTEGER n
   REAL(dpreal) x(n)
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id 
   INTEGER(HSIZE_T) dims(1) ! one dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=n   ! only one data item 
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE,x,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(x,n, MPI_REAL8,0,comm,ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_simple_dpreal

  
! writes a 1 dimensional data set in single precision 
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_simple_spreal(group_id,tag,n,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(1)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   REAL(spreal) x(n)
   INTEGER n
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id,dsp_simple_id
   hdftime=hdftime-MPI_WTIME()
   dims(1)=n               !  write needs dimensions of arrays 
   ! 
   CALL h5screate_simple_f(1,dims,dsp_simple_id,ierr);
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_REAL,dsp_simple_id,dset_id,ierr) 
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,x,dims,ierr)
   END IF
   CALL h5dclose_f(dset_id, ierr)
   CALL h5sclose_f(dsp_simple_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_simple_spreal
  
! reads a simple single precision array
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_simple_spreal(group_id,dsetname,n,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   INTEGER n
   REAL(spreal) x(n)
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id 
   INTEGER(HSIZE_T) dims(1) ! one dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=n   ! only one data item 
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_REAL,x,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(x,n, MPI_REAL4,0,comm,ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_simple_spreal
   
  
! writes a 1 dimensional data set in integer
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_simple_int(group_id,tag,n,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(1)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   INTEGER x(n)
   INTEGER n
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id,dsp_simple_id
   hdftime=hdftime-MPI_WTIME()
   dims(1)=n               !  write needs dimensions of arrays 
   ! 
   CALL h5screate_simple_f(1,dims,dsp_simple_id,ierr);
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_INTEGER,dsp_simple_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,x,dims,ierr)
   END IF 
   CALL h5dclose_f(dset_id, ierr)
   CALL h5sclose_f(dsp_simple_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_simple_int

! reads a simple integer array
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_simple_int(group_id,dsetname,n,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   INTEGER n
   INTEGER x(n)
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id 
   INTEGER(HSIZE_T) dims(1) ! one dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=n   ! only one data item 
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER,x,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(x,n, MPI_INTEGER,0,comm,ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_simple_int
   
   
  SUBROUTINE write_scalar_dpreal(group_id,tag,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(1)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   REAL(dpreal) x
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id
   hdftime=hdftime-MPI_WTIME()
   dims(1)=1               !  write needs dimensions of arrays (here single value)
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_DOUBLE,dsp_scalar_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,x,dims,ierr)
   END IF
   CALL h5dclose_f(dset_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_scalar_dpreal

! reads a scalar double precision
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_scalar_dpreal(group_id,dsetname,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   REAL(dpreal) x
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id
   REAL(dpreal) dset(1)     ! need array for read subroutine 
   INTEGER(HSIZE_T) dims(1) ! one dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=1   ! only one data item 
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dset,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(dset,1, MPI_REAL8,0,comm,ierr)

   x=dset(1)   
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_scalar_dpreal

  
! writes a scalar single precision
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_scalar_spreal(group_id,tag,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(1)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   REAL(spreal) x
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id
   hdftime=hdftime-MPI_WTIME()
   dims(1)=1               !  write needs dimensions of arrays (here single value)
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_REAL,dsp_scalar_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL,x,dims,ierr)
   END IF
   CALL h5dclose_f(dset_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_scalar_spreal
  
! reads a scalar single precision
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_scalar_spreal(group_id,dsetname,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   REAL(spreal) x
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id
   REAL(spreal) dset(1)     ! need array for read subroutine 
   INTEGER(HSIZE_T) dims(1) ! one dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=1   ! only one data item 
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_REAL,dset,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(dset,1, MPI_REAL4,0,comm,ierr)

   x=dset(1)  
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_scalar_spreal

  
! writes a scalar integer
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_scalar_int(group_id,tag,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(1)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   INTEGER x
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id
   hdftime=hdftime-MPI_WTIME()
   dims(1)=1               !  write needs dimensions of arrays (here single value)
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_INTEGER,dsp_scalar_id,dset_id,ierr) 
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,x,dims,ierr)
   END IF 
   CALL h5dclose_f(dset_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_scalar_int

! reads a scalar integer
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_scalar_int(group_id,dsetname,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   INTEGER x
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id
   INTEGER dset(1)     ! need array for read subroutine 
   INTEGER(HSIZE_T) dims(1) ! one dim data set  
   INTEGER myid
   INTEGER comm
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   dims(1)=1   ! only one data item 
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(dset,1, MPI_INTEGER,0,comm,ierr)

   x=dset(1)   
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_scalar_int
 
! writes a string 
! non-parallel, only master=.true. writes the data
! the dataset is named tag and written to group group_id
  
  SUBROUTINE write_character(group_id,tag,x,master)
   IMPLICIT NONE
   LOGICAL master
   INTEGER(HSIZE_T) dims(1)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
   CHARACTER(LEN=*) tag
   CHARACTER(LEN=*) x
   INTEGER n
   INTEGER ierr
   INTEGER(HID_T) :: dset_id,group_id,dsp_simple_id
   hdftime=hdftime-MPI_WTIME()
   n=len_trim(x)           !  length without trailing spaces
   dims(1)=n               !  write needs dimensions of arrays 
   ! 
   CALL h5screate_simple_f(1,dims,dsp_simple_id,ierr)
   CALL h5dcreate_f(group_id,trim(tag),H5T_NATIVE_CHARACTER,dsp_simple_id,dset_id,ierr)  
   IF(master) THEN
    CALL h5dwrite_f(dset_id, H5T_NATIVE_CHARACTER,trim(x),dims,ierr)
   END IF 
   CALL h5dclose_f(dset_id, ierr)
   CALL h5sclose_f(dsp_simple_id, ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE write_character

! reads a string
! from dataset with name dsetname in group group_id
! non-parallel, id=0 reads the data
! the data is then broadcasted to the rest of the 
! processes in communicator comm

  SUBROUTINE read_character(group_id,dsetname,x,comm)
   IMPLICIT NONE
   CHARACTER(LEN=*) dsetname
   INTEGER n
   CHARACTER(LEN=*) x
   INTEGER ierr
   INTEGER(HID_T) dset_id,group_id,space_id 
   INTEGER myid
   INTEGER comm
   INTEGER(HSIZE_T) dims(1),maxdims(1)
   hdftime=hdftime-MPI_WTIME()
   ! find rank 
   CALL MPI_COMM_RANK(comm,myid,ierr)
   
   ! open the data set in current group 
   CALL h5dopen_f(group_id,trim(dsetname), dset_id, ierr)
   ! find length of the string 
   CALL h5dget_space_f(dset_id, space_id, ierr) 
   CALL h5sget_simple_extent_dims_f(space_id, dims, maxdims, ierr)
   n=dims(1)
   
   ! read the data on id = 0  
   IF(myid.EQ.0) THEN
    CALL h5dread_f(dset_id, H5T_NATIVE_CHARACTER,x,dims,ierr)
   END IF
   ! close data set again 
   CALL h5dclose_f(dset_id, ierr)

   ! broadcats to rest of processes
   CALL MPI_BCAST(x,n, MPI_CHARACTER,0,comm,ierr)
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE read_character
   
! next set of subroutines 
! opens the group in an existing file 
! writes/reads an additional datum to/from this group
!  (maybe add energy to wave function) 
  
  SUBROUTINE writedata_dpreal(filename,group,tag,comm,x)
   IMPLICIT NONE
   CHARACTER(LEN=*) filename,group,tag
   INTEGER comm,ierr,peid
   REAL(dpreal) :: x  
   INTEGER(HID_T) f_id,group_id
   hdftime=hdftime-MPI_WTIME()
!  get rank 
   
   CALL MPI_COMM_RANK(comm,peid,ierr)  
   
!  open file    
   CALL open_write_h5file(filename,comm,f_id)

!  open group 
   CALL h5gopen_f(f_id,trim(group),group_id,ierr)
   
!  write datum 
   
   CALL write_scalar_dpreal(group_id,trim(tag),x,peid.eq.0)
     
!  close group    
   CALL h5gclose_f(group_id, ierr)
   
   
   CALL close_h5file(f_id)
   
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE writedata_dpreal
  
  SUBROUTINE readdata_dpreal(filename,group,tag,comm,x)
   IMPLICIT NONE
   CHARACTER(LEN=*) filename,group,tag
   INTEGER comm,ierr
   REAL(dpreal) :: x  
   INTEGER(HID_T) f_id,group_id
   hdftime=hdftime-MPI_WTIME() 
   
!  open file    
   CALL open_read_h5file(filename,comm,f_id)

!  open group 
   CALL h5gopen_f(f_id,trim(group),group_id,ierr)
   
!  write datum 
   
   CALL read_scalar_dpreal(group_id,trim(tag),x,comm)
     
!  close group    
   CALL h5gclose_f(group_id, ierr)
   
   
   CALL close_h5file(f_id)
   
   hdftime=hdftime+MPI_WTIME()
  END SUBROUTINE readdata_dpreal
  
  
!! do complex later ....  
!! writes a scalar double precision complex
!! non-parallel, only master=.true. writes the data
!! the dataset is named tag and written to group group_id
!  
!  SUBROUTINE write_scalar_dpcompl(group_id,tag,x,master)
!   IMPLICIT NONE
!   LOGICAL master
!   INTEGER(HSIZE_T) dims(1)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
!   CHARACTER(LEN=*) tag
!   COMPLEX(dpreal) x
!   INTEGER ierr
!   INTEGER(HID_T) :: dset_id,group_id
!   hdftime=hdftime-MPI_WTIME()
!   dims(1)=1               !  write needs dimensions of arrays (here single value)
!   CALL h5dcreate_f(group_id,trim(tag),compl_dpreal_id,dsp_scalar_id,dset_id,ierr)  
!   IF(master) THEN
!    CALL h5dwrite_f(dset_id, compl_dpreal_id,x,dims,ierr)
!   END IF 
!   CALL h5dclose_f(dset_id, ierr)
!   hdftime=hdftime+MPI_WTIME()
!  END SUBROUTINE write_scalar_dpcompl
!
!! writes a scalar single precision complex
!! non-parallel, only master=.true. writes the data
!! the dataset is named tag and written to group group_id
!  
!  SUBROUTINE write_scalar_spcompl(group_id,tag,x,master)
!   IMPLICIT NONE
!   LOGICAL master
!   INTEGER(HSIZE_T) dims(1)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
!   CHARACTER(LEN=*) tag
!   COMPLEX(spreal) x
!   INTEGER ierr
!   INTEGER(HID_T) :: dset_id,group_id
!   TYPE complex_struct
!    REAL(spreal) :: re_part
!    REAL(spreal) :: im_part
!   END TYPE 
!   TYPE(complex_struct) :: y
!   hdftime=hdftime-MPI_WTIME()
!   dims(1)=1               !  write needs dimensions of arrays (here single value)
!   CALL h5dcreate_f(group_id,trim(tag),compl_spreal_id,dsp_scalar_id,dset_id,ierr)  
!   IF(master) THEN
!    CALL h5dwrite_f(dset_id, compl_spreal_id,x,dims,ierr)
!   END IF
!   CALL h5dclose_f(dset_id, ierr)
!   hdftime=hdftime+MPI_WTIME()
!  END SUBROUTINE write_scalar_spcompl
!  
END MODULE hdf_tool
 
