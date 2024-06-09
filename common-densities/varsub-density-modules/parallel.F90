#include "fdefs.h"

!! module parallel reads distribution of processors
!! starts up MPI and defines communicators

MODULE parallel
 USE precision
 USE mpi_const

 PRIVATE 
 
 ! npe parameters for public use 
 INTEGER,PUBLIC :: npe_p12,npe_q4,npe_p3,npe_alpha,npe_ener,npe
 
 ! communicators defined 
 INTEGER,PUBLIC :: commp12,commq4,commp3,commalpha,commener
 INTEGER,PUBLIC :: commmesh,commampnn,commamp3n,commamp4n,commall
 INTEGER,PUBLIC :: commp12p3,commp12alpha,commp3alpha,commp3q4ener,&
                 & commq4ener,commp3q4,commp12p3alphaener,&
                 & commp3q4alphaener,commp3q4alpha,commp12q4alphaener

 ! id's within communicators 
 INTEGER,PUBLIC :: myid,myp12id,myq4id,myp3id,myalphaid,myenerid, &
             &      mymeshid,myampnnid,myamp3nid,myamp4nid,&
             &      myp12p3id,myp12alphaid,myp3alphaid,myp3q4enerid,&
             &      myq4enerid,myp3q4id,myp3q4alphaenerid,myp12q4alphaenerid,&
             &      myp3q4alphaid,myp12p3alphaenerid

 ! is true for PE=0 only 
 ! on this processors prints to stdout shall be performed 

 ! subroutines provided by the module 
 PUBLIC initparallel,printparallel
 

CONTAINS 

 ! routine initializes MPI and defines basic 
 ! variables used to parallize the code 
 ! mpi_finalize needs to be called at the end of the main program 
 ! in order to have CPU's finish their tasks 

 SUBROUTINE initparallel
  IMPLICIT NONE
  INTEGER ierr

  ! start up MPI library

  CALL MPI_INIT(ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  ! read in npe parameters  
  CALL readpara
 
  ! check consistency of read in npe_* with total number of npe

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  IF(npe_p12*npe_q4*npe_p3*npe_alpha*npe_ener.NE.npe) &
      STOP 'npe inconsistent'

  ! set up the communicators 

  CALL initcomm  

  ! fix one PE for printout
  master=(myid.EQ.0) 
  
  ! if possible redirect output to nul device except on 
  ! master to get potential output only on master processors
#ifdef SUPPRESSOUT
  IF(.NOT.master) THEN
   OPEN(unit=stdout, file="/dev/null", status="OLD")     ! nullfile="/dev/null"
#ifdef BGQ   
   OPEN(unit=stderr, file="error.out", status="OLD")     ! nullfile="error.out" on BGQ
#else
   OPEN(unit=stderr, file="/dev/null", status="OLD")     ! nullfile="/dev/null"
#endif   
  END IF
  IF(master) THEN
   WRITE(*,*) 'Redirect output to NULL device on non-master processors'
  END IF
#endif  
 END SUBROUTINE initparallel

 ! prints out parameters used for parallelization 
 ! should be called only from master processor 
 ! in order to have only one copy in the output 

 SUBROUTINE printparallel
  IMPLICIT NONE
  INTEGER ncomm,ierr
  
  IF(master) THEN   
   WRITE(*,*) 
   WRITE(*,*) 'actual size of communicators:'

   CALL MPI_COMM_SIZE(commp12,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commp12:     ',ncomm

   CALL MPI_COMM_SIZE(commq4,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commq4:      ',ncomm

   CALL MPI_COMM_SIZE(commp3,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commp3:      ',ncomm

   CALL MPI_COMM_SIZE(commmesh,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commmesh:    ',ncomm

   CALL MPI_COMM_SIZE(commampnn,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commampnn:    ',ncomm

   CALL MPI_COMM_SIZE(commamp3n,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commamp3n:    ',ncomm

   CALL MPI_COMM_SIZE(commamp4n,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commamp4n:    ',ncomm

   CALL MPI_COMM_SIZE(commalpha,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commalpha:   ',ncomm

   CALL MPI_COMM_SIZE(commener,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commener:    ',ncomm   

   CALL MPI_COMM_SIZE(commall,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commall:     ',ncomm   

   CALL MPI_COMM_SIZE(commp12p3,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commp12p3:   ',ncomm

   CALL MPI_COMM_SIZE(commp12alpha,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commp12alpha:',ncomm
   
   CALL MPI_COMM_SIZE(commp3alpha,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commp3alpha:',ncomm
   
   CALL MPI_COMM_SIZE(commp3q4ener,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commp3q4ener:',ncomm

   CALL MPI_COMM_SIZE(commp3q4,ncomm,ierr)
   IF(ierr.NE.0) STOP 'mpi failure'
   WRITE(*,*) 'commp3q4:',ncomm

  END IF ! master

#ifdef DEBUG
  ! if debugging print independent info on all processors 

  CALL MPI_COMM_SIZE(commp12,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp12:     ',myid,ncomm,myp12id

  CALL MPI_COMM_SIZE(commq4,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commq4:      ',myid,ncomm,myq4id

  CALL MPI_COMM_SIZE(commp3,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp3:      ',myid,ncomm,myp3id

  CALL MPI_COMM_SIZE(commmesh,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commmesh:    ',myid,ncomm,mymeshid

  CALL MPI_COMM_SIZE(commampnn,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commampnn:    ',myid,ncomm,myampnnid

  CALL MPI_COMM_SIZE(commamp3n,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commamp3n:    ',myid,ncomm,myamp3nid

  CALL MPI_COMM_SIZE(commamp4n,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commamp4n:    ',myid,ncomm,myamp4nid

  CALL MPI_COMM_SIZE(commalpha,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commalpha:   ',myid,ncomm,myalphaid

  CALL MPI_COMM_SIZE(commener,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commener:    ',myid,ncomm,myenerid

  CALL MPI_COMM_SIZE(commall,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commall:     ',myid,ncomm,myid

  CALL MPI_COMM_SIZE(commp12p3,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp12p3:     ',myid,ncomm,myp12p3id

  CALL MPI_COMM_SIZE(commp12alpha,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp12alpha:  ',myid,ncomm,myp12alphaid

  CALL MPI_COMM_SIZE(commp3alpha,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp3alpha:  ',myid,ncomm,myp3alphaid
  
  CALL MPI_COMM_SIZE(commp3q4ener,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp3q4ener:  ',myid,ncomm,myp3q4enerid
  
  CALL MPI_COMM_SIZE(commp3q4,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp3q4:  ',myid,ncomm,myp3q4id
  
  CALL MPI_COMM_SIZE(commq4ener,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commq4ener:  ',myid,ncomm,myq4enerid
  
  CALL MPI_COMM_SIZE(commp3q4alphaener,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp3q4alphaener:  ',myid,ncomm,myp3q4alphaenerid
  
  CALL MPI_COMM_SIZE(commp3q4alpha,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp3q4alpha:  ',myid,ncomm,myp3q4alphaid
  
  CALL MPI_COMM_SIZE(commp12q4alphaener,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp12q4alphaener:  ',myid,ncomm,myp12q4alphaenerid
  
  CALL MPI_COMM_SIZE(commp12p3alphaener,ncomm,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  WRITE(*,'(A,3I8)') 'commp12p3alphaener:  ',myid,ncomm,myp12p3alphaenerid
  
#endif

 END SUBROUTINE printparallel

 ! subroutine builds cartesian communicator 
 ! and 
 SUBROUTINE initcomm
  IMPLICIT NONE
  INTEGER ierr
  INTEGER dims(5)
  LOGICAL per(5),ret(5) 


  ! first all PE's as cartesian communicator 
  
  dims(1)=npe_p12
  dims(2)=npe_p3
  dims(3)=npe_q4
  dims(4)=npe_alpha
  dims(5)=npe_ener
  per(1)=.FALSE.
  per(2)=.FALSE.
  per(3)=.FALSE.
  per(4)=.FALSE.
  per(5)=.FALSE.
  
  CALL MPI_CART_CREATE(MPI_COMM_WORLD,5,dims,per,.false.,commall,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commall,myid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'


  ! set up subcommunicators 
  ! ret variable=.true. implies that 
  ! PE's with different id's in this dimension remain part of the communicator
  ! ret variable=.false. means that 
  ! PE's with different id's in this dimension are in different communicators 

  ! set up subcommunicator for mesh 

  ret(1)=.true.    
  ret(2)=.true.
  ret(3)=.true.
  ret(4)=.false.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commmesh,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commmesh,mymeshid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  ! set up subcommunicators for mesh & channels  (amplitudes)
  ! nn case 
  ret(1)=.true.    
  ret(2)=.false.
  ret(3)=.false.
  ret(4)=.true.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commampnn,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commampnn,myampnnid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  ! 3n case 
  ret(1)=.true.    
  ret(2)=.true.
  ret(3)=.false.
  ret(4)=.true.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commamp3n,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commamp3n,myamp3nid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  ! nn case 
  ret(1)=.true.    
  ret(2)=.true.
  ret(3)=.true.
  ret(4)=.true.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commamp4n,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commamp4n,myamp4nid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  ! set up subcommunicator for p12 

  ret(1)=.true.    
  ret(2)=.false.
  ret(3)=.false.
  ret(4)=.false.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commp12,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp12,myp12id,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'


  ! set up subcommunicator for q4

  ret(1)=.false.    
  ret(2)=.false.
  ret(3)=.true.
  ret(4)=.false.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commq4,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commq4,myq4id,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  ! set up subcommunicator for p3

  ret(1)=.false.    
  ret(2)=.true.
  ret(3)=.false.
  ret(4)=.false.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commp3,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp3,myp3id,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'


  ! set up subcommunicator for alpha

  ret(1)=.false.    
  ret(2)=.false.
  ret(3)=.false.
  ret(4)=.true.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commalpha,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commalpha,myalphaid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  ! set up subcommunicator for energies

  ret(1)=.false.    
  ret(2)=.false.
  ret(3)=.false.
  ret(4)=.false.
  ret(5)=.true.
  
  CALL MPI_CART_SUB(commall,ret,commener,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commener,myenerid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'


  ! set up subcommunicator for p12 and p3 for G_A

  ret(1)=.true.    
  ret(2)=.true.
  ret(3)=.false.
  ret(4)=.false.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commp12p3,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp12p3,myp12p3id,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  ! set up subcommunicator for p12alpha

  ret(1)=.true.    
  ret(2)=.false.
  ret(3)=.false.
  ret(4)=.true.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commp12alpha,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp12alpha,myp12alphaid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  
  ! set up subcommunicator for p3alpha

  ret(1)=.false.    
  ret(2)=.true.
  ret(3)=.false.
  ret(4)=.true.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commp3alpha,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp3alpha,myp3alphaid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  
  
  ret(1)=.false.    
  ret(2)=.true.
  ret(3)=.true.
  ret(4)=.false.
  ret(5)=.true.
  
  CALL MPI_CART_SUB(commall,ret,commp3q4ener,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp3q4ener,myp3q4enerid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'

  ret(1)=.false.    
  ret(2)=.true.
  ret(3)=.true.
  ret(4)=.false.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commp3q4,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp3q4,myp3q4id,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  
  ret(1)=.false.    
  ret(2)=.false.
  ret(3)=.true.
  ret(4)=.false.
  ret(5)=.true.
  
  CALL MPI_CART_SUB(commall,ret,commq4ener,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commq4ener,myq4enerid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  
  ret(1)=.false.    
  ret(2)=.true.
  ret(3)=.true.
  ret(4)=.true.
  ret(5)=.true.
  
  CALL MPI_CART_SUB(commall,ret,commp3q4alphaener,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp3q4alphaener,myp3q4alphaenerid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  
  ret(1)=.false.    
  ret(2)=.true.
  ret(3)=.true.
  ret(4)=.true.
  ret(5)=.false.
  
  CALL MPI_CART_SUB(commall,ret,commp3q4alpha,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp3q4alpha,myp3q4alphaid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  
  ret(1)=.true.    
  ret(2)=.false.
  ret(3)=.true.
  ret(4)=.true.
  ret(5)=.true.
  
  CALL MPI_CART_SUB(commall,ret,commp12q4alphaener,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp12q4alphaener,myp12q4alphaenerid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  
  ret(1)=.true.    
  ret(2)=.true.
  ret(3)=.false.
  ret(4)=.true.
  ret(5)=.true.
  
  CALL MPI_CART_SUB(commall,ret,commp12p3alphaener,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  CALL MPI_COMM_RANK(commp12p3alphaener,myp12p3alphaenerid,ierr)
  IF(ierr.NE.0) STOP 'mpi failure'
  

 END SUBROUTINE initcomm
 
!! subroutine reads in the file parallel.dat
!! that contains the distribution of processes  
!! over the CPUs 

 SUBROUTINE readpara
  IMPLICIT NONE 

  OPEN(UNIT=75,FILE='parallel.dat',FORM='FORMATTED',STATUS='OLD')
  READ(75,*,END=100,ERR=100) 
  READ(75,*,END=100,ERR=100) npe_p12,npe_p3,npe_q4,npe_alpha,npe_ener
  CLOSE(75)

  RETURN

100 CONTINUE 
  !     treatment of I/O errors 
  STOP 'problem with parallel.dat'
 END SUBROUTINE readpara

END MODULE parallel
