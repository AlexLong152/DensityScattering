module meshpoints
 use precision
 use constants
 use gauss
 use parallel
 use mpi_const
 use hdf5
 use hdf_tool

 implicit none
 private
 public meshinit,meshprint,writep12p,readp12p,&
      & writep3p,readp3p,writeq4p,readq4p, &
      & writep34p,readp34p,writeqp,readqp,&
      & writephysconst,readphysconst
 
 integer :: p12n1,p12n2,p3n1,p3n2,q4n1,q4n2,qn1,qn2,p34n1,p34n2,ip12,ip3,iq4,ip34,iq
 integer,public :: p12n,p3n,q4n,qn,p34n,xintn,phiintn,xintpoln
 real(dpreal),allocatable :: p12gp(:,:),p12gpp(:,:),p3gp(:,:),p3gpp(:,:),q4gp(:,:),q4gpp(:,:)
 real(dpreal),allocatable :: qgp(:,:),qgpp(:,:),p34gp(:,:),p34gpp(:,:)
 real(dpreal),public :: p12p1,p12p2,p12p3,p3p1,p3p2,p3p3,q4p1,q4p2,q4p3,qp1,qp2,qp3,p34p1,p34p2,p34p3
 real(dpreal),allocatable,public :: p12p(:),p3p(:),q4p(:),qp(:),p34p(:),p12w(:),p3w(:),q4w(:),qw(:),p34w(:)
 real(dpreal),allocatable,public :: xintp(:),xintw(:),phiintp(:),phiintw(:),xintpolp(:),xintpolw(:)
 
 ! sets of parameters for TNF files, meshtype = 'TNF'
 real(dpreal) :: PG_tnf(5),QG_tnf(5)
 integer      :: NP_tnf(4),NQ_tnf(4)
 real(dpreal) :: QMIN,QROM1,QROM2,QROM3,QROM4,QROM5,QROM6
 real(dpreal) :: PMIN,PROM1,PROM2,PROM3,PROM4,PROM5,PROM6
 integer :: NP1ROM,NPINT1ROM,NPINT2ROM,NQ1ROM,NQINT1ROM,NQINT2ROM
 
 ! r-space grids 
 integer :: r12n1,r12n2,r3n1,r3n2,r4n1,r4n2,rn1,rn2,r34n1,r34n2,ir12,ir3,ir4,ir34,ir
 real(dpreal),public :: r12p1,r12p2,r12p3,r3p1,r3p2,r3p3,r4p1,r4p2,r4p3,rp1,rp2,rp3,r34p1,r34p2,r34p3
 integer,public :: r12n,r3n,r4n,rn,r34n
 real(dpreal),allocatable,public :: r12p(:),r3p(:),r4p(:),rp(:),r34p(:),r12w(:),r3w(:),r4w(:),rw(:),r34w(:)

 
 logical,public :: scatt
 character(4) :: meshtype,rmeshtype
 
contains
 
 subroutine meshinit
  implicit none
  real(dpreal) qbar,pbar 
  integer ip12,ip3
  
  ! start with p-space mesh 
  ! Parameters are read from the file 'meshpoints.dat'
  open(unit=75,file='meshpoints.dat',form='formatted',status='old')

  ! Read mesh point distribution type
  read(75,*,end=100,err=100) meshtype

  if(trim(meshtype).NE.'FLEX') then
    ! Read mesh parameters
    read(75,*,end=100,err=100) p12p1,p12p2,p12p3,p12n1,p12n2
    read(75,*,end=100,err=100) p3p1,p3p2,p3p3,p3n1,p3n2
    read(75,*,end=100,err=100) q4p1,q4p2,q4p3,q4n1,q4n2
    read(75,*,end=100,err=100) qp1,qp2,qp3,qn1,qn2
    read(75,*,end=100,err=100) p34p1,p34p2,p34p3,p34n1,p34n2
    read(75,*,end=100,err=100) xintn,phiintn,xintpoln
    if(trim(meshtype).EQ.'TNF') then
     read(75,*,end=100,err=100) PG_tnf(1:5)
     read(75,*,end=100,err=100) NP_tnf(1:4)
     read(75,*,end=100,err=100) QG_tnf(1:5)
     read(75,*,end=100,err=100) NQ_tnf(1:4)
    end if
    if(trim(meshtype).EQ.'ROM') then
     read(75,*,end=100,err=100) PMIN,PROM1,PROM2,PROM3,PROM4,PROM5,PROM6  
     read(75,*,end=100,err=100) NP1ROM,NPINT1ROM,NPINT2ROM
     read(75,*,end=100,err=100) QMIN,QROM1,QROM2,QROM3,QROM4,QROM5,QROM6
     read(75,*,end=100,err=100) NQ1ROM,NQINT1ROM,NQINT2ROM
    end if
    if(trim(meshtype).EQ.'READ') then ! read P12 and P3 grids 
       read(75,*,end=100,err=100) p12n,p3n    ! read number of points first
     allocate(p12p(p12n),p12w(p12n))
     allocate(p3p(p3n),p3w(p3n))
     DO ip12=1,P12N 
      read(75,*,end=100,err=100) p12p(ip12),p12w(ip12)
     END DO
     DO ip3=1,P3N 
      read(75,*,end=100,err=100) p3p(ip3),p3w(ip3)
     END DO
    endif
    close(75)
  
  ! Definition of meshpoint parameters
    if(trim(meshtype).eq.'TNF') then
     p12n=np_tnf(1)+np_tnf(2)+np_tnf(3)+np_tnf(4)
     p3n=nq_tnf(1)+nq_tnf(2)+nq_tnf(3)+nq_tnf(4)
     q4n=q4n1+q4n2
     qn=qn1+qn2
     p34n=p34n1+p34n2
    else if(trim(meshtype).eq.'ROM') then
     p12n=NP1ROM+NPINT1ROM+NPINT2ROM+2
     p3n=NQ1ROM+NQINT1ROM+NQINT2ROM+2
     q4n=q4n1+q4n2
     qn=qn1+qn2
     p34n=p34n1+p34n2
    else if(trim(meshtype).eq.'READ') then
     q4n=q4n1+q4n2
     qn=qn1+qn2
     p34n=p34n1+p34n2 
    else
     p12n=p12n1+p12n2
     p3n=p3n1+p3n2
     q4n=q4n1+q4n2
     qn=qn1+qn2
     p34n=p34n1+p34n2
    end if
  
  ! Momentum vectors and weights are allocated
    if (trim(meshtype).ne.'READ') then
     allocate(p12p(p12n))
     allocate(p3p(p3n))
     allocate(p12w(p12n))
     allocate(p3w(p3n))
    end if
  
    allocate(q4p(q4n))
    allocate(qp(qn))
    allocate(p34p(p34n))
  
  
    allocate(q4w(q4n))
    allocate(qw(qn))
    allocate(p34w(p34n))
  
    if (trim(meshtype).eq.'TRNS') then
   
    ! Integration points are calculated with subroutine "trns"
     call trns(p12n1,p12n2,p12n,p12p1,p12p2,p12p3,p12p,p12w)
     call trns(p3n1,p3n2,p3n,p3p1,p3p2,p3p3,p3p,p3w)
     call trns(q4n1,q4n2,q4n,q4p1,q4p2,q4p3,q4p,q4w)
     call trns(qn1,qn2,qn,qp1,qp2,qp3,qp,qw)
     call trns(p34n1,p34n2,p34n,p34p1,p34p2,p34p3,p34p,p34w)

    else if (trim(meshtype).eq.'CHEB') then

     allocate(p12gp(p12n,p12n))
     allocate(p12gpp(p12n,p12n))

     allocate(p3gp(p3n,p3n))
     allocate(p3gpp(p3n,p3n))

     allocate(q4gp(q4n,q4n))
     allocate(q4gpp(q4n,q4n))

     allocate(qgp(qn,qn))
     allocate(qgpp(qn,qn))

     allocate(p34gp(p34n,p34n))
     allocate(p34gpp(p34n,p34n))
    
     call chebyshev(p12n,p12p2,p12p3,p12p,p12w,p12gp,p12gpp)
     call chebyshev(p3n,p3p2,p3p3,p3p,p3w,p3gp,p3gpp)
     call chebyshev(q4n,q4p2,q4p3,q4p,q4w,q4gp,q4gpp)
     call chebyshev(qn,qp2,qp3,qp,qw,qgp,qgpp)
     call chebyshev(p34n,p34p2,p34p3,p34p,p34w,p34gp,p34gpp)

     deallocate(p12gp)
     deallocate(p12gpp)

     deallocate(p3gp)
     deallocate(p3gpp)

     deallocate(q4gp)
     deallocate(q4gpp)

     deallocate(qgp)
     deallocate(qgpp)

     deallocate(p34gp)
     deallocate(p34gpp)

    else if(trim(meshtype).eq.'TNF') then
     if(Q4N.NE.1) stop 'TNF mesh is only 3N!!!'
     call trns(q4n1,q4n2,q4n,q4p1,q4p2,q4p3,q4p,q4w)
     call trns(qn1,qn2,qn,qp1,qp2,qp3,qp,qw)
     call trns(p34n1,p34n2,p34n,p34p1,p34p2,p34p3,p34p,p34w)
     CALL IMPULSE (NP_tnf,PG_tnf,P12P,P12W,PBAR)
     CALL IMPULSE (NQ_tnf,QG_tnf,P3P,P3W,QBAR)
    else if(trim(meshtype).eq.'ROM') then
     if(Q4N.NE.1) stop 'ROM mesh is only 3N!!!'
     call trns(q4n1,q4n2,q4n,q4p1,q4p2,q4p3,q4p,q4w)
     call trns(qn1,qn2,qn,qp1,qp2,qp3,qp,qw)
     call trns(p34n1,p34n2,p34n,p34p1,p34p2,p34p3,p34p,p34w)
     call QV4SETN(PMIN,PROM1,PROM2,PROM3,PROM4,PROM5,PROM6,NP1ROM,NPINT1ROM,NPINT2ROM,P12P,P12W)
     call QV4SETN(QMIN,QROM1,QROM2,QROM3,QROM4,QROM5,QROM6,NQ1ROM,NQINT1ROM,NQINT2ROM,P3P,P3W)
    else if(trim(meshtype).eq.'READ') then
     if(Q4N.NE.1) stop 'READ mesh is only 3N!!!'
     call trns(q4n1,q4n2,q4n,q4p1,q4p2,q4p3,q4p,q4w)
     call trns(qn1,qn2,qn,qp1,qp2,qp3,qp,qw)
     call trns(p34n1,p34n2,p34n,p34p1,p34p2,p34p3,p34p,p34w)
     ! mesh points have been fixed already
    else
     stop 'meshtype unknown!'

    end if
  else    ! FLEX points 
   read(75,*,end=100,err=100) xintn,phiintn,xintpoln  ! first read angular grid points 

   ! define P12 grid 
   read(75,*,end=100,err=100) meshtype   ! reuse meshtype now for individual momenta 
   SELECT CASE(trim(meshtype))           ! meshtype for P12 points
    CASE('TRNS')
      read(75,*,end=100,err=100) p12p1,p12p2,p12p3,p12n1,p12n2
      p12n=p12n1+p12n2
      allocate(p12p(p12n),p12w(p12n))
      call trns(p12n1,p12n2,p12n,p12p1,p12p2,p12p3,p12p,p12w)
    CASE('CHEB')
      read(75,*,end=100,err=100) p12p2,p12p3,p12n
      p12p1=0.0
      p12n1=0
      p12n2=0
      allocate(p12p(p12n),p12w(p12n))
      allocate(p12gp(p12n,p12n),p12gpp(p12n,p12n))
      call chebyshev(p12n,p12p2,p12p3,p12p,p12w,p12gp,p12gpp)
      deallocate(p12gp,p12gpp)
    CASE('READ')  
      read(75,*,end=100,err=100) p12n
      p12p1=0.0
      p12p2=0.0
      p12p3=0.0
      p12n1=0
      p12n2=0
      allocate(p12p(p12n),p12w(p12n))
      do ip12=1,P12N 
        read(75,*,end=100,err=100) p12p(ip12),p12w(ip12)
      end do
    CASE DEFAULT
      stop 'meshtype unknown!'
   END SELECT
      

   ! define P3 grid 
   read(75,*,end=100,err=100) meshtype   ! reuse meshtype now for individual momenta 
   SELECT CASE(trim(meshtype))           ! meshtype for P3 points
    CASE('TRNS')
      read(75,*,end=100,err=100) p3p1,p3p2,p3p3,p3n1,p3n2
      p3n=p3n1+p3n2
      allocate(p3p(p3n),p3w(p3n))
      call trns(p3n1,p3n2,p3n,p3p1,p3p2,p3p3,p3p,p3w)
    CASE('CHEB')
      read(75,*,end=100,err=100) p3p2,p3p3,p3n
      p3p1=0.0
      p3n1=0
      p3n2=0
      allocate(p3p(p3n),p3w(p3n))
      allocate(p3gp(p3n,p3n),p3gpp(p3n,p3n))
      call chebyshev(p3n,p3p2,p3p3,p3p,p3w,p3gp,p3gpp)
      deallocate(p3gp,p3gpp)
    CASE('READ')  
      read(75,*,end=100,err=100) p3n
      p3p1=0.0
      p3p2=0.0
      p3p3=0.0
      p3n1=0
      p3n2=0
      allocate(p3p(p3n),p3w(p3n))
      do ip3=1,P3N 
        read(75,*,end=100,err=100) p3p(ip3),p3w(ip3)
      end do
    CASE DEFAULT
      stop 'meshtype unknown!'
   END SELECT
      

   ! define Q4 grid 
   read(75,*,end=100,err=100) meshtype   ! reuse meshtype now for individual momenta 
   SELECT CASE(trim(meshtype))           ! meshtype for Q4 points
    CASE('TRNS')
      read(75,*,end=100,err=100) q4p1,q4p2,q4p3,q4n1,q4n2
      q4n=q4n1+q4n2
      allocate(q4p(q4n),q4w(q4n))
      call trns(q4n1,q4n2,q4n,q4p1,q4p2,q4p3,q4p,q4w)
    CASE('CHEB')
      read(75,*,end=100,err=100) q4p2,q4p3,q4n
      q4p1=0.0
      q4n1=0
      q4n2=0
      allocate(q4p(q4n),q4w(q4n))
      allocate(q4gp(q4n,q4n),q4gpp(q4n,q4n))
      call chebyshev(q4n,q4p2,q4p3,q4p,q4w,q4gp,q4gpp)
      deallocate(q4gp,q4gpp)
    CASE('READ')  
      read(75,*,end=100,err=100) q4n
      q4p1=0.0
      q4p2=0.0
      q4p3=0.0
      q4n1=0
      q4n2=0
      allocate(q4p(q4n),q4w(q4n))
      do iq4=1,Q4N 
        read(75,*,end=100,err=100) q4p(iq4),q4w(iq4)
      end do
    CASE DEFAULT
      stop 'meshtype unknown!'
   END SELECT
      

   ! define P34 grid 
   read(75,*,end=100,err=100) meshtype   ! reuse meshtype now for individual momenta 
   SELECT CASE(trim(meshtype))           ! meshtype for P34 points
    CASE('TRNS')
      read(75,*,end=100,err=100) p34p1,p34p2,p34p3,p34n1,p34n2
      p34n=p34n1+p34n2
      allocate(p34p(p34n),p34w(p34n))
      call trns(p34n1,p34n2,p34n,p34p1,p34p2,p34p3,p34p,p34w)
    CASE('CHEB')
      read(75,*,end=100,err=100) p34p2,p34p3,p34n
      p34p1=0.0
      p34n1=0
      p34n2=0
      allocate(p34p(p34n),p34w(p34n))
      allocate(p34gp(p34n,p34n),p34gpp(p34n,p34n))
      call chebyshev(p34n,p34p2,p34p3,p34p,p34w,p34gp,p34gpp)
      deallocate(p34gp,p34gpp)
    CASE('READ')  
      read(75,*,end=100,err=100) p34n
      p34p1=0.0
      p34p2=0.0
      p34p3=0.0
      p34n1=0
      p34n2=0
      allocate(p34p(p34n),p34w(p34n))
      do ip34=1,P34N 
        read(75,*,end=100,err=100) p34p(ip34),p34w(ip34)
      end do
    CASE DEFAULT
      stop 'meshtype unknown!'
   END SELECT
      

   ! define Q grid 
   read(75,*,end=100,err=100) meshtype   ! reuse meshtype now for individual momenta 
   SELECT CASE(trim(meshtype))           ! meshtype for Q points
    CASE('TRNS')
      read(75,*,end=100,err=100) qp1,qp2,qp3,qn1,qn2
      qn=qn1+qn2
      allocate(qp(qn),qw(qn))
      call trns(qn1,qn2,qn,qp1,qp2,qp3,qp,qw)
    CASE('CHEB')
      read(75,*,end=100,err=100) qp2,qp3,qn
      qp1=0.0
      qn1=0
      qn2=0
      allocate(qp(qn),qw(qn))
      allocate(qgp(qn,qn),qgpp(qn,qn))
      call chebyshev(qn,qp2,qp3,qp,qw,qgp,qgpp)
      deallocate(qgp,qgpp)
    CASE('READ')  
      read(75,*,end=100,err=100) qn
      qp1=0.0
      qp2=0.0
      qp3=0.0
      qn1=0
      qn2=0
      allocate(qp(qn),qw(qn))
      do iq=1,QN 
        read(75,*,end=100,err=100) qp(iq),qw(iq)
      end do
    CASE DEFAULT
      stop 'meshtype unknown!'
   END SELECT
      
  end if  ! FLEX points 
  ! Grid points for anglular integration and interpolation are computed
  
  allocate(xintp(xintn))
  allocate(xintw(xintn))
  call gauleg(xintn,-1.0_dpreal,1.0_dpreal,xintp,xintw)
  
  allocate(phiintp(phiintn))
  allocate(phiintw(phiintn))
  call gauleg(phiintn,0.0_dpreal,2*pi,phiintp,phiintw)

  allocate(xintpolp(xintpoln))
  allocate(xintpolw(xintpoln))
  call gauleg(xintpoln,-1.0_dpreal,1.0_dpreal,xintpolp,xintpolw)
  
 ! now r-space mesh 
  
  open(unit=75,file='rmeshpoints.dat',form='formatted',status='old')

  ! Read mesh point distribution type
  read(75,*,end=100,err=100) rmeshtype

  SELECT CASE(trim(rmeshtype))
   CASE('TRNS')  ! Read mesh parameters for TRNS 
    read(75,*,end=100,err=100) r12p1,r12p2,r12p3,r12n1,r12n2
    read(75,*,end=100,err=100) r3p1,r3p2,r3p3,r3n1,r3n2
    read(75,*,end=100,err=100) r4p1,r4p2,r4p3,r4n1,r4n2
    read(75,*,end=100,err=100) rp1,rp2,rp3,rn1,rn2
    read(75,*,end=100,err=100) r34p1,r34p2,r34p3,r34n1,r34n2
    r12n=r12n1+r12n2
    r3n=r3n1+r3n2
    r4n=r4n1+r4n2
    rn=rn1+rn2
    r34n=r34n1+r34n2
    ! allocate grid point arrays 
    allocate(r12p(r12n))
    allocate(r3p(r3n))
    allocate(r12w(r12n))
    allocate(r3w(r3n))
    allocate(r4p(r4n))
    allocate(rp(rn))
    allocate(r34p(r34n))
    allocate(r4w(r4n))
    allocate(rw(rn))
    allocate(r34w(r34n))
    
    ! Integration points are calculated with subroutine "trns"
    call trns(r12n1,r12n2,r12n,r12p1,r12p2,r12p3,r12p,r12w)
    call trns(r3n1,r3n2,r3n,r3p1,r3p2,r3p3,r3p,r3w)
    call trns(r4n1,r4n2,r4n,r4p1,r4p2,r4p3,r4p,r4w)
    call trns(rn1,rn2,rn,rp1,rp2,rp3,rp,rw)
    call trns(r34n1,r34n2,r34n,r34p1,r34p2,r34p3,r34p,r34w)
   
   CASE('READ') ! read P12 and P3 grids 
    read(75,*,end=100,err=100) r12n,r3n,r4n,r34n,rn    ! read number of points first    
    allocate(r12p(r12n),r12w(r12n))
    allocate(r3p(r3n),r3w(r3n))
    allocate(r4p(r4n),r4w(r4n))
    allocate(r34p(r34n),r34w(r34n))
    allocate(rp(rn),rw(rn))
    DO ir12=1,R12N 
     read(75,*,end=100,err=100) r12p(ir12),r12w(ir12)
    END DO
    DO ir3=1,R3N 
     read(75,*,end=100,err=100) r3p(ir3),r3w(ir3)
    END DO
    DO ir4=1,R4N 
     read(75,*,end=100,err=100) r4p(ir4),r4w(ir4)
    END DO
    DO ir34=1,R34N 
     read(75,*,end=100,err=100) r34p(ir34),r34w(ir34)
    END DO
    DO ir=1,RN 
     read(75,*,end=100,err=100) rp(ir),rw(ir)
    END DO
   CASE DEFAULT
    STOP 'meshtype not defined'
  END SELECT   
    
  close(75)


  return
  ! Regular end of subroutine
100 continue
  ! Treatment of i/o errors
  stop 'problem with meshpoints.dat or rmeshpoints.dat'
 end subroutine meshinit
 
 subroutine meshprint
  implicit none
  real(dpreal) sum
  
  ! Output of the meshpoint parameters
  if (master) then
   IF(meshtype.eq.'TNF') THEN
    write(*,'(A,5E15.6)') 'pg parameters: ',pg_tnf(1:5)
    write(*,'(A,4I5)')    'np parameters: ',np_tnf(1:4)
    
    do ip12=1,p12n
     write(*,*) 'momentum point p12',ip12,'(in fm^-1)=',p12p(ip12)
    end do
   
    write(*,*)

    write(*,'(A,5E15.6)') 'pg parameters: ',pg_tnf(1:5)
    write(*,'(A,4I5)')    'np parameters: ',np_tnf(1:4)
    
    do ip3=1,p3n
     write(*,*) 'momentum point p3',ip3,'(in fm^-1)=',p3p(ip3)
    end do

   ELSE IF(trim(meshtype).eq.'ROM') THEN
    WRITE(*,'(A,7E15.6)') 'P parameters:', PMIN,PROM1,PROM2,PROM3,PROM4,PROM5,PROM6  
    WRITE(*,'(A,3I5)')    'NP parameters:', NP1ROM,NPINT1ROM,NPINT2ROM
    WRITE(*,'(A,7E15.6)') 'Q parameters:', QMIN,QROM1,QROM2,QROM3,QROM4,QROM5,QROM6
    WRITE(*,'(A,3I5)')    'NQ parameters:', NQ1ROM,NQINT1ROM,NQINT2ROM
    
   ELSE
    write(*,'(A,A,E15.6,A,E15.6,A,E15.6,A,I4,A,I4,A,I4)') & 
    & 'p12 parameters: ',' p1=',p12p1,' p2=',p12p2,' p3=',p12p3,' n1=',p12n1,' n2=',p12n2,' n_ges=',p12n
    do ip12=1,p12n
     write(*,'(A,I6,2E20.8)') 'P12MESH: ',ip12,p12p(ip12),p12w(ip12)
    end do
   
    write(*,*)

    write(*,'(A,A,E15.6,A,E15.6,A,E15.6,A,I4,A,I4,A,I4)') &
    & 'p3 parameters:  ',' p1=',p3p1,' p2=',p3p2,' p3=',p3p3,' n1=',p3n1,' n2=',p3n2,' n_ges=',p3n
    do ip3=1,p3n
     write(*,'(A,I6,2E20.8)') 'P3MESH:  ',ip3,p3p(ip3),p3w(ip3)
    end do
   END IF 
   write(*,*)

   write(*,'(A,A,E15.6,A,E15.6,A,E15.6,A,I4,A,I4,A,I4)') & 
   & 'q4 parameters:  ',' p1=',q4p1,' p2=',q4p2,' p3=',q4p3,' n1=',q4n1,' n2=',q4n2,' n_ges=',q4n
   do iq4=1,q4n
    write(*,'(A,I6,2E20.8)') 'Q4MESH:  ',iq4,q4p(iq4),q4w(iq4)
   end do

   write(*,*)

   write(*,'(A,A,E15.6,A,E15.6,A,E15.6,A,I4,A,I4,A,I4)') & 
   & 'q parameters:   ',' p1=',qp1,' p2=',qp2,' p3=',qp3,' n1=',qn1,' n2=',qn2,' n_ges=',qn
   do iq=1,qn
    write(*,'(A,I6,2E20.8)') 'QMESH:   ',iq,qp(iq),qw(iq)
   end do

   write(*,*)

   write(*,'(A,A,E15.6,A,E15.6,A,E15.6,A,I4,A,I4,A,I4)') &
   & 'p34 parameters: ',' p1=',p34p1,' p2=',p34p2,' p3=',p34p3,' n1=',p34n1,' n2=',p34n2,' n_ges=',p34n
   do ip34=1,p34n
    write(*,'(A,I6,2E20.8)') 'P34MESH: ',ip34,p34p(ip34),p34w(ip34)
   end do
   
  end if
  
  ! basic test of weights 
  IF(master) THEN
   sum=0.0
   DO ip12=1,P12N
    sum=sum+P12W(ip12)
   END DO
   WRITE(*,'(A,E20.10)') 'SUMP12WEIGHT: ',sum

   sum=0.0
   DO ip3=1,P3N
    sum=sum+P3W(ip3)
   END DO
   WRITE(*,'(A,E20.10)') 'SUMP3WEIGHT: ',sum
   sum=0.0
   DO iq4=1,Q4N
    sum=sum+Q4W(iq4)
   END DO
   WRITE(*,'(A,E20.10)') 'SUMQ4WEIGHT: ',sum
   sum=0.0
   DO ip34=1,P34N
    sum=sum+P34W(ip34)
   END DO
   WRITE(*,'(A,E20.10)') 'SUMP34WEIGHT: ',sum
   sum=0.0
   DO iQ=1,QN
    sum=sum+QW(iq)
   END DO
   WRITE(*,'(A,E20.10)') 'SUMQWEIGHT: ',sum
  END IF
  
 end subroutine meshprint
 
 ! output routine for meshpoints 
 ! it assumed that a group id opened into which 
 ! the information is written 
 ! 
 
 ! 1. p12 grid 
 ! info: meshtype,p12p1,p12p2,p12p3,
 !       p12n1,p12n2,p12n,p12p,p12w
 ! maybe other meshtypes than TRNS require more info
 ! this is not implemented yet 
 
 SUBROUTINE writep12p(group_amp_id)
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  
  ! first real parameters 
  
  CALL write_scalar_dpreal(group_amp_id,'p12p1',p12p1,master)
  CALL write_scalar_dpreal(group_amp_id,'p12p2',p12p2,master)
  CALL write_scalar_dpreal(group_amp_id,'p12p3',p12p3,master)
  
  ! then type 
  CALL write_character(group_amp_id,'meshtype',meshtype,master)
  
  ! then integer parameters 
  CALL write_scalar_int(group_amp_id,'p12n1',p12n1,master)
  CALL write_scalar_int(group_amp_id,'p12n2',p12n2,master)
  CALL write_scalar_int(group_amp_id,'p12ntot',p12n,master)
  
  ! then meshpoints and weights 
  CALL write_simple_dpreal(group_amp_id,'p12p',p12n,p12p,master)
  CALL write_simple_dpreal(group_amp_id,'p12w',p12n,p12w,master)
  
 END SUBROUTINE writep12p
 
 ! 2. p3 grid 
 ! info: meshtype,p3p1,p3p2,p3p3,
 !       p3n1,p3n2,p3n,p3p,p3w
 ! maybe other meshtypes than TRNS require more info
 ! this is not implemented yet 
 
 SUBROUTINE writep3p(group_amp_id)
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  
  ! first real parameters 
  
  CALL write_scalar_dpreal(group_amp_id,'p3p1',p3p1,master)
  CALL write_scalar_dpreal(group_amp_id,'p3p2',p3p2,master)
  CALL write_scalar_dpreal(group_amp_id,'p3p3',p3p3,master)
  
  ! then type 
  ! assume that this is stored with p12 
  ! CALL write_character(group_amp_id,'meshtype',meshtype,master)
  
  ! then integer parameters 
  CALL write_scalar_int(group_amp_id,'p3n1',p3n1,master)
  CALL write_scalar_int(group_amp_id,'p3n2',p3n2,master)
  CALL write_scalar_int(group_amp_id,'p3ntot',p3n,master)
  
  ! then meshpoints and weights 
  CALL write_simple_dpreal(group_amp_id,'p3p',p3n,p3p,master)
  CALL write_simple_dpreal(group_amp_id,'p3w',p3n,p3w,master)
  
 END SUBROUTINE writep3p
 
 ! 3. q4 grid 
 ! info: meshtype,q4p1,q4p2,q4p3,
 !       q4n1,q4n2,q4n,q4p,q4w
 ! maybe other meshtypes than TRNS require more info
 ! this is not implemented yet 
 
 SUBROUTINE writeq4p(group_amp_id)
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  
  ! first real parameters 
  
  CALL write_scalar_dpreal(group_amp_id,'q4p1',q4p1,master)
  CALL write_scalar_dpreal(group_amp_id,'q4p2',q4p2,master)
  CALL write_scalar_dpreal(group_amp_id,'q4p3',q4p3,master)
  
  ! then type 
  ! assume that this was written with p12 
  ! CALL write_character(group_amp_id,'meshtype',meshtype,master)
  
  ! then integer parameters 
  CALL write_scalar_int(group_amp_id,'q4n1',q4n1,master)
  CALL write_scalar_int(group_amp_id,'q4n2',q4n2,master)
  CALL write_scalar_int(group_amp_id,'q4ntot',q4n,master)
  
  ! then meshpoints and weights 
  CALL write_simple_dpreal(group_amp_id,'q4p',q4n,q4p,master)
  CALL write_simple_dpreal(group_amp_id,'q4w',q4n,q4w,master)
  
 END SUBROUTINE writeq4p
 
 ! 4. p34 grid 
 ! info: meshtype,p34p1,p34p2,p34p3,
 !       p34n1,p34n2,p34n,p34p,p34w
 ! maybe other meshtypes than TRNS require more info
 ! this is not implemented yet 
 
 SUBROUTINE writep34p(group_amp_id)
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  
  ! first real parameters 
  
  CALL write_scalar_dpreal(group_amp_id,'p34p1',p34p1,master)
  CALL write_scalar_dpreal(group_amp_id,'p34p2',p34p2,master)
  CALL write_scalar_dpreal(group_amp_id,'p34p3',p34p3,master)
  
  ! then type 
  ! CALL write_character(group_amp_id,'meshtype',meshtype,master)
  ! assume that this was written with p12 grid 
  
  ! then integer parameters 
  CALL write_scalar_int(group_amp_id,'p34n1',p34n1,master)
  CALL write_scalar_int(group_amp_id,'p34n2',p34n2,master)
  CALL write_scalar_int(group_amp_id,'p34ntot',p34n,master)
  
  ! then meshpoints and weights 
  CALL write_simple_dpreal(group_amp_id,'p34p',p34n,p34p,master)
  CALL write_simple_dpreal(group_amp_id,'p34w',p34n,p34w,master)
  
 END SUBROUTINE writep34p
 
 ! 5. q grid 
 ! info: meshtype,qp1,qp2,qp3,
 !       qn1,qn2,qn,qp,qw
 ! maybe other meshtypes than TRNS require more info
 ! this is not implemented yet 
 
 SUBROUTINE writeqp(group_amp_id)
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  
  ! first real parameters 
  
  CALL write_scalar_dpreal(group_amp_id,'qp1',qp1,master)
  CALL write_scalar_dpreal(group_amp_id,'qp2',qp2,master)
  CALL write_scalar_dpreal(group_amp_id,'qp3',qp3,master)
  
  ! then type 
  ! CALL write_character(group_amp_id,'meshtype',meshtype,master)
  ! assume that this was written with p12 grid 
  
  ! then integer parameters 
  CALL write_scalar_int(group_amp_id,'qn1',qn1,master)
  CALL write_scalar_int(group_amp_id,'qn2',qn2,master)
  CALL write_scalar_int(group_amp_id,'qntot',qn,master)
  
  ! then meshpoints and weights 
  CALL write_simple_dpreal(group_amp_id,'qp',qn,qp,master)
  CALL write_simple_dpreal(group_amp_id,'qw',qn,qw,master)
  
 END SUBROUTINE writeqp
 
 
 ! input routine for meshpoints 
 ! it assumed that a group id opened from which 
 ! the information is read 
 ! supplemental info is written out to stdout
 ! meshpoints are read to given array
 ! these arrays are automatically allocated 
 
 ! 1. p12 grid 
 ! info printed out: meshtype,p12p1,p12p2,p12p3,
 !       p12n1,p12n2,p12n
 ! maybe other meshtypes than TRNS require more info
 ! this is not implemented yet 
 
 SUBROUTINE readp12p(group_amp_id,meshp,meshw,nsize)
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  REAL(dpreal),ALLOCATABLE :: meshp(:),meshw(:)
  INTEGER nsize 
  REAL(dpreal) :: p12p1_temp,p12p2_temp,p12p3_temp
  INTEGER :: p12n1_temp,p12n2_temp,ip
  CHARACTER(LEN=10) :: meshtype_temp
  
  ! first read the real parameters 
  
  CALL read_scalar_dpreal(group_amp_id,'p12p1',p12p1_temp,commall)
  CALL read_scalar_dpreal(group_amp_id,'p12p2',p12p2_temp,commall)
  CALL read_scalar_dpreal(group_amp_id,'p12p3',p12p3_temp,commall)
  
  ! then type 
  CALL read_character(group_amp_id,'meshtype',meshtype_temp,commall)
  
  ! then integer parameters 
  CALL read_scalar_int(group_amp_id,'p12n1',p12n1_temp,commall)
  CALL read_scalar_int(group_amp_id,'p12n2',p12n2_temp,commall)
  CALL read_scalar_int(group_amp_id,'p12ntot',nsize,commall)
  
  ! then meshpoints and weights 
  ALLOCATE(meshp(nsize),meshw(nsize))
  CALL read_simple_dpreal(group_amp_id,'p12p',nsize,meshp,commall)
  CALL read_simple_dpreal(group_amp_id,'p12w',nsize,meshw,commall)
  
  ! print out the parameters on master processor
  IF(master) THEN
   WRITE(*,*) 'P12 grid read in from file'
   
   write(*,'(A,A,E15.6,A,E15.6,A,E15.6,A,I5,A,I5,A,I5)') 'p12 read parameters: ',&
              ' p1=',p12p1_temp,&
   &          ' p2=',p12p2_temp,' p3=',p12p3_temp,&
   &          ' n1=',p12n1_temp,' n2=',p12n2_temp,' n_ges=',nsize
!  do ip=1,nsize
!   write(*,'(A,I4,2E15.6)') 'p12_read:',ip,meshp(ip),meshw(ip)
!  end do
   
  END IF 
  
 END SUBROUTINE readp12p

 ! 2. p3 grid 
 ! info printed out: meshtype,p3p1,p3p2,p3p3,
 !       p3n1,p3n2,p3n
 ! maybe other meshtypes than TRNS require more info
 ! this is not implemented yet 
 
 SUBROUTINE readp3p(group_amp_id,meshp,meshw,nsize)
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  REAL(dpreal),ALLOCATABLE :: meshp(:),meshw(:)
  INTEGER nsize 
  REAL(dpreal) :: p3p1_temp,p3p2_temp,p3p3_temp
  INTEGER :: p3n1_temp,p3n2_temp,ip
  CHARACTER(LEN=10) :: meshtype_temp
  
  ! first read the real parameters 
  
  CALL read_scalar_dpreal(group_amp_id,'p3p1',p3p1_temp,commall)
  CALL read_scalar_dpreal(group_amp_id,'p3p2',p3p2_temp,commall)
  CALL read_scalar_dpreal(group_amp_id,'p3p3',p3p3_temp,commall)
  
  ! then type 
  CALL read_character(group_amp_id,'meshtype',meshtype_temp,commall)
  
  ! then integer parameters 
  CALL read_scalar_int(group_amp_id,'p3n1',p3n1_temp,commall)
  CALL read_scalar_int(group_amp_id,'p3n2',p3n2_temp,commall)
  CALL read_scalar_int(group_amp_id,'p3ntot',nsize,commall)
  
  ! then meshpoints and weights 
  ALLOCATE(meshp(nsize),meshw(nsize))
  CALL read_simple_dpreal(group_amp_id,'p3p',nsize,meshp,commall)
  CALL read_simple_dpreal(group_amp_id,'p3w',nsize,meshw,commall)
  
  ! print out the parameters on master processor
  IF(master) THEN
   WRITE(*,*) 'P3 grid read in from file'
   
   write(*,'(A,A,E15.6,A,E15.6,A,E15.6,A,I5,A,I5,A,I5)') 'p3 read parameters: ',&
              ' p1=',p3p1_temp,&
   &          ' p2=',p3p2_temp,' p3=',p3p3_temp,&
   &          ' n1=',p3n1_temp,' n2=',p3n2_temp,' n_ges=',nsize
   do ip=1,nsize
    write(*,'(A,I4,2E15.6)') 'p3_read:',ip,meshp(ip),meshw(ip)
   end do
   
  END IF 
  
 END SUBROUTINE readp3p

 ! 3. q4 grid 
 ! info printed out: meshtype,q4p1,q4p2,q4p3,
 !       q4n1,q4n2,q4n
 ! maybe other meshtypes than TRNS require more info
 ! this is not implemented yet 
 
 SUBROUTINE readq4p(group_amp_id,meshp,meshw,nsize)
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  REAL(dpreal),ALLOCATABLE :: meshp(:),meshw(:)
  INTEGER nsize 
  REAL(dpreal) :: q4p1_temp,q4p2_temp,q4p3_temp
  INTEGER :: q4n1_temp,q4n2_temp,ip
  CHARACTER(LEN=10) :: meshtype_temp
  
  ! first read the real parameters 
  
  CALL read_scalar_dpreal(group_amp_id,'q4p1',q4p1_temp,commall)
  CALL read_scalar_dpreal(group_amp_id,'q4p2',q4p2_temp,commall)
  CALL read_scalar_dpreal(group_amp_id,'q4p3',q4p3_temp,commall)
  
  ! then type 
  CALL read_character(group_amp_id,'meshtype',meshtype_temp,commall)
  
  ! then integer parameters 
  CALL read_scalar_int(group_amp_id,'q4n1',q4n1_temp,commall)
  CALL read_scalar_int(group_amp_id,'q4n2',q4n2_temp,commall)
  CALL read_scalar_int(group_amp_id,'q4ntot',nsize,commall)
  
  ! then meshpoints and weights 
  ALLOCATE(meshp(nsize),meshw(nsize))
  CALL read_simple_dpreal(group_amp_id,'q4p',nsize,meshp,commall)
  CALL read_simple_dpreal(group_amp_id,'q4w',nsize,meshw,commall)
  
  ! print out the parameters on master processor
  IF(master) THEN
   WRITE(*,*) 'Q4 grid read in from file'
   
   write(*,'(A,A,E15.6,A,E15.6,A,E15.6,A,I5,A,I5,A,I5)') 'q4 read parameters: ',&
              ' p1=',q4p1_temp,&
   &          ' p2=',q4p2_temp,' p3=',q4p3_temp,&
   &          ' n1=',q4n1_temp,' n2=',q4n2_temp,' n_ges=',nsize
   do ip=1,nsize
    write(*,'(A,I4,2E15.6)') 'q4_read:',ip,meshp(ip),meshw(ip)
   end do
   
  END IF 
  
 END SUBROUTINE readq4p

 ! 4. p34 grid 
 ! info printed out: meshtype,p34p1,p34p2,p34p3,
 !       p34n1,p34n2,p34n
 ! maybe other meshtypes than TRNS require more info
 ! this is not implemented yet 
 
 SUBROUTINE readp34p(group_amp_id,meshp,meshw,nsize)
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  REAL(dpreal),ALLOCATABLE :: meshp(:),meshw(:)
  INTEGER nsize 
  REAL(dpreal) :: p34p1_temp,p34p2_temp,p34p3_temp
  INTEGER :: p34n1_temp,p34n2_temp,ip
  CHARACTER(LEN=10) :: meshtype_temp
  
  ! first read the real parameters 
  
  CALL read_scalar_dpreal(group_amp_id,'p34p1',p34p1_temp,commall)
  CALL read_scalar_dpreal(group_amp_id,'p34p2',p34p2_temp,commall)
  CALL read_scalar_dpreal(group_amp_id,'p34p3',p34p3_temp,commall)
  
  ! then type 
  CALL read_character(group_amp_id,'meshtype',meshtype_temp,commall)
  
  ! then integer parameters 
  CALL read_scalar_int(group_amp_id,'p34n1',p34n1_temp,commall)
  CALL read_scalar_int(group_amp_id,'p34n2',p34n2_temp,commall)
  CALL read_scalar_int(group_amp_id,'p34ntot',nsize,commall)
  
  ! then meshpoints and weights 
  ALLOCATE(meshp(nsize),meshw(nsize))
  CALL read_simple_dpreal(group_amp_id,'p34p',nsize,meshp,commall)
  CALL read_simple_dpreal(group_amp_id,'p34w',nsize,meshw,commall)
  
  ! print out the parameters on master processor
  IF(master) THEN
   WRITE(*,*) 'P34 grid read in from file'
   
   write(*,'(A,A,E15.6,A,E15.6,A,E15.6,A,I5,A,I5,A,I5)') 'p34 read parameters: ',&
              ' p1=',p34p1_temp,&
   &          ' p2=',p34p2_temp,' p3=',p34p3_temp,&
   &          ' n1=',p34n1_temp,' n2=',p34n2_temp,' n_ges=',nsize
   do ip=1,nsize
    write(*,'(A,I4,2E15.6)') 'p34_read:',ip,meshp(ip),meshw(ip)
   end do
   
  END IF 
  
 END SUBROUTINE readp34p

 ! 5. q grid 
 ! info printed out: meshtype,qp1,qp2,qp3,
 !       qn1,qn2,qn
 ! maybe other meshtypes than TRNS require more info
 ! this is not implemented yet 
 
 SUBROUTINE readqp(group_amp_id,meshp,meshw,nsize)
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  REAL(dpreal),ALLOCATABLE :: meshp(:),meshw(:)
  INTEGER nsize 
  REAL(dpreal) :: qp1_temp,qp2_temp,qp3_temp
  INTEGER :: qn1_temp,qn2_temp,ip
  CHARACTER(LEN=10) :: meshtype_temp
  
  ! first read the real parameters 
  
  CALL read_scalar_dpreal(group_amp_id,'qp1',qp1_temp,commall)
  CALL read_scalar_dpreal(group_amp_id,'qp2',qp2_temp,commall)
  CALL read_scalar_dpreal(group_amp_id,'qp3',qp3_temp,commall)
  
  ! then type 
  CALL read_character(group_amp_id,'meshtype',meshtype_temp,commall)
  
  ! then integer parameters 
  CALL read_scalar_int(group_amp_id,'qn1',qn1_temp,commall)
  CALL read_scalar_int(group_amp_id,'qn2',qn2_temp,commall)
  CALL read_scalar_int(group_amp_id,'qntot',nsize,commall)
  
  ! then meshpoints and weights 
  ALLOCATE(meshp(nsize),meshw(nsize))
  CALL read_simple_dpreal(group_amp_id,'qp',nsize,meshp,commall)
  CALL read_simple_dpreal(group_amp_id,'qw',nsize,meshw,commall)
  
  ! print out the parameters on master processor
  IF(master) THEN
   WRITE(*,*) 'Q grid read in from file'
   
   write(*,'(A,A,E15.6,A,E15.6,A,E15.6,A,I5,A,I5,A,I5)') 'q read parameters: ',&
              ' p1=',qp1_temp,&
   &          ' p2=',qp2_temp,' p3=',qp3_temp,&
   &          ' n1=',qn1_temp,' n2=',qn2_temp,' n_ges=',nsize
   do ip=1,nsize
    write(*,'(A,I4,2E15.6)') 'q_read:',ip,meshp(ip),meshw(ip)
   end do
   
  END IF 
  
 END SUBROUTINE readqp

 ! write out the physical constants to file maybe this should be moved elsewhere later 
 
 SUBROUTINE writephysconst(group_amp_id)   ! put to phys const ????
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  
  CALL write_scalar_dpreal(group_amp_id,'hbarc',hbarc,master)
  CALL write_scalar_dpreal(group_amp_id,'mprot',mprot,master)
  CALL write_scalar_dpreal(group_amp_id,'mneu',mneu,master)
  CALL write_scalar_dpreal(group_amp_id,'mlam',mlam,master)
  CALL write_scalar_dpreal(group_amp_id,'msigp',msigp,master)
  CALL write_scalar_dpreal(group_amp_id,'msig0',msig0,master)
  CALL write_scalar_dpreal(group_amp_id,'msigm',msigm,master)
  CALL write_scalar_dpreal(group_amp_id,'m_pip',m_pip,master)
  CALL write_scalar_dpreal(group_amp_id,'m_pi0',m_pi0,master)
  CALL write_scalar_dpreal(group_amp_id,'m_pi',m_pi,master)
  
  CALL write_scalar_dpreal(group_amp_id,'mxi0',mxi0,master)
  CALL write_scalar_dpreal(group_amp_id,'mxim',mxim,master)
  CALL write_scalar_dpreal(group_amp_id,'mxiave',mxiave,master)
  
  CALL write_scalar_dpreal(group_amp_id,'msigave',msigave,master)
  CALL write_scalar_dpreal(group_amp_id,'mnuc',mnuc,master)
  
  CALL write_scalar_dpreal(group_amp_id,'m_eta',m_eta,master)
  CALL write_scalar_dpreal(group_amp_id,'m_rho',m_rho,master)
  CALL write_scalar_dpreal(group_amp_id,'m_omega',m_omega,master)
  
  CALL write_scalar_dpreal(group_amp_id,'m_kaonp',m_kaonp,master)
  CALL write_scalar_dpreal(group_amp_id,'m_kaon0',m_kaon0,master)
  CALL write_scalar_dpreal(group_amp_id,'m_kaon',m_kaon,master)
  
  CALL write_scalar_dpreal(group_amp_id,'alpha_fein',alpha_fein,master)
  
 END SUBROUTINE writephysconst
 
 ! read the physical constants to file maybe this should be moved elsewhere later 
 
 SUBROUTINE readphysconst(group_amp_id)   ! put to phys const ????
  IMPLICIT NONE
  INTEGER(HID_T) group_amp_id
  REAL(dpreal) :: hbarc,mprot,mneu,mlam,msigp,msig0,msigm,m_pip,m_pi0, &
  &   m_pi,mxi0,mxim,mxiave,msigave,mnuc,m_eta,m_rho,m_omega, &
  &   m_kaonp,m_kaon0,m_kaon,alpha_fein
  LOGICAL newconst_exist
  INTEGER ierr
  
  CALL read_scalar_dpreal(group_amp_id,'hbarc',hbarc,commall)
  CALL read_scalar_dpreal(group_amp_id,'mprot',mprot,commall)
  CALL read_scalar_dpreal(group_amp_id,'mneu',mneu,commall)
  CALL read_scalar_dpreal(group_amp_id,'mlam',mlam,commall)
  CALL read_scalar_dpreal(group_amp_id,'msigp',msigp,commall)
  CALL read_scalar_dpreal(group_amp_id,'msig0',msig0,commall)
  CALL read_scalar_dpreal(group_amp_id,'msigm',msigm,commall)
  CALL read_scalar_dpreal(group_amp_id,'m_pip',m_pip,commall)
  CALL read_scalar_dpreal(group_amp_id,'m_pi0',m_pi0,commall)
  CALL read_scalar_dpreal(group_amp_id,'m_pi',m_pi,commall)
  
  ! the following lines also work for older versions of the read/write constants 
  ! routines by checking whether constants are stored in the file 
  
  CALL h5lexists_f(group_amp_id,trim('mxi0'),newconst_exist,ierr)
  
  if(newconst_exist) then
   CALL read_scalar_dpreal(group_amp_id,'mxi0',mxi0,commall)
   CALL read_scalar_dpreal(group_amp_id,'mxim',mxim,commall)
   CALL read_scalar_dpreal(group_amp_id,'mxiave',mxiave,commall)
  
   CALL read_scalar_dpreal(group_amp_id,'msigave',msigave,commall)
   CALL read_scalar_dpreal(group_amp_id,'mnuc',mnuc,commall)
  
   CALL read_scalar_dpreal(group_amp_id,'m_eta',m_eta,commall)
   CALL read_scalar_dpreal(group_amp_id,'m_rho',m_rho,commall)
   CALL read_scalar_dpreal(group_amp_id,'m_omega',m_omega,commall)
   
   CALL read_scalar_dpreal(group_amp_id,'m_kaonp',m_kaonp,commall)
   CALL read_scalar_dpreal(group_amp_id,'m_kaon0',m_kaon0,commall)
   CALL read_scalar_dpreal(group_amp_id,'m_kaon',m_kaon,commall)
   
   CALL read_scalar_dpreal(group_amp_id,'alpha_fein',alpha_fein,commall)
  end if   
  
  IF(master) THEN
   WRITE(*,*) 'Phyiscal constants read in from file'
  
   WRITE(*,*) 'hbarc = ',hbarc 
   WRITE(*,*) 'mprot = ',mprot*hbarc
   WRITE(*,*) 'mneu  = ', mneu*hbarc
   WRITE(*,*) 'mlam  = ', mlam*hbarc
   WRITE(*,*) 'msigp = ',msigp*hbarc
   WRITE(*,*) 'msig0 = ',msig0*hbarc
   WRITE(*,*) 'msigm = ',msigm*hbarc
   WRITE(*,*) 'm_pip = ',m_pip*hbarc
   WRITE(*,*) 'm_pi0 = ',m_pi0*hbarc
   WRITE(*,*) 'm_pi  = ',m_pi*hbarc
   if(newconst_exist) then
    WRITE(*,*) 'mxi0      = ',mxi0*hbarc 
    WRITE(*,*) 'mxim      = ',mxim*hbarc 
    WRITE(*,*) 'mxiave    = ',mxiave*hbarc 
    
    WRITE(*,*) 'msigave   = ',msigave*hbarc 
    WRITE(*,*) 'mnuc      = ',mnuc*hbarc 
    
    WRITE(*,*) 'm_eta     = ',m_eta*hbarc 
    WRITE(*,*) 'm_rho     = ',m_rho*hbarc 
    WRITE(*,*) 'm_omega   = ',m_omega*hbarc 
    
    WRITE(*,*) 'm_kaonp   = ',m_kaonp*hbarc
    WRITE(*,*) 'm_kaon0   = ',m_kaon0*hbarc 
    WRITE(*,*) 'm_kaon    = ',m_kaon*hbarc 
    
    WRITE(*,*) 'alpha_fein= ',alpha_fein
    
   else
    WRITE(*,*) '...... no further constants in this file'   
   endif    
  END IF 
   
 END SUBROUTINE readphysconst

! the next subroutine implements Roman Skibinskis meshpoints for the 3NF
! this routine is used when ROM is set as the meshpoint flag
!    NQ1 gauss points are put in [0,Q1] 
!    NQINT1/2 points in [Q2,Q3] and [Q3,Q4]
!    NQINT2   points in [Q4,Q5]
!    additionally the first points is changed to QMIN 
!    and two further points are added at the end between Q5 and Q6 
!       and at Q6 
!  !!! these grid points should not be used for integration, since the weights 
!  are partly distorted !!!  
 SUBROUTINE QV4SETN(QMIN,Q1,Q2,Q3,Q4,Q5,Q6,NQ1,NQINT1,NQINT2,Q,GQ)
  IMPLICIT NONE
  REAL(dpreal) :: Q(NQ1+NQINT1+NQINT2+2),GQ(NQ1+NQINT1+NQINT2+2)
  REAL(dpreal) :: QMIN,Q1,Q2,Q3,Q4,Q5,Q6
  INTEGER NQ1,NQINT1,NQINT2,NQ2,NQ
!
!*****Q-points for calculation of V4 and V4(1+P). While using this
!     subroutine put NQ=16 in PARAMETER STATEMENTS.
!     (Obtained from QV4SET)
!

  REAL(dpreal) :: X(NQ1),GX(NQ1),X2(NQ1),GX2(NQ1),X1(NQINT1+NQINT2),GX1(NQINT1+NQINT2)
  INTEGER I
  NQ2=NQINT1+NQINT2
  NQ=NQ1+NQINT1+NQINT2+2
  
  CALL gauleg(NQ1,0.0_dpreal,Q1,X2,GX2)
      
  DO I=1,NQ1
   Q(I)=X2(I)
   GQ(I)=GX2(I)
  END DO 
  Q(1)=QMIN
  GQ(1)=0.0
  
  CALL TRNSMIN(NQINT1,NQINT2,NQ2,Q2,Q3,Q4,Q5,X1,GX1)
  DO I=1,NQ2
   Q(I+NQ1)=X1(I)
   GQ(I+NQ1)=GX1(I)
  END DO
  GQ(NQ-1)=0.0
  GQ(NQ)=0.0
  Q(NQ-1)=Q(NQ-2)+0.4*(Q6-Q(NQ-2))
  Q(NQ)=Q6
  RETURN
 END SUBROUTINE QV4SETN
END MODULE meshpoints
