#include "fdefs.h"

module meshscatt
 use precision
 use parallel
 use mpi_const
 use constants
 use meshpoints
 use gauss
        
 implicit none
 private
 public meshscatt_init
 
 integer,public :: scatt_mode,enern=0,q0n=0
 real(dpreal),allocatable,public :: q0p(:),q0w(:)
 complex(dpreal),allocatable,public :: ener_kin(:),ener(:)
 
contains

 subroutine meshscatt_init(edeut)
  implicit none
 
  integer :: ip12,ip3,iq0,iener
  logical :: p12_custom,p3_custom
  real(dpreal) :: edeut,etemp,q0tilde

  ! make the routine reentry capable to be able to readjust energies 
  ! after deuteron has been calculated 
  
  if(allocated(ener_kin)) deallocate(ener_kin)
  if(allocated(ener)) deallocate(ener)
  if(allocated(q0w)) deallocate(q0w)
  if(allocated(q0p)) deallocate(q0p)
  
  ! Parameters are read from the file 'ener.dat'
  scatt=.true.
  open(unit=75,file='ener.dat',form='formatted',status='old')
  read(75,*,end=100,err=100) scatt_mode
  read(75,*,end=100,err=100) p12_custom,p3_custom
  read(75,*,end=100,err=100) enern
  allocate(ener_kin(enern))
  allocate(ener(enern))
  do iener=1,enern
   read(75,*,end=100,err=100) etemp
   ener_kin(iener)=etemp/hbarc ! energies read are in MeV
   if (etemp.le.0.0_dpreal) then
    scatt=.false.
   end if
  end do
  close(75)
  ! The number of q0 and ener points is equal
  q0n=enern
 
  ! set total energy to the kinetic energy plus the deuteron binding energy
  ener=ener_kin+edeut
 
  ! Momentum points for q0 are calculated
  ! q0w is now a dummy field required for testing routines only and set to zero here
  allocate(q0w(q0n))
  q0w=0.0_dpreal
  allocate(q0p(q0n))
  q0p=0.0_dpreal
  if (scatt) then
      
   if (scatt_mode.eq.3) then
    do iq0=1,q0n
     q0p(iq0)=sqrt(4.0_dpreal*ener_kin(iq0)*mnuc/3.0_dpreal)
    end do
   else if (scatt_mode.eq.431) then
    do iq0=1,q0n
     q0p(iq0)=sqrt(3.0_dpreal*ener_kin(iq0)*mnuc/2.0_dpreal)
    end do
   else if (scatt_mode.eq.422) then
    do iq0=1,q0n
     q0p(iq0)=sqrt(2.0_dpreal*ener_kin(iq0)*mnuc)
    end do
   else
    stop 'Scattering mode unknown'
   end if
  else
   write(*,*) 'Negative energies are present, therefore no scattering observables can be computed'
  end if ! scatt
  
  ! this variable is not used elsewhere except in p3_custom branch
  ! just define it work nicer printout later 
  q0tilde=0.0_dpreal
  
  ! DEFINITION OF CUSTOMIZED P3 GRID POINTS
  if (p3_custom) then
    ! Definition of q0tilde (uses first energy)
    q0tilde=sqrt(q0p(1)**2+4.0_dpreal*mnuc*edeut/3.0_dpreal)*0.999999

   ! Exact replica of previous 64 point p3 grid (variable number)
!   call gauleg(p3n/8-2,0.0_dpreal,0.4_dpreal*q0tilde,p3p(1:p3n/8-2),p3w(1:p3n/8-2))
!   call gauleg(p3n/8-1,0.4_dpreal*q0tilde,0.5_dpreal*q0tilde,p3p(p3n/8-1:p3n/4-3),p3w(p3n/8-1:p3n/4-3))
!   call gauleg(p3n/8,0.5_dpreal*q0tilde,0.6_dpreal*q0tilde,p3p(p3n/4-2:3*p3n/8-3),p3w(p3n/4-2:3*p3n/8-3))
!   call gauleg(p3n/8-1,0.6_dpreal*q0tilde,0.9_dpreal*q0tilde,p3p(3*p3n/8-2:p3n/2-4),p3w(3*p3n/8-2:p3n/2-4))
!   call gauleg(p3n/8-1,0.9_dpreal*q0tilde,q0tilde,p3p(p3n/2-3:5*p3n/8-5),p3w(p3n/2-3:5*p3n/8-5))
!   p3p(5*p3n/8-4)=q0tilde
!   p3w(5*p3n/8-4)=0.0_dpreal
!   call gauleg(p3n/8-1,q0tilde,1.1_dpreal*q0tilde,p3p(5*p3n/8-3:3*p3n/4-5),p3w(5*p3n/8-3:3*p3n/4-5))
!   call gauleg(p3n/8-1,1.1_dpreal*q0tilde,2.0_dpreal*q0tilde,p3p(3*p3n/4-4:7*p3n/8-6),p3w(3*p3n/4-4:7*p3n/8-6))
!   call gauleg(p3n/8+6,2.0_dpreal*q0tilde,p3p3,p3p(7*p3n/8-5:p3n),p3w(7*p3n/8-5:p3n))
   
   ! More general grid which is still close to the 64 point grid but more flexible for different numbers of p3 grid points
   call gauleg(p3n/8-1,0.0_dpreal,0.4_dpreal*q0tilde,p3p(1:p3n/8-1),p3w(1:p3n/8-1))
   call gauleg(p3n/8,0.4_dpreal*q0tilde,0.5_dpreal*q0tilde,p3p(p3n/8:p3n/4-1),p3w(p3n/8:p3n/4-1))
   call gauleg(p3n/8-1,0.5_dpreal*q0tilde,0.6_dpreal*q0tilde,p3p(p3n/4:3*p3n/8-2),p3w(p3n/4:3*p3n/8-2))
   call gauleg(p3n/8-1,0.6_dpreal*q0tilde,0.9_dpreal*q0tilde,p3p(3*p3n/8-1:p3n/2-3),p3w(3*p3n/8-1:p3n/2-3))
   call gauleg(p3n/8-1,0.9_dpreal*q0tilde,q0tilde,p3p(p3n/2-2:5*p3n/8-4),p3w(p3n/2-2:5*p3n/8-4))
   p3p(5*p3n/8-3)=q0tilde
   p3w(5*p3n/8-3)=0.0_dpreal
   call gauleg(p3n/8-1,q0tilde,1.1_dpreal*q0tilde,p3p(5*p3n/8-2:3*p3n/4-4),p3w(5*p3n/8-2:3*p3n/4-4))
   call gauleg(p3n/8,1.1_dpreal*q0tilde,2.0_dpreal*q0tilde,p3p(3*p3n/4-3:7*p3n/8-4),p3w(3*p3n/4-3:7*p3n/8-4))
   call gauleg(p3n/8+4,2.0_dpreal*q0tilde,p3p3,p3p(7*p3n/8-3:p3n),p3w(7*p3n/8-3:p3n))
 
   ! 70 points (fixed number)

!   call gauleg(19,0.0_dpreal,0.800511_dpreal,p3p(1:19),p3w(1:19))
!   call gauleg(29,0.800511_dpreal,1.200511_dpreal,p3p(20:48),p3w(20:48))
!   call gauleg(12,1.200511_dpreal,3.0_dpreal,p3p(49:60),p3w(49:60))
!   call gauleg(10,3.0_dpreal,20.0_dpreal,p3p(60:70),p3w(60:70))

   ! 80 points (fixed number)

!   call gauleg(19,0.0_dpreal,0.800511_dpreal,p3p(1:19),p3w(1:19))
!   call gauleg(31,0.800511_dpreal,1.200511_dpreal,p3p(20:50),p3w(20:50))
!   call gauleg(16,1.200511_dpreal,3.0_dpreal,p3p(51:66),p3w(51:66))
!   call gauleg(14,3.0_dpreal,20.0_dpreal,p3p(67:80),p3w(67:80))

   ! 80 points with \tilde{q}_0 as point with weight 0 (fixed number)

!   call gauleg(34,0.0_dpreal,1.000511_dpreal,p3p(1:34),p3w(1:34))
!   p3p(35)=1.000511_dpreal
!   p3w(35)=0.0_dpreal
!   call gauleg(31,1.000511_dpreal,3.0_dpreal,p3p(36:66),p3w(36:66))
!   call gauleg(14,3.0_dpreal,20.0_dpreal,p3p(67:80),p3w(67:80))

   ! 40 points with \tilde{q}_0 as point with weight 0 (fixed number)

!   call gauleg(9,0.0_dpreal,0.5_dpreal,p3p(1:9),p3w(1:9))
!   call gauleg(11,0.5_dpreal,1.000511_dpreal,p3p(10:20),p3w(10:20))
!   p3p(21)=1.000511_dpreal
!   p3w(21)=0.0_dpreal
!   call gauleg(11,1.000511_dpreal,2.0_dpreal,p3p(22:32),p3w(22:32))
!   call gauleg(8,2.0_dpreal,10.0_dpreal,p3p(33:40),p3w(33:40))

   ! 40 points with \tilde{q}_0 as point with weight 0 and more points at 0.5 and 1.0005 fm^-1 (fixed number)

!   call gauleg(4,0.0_dpreal,0.4_dpreal,p3p(1:4),p3w(1:4))
!   call gauleg(5,0.4_dpreal,0.5_dpreal,p3p(5:9),p3w(5:9))
!   call gauleg(6,0.5_dpreal,0.6_dpreal,p3p(10:15),p3w(10:15))
!   call gauleg(5,0.6_dpreal,0.9_dpreal,p3p(16:20),p3w(16:20))
!   call gauleg(5,0.9_dpreal,1.000511_dpreal,p3p(21:25),p3w(21:25))
!   p3p(26)=1.000511_dpreal
!   p3w(26)=0.0_dpreal
!   call gauleg(5,1.000511_dpreal,1.1_dpreal,p3p(27:31),p3w(27:31))
!   call gauleg(4,1.1_dpreal,2.0_dpreal,p3p(32:35),p3w(32:35))
!   call gauleg(5,2.0_dpreal,10.0_dpreal,p3p(36:40),p3w(36:40))

   ! 48 points with \tilde{q}_0 as point with weight 0 and more points at 0.5 and 1.0005 fm^-1 (fixed number)

!   call gauleg(4,0.0_dpreal,0.4_dpreal,p3p(1:4),p3w(1:4))
!   call gauleg(5,0.4_dpreal,0.5_dpreal,p3p(5:9),p3w(5:9))
!   call gauleg(6,0.5_dpreal,0.6_dpreal,p3p(10:15),p3w(10:15))
!   call gauleg(5,0.6_dpreal,0.9_dpreal,p3p(16:20),p3w(16:20))
!   call gauleg(5,0.9_dpreal,1.000511_dpreal,p3p(21:25),p3w(21:25))
!   p3p(26)=1.000511_dpreal
!   p3w(26)=0.0_dpreal
!   call gauleg(5,1.000511_dpreal,1.1_dpreal,p3p(27:31),p3w(27:31))
!   call gauleg(7,1.1_dpreal,2.0_dpreal,p3p(32:38),p3w(32:38))
!   call gauleg(10,2.0_dpreal,10.0_dpreal,p3p(39:48),p3w(39:48))
  
   ! 64 points with \tilde{q}_0 as point with weight 0 and more points at 0.5 and 1.0005 fm^-1 (fixed number)
  
!   q0tilde=1.00048_dpreal
!   call gauleg(6,0.0_dpreal,0.4_dpreal,p3p(1:6),p3w(1:6))
!   call gauleg(7,0.4_dpreal,0.5_dpreal,p3p(7:13),p3w(7:13))
!   call gauleg(8,0.5_dpreal,0.6_dpreal,p3p(14:21),p3w(14:21))
!   call gauleg(7,0.6_dpreal,0.9_dpreal,p3p(22:28),p3w(22:28))
!   call gauleg(7,0.9_dpreal,q0tilde,p3p(29:35),p3w(29:35))
!   p3p(36)=q0tilde
!   p3w(36)=0.0_dpreal
!   call gauleg(7,q0tilde,1.1_dpreal,p3p(37:43),p3w(37:43))
!   call gauleg(7,1.1_dpreal,2.0_dpreal,p3p(44:50),p3w(44:50))
!   call gauleg(14,2.0_dpreal,15.0_dpreal,p3p(51:64),p3w(51:64))
  
   ! 256 points with \tilde{q}_0 as point with weight 0 and more points at 0.5 and 1.0005 fm^-1 (fixed number)
  
!   q0tilde=1.00048_dpreal
!   call gauleg(24,0.0_dpreal,0.4_dpreal,p3p(1:24),p3w(1:24))
!   call gauleg(28,0.4_dpreal,0.5_dpreal,p3p(25:52),p3w(25:52))
!   call gauleg(32,0.5_dpreal,0.6_dpreal,p3p(53:84),p3w(53:84))
!   call gauleg(28,0.6_dpreal,0.9_dpreal,p3p(85:112),p3w(85:112))
!   call gauleg(30,0.9_dpreal,q0tilde,p3p(113:142),p3w(113:142))
!   p3p(143)=q0tilde
!   p3w(143)=0.0_dpreal
!   call gauleg(29,q0tilde,1.1_dpreal,p3p(144:172),p3w(144:172))
!   call gauleg(28,1.1_dpreal,2.0_dpreal,p3p(173:200),p3w(173:200))
!   call gauleg(56,2.0_dpreal,15.0_dpreal,p3p(201:256),p3w(201:256))

   ! 80 points with \tilde{q}_0 as point with weight 0 and more points near p3=0.5 fm^-1 (STANDARD) (fixed number)

!   call gauleg(22,0.0_dpreal,0.5_dpreal,p3p(1:22),p3w(1:22))
!   call gauleg(22,0.5_dpreal,1.000511_dpreal,p3p(23:44),p3w(23:44))
!   p3p(45)=1.000511_dpreal
!   p3w(45)=0.0_dpreal
!   call gauleg(16,1.000511_dpreal,2.0_dpreal,p3p(46:61),p3w(46:61))
!   call gauleg(19,2.0_dpreal,20.0_dpreal,p3p(62:80),p3w(62:80))

   ! 90 points with \tilde{q}_0 as point with weight 0 and more points near p3=0.5 fm^-1 (fixed number)

!   call gauleg(22,0.0_dpreal,0.5_dpreal,p3p(1:22),p3w(1:22))
!   call gauleg(22,0.5_dpreal,1.000511_dpreal,p3p(23:44),p3w(23:44))
!   p3p(45)=1.000511_dpreal
!   p3w(45)=0.0_dpreal
!   call gauleg(16,1.000511_dpreal,2.0_dpreal,p3p(46:61),p3w(46:61))
!   call gauleg(19,2.0_dpreal,8.0_dpreal,p3p(62:80),p3w(62:80))
!   call gauleg(10,8.0_dpreal,20.0_dpreal,p3p(81:90),p3w(81:90))

   ! 100 points with \tilde{q}_0 as point with weight 0 and more points near p3=0.5 fm^-1 (fixed number)

!   call gauleg(27,0.0_dpreal,0.5_dpreal,p3p(1:27),p3w(1:27))
!   call gauleg(27,0.5_dpreal,1.000511_dpreal,p3p(28:54),p3w(28:54))
!   p3p(55)=1.000511_dpreal
!   p3w(55)=0.0_dpreal
!   call gauleg(16,1.000511_dpreal,2.0_dpreal,p3p(56:71),p3w(56:71))
!   call gauleg(19,2.0_dpreal,8.0_dpreal,p3p(72:90),p3w(72:90))
!   call gauleg(10,8.0_dpreal,20.0_dpreal,p3p(91:100),p3w(91:100))

   ! 40 points only around p3=1.00051 fm^-1 and width 0.04 (fixed number)
  
!   call gauleg(19,0.934947_dpreal,1.00051_dpreal,p3p(1:19),p3w(1:19))
!   p3p(20)=1.00051_dpreal
!   p3w(20)=0.0_dpreal
!   call gauleg(20,1.00051_dpreal,1.06204_dpreal,p3p(21:40),p3w(21:40))
  
   ! 40 points only around p3=1.00051 fm^-1 and width 0.01 (fixed number)
  
!   call gauleg(19,0.98453_dpreal,1.0005114778209308_dpreal,p3p(1:19),p3w(1:19))
!   p3p(20)=1.0005114778209308_dpreal
!   p3w(20)=0.0_dpreal
!   call gauleg(20,1.0005114778209308_dpreal,1.01624_dpreal,p3p(21:40),p3w(21:40))
  
   ! 40 points only around p3=1.00051 fm^-1 and width 0.002 (fixed number)
  
!   call gauleg(19,0.997336_dpreal,1.00051_dpreal,p3p(1:19),p3w(1:19))
!   p3p(20)=1.00051_dpreal
!   p3w(20)=0.0_dpreal
!   call gauleg(20,1.00051_dpreal,1.00368_dpreal,p3p(21:40),p3w(21:40))
  end if ! p3_custom
  
  
  ! DEFINITION OF CUSTOMIZED P12 GRID POINTS
  if (p12_custom) then
  
  ! 70 points (fixed number)

!  call gauleg(19,0.0_dpreal,0.800511_dpreal,p12p(1:19),p12w(1:19))
!  call gauleg(29,0.800511_dpreal,1.200511_dpreal,p12p(20:48),p12w(20:48))
!  call gauleg(12,1.200511_dpreal,3.0_dpreal,p12p(49:60),p12w(49:60))
!  call gauleg(10,3.0_dpreal,20.0_dpreal,p12p(60:70),p12w(60:70))
  
  end if ! p12_custom
 
  ! START TEST: MANUALLY SETTING POINTS FOR Q4, P34 AND Q

!  q4p=0.0_dpreal
!  q4w=0.0_dpreal
  
!  p34p=0.0_dpreal
!  p34w=0.0_dpreal
  
!  qp=0.0_dpreal
!  qw=0.0_dpreal
  
  ! END TEST: MANUALLY SETTING POINTS FOR Q4, P34 AND Q
  
  if (master) then
   if (scatt_mode.eq.3) then
    write(*,*) 'scattering mode: 3N'
   else if (scatt_mode.eq.431) then
    write(*,*) 'scattering mode: 4N 3+1'
   else if (scatt_mode.eq.422) then
    write(*,*) 'scattering mode: 4N 2+2'
   end if
   write(*,*)
   do iq0=1,q0n
    write(*,*) 'momentum point q0',iq0,'(in fm^-1)=',q0p(iq0)
   end do
   write(*,*)
   do iener=1,enern
    write(*,*) 'four-body kinetic energy E',iener,'(in MeV)=',ener_kin(iener)*hbarc
   end do
   write(*,*)
   do iq0=1,q0n
    write(*,*) 'four-body total energy E',iq0,'(in MeV)=',ener(iq0)*hbarc
   end do
   write(*,*)
   write(*,*) 'momentum point q0tilde (in fm^-1)=',q0tilde
   write(*,*)
   if (p12_custom) then
    write(*,*) 'custom p12 points were set:'
    do ip12=1,p12n
     write(*,*) 'new momentum point p12',ip12,'(in fm^-1)=',p12p(ip12)
    end do
    write(*,*)
   end if
   if (p3_custom) then
    write(*,*) 'custom p3 points were set:'
    do ip3=1,p3n
     write(*,*) 'new momentum point p3',ip3,'(in fm^-1)=',p3p(ip3)
    end do
   end if
  end if
  
  return
  ! Regular end of subroutine
100 continue
  ! Treatment of i/o errors
  stop 'problem with ener.dat'
 end subroutine meshscatt_init
    
end module meshscatt