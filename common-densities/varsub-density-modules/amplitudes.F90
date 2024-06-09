#include "fdefs.h"

!! module amplitudes requires initialized meshpoints,parallel and pwbook
!! and provides types and basic operations for the NN,3N, 4N31 and 4N22 amplitudes  
!! module assumes that momenta and channels are distributed over the pertinent communicators

MODULE amplitudes
 USE precision
 USE pwbook
 USE parallel 
 USE meshpoints
 use meshscatt
 USE spline
 USE mpi_const
 USE hdf_tool
 USE HDF5
 USE constants
 USE clebsch
 USE besselfu
 USE gauss
 
! use threebodyscatt

 PRIVATE 

 ! definition of types for the amplitudes 
 TYPE AMPNN
  SEQUENCE
  COMPLEX(spreal),ALLOCATABLE :: amp(:,:)
 END TYPE AMPNN

 TYPE AMP3N
  SEQUENCE
  COMPLEX(spreal),ALLOCATABLE :: amp(:,:,:)
 END TYPE AMP3N

 TYPE AMP4N31
  SEQUENCE
  COMPLEX(spreal),ALLOCATABLE :: amp(:,:,:,:)
 END TYPE AMP4N31

 TYPE AMP4N22
  SEQUENCE
  COMPLEX(spreal),ALLOCATABLE :: amp(:,:,:,:)
 END TYPE AMP4N22

 TYPE RAMPNN
  SEQUENCE
  COMPLEX(spreal),ALLOCATABLE :: amp(:,:)
 END TYPE RAMPNN

 TYPE RAMP3N
  SEQUENCE
  COMPLEX(spreal),ALLOCATABLE :: amp(:,:,:)
 END TYPE RAMP3N

 TYPE RAMP4N31
  SEQUENCE
  COMPLEX(spreal),ALLOCATABLE :: amp(:,:,:,:)
 END TYPE RAMP4N31

 TYPE RAMP4N22
  SEQUENCE
  COMPLEX(spreal),ALLOCATABLE :: amp(:,:,:,:)
 END TYPE RAMP4N22

 PUBLIC AMPNN,AMP3N,AMP4N31,AMP4N22
 PUBLIC RAMPNN,RAMP3N,RAMP4N31,RAMP4N22
 
 ! definition of interface routines for the standard operators 

 ! product is either skalar * amp 
 ! or amp * amp 
 ! scalars can be spreal,dpreal,spcmplx,dpcmplx

 INTERFACE normalize
   MODULE PROCEDURE norm_nn,norm_3n,norm_4n31,norm_4n22,norm_4nboth
 END INTERFACE

 INTERFACE OPERATOR (*) 
  MODULE PROCEDURE dotprod_nn,prodsprealnn,proddprealnn,prodspcmplxnn,proddpcmplxnn
  MODULE PROCEDURE dotprod_3n,prodspreal3n,proddpreal3n,prodspcmplx3n,proddpcmplx3n
  MODULE PROCEDURE dotprod_4n31,prodspreal4n31,proddpreal4n31,prodspcmplx4n31,proddpcmplx4n31
  MODULE PROCEDURE dotprod_4n22,prodspreal4n22,proddpreal4n22,prodspcmplx4n22,proddpcmplx4n22
  MODULE PROCEDURE rdotprod_nn,rprodsprealnn,rproddprealnn,rprodspcmplxnn,rproddpcmplxnn
  MODULE PROCEDURE rdotprod_3n,rprodspreal3n,rproddpreal3n,rprodspcmplx3n,rproddpcmplx3n
  MODULE PROCEDURE rdotprod_4n31,rprodspreal4n31,rproddpreal4n31,rprodspcmplx4n31,rproddpcmplx4n31
  MODULE PROCEDURE rdotprod_4n22,rprodspreal4n22,rproddpreal4n22,rprodspcmplx4n22,rproddpcmplx4n22
 END INTERFACE
 
 ! assignment of amp = amp or amp=scalar 
 ! scalars can be spreal,dpreal,spcmplx,dpcmplx

 INTERFACE ASSIGNMENT (=) 
  MODULE PROCEDURE set_psipsinn,set_psisprealnn,set_psidprealnn,set_psispcmplxnn,set_psidpcmplxnn
  MODULE PROCEDURE set_psipsi3n,set_psispreal3n,set_psidpreal3n,set_psispcmplx3n,set_psidpcmplx3n
  MODULE PROCEDURE set_psipsi4n31,set_psispreal4n31,set_psidpreal4n31,set_psispcmplx4n31,set_psidpcmplx4n31
  MODULE PROCEDURE set_psipsi4n22,set_psispreal4n22,set_psidpreal4n22,set_psispcmplx4n22,set_psidpcmplx4n22
  MODULE PROCEDURE rset_psipsinn,rset_psisprealnn,rset_psidprealnn,rset_psispcmplxnn,rset_psidpcmplxnn
  MODULE PROCEDURE rset_psipsi3n,rset_psispreal3n,rset_psidpreal3n,rset_psispcmplx3n,rset_psidpcmplx3n
  MODULE PROCEDURE rset_psipsi4n31,rset_psispreal4n31,rset_psidpreal4n31,rset_psispcmplx4n31,rset_psidpcmplx4n31
  MODULE PROCEDURE rset_psipsi4n22,rset_psispreal4n22,rset_psidpreal4n22,rset_psispcmplx4n22,rset_psidpcmplx4n22
  MODULE PROCEDURE fourier_r2p_psinn,fourier_r2p_psi3n,fourier_r2p_psi4n31,fourier_r2p_psi4n22
  MODULE PROCEDURE fourier_p2r_psinn,fourier_p2r_psi3n,fourier_p2r_psi4n31,fourier_p2r_psi4n22
 END INTERFACE

 ! addition and subtraction of amplitude
 INTERFACE OPERATOR (+) 
  MODULE PROCEDURE plus_psinn,plus_psi3n,plus_psi4n31,plus_psi4n22
  MODULE PROCEDURE rplus_psinn,rplus_psi3n,rplus_psi4n31,rplus_psi4n22
 END INTERFACE

 INTERFACE OPERATOR (-) 
  MODULE PROCEDURE minus_psinn,minus_psi3n,minus_psi4n31,minus_psi4n22
  MODULE PROCEDURE rminus_psinn,rminus_psi3n,rminus_psi4n31,rminus_psi4n22
 END INTERFACE

 
 ! to collect all points of one channel
 INTERFACE collect_channel
  MODULE PROCEDURE collect_channel_nn,collect_channel_3n
 END INTERFACE  
 
 ! to write the amplitudes 
 INTERFACE writeamp
  MODULE PROCEDURE writeamp_nn,writeamp_3n,writeamp_4n31,writeamp_4n22
 END INTERFACE 
 
 
 ! to read the amplitudes 
 INTERFACE readamp
  MODULE PROCEDURE readamp_nn,readamp_3n,readamp_4n31,readamp_4n22
 END INTERFACE 
 
 INTERFACE printinfo
  MODULE PROCEDURE printinfo_3n_wave,printinfo_3n,print_info_nn,&
  &                printinfo_4n,printinfo_4n31,printinfo_4n22
 END INTERFACE
 
 INTERFACE project_LS
  MODULE PROCEDURE project_LS_3n 
 END INTERFACE
 
 INTERFACE project_principalS
  MODULE PROCEDURE project_principalS_3n 
 END INTERFACE
 
 INTERFACE printmomdist
  MODULE PROCEDURE printmomdist_3n
 END INTERFACE
 
 INTERFACE printmomcorr
  MODULE PROCEDURE printmomcorr_3n
 END INTERFACE
 
 INTERFACE testamp
  MODULE PROCEDURE testamp_nn,testamp_3n,testamp_4n31,testamp_4n22 
 END INTERFACE 
 
 INTERFACE printrms
  MODULE PROCEDURE printrms_nn,printrms_3n,printrms_4n31,printrms_4n22
 END INTERFACE printrms
 
 PUBLIC OPERATOR(*),OPERATOR(+),OPERATOR(-),ASSIGNMENT(=),initamp,&
 &      collect_channel_3n,normalize,printmomcorr,printmomdist,testamp,printrms
 public write_thatphi_3N,collect_channel,writeamp,readamp,printinfo,project_LS,project_principalS,printamp
 
 ! variables needed for read and write subroutines

 integer :: p12n_read,p3n_read,alpha3ncdepmax_read
 integer,allocatable :: alpha3n_index(:),indx_p12_read_p12(:,:),indx_p3_read_p3(:,:)
 real(dpreal),allocatable :: spl_p12_read_p12(:,:),spl_p3_read_p3(:,:)

 public open_writefile_thatphi_3N,writefile_thatphi_3N,open_readfile_thatphi_3N,readfile_thatphi_3N

 ! integers that define the local components of the amplitudes 

 INTEGER,PUBLIC :: myp12,mynp12          ! first element of local block and number of elements 
 INTEGER,PUBLIC :: myp3,mynp3
 INTEGER,PUBLIC :: myq4,mynq4
 INTEGER,PUBLIC :: myp34,mynp34
 INTEGER,PUBLIC :: myq,mynq
 INTEGER,PUBLIC :: myq0,mynq0
 
 INTEGER,PUBLIC :: mynalphaNN    ! for (2 body) - number of elements, each npe_alpha element is local starting from myalphaid+1
 INTEGER,PUBLIC :: myalphaNN     ! newly defined and, as stated above, always equal to myalphaid+1
 INTEGER,PUBLIC :: mynalpha3N    ! for (3 body) - number of elements, each npe_alpha element is local starting from myalphaid+1
 INTEGER,PUBLIC :: myalpha3N     ! newly defined and, as stated above, always equal to myalphaid+1
 INTEGER,PUBLIC :: mynalpha4N31  ! for (3+1 body) - number of elements, each npe_alpha element is local starting from myalphaid+1
 INTEGER,PUBLIC :: myalpha4N31   ! newly defined and, as stated above, always equal to myalphaid+1
 INTEGER,PUBLIC :: mynbeta4N22   ! for (2+2 body) - number of elements, each npe_alpha element is local starting from myalphaid+1
 INTEGER,PUBLIC :: mybeta4N22    ! newly defined and, as stated above, always equal to myalphaid+1

 INTEGER,ALLOCATABLE,PUBLIC :: firstp12(:),lastp12(:),nump12(:) ! same as above for pe=0,npe_p12-1 
 INTEGER,ALLOCATABLE,PUBLIC :: firstp3(:),lastp3(:),nump3(:) ! same as above for pe=0,npe_p3-1
 INTEGER,ALLOCATABLE,PUBLIC :: firstq4(:),lastq4(:),numq4(:) ! same as above for pe=0,npe_q4-1

 INTEGER,ALLOCATABLE,PUBLIC :: firstp34(:),lastp34(:),nump34(:) ! same as above for pe=0,npe_p3-1
 INTEGER,ALLOCATABLE,PUBLIC :: firstq(:),lastq(:),numq(:) ! same as above for pe=0,npe_q4-1

 INTEGER,ALLOCATABLE,PUBLIC :: firstq0(:),lastq0(:),numq0(:) ! same as above for pe=0,npe_ener-1

 INTEGER,ALLOCATABLE,PUBLIC :: numalphaNN(:),numalpha3N(:),numalpha4N31(:),numbeta4N22(:) ! same for pe=0,npe_alpha-1

 
 ! integers that define the local components of the amplitudes 

 INTEGER,PUBLIC :: myr12,mynr12          ! first element of local block and number of elements 
 INTEGER,PUBLIC :: myr3,mynr3
 INTEGER,PUBLIC :: myr4,mynr4
 INTEGER,PUBLIC :: myr34,mynr34
 INTEGER,PUBLIC :: myr,mynr

 
 INTEGER,ALLOCATABLE,PUBLIC :: firstr12(:),lastr12(:),numr12(:) ! same as above for pe=0,npe_p12-1 
 INTEGER,ALLOCATABLE,PUBLIC :: firstr3(:),lastr3(:),numr3(:) ! same as above for pe=0,npe_p3-1
 INTEGER,ALLOCATABLE,PUBLIC :: firstr4(:),lastr4(:),numr4(:) ! same as above for pe=0,npe_q4-1

 INTEGER,ALLOCATABLE,PUBLIC :: firstr34(:),lastr34(:),numr34(:) ! same as above for pe=0,npe_p3-1
 INTEGER,ALLOCATABLE,PUBLIC :: firstr(:),lastr(:),numr(:) ! same as above for pe=0,npe_q4-1

 ! precomputed integration weights 
 
 REAL(spreal),ALLOCATABLE :: p12weight(:),p12p3weight(:,:),p12p3q4weight(:,:,:),p12p34qweight(:,:,:)
 REAL(spreal),ALLOCATABLE :: r12weight(:),r12r3weight(:,:),r12r3r4weight(:,:,:),r12r34rweight(:,:,:)
 
 ! grid points for Fourier trafo with  Bessel functions 
 REAL(dpreal) :: rintp1,rintp2,rintp3
 REAL(dpreal),ALLOCATABLE :: rintp(:),rintw(:) 
 INTEGER  :: rintn1,rintn2,rintn
 
 ! grid points for Fourier trafo with  Bessel functions 
 REAL(dpreal) :: pintp1,pintp2,pintp3
 REAL(dpreal),ALLOCATABLE :: pintp(:),pintw(:) 
 INTEGER  :: pintn1,pintn2,pintn
 
 ! arrays that contain the splines for the calculations of Fourier trafos 
 REAL(dpreal),ALLOCATABLE :: r12spl(:,:),r3spl(:,:),r4spl(:,:),r34spl(:,:),rspl(:,:)
 INTEGER,ALLOCATABLE :: r12indx(:,:),r3indx(:,:),r4indx(:,:),r34indx(:,:),rindx(:,:)
 REAL(dpreal),ALLOCATABLE :: p12spl(:,:),p3spl(:,:),q4spl(:,:),p34spl(:,:),qspl(:,:)
 INTEGER,ALLOCATABLE :: p12indx(:,:),p3indx(:,:),q4indx(:,:),p34indx(:,:),qindx(:,:)
 
 ! arrays that contain integrated bessel functions 
 REAL(spreal),ALLOCATABLE  :: besslr2p12(:,:,:)    !   (mynp12,r12n,0:l12max)
 REAL(spreal),ALLOCATABLE  :: besslr2p3(:,:,:)     !   (mynp3,r3n,0:l3max)
 REAL(spreal),ALLOCATABLE  :: besslr2p4(:,:,:)     !   (mynq4,r4n,0:l4max)
 REAL(spreal),ALLOCATABLE  :: besslr2p34(:,:,:)    !   (mynp34,r34n,0:l12max)
 REAL(spreal),ALLOCATABLE  :: besslr2p(:,:,:)      !   (mynq,rn,0:lammax)
 
 REAL(spreal),ALLOCATABLE  :: besslp2r12(:,:,:)    !   (mynr12,p12n,0:l12max)
 REAL(spreal),ALLOCATABLE  :: besslp2r3(:,:,:)     !   (mynr3,p3n,0:l3max)
 REAL(spreal),ALLOCATABLE  :: besslp2r4(:,:,:)     !   (mynr4,q4n,0:l4max)
 REAL(spreal),ALLOCATABLE  :: besslp2r34(:,:,:)    !   (mynr34,p34n,0:l12max)
 REAL(spreal),ALLOCATABLE  :: besslp2r(:,:,:)      !   (mynr,qn,0:lammax)
 
 ! timing of Bessel functions 
 REAL(dpreal),PUBLIC :: besstime=0.0,bessinttime=0.0,bessredtime=0.0,bessspltime=0.0
 REAL(dpreal),PUBLIC :: readamptime=0.0,writeamptime=0.0
 
 CHARACTER(LEN=4) :: rmeshtype
 
 REAL(dpreal) :: tolerance=1E-4_dpreal
 
CONTAINS

 ! initialization routine needs to be called before using amplitudes 
 ! the routine determines the size of amplitude arrays based on meshpoints 
 ! and number of processors
 
 SUBROUTINE initamp
  IMPLICIT NONE 
  INTEGER id
  INTEGER ip12,ip3,iq4,ip34,iq
  INTEGER ir12,ir3,ir4,ir34,ir
  
  ! first printout version running
  IF(master) THEN
   WRITE(*,*) 'Initamp version: ',VERREV
!   WRITE(*,*) 'Date           : ',VERDATE
  END IF
  
  ! first prepare p-space amplitudes 
  ! distribution of P12 mesh 
  
  IF(allocated(firstp12)) deallocate(firstp12)
  IF(allocated(lastp12)) deallocate(lastp12)
  IF(allocated(nump12)) deallocate(nump12)
  
  ALLOCATE(firstp12(0:npe_p12-1),lastp12(0:npe_p12-1),nump12(0:npe_p12-1))
  
  DO id=0,npe_p12-1
   CALL distr_block(P12N,npe_p12,id,firstp12(id),lastp12(id),nump12(id))   
  END DO
  myp12=firstp12(myp12id)
  mynp12=nump12(myp12id)
  
  ! distribution of P3 mesh 
  
  IF(allocated(firstp3)) deallocate(firstp3)
  IF(allocated(lastp3)) deallocate(lastp3)
  IF(allocated(nump3)) deallocate(nump3)
  
  ALLOCATE(firstp3(0:npe_p3-1),lastp3(0:npe_p3-1),nump3(0:npe_p3-1))
  
  DO id=0,npe_p3-1
   CALL distr_block(P3N,npe_p3,id,firstp3(id),lastp3(id),nump3(id))   
  END DO
  myp3=firstp3(myp3id)
  mynp3=nump3(myp3id)
  
  ! distribution of Q4 mesh 
  
  IF(allocated(firstq4)) deallocate(firstq4)
  IF(allocated(lastq4)) deallocate(lastq4)
  IF(allocated(numq4)) deallocate(numq4)
  
  ALLOCATE(firstq4(0:npe_q4-1),lastq4(0:npe_q4-1),numq4(0:npe_q4-1))
  
  DO id=0,npe_q4-1
   CALL distr_block(Q4N,npe_q4,id,firstq4(id),lastq4(id),numq4(id))   
  END DO
  myq4=firstq4(myq4id)
  mynq4=numq4(myq4id)
  
  ! distribution of P34 mesh 
  
  IF(allocated(firstp34)) deallocate(firstp34)
  IF(allocated(lastp34)) deallocate(lastp34)
  IF(allocated(nump34)) deallocate(nump34)
  
  ALLOCATE(firstp34(0:npe_p3-1),lastp34(0:npe_p3-1),nump34(0:npe_p3-1))
  
  DO id=0,npe_p3-1
   CALL distr_block(P34N,npe_p3,id,firstp34(id),lastp34(id),nump34(id))   
  END DO
  myp34=firstp34(myp3id)
  mynp34=nump34(myp3id)
  
  ! distribution of Q mesh 
  
  IF(allocated(firstq)) deallocate(firstq)
  IF(allocated(lastq)) deallocate(lastq)
  IF(allocated(numq)) deallocate(numq)
  
  ALLOCATE(firstq(0:npe_q4-1),lastq(0:npe_q4-1),numq(0:npe_q4-1))
  
  DO id=0,npe_q4-1
   CALL distr_block(QN,npe_q4,id,firstq(id),lastq(id),numq(id))   
  END DO
  myq=firstq(myq4id)
  mynq=numq(myq4id)
  
  ! distribution of Q0 mesh for energies (set of read in in meshpoints)
  
  IF(allocated(firstq0)) deallocate(firstq0)
  IF(allocated(lastq0)) deallocate(lastq0)
  IF(allocated(numq0)) deallocate(numq0)
  
  ALLOCATE(firstq0(0:npe_ener-1),lastq0(0:npe_ener-1),numq0(0:npe_ener-1))
  
  DO id=0,npe_ener-1
   CALL distr_block(Q0N,npe_ener,id,firstq0(id),lastq0(id),numq0(id))
  END DO
  myq0=firstq0(myenerid)
  mynq0=numq0(myenerid)
  
  ! distribution of pw channels 
  ! NN case 
  IF(allocated(numalphaNN)) deallocate(numalphaNN)
  ALLOCATE(numalphaNN(0:npe_alpha-1))
  DO id=0,npe_alpha-1
   CALL distr_piece(alphaNNcdepmax,npe_alpha,id,numalphaNN(id))
  END DO
  myalphaNN=myalphaid+1
  mynalphaNN=numalphaNN(myalphaid)
  
  ! 3N case 
  IF(allocated(numalpha3N)) deallocate(numalpha3N)
  ALLOCATE(numalpha3N(0:npe_alpha-1))
  DO id=0,npe_alpha-1
   CALL distr_piece(alpha3Ncdepmax,npe_alpha,id,numalpha3N(id))
  END DO
  myalpha3N=myalphaid+1
  mynalpha3N=numalpha3N(myalphaid)
  
  ! 4N31 case 
  IF(allocated(numalpha4N31)) deallocate(numalpha4N31)
  ALLOCATE(numalpha4N31(0:npe_alpha-1))
  DO id=0,npe_alpha-1
   CALL distr_piece(alpha4N31cdepmax,npe_alpha,id,numalpha4N31(id))
  END DO
  myalpha4N31=myalphaid+1
  mynalpha4N31=numalpha4N31(myalphaid)
  
  ! 4N22 case 
  IF(allocated(numbeta4N22)) deallocate(numbeta4N22)
  ALLOCATE(numbeta4N22(0:npe_alpha-1))
  DO id=0,npe_alpha-1
   CALL distr_piece(beta4N22cdepmax,npe_alpha,id,numbeta4N22(id))
  END DO
  mybeta4N22=myalphaid+1
  mynbeta4N22=numbeta4N22(myalphaid)
  
  ! then prepare r-space amplitudes 
  
  IF(allocated(firstr12)) deallocate(firstr12)
  IF(allocated(lastr12)) deallocate(lastr12)
  IF(allocated(numr12)) deallocate(numr12)
  
  ALLOCATE(firstr12(0:npe_p12-1),lastr12(0:npe_p12-1),numr12(0:npe_p12-1))
  
  DO id=0,npe_p12-1
   CALL distr_block(R12N,npe_p12,id,firstr12(id),lastr12(id),numr12(id))   
  END DO
  myr12=firstr12(myp12id)
  mynr12=numr12(myp12id)
  
  ! distribution of R3 mesh 
  
  IF(allocated(firstr3)) deallocate(firstr3)
  IF(allocated(lastr3)) deallocate(lastr3)
  IF(allocated(numr3)) deallocate(numr3)
  
  ALLOCATE(firstr3(0:npe_p3-1),lastr3(0:npe_p3-1),numr3(0:npe_p3-1))
  
  DO id=0,npe_p3-1
   CALL distr_block(R3N,npe_p3,id,firstr3(id),lastr3(id),numr3(id))   
  END DO
  myr3=firstr3(myp3id)
  mynr3=numr3(myp3id)
  
  ! distribution of Q4 mesh 
  
  IF(allocated(firstr4)) deallocate(firstr4)
  IF(allocated(lastr4)) deallocate(lastr4)
  IF(allocated(numr4)) deallocate(numr4)
  
  ALLOCATE(firstr4(0:npe_q4-1),lastr4(0:npe_q4-1),numr4(0:npe_q4-1))
  
  DO id=0,npe_q4-1
   CALL distr_block(R4N,npe_q4,id,firstr4(id),lastr4(id),numr4(id))   
  END DO
  myr4=firstr4(myq4id)
  mynr4=numr4(myq4id)
  
  ! distribution of P34 mesh 
  
  IF(allocated(firstr34)) deallocate(firstr34)
  IF(allocated(lastr34)) deallocate(lastr34)
  IF(allocated(numr34)) deallocate(numr34)
  
  ALLOCATE(firstr34(0:npe_p3-1),lastr34(0:npe_p3-1),numr34(0:npe_p3-1))
  
  DO id=0,npe_p3-1
   CALL distr_block(R34N,npe_p3,id,firstr34(id),lastr34(id),numr34(id))   
  END DO
  myr34=firstr34(myp3id)
  mynr34=numr34(myp3id)
  
  ! distribution of Q mesh 
  
  IF(allocated(firstr)) deallocate(firstr)
  IF(allocated(lastr)) deallocate(lastr)
  IF(allocated(numr)) deallocate(numr)
  
  ALLOCATE(firstr(0:npe_q4-1),lastr(0:npe_q4-1),numr(0:npe_q4-1))
  
  DO id=0,npe_q4-1
   CALL distr_block(RN,npe_q4,id,firstr(id),lastr(id),numr(id))   
  END DO
  myr=firstr(myq4id)
  mynr=numr(myq4id)
  
  
  ! define weights for different integrations (locally only)
  
  IF(allocated(p12weight)) DEALLOCATE(p12weight)
  IF(allocated(p12p3weight)) DEALLOCATE(p12p3weight)
  IF(allocated(p12p3q4weight)) DEALLOCATE(p12p3q4weight)
  IF(allocated(p12p34qweight)) DEALLOCATE(p12p34qweight)
  ALLOCATE(p12weight(mynp12))
  ALLOCATE(p12p3weight(mynp12,mynp3))
  ALLOCATE(p12p3q4weight(mynp12,mynp3,mynq4))
  ALLOCATE(p12p34qweight(mynp12,mynp34,mynq))
  
  DO ip12=1,mynp12
   p12weight(ip12)=P12W(ip12+myp12-1)*P12P(ip12+myp12-1)**2
   DO ip3=1,mynp3
    p12p3weight(ip12,ip3)=p12weight(ip12)*P3W(ip3+myp3-1)*P3P(ip3+myp3-1)**2
    DO iq4=1,mynq4
     p12p3q4weight(ip12,ip3,iq4)=p12p3weight(ip12,ip3)*Q4W(iq4+myq4-1)*Q4P(iq4+myq4-1)**2
    END DO
   END DO
   
   DO ip34=1,mynp34
    DO iq=1,mynq
     p12p34qweight(ip12,ip34,iq)=p12weight(ip12) &
          *P34W(ip34+myp34-1)*P34P(ip34+myp34-1)**2 &
          *QW(iq+myq-1)*QP(iq+myq-1)**2 
    END DO
   END DO
  END DO
  
! define weights for different r-integrations  (locally only)
  
  IF(allocated(r12weight)) DEALLOCATE(r12weight)
  IF(allocated(r12r3weight)) DEALLOCATE(r12r3weight)
  IF(allocated(r12r3r4weight)) DEALLOCATE(r12r3r4weight)
  IF(allocated(r12r34rweight)) DEALLOCATE(r12r34rweight)
  ALLOCATE(r12weight(mynr12))
  ALLOCATE(r12r3weight(mynr12,mynr3))
  ALLOCATE(r12r3r4weight(mynr12,mynr3,mynr4))
  ALLOCATE(r12r34rweight(mynr12,mynr34,mynr))
  
  DO ir12=1,mynr12
   r12weight(ir12)=R12W(ir12+myr12-1)*R12P(ir12+myr12-1)**2
   DO ir3=1,mynr3
    r12r3weight(ir12,ir3)=r12weight(ir12)*R3W(ir3+myr3-1)*R3P(ir3+myr3-1)**2
    DO ir4=1,mynr4
     r12r3r4weight(ir12,ir3,ir4)=r12r3weight(ir12,ir3)*R4W(ir4+myr4-1)*R4P(ir4+myr4-1)**2
    END DO
   END DO
   
   DO ir34=1,mynr34
    DO ir=1,mynr
     r12r34rweight(ir12,ir34,ir)=r12weight(ir12) &
          *R34W(ir34+myr34-1)*R34P(ir34+myr34-1)**2 &
          *RW(ir+myr-1)*RP(ir+myr-1)**2 
    END DO
   END DO
  END DO

  
  ! read in amp parameters
  ! this is here only the mesh points for the integration of besselfunctions 
  
  open(unit=75,file='amplitudes.dat',form='formatted',status='old')

  ! Read mesh point distribution type
  read(75,*,end=100,err=100) rmeshtype

  SELECT CASE(trim(rmeshtype))
   CASE('TRNS')  ! Read mesh parameters for TRNS 
    read(75,*,end=100,err=100) rintp1,rintp2,rintp3,rintn1,rintn2
    read(75,*,end=100,err=100) pintp1,pintp2,pintp3,pintn1,pintn2
    
    rintn=rintn1+rintn2
    pintn=pintn1+pintn2
    ! allocate grid point arrays 
    allocate(rintp(rintn),rintw(rintn))
    allocate(pintp(pintn),pintw(pintn))
    
    ! Integration points are calculated with subroutine "trns"
    call trns(rintn1,rintn2,rintn,rintp1,rintp2,rintp3,rintp,rintw)
    call trns(pintn1,pintn2,pintn,pintp1,pintp2,pintp3,pintp,pintw)
    
   CASE DEFAULT
    STOP 'meshtype not defined'
  END SELECT   
  
  read(75,*,end=100,err=100) tolerance   !  bound for interpolation of low mom part 
  close(75)
  
  return
  ! Regular end of subroutine
100 continue
  ! Treatment of i/o errors
  stop 'problem with amplitudes.dat'
  
  
 END SUBROUTINE initamp
 
 ! printing routine prints out the distribution of mesh points for 
 ! the amplitudes 
 ! is called by initamp on first processor only 
 
 SUBROUTINE printamp
  IMPLICIT NONE 
  INTEGER id12,id3,id4,idalpha,irint,ipint
  
  WRITE(*,*) 
  WRITE(*,*) 'Initialization of amplitudes module' 
  WRITE(*,*)
        ! first printout version running
  WRITE(*,*) 'Amplitudes version: ',VERREV
!  WRITE(*,*) 'Date              : ',VERDATE
  WRITE(*,*)
  WRITE(*,*) 'Distribution of mesh points' 
  DO id12=0,npe_p12-1
   WRITE(*,'(A,4I5)') 'ID12MESH: ',id12,firstp12(id12),lastp12(id12),nump12(id12)
  END DO
  DO id3=0,npe_p3-1
   WRITE(*,'(A,4I5)') 'ID3MESH:  ',id3,firstp3(id3),lastp3(id3),nump3(id3)
  END DO
  DO id4=0,npe_q4-1
   WRITE(*,'(A,4I5)') 'ID4MESH:  ',id4,firstq4(id4),lastq4(id4),numq4(id4)
  END DO
  DO id3=0,npe_p3-1
   WRITE(*,'(A,4I5)') 'ID34MESH: ',id3,firstp34(id3),lastp34(id3),nump34(id3)
  END DO
  DO id4=0,npe_q4-1 
   WRITE(*,'(A,4I5)') 'IDMESH:   ',id4,firstq(id4),lastq(id4),numq(id4)
  END DO
  
  DO id12=0,npe_p12-1
   WRITE(*,'(A,4I5)') 'IDR12MESH: ',id12,firstr12(id12),lastr12(id12),numr12(id12)
  END DO
  DO id3=0,npe_p3-1
   WRITE(*,'(A,4I5)') 'IDR3MESH:  ',id3,firstr3(id3),lastr3(id3),numr3(id3)
  END DO
  DO id4=0,npe_q4-1
   WRITE(*,'(A,4I5)') 'IDR4MESH:  ',id4,firstr4(id4),lastr4(id4),numr4(id4)
  END DO
  DO id3=0,npe_p3-1
   WRITE(*,'(A,4I5)') 'IDR34MESH:  ',id3,firstr3(id3),lastr3(id3),numr3(id3)
  END DO
  DO id4=0,npe_q4-1
   WRITE(*,'(A,4I5)') 'IDRMESH:  ',id4,firstr4(id4),lastr4(id4),numr4(id4)
  END DO
  
  WRITE(*,*) 'Distribution of alpha channels'
  DO idalpha=0,npe_alpha-1
   WRITE(*,'(A,5I5)') 'IDALPHA:  ',idalpha,numalphaNN(idalpha),numalpha3N(idalpha),& 
        numalpha4N31(idalpha),numbeta4N22(idalpha)
  END DO
  
  WRITE(*,*) 
  WRITE(*,*) 
  WRITE(*,*) 'Mesh points for Fourier trafo'
  WRITE(*,'(A,3I5)')    'rintn1,rintn2,rintn  =',rintn1,rintn2,rintn
  WRITE(*,'(A,3E15.6)') 'rintp1,rintp2,rintp3 =',rintp1,rintp2,rintp3
  WRITE(*,'(A6,2A15)') 'irint','rintp','rintw'
  DO irint=1,rintn
   WRITE(*,'(I6,2E15.6)') irint,rintp(irint),rintw(irint)
  END DO
  WRITE(*,*) 
  
  WRITE(*,'(A,3I5)')    'pintn1,pintn2,pintn  =',pintn1,pintn2,pintn
  WRITE(*,'(A,3E15.6)') 'pintp1,pintp2,pintp3 =',pintp1,pintp2,pintp3
  WRITE(*,'(A6,2A15)') 'ipint','pintp','rintw'
  DO ipint=1,pintn
   WRITE(*,'(I6,2E15.6)') ipint,pintp(ipint),pintw(ipint)
  END DO
  WRITE(*,*) 
  WRITE(*,*) 'tolerance for interpolations in amplitudes: ',tolerance
  WRITE(*,*)
  
 END SUBROUTINE printamp
 
 ! now define functions for products of amplitudes
 ! with scalars 
 ! the product results in corresponding type of amplitude
 ! communication is not required 
 ! 1. case:  nn amplitude 
 ! 1a. scalar = REAL(spreal)
 FUNCTION prodsprealnn(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMPNN) :: prodsprealnn
  TYPE (AMPNN),INTENT(IN)    :: PSI1
  REAL(spreal),INTENT(IN)    :: lambda
  INTEGER ip12,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in prodsprealnn !!'
  END IF
#endif
  prodsprealnn=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalphaNN   ! set the amplitude 
   DO ip12=1,mynp12
    prodsprealnn%amp(ip12,alpha)=lambda*PSI1%amp(ip12,alpha)
   END DO
  END DO
  
 END FUNCTION prodsprealnn
 
 ! 1b. scalar = REAL(dpreal)
 FUNCTION proddprealnn(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMPNN) :: proddprealnn
  TYPE (AMPNN),INTENT(IN)    :: PSI1
  REAL(dpreal),INTENT(IN)    :: lambda
  INTEGER ip12,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddprealnn !!'
  END IF
#endif
  proddprealnn=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalphaNN   ! set the amplitude 
   DO ip12=1,mynp12
    proddprealnn%amp(ip12,alpha)=lambda*PSI1%amp(ip12,alpha)
   END DO
  END DO
  
 END FUNCTION proddprealnn
 
 ! 1c. scalar = COMPLEX(spreal)
 FUNCTION prodspcmplxnn(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMPNN) :: prodspcmplxnn
  TYPE (AMPNN),INTENT(IN)    :: PSI1
  COMPLEX(spreal),INTENT(IN)    :: lambda
  INTEGER ip12,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in prodspcmplxnn !!'
  END IF
#endif
  prodspcmplxnn=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalphaNN   ! set the amplitude 
   DO ip12=1,mynp12
    prodspcmplxnn%amp(ip12,alpha)=lambda*PSI1%amp(ip12,alpha)
   END DO
  END DO
  
 END FUNCTION prodspcmplxnn
 
 ! 1d. scalar = COMPLEX(dpreal)
 FUNCTION proddpcmplxnn(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMPNN) :: proddpcmplxnn
  TYPE (AMPNN),INTENT(IN)    :: PSI1
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  INTEGER ip12,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddpcmplxnn !!'
  END IF
#endif
  proddpcmplxnn=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalphaNN   ! set the amplitude 
   DO ip12=1,mynp12
    proddpcmplxnn%amp(ip12,alpha)=lambda*PSI1%amp(ip12,alpha)
   END DO
  END DO
  
 END FUNCTION proddpcmplxnn
 
 
 ! 2. case:  3n amplitude 
 ! 2a. scalar = REAL(spreal)
 FUNCTION prodspreal3n(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP3N) :: prodspreal3n
  TYPE (AMP3N),INTENT(IN)    :: PSI1
  REAL(spreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip3,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in prodspreal3n !!'
  END IF
#endif
  prodspreal3n=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha3N   ! set the amplitude 
   DO ip3=1,mynp3
    DO ip12=1,mynp12
     prodspreal3n%amp(ip12,ip3,alpha)=lambda*PSI1%amp(ip12,ip3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION prodspreal3n
 
 ! 2b. scalar = REAL(dpreal)
 FUNCTION proddpreal3n(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP3N) :: proddpreal3n
  TYPE (AMP3N),INTENT(IN)    :: PSI1
  REAL(dpreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip3,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddpreal3n !!'
  END IF
#endif
  proddpreal3n=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha3N   ! set the amplitude 
   DO ip3=1,mynp3
    DO ip12=1,mynp12
     proddpreal3n%amp(ip12,ip3,alpha)=lambda*PSI1%amp(ip12,ip3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION proddpreal3n
 
 ! 2c. scalar = COMPLEX(spreal)
 FUNCTION prodspcmplx3n(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP3N) :: prodspcmplx3n
  TYPE (AMP3N),INTENT(IN)    :: PSI1
  COMPLEX(spreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip3,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in prodspcmplx3n !!'
  END IF
#endif
  prodspcmplx3n=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha3N   ! set the amplitude 
   DO ip3=1,mynp3
    DO ip12=1,mynp12
     prodspcmplx3n%amp(ip12,ip3,alpha)=lambda*PSI1%amp(ip12,ip3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION prodspcmplx3n
 
 ! 2d. scalar = COMPLEX(dpreal)
 FUNCTION proddpcmplx3n(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP3N) :: proddpcmplx3n
  TYPE (AMP3N),INTENT(IN)    :: PSI1
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip3,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddpcmplx3n !!'
  END IF
#endif
  proddpcmplx3n=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha3N   ! set the amplitude 
   DO ip3=1,mynp3
    DO ip12=1,mynp12
     proddpcmplx3n%amp(ip12,ip3,alpha)=lambda*PSI1%amp(ip12,ip3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION proddpcmplx3n
 
 ! 3. case:  4n31 amplitude 
 ! 3a. scalar = REAL(spreal)
 FUNCTION prodspreal4n31(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP4N31) :: prodspreal4n31
  TYPE (AMP4N31),INTENT(IN)    :: PSI1
  REAL(spreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip3,iq4,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in prodspreal4n31 !!'
  END IF
#endif
  prodspreal4n31=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha4N31   ! set the amplitude 
   DO iq4=1,mynq4
    DO ip3=1,mynp3
     DO ip12=1,mynp12
      prodspreal4n31%amp(ip12,ip3,iq4,alpha)=lambda*PSI1%amp(ip12,ip3,iq4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION prodspreal4n31
 
 ! 3b. scalar = REAL(dpreal)
 FUNCTION proddpreal4n31(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP4N31) :: proddpreal4n31
  TYPE (AMP4N31),INTENT(IN)    :: PSI1
  REAL(dpreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip3,iq4,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddpreal4n31 !!'
  END IF
#endif
  proddpreal4n31=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha4N31   ! set the amplitude 
   DO iq4=1,mynq4
    DO ip3=1,mynp3
     DO ip12=1,mynp12
      proddpreal4n31%amp(ip12,ip3,iq4,alpha)=lambda*PSI1%amp(ip12,ip3,iq4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION proddpreal4n31
 
 ! 3c. scalar = COMPLEX(spreal)
 FUNCTION prodspcmplx4n31(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP4N31) :: prodspcmplx4n31
  TYPE (AMP4N31),INTENT(IN)    :: PSI1
  COMPLEX(spreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip3,iq4,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in prodspcmplx4n31 !!'
  END IF
#endif
  prodspcmplx4n31=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha4N31   ! set the amplitude 
   DO iq4=1,mynq4
    DO ip3=1,mynp3
     DO ip12=1,mynp12
      prodspcmplx4n31%amp(ip12,ip3,iq4,alpha)=lambda*PSI1%amp(ip12,ip3,iq4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION prodspcmplx4n31
 
 ! 3d. scalar = COMPLEX(dpreal)
 FUNCTION proddpcmplx4n31(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP4N31) :: proddpcmplx4n31
  TYPE (AMP4N31),INTENT(IN)    :: PSI1
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip3,iq4,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddpcmplx4n31 !!'
  END IF
#endif
  proddpcmplx4n31=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha4N31   ! set the amplitude 
   DO iq4=1,mynq4
    DO ip3=1,mynp3
     DO ip12=1,mynp12
      proddpcmplx4n31%amp(ip12,ip3,iq4,alpha)=lambda*PSI1%amp(ip12,ip3,iq4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION proddpcmplx4n31
 
 ! 4. case:  4n22 amplitude 
 ! 4a. scalar = REAL(spreal)
 FUNCTION prodspreal4n22(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP4N22) :: prodspreal4n22
  TYPE (AMP4N22),INTENT(IN)    :: PSI1
  REAL(spreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip34,iq,beta
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in prodspreal4n22 !!'
  END IF
#endif
  prodspreal4n22=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO beta=1,mynbeta4N22   ! set the amplitude 
   DO iq=1,mynq
    DO ip34=1,mynp34
     DO ip12=1,mynp12
      prodspreal4n22%amp(ip12,ip34,iq,beta)=lambda*PSI1%amp(ip12,ip34,iq,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION prodspreal4n22
 
 ! 4b. scalar = REAL(dpreal)
 FUNCTION proddpreal4n22(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP4N22) :: proddpreal4n22
  TYPE (AMP4N22),INTENT(IN)    :: PSI1
  REAL(dpreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip34,iq,beta
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddpreal4n22 !!'
  END IF
#endif
  proddpreal4n22=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO beta=1,mynbeta4N22   ! set the amplitude 
   DO iq=1,mynq
    DO ip34=1,mynp34
     DO ip12=1,mynp12
      proddpreal4n22%amp(ip12,ip34,iq,beta)=lambda*PSI1%amp(ip12,ip34,iq,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION proddpreal4n22
 
 ! 4c. scalar = COMPLEX(spreal)
 FUNCTION prodspcmplx4n22(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP4N22) :: prodspcmplx4n22
  TYPE (AMP4N22),INTENT(IN)    :: PSI1
  COMPLEX(spreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip34,iq,beta
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in prodspcmplx4n22 !!'
  END IF
#endif
  prodspcmplx4n22=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO beta=1,mynbeta4N22   ! set the amplitude 
   DO iq=1,mynq
    DO ip34=1,mynp34
     DO ip12=1,mynp12
      prodspcmplx4n22%amp(ip12,ip34,iq,beta)=lambda*PSI1%amp(ip12,ip34,iq,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION prodspcmplx4n22
 
 ! 4d. scalar = COMPLEX(dpreal)
 FUNCTION proddpcmplx4n22(lambda,PSI1)
  IMPLICIT NONE
  TYPE (AMP4N22) :: proddpcmplx4n22
  TYPE (AMP4N22),INTENT(IN)    :: PSI1
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  INTEGER ip12,ip34,iq,beta
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddpcmplx4n22 !!'
  END IF
#endif
  proddpcmplx4n22=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO beta=1,mynbeta4N22   ! set the amplitude 
   DO iq=1,mynq
    DO ip34=1,mynp34
     DO ip12=1,mynp12
      proddpcmplx4n22%amp(ip12,ip34,iq,beta)=lambda*PSI1%amp(ip12,ip34,iq,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION proddpcmplx4n22
 
 
 ! now define functions for scalar products of amplitudes 
 ! the scalar product results in COMPLEX(spreal)
 ! communication is required, therefore scalarproduct has to be 
 ! call by all PE's in commamp 
 
 ! 1. case:  nn amplitude 
 
 FUNCTION dotprod_nn(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMPNN),INTENT(IN)    :: PSI2
  TYPE (AMPNN),INTENT(IN)    :: PSI1
  COMPLEX(dpreal) :: dotprod_nn
  COMPLEX(dpreal) :: locprod
  INTEGER ip12,alpha,ierr
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in dotprod !!'
  END IF
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in dotprod !!'
  END IF
#endif
  
  ! do local summation with p12 weights
  
  locprod=0.0_dpreal    
  DO alpha=1,mynalphaNN
   DO ip12=1,mynp12
    locprod=locprod &
         + conjg(PSI1%amp(ip12,alpha))*PSI2%amp(ip12,alpha) &
         *p12weight(ip12)
   END DO
  END DO
  
  ! sum local results using MPI 
  CALL MPI_ALLREDUCE(locprod,dotprod_nn,1,MPI_COMPLEX16,MPI_SUM,commampnn,ierr) ! here it is assumed that product is COMPLEX(dpreal)
  
 END FUNCTION dotprod_nn
 
 ! 2. case:  3n amplitude 
 
 FUNCTION dotprod_3n(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP3N),INTENT(IN)    :: PSI2
  TYPE (AMP3N),INTENT(IN)    :: PSI1
  COMPLEX(dpreal) :: dotprod_3n
  COMPLEX(dpreal) :: locprod
  INTEGER ip12,ip3,alpha,ierr
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in dotprod !!'
  END IF
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in dotprod !!'
  END IF
#endif
  
  ! do local summation with p12,p3 weights
  
  locprod=0.0_dpreal
  DO alpha=1,mynalpha3N
   DO ip3=1,mynp3
    DO ip12=1,mynp12
     locprod=locprod &
          + conjg(PSI1%amp(ip12,ip3,alpha))*PSI2%amp(ip12,ip3,alpha) &
          *p12p3weight(ip12,ip3)
    END DO
   END DO
  END DO
  
  ! sum local results using MPI 
  CALL MPI_ALLREDUCE(locprod,dotprod_3n,1,MPI_COMPLEX16,MPI_SUM,commamp3n,ierr) ! here it is assumed that product is COMPLEX(dpreal)
  
 END FUNCTION dotprod_3n
 
 ! 3. case:  4n31 amplitude 
 
 FUNCTION dotprod_4n31(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP4N31),INTENT(IN)    :: PSI2
  TYPE (AMP4N31),INTENT(IN)    :: PSI1
  COMPLEX(dpreal) :: dotprod_4n31
  COMPLEX(dpreal) :: locprod
  INTEGER ip12,ip3,iq4,alpha,ierr
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in dotprod !!'
  END IF
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in dotprod !!'
  END IF
#endif
  
  ! do local summation with p12,p3,q4 weights
  
  locprod=0.0_dpreal
  DO alpha=1,mynalpha4N31
   DO iq4=1,mynq4
    DO ip3=1,mynp3
     DO ip12=1,mynp12
      locprod=locprod &
           + conjg(PSI1%amp(ip12,ip3,iq4,alpha))*PSI2%amp(ip12,ip3,iq4,alpha) &
           *p12p3q4weight(ip12,ip3,iq4)
     END DO
    END DO
   END DO
  END DO
  
  ! sum local results using MPI 
  CALL MPI_ALLREDUCE(locprod,dotprod_4n31,1,MPI_COMPLEX16,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is COMPLEX(dpreal)
  
 END FUNCTION dotprod_4n31
 
 ! 4. case:  4n22 amplitude 
 
 FUNCTION dotprod_4n22(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP4N22),INTENT(IN)    :: PSI2
  TYPE (AMP4N22),INTENT(IN)    :: PSI1
  COMPLEX(dpreal) :: dotprod_4n22
  COMPLEX(dpreal) :: locprod
  INTEGER ip12,ip34,iq,beta,ierr
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in dotprod !!'
  END IF
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in dotprod !!'
  END IF
#endif
  
  ! do local summation with p12,p34,q weights
  
  locprod=0.0_dpreal
  DO beta=1,mynbeta4N22
   DO iq=1,mynq
    DO ip34=1,mynp34
     DO ip12=1,mynp12
      locprod=locprod &
           + conjg(PSI1%amp(ip12,ip34,iq,beta))*PSI2%amp(ip12,ip34,iq,beta) &
           *p12p34qweight(ip12,ip34,iq)
     END DO
    END DO
   END DO
  END DO
  
  ! sum local results using MPI 
  CALL MPI_ALLREDUCE(locprod,dotprod_4n22,1,MPI_COMPLEX16,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is COMPLEX(dpreal)
  
 END FUNCTION dotprod_4n22
 
 
 ! module procedures for assigment statement 
 
 ! first part: assigment of amplitude to amplitude 
 
 SUBROUTINE set_psipsinn(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMPNN),INTENT(IN)    :: PSI2
  TYPE (AMPNN),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynalphann))
  END IF
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  
  PSI1%amp=PSI2%amp   ! copy amplitude element by element 
  
 END SUBROUTINE set_psipsinn
 
 SUBROUTINE set_psipsi3n(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP3N),INTENT(IN)    :: PSI2
  TYPE (AMP3N),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp3,mynalpha3n))
  END IF
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  
  PSI1%amp=PSI2%amp   ! copy amplitude element by element 
  
 END SUBROUTINE set_psipsi3n
 
 SUBROUTINE set_psipsi4n31(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP4N31),INTENT(IN)    :: PSI2
  TYPE (AMP4N31),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp3,mynq4,mynalpha4n31))
  END IF
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  
  PSI1%amp=PSI2%amp   ! copy amplitude element by element 
  
 END SUBROUTINE set_psipsi4n31
 
 SUBROUTINE set_psipsi4n22(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP4N22),INTENT(IN)    :: PSI2
  TYPE (AMP4N22),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp34,mynq,mynbeta4n22))
  END IF
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  
  PSI1%amp=PSI2%amp   ! copy amplitude element by element 
  
 END SUBROUTINE set_psipsi4n22
 
 ! assignment of spcmplx to amplitude 
 
 SUBROUTINE set_psispcmplxnn(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(spreal),INTENT(IN) :: lambda
  TYPE (AMPNN),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynalphann))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psispcmplxnn
 
 SUBROUTINE set_psispcmplx3n(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(spreal),INTENT(IN)    :: lambda
  TYPE (AMP3N),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp3,mynalpha3n))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psispcmplx3n
 
 SUBROUTINE set_psispcmplx4n31(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(spreal),INTENT(IN)    :: lambda
  TYPE (AMP4N31),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp3,mynq4,mynalpha4n31))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psispcmplx4n31
 
 SUBROUTINE set_psispcmplx4n22(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(spreal),INTENT(IN)    :: lambda
  TYPE (AMP4N22),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp34,mynq,mynbeta4n22))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psispcmplx4n22
 
 
 ! assignment of spreal to amplitude 
 
 SUBROUTINE set_psisprealnn(PSI1,lambda)
  IMPLICIT NONE
  REAL(spreal),INTENT(IN) :: lambda
  TYPE (AMPNN),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynalphann))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psisprealnn
 
 SUBROUTINE set_psispreal3n(PSI1,lambda)
  IMPLICIT NONE
  REAL(spreal),INTENT(IN)    :: lambda
  TYPE (AMP3N),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp3,mynalpha3n))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psispreal3n
 
 SUBROUTINE set_psispreal4n31(PSI1,lambda)
  IMPLICIT NONE
  REAL(spreal),INTENT(IN)    :: lambda
  TYPE (AMP4N31),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp3,mynq4,mynalpha4n31))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psispreal4n31
 
 SUBROUTINE set_psispreal4n22(PSI1,lambda)
  IMPLICIT NONE
  REAL(spreal),INTENT(IN)    :: lambda
  TYPE (AMP4N22),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp34,mynq,mynbeta4n22))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psispreal4n22
 
 
 ! assignment of dpcmplx to amplitude 
 
 SUBROUTINE set_psidpcmplxnn(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(dpreal),INTENT(IN) :: lambda
  TYPE (AMPNN),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynalphann))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psidpcmplxnn
 
 SUBROUTINE set_psidpcmplx3n(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  TYPE (AMP3N),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp3,mynalpha3n))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psidpcmplx3n
 
 SUBROUTINE set_psidpcmplx4n31(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  TYPE (AMP4N31),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp3,mynq4,mynalpha4n31))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psidpcmplx4n31
 
 SUBROUTINE set_psidpcmplx4n22(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  TYPE (AMP4N22),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp34,mynq,mynbeta4n22))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psidpcmplx4n22
 
 
 ! assignment of dpreal to amplitude 
 
 SUBROUTINE set_psidprealnn(PSI1,lambda)
  IMPLICIT NONE
  REAL(dpreal),INTENT(IN) :: lambda
  TYPE (AMPNN),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynalphann))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psidprealnn
 
 SUBROUTINE set_psidpreal3n(PSI1,lambda)
  IMPLICIT NONE
  REAL(dpreal),INTENT(IN)    :: lambda
  TYPE (AMP3N),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp3,mynalpha3n))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psidpreal3n
 
 SUBROUTINE set_psidpreal4n31(PSI1,lambda)
  IMPLICIT NONE
  REAL(dpreal),INTENT(IN)    :: lambda
  TYPE (AMP4N31),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp3,mynq4,mynalpha4n31))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psidpreal4n31
 
 SUBROUTINE set_psidpreal4n22(PSI1,lambda)
  IMPLICIT NONE
  REAL(dpreal),INTENT(IN)    :: lambda
  TYPE (AMP4N22),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynp12,mynp34,mynq,mynbeta4n22))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE set_psidpreal4n22
 
 ! now the set of procedure for the OPERATOR(+)
 ! 1. case nn amplitude 
 
 FUNCTION plus_psinn(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMPNN),INTENT(IN)    :: PSI2
  TYPE (AMPNN),INTENT(IN)    :: PSI1
  TYPE (AMPNN)  :: plus_psinn
  INTEGER alpha,ip12
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in plus_psinn !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in plus_psinn !!'
  END IF
#endif
  
  plus_psinn=0.0_spreal   ! allocate & predefine psinn 
  
  DO alpha=1,mynalphaNN
   DO ip12=1,mynp12
    plus_psinn%amp(ip12,alpha)=PSI1%amp(ip12,alpha)+PSI2%amp(ip12,alpha)
   END DO
  END DO
  
 END FUNCTION plus_psinn
 
 ! 2. case 3n amplitude 
 
 FUNCTION plus_psi3n(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP3N),INTENT(IN)    :: PSI2
  TYPE (AMP3N),INTENT(IN)    :: PSI1
  TYPE (AMP3N)   :: plus_psi3n
  INTEGER alpha,ip12,ip3
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in plus_psi3n !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in plus_psi3n !!'
  END IF
#endif
  
  plus_psi3n=0.0_spreal   ! allocate & predefine psi3n 
  
  DO alpha=1,mynalpha3N
   DO ip3=1,mynp3
    DO ip12=1,mynp12
     plus_psi3n%amp(ip12,ip3,alpha)=PSI1%amp(ip12,ip3,alpha)+PSI2%amp(ip12,ip3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION plus_psi3n
 
 ! 3. case 4n31 amplitude 
 
 FUNCTION plus_psi4n31(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP4N31),INTENT(IN)    :: PSI2
  TYPE (AMP4N31),INTENT(IN)    :: PSI1
  TYPE (AMP4N31)   :: plus_psi4n31
  INTEGER alpha,ip12,ip3,iq4
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in plus_psi4n31 !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in plus_psi4n31 !!'
  END IF
#endif
  
  plus_psi4n31=0.0_spreal   ! allocate & predefine psi4n31 
  
  DO alpha=1,mynalpha4N31
   DO iq4=1,mynq4
    DO ip3=1,mynp3
     DO ip12=1,mynp12
      plus_psi4n31%amp(ip12,ip3,iq4,alpha)=PSI1%amp(ip12,ip3,iq4,alpha)&
           +PSI2%amp(ip12,ip3,iq4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION plus_psi4n31
 
 ! 4. case 4n22 amplitude 
 
 FUNCTION plus_psi4n22(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP4N22),INTENT(IN)    :: PSI2
  TYPE (AMP4N22),INTENT(IN)    :: PSI1
  TYPE (AMP4N22)   :: plus_psi4n22
  INTEGER beta,ip12,ip34,iq
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in plus_psi4n22 !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in plus_psi4n22 !!'
  END IF
#endif
  
  plus_psi4n22=0.0_spreal   ! allocate & predefine psi4n22 
  
  DO beta=1,mynbeta4N22
   DO iq=1,mynq
    DO ip34=1,mynp34
     DO ip12=1,mynp12
      plus_psi4n22%amp(ip12,ip34,iq,beta)=PSI1%amp(ip12,ip34,iq,beta)&
           +PSI2%amp(ip12,ip34,iq,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION plus_psi4n22
 
 ! now the set of procedure for the OPERATOR(-)
 ! 1. case nn amplitude 
 
 FUNCTION minus_psinn(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMPNN),INTENT(IN)    :: PSI2
  TYPE (AMPNN),INTENT(IN)    :: PSI1
  TYPE (AMPNN)  :: minus_psinn
  INTEGER alpha,ip12
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in minus_psinn !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in minus_psinn !!'
  END IF
#endif
  
  minus_psinn=0.0_spreal   ! allocate & predefine psinn 
  
  DO alpha=1,mynalphaNN
   DO ip12=1,mynp12
    minus_psinn%amp(ip12,alpha)=PSI1%amp(ip12,alpha)-PSI2%amp(ip12,alpha)
   END DO
  END DO
  
 END FUNCTION minus_psinn
 
 ! 2. case 3n amplitude 
 
 FUNCTION minus_psi3n(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP3N),INTENT(IN)    :: PSI2
  TYPE (AMP3N),INTENT(IN)    :: PSI1
  TYPE (AMP3N)   :: minus_psi3n
  INTEGER alpha,ip12,ip3
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in minus_psi3n !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in minus_psi3n !!'
  END IF
#endif
  
  minus_psi3n=0.0_spreal   ! allocate & predefine psi3n 
  
  DO alpha=1,mynalpha3N
   DO ip3=1,mynp3
    DO ip12=1,mynp12
     minus_psi3n%amp(ip12,ip3,alpha)=PSI1%amp(ip12,ip3,alpha)-PSI2%amp(ip12,ip3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION minus_psi3n
 
 ! 3. case 4n31 amplitude 
 
 FUNCTION minus_psi4n31(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP4N31),INTENT(IN)    :: PSI2
  TYPE (AMP4N31),INTENT(IN)    :: PSI1
  TYPE (AMP4N31)   :: minus_psi4n31
  INTEGER alpha,ip12,ip3,iq4
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in minus_psi4n31 !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in minus_psi4n31 !!'
  END IF
#endif
  
  minus_psi4n31=0.0_spreal   ! allocate & predefine psi4n31 
  
  DO alpha=1,mynalpha4N31
   DO iq4=1,mynq4
    DO ip3=1,mynp3
     DO ip12=1,mynp12
      minus_psi4n31%amp(ip12,ip3,iq4,alpha)=PSI1%amp(ip12,ip3,iq4,alpha)&
           -PSI2%amp(ip12,ip3,iq4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION minus_psi4n31
 
 ! 4. case 4n22 amplitude 
 
 FUNCTION minus_psi4n22(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP4N22),INTENT(IN)    :: PSI2
  TYPE (AMP4N22),INTENT(IN)    :: PSI1
  TYPE (AMP4N22)   :: minus_psi4n22
  INTEGER beta,ip12,ip34,iq
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in minus_psi4n22 !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in minus_psi4n22 !!'
  END IF
#endif
  
  minus_psi4n22=0.0_spreal   ! allocate & predefine psi4n22 
  
  DO beta=1,mynbeta4N22
   DO iq=1,mynq
    DO ip34=1,mynp34
     DO ip12=1,mynp12
      minus_psi4n22%amp(ip12,ip34,iq,beta)=PSI1%amp(ip12,ip34,iq,beta)&
           -PSI2%amp(ip12,ip34,iq,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION minus_psi4n22
 
 
 !! collect channel subroutines for easier access to full amplitudes on one processors 
 !! mostly required for I/O of amplitudes 
 
 SUBROUTINE collect_channel_NN(phi,ampch,alpha)
  IMPLICIT NONE 
  TYPE (AMPNN)   :: phi
  COMPLEX(spreal) :: ampch(p12n)
  COMPLEX(spreal),ALLOCATABLE :: step1(:)
  INTEGER pe,alpha,localpha,ierr
  
#ifdef DEBUG
   IF(alpha.LE.0.OR.alpha.GT.alphaNNcdepmax) THEN
    WRITE(*,*) 'alpha in collect_channel_NN',alpha
    STOP 'problem in collect'
   END IF
#endif 
  
  pe=mod(alpha-1,npe_alpha)
  localpha=(alpha-1-pe)/npe_alpha+1
  
  ALLOCATE(step1(mynp12))
  
  IF(myalphaid.EQ.pe) THEN   
   step1=phi%amp(:,localpha)
  END IF
  
  CALL MPI_BCAST(step1,mynp12,MPI_COMPLEX8, &
 &   pe,commalpha,ierr)
  CALL collect_part_spcmplx(step1,ampch,1,p12n,1,commp12)
  
  DEALLOCATE(step1)
    
 END SUBROUTINE collect_channel_NN
 
 !! collect channel subroutines for easier access to full amplitudes on one processors 
 !! mostly required for I/O of amplitudes 
 
 SUBROUTINE collect_channel_3n(phi,ampch,alpha)
  IMPLICIT NONE 
  TYPE (AMP3N)   :: phi
  COMPLEX(spreal) :: ampch(p12n,p3n)
  COMPLEX(spreal),ALLOCATABLE :: step1(:,:),step2(:,:)
  INTEGER pe,alpha,localpha,ierr
  
#ifdef DEBUG
  IF(alpha.LE.0.OR.alpha.GT.alpha3Ncdepmax) THEN
   WRITE(*,*) 'alpha in collect_channel_3N',alpha
   STOP 'problem in collect'
  END IF
#endif 
  
  pe=mod(alpha-1,npe_alpha)
  localpha=(alpha-1-pe)/npe_alpha+1
  
  ALLOCATE(step1(mynp12,mynp3))
  ALLOCATE(step2(mynp12,p3n))
  
  IF(myalphaid.EQ.pe) THEN   
   step1=phi%amp(:,:,localpha)
  END IF
  
  CALL MPI_BCAST(step1,mynp12*mynp3,MPI_COMPLEX8,pe,commalpha,ierr)
  CALL collect_part_spcmplx(step1,step2,mynp12,p3n,1,commp3)
  CALL collect_part_spcmplx(step2,ampch,1,p12n,p3n,commp12)
  
  DEALLOCATE(step1,step2)
  
 END SUBROUTINE collect_channel_3n
 
 !!    subroutines to normalize amplitudes 
 SUBROUTINE norm_nn(phi)
  IMPLICIT NONE
  TYPE(ampnn) phi
  TYPE(ampnn) phih
  REAL(dpreal) fakt
  fakt=phi*phi
  fakt=1.0_dpreal/sqrt(fakt)
  phih=fakt*phi
  phi=phih
  DEALLOCATE(phih%amp)
 END SUBROUTINE norm_nn
 
 SUBROUTINE norm_3n(phi)
  IMPLICIT NONE
  TYPE(amp3n) phi
  TYPE(amp3n) phih
  REAL(dpreal) fakt
  fakt=phi*phi
  fakt=1.0_dpreal/sqrt(fakt)
  phih=fakt*phi
  phi=phih
  DEALLOCATE(phih%amp)
 END SUBROUTINE norm_3n
 
 SUBROUTINE norm_4n31(phi)
  IMPLICIT NONE
  TYPE(amp4n31) phi
  TYPE(amp4n31) phih
  REAL(dpreal) fakt
  fakt=phi*phi
  fakt=1.0_dpreal/sqrt(fakt)
  phih=fakt*phi
  phi=phih
  DEALLOCATE(phih%amp)
 END SUBROUTINE norm_4n31
 
 SUBROUTINE norm_4n22(phi)
  IMPLICIT NONE
  TYPE(amp4n22) phi
  TYPE(amp4n22) phih
  REAL(dpreal) fakt
  fakt=phi*phi
  fakt=1.0_dpreal/sqrt(fakt)
  phih=fakt*phi
  phi=phih
  DEALLOCATE(phih%amp)
 END SUBROUTINE norm_4n22
 
 SUBROUTINE norm_4nboth(phi1,phi2)
  IMPLICIT NONE
  TYPE(amp4n31) phi1  
  TYPE(amp4n22) phi2
  TYPE(amp4n31) phih1  
  TYPE(amp4n22) phih2
  REAL(dpreal) fakt
  fakt=phi1*phi1+phi2*phi2
  fakt=1.0_dpreal/sqrt(fakt)
  phih1=fakt*phi1
  phih2=fakt*phi2
  phi1=phih1
  phi2=phih2  
  DEALLOCATE(phih1%amp)
  DEALLOCATE(phih2%amp)
 END SUBROUTINE norm_4nboth

 ! This subroutine creates and prepares a file 'filename' where the
 ! amplitude is written to. It also writes the header of the file which 
 ! contains information about grid points and alpha-bookkeeping.
 ! The 'tag' is written to the file as an identifier at the top
 ! of the file. The 'filenum' is used as an internal reference
 ! and is needed when an amplitude is actually written to the file

 subroutine open_writefile_thatphi_3N(filenum,filename,tag)
  implicit none
  
  integer :: filenum,ip12,ip3,iq0
  integer :: al3n,l12,s12,j12,t12,l3,i3,j3,tau3,mtau3,alphann,pari
  character(len=*) :: filename,tag
  
  ! MASTER ONLY: write header of file

  if (master) then

   open(unit=filenum,file=filename,form='formatted',status='unknown')
  
   write(filenum,*) tag
  
   write(filenum,*) 'number of p12 grid points:'
   write(filenum,*) p12n
   write(filenum,*) 'all p12 grid points:'
   do ip12=1,p12n
    write(filenum,*) p12p(ip12)
   end do
  
   write(filenum,*) 'number of p3 grid points:'
   write(filenum,*) p3n
   write(filenum,*) 'all p3 grid points:'
   do ip3=1,p3n
    write(filenum,*) p3p(ip3)
   end do
   
   write(filenum,*) 'number of alpha3n channels:'
   write(filenum,*) alpha3ncdepmax
   write(filenum,*) 'quantum numbers of all alpha3n channels:'
   do al3n=1,alpha3ncdepmax
    call get3nqn(al3n,l12,s12,j12,t12,l3,i3,j3,tau3,mtau3,alphann,pari)
    write(filenum,*) l12,s12,j12,t12,l3,i3,j3,tau3,mtau3,alphann,pari
   end do
  
   write(filenum,*) 'number of energies:'
   write(filenum,*) q0n
   write(filenum,*) 'all energies:'
   do iq0=1,q0n
    write(filenum,*) ener(iq0)
   end do
   
  end if ! master

 end subroutine open_writefile_thatphi_3N
 

 ! This subroutine writes a specific amplitude to the file with the
 ! internal reference number 'filenum' which was set in the routine
 ! 'open_writefile_thatphi_3N'. The amplitude 'amp_in' and the parameters
 ! 'ener', 'md' and 'mns' are input parameters which are written
 ! to the file.
 
 subroutine writefile_thatphi_3N(filenum,amp_in,ener,md,mns)
  implicit none
  
  integer :: filenum,md,mns,al3n,local3n,ierr
  complex(dpreal) :: ener
  complex(dpreal),allocatable :: amp_temp1(:,:,:),amp_temp2(:,:,:),amp_temp_write(:,:)
  type(amp3n) :: amp_in
  
  if (master) then
   
   write(filenum,*) 'Amplitude for energy (fm^-1) ='
   write(filenum,*) ener
   write(filenum,*) 'md ='
   write(filenum,*) md
   write(filenum,*) 'mns ='
   write(filenum,*) mns
   
  end if

  allocate(amp_temp1(mynp12,mynp3,mynalpha3n))
  amp_temp1=amp_in%amp

  allocate(amp_temp2(p12n,mynp3,mynalpha3n))
  call collect_part_dpcmplx(amp_temp1,amp_temp2,1,p12n,mynp3*mynalpha3n,commp12)
  deallocate(amp_temp1)
  
  allocate(amp_temp1(p12n,p3n,mynalpha3n))
  call collect_part_dpcmplx(amp_temp2,amp_temp1,p12n,p3n,mynalpha3n,commp3)
  deallocate(amp_temp2)

  allocate(amp_temp_write(p12n,p3n))
  
  ! ALL: broadcast all local alpha3n channels to all processes
  ! MASTER: write the amplitude to the file channel by channel
  
  do al3n=1,alpha3ncdepmax ! al3n is the global alpha3n index
   if (mod(al3n-1,npe_alpha).eq.myalphaid) then
    local3n=(al3n-myalpha3n)/npe_alpha+1
    amp_temp_write=amp_temp1(:,:,local3n)
   end if
   call mpi_bcast(amp_temp_write,p12n*p3n,mpi_complex16,mod(al3n-1,npe_alpha),commalpha,ierr)
   if (master) then
    write(filenum,'(6E15.6)') amp_temp_write

!!$    do ip3=1,p3n
!!$     do ip12=1,p12n
!!$      write(filenum,*) amp_temp_write(ip12,ip3)
!!$     end do
!!$    end do
      
   end if
  end do
  
  deallocate(amp_temp1)
  deallocate(amp_temp_write)
  
 end subroutine writefile_thatphi_3N
 

 ! This subroutine reads the header of the file 'filename' where the
 ! amplitude is read from. It also prepares interpolation and transition
 ! from the bookkeeping in the file to the current bookkeeping. Again,
 ! 'filenum' is an internal reference number needed again in the subsequent
 ! application of 'readfile_thatphi_3N'.

 subroutine open_readfile_thatphi_3N(filenum,filename)
  implicit none
  
  integer :: filenum,ip12,ip3,q0n_read,iq0,ierr
  integer :: al3n,l12,s12,j12,t12,l3,i3,j3,tau3,mtau3,alphann,pari
  integer :: al3n_,l12_,s12_,j12_,t12_,l3_,i3_,j3_,tau3_,mtau3_,alphann_,pari_
  real(dpreal),allocatable :: p12p_read(:),p3p_read(:)
  character(len=*) :: filename
  
  ! read parameters of scattering amplitude by master process
  ! and broadcast to all processes
  
  if (master) then
   open(unit=filenum,file=filename,form='formatted',status='unknown')

   read(filenum,*) ! read over tag in the first line

   ! read p12n and all p12 grid points from the file
   read(filenum,*) ! read over comment
   read(filenum,*) p12n_read
  end if

  call mpi_bcast(p12n_read,1,mpi_integer,0,commall,ierr)
  if(.not.allocated(p12p_read)) allocate(p12p_read(p12n_read))

  if (master) then
   read(filenum,*) ! read over comment
   do ip12=1,p12n_read
    read(filenum,*) p12p_read(ip12)
   end do
  end if

  call mpi_bcast(p12p_read,p12n_read,mpi_real8,0,commall,ierr)

  if (master) then
   ! read p3n and all p3 grid points from the file
   read(filenum,*) ! read over comment
   read(filenum,*) p3n_read
  end if

  call mpi_bcast(p3n_read,1,mpi_integer,0,commall,ierr)
  if(.not.allocated(p3p_read)) allocate(p3p_read(p3n_read))
   
  if (master) then
   read(filenum,*) ! read over comment
   do ip3=1,p3n_read
    read(filenum,*) p3p_read(ip3)
   end do
  end if

  call mpi_bcast(p3p_read,p3n_read,mpi_real8,0,commall,ierr)

  if (master) then   
   ! read alpha3n channels from the file
   read(filenum,*) ! read over comment
   read(filenum,*) alpha3ncdepmax_read
   read(filenum,*) ! read over comment
  end if
   
  call mpi_bcast(alpha3ncdepmax_read,1,mpi_integer,0,commall,ierr)
  if(.not.allocated(alpha3n_index)) allocate(alpha3n_index(alpha3ncdepmax_read))
   
  if (master) then
   do al3n_=1,alpha3ncdepmax_read
    
    read(filenum,*) l12_,s12_,j12_,t12_,l3_,i3_,j3_,tau3_,mtau3_,alphann_,pari_
    
    ! create index file for alpha3n indices
    
    alpha3n_index(al3n_)=0
    do al3n=1,alpha3ncdepmax
     call get3nqn(al3n,l12,s12,j12,t12,l3,i3,j3,tau3,mtau3,alphann,pari)
     if ((l12.eq.l12_).and.(s12.eq.s12_).and.(j12.eq.j12_).and.(t12.eq.t12_)&
          .and.(l3.eq.l3_).and.(i3.eq.i3_).and.(j3.eq.j3_).and.(tau3.eq.tau3_)&
          .and.(mtau3.eq.mtau3_).and.(alphann.eq.alphann_).and.(pari.eq.pari_)) then
      alpha3n_index(al3n_)=al3n
     end if
    end do
    
   end do ! al3n_
  end if ! master
  
  call mpi_bcast(alpha3n_index,alpha3ncdepmax_read,mpi_integer,0,commall,ierr)

  if (master) then
   read(filenum,*) ! read over comment
   read(filenum,*) q0n_read
   read(filenum,*) ! read over comment
   do iq0=1,q0n_read
    read(filenum,*) ! read over energies
   end do  
  end if ! master

  call mpi_bcast(q0n_read,1,mpi_integer,0,commall,ierr)

  ! ALL: prepare interpolation to local p12 and p3

  if(.not.allocated(indx_p12_read_p12)) allocate(indx_p12_read_p12(p12n,4))
  if(.not.allocated(spl_p12_read_p12)) allocate(spl_p12_read_p12(p12n,4))
  call cubherm_dp(p12p_read,p12n_read,p12p,p12n,spl_p12_read_p12,indx_p12_read_p12)
  
  if(.not.allocated(indx_p3_read_p3)) allocate(indx_p3_read_p3(p3n,4))
  if(.not.allocated(spl_p3_read_p3)) allocate(spl_p3_read_p3(p3n,4))
  call cubherm_dp(p3p_read,p3n_read,p3p,p3n,spl_p3_read_p3,indx_p3_read_p3)

 end subroutine open_readfile_thatphi_3N
 

 ! This subroutine reads the amplitudes stored in the file with
 ! internal reference number 'filenum' subsequently. The amplitude
 ! is given in 'amp_out' and is stored using current grid points and
 ! alpha-bookkeeping according to the results obtained in the
 ! subroutine 'open_readfile_thatphi_3N'. The parameters 'ener',
 ! 'md', and 'mns' are read from the file and given as output from
 ! the subroutine
 
 subroutine readfile_thatphi_3N(filenum,amp_out,ener,md,mns)
  implicit none
  
  integer :: filenum,md,mns,indx_p12,indx_p3,ip12,ip3,al3n_,al3n,r,s,local3n,ierr
  complex(dpreal) :: ener
  complex(dpreal) :: spl_p12,spl_p3
  complex(dpreal),allocatable :: amp_temp_read(:,:)
  type(amp3n) :: amp_out
  
  amp_out=0.0_dpreal
  
  ! MASTER ONLY: read amplitude header

  if (master) then
   
   read(filenum,*) ! read over comment
   read(filenum,*) ener
   read(filenum,*) ! read over comment
   read(filenum,*) md
   read(filenum,*) ! read over comment
   read(filenum,*) mns

  end if

  ! ALL: broadcast header parameters

  call mpi_bcast(ener,1,mpi_real8,0,commall,ierr)
  call mpi_bcast(md,1,mpi_integer,0,commall,ierr)
  call mpi_bcast(mns,1,mpi_integer,0,commall,ierr)
  
  allocate(amp_temp_read(p12n_read,p3n_read))

  ! ALL: go through all amplitudes available in file

  do al3n_=1,alpha3ncdepmax_read
    
   if (master) then
    read(filenum,'(6E15.6)') amp_temp_read
   end if

   al3n=alpha3n_index(al3n_) ! al3n is the global alpha3n index
!   write(*,*) 'al3n =',al3n
    
   ! send out amplitudes to all processes with alpha3n channels present in current bookkeeping
    
   if (al3n.ne.0) then
    call mpi_bcast(amp_temp_read,p12n_read*p3n_read,mpi_complex16,0,commall,ierr)
    if (mod(al3n-1,npe_alpha).eq.myalphaid) then

     do ip12=1,mynp12
      do ip3=1,mynp3
       
       do r=1,4
        do s=1,4
         indx_p12=indx_p12_read_p12(myp12+ip12-1,r)
         spl_p12=spl_p12_read_p12(myp12+ip12-1,r)
         indx_p3=indx_p3_read_p3(myp3+ip3-1,s)
         spl_p3=spl_p3_read_p3(myp3+ip3-1,s)

         local3n=(al3n-myalpha3n)/npe_alpha+1
         amp_out%amp(ip12,ip3,local3n)=amp_out%amp(ip12,ip3,local3n)+spl_p12*spl_p3*amp_temp_read(indx_p12,indx_p3)
      
        end do ! s
       end do ! r

      end do ! ip3
     end do ! ip12
     
    end if ! myalphaid
    
   end if ! amplitude present condition

  end do ! al3n_

  deallocate(amp_temp_read)
   
 end subroutine readfile_thatphi_3N

 ! now a set of subroutines to read and write an amplitude
 ! it is assumed that HDF5 is started 
 
 
 ! NN amplitude
 ! open file <filename>
 ! write amplitude PSI to file <filename>
 ! in the group ampname
 ! file is then closed 
 
 SUBROUTINE writeamp_nn(PSI,filename,ampname)
  IMPLICIT NONE 
  TYPE(AMPNN),TARGET :: PSI 
  CHARACTER(LEN=*) filename
  CHARACTER(LEN=*) ampname
  INTEGER(HID_T) writef_id
  INTEGER(HID_T) group_amp_id
  INTEGER(HID_T) dset_nn_id,dsp_nn_id,dsp_hyper_id
  INTEGER(HID_T) msp_nn_id,msp_hyper_id
  INTEGER(HSIZE_T) dims(4)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
  INTEGER(HSIZE_T) start(4)     ! storage for offset in output arrays
  INTEGER(HSIZE_T) block(4)     ! storage for blocksize in output arrays
  INTEGER(HSIZE_T) stride(4)    ! storage for stride in output arrays
  INTEGER(HSIZE_T) count(4)     ! storage for count in output arrays
  INTEGER ierr
  REAL(spreal),ALLOCATABLE  :: buf(:,:,:)
  INTEGER(lint) blocksize

  writeamptime=writeamptime-MPI_WTIME()
  
  ! first open the hdf file
  ! is always assumed that the file exits
  ! copy empty file if necessary
  
  CALL open_write_h5file(filename,commall,writef_id)
  
  ! then create the group <ampname>
  
  CALL h5gcreate_f(writef_id,trim(ampname),group_amp_id, ierr)
  
  ! file should contain all info to interprete the data 
  !    1. meshpoints for P12 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  CALL writep12p(group_amp_id)
  
  !    2. NN channels qnalphaNN(:,:) including number of channels 
  CALL writennchannels(group_amp_id)
  
  !    3. physical constants from constant mod are also written 
  CALL writephysconst(group_amp_id) 
  !    4. NN amplitude 
  hdftime=hdftime-MPI_WTIME()
  ! this needs to be written in parallel since the data 
  ! is distributed 
  ! make sure that that only one processor writes the data 
  ! in case the distribution is not complete 
  
  ! here data is distributed in in commp12 and commalpha 
  ! that writing processors should have myp3id=myq4id=myenerid=0
  ! put the hyperslabs for other processors to empty 
  
  !     need to create a data set in the file 
  !     with complete number amp-elements 
  !     outline shape of global set of data 
  dims(1)=2                  !  complex = 2 * real
  dims(2)=P12N          !  write needs dimensions of arrays 
  dims(3)=alphaNNcdepmax     !  here array size in data file  
  !     create corresponding data space            
  CALL h5screate_simple_f(3,dims,dsp_nn_id,ierr);

  !     describe the memory layout of the data
  dims(1)=2                  ! complex = 2 * real
  dims(2)=mynp12        ! write needs dimension of arrays
  dims(3)=mynalphaNN         ! here size of array in memory 
  !     create corresponding mem space    
  CALL h5screate_simple_f(3,dims,msp_nn_id,ierr);
  
  ! select the hyperslap in the file  that corresponds to the local 
  ! NN data 
  IF(myp3id.EQ.0 .AND. myq4id.EQ.0 .AND. myenerid.EQ.0) THEN
   start(1)=0
   start(2)=myp12-1             
   start(3)=myalphaid 
   block(1)=2
   block(2)=mynp12
   block(3)=1
   stride(1)=1
   stride(2)=1
   stride(3)=npe_alpha
   count(1)=1
   count(2)=1
   count(3)=mynalphaNN
  ELSE
   start=0
   block=1
   stride=1
   count=0
  END IF
  
  CALL h5scopy_f(dsp_nn_id,dsp_hyper_id,ierr)
  CALL h5sselect_hyperslab_f(dsp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)

  ! select the hyperslap in memorythat corresponds to the local 
  ! NN data 
  IF(myp3id.EQ.0 .AND. myq4id.EQ.0 .AND. myenerid.EQ.0) THEN
   start(1)=0
   start(2)=0             
   start(3)=0
   block(1)=2
   block(2)=mynp12
   block(3)=mynalphaNN
   stride(1)=1
   stride(2)=1
   stride(3)=1
   count(1)=1
   count(2)=1
   count(3)=1
  ELSE
   start=0
   block=1
   stride=1
   count=0
  END IF   
     
  CALL h5scopy_f(msp_nn_id,msp_hyper_id,ierr)
  CALL h5sselect_hyperslab_f(msp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)
         
!    create the dataset in file 
  CALL h5dcreate_f(group_amp_id,trim('nnamp'),H5T_NATIVE_REAL,&
     &             dsp_nn_id,dset_nn_id,ierr)
    
     
! resolve compatiblity issue by copying to REAL buffer 
  ALLOCATE(buf(2,mynp12,mynalphaNN))
  buf(1,:,:)=REAL(PSI%amp(:,:))
  buf(2,:,:)=AIMAG(PSI%amp(:,:))

#ifdef HDF_MPIO
      ! check for consistency with MPI
  blocksize=4_lint*2*mynp12*mynalphaNN  
  IF(blocksize.GE.blocklimit) THEN
    WRITE(*,*) 'block too large to be handled with MPI-IO'
    WRITE(*,*) 'blocksize = ',blocksize
    WRITE(*,*) 'blocklimit = ',blocklimit
    CALL abort
  END IF
#endif
 
! write data of amplitude (use here collective communication)       
  CALL h5dwrite_f(dset_nn_id, H5T_NATIVE_REAL,& 
     &      buf,dims,ierr,&
     &      mem_space_id = msp_hyper_id,&
     &      file_space_id = dsp_hyper_id,& 
     &      xfer_prp = pcollectwrite_id)

  DEALLOCATE(buf)

  hdftime=hdftime+MPI_WTIME()
  
! close data set   
  CALL h5dclose_f(dset_nn_id,ierr)
  
! do not need mem_space anymore      
  CALL h5sclose_f(msp_hyper_id,ierr)
  CALL h5sclose_f(msp_nn_id,ierr)
! and not data space anymore    
  CALL h5sclose_f(dsp_hyper_id,ierr)
  CALL h5sclose_f(dsp_nn_id,ierr) 
! close the group for this amplitude 
  CALL h5gclose_f(group_amp_id,ierr) 
! finally close the file 
  CALL close_h5file(writef_id)

  writeamptime=writeamptime+MPI_WTIME()
  
 END SUBROUTINE writeamp_nn

 ! NN amplitude
 ! open file <filename>
 ! read amplitude PSI from group ampname 
 ! of file <filename>
 ! file is then closed
 ! print the parameters and give back the amplitude
 
 SUBROUTINE readamp_nn(PSI,filename,ampname)
  IMPLICIT NONE 
  TYPE(AMPNN),TARGET :: PSI 
  CHARACTER(LEN=*) filename
  CHARACTER(LEN=*) ampname
  INTEGER(HID_T) readf_id
  INTEGER(HID_T) group_amp_id
  INTEGER(HID_T) dset_nn_id,dsp_nn_id,dsp_hyper_id
  INTEGER(HID_T) msp_nn_id,msp_hyper_id
  INTEGER(HSIZE_T) dims(4)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
  INTEGER(HSIZE_T) start(4)     ! storage for offset in output arrays
  INTEGER(HSIZE_T) block(4)     ! storage for blocksize in output arrays
  INTEGER(HSIZE_T) stride(4)    ! storage for stride in output arrays
  INTEGER(HSIZE_T) count(4)     ! storage for count in output arrays
  INTEGER ierr,p12n_read,alphacdepmax_read,alphamax_read,ip12,alpha
  ! buffer to read the complex data 
  REAL(spreal),ALLOCATABLE  :: buf(:,:,:),&
     &                         psi_tmp_p12(:,:,:), &
     &                         psi_inter_p12(:,:,:),psi_tmp_alpha(:,:,:)
  REAL(dpreal),ALLOCATABLE  :: p12p_read(:),p12w_read(:)
  INTEGER,ALLOCATABLE  :: qn_read(:,:)
  INTEGER mynp12_read,myp12_read,myendp12_read,mynalpha_read
  REAL(dpreal),ALLOCATABLE :: spl12(:,:)
  INTEGER,ALLOCATABLE :: indx12(:,:),alphap(:)
  REAL(dpreal) :: fakt 
  INTEGER alphaloc,l12,s12,j12,t12,mt12
  INTEGER(lint) blocksize
  
  readamptime=readamptime-MPI_WTIME()
  ! first open the hdf file
  ! is always assumed that the file exits
  ! copy empty file if necessary
  
  CALL open_read_h5file(filename,commall,readf_id)
  
  ! then create the group <ampname>
  
  CALL h5gopen_f(readf_id,trim(ampname),group_amp_id, ierr)
  
  ! file should contain all info to interprete the data 
  !    1. meshpoints for P12 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  !       the meshpoints and weights are allocated within the routine
  CALL readp12p(group_amp_id,p12p_read,p12w_read,p12n_read)
  
  !    2. NN channels qn_read(:,:) including number of channels
  !       the quantum numbers are allocated within the routine 
  CALL readnnchannels(group_amp_id,qn_read,alphamax_read,alphacdepmax_read)
  
  !    3. physical constants from constant mod are also read 
  !       and printed. No values are used. 
  CALL readphysconst(group_amp_id) 
  
  !    4. NN amplitude 
  
  ! first read the amplitude using the mesh stored 
  
  ! distribution of mesh and channels 
  
  CALL distr_block(p12n_read,npe_p12,myp12id,myp12_read,myendp12_read,mynp12_read)
  CALL distr_piece(alphacdepmax_read,npe_alpha,myalphaid,mynalpha_read)
  
  ! prepare buffer for real and imaginary part of the amplitude 
  
  hdftime=hdftime-MPI_WTIME()
  
  ALLOCATE(buf(2,mynp12_read,mynalpha_read))

  
  ! here data is distributed in in commp12 and commalpha 
  ! the reading processors should have myp3id=myq4id=myenerid=0
  ! this is element 0 in commp3q4ener
  ! put the hyperslabs for other processors to empty 
  
  !     need to open the data set in the file 
  CALL h5dopen_f(group_amp_id,trim('nnamp'),dset_nn_id,ierr)
  
  ! get the data space descriptor 
  CALL h5dget_space_f(dset_nn_id, dsp_nn_id, ierr)
     
  !     describe the memory layout of the data
  dims(1)=2                  ! complex = 2 * real
  dims(2)=mynp12_read        ! read needs dimension of arrays
  dims(3)=mynalpha_read      ! here size of array in memory 
  !     create corresponding mem space    
  CALL h5screate_simple_f(3,dims,msp_nn_id,ierr);
  
  ! select the hyperslap in the file  that corresponds to the local 
  ! NN data 
  IF(myp3q4enerid.EQ.0) THEN
   start(1)=0
   start(2)=myp12_read-1             
   start(3)=myalphaid 
   block(1)=2
   block(2)=mynp12_read
   block(3)=1
   stride(1)=1
   stride(2)=1
   stride(3)=npe_alpha
   count(1)=1
   count(2)=1
   count(3)=mynalpha_read
  ELSE
   start=0
   block=1
   stride=1
   count=0
  END IF
  
  CALL h5scopy_f(dsp_nn_id,dsp_hyper_id,ierr)
  CALL h5sselect_hyperslab_f(dsp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)

  ! select the hyperslap in memorythat corresponds to the local 
  ! NN data 
  IF(myp3q4enerid.EQ.0) THEN
   start(1)=0
   start(2)=0             
   start(3)=0
   block(1)=2
   block(2)=mynp12_read
   block(3)=mynalpha_read
   stride(1)=1
   stride(2)=1
   stride(3)=1
   count(1)=1
   count(2)=1
   count(3)=1
  ELSE
   start=0
   block=1
   stride=1
   count=0
  END IF   
     
  CALL h5scopy_f(msp_nn_id,msp_hyper_id,ierr)
  CALL h5sselect_hyperslab_f(msp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)
         
!     describe the memory layout of the data again 
  dims(1)=2                  ! complex = 2 * real
  dims(2)=mynp12_read        ! read needs dimension of arrays
  dims(3)=mynalpha_read      ! here size of array in memory 

! read data 

#ifdef HDF_MPIO
      ! check for consistency with MPI
  blocksize=4_lint*2*mynp12_read*mynalpha_read  
  IF(blocksize.GE.blocklimit) THEN
    WRITE(*,*) 'block too large to be handled with MPI-IO'
    WRITE(*,*) 'blocksize = ',blocksize
    WRITE(*,*) 'blocklimit = ',blocklimit
    CALL abort
  END IF
#endif
  
  CALL h5dread_f(dset_nn_id, H5T_NATIVE_REAL,& 
     &      buf,dims,ierr,&
     &      mem_space_id = msp_hyper_id,&
     &      file_space_id = dsp_hyper_id,& 
     &      xfer_prp = pcollectread_id)
     
     
! here all data is read into buffers on myp3q4enerid=0      
  hdftime=hdftime+MPI_WTIME()   
! close hdf handles 
     ! close data set   
  CALL h5dclose_f(dset_nn_id,ierr)
  
! do not need mem_space anymore      
  CALL h5sclose_f(msp_hyper_id,ierr)
  CALL h5sclose_f(msp_nn_id,ierr)
! and not data space anymore    
  CALL h5sclose_f(dsp_hyper_id,ierr)
  CALL h5sclose_f(dsp_nn_id,ierr) 
! close the group for this amplitude 
  CALL h5gclose_f(group_amp_id,ierr) 
! finally close the file 
  CALL close_h5file(readf_id)

! need data on all PEs in commp3q4ener 
  CALL MPI_BCAST(buf,2*mynp12_read*mynalpha_read,MPI_REAL4,&
    &            0,commp3q4ener,ierr)
  
! set up the relation between current channels and the channels 
! read in 
  ALLOCATE(alphap(mynalphaNN))
  IF(cdepNN) THEN  ! use all channels 
    CALL set_alphann(qn_read,alphacdepmax_read,alphap,alphacdepmax_read,cdepNN)
  ELSE ! only use the charge independent ones 
    CALL set_alphann(qn_read,alphacdepmax_read,alphap,alphamax_read,cdepNN)
  END IF

!  devide buffer by proper l factor for interpolation
  
  DO alphaloc=1,mynalpha_read   
    alpha=myalphaid+1+npe_alpha*(alphaloc-1)
    l12=qn_read(1,alpha)
    
    DO ip12=1,mynp12_read
     fakt=max(p12p_read(myp12_read-1+ip12)**l12,tolerance)   
     buf(1:2,ip12,alphaloc)=buf(1:2,ip12,alphaloc)/fakt
    END DO ! ip12 
  END DO   ! alphaloc 
  
! prepare the spline elements to interpolate
! from p12p_read to p12p  
  ALLOCATE(spl12(4,mynp12),indx12(4,mynp12))
  CALL cubfast_dp(p12p_read,p12n_read,&
        &    p12p(myp12:myp12+mynp12-1),mynp12, &
        &    spl12,indx12)
  
! continue with interpolation of the data 
! collect all p12 mesh points to the processors 
  ALLOCATE(psi_tmp_p12(2,p12n_read,mynalpha_read))
  CALL collect_part_spreal(buf,psi_tmp_p12,2,p12n_read,mynalpha_read,commp12)
  DEALLOCATE(buf)
! then interpolate  
  ALLOCATE(psi_inter_p12(2,mynp12,mynalpha_read))
 
! now interpolation of psi/p12**l12
  
  psi_inter_p12=0.0
  DO alpha=1,mynalpha_read
   DO ip12=1,mynp12     
    psi_inter_p12(1:2,ip12,alpha) = &
     &    spl12(1,ip12)*psi_tmp_p12(1:2,indx12(1,ip12),alpha) &
     &  + spl12(2,ip12)*psi_tmp_p12(1:2,indx12(2,ip12),alpha) &
     &  + spl12(3,ip12)*psi_tmp_p12(1:2,indx12(3,ip12),alpha) &
     &  + spl12(4,ip12)*psi_tmp_p12(1:2,indx12(4,ip12),alpha)
   END DO
  END DO 
  DEALLOCATE(spl12,indx12)
  DEALLOCATE(psi_tmp_p12) 

! now collect all channels to reorder the channels 
! to current bookkeeping, drop channels that cannot be assigned 
! leave components zero that are not found in file   
! start with collecting all channels  
  
  ALLOCATE(psi_tmp_alpha(2,mynp12,alphacdepmax_read))
  CALL collect_piece_spreal(psi_inter_p12,psi_tmp_alpha, &
  &                  2*mynp12,alphacdepmax_read,1,commalpha)
  
  DEALLOCATE(psi_inter_p12)
  
  IF(.NOT. allocated(PSI%amp)) &
    &   ALLOCATE(PSI%amp(mynp12,mynalphaNN))
    
  PSI%amp=0.0   ! 
  DO alphaloc=1,mynalphaNN     
   IF(alphap(alphaloc).NE.0) THEN
    alpha=myalphaid+1+npe_alpha*(alphaloc-1)   
    call getNNqn(alpha,l12,s12,j12,t12,mt12)     
    DO ip12=1,mynp12
     fakt=p12p(myp12-1+ip12)**l12   
     PSI%amp(ip12,alphaloc) &
     &  = fakt*psi_tmp_alpha(1,ip12,alphap(alphaloc)) &
         + fakt*IU*psi_tmp_alpha(2,ip12,alphap(alphaloc))
    END DO
   END IF
  END DO
  
  DEALLOCATE(psi_tmp_alpha)
  DEALLOCATE(alphap)
  DEALLOCATE(p12p_read,p12w_read)
  DEALLOCATE(qn_read)
  
  readamptime=readamptime+MPI_WTIME()

 END SUBROUTINE readamp_nn

 SUBROUTINE set_alphann(qn,alphamaxcdep_read,alphap,alphamax_read,cdep)
  IMPLICIT NONE
  INTEGER alphamax_read,alphamaxcdep_read
  INTEGER qn(5,alphamaxcdep_read),alphap(mynalphaNN)
  LOGICAL cdep
  INTEGER alphaneu,alphaloc,alphaalt
  INTEGER l12,s12,j12,t12,mt12
  INTEGER l12p,s12p,j12p,t1p,t12p,mt12p
  
  alphap=0
  
  DO alphaloc=1,mynalphaNN
   alphaneu=myalphaid+1+npe_alpha*(alphaloc-1)
   call getNNqn(alphaneu,l12,s12,j12,t12,mt12)
   DO alphaalt=1,alphamax_read
    l12p =qn(1,alphaalt)
    s12p =qn(2,alphaalt)
    j12p =qn(3,alphaalt)
    t12p =qn(4,alphaalt)
    mt12p=qn(5,alphaalt)
    IF((l12p.eq.l12).AND.&
    &  (s12p.eq.s12).AND.&
    &  (j12p.eq.j12).AND.&
    &  (t12p.eq.t12).AND.&
    &  ((mt12p.eq.mt12).OR. (.NOT. cdep))) THEN
     alphap(alphaloc)=alphaalt
     EXIT
    END IF
    IF(alphaalt.EQ.alphamax_read &
   &  .and.myenerid.EQ.0  &
   &  .and.mymeshid.EQ.0) &
     &  WRITE(*,*) 'NN Channel not found:',alphaneu
   END DO
  END DO ! alphaloc
  
  
 END SUBROUTINE set_alphann
 
 ! 3N amplitude
 ! open file <filename>
 ! write amplitude PSI to file <filename>
 ! in the group ampname
 ! file is then closed 
 
 SUBROUTINE writeamp_3n(PSI,filename,ampname)
  IMPLICIT NONE 
  TYPE(AMP3N),TARGET :: PSI 
  CHARACTER(LEN=*) filename
  CHARACTER(LEN=*) ampname
  INTEGER(HID_T) writef_id
  INTEGER(HID_T) group_amp_id
  INTEGER(HID_T) dset_nn_id,dsp_nn_id,dsp_hyper_id
  INTEGER(HID_T) msp_nn_id,msp_hyper_id
  INTEGER(HSIZE_T) dims(4)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
  INTEGER(HSIZE_T) start(4)     ! storage for offset in output arrays
  INTEGER(HSIZE_T) block(4)     ! storage for blocksize in output arrays
  INTEGER(HSIZE_T) stride(4)    ! storage for stride in output arrays
  INTEGER(HSIZE_T) count(4)     ! storage for count in output arrays
  INTEGER ierr
  REAL(spreal),ALLOCATABLE  :: buf(:,:,:,:)
  INTEGER(lint) blocksize

  writeamptime=writeamptime-MPI_WTIME()
  
  ! first open the hdf file
  ! is always assumed that the file exits
  ! copy empty file if necessary
  
  CALL open_write_h5file(filename,commall,writef_id)
  
  ! then create the group <ampname>

  CALL h5gcreate_f(writef_id,trim(ampname),group_amp_id, ierr)

  ! file should contain all info to interprete the data 
  !    1. meshpoints for P12 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  CALL writep12p(group_amp_id)
   
  !    2. meshpoints for P3 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  CALL writep3p(group_amp_id)

  !    3. 3N channels qnalpha3N(:,:) including number of channels 
  CALL write_3n_channels(group_amp_id)
  
  !    4. physical constants from constant mod are also written 
  CALL writephysconst(group_amp_id) 
  
  !    5. 3N amplitude if it exists 
  
  hdftime=hdftime-MPI_WTIME()
  
  IF(allocated(PSI%amp)) THEN
  ! this needs to be written in parallel since the data 
  ! is distributed 
  ! make sure that that only one processor writes the data 
  ! in case the distribution is not complete 
  
  ! here data is distributed in in commp12 and commalpha 
  ! that writing processors should have myp3id=myq4id=myenerid=0
  ! put the hyperslabs for other processors to empty 
  
  !     need to create a data set in the file 
  !     with complete number amp-elements 
  !     outline shape of global set of data 
   dims(1)=2                  !  complex = 2 * real
   dims(2)=P12N          !  write needs dimensions of arrays 
   dims(3)=P3N
   dims(4)=alpha3Ncdepmax     !  here array size in data file  
  !     create corresponding data space            
   CALL h5screate_simple_f(4,dims,dsp_nn_id,ierr);

  !     describe the memory layout of the data
   dims(1)=2                  ! complex = 2 * real
   dims(2)=mynp12        ! write needs dimension of arrays
   dims(3)=mynp3
   dims(4)=mynalpha3N         ! here size of array in memory 
  !     create corresponding mem space    
   CALL h5screate_simple_f(4,dims,msp_nn_id,ierr);
  
  ! select the hyperslap in the file  that corresponds to the local 
  ! 3N data 
   IF(myq4id.EQ.0 .AND. myenerid.EQ.0) THEN
    start(1)=0
    start(2)=myp12-1
    start(3)=myp3-1
    start(4)=myalphaid 
    block(1)=2
    block(2)=mynp12
    block(3)=mynp3
    block(4)=1
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=npe_alpha
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=mynalpha3N
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF
  
   CALL h5scopy_f(dsp_nn_id,dsp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(dsp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)

  ! select the hyperslap in memorythat corresponds to the local 
  ! 3N data 

   IF(myq4id.EQ.0 .AND. myenerid.EQ.0) THEN
    start(1)=0
    start(2)=0             
    start(3)=0
    start(4)=0
    block(1)=2
    block(2)=mynp12
    block(3)=mynp3
    block(4)=mynalpha3N
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=1
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=1
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF   
     
   CALL h5scopy_f(msp_nn_id,msp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(msp_hyper_id,H5S_SELECT_SET_F,&
      &                      start, count,ierr,stride,block)
       
      
!    create the dataset in file 
   CALL h5dcreate_f(group_amp_id,trim('3namp'),H5T_NATIVE_REAL,&
      &             dsp_nn_id,dset_nn_id,ierr)
    

      
! resolve compatiblity issue by copying to REAL buffer 
   ALLOCATE(buf(2,mynp12,mynp3,mynalpha3N))
   buf(1,:,:,:)=REAL(PSI%amp(:,:,:))
   buf(2,:,:,:)=AIMAG(PSI%amp(:,:,:))

#ifdef HDF_MPIO
      ! check for consistency with MPI
  blocksize=4_lint*2*mynp12*mynp3*mynalpha3N  
  IF(blocksize.GE.blocklimit) THEN
    WRITE(*,*) 'block too large to be handled with MPI-IO'
    WRITE(*,*) 'blocksize = ',blocksize
    WRITE(*,*) 'blocklimit = ',blocklimit
    CALL abort
  END IF
#endif
 
      
! write data of amplitude (use here collective communication)       
   CALL h5dwrite_f(dset_nn_id, H5T_NATIVE_REAL,& 
     &      buf,dims,ierr,&
     &      mem_space_id = msp_hyper_id,&
     &      file_space_id = dsp_hyper_id,& 
     &      xfer_prp = pcollectwrite_id)

   DEALLOCATE(buf)

! close data set   
   CALL h5dclose_f(dset_nn_id,ierr)
  
! do not need mem_space anymore      
   CALL h5sclose_f(msp_hyper_id,ierr)
   CALL h5sclose_f(msp_nn_id,ierr)
! and not data space anymore    
   CALL h5sclose_f(dsp_hyper_id,ierr)
   CALL h5sclose_f(dsp_nn_id,ierr) 
  END IF  ! 3N exists 
  
  hdftime=hdftime+MPI_WTIME()
  
! close the group for this amplitude 
   CALL h5gclose_f(group_amp_id,ierr) 
! finally close the file 
   CALL close_h5file(writef_id)  
  
  writeamptime=writeamptime+MPI_WTIME()
  
 END SUBROUTINE writeamp_3n

 ! 3N  amplitude
 ! open file <filename>
 ! read amplitude PSI from group ampname 
 ! of file <filename>
 ! file is then closed
 ! print the parameters and give back the amplitude
 
 SUBROUTINE readamp_3n(PSI,filename,ampname)
  IMPLICIT NONE 
  TYPE(AMP3N),TARGET :: PSI 
  CHARACTER(LEN=*) filename
  CHARACTER(LEN=*) ampname
  INTEGER(HID_T) readf_id
  INTEGER(HID_T) group_amp_id
  INTEGER(HID_T) dset_nn_id,dsp_nn_id,dsp_hyper_id
  INTEGER(HID_T) msp_nn_id,msp_hyper_id
  INTEGER(HSIZE_T) dims(4)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
  INTEGER(HSIZE_T) start(4)     ! storage for offset in output arrays
  INTEGER(HSIZE_T) block(4)     ! storage for blocksize in output arrays
  INTEGER(HSIZE_T) stride(4)    ! storage for stride in output arrays
  INTEGER(HSIZE_T) count(4)     ! storage for count in output arrays
  INTEGER ierr,p12n_read,p3n_read
  INTEGER alpha3Ncdepmax_read,alpha3Nmax_read
  INTEGER ip12,ip3,alpha
  ! buffer to read the complex data 
  REAL(spreal),ALLOCATABLE  :: buf(:,:,:,:),&
     &                         psi_tmp_p12(:,:,:,:), &
     &                         psi_inter_p12(:,:,:,:),&
     &                         psi_tmp_p3(:,:,:,:), &
     &                         psi_inter_p3(:,:,:,:),&
     &                         psi_tmp_alpha(:,:,:,:)
  REAL(dpreal),ALLOCATABLE  :: p12p_read(:),p12w_read(:),&
                       &       p3p_read(:),p3w_read(:)
  INTEGER,ALLOCATABLE  :: qnalpha3N_read(:,:)
  INTEGER mynp12_read,myp12_read,myendp12_read
  INTEGER mynp3_read,myp3_read,myendp3_read
  INTEGER mynalpha3N_read
  REAL(dpreal),ALLOCATABLE :: spl12(:,:),spl3(:,:)
  INTEGER,ALLOCATABLE :: indx12(:,:),indx3(:,:),alphap(:)
  LOGICAL amp_exists
  REAL(dpreal) :: fakt 
  INTEGER alphaloc,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  INTEGER(lint) blocksize

  readamptime=readamptime-MPI_WTIME()
  
  ! first open the hdf file
  ! is always assumed that the file exits
  ! copy empty file if necessary
  CALL open_read_h5file(filename,commall,readf_id)
  
  ! then create the group <ampname>
  
  CALL h5gopen_f(readf_id,trim(ampname),group_amp_id, ierr)

  ! file should contain all info to interprete the data 
  !    1. meshpoints for P12/P3 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  !       the meshpoints and weights are allocated within the routine
  CALL readp12p(group_amp_id,p12p_read,p12w_read,p12n_read)
  CALL readp3p(group_amp_id,p3p_read,p3w_read,p3n_read)

  !    2. 3N/3N channels qn_read(:,:) including number of channels
  !       the quantum numbers are allocated within the routine 
  CALL read_3n_channels(group_amp_id,qnalpha3N_read,&
     &                   alpha3Nmax_read,alpha3Ncdepmax_read)
  
  !    3. physical constants from constant mod are also read 
  !       and printed. No values are used. 
  CALL readphysconst(group_amp_id) 

! prepare the spline elements to interpolate
! from p12p_read to p12p 
! and p3p_read to p3p  
   
  ALLOCATE(spl12(4,mynp12),indx12(4,mynp12))
  ALLOCATE(spl3(4,mynp3),indx3(4,mynp3))
  CALL cubfast_dp(p12p_read,p12n_read,&
        &    p12p(myp12:myp12+mynp12-1),mynp12, &
        &    spl12,indx12)
  CALL cubfast_dp(p3p_read,p3n_read,&
        &    p3p(myp3:myp3+mynp3-1),mynp3, &
        &    spl3,indx3)
        
  ! distribution of mesh and channels 
  CALL distr_block(p12n_read,npe_p12,myp12id,myp12_read,myendp12_read,mynp12_read)
  CALL distr_block(p3n_read,npe_p3,myp3id,myp3_read,myendp3_read,mynp3_read)
  CALL distr_piece(alpha3Ncdepmax_read,npe_alpha,myalphaid,mynalpha3N_read)
  
  !    4. 3N amplitude if it exists 
  
  CALL h5lexists_f(group_amp_id, trim('3namp'),amp_exists,ierr)
  
  hdftime=hdftime-MPI_WTIME()
  
  IF(amp_exists) THEN
  ! first read the amplitude using the mesh stored 
  ! prepare buffer for real and imaginary part of the amplitude 
  
   ALLOCATE(buf(2,mynp12_read,mynp3_read,mynalpha3N_read))

  
  ! here data is distributed in in commp12 and commalpha 
  ! the reading processors should have myp3id=myq4id=myenerid=0
  ! this is element 0 in commp3q4ener
  ! put the hyperslabs for other processors to empty 

  !     need to open the data set in the file 
   CALL h5dopen_f(group_amp_id,trim('3namp'),dset_nn_id,ierr)
  
  ! get the data space descriptor 
   CALL h5dget_space_f(dset_nn_id, dsp_nn_id, ierr)
     
  !     describe the memory layout of the data
   dims(1)=2                  ! complex = 2 * real
   dims(2)=mynp12_read        ! read needs dimension of arrays
   dims(3)=mynp3_read
   dims(4)=mynalpha3N_read      ! here size of array in memory 
  !     create corresponding mem space    
   CALL h5screate_simple_f(4,dims,msp_nn_id,ierr);
  
  ! select the hyperslap in the file  that corresponds to the local 
  ! 3N data 
   IF(myq4enerid.EQ.0) THEN
    start(1)=0
    start(2)=myp12_read-1
    start(3)=myp3_read-1
    start(4)=myalphaid 
    block(1)=2
    block(2)=mynp12_read
    block(3)=mynp3_read
    block(4)=1
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=npe_alpha
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=mynalpha3N_read
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF
  
   CALL h5scopy_f(dsp_nn_id,dsp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(dsp_hyper_id,H5S_SELECT_SET_F,&
      &                      start, count,ierr,stride,block)

  ! select the hyperslap in memorythat corresponds to the local 
  ! 3N data 
   IF(myq4enerid.EQ.0) THEN
    start(1)=0
    start(2)=0             
    start(3)=0
    start(4)=0
    block(1)=2
    block(2)=mynp12_read
    block(3)=mynp3_read 
    block(4)=mynalpha3N_read
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=1
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=1
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF   
     
   CALL h5scopy_f(msp_nn_id,msp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(msp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)
         
!     describe the memory layout of the data again 
   dims(1)=2                  ! complex = 2 * real
   dims(2)=mynp12_read        ! read needs dimension of arrays
   dims(3)=mynp3_read
   dims(4)=mynalpha3N_read      ! here size of array in memory  

#ifdef HDF_MPIO
      ! check for consistency with MPI
   blocksize=4_lint*2*mynp12_read*mynp3_read*mynalpha3N_read  
   IF(blocksize.GE.blocklimit) THEN
    WRITE(*,*) 'block too large to be handled with MPI-IO'
    WRITE(*,*) 'blocksize = ',blocksize
    WRITE(*,*) 'blocklimit = ',blocklimit
    CALL abort
   END IF
#endif
  
   
! read data 

   CALL h5dread_f(dset_nn_id, H5T_NATIVE_REAL,& 
     &      buf,dims,ierr,&
     &      mem_space_id = msp_hyper_id,&
     &      file_space_id = dsp_hyper_id,& 
     &      xfer_prp = pcollectread_id)
     
     
! here all data is read into buffers on myq4enerid=0      
     
! close hdf handles 
     ! close data set   
   CALL h5dclose_f(dset_nn_id,ierr)
  
! do not need mem_space anymore      
   CALL h5sclose_f(msp_hyper_id,ierr)
   CALL h5sclose_f(msp_nn_id,ierr)
! and not data space anymore    
   CALL h5sclose_f(dsp_hyper_id,ierr)
   CALL h5sclose_f(dsp_nn_id,ierr) 

! need data on all PEs in commp3q4ener 
   CALL MPI_BCAST(buf,2*mynp12_read*mynp3_read*mynalpha3N_read,&
    &            MPI_REAL4,0,commq4ener,ierr)
 
! set up the relation between current channels and the channels 
! read in 
   ALLOCATE(alphap(mynalpha3N))
   IF(cdep3N) THEN  ! use all channels 
    CALL set_alpha_3n(qnalpha3N_read,alpha3Ncdepmax_read,alphap,alpha3Ncdepmax_read,cdep3N)
   ELSE ! only use the charge independent ones 
    CALL set_alpha_3n(qnalpha3N_read,alpha3Ncdepmax_read,alphap,alpha3Nmax_read,cdep3N)
   END IF
   

   !  devide buffer by proper l factor for interpolation
  
  DO alphaloc=1,mynalpha3N_read   
    alpha=myalphaid+1+npe_alpha*(alphaloc-1)
    l12=qnalpha3N_read(1,alpha)
    l3=qnalpha3N_read(5,alpha)
    
    DO ip3=1,mynp3_read
     DO ip12=1,mynp12_read   
      fakt=max(p12p_read(myp12_read-1+ip12)**l12*p3p_read(myp3_read-1+ip3)**l3,tolerance)   
      buf(1:2,ip12,ip3,alphaloc)=buf(1:2,ip12,ip3,alphaloc)/fakt
     END DO ! ip12 
    END DO  ! ip3 
  END DO   ! alphaloc 
  
! continue with interpolation of the data 
! collect all p12 mesh points to the processors 
   ALLOCATE(psi_tmp_p12(2,p12n_read,mynp3_read,mynalpha3N_read))
   CALL collect_part_spreal(buf,psi_tmp_p12,2,p12n_read,&
        &                   mynalpha3N_read*mynp3_read,commp12)
   DEALLOCATE(buf)
! then interpolate  
   ALLOCATE(psi_inter_p12(2,mynp12,mynp3_read,mynalpha3N_read))
 
! now interpolation of psi/p12**l12/p3**l3
     
   psi_inter_p12=0.0
   DO alpha=1,mynalpha3N_read
    DO ip3=1,mynp3_read
     DO ip12=1,mynp12     
      psi_inter_p12(1:2,ip12,ip3,alpha) = &
     &    spl12(1,ip12)*psi_tmp_p12(1:2,indx12(1,ip12),ip3,alpha) &
     &  + spl12(2,ip12)*psi_tmp_p12(1:2,indx12(2,ip12),ip3,alpha) &
     &  + spl12(3,ip12)*psi_tmp_p12(1:2,indx12(3,ip12),ip3,alpha) &
     &  + spl12(4,ip12)*psi_tmp_p12(1:2,indx12(4,ip12),ip3,alpha)
     END DO
    END DO
   END DO 
   DEALLOCATE(psi_tmp_p12) 

! continue with interpolation of the p3 data 
! collect all p3 mesh points to the processors 
   ALLOCATE(psi_tmp_p3(2,mynp12,p3n_read,mynalpha3N_read))
   CALL collect_part_spreal(psi_inter_p12,psi_tmp_p3,&
        &                   2*mynp12,p3n_read,&
        &                   mynalpha3N_read,commp3)
   DEALLOCATE(psi_inter_p12)
! then interpolate  
   ALLOCATE(psi_inter_p3(2,mynp12,mynp3,mynalpha3N_read))
 
! now interpolation of psi/p12**l12/p3**l3
     
   psi_inter_p3=0.0
   DO alpha=1,mynalpha3N_read
    DO ip3=1,mynp3
     psi_inter_p3(1:2,1:mynp12,ip3,alpha) = &
     &    spl3(1,ip3)*psi_tmp_p3(1:2,1:mynp12,indx3(1,ip3),alpha) &
     &  + spl3(2,ip3)*psi_tmp_p3(1:2,1:mynp12,indx3(2,ip3),alpha) &
     &  + spl3(3,ip3)*psi_tmp_p3(1:2,1:mynp12,indx3(3,ip3),alpha) &
     &  + spl3(4,ip3)*psi_tmp_p3(1:2,1:mynp12,indx3(4,ip3),alpha)
    END DO
   END DO 
   DEALLOCATE(psi_tmp_p3) 
  
! now collect all channels to reorder the channels 
! to current bookkeeping, drop channels that cannot be assigned 
! leave components zero that are not found in file   
! start with collecting all channels  
  
   ALLOCATE(psi_tmp_alpha(2,mynp12,mynp3,alpha3Ncdepmax_read))
   CALL collect_piece_spreal(psi_inter_p3,psi_tmp_alpha, &
  &                  2*mynp12*mynp3,alpha3Ncdepmax_read,1,commalpha)
  
   DEALLOCATE(psi_inter_p3)
  
   IF(.NOT. allocated(PSI%amp)) &
    &   ALLOCATE(PSI%amp(mynp12,mynp3,mynalpha3N))
    
   PSI%amp=0.0   ! 
   DO alphaloc=1,mynalpha3N     
    IF(alphap(alphaloc).NE.0) THEN
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)   
     CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)
     DO ip3=1,mynp3
      DO ip12=1,mynp12
      fakt=p12p(myp12-1+ip12)**l12*p3p(myp3-1+ip3)**l3   
      PSI%amp(ip12,ip3,alphaloc) &
     &  = fakt*psi_tmp_alpha(1,ip12,ip3,alphap(alphaloc)) &
         + fakt*IU*psi_tmp_alpha(2,ip12,ip3,alphap(alphaloc))
      END DO
     END DO
    END IF
   END DO
  
   DEALLOCATE(psi_tmp_alpha)
   DEALLOCATE(alphap)
  END IF  ! 3N exists 
  
  hdftime=hdftime+MPI_WTIME()   
  
  DEALLOCATE(spl12,indx12)
  DEALLOCATE(spl3,indx3)
  DEALLOCATE(p12p_read,p12w_read)
  DEALLOCATE(p3p_read,p3w_read)
  DEALLOCATE(qnalpha3N_read)
 
! close the group for this amplitude 
  CALL h5gclose_f(group_amp_id,ierr) 
! finally close the file 
  CALL close_h5file(readf_id)  
  readamptime=readamptime+MPI_WTIME()
   
 END SUBROUTINE readamp_3n

 SUBROUTINE set_alpha_3n(qn,alphamaxcdep_read,alphap,alphamax_read,cdep)
  IMPLICIT NONE
  INTEGER alphamax_read,alphamaxcdep_read
  INTEGER qn(11,alphamaxcdep_read),alphap(mynalpha3N)
  LOGICAL cdep
  INTEGER alphaneu,alphaloc,alphaalt
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  INTEGER l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,mtau3p
  
  alphap=0
  
  DO alphaloc=1,mynalpha3N
   alphaneu=myalphaid+1+npe_alpha*(alphaloc-1)
   call get3Nqn(alphaneu,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)
   DO alphaalt=1,alphamax_read
    l12p =qn(1,alphaalt)
    s12p =qn(2,alphaalt)
    j12p =qn(3,alphaalt)
    t12p =qn(4,alphaalt)
    l3p =qn(5,alphaalt)
    I3p =qn(6,alphaalt)
    j3p =qn(7,alphaalt)
    tau3p =qn(8,alphaalt)
    mtau3p =qn(9,alphaalt)

    IF((l12p.eq.l12).AND.&
    &  (s12p.eq.s12).AND.&
    &  (j12p.eq.j12).AND.&
    &  (t12p.eq.t12).AND.&
    &  (l3p.eq.l3).AND.&
    &  (I3p.eq.I3).AND.&
    &  (j3p.eq.j3).AND.&
    &  (tau3p.eq.tau3).AND.&
    &  ((mtau3p.eq.mtau3).OR. (.NOT. cdep))) THEN
     alphap(alphaloc)=alphaalt
     EXIT
    END IF
    IF(alphaalt.EQ.alphamax_read &
   &  .and.myenerid.EQ.0  &
   &  .and.mymeshid.EQ.0) THEN
     WRITE(*,*) '3N Channel not found:',alphaneu
    ENDIF
   END DO
  END DO ! alphaloc
  
 END SUBROUTINE set_alpha_3n
 
 ! for testing print out a few simple expectation values 
 ! 1) kinetic energy for each component and  norm for each component 
 ! 2) LS probabilities for each component 
 ! first case: 3N wave function
 
 
 SUBROUTINE printekin_3n_wave(PSI)
  IMPLICIT NONE
  TYPE(AMP3N) PSI
  REAL(spreal) :: locekin,ekin,locnorm,norm
  INTEGER ip12,ip3,alphaloc,alpha,ierr
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  INTEGER alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p, &
 &        mtau3p,alpha2bp,parip,mt12,mt1
  COMPLEX(spreal),ALLOCATABLE :: psiamp3N(:,:,:)
  REAL(spreal),ALLOCATABLE :: isomass(:,:,:)
  REAL(spreal),ALLOCATABLE :: mbaryon(:,:),mthres(:)
  REAL(spreal) :: isofakt 
  
!  prepare masses depending on isospin 
  ALLOCATE(mbaryon(0:2,-2:2),mthres(-4:4))
  mbaryon=0.0_spreal
! set all baryons, only nucleons are used   
  mbaryon(0,0)=mlam
  mbaryon(1,-1)=mneu
  mbaryon(1,1)=mprot
  mbaryon(2,-2)=msigm
  mbaryon(2,0)=msig0
  mbaryon(2,2)=msigp
  mthres=0.0_spreal
  mthres(1)=2.0*mprot+mneu
  mthres(-1)=2.0*mneu+mprot
  
  
  IF(allocated(PSI%amp)) THEN ! calculate based on NNY amplitude 
! start with kinetic energy
! since there is a possibility to have isospin transitions
! I need to collect all channels 
   ALLOCATE(psiamp3N(mynp12,mynp3,alpha3Ncdepmax))
   CALL collect_piece_spcmplx(PSI%amp,psiamp3N,&
          &     mynp12*mynp3,alpha3Ncdepmax,1,commalpha)

! prepare the isospin masses for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    isomass(i,alphap,alpha)
!       i=1   total mass 
!       i=2   reduced mass average particles 12 
!       i=3   reduced mass average particle 3          
          
          
   ALLOCATE(isomass(1:3,alpha3Ncdepmax,mynalpha3N)) 
   isomass=0.0_spreal
   
   DO alphap=1,alpha3Ncdepmax
    CALL get3Nqn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,&
       &          j3p,tau3p,mtau3p,alpha2bp,parip) 
    t12p=2*t12p   !   use convention that isospins are twice    
    DO alphaloc=1,mynalpha3N
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)
     t12=2*t12   !   use convention that isospins are twice
     IF(l12.eq.l12p .and. s12.eq.s12p .and.&
      & j12.eq.j12p .and.  l3.eq.l3p  .and.&
      &  I3.eq.I3p  .and.  j3.eq.j3p  .and.&
      &  mtau3.eq.mtau3p) THEN
      DO mt12=-t12,t12,2
       IF(abs(mtau3-mt12).GT.1) cycle
       DO mt1=-1,1,2 
        IF(abs(mt12-mt1).GT.1) cycle
        isofakt = CG(t12,1,tau3,mt12,mtau3-mt12,mtau3) &
        & *CG(1,1,t12,mt1,mt12-mt1,mt12) &
        & *CG(t12p,1,tau3p,mt12,mtau3-mt12,mtau3) &
        & *CG(1,1,t12p,mt1,mt12-mt1,mt12)
        
        isomass(1,alphap,alphaloc) = isomass(1,alphap,alphaloc) &
        & + isofakt &
        &   *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12) &
        &    - mthres(mtau3))
        isomass(2,alphap,alphaloc) = isomass(2,alphap,alphaloc) &
        & + isofakt &
        &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  & 
        &      /(2.0*mbaryon(1,mt1)*mbaryon(1,mt12-mt1))
        isomass(3,alphap,alphaloc) = isomass(3,alphap,alphaloc) &
        & + isofakt &
        &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12))  & 
        &      /(2.0*(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  &
                    *mbaryon(1,mtau3-mt12))
       END DO ! mt1 
      END DO  ! mt12   
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   locekin=0.0_dpreal 
   DO alphap=1,alpha3Ncdepmax
    DO alphaloc=1,mynalpha3N
     IF(isomass(2,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
      DO ip3=1,mynp3
       DO ip12=1,mynp12
        locekin=locekin &
       & + conjg(psiamp3N(ip12,ip3,alphap)) &
            *PSI%amp(ip12,ip3,alphaloc) &
       &    *p12p3weight(ip12,ip3) &
       &    *(isomass(1,alphap,alphaloc) &
       &    + isomass(2,alphap,alphaloc)*P12P(myp12+ip12-1)**2 & 
       &    + isomass(3,alphap,alphaloc)*P3P(myp3+ip3-1)**2)
       
       END DO ! ip12 
      END DO ! ip3 
     END IF
    END DO ! alphaloc
   END DO ! alphap

   locnorm=0.0_dpreal 
   DO alphaloc=1,mynalpha3N
    DO ip3=1,mynp3
     DO ip12=1,mynp12
      locnorm=locnorm &
     & + conjg(PSI%amp(ip12,ip3,alphaloc)) &
          *PSI%amp(ip12,ip3,alphaloc) &
     &    *p12p3weight(ip12,ip3) 
     END DO ! ip12 
    END DO ! ip3 
   END DO ! alphaloc
   
   DEALLOCATE(psiamp3N,isomass) 
   
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locekin,ekin,1,MPI_REAL4,MPI_SUM,commamp3n,ierr) ! here it is assumed that product is REAL(spreal)
   CALL MPI_ALLREDUCE(locnorm,norm,1,MPI_REAL4,MPI_SUM,commamp3n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A,E15.6)') 'kinetic energy from 3N wave [MeV]: ',ekin*hbarc
    WRITE(*,'(A,E15.6)') 'Norm of 3N wave amplitude:         ',norm
   END IF
   
  END IF   ! 3N wave 
  
  DEALLOCATE(mbaryon,mthres) 
  
 END SUBROUTINE printekin_3n_wave
 
 SUBROUTINE printls_3n_wave(PSI)
  IMPLICIT NONE
  TYPE(AMP3N) PSI
  REAL(spreal),ALLOCATABLE :: locpls(:,:),pls(:,:)
  INTEGER ip12,ip3,alphaloc,alpha,ierr
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  INTEGER alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p, &
 &        mtau3p,alpha2bp,parip,mt12,mt1,bl,bs
  COMPLEX(spreal),ALLOCATABLE :: psiamp3N(:,:,:)
  REAL(spreal),ALLOCATABLE :: recoupl(:,:,:,:) 
  REAL(spreal) :: norm 
  

  
  IF(allocated(PSI%amp)) THEN ! calculate based on 3N wave 
! start with kinetic energy
! since there is a possibility to have isospin transitions
! I need to collect all channels 
   ALLOCATE(psiamp3N(mynp12,mynp3,alpha3Ncdepmax))
   CALL collect_piece_spcmplx(PSI%amp,psiamp3N,&
          &     mynp12*mynp3,alpha3Ncdepmax,1,commalpha)

! prepare the recoupl coeffs for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    recoupl(L,S,alphap,alpha)
!       L total angular momentum 0:j3_max+3/2 
!       S total spin 0  =  1/2  1 = 3/2           
          
          
   ALLOCATE(recoupl(0:(j3max+3)/2,0:1,alpha3Ncdepmax,mynalpha3N))
   
   recoupl=0.0_spreal
   
   DO alphap=1,alpha3Ncdepmax
    call get3Nqn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,&
       &          j3p,tau3p,mtau3p,alpha2bp,parip)
    DO alphaloc=1,mynalpha3N
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)
     IF(t12.eq.t12p .and. &
      & tau3.eq.tau3p .and.  mtau3.eq.mtau3p .and. &
      & l12.eq.l12p .and. l3.eq.l3p .and. &
      & s12.eq.s12p .and. j3.eq.j3p) THEN
      DO bs=1,3,2 
       DO bl=abs(j3-bs)/2,(j3+bs)/2        
        recoupl(bl,bs/2,alphap,alphaloc) = & 
        &  SDCH(bl)*SDCH(j12)*SD2CH(I3)*SD2CH(bs)  &
        &   *C9J(2*l12,2*s12,2*j12,2*l3,1,I3,2*bl,bs,j3) &
        &   *SDCH(bl)*SDCH(j12p)*SD2CH(I3p)*SD2CH(bs)  &
        &   *C9J(2*l12,2*s12,2*j12p,2*l3,1,I3p,2*bl,bs,j3)
       END DO ! bs
      END DO  ! bl   
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   ALLOCATE(locpls(0:(j3max+3)/2,0:1),pls(0:(j3max+3)/2,0:1))
   
   locpls=0.0_dpreal 
   DO bl=0,(j3max+3)/2
    DO bs=1,3,2 
     DO alphap=1,alpha3Ncdepmax
      DO alphaloc=1,mynalpha3N
       IF(recoupl(bl,bs/2,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
        DO ip3=1,mynp3
         DO ip12=1,mynp12
          locpls(bl,bs/2)=locpls(bl,bs/2) &
       & + conjg(psiamp3N(ip12,ip3,alphap)) &
            *PSI%amp(ip12,ip3,alphaloc) &
       &    *p12p3weight(ip12,ip3) &
       &    *recoupl(bl,bs/2,alphap,alphaloc)
         END DO ! ip12 
        END DO ! ip3 
       END IF
      END DO ! alphaloc
     END DO ! alphap
    END DO ! bs 
   END DO ! bl 
   
   
   
   DEALLOCATE(psiamp3N,recoupl) 
   
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locpls,pls,j3max+5,MPI_REAL4,MPI_SUM,commamp3n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A)') 'LS probabilities from 3N wave: '
    norm=0.0_spreal
    DO bl=0,(j3max+3)/2
     DO bs=1,3,2 
      IF(pls(bl,bs/2).NE.0.0_spreal) THEN
       WRITE(*,'(A,2I4,2X,E15.6)') 'PLS-3N: ',bl,bs,pls(bl,bs/2)
       norm=norm+pls(bl,bs/2)
      END IF
     END DO
    END DO
    WRITE(*,'(A,8X,2X,E15.6)') 'PLS-3N: ',norm
   END IF
   
   DEALLOCATE(locpls,pls)
   
  END IF   ! 3N ampl  
  
 END SUBROUTINE printls_3n_wave
 
 SUBROUTINE printiso_3n_wave(PSI)
  IMPLICIT NONE
  TYPE(AMP3N) PSI
  REAL(spreal),ALLOCATABLE :: locpiso(:),piso(:)
  INTEGER ip12,ip3,alphaloc,alpha,ierr
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  INTEGER alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p, &
 &        mtau3p,alpha2bp,parip,pncase
  COMPLEX(spreal),ALLOCATABLE :: psiamp3N(:,:,:)
  REAL(spreal),ALLOCATABLE :: recoupl(:,:,:) 
  REAL(spreal) :: norm 
  
! start with projection on spectator proton and neutron in 
! case 1 and 2 of recoupl
  
  IF(allocated(PSI%amp)) THEN ! calculate based on 3N wave 
! start with kinetic energy
! since there is a possibility to have isospin transitions
! I need to collect all channels 
   ALLOCATE(psiamp3N(mynp12,mynp3,alpha3Ncdepmax))
   CALL collect_piece_spcmplx(PSI%amp,psiamp3N,&
          &     mynp12*mynp3,alpha3Ncdepmax,1,commalpha)

! prepare the recoupl coeffs for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    recoupl(L,S,alphap,alpha)
!       L total angular momentum 0:j3_max+3/2 
!       S total spin 0  =  1/2  1 = 3/2           
          
          
   ALLOCATE(recoupl(2,alpha3Ncdepmax,mynalpha3N))
   
   recoupl=0.0_spreal
   
   DO alphap=1,alpha3Ncdepmax
    call get3Nqn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,&
       &          j3p,tau3p,mtau3p,alpha2bp,parip)
    DO alphaloc=1,mynalpha3N
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)
     IF(t12.eq.t12p .and. &
      & mtau3.eq.mtau3p .and. &
      & l12.eq.l12p .and. l3.eq.l3p .and. &
      & s12.eq.s12p .and. j3.eq.j3p .and. &
      &  I3.eq.I3p  .and.j12.eq.j12p) THEN
      DO pncase=1,2
       IF(tau3.eq.tau3p) THEN
        recoupl(pncase,alphap,alphaloc) = 0.5_dpreal 
       END IF
       recoupl(pncase,alphap,alphaloc) &
        &   = recoupl(pncase,alphap,alphaloc) &
       &   + (-1.0)**(pncase+1)*1.5*sqrt(2.0) &
        &   *SDCH(t12p)*SD2CH(tau3) &
        &   * C9J(2*t12p,2*t12,0,1,1,2,tau3p,tau3,2) &
        &   * CG(tau3,2,tau3p,mtau3,0,mtau3p)
      END DO  ! pncase   
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   ALLOCATE(locpiso(2),piso(2))
   
   locpiso=0.0_dpreal 
   DO pncase=1,2
    DO alphap=1,alpha3Ncdepmax
     DO alphaloc=1,mynalpha3N
      IF(recoupl(pncase,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
       DO ip3=1,mynp3
        DO ip12=1,mynp12
         locpiso(pncase)=locpiso(pncase) &
      & + conjg(psiamp3N(ip12,ip3,alphap)) &
      &    *PSI%amp(ip12,ip3,alphaloc) &
      &    *p12p3weight(ip12,ip3) &
      &    *recoupl(pncase,alphap,alphaloc)
        END DO ! ip12 
       END DO ! ip3 
      END IF
     END DO ! alphaloc
    END DO ! alphap
   END DO ! pncase
   
   
   
   DEALLOCATE(psiamp3N,recoupl) 
   
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locpiso,piso,2,MPI_REAL4,MPI_SUM,commamp3n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A)') 'ISO probabilities from 3N wave: '
    norm=0.0_spreal
    IF(piso(1).NE.0.0_spreal) THEN
     WRITE(*,'(A,2X,E15.6)') 'Spectator proton (WAVE):  ',piso(1)
     norm=norm+piso(1)
    END IF
     
    IF(piso(2).NE.0.0_spreal) THEN
     WRITE(*,'(A,2X,E15.6)') 'Spectator neutron (WAVE): ',piso(2)
     norm=norm+piso(2)
    END IF
    
    WRITE(*,'(A,2X,E15.6)')  'Sum of Spectator(WAVE):   ',norm
   END IF
   
   DEALLOCATE(locpiso,piso)
   
  END IF   ! 3N ampl  
  
 END SUBROUTINE printiso_3n_wave

 
 
 SUBROUTINE printinfo_3n_wave(PSI)
  IMPLICIT NONE
  TYPE(AMP3N) PSI 
  
  CALL printekin_3n_wave(PSI)
  CALL printls_3n_wave(PSI)
  CALL printiso_3n_wave(PSI)
  
 END SUBROUTINE printinfo_3n_wave
 
! first case: 3N wave function & fad component 
 
 
 SUBROUTINE printekin_3n(PSI,PSIFAD)
  IMPLICIT NONE
  TYPE(AMP3N) PSI,PSIFAD
  REAL(spreal) :: locekin,ekin,locnorm,norm
  INTEGER ip12,ip3,alphaloc,alpha,ierr
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  INTEGER alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p, &
 &        mtau3p,alpha2bp,parip,mt12,mt1
  COMPLEX(spreal),ALLOCATABLE :: psiamp3N(:,:,:)
  REAL(spreal),ALLOCATABLE :: isomass(:,:,:)
  REAL(spreal),ALLOCATABLE :: mbaryon(:,:),mthres(:)
  REAL(spreal) :: isofakt 
  
!  prepare masses depending on isospin 
  ALLOCATE(mbaryon(0:2,-2:2),mthres(-4:4))
  mbaryon=0.0_spreal
! set all baryons, only nucleons are used   
  mbaryon(0,0)=mlam
  mbaryon(1,-1)=mneu
  mbaryon(1,1)=mprot
  mbaryon(2,-2)=msigm
  mbaryon(2,0)=msig0
  mbaryon(2,2)=msigp
  mthres=0.0_spreal
  mthres(1)=2.0*mprot+mneu
  mthres(-1)=2.0*mneu+mprot
  
  
  IF(allocated(PSI%amp).AND. allocated(PSIFAD%amp)) THEN ! calculate based on 3N wave and fadcomp 
! start with kinetic energy
! since there is a possibility to have isospin transitions
! I need to collect all channels 
   ALLOCATE(psiamp3N(mynp12,mynp3,alpha3Ncdepmax))
   CALL collect_piece_spcmplx(PSIFAD%amp,psiamp3N,&
          &     mynp12*mynp3,alpha3Ncdepmax,1,commalpha)

! prepare the isospin masses for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    isomass(i,alphap,alpha)
!       i=1   total mass 
!       i=2   reduced mass average particles 12 
!       i=3   reduced mass average particle 3          
          
          
   ALLOCATE(isomass(1:3,alpha3Ncdepmax,mynalpha3N)) 
   isomass=0.0_spreal
   
   DO alphap=1,alpha3Ncdepmax
    CALL get3Nqn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,&
       &          j3p,tau3p,mtau3p,alpha2bp,parip) 
    t12p=2*t12p   !   use convention that isospins are twice    
    DO alphaloc=1,mynalpha3N
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)
     t12=2*t12   !   use convention that isospins are twice
     IF(l12.eq.l12p .and. s12.eq.s12p .and.&
      & j12.eq.j12p .and.  l3.eq.l3p  .and.&
      &  I3.eq.I3p  .and.  j3.eq.j3p  .and.&
      &  mtau3.eq.mtau3p) THEN
      DO mt12=-t12,t12,2
       IF(abs(mtau3-mt12).GT.1) cycle
       DO mt1=-1,1,2 
        IF(abs(mt12-mt1).GT.1) cycle
        isofakt = CG(t12,1,tau3,mt12,mtau3-mt12,mtau3) &
        & *CG(1,1,t12,mt1,mt12-mt1,mt12) &
        & *CG(t12p,1,tau3p,mt12,mtau3-mt12,mtau3) &
        & *CG(1,1,t12p,mt1,mt12-mt1,mt12)
        
        isomass(1,alphap,alphaloc) = isomass(1,alphap,alphaloc) &
        & + isofakt &
        &   *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12) &
        &    - mthres(mtau3))
        isomass(2,alphap,alphaloc) = isomass(2,alphap,alphaloc) &
        & + isofakt &
        &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  & 
        &      /(2.0*mbaryon(1,mt1)*mbaryon(1,mt12-mt1))
        isomass(3,alphap,alphaloc) = isomass(3,alphap,alphaloc) &
        & + isofakt &
        &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12))  & 
        &      /(2.0*(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  &
                    *mbaryon(1,mtau3-mt12))
       END DO ! mt1 
      END DO  ! mt12   
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   locekin=0.0_dpreal 
   DO alphap=1,alpha3Ncdepmax
    DO alphaloc=1,mynalpha3N
     IF(isomass(2,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
      DO ip3=1,mynp3
       DO ip12=1,mynp12
        locekin=locekin &
       & + conjg(psiamp3N(ip12,ip3,alphap)) &
            *PSI%amp(ip12,ip3,alphaloc) &
       &    *p12p3weight(ip12,ip3) &
       &    *(isomass(1,alphap,alphaloc) &
       &    + isomass(2,alphap,alphaloc)*P12P(myp12+ip12-1)**2 & 
       &    + isomass(3,alphap,alphaloc)*P3P(myp3+ip3-1)**2)
       
       END DO ! ip12 
      END DO ! ip3 
     END IF
    END DO ! alphaloc
   END DO ! alphap

   locnorm=0.0_dpreal 
   DO alphaloc=1,mynalpha3N
    DO ip3=1,mynp3
     DO ip12=1,mynp12
      locnorm=locnorm &
     & + conjg(PSIFAD%amp(ip12,ip3,alphaloc)) &
          *PSI%amp(ip12,ip3,alphaloc) &
     &    *p12p3weight(ip12,ip3) 
     END DO ! ip12 
    END DO ! ip3 
   END DO ! alphaloc
   
   DEALLOCATE(psiamp3N,isomass) 
   
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locekin,ekin,1,MPI_REAL4,MPI_SUM,commamp3n,ierr) ! here it is assumed that product is REAL(spreal)
   CALL MPI_ALLREDUCE(locnorm,norm,1,MPI_REAL4,MPI_SUM,commamp3n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A,E15.6)') &
    &  'kinetic energy from 3N wave and fadcomp [MeV]: ',3.0*ekin*hbarc
    WRITE(*,'(A,E15.6)') &
    &  'Norm of 3N wave and fadcomp amplitude:         ',3.0*norm
   END IF
   
  END IF   ! 3N wave 
  
  DEALLOCATE(mbaryon,mthres) 
  
 END SUBROUTINE printekin_3n
 
 SUBROUTINE printls_3n(PSI,PSIFAD)
  IMPLICIT NONE
  TYPE(AMP3N) PSI,PSIFAD
  REAL(spreal),ALLOCATABLE :: locpls(:,:),pls(:,:)
  INTEGER ip12,ip3,alphaloc,alpha,ierr
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  INTEGER alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p, &
 &        mtau3p,alpha2bp,parip,mt12,mt1,bl,bs
  COMPLEX(spreal),ALLOCATABLE :: psiamp3N(:,:,:)
  REAL(spreal),ALLOCATABLE :: recoupl(:,:,:,:) 
  REAL(spreal) :: norm 
  

  
  IF(allocated(PSI%amp) .AND. allocated(PSIFAD%amp)) THEN ! calculate based on 3N wave and fadcomp
! start with kinetic energy
! since there is a possibility to have isospin transitions
! I need to collect all channels 
   ALLOCATE(psiamp3N(mynp12,mynp3,alpha3Ncdepmax))
   CALL collect_piece_spcmplx(PSIFAD%amp,psiamp3N,&
          &     mynp12*mynp3,alpha3Ncdepmax,1,commalpha)

! prepare the recoupl coeffs for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    recoupl(L,S,alphap,alpha)
!       L total angular momentum 0:j3_max+3/2 
!       S total spin 0  =  1/2  1 = 3/2           
          
          
   ALLOCATE(recoupl(0:(j3max+3)/2,0:1,alpha3Ncdepmax,mynalpha3N))
   
   recoupl=0.0_spreal
   
   DO alphap=1,alpha3Ncdepmax
    call get3Nqn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,&
       &          j3p,tau3p,mtau3p,alpha2bp,parip)
    DO alphaloc=1,mynalpha3N
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)
     IF(t12.eq.t12p .and. &
      & tau3.eq.tau3p .and.  mtau3.eq.mtau3p .and. &
      & l12.eq.l12p .and. l3.eq.l3p .and. &
      & s12.eq.s12p .and. j3.eq.j3p) THEN
      DO bs=1,3,2 
       DO bl=abs(j3-bs)/2,(j3+bs)/2        
        recoupl(bl,bs/2,alphap,alphaloc) = & 
        &  SDCH(bl)*SDCH(j12)*SD2CH(I3)*SD2CH(bs)  &
        &   *C9J(2*l12,2*s12,2*j12,2*l3,1,I3,2*bl,bs,j3) &
        &   *SDCH(bl)*SDCH(j12p)*SD2CH(I3p)*SD2CH(bs)  &
        &   *C9J(2*l12,2*s12,2*j12p,2*l3,1,I3p,2*bl,bs,j3)
       END DO ! bs
      END DO  ! bl   
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   ALLOCATE(locpls(0:(j3max+3)/2,0:1),pls(0:(j3max+3)/2,0:1))
   
   locpls=0.0_dpreal 
   DO bl=0,(j3max+3)/2
    DO bs=1,3,2 
     DO alphap=1,alpha3Ncdepmax
      DO alphaloc=1,mynalpha3N
       IF(recoupl(bl,bs/2,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
        DO ip3=1,mynp3
         DO ip12=1,mynp12
          locpls(bl,bs/2)=locpls(bl,bs/2) &
       & + conjg(psiamp3N(ip12,ip3,alphap)) &
            *PSI%amp(ip12,ip3,alphaloc) &
       &    *p12p3weight(ip12,ip3) &
       &    *recoupl(bl,bs/2,alphap,alphaloc)
         END DO ! ip12 
        END DO ! ip3 
       END IF
      END DO ! alphaloc
     END DO ! alphap
    END DO ! bs 
   END DO ! bl 
   
   
   
   DEALLOCATE(psiamp3N,recoupl) 
   
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locpls,pls,j3max+5,MPI_REAL4,MPI_SUM,commamp3n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A)') 'LS probabilities from 3N wave: '
    norm=0.0_spreal
    DO bl=0,(j3max+3)/2
     DO bs=1,3,2 
      IF(pls(bl,bs/2).NE.0.0_spreal) THEN
       WRITE(*,'(A,2I4,2X,E15.6)') 'PLS-3N(WAVE-FAD): ',bl,bs,3.0*pls(bl,bs/2)
       norm=norm+3.0*pls(bl,bs/2)
      END IF
     END DO
    END DO
    WRITE(*,'(A,8X,2X,E15.6)') 'PLS-3N(WAVE-FAD): ',norm
   END IF
   
   DEALLOCATE(locpls,pls)
   
  END IF   ! 3N ampl  
  
 END SUBROUTINE printls_3n

 
 SUBROUTINE printinfo_3n(PSI,PSIFAD)
  IMPLICIT NONE
  TYPE(AMP3N) PSI,PSIFAD 
  
  CALL printekin_3n(PSI,PSIFAD)
  CALL printls_3n(PSI,PSIFAD)
  
 END SUBROUTINE printinfo_3n

 ! first case: 4N wave function 3+1,2+2  & yak components 
 
 SUBROUTINE printekin_4n(PSI4N31,PSI4N22,PHI4N31,PHI4N22)
  IMPLICIT NONE
  TYPE(AMP4N31) PSI4N31,PHI4N31 
  TYPE(AMP4N22) PSI4N22,PHI4N22
  
  REAL(spreal) :: locekin,ekin,locnorm,norm
  INTEGER ip12,ip3,iq4,ip34,iq,alphaloc,alpha,alphap,ierr
  INTEGER l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p,alpha2bp,alpha3bp,parip
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari,mtau3,mt12,mt1,mt34,mt3
  INTEGER l34p,s34p,j34p,t34p,lamp,Ip,l34,s34,j34,t34,lam,I,alpha12,alpha34,alpha12p,alpha34p
  COMPLEX(spreal),ALLOCATABLE :: psiamp4N(:,:,:,:)
  REAL(spreal),ALLOCATABLE :: isomass(:,:,:)
  REAL(spreal),ALLOCATABLE :: mbaryon(:,:),mthres(:)
  REAL(spreal) :: isofakt 
  
!  prepare masses depending on isospin 
  ALLOCATE(mbaryon(0:2,-2:2),mthres(-4:4))
  mbaryon=0.0_spreal
! set all baryons, only nucleons are used   
  mbaryon(0,0)=mlam
  mbaryon(1,-1)=mneu
  mbaryon(1,1)=mprot
  mbaryon(2,-2)=msigm
  mbaryon(2,0)=msig0
  mbaryon(2,2)=msigp
  mthres=0.0_spreal
  mthres(4)=4.0*mprot
  mthres(2)=3.0*mprot+mneu
  mthres(0)=2.0*mneu+2.0*mprot
  mthres(-2)=3.0*mneu+mprot
  mthres(-4)=4.0*mneu
  
  IF(allocated(PSI4N31%amp).AND. allocated(PSI4N22%amp) .AND. &
    & allocated(PHI4N31%amp).AND. allocated(PHI4N22%amp)  ) THEN ! calculate based on 4N wave and Yakubovsky components  
! start with kinetic energy part in 3+1 coordinates 
! since there is a possibility to have isospin transitions
! I need to collect all channels
    
! start with 3+1 part 
    
   ALLOCATE(psiamp4N(mynp12,mynp3,mynq4,alpha4N31cdepmax))
   CALL collect_piece_spcmplx(PHI4N31%amp,psiamp4N,&
          &     mynp12*mynp3*mynq4,alpha4N31cdepmax,1,commalpha)

! prepare the isospin masses for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    isomass(i,alphap,alpha)
!       i=1   total mass 
!       i=2   reduced mass average particles 12 
!       i=3   reduced mass average particle 3
!       i=4   reduced mass average particle 4           
          
          
   ALLOCATE(isomass(1:4,alpha4N31cdepmax,mynalpha4N31)) 
   isomass=0.0_spreal
   
   DO alphap=1,alpha4N31cdepmax
    call get4N31qn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p,alpha2bp,alpha3bp,parip)
    t12p=2*t12p   !   use convention that isospins are twice 
    tau4p=2*tau4p
    mtau4p=2*mtau4p
    DO alphaloc=1,mynalpha4N31
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     call get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari)
     t12=2*t12   !   use convention that isospins are twice
     tau4=2*tau4
     mtau4=2*mtau4
     IF(l12.eq.l12p .and. s12.eq.s12p .and.&
      & j12.eq.j12p .and.  l3.eq.l3p  .and.&
      &  I3.eq.I3p  .and.  j3.eq.j3p  .and.&
      &  l4.eq.l4p  .and.  I4.eq.I4p  .and.&
      &  j4.eq.j4p  .and.  mtau4.eq.mtau4p) THEN
      DO mtau3=-tau3,tau3,2
       IF(abs(mtau4-mtau3).GT.1) cycle
       DO mt12=-t12,t12,2
        IF(abs(mtau3-mt12).GT.1) cycle
        DO mt1=-1,1,2 
         IF(abs(mt12-mt1).GT.1) cycle
         isofakt = CG(tau3,1,tau4,mtau3,mtau4-mtau3,mtau4) &
         & *CG(t12,1,tau3,mt12,mtau3-mt12,mtau3) &
         & *CG(1,1,t12,mt1,mt12-mt1,mt12) &
         & *CG(tau3p,1,tau4p,mtau3,mtau4-mtau3,mtau4) &
         & *CG(t12p,1,tau3p,mt12,mtau3-mt12,mtau3) &
         & *CG(1,1,t12p,mt1,mt12-mt1,mt12)

         isomass(1,alphap,alphaloc) = isomass(1,alphap,alphaloc) &
         & + isofakt &
         &   *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12)+mbaryon(1,mtau4-mtau3) &
         &    - mthres(mtau4))
         isomass(2,alphap,alphaloc) = isomass(2,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  & 
         &      /(2.0*mbaryon(1,mt1)*mbaryon(1,mt12-mt1))
         isomass(3,alphap,alphaloc) = isomass(3,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12))  & 
         &      /(2.0*(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  &
                     *mbaryon(1,mtau3-mt12))
         isomass(4,alphap,alphaloc) = isomass(4,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12)+mbaryon(1,mtau4-mtau3))  & 
         &      /(2.0*(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12))  &
                     *mbaryon(1,mtau4-mtau3))            
        END DO ! mt1 
       END DO  ! mt12 
      END DO   ! mtau3
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   locekin=0.0_dpreal 
   DO alphap=1,alpha4N31cdepmax
    DO alphaloc=1,mynalpha4N31
     IF(isomass(2,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
      DO iq4=1,mynq4
       DO ip3=1,mynp3
        DO ip12=1,mynp12
         locekin=locekin &
        & + 12.0*conjg(psiamp4N(ip12,ip3,iq4,alphap)) &
             *PSI4N31%amp(ip12,ip3,iq4,alphaloc) &
        &    *p12p3q4weight(ip12,ip3,iq4) &
        &    *(isomass(1,alphap,alphaloc) &
        &    + isomass(2,alphap,alphaloc)*P12P(myp12+ip12-1)**2 & 
        &    + isomass(3,alphap,alphaloc)*P3P(myp3+ip3-1)**2  &
        &    + isomass(4,alphap,alphaloc)*Q4P(myq4+iq4-1)**2)

        END DO ! ip12 
       END DO ! ip3
      END DO  ! iq4 
     END IF
    END DO ! alphaloc
   END DO ! alphap

   locnorm=0.0_dpreal 
   DO alphaloc=1,mynalpha4N31
    DO iq4=1,mynq4
     DO ip3=1,mynp3
      DO ip12=1,mynp12
       locnorm=locnorm &
     & + 12.0*conjg(PHI4N31%amp(ip12,ip3,iq4,alphaloc)) &
          *PSI4N31%amp(ip12,ip3,iq4,alphaloc) &
     &    *p12p3q4weight(ip12,ip3,iq4)   
      END DO ! ip12 
     END DO ! ip3  
    END DO ! iq4
   END DO ! alphaloc
   
   DEALLOCATE(psiamp4N,isomass) 
   
   
   
! start with 2+2 part 
    
   ALLOCATE(psiamp4N(mynp12,mynp34,mynq,beta4N22cdepmax))
   CALL collect_piece_spcmplx(PHI4N22%amp,psiamp4N,&
          &     mynp12*mynp34*mynq,beta4N22cdepmax,1,commalpha)

! prepare the isospin masses for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    isomass(i,alphap,alpha)
!       i=1   total mass 
!       i=2   reduced mass average particles 12 
!       i=3   reduced mass average particle 34
!       i=4   reduced mass average particle  12+34           
          
          
   ALLOCATE(isomass(1:4,beta4N22cdepmax,mynbeta4N22)) 
   isomass=0.0_spreal
   
   DO alphap=1,beta4N22cdepmax
    call get4N22qn(alphap,l12p,s12p,j12p,t12p,l34p,s34p,j34p,t34p,lamp,Ip,j4p,tau4p,mtau4p,alpha12p,alpha34p,parip)
    t12p=2*t12p   !   use convention that isospins are twice 
    t34p=2*t34p 
    tau4p=2*tau4p
    mtau4p=2*mtau4p
    DO alphaloc=1,mynbeta4N22
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     call get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari)
     t12=2*t12   !   use convention that isospins are twice
     t34=2*t34 
     tau4=2*tau4
     mtau4=2*mtau4
     IF(l12.eq.l12p .and. s12.eq.s12p .and.&
      & j12.eq.j12p .and.  l34.eq.l34p  .and.&
      &  s34.eq.s34p  .and.  j34.eq.j34p  .and.&
      &  lam.eq.lamp  .and.  I.eq.Ip  .and.&
      &  j4.eq.j4p  .and.  mtau4.eq.mtau4p) THEN
      DO mt12=-t12,t12,2
       IF(abs(mtau4-mt12).GT.t34) cycle
       DO mt3=-1,1,2
        IF(abs(mtau4-mt12-mt3).GT.1) cycle
        DO mt1=-1,1,2 
         IF(abs(mt12-mt1).GT.1) cycle
         isofakt = CG(t12,t34,tau4,mt12,mtau4-mt12,mtau4) &
         & *CG(1,1,t34,mt3,mtau4-mt12-mt3,mtau4-mt12) &
         & *CG(1,1,t12,mt1,mt12-mt1,mt12) &
         & *CG(t12p,t34p,tau4p,mt12,mtau4-mt12,mtau4)  &
         & *CG(1,1,t34p,mt3,mtau4-mt12-mt3,mtau4-mt12) &
         & *CG(1,1,t12p,mt1,mt12-mt1,mt12)

         isomass(1,alphap,alphaloc) = isomass(1,alphap,alphaloc) &
         & + isofakt &
         &   *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mt3)+mbaryon(1,mtau4-mt12-mt3) &
         &    - mthres(mtau4))
         isomass(2,alphap,alphaloc) = isomass(2,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  & 
         &      /(2.0*mbaryon(1,mt1)*mbaryon(1,mt12-mt1))
         isomass(3,alphap,alphaloc) = isomass(3,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt3)+mbaryon(1,mtau4-mt12-mt3))  & 
         &      /(2.0*(mbaryon(1,mt3)*mbaryon(1,mtau4-mt12-mt3)))
         isomass(4,alphap,alphaloc) = isomass(4,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mt3)+mbaryon(1,mtau4-mt12-mt3))  & 
         &      /(2.0*(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))   &
         &            *(mbaryon(1,mt3)+mbaryon(1,mtau4-mt12-mt3)))            
        END DO ! mt1 
       END DO  ! mt3 
      END DO   ! mt12
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   ! do not set locekin to zero since the 3+1 part needs to be included 
   DO alphap=1,beta4N22cdepmax
    DO alphaloc=1,mynbeta4N22
     IF(isomass(2,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
      DO iq=1,mynq
       DO ip34=1,mynp34
        DO ip12=1,mynp12
         locekin=locekin &
        & + 6.0*conjg(psiamp4N(ip12,ip34,iq,alphap)) &
        &    *PSI4N22%amp(ip12,ip34,iq,alphaloc) &
        &    *p12p34qweight(ip12,ip34,iq) &
        &    *(isomass(1,alphap,alphaloc) &
        &    + isomass(2,alphap,alphaloc)*P12P(myp12+ip12-1)**2 & 
        &    + isomass(3,alphap,alphaloc)*P34P(myp34+ip34-1)**2  &
        &    + isomass(4,alphap,alphaloc)*QP(myq+iq-1)**2)

        END DO ! ip12 
       END DO ! ip3
      END DO  ! iq4 
     END IF
    END DO ! alphaloc
   END DO ! alphap

   ! do not set locnorm to zero since 3+1 part needs to be kept 
   DO alphaloc=1,mynbeta4N22
    DO iq=1,mynq
     DO ip34=1,mynp34
      DO ip12=1,mynp12
       locnorm=locnorm &
     & + 6.0*conjg(PHI4N22%amp(ip12,ip34,iq,alphaloc)) &
          *PSI4N22%amp(ip12,ip34,iq,alphaloc) &
     &    *p12p34qweight(ip12,ip34,iq)   
      END DO ! ip12 
     END DO ! ip3  
    END DO ! iq4
   END DO ! alphaloc
   
   DEALLOCATE(psiamp4N,isomass) 
   
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locekin,ekin,1,MPI_REAL4,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is REAL(spreal)
   CALL MPI_ALLREDUCE(locnorm,norm,1,MPI_REAL4,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A,E15.6)') &
    &  'kinetic energy from 4N wave and Yakubovsky components [MeV]: ',ekin*hbarc
    WRITE(*,'(A,E15.6)') &
    &  'Norm of 4N wave and Yakubovsky components:         ',norm
   END IF
   
  END IF   ! 4N wave 
  
  DEALLOCATE(mbaryon,mthres) 
  
 END SUBROUTINE printekin_4n
 
 
  SUBROUTINE printekin_4n31(PSI4N31)
  IMPLICIT NONE
  TYPE(AMP4N31) PSI4N31 
  
  REAL(spreal) :: locekin,ekin,locnorm,norm
  INTEGER ip12,ip3,iq4,ip34,iq,alphaloc,alpha,alphap,ierr
  INTEGER l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p,alpha2bp,alpha3bp,parip
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari,mtau3,mt12,mt1,mt34,mt3
  INTEGER alpha12,alpha34,alpha12p,alpha34p
  COMPLEX(spreal),ALLOCATABLE :: psiamp4N(:,:,:,:)
  REAL(spreal),ALLOCATABLE :: isomass(:,:,:)
  REAL(spreal),ALLOCATABLE :: mbaryon(:,:),mthres(:)
  REAL(spreal) :: isofakt 
  
!  prepare masses depending on isospin 
  ALLOCATE(mbaryon(0:2,-2:2),mthres(-4:4))
  mbaryon=0.0_spreal
! set all baryons, only nucleons are used   
  mbaryon(0,0)=mlam
  mbaryon(1,-1)=mneu
  mbaryon(1,1)=mprot
  mbaryon(2,-2)=msigm
  mbaryon(2,0)=msig0
  mbaryon(2,2)=msigp
  mthres=0.0_spreal
  mthres(4)=4.0*mprot
  mthres(2)=3.0*mprot+mneu
  mthres(0)=2.0*mneu+2.0*mprot
  mthres(-2)=3.0*mneu+mprot
  mthres(-4)=4.0*mneu
  
  IF(allocated(PSI4N31%amp)) THEN ! calculate based on 4N31 wave 
! start with kinetic energy part in 3+1 coordinates 
! since there is a possibility to have isospin transitions
! I need to collect all channels
    
! start with 3+1 part 
    
   ALLOCATE(psiamp4N(mynp12,mynp3,mynq4,alpha4N31cdepmax))
   CALL collect_piece_spcmplx(PSI4N31%amp,psiamp4N,&
          &     mynp12*mynp3*mynq4,alpha4N31cdepmax,1,commalpha)

! prepare the isospin masses for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    isomass(i,alphap,alpha)
!       i=1   total mass 
!       i=2   reduced mass average particles 12 
!       i=3   reduced mass average particle 3
!       i=4   reduced mass average particle 4           
          
          
   ALLOCATE(isomass(1:4,alpha4N31cdepmax,mynalpha4N31)) 
   isomass=0.0_spreal
   
   DO alphap=1,alpha4N31cdepmax
    call get4N31qn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p,alpha2bp,alpha3bp,parip)
    t12p=2*t12p   !   use convention that isospins are twice 
    tau4p=2*tau4p
    mtau4p=2*mtau4p
    DO alphaloc=1,mynalpha4N31
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     call get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari)
     t12=2*t12   !   use convention that isospins are twice
     tau4=2*tau4
     mtau4=2*mtau4
     IF(l12.eq.l12p .and. s12.eq.s12p .and.&
      & j12.eq.j12p .and.  l3.eq.l3p  .and.&
      &  I3.eq.I3p  .and.  j3.eq.j3p  .and.&
      &  l4.eq.l4p  .and.  I4.eq.I4p  .and.&
      &  j4.eq.j4p  .and.  mtau4.eq.mtau4p) THEN
      DO mtau3=-tau3,tau3,2
       IF(abs(mtau4-mtau3).GT.1) cycle
       DO mt12=-t12,t12,2
        IF(abs(mtau3-mt12).GT.1) cycle
        DO mt1=-1,1,2 
         IF(abs(mt12-mt1).GT.1) cycle
         isofakt = CG(tau3,1,tau4,mtau3,mtau4-mtau3,mtau4) &
         & *CG(t12,1,tau3,mt12,mtau3-mt12,mtau3) &
         & *CG(1,1,t12,mt1,mt12-mt1,mt12) &
         & *CG(tau3p,1,tau4p,mtau3,mtau4-mtau3,mtau4) &
         & *CG(t12p,1,tau3p,mt12,mtau3-mt12,mtau3) &
         & *CG(1,1,t12p,mt1,mt12-mt1,mt12)

         isomass(1,alphap,alphaloc) = isomass(1,alphap,alphaloc) &
         & + isofakt &
         &   *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12)+mbaryon(1,mtau4-mtau3) &
         &    - mthres(mtau4))
         isomass(2,alphap,alphaloc) = isomass(2,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  & 
         &      /(2.0*mbaryon(1,mt1)*mbaryon(1,mt12-mt1))
         isomass(3,alphap,alphaloc) = isomass(3,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12))  & 
         &      /(2.0*(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  &
                     *mbaryon(1,mtau3-mt12))
         isomass(4,alphap,alphaloc) = isomass(4,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12)+mbaryon(1,mtau4-mtau3))  & 
         &      /(2.0*(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mtau3-mt12))  &
                     *mbaryon(1,mtau4-mtau3))            
        END DO ! mt1 
       END DO  ! mt12 
      END DO   ! mtau3
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   locekin=0.0_dpreal 
   DO alphap=1,alpha4N31cdepmax
    DO alphaloc=1,mynalpha4N31
     IF(isomass(2,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
      DO iq4=1,mynq4
       DO ip3=1,mynp3
        DO ip12=1,mynp12
         locekin=locekin &
        & + conjg(psiamp4N(ip12,ip3,iq4,alphap)) &
             *PSI4N31%amp(ip12,ip3,iq4,alphaloc) &
        &    *p12p3q4weight(ip12,ip3,iq4) &
        &    *(isomass(1,alphap,alphaloc) &
        &    + isomass(2,alphap,alphaloc)*P12P(myp12+ip12-1)**2 & 
        &    + isomass(3,alphap,alphaloc)*P3P(myp3+ip3-1)**2  &
        &    + isomass(4,alphap,alphaloc)*Q4P(myq4+iq4-1)**2)

        END DO ! ip12 
       END DO ! ip3
      END DO  ! iq4 
     END IF
    END DO ! alphaloc
   END DO ! alphap

   locnorm=0.0_dpreal 
   DO alphaloc=1,mynalpha4N31
    DO iq4=1,mynq4
     DO ip3=1,mynp3
      DO ip12=1,mynp12
       locnorm=locnorm &
     & + conjg(PSI4N31%amp(ip12,ip3,iq4,alphaloc)) &
          *PSI4N31%amp(ip12,ip3,iq4,alphaloc) &
     &    *p12p3q4weight(ip12,ip3,iq4)   
      END DO ! ip12 
     END DO ! ip3  
    END DO ! iq4
   END DO ! alphaloc
   
   DEALLOCATE(psiamp4N,isomass) 
   
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locekin,ekin,1,MPI_REAL4,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is REAL(spreal)
   CALL MPI_ALLREDUCE(locnorm,norm,1,MPI_REAL4,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A,E15.6)') &
    &  'kinetic energy from 4N31 wave [MeV]: ',ekin*hbarc
    WRITE(*,'(A,E15.6)') &
    &  'Norm of 4N31 wave:         ',norm
   END IF
   
  END IF   ! amp4N31   
  
  DEALLOCATE(mbaryon,mthres) 
  
 END SUBROUTINE printekin_4n31
 
  
 SUBROUTINE printekin_4n22(PSI4N22)
  IMPLICIT NONE
  TYPE(AMP4N22) PSI4N22
  
  REAL(spreal) :: locekin,ekin,locnorm,norm
  INTEGER ip12,ip3,iq4,ip34,iq,alphaloc,alpha,alphap,ierr
  INTEGER l12p,s12p,j12p,t12p,j4p,tau4p,mtau4p,alpha2bp,parip
  INTEGER l12,s12,j12,t12,j4,tau4,mtau4,alpha2b,pari,mt12,mt1,mt34,mt3
  INTEGER l34p,s34p,j34p,t34p,lamp,Ip,l34,s34,j34,t34,lam,I,alpha12,alpha34,alpha12p,alpha34p
  COMPLEX(spreal),ALLOCATABLE :: psiamp4N(:,:,:,:)
  REAL(spreal),ALLOCATABLE :: isomass(:,:,:)
  REAL(spreal),ALLOCATABLE :: mbaryon(:,:),mthres(:)
  REAL(spreal) :: isofakt 
  
!  prepare masses depending on isospin 
  ALLOCATE(mbaryon(0:2,-2:2),mthres(-4:4))
  mbaryon=0.0_spreal
! set all baryons, only nucleons are used   
  mbaryon(0,0)=mlam
  mbaryon(1,-1)=mneu
  mbaryon(1,1)=mprot
  mbaryon(2,-2)=msigm
  mbaryon(2,0)=msig0
  mbaryon(2,2)=msigp
  mthres=0.0_spreal
  mthres(4)=4.0*mprot
  mthres(2)=3.0*mprot+mneu
  mthres(0)=2.0*mneu+2.0*mprot
  mthres(-2)=3.0*mneu+mprot
  mthres(-4)=4.0*mneu
  
  
  IF(allocated(PSI4N22%amp)) THEN ! calculate based on 4N22 wave 
! start with 2+2 part 
    
   ALLOCATE(psiamp4N(mynp12,mynp34,mynq,beta4N22cdepmax))
   CALL collect_piece_spcmplx(PSI4N22%amp,psiamp4N,&
          &     mynp12*mynp34*mynq,beta4N22cdepmax,1,commalpha)

! prepare the isospin masses for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    isomass(i,alphap,alpha)
!       i=1   total mass 
!       i=2   reduced mass average particles 12 
!       i=3   reduced mass average particle 34
!       i=4   reduced mass average particle  12+34           
          
          
   ALLOCATE(isomass(1:4,beta4N22cdepmax,mynbeta4N22)) 
   isomass=0.0_spreal
   
   DO alphap=1,beta4N22cdepmax
    call get4N22qn(alphap,l12p,s12p,j12p,t12p,l34p,s34p,j34p,t34p,lamp,Ip,j4p,tau4p,mtau4p,alpha12p,alpha34p,parip)
    t12p=2*t12p   !   use convention that isospins are twice 
    t34p=2*t34p 
    tau4p=2*tau4p
    mtau4p=2*mtau4p
    DO alphaloc=1,mynbeta4N22
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     call get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari)
     t12=2*t12   !   use convention that isospins are twice
     t34=2*t34 
     tau4=2*tau4
     mtau4=2*mtau4
     IF(l12.eq.l12p .and. s12.eq.s12p .and.&
      & j12.eq.j12p .and.  l34.eq.l34p  .and.&
      &  s34.eq.s34p  .and.  j34.eq.j34p  .and.&
      &  lam.eq.lamp  .and.  I.eq.Ip  .and.&
      &  j4.eq.j4p  .and.  mtau4.eq.mtau4p) THEN
      DO mt12=-t12,t12,2
       IF(abs(mtau4-mt12).GT.t34) cycle
       DO mt3=-1,1,2
        IF(abs(mtau4-mt12-mt3).GT.1) cycle
        DO mt1=-1,1,2 
         IF(abs(mt12-mt1).GT.1) cycle
         isofakt = CG(t12,t34,tau4,mt12,mtau4-mt12,mtau4) &
         & *CG(1,1,t34,mt3,mtau4-mt12-mt3,mtau4-mt12) &
         & *CG(1,1,t12,mt1,mt12-mt1,mt12) &
         & *CG(t12p,t34p,tau4p,mt12,mtau4-mt12,mtau4)  &
         & *CG(1,1,t34p,mt3,mtau4-mt12-mt3,mtau4-mt12) &
         & *CG(1,1,t12p,mt1,mt12-mt1,mt12)

         isomass(1,alphap,alphaloc) = isomass(1,alphap,alphaloc) &
         & + isofakt &
         &   *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mt3)+mbaryon(1,mtau4-mt12-mt3) &
         &    - mthres(mtau4))
         isomass(2,alphap,alphaloc) = isomass(2,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  & 
         &      /(2.0*mbaryon(1,mt1)*mbaryon(1,mt12-mt1))
         isomass(3,alphap,alphaloc) = isomass(3,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt3)+mbaryon(1,mtau4-mt12-mt3))  & 
         &      /(2.0*(mbaryon(1,mt3)*mbaryon(1,mtau4-mt12-mt3)))
         isomass(4,alphap,alphaloc) = isomass(4,alphap,alphaloc) &
         & + isofakt &
         &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1)+mbaryon(1,mt3)+mbaryon(1,mtau4-mt12-mt3))  & 
         &      /(2.0*(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))   &
         &            *(mbaryon(1,mt3)+mbaryon(1,mtau4-mt12-mt3)))            
        END DO ! mt1 
       END DO  ! mt3 
      END DO   ! mt12
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   locekin=0.0_dpreal 
   DO alphap=1,beta4N22cdepmax
    DO alphaloc=1,mynbeta4N22
     IF(isomass(2,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
      DO iq=1,mynq
       DO ip34=1,mynp34
        DO ip12=1,mynp12
         locekin=locekin &
        & + conjg(psiamp4N(ip12,ip34,iq,alphap)) &
        &    *PSI4N22%amp(ip12,ip34,iq,alphaloc) &
        &    *p12p34qweight(ip12,ip34,iq) &
        &    *(isomass(1,alphap,alphaloc) &
        &    + isomass(2,alphap,alphaloc)*P12P(myp12+ip12-1)**2 & 
        &    + isomass(3,alphap,alphaloc)*P34P(myp34+ip34-1)**2  &
        &    + isomass(4,alphap,alphaloc)*QP(myq+iq-1)**2)

        END DO ! ip12 
       END DO ! ip3
      END DO  ! iq4 
     END IF
    END DO ! alphaloc
   END DO ! alphap

   locnorm=0.0_dpreal 
   DO alphaloc=1,mynbeta4N22
    DO iq=1,mynq
     DO ip34=1,mynp34
      DO ip12=1,mynp12
       locnorm=locnorm &
     & + conjg(PSI4N22%amp(ip12,ip34,iq,alphaloc)) &
          *PSI4N22%amp(ip12,ip34,iq,alphaloc) &
     &    *p12p34qweight(ip12,ip34,iq)   
      END DO ! ip12 
     END DO ! ip3  
    END DO ! iq4
   END DO ! alphaloc
   
   DEALLOCATE(psiamp4N,isomass) 
   
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locekin,ekin,1,MPI_REAL4,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is REAL(spreal)
   CALL MPI_ALLREDUCE(locnorm,norm,1,MPI_REAL4,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A,E15.6)') &
    &  'kinetic energy from 4N22 wave [MeV]: ',ekin*hbarc
    WRITE(*,'(A,E15.6)') &
    &  'Norm of 4N22 wave:         ',norm
   END IF
   
  END IF   ! 4N22 wave 
  
  DEALLOCATE(mbaryon,mthres) 
  
 END SUBROUTINE printekin_4n22
 
 
 
 SUBROUTINE printls_4n(PSI4N31,PSI4N22,PHI4N31,PHI4N22)
  IMPLICIT NONE
  TYPE(AMP4N31) PSI4N31,PHI4N31 
  TYPE(AMP4N22) PSI4N22,PHI4N22
  
  REAL(spreal),ALLOCATABLE :: locpls(:,:),pls(:,:)
  INTEGER ip12,ip3,iq4,ip34,iq,alphaloc,alpha,alphap,ierr
  INTEGER l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p,alpha2bp,alpha3bp,parip
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari,mtau3,mt12,mt1,mt34,mt3
  INTEGER l34p,s34p,j34p,t34p,lamp,Ip,l34,s34,j34,t34,lam,I,alpha12,alpha34,alpha12p,alpha34p,bl,bs,bl3,bs3
  COMPLEX(spreal),ALLOCATABLE :: psiamp4N(:,:,:,:)
  REAL(spreal),ALLOCATABLE :: recoupl(:,:,:,:) 
  REAL(spreal) :: norm 
  

  
  IF(allocated(PSI4N31%amp).AND. allocated(PSI4N22%amp) .AND. &
    & allocated(PHI4N31%amp).AND. allocated(PHI4N22%amp)  ) THEN ! calculate based on 4N wave and Yakubovsky components
! start with kinetic energy
! since there is a possibility to have isospin transitions
! I need to collect all channels 
! start with 3+1 part 
    
   ALLOCATE(psiamp4N(mynp12,mynp3,mynq4,alpha4N31cdepmax))
   CALL collect_piece_spcmplx(PHI4N31%amp,psiamp4N,&
          &     mynp12*mynp3*mynq4,alpha4N31cdepmax,1,commalpha)

! prepare the recoupl coeffs for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    recoupl(L,S,alphap,alpha)
!       L total angular momentum 0:j4_max+2 
!       S total spin 0 : 2        
          
          
   ALLOCATE(recoupl(0:(j4max+2),0:2,alpha4N31cdepmax,mynalpha4N31))
   
   recoupl=0.0_spreal
   
   DO alphap=1,alpha4N31cdepmax
    call get4N31qn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p,alpha2bp,alpha3bp,parip)
    DO alphaloc=1,mynalpha4N31
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari) 
     IF(t12.eq.t12p .and. &
      & tau3.eq.tau3p .and. tau4.eq.tau4p .and. mtau4.eq.mtau4p .and. &
      & l12.eq.l12p .and. l3.eq.l3p .and. l4.eq.l4p .and. &
      & s12.eq.s12p .and. j4.eq.j4p) THEN
      DO bs=0,2 
       DO bl=abs(j4-bs),j4+bs 
        DO bl3=abs(bl-l4),bl+l4
         DO bs3=abs(j3-2*bl3),j3+2*bl3,2  
          recoupl(bl,bs,alphap,alphaloc) = recoupl(bl,bs,alphap,alphaloc) + & 
        &  SDCH(bl3)*SDCH(j12)*SD2CH(I3)*SD2CH(bs3)  &
        &   *C9J(2*l12,2*s12,2*j12,2*l3,1,I3,2*bl3,bs3,j3) &
        &   *SD2CH(j3)*SD2CH(I4)*SDCH(bl)*SDCH(bs)  &
        &   *C9J(2*bl3,bs3,j3,2*l4,1,I4,2*bl,2*bs,2*j4) &        
        &   *SDCH(bl3)*SDCH(j12p)*SD2CH(I3p)*SD2CH(bs3)  &
        &   *C9J(2*l12,2*s12,2*j12p,2*l3,1,I3p,2*bl3,bs3,j3p) &
        &   *SD2CH(j3p)*SD2CH(I4p)*SDCH(bl)*SDCH(bs)  &
        &   *C9J(2*bl3,bs3,j3p,2*l4,1,I4p,2*bl,2*bs,2*j4)
         END DO ! bs3
        END DO  ! bl3 
       END DO ! bs
      END DO  ! bl   
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   ALLOCATE(locpls(0:(j4max+2),0:2),pls(0:(j4max+2),0:2))
   
   locpls=0.0_dpreal 
   DO bl=0,j4max+2
    DO bs=0,2 
     DO alphap=1,alpha4N31cdepmax
      DO alphaloc=1,mynalpha4N31
       IF(recoupl(bl,bs,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
        DO iq4=1,mynq4
         DO ip3=1,mynp3
          DO ip12=1,mynp12
           locpls(bl,bs)=locpls(bl,bs) &
       & + 12.0*conjg(psiamp4N(ip12,ip3,iq4,alphap)) &
       &    *PSI4N31%amp(ip12,ip3,iq4,alphaloc) &
       &    *p12p3q4weight(ip12,ip3,iq4) &
       &    *recoupl(bl,bs,alphap,alphaloc)
          END DO ! ip12 
         END DO ! ip3 
        END DO  ! iq4
       END IF
      END DO ! alphaloc
     END DO ! alphap
    END DO ! bs 
   END DO ! bl 
   
   DEALLOCATE(psiamp4N,recoupl) 
   
  ! continue with 2+2 part 
    
   ALLOCATE(psiamp4N(mynp12,mynp34,mynq,beta4N22cdepmax))
   CALL collect_piece_spcmplx(PHI4N22%amp,psiamp4N,&
          &     mynp12*mynp34*mynq,beta4N22cdepmax,1,commalpha)

! prepare the recoupl coeffs for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    recoupl(L,S,alphap,alpha)
!       L total angular momentum 0:j4_max+2 
!       S total spin 0 : 2           
          
          
   ALLOCATE(recoupl(0:(j4max+2),0:2,beta4N22cdepmax,mynbeta4N22))
   
   recoupl=0.0_spreal
   
   DO alphap=1,beta4N22cdepmax
    CALL get4N22qn(alphap,l12p,s12p,j12p,t12p,l34p,s34p,j34p,t34p,lamp,Ip,j4p,tau4p,mtau4p,alpha12p,alpha34p,parip) 
    DO alphaloc=1,mynbeta4N22
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari) 
     IF(t12.eq.t12p .and. &
      & t34.eq.t34p .and. tau4.eq.tau4p .and. mtau4.eq.mtau4p .and. &
      & l12.eq.l12p .and. l34.eq.l34p .and. lam.eq.lamp .and. &
      & s12.eq.s12p .and. s34.eq.s34p .and. j4.eq.j4p) THEN
      DO bs=0,2 
       DO bl=abs(j4-bs),j4+bs 
        DO bl3=abs(I-s12),I+s12  
         recoupl(bl,bs,alphap,alphaloc) = recoupl(bl,bs,alphap,alphaloc) + & 
        &  (-1)**(l12+s12+j12)*SDCH(bl3)*SDCH(j12)  &
        &   *C6J(2*lam,2*l12,2*bl3,2*s12,2*I,2*j12) &
        &   *SDCH(I)*SDCH(j34)*SDCH(bl)*SDCH(bs)  &
        &   *C9J(2*bl3,2*s12,2*I,2*l34,2*s34,2*j34,2*bl,2*bs,2*j4) &  
        &   *(-1)**(l12+s12+j12p)*SDCH(bl3)*SDCH(j12p)  &
        &   *C6J(2*lam,2*l12,2*bl3,2*s12,2*Ip,2*j12p) &
        &   *SDCH(Ip)*SDCH(j34p)*SDCH(bl)*SDCH(bs)  &
        &   *C9J(2*bl3,2*s12,2*Ip,2*l34,2*s34,2*j34p,2*bl,2*bs,2*j4) 
        END DO  ! bl3 
       END DO ! bs
      END DO  ! bl   
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   ! do not set locpls to zero or allocate since the 3+1 part needs to be kept 
   DO bl=0,j4max+2
    DO bs=0,2 
     DO alphap=1,beta4N22cdepmax
      DO alphaloc=1,mynbeta4N22
       IF(recoupl(bl,bs,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
        DO iq=1,mynq
         DO ip34=1,mynp34
          DO ip12=1,mynp12
           locpls(bl,bs)=locpls(bl,bs) &
       & + 6.0*conjg(psiamp4N(ip12,ip34,iq,alphap)) &
            *PSI4N22%amp(ip12,ip34,iq,alphaloc) &
       &    *p12p34qweight(ip12,ip34,iq) &
       &    *recoupl(bl,bs,alphap,alphaloc)
          END DO ! ip12 
         END DO ! ip34
        END DO  ! iq
       END IF
      END DO ! alphaloc
     END DO ! alphap
    END DO ! bs 
   END DO ! bl 
   
   DEALLOCATE(psiamp4N,recoupl) 
    
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locpls,pls,3*(j4max+3),MPI_REAL4,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A)') 'LS probabilities from 4N wave and Yak-component: '
    norm=0.0_spreal
    DO bl=0,j4max+2
     DO bs=0,2 
      IF(pls(bl,bs).NE.0.0_spreal) THEN
       WRITE(*,'(A,2I4,2X,E15.6)') 'PLS-4N(WAVE-Yak): ',bl,bs,pls(bl,bs)
       norm=norm+pls(bl,bs)
      END IF
     END DO
    END DO
    WRITE(*,'(A,8X,2X,E15.6)') 'PLS-4N(WAVE-Yak): ',norm
   END IF
   
   DEALLOCATE(locpls,pls)
   
  END IF   ! 4N ampl  
  
 END SUBROUTINE printls_4n

SUBROUTINE printls_4n31(PSI4N31)
  IMPLICIT NONE
  TYPE(AMP4N31) PSI4N31 
  
  REAL(spreal),ALLOCATABLE :: locpls(:,:),pls(:,:)
  INTEGER ip12,ip3,iq4,alphaloc,alpha,alphap,ierr
  INTEGER l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p,alpha2bp,alpha3bp,parip
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari,mtau3,mt12,mt1,mt34,mt3
  INTEGER bl,bs,bl3,bs3
  COMPLEX(spreal),ALLOCATABLE :: psiamp4N(:,:,:,:)
  REAL(spreal),ALLOCATABLE :: recoupl(:,:,:,:) 
  REAL(spreal) :: norm 
  

  
  IF(allocated(PSI4N31%amp)) THEN ! calculate based on 4N wave and Yakubovsky components
! start with kinetic energy
! since there is a possibility to have isospin transitions
! I need to collect all channels 
! start with 3+1 part 
    
   ALLOCATE(psiamp4N(mynp12,mynp3,mynq4,alpha4N31cdepmax))
   CALL collect_piece_spcmplx(PSI4N31%amp,psiamp4N,&
          &     mynp12*mynp3*mynq4,alpha4N31cdepmax,1,commalpha)

! prepare the recoupl coeffs for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    recoupl(L,S,alphap,alpha)
!       L total angular momentum 0:j4_max+2 
!       S total spin 0 : 2        
          
          
   ALLOCATE(recoupl(0:(j4max+2),0:2,alpha4N31cdepmax,mynalpha4N31))
   
   recoupl=0.0_spreal
   
   DO alphap=1,alpha4N31cdepmax
    call get4N31qn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p,alpha2bp,alpha3bp,parip)
    DO alphaloc=1,mynalpha4N31
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari) 
     IF(t12.eq.t12p .and. &
      & tau3.eq.tau3p .and. tau4.eq.tau4p .and. mtau4.eq.mtau4p .and. &
      & l12.eq.l12p .and. l3.eq.l3p .and. l4.eq.l4p .and. &
      & s12.eq.s12p .and. j4.eq.j4p) THEN
      DO bs=0,2 
       DO bl=abs(j4-bs),j4+bs 
        DO bl3=abs(bl-l4),bl+l4
         DO bs3=abs(j3-2*bl3),j3+2*bl3,2  
          recoupl(bl,bs,alphap,alphaloc) = recoupl(bl,bs,alphap,alphaloc) + & 
        &  SDCH(bl3)*SDCH(j12)*SD2CH(I3)*SD2CH(bs3)  &
        &   *C9J(2*l12,2*s12,2*j12,2*l3,1,I3,2*bl3,bs3,j3) &
        &   *SD2CH(j3)*SD2CH(I4)*SDCH(bl)*SDCH(bs)  &
        &   *C9J(2*bl3,bs3,j3,2*l4,1,I4,2*bl,2*bs,2*j4) &        
        &   *SDCH(bl3)*SDCH(j12p)*SD2CH(I3p)*SD2CH(bs3)  &
        &   *C9J(2*l12,2*s12,2*j12p,2*l3,1,I3p,2*bl3,bs3,j3p) &
        &   *SD2CH(j3p)*SD2CH(I4p)*SDCH(bl)*SDCH(bs)  &
        &   *C9J(2*bl3,bs3,j3p,2*l4,1,I4p,2*bl,2*bs,2*j4)
         END DO ! bs3
        END DO  ! bl3 
       END DO ! bs
      END DO  ! bl   
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   ALLOCATE(locpls(0:(j4max+2),0:2),pls(0:(j4max+2),0:2))
   
   locpls=0.0_dpreal 
   DO bl=0,j4max+2
    DO bs=0,2 
     DO alphap=1,alpha4N31cdepmax
      DO alphaloc=1,mynalpha4N31
       IF(recoupl(bl,bs,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
        DO iq4=1,mynq4
         DO ip3=1,mynp3
          DO ip12=1,mynp12
           locpls(bl,bs)=locpls(bl,bs) &
       & +   conjg(psiamp4N(ip12,ip3,iq4,alphap)) &
       &    *PSI4N31%amp(ip12,ip3,iq4,alphaloc) &
       &    *p12p3q4weight(ip12,ip3,iq4) &
       &    *recoupl(bl,bs,alphap,alphaloc)
          END DO ! ip12 
         END DO ! ip3 
        END DO  ! iq4
       END IF
      END DO ! alphaloc
     END DO ! alphap
    END DO ! bs 
   END DO ! bl 
   
   DEALLOCATE(psiamp4N,recoupl) 
   
   
  ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locpls,pls,3*(j4max+3),MPI_REAL4,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A)') 'LS probabilities from 4N wave and Yak-component: '
    norm=0.0_spreal
    DO bl=0,j4max+2
     DO bs=0,2 
      IF(pls(bl,bs).NE.0.0_spreal) THEN
       WRITE(*,'(A,2I4,2X,E15.6)') 'PLS-4N(WAVE4N31): ',bl,bs,pls(bl,bs)
       norm=norm+pls(bl,bs)
      END IF
     END DO
    END DO
    WRITE(*,'(A,8X,2X,E15.6)') 'PLS-4N(WAVE4N31): ',norm
   END IF
   
   DEALLOCATE(locpls,pls)   

   END IF  ! 4N31 part   
 END SUBROUTINE printls_4n31

 SUBROUTINE printls_4n22(PSI4N22)
  IMPLICIT NONE
  TYPE(AMP4N22) PSI4N22
  
  REAL(spreal),ALLOCATABLE :: locpls(:,:),pls(:,:)
  INTEGER ip12,ip34,iq,alphaloc,alpha,alphap,ierr
  INTEGER l12p,s12p,j12p,t12p,j4p,tau4p,mtau4p,alpha2bp,parip
  INTEGER l12,s12,j12,t12,j4,tau4,mtau4,alpha2b,pari,mt12,mt1,mt34,mt3
  INTEGER l34p,s34p,j34p,t34p,lamp,Ip,l34,s34,j34,t34,lam,I,alpha12,alpha34,alpha12p,alpha34p,bl,bs,bl3,bs3
  COMPLEX(spreal),ALLOCATABLE :: psiamp4N(:,:,:,:)
  REAL(spreal),ALLOCATABLE :: recoupl(:,:,:,:) 
  REAL(spreal) :: norm   
  
   ! continue with 2+2 part 
  IF(allocated(PSI4N22%amp)) THEN ! calculate based on 4N wave and Yakubovsky components 
 
    
   ALLOCATE(psiamp4N(mynp12,mynp34,mynq,beta4N22cdepmax))
   CALL collect_piece_spcmplx(PSI4N22%amp,psiamp4N,&
          &     mynp12*mynp34*mynq,beta4N22cdepmax,1,commalpha)

! prepare the recoupl coeffs for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    recoupl(L,S,alphap,alpha)
!       L total angular momentum 0:j4_max+2 
!       S total spin 0 : 2           
          
          
   ALLOCATE(recoupl(0:(j4max+2),0:2,beta4N22cdepmax,mynbeta4N22))
   
   recoupl=0.0_spreal
   
   DO alphap=1,beta4N22cdepmax
    CALL get4N22qn(alphap,l12p,s12p,j12p,t12p,l34p,s34p,j34p,t34p,lamp,Ip,j4p,tau4p,mtau4p,alpha12p,alpha34p,parip) 
    DO alphaloc=1,mynbeta4N22
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari) 
     IF(t12.eq.t12p .and. &
      & t34.eq.t34p .and. tau4.eq.tau4p .and. mtau4.eq.mtau4p .and. &
      & l12.eq.l12p .and. l34.eq.l34p .and. lam.eq.lamp .and. &
      & s12.eq.s12p .and. s34.eq.s34p .and. j4.eq.j4p) THEN
      DO bs=0,2 
       DO bl=abs(j4-bs),j4+bs 
        DO bl3=abs(I-s12),I+s12  
         recoupl(bl,bs,alphap,alphaloc) = recoupl(bl,bs,alphap,alphaloc) + & 
        &  (-1)**(l12+s12+j12)*SDCH(bl3)*SDCH(j12)  &
        &   *C6J(2*lam,2*l12,2*bl3,2*s12,2*I,2*j12) &
        &   *SDCH(I)*SDCH(j34)*SDCH(bl)*SDCH(bs)  &
        &   *C9J(2*bl3,2*s12,2*I,2*l34,2*s34,2*j34,2*bl,2*bs,2*j4) &  
        &   *(-1)**(l12+s12+j12p)*SDCH(bl3)*SDCH(j12p)  &
        &   *C6J(2*lam,2*l12,2*bl3,2*s12,2*Ip,2*j12p) &
        &   *SDCH(Ip)*SDCH(j34p)*SDCH(bl)*SDCH(bs)  &
        &   *C9J(2*bl3,2*s12,2*Ip,2*l34,2*s34,2*j34p,2*bl,2*bs,2*j4) 
        END DO  ! bl3 
       END DO ! bs
      END DO  ! bl   
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   ALLOCATE(locpls(0:(j4max+2),0:2),pls(0:(j4max+2),0:2)) 
   locpls=0.0_dpreal  
   
   DO bl=0,j4max+2
    DO bs=0,2 
     DO alphap=1,beta4N22cdepmax
      DO alphaloc=1,mynbeta4N22
       IF(recoupl(bl,bs,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
        DO iq=1,mynq
         DO ip34=1,mynp34
          DO ip12=1,mynp12
           locpls(bl,bs)=locpls(bl,bs) &
       & +   conjg(psiamp4N(ip12,ip34,iq,alphap)) &
            *PSI4N22%amp(ip12,ip34,iq,alphaloc) &
       &    *p12p34qweight(ip12,ip34,iq) &
       &    *recoupl(bl,bs,alphap,alphaloc)
          END DO ! ip12 
         END DO ! ip34
        END DO  ! iq
       END IF
      END DO ! alphaloc
     END DO ! alphap
    END DO ! bs 
   END DO ! bl 
   
   DEALLOCATE(psiamp4N,recoupl) 
    
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locpls,pls,3*(j4max+3),MPI_REAL4,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A)') 'LS probabilities from 4N wave and Yak-component: '
    norm=0.0_spreal
    DO bl=0,j4max+2
     DO bs=0,2 
      IF(pls(bl,bs).NE.0.0_spreal) THEN
       WRITE(*,'(A,2I4,2X,E15.6)') 'PLS-4N(WAVE4N22): ',bl,bs,pls(bl,bs)
       norm=norm+pls(bl,bs)
      END IF
     END DO
    END DO
    WRITE(*,'(A,8X,2X,E15.6)') 'PLS-4N(WAVE4N22): ',norm
   END IF
   
   DEALLOCATE(locpls,pls)
   
  END IF   ! 4N22 ampl  
  
 END SUBROUTINE printls_4n22

 
 SUBROUTINE printinfo_4n(PSI4N31,PSI4N22,PHI4N31,PHI4N22)
  IMPLICIT NONE
  TYPE(AMP4N31) PSI4N31,PHI4N31 
  TYPE(AMP4N22) PSI4N22,PHI4N22 
  
  CALL printekin_4n(PSI4N31,PSI4N22,PHI4N31,PHI4N22)
  CALL printls_4n(PSI4N31,PSI4N22,PHI4N31,PHI4N22)
  
 END SUBROUTINE printinfo_4n

 
 SUBROUTINE printinfo_4n31(PSI4N31)
  IMPLICIT NONE
  TYPE(AMP4N31) PSI4N31 
  
  CALL printekin_4n31(PSI4N31)
  CALL printls_4n31(PSI4N31)
  
 END SUBROUTINE printinfo_4n31
 
 SUBROUTINE printinfo_4n22(PSI4N22)
  IMPLICIT NONE
  TYPE(AMP4N22) PSI4N22 
  
  CALL printekin_4n22(PSI4N22)
  CALL printls_4n22(PSI4N22)
  
 END SUBROUTINE printinfo_4n22
 ! for testing print out a few simple expectation values 
 ! 1) kinetic energy for each component and  norm for each component 
 ! 2) LS probabilities for each component 
 ! next case: 2N wave function
 
 
 SUBROUTINE printekin_nn(PSI)
  IMPLICIT NONE
  TYPE(AMPNN) PSI
  REAL(spreal) :: locekin,ekin,locnorm,norm
  INTEGER ip12,alphaloc,alpha,ierr
  INTEGER l12,s12,j12,t12,mt12
  INTEGER alphap,l12p,s12p,j12p,t12p,mt12p,mt1
  COMPLEX(spreal),ALLOCATABLE :: psiampNN(:,:)
  REAL(spreal),ALLOCATABLE :: isomass(:,:,:)
  REAL(spreal),ALLOCATABLE :: mbaryon(:,:),mthres(:)
  REAL(spreal) :: isofakt 
  
!  prepare masses depending on isospin 
  ALLOCATE(mbaryon(0:2,-2:2),mthres(-2:2))
  mbaryon=0.0_spreal
! set all baryons, only nucleons are used   
  mbaryon(0,0)=mlam
  mbaryon(1,-1)=mneu
  mbaryon(1,1)=mprot
  mbaryon(2,-2)=msigm
  mbaryon(2,0)=msig0
  mbaryon(2,2)=msigp
  mthres=0.0_spreal
  mthres(2)=2.0*mprot
  mthres(0)=mprot+mneu
  mthres(-2)=2.0*mneu
  
  
  IF(allocated(PSI%amp)) THEN ! calculate based on NN amplitude 
! start with kinetic energy
! since there is a possibility to have isospin transitions
! I need to collect all channels 
   ALLOCATE(psiampNN(mynp12,alphaNNcdepmax))
   CALL collect_piece_spcmplx(PSI%amp,psiampNN,&
          &     mynp12,alphaNNcdepmax,1,commalpha)

! prepare the isospin masses for each channel combination 
! very simple, just put zeros when it is zero, no bookkeeping here 
!    isomass(i,alphap,alpha)
!       i=1   total mass 
!       i=2   reduced mass average particles 12 
!       i=3   reduced mass average particle 3          
          
          
   ALLOCATE(isomass(1:2,alphaNNcdepmax,mynalphaNN)) 
   isomass=0.0_spreal
   
   DO alphap=1,alphaNNcdepmax
    CALL getNNqn(alphap,l12p,s12p,j12p,t12p,mt12p) 
    t12p=2*t12p   !   use convention that isospins are twice  
    mt12p=2*mt12p
    DO alphaloc=1,mynalphaNN
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL getNNqn(alpha,l12,s12,j12,t12,mt12)
     t12=2*t12   !   use convention that isospins are twice
     mt12=2*mt12
     IF(l12.eq.l12p .and. s12.eq.s12p .and.&
      & j12.eq.j12p .and.&
      &  mt12.eq.mt12p) THEN
       DO mt1=-1,1,2 
        IF(abs(mt12-mt1).GT.1) cycle
        isofakt = CG(1,1,t12,mt1,mt12-mt1,mt12) &
        &        *CG(1,1,t12p,mt1,mt12-mt1,mt12)
        
        isomass(1,alphap,alphaloc) = isomass(1,alphap,alphaloc) &
        & + isofakt &
        &   *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1) &
        &    - mthres(mt12))
        isomass(2,alphap,alphaloc) = isomass(2,alphap,alphaloc) &
        & + isofakt &
        &  *(mbaryon(1,mt1)+mbaryon(1,mt12-mt1))  & 
        &      /(2.0*mbaryon(1,mt1)*mbaryon(1,mt12-mt1))
       END DO  ! mt1 
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   locekin=0.0_dpreal 
   DO alphap=1,alphaNNcdepmax
    DO alphaloc=1,mynalphaNN
     IF(isomass(2,alphap,alphaloc).NE.0.0_spreal) THEN  ! use channel if coupled 
      DO ip12=1,mynp12
       locekin=locekin &
       & + conjg(psiampNN(ip12,alphap)) &
            *PSI%amp(ip12,alphaloc) &
       &    *p12weight(ip12) &
       &    *(isomass(1,alphap,alphaloc) &
       &    + isomass(2,alphap,alphaloc)*P12P(myp12+ip12-1)**2)   
      END DO ! ip12 
     END IF
    END DO ! alphaloc
   END DO ! alphap

   locnorm=0.0_dpreal 
   DO alphaloc=1,mynalphaNN
    DO ip12=1,mynp12
     locnorm=locnorm &
     & + conjg(PSI%amp(ip12,alphaloc)) &
          *PSI%amp(ip12,alphaloc) &
     &    *p12weight(ip12) 
    END DO ! ip12 
   END DO ! alphaloc
   
   DEALLOCATE(psiampNN,isomass) 
   
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locekin,ekin,1,MPI_REAL4,MPI_SUM,commampnn,ierr) ! here it is assumed that product is REAL(spreal)
   CALL MPI_ALLREDUCE(locnorm,norm,1,MPI_REAL4,MPI_SUM,commampnn,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A,E15.6)') 'kinetic energy of NN amplitude [MeV]: ',ekin*hbarc
    WRITE(*,'(A,E15.6)') 'Norm of NN  amplitude:         ',norm
   END IF
   
  END IF   ! 3N wave 
  
  DEALLOCATE(mbaryon,mthres) 
  
 END SUBROUTINE printekin_nn
 
 SUBROUTINE printls_nn(PSI)
  IMPLICIT NONE
  TYPE(AMPNN) PSI
  REAL(spreal),ALLOCATABLE :: locpls(:),pls(:)
  INTEGER ip12,ip3,alphaloc,alpha,ierr
  INTEGER l12,s12,j12,t12,mt12
  REAL(spreal) :: norm 
  
 
  IF(allocated(PSI%amp)) THEN ! calculate based on NN wave 
   ALLOCATE(locpls(alphaNNcdepmax),pls(alphaNNcdepmax))
 
 ! calculate the local sums 
   locpls=0.0_dpreal 
   DO alphaloc=1,mynalphaNN
    alpha=myalphaid+1+(alphaloc-1)*npe_alpha
    DO ip12=1,mynp12
     locpls(alpha)=locpls(alpha) &
     & + conjg(PSI%amp(ip12,alphaloc)) &
          *PSI%amp(ip12,alphaloc) &
     &    *p12weight(ip12) 
    END DO ! ip12 
   END DO ! alphaloc
   
   ! sum local results using MPI 
   CALL MPI_ALLREDUCE(locpls,pls,alphaNNcdepmax,MPI_REAL4,MPI_SUM,commampnn,ierr) ! here it is assumed that product is REAL(spreal)
   
   IF(master) THEN
    WRITE(*,'(A)') 'Channel probabilities from NN wave: '
    norm=0.0_spreal
    DO alpha=1,alphaNNcdepmax 
     CALL getnnqn(alpha,l12,s12,j12,t12,mt12)
     IF(abs(pls(alpha)).GT.50.0*eps_sp) THEN
      WRITE(*,'(A,5I4,2X,E15.6)') 'PCHAN-NN: ',l12,s12,j12,t12,mt12,pls(alpha)
      norm=norm+pls(alpha)
     END IF
    END DO
    WRITE(*,'(A,8X,2X,E15.6)') 'PCHAN-NN-NORM: ',norm
   END IF
   
   DEALLOCATE(locpls,pls)
   
  END IF   ! NN ampl  
  
 END SUBROUTINE printls_nn
 
 
 SUBROUTINE print_info_nn(PSI)
  IMPLICIT NONE
  TYPE(AMPNN) PSI 
  
  CALL printekin_nn(PSI)
  CALL printls_nn(PSI) 
  
 END SUBROUTINE print_info_nn

 
 ! printmomdist_3n calculates and prints the momentum distribition 
 ! for a given amplitude
 ! the calculation is done separately for protons and neutrons 
 ! and for all available charges of the amplitude
 ! also the overall normalization is calculated 
 
 SUBROUTINE printmomdist_3n(PSI)
  IMPLICIT NONE
  TYPE(AMP3N) PSI
  REAL(spreal) :: norm
  INTEGER ip12,ip3,alphaloc,alpha,ierr
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  INTEGER alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p, &
 &        mtau3p,alpha2bp,parip,mt3
  COMPLEX(spreal),ALLOCATABLE :: psiamp3N(:,:,:)
  REAL(spreal),ALLOCATABLE :: isofakt(:,:,:,:) ! (alpha,alphaloc,mt3,mtau3) 
  REAL(spreal),ALLOCATABLE :: locmomdist(:,:,:),redmomdist(:,:,:),&
  &                           momdist(:,:,:)  ! (ip3,mt3,mtau3) 
  REAL(spreal),ALLOCATABLE :: isonorm(:,:)  !  (ip3,mt3,mtau3) 
  

  
  IF(allocated(PSI%amp)) THEN ! calculate based on 3N amplitude 
! since there is a possibility to have isospin transitions
! I need to collect all channels 
   ALLOCATE(psiamp3N(mynp12,mynp3,alpha3Ncdepmax))
   CALL collect_piece_spcmplx(PSI%amp,psiamp3N,&
          &     mynp12*mynp3,alpha3Ncdepmax,1,commalpha)       
          
          
   ALLOCATE(isofakt(mynalpha3N,alpha3Ncdepmax,-1:1,mtau3min:mtau3max)) 
   isofakt=0.0
   DO alphap=1,alpha3Ncdepmax
    CALL get3Nqn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,&
       &          j3p,tau3p,mtau3p,alpha2bp,parip) 
    t12p=2*t12p   !   use convention that isospins are twice    
    DO alphaloc=1,mynalpha3N
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)
     t12=2*t12   !   use convention that isospins are twice
     IF(l12.eq.l12p .and. s12.eq.s12p .and.&
      & j12.eq.j12p .and.  l3.eq.l3p  .and.&
      &  I3.eq.I3p  .and.  j3.eq.j3p  .and.&
      & t12.eq.t12p .and. &
      &  mtau3.eq.mtau3p) THEN
      DO mt3=-1,1,2 
       IF(abs(mtau3-mt3).GT.t12) cycle
       isofakt(alphaloc,alphap,mt3,mtau3) = CG(t12,1,tau3,mtau3-mt3,mt3,mtau3) &
       &                         *CG(t12p,1,tau3p,mtau3-mt3,mt3,mtau3)
      END DO ! mt3   
     END IF ! qn numbers 
    END DO ! alphaloc
   END DO  ! alphap 
   
   ALLOCATE(isonorm(-1:1,mtau3min:mtau3max),&
     & locmomdist(mynp3,-1:1,mtau3min:mtau3max),&
     & momdist(P3N,-1:1,mtau3min:mtau3max),&
     & redmomdist(mynp3,-1:1,mtau3min:mtau3max))
   isonorm=0.0
   locmomdist=0.0
   
   DO mtau3=mtau3min,mtau3max,2
    DO mt3=-1,1,2
     DO alphap=1,alpha3Ncdepmax
      DO alphaloc=1,mynalpha3N
       IF(abs(isofakt(alphaloc,alphap,mt3,mtau3)).GT.50.0*eps_sp) THEN  ! use channel if coupled 
        DO ip3=1,mynp3
         DO ip12=1,mynp12
          locmomdist(ip3,mt3,mtau3)=locmomdist(ip3,mt3,mtau3) &
       &   + conjg(psiamp3N(ip12,ip3,alphap)) &
       &    *PSI%amp(ip12,ip3,alphaloc) &
       &    *p12weight(ip12) &
       &    *isofakt(alphaloc,alphap,mt3,mtau3)
         END DO ! ip12 
        END DO ! ip3 
       END IF
      END DO ! alphaloc
     END DO ! alphap
    END DO ! mt3
   END DO  ! mtau3 

   DEALLOCATE(psiamp3N,isofakt) 
     
   ! global sums over ip12 and alpha
   CALL MPI_ALLREDUCE(locmomdist,redmomdist,mynp3*3*(mtau3max-mtau3min+1),&
     &                MPI_REAL4,MPI_SUM,commp12alpha,ierr) ! here it is assumed that product is REAL(spreal)
   
   ! then collect the result on all processors 
     
   CALL collect_part_spreal(redmomdist,momdist,1,p3n,&
              &             3*(mtau3max-mtau3min+1),commp3)  
              
   ! calc the norm to check what has to be printed 
              
   DO mtau3=mtau3min,mtau3max,2
    DO mt3=-1,1,2 
     DO ip3=1,P3N
      isonorm(mt3,mtau3)=isonorm(mt3,mtau3) &
     & + momdist(ip3,mt3,mtau3) &
     &    *P3W(ip3)*P3P(ip3)**2
     END DO ! ip3 
    END DO ! mt3 
   END DO ! mtau3
   
   norm=0.0
   DO mtau3=mtau3min,mtau3max,2
    DO mt3=-1,1,2 
     IF(isonorm(mt3,mtau3).GT.50.0*eps_sp) THEN
      IF(master) THEN 
       WRITE(*, &
       &'("Probability for isospin state for mt3= ",I4,&
       &  "/2  mtau3= ",I4,"/2: ",E15.6)') &
       &  mt3,mtau3,isonorm(mt3,mtau3)
       DO ip3=1,P3N
        WRITE(*,'(A,2I4,2E15.6)') 'PROPP1N: ',mt3,mtau3,P3P(ip3),momdist(ip3,mt3,mtau3)
       END DO
       
      END IF
     END IF
     norm = norm + isonorm(mt3,mtau3)
    END DO ! mt3
   END DO ! mtau3 
   
   If(master) THEN
    WRITE(*,'(A,E15.6)') 'Norm in momentum distribution: ',norm
   END IF
 
   DEALLOCATE(isonorm,locmomdist,momdist,redmomdist)
  END IF   ! 3N wave 
  
  
 END SUBROUTINE printmomdist_3n
 
 ! printmomcorr_3n calculates and prints the momentum distribition 
 ! for a given amplitude
 ! the calculation is done separately for pp, np, and nn pairs 
 ! and for all available charges of the amplitude
 ! also the overall normalization is calculated 
 
 SUBROUTINE printmomcorr_3n(PSI)
  IMPLICIT NONE
  TYPE(AMP3N) PSI
  REAL(spreal) :: norm
  INTEGER ip12,ip3,alphaloc,alpha,ierr
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  INTEGER alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p, &
 &        mtau3p,alpha2bp,parip,mt12
  COMPLEX(spreal),ALLOCATABLE :: psiamp3N(:,:,:)
  REAL(spreal),ALLOCATABLE :: isofakt(:,:,:,:) ! (alpha,alphaloc,mt12,mtau3) 
  REAL(spreal),ALLOCATABLE :: locmomdist(:,:,:),redmomdist(:,:,:),&
  &                           momdist(:,:,:)  ! (ip12,mt12,mtau3) 
  REAL(spreal),ALLOCATABLE :: isonorm(:,:)  !  (ip12,mt12,mtau3) 

  

  
  IF(allocated(PSI%amp)) THEN ! calculate based on NNY amplitude 
! since there is a possibility to have isospin transitions
! I need to collect all channels 
   ALLOCATE(psiamp3N(mynp12,mynp3,alpha3Ncdepmax))
   CALL collect_piece_spcmplx(PSI%amp,psiamp3N,&
          &     mynp12*mynp3,alpha3Ncdepmax,1,commalpha)    
          
          
   ALLOCATE(isofakt(mynalpha3N,alpha3Ncdepmax,-2:2,mtau3min:mtau3max)) 
   isofakt=0.0
   DO alphap=1,alpha3Ncdepmax
    CALL get3Nqn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,&
       &          j3p,tau3p,mtau3p,alpha2bp,parip) 
    t12p=2*t12p   !   use convention that isospins are twice    
    DO alphaloc=1,mynalpha3N
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)
     CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)
     t12=2*t12   !   use convention that isospins are twice
     IF(l12.eq.l12p .and. s12.eq.s12p .and.&
      & j12.eq.j12p .and.  l3.eq.l3p  .and.&
      &  I3.eq.I3p  .and.  j3.eq.j3p  .and.&
      & t12.eq.t12p .and. &
      &  mtau3.eq.mtau3p) THEN
      DO mt12=-2,2,2 
       IF(abs(mt12).GT.t12) cycle
       isofakt(alphaloc,alphap,mt12,mtau3) = CG(t12,1,tau3,mt12,mtau3-mt12,mtau3) &
       &                         *CG(t12p,1,tau3p,mt12,mtau3-mt12,mtau3)
      END DO ! mt12  
     END IF
    END DO ! alphaloc
   END DO  ! alphap 
   
   ALLOCATE(isonorm(-2:2,mtau3min:mtau3max),&
     & locmomdist(mynp12,-2:2,mtau3min:mtau3max),&
     & momdist(P12N,-2:2,mtau3min:mtau3max),&
     & redmomdist(mynp12,-2:2,mtau3min:mtau3max))
   isonorm=0.0
   locmomdist=0.0
   
   DO mtau3=mtau3min,mtau3max,2
    DO mt12=-2,2,2
     DO alphap=1,alpha3Ncdepmax
      DO alphaloc=1,mynalpha3N
       IF(abs(isofakt(alphaloc,alphap,mt12,mtau3)).GT.50.0*eps_sp) THEN  ! use channel if coupled 
        DO ip3=1,mynp3
         DO ip12=1,mynp12
          locmomdist(ip12,mt12,mtau3)=locmomdist(ip12,mt12,mtau3) &
       &   + conjg(psiamp3N(ip12,ip3,alphap)) &
       &    *PSI%amp(ip12,ip3,alphaloc) &
       &    *P3P(myp3+ip3-1)**2*P3W(myp3+ip3-1) &
       &    *isofakt(alphaloc,alphap,mt12,mtau3)
         END DO ! ip12 
        END DO ! ip3 
       END IF
      END DO ! alphaloc
     END DO ! alphap
    END DO ! mt3
   END DO  ! mtau3 

   DEALLOCATE(psiamp3N,isofakt) 
     
   ! global sums over ip12 and alpha
   CALL MPI_ALLREDUCE(locmomdist,redmomdist,mynp12*5*(mtau3max-mtau3min+1),&
     &                MPI_REAL4,MPI_SUM,commp3alpha,ierr) ! here it is assumed that product is REAL(spreal)
   
   ! then collect the result on all processors 
     
   CALL collect_part_spreal(redmomdist,momdist,1,p12n,&
              &             5*(mtau3max-mtau3min+1),commp12)  
              
   ! calc the norm to check what has to be printed 
              
   DO mtau3=mtau3min,mtau3max,2
    DO mt12=-2,2,2 
     DO ip12=1,P12N
      isonorm(mt12,mtau3)=isonorm(mt12,mtau3) &
     & + momdist(ip12,mt12,mtau3) &
     &    *P12W(ip12)*P12P(ip12)**2
     END DO ! ip3 
    END DO ! mt3 
   END DO ! mtau3
   
   norm=0.0
   DO mtau3=mtau3min,mtau3max,2
    DO mt12=-2,2,2 
     IF(isonorm(mt12,mtau3).GT.50.0*eps_sp) THEN
      IF(master) THEN 
       WRITE(*, &
       &'("Probability for isospin state for mt12= ",I4,&
       &  "/2  mtau3= ",I4,"/2: ",E15.6)') &
       &  mt12,mtau3,isonorm(mt12,mtau3)
       DO ip12=1,P12N
        WRITE(*,'(A,2I4,2E15.6)') 'PROPC2N: ',mt12,mtau3,P12P(ip12),momdist(ip12,mt12,mtau3)
       END DO
       
      END IF
     END IF
     norm = norm + isonorm(mt12,mtau3)
    END DO ! mt3
   END DO ! mtau3 
   
   If(master) THEN
    WRITE(*,'(A,E15.6)') 'Norm in momentum correlation: ',norm
   END IF
 
   DEALLOCATE(isonorm,locmomdist,momdist,redmomdist)
  END IF   ! 3N wave 
  
  
 END SUBROUTINE printmomcorr_3n
 
 SUBROUTINE project_LS_3n(bl,bs,PSIIN,PSIOUT)
  IMPLICIT NONE
  INTEGER bl,bs 
  TYPE(AMP3N) :: PSIIN,PSIOUT
  COMPLEX(spreal),ALLOCATABLE :: PSILS(:,:,:),PSILSLOC(:,:,:)
  REAL(spreal) :: faktls 
  INTEGER localpha,alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,pari,alpha2b
  INTEGER alphap,l12p,l3p,blp,s12p,bsp,j3p,t12p,tau3p,mtau3p,parip
  INTEGER ip12,ip3,ierr
  
  
  !  trafo to LS channels 
  
  ALLOCATE(PSILSLOC(mynp12,mynp3,alpha3NLScdepmax))
  ALLOCATE(PSILS(mynp12,mynp3,alpha3NLScdepmax))
  PSILSLOC=0.0
  
  DO alphap=1,alpha3NLScdepmax
   CALL get3NLSqn(alphap,l12p,l3p,blp,s12p,bsp,j3p,t12p,tau3p,mtau3p,parip)
   DO localpha=1,mynalpha3N   
    alpha=myalphaid+1+(localpha-1)*npe_alpha
    CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)  
    IF(l12.eq.l12p.AND.l3.eq.l3p.AND.s12.eq.s12p.AND.t12.EQ.t12p  &
    &   .AND.j3.EQ.j3p.AND.pari.EQ.parip.AND.tau3.EQ.tau3p.AND.mtau3.EQ.mtau3p) THEN
 
     faktls=SDCH(blp)*SDCH(j12)*SD2CH(I3)*SD2CH(bsp)  &
        &   *C9J(2*l12,2*s12,2*j12,2*l3,1,I3,2*blp,bsp,j3)
        
     DO ip3=1,mynp3
      DO ip12=1,mynp12
       PSILSLOC(ip12,ip3,alphap)=PSILSLOC(ip12,ip3,alphap)  &
    &   +faktls*PSIIN%amp(ip12,ip3,localpha)
      END DO
     END DO
       
    END IF ! qn numbers 
     
   END DO ! localpha
  END DO  ! alphap
  
  CALL MPI_ALLREDUCE(PSILSLOC,PSILS,mynp12*mynp3*alpha3NLScdepmax,MPI_COMPLEX8,MPI_SUM,commalpha,ierr) 
  DEALLOCATE(PSILSLOC)
  
  
  ! then trafo back to orginal jj including projection on bl and bs 
  PSIOUT=0.0
  
  DO alphap=1,alpha3NLScdepmax
   CALL get3NLSqn(alphap,l12p,l3p,blp,s12p,bsp,j3p,t12p,tau3p,mtau3p,parip)
   DO localpha=1,mynalpha3N   
    alpha=myalphaid+1+(localpha-1)*npe_alpha
    CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)  
    IF(l12.eq.l12p.AND.l3.eq.l3p.AND.s12.eq.s12p.AND.t12.EQ.t12p  &
    &   .AND.j3.EQ.j3p.AND.pari.EQ.parip.AND.tau3.EQ.tau3p.AND.mtau3.EQ.mtau3p &
    &   .AND. blp.EQ.bl .AND. bsp.EQ.bs) THEN  ! this line includes the projection 
 
     faktls=SDCH(blp)*SDCH(j12)*SD2CH(I3)*SD2CH(bsp)  &
        &   *C9J(2*l12,2*s12,2*j12,2*l3,1,I3,2*blp,bsp,j3)
        
     DO ip3=1,mynp3
      DO ip12=1,mynp12
       PSIOUT%amp(ip12,ip3,localpha)=PSIOUT%amp(ip12,ip3,localpha)  &
    &   +faktls*PSILS(ip12,ip3,alphap)
      END DO
     END DO
       
    END IF ! qn numbers 
     
   END DO ! localpha
  END DO  ! alphap
  
  DEALLOCATE(PSILS)

 END SUBROUTINE  project_LS_3n
 
 
 
 SUBROUTINE project_principalS_3n(PSIIN,PSIOUT)
  IMPLICIT NONE
  TYPE(AMP3N) :: PSIIN,PSIOUT
  COMPLEX(spreal),ALLOCATABLE :: PSILS(:,:,:),PSILSLOC(:,:,:)
  REAL(spreal) :: faktls 
  INTEGER localpha,alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,pari,alpha2b
  INTEGER alphap,l12p,l3p,blp,s12p,bsp,j3p,t12p,tau3p,mtau3p,parip
  INTEGER alphapp,l12pp,l3pp,blpp,s12pp,bspp,j3pp,t12pp,tau3pp,mtau3pp,paripp
  INTEGER ip12,ip3,ierr
  
  
  !  trafo to LS channels 
  
  ALLOCATE(PSILSLOC(mynp12,mynp3,alpha3NLScdepmax))
  ALLOCATE(PSILS(mynp12,mynp3,alpha3NLScdepmax))
  PSILSLOC=0.0
  
  DO alphap=1,alpha3NLScdepmax
   CALL get3NLSqn(alphap,l12p,l3p,blp,s12p,bsp,j3p,t12p,tau3p,mtau3p,parip)
   DO localpha=1,mynalpha3N   
    alpha=myalphaid+1+(localpha-1)*npe_alpha
    CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)  
    IF(l12.eq.l12p.AND.l3.eq.l3p.AND.s12.eq.s12p.AND.t12.EQ.t12p  &
    &   .AND.j3.EQ.j3p.AND.pari.EQ.parip.AND.tau3.EQ.tau3p.AND.mtau3.EQ.mtau3p) THEN
 
     faktls=SDCH(blp)*SDCH(j12)*SD2CH(I3)*SD2CH(bsp)  &
        &   *C9J(2*l12,2*s12,2*j12,2*l3,1,I3,2*blp,bsp,j3)
        
     DO ip3=1,mynp3
      DO ip12=1,mynp12
       PSILSLOC(ip12,ip3,alphap)=PSILSLOC(ip12,ip3,alphap)  &
    &   +faktls*PSIIN%amp(ip12,ip3,localpha)
      END DO
     END DO
       
    END IF ! qn numbers 
     
   END DO ! localpha
  END DO  ! alphap
  
  CALL MPI_ALLREDUCE(PSILSLOC,PSILS,mynp12*mynp3*alpha3NLScdepmax,MPI_COMPLEX8,MPI_SUM,commalpha,ierr) 
  DEALLOCATE(PSILSLOC)
  
  
  ! then trafo back to orginal jj including projection on bl and bs 
  PSIOUT=0.0
  
  DO alphapp=1,alpha3NLScdepmax
   CALL get3NLSqn(alphapp,l12pp,l3pp,blpp,s12pp,bspp,j3pp,t12pp,tau3pp,mtau3pp,paripp)
    DO alphap=1,alpha3NLScdepmax
    CALL get3NLSqn(alphap,l12p,l3p,blp,s12p,bsp,j3p,t12p,tau3p,mtau3p,parip)
    IF(l12p.eq.l12pp.AND.l3p.EQ.l3pp.AND.blp.EQ.blpp.AND.bsp.EQ.bspp  &
    &   .AND.j3p.EQ.j3pp.AND.tau3p.EQ.tau3pp.AND.mtau3p.EQ.mtau3pp.AND.parip.EQ.paripp &
    &   .AND.bsp.EQ.1.AND.blp.EQ.0.AND.j3p.EQ.1.AND.tau3p.EQ.1) THEN   ! projection on L=0 and S=1/2 (J3=1/2) and T3=1/2 
     DO localpha=1,mynalpha3N   
      alpha=myalphaid+1+(localpha-1)*npe_alpha
      CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)  
      IF(l12.eq.l12p.AND.l3.eq.l3p.AND.s12.eq.s12p.AND.t12.EQ.t12p  &
      &   .AND.j3.EQ.j3p.AND.pari.EQ.parip.AND.tau3.EQ.tau3p.AND.mtau3.EQ.mtau3p ) THEN 

       faktls=SDCH(blp)*SDCH(j12)*SD2CH(I3)*SD2CH(bsp)  &
          &   *C9J(2*l12,2*s12,2*j12,2*l3,1,I3,2*blp,bsp,j3)  &
          &   *0.5*(-1)**(s12pp-s12p)
          
       DO ip3=1,mynp3
        DO ip12=1,mynp12
         PSIOUT%amp(ip12,ip3,localpha)=PSIOUT%amp(ip12,ip3,localpha)  &
      &   +faktls*PSILS(ip12,ip3,alphapp)
        END DO
       END DO

      END IF ! qn numbers 

     END DO ! localpha
    END IF ! alphapp / alphap qn
   END DO  ! alphap
  END DO   ! alphapp
  
  DEALLOCATE(PSILS)
 END SUBROUTINE project_principalS_3n
 
 ! 4N31 amplitude
 ! open file <filename>
 ! write amplitude PSI to file <filename>
 ! in the group ampname
 ! file is then closed 
 
 SUBROUTINE writeamp_4n31(PSI,filename,ampname)
  IMPLICIT NONE 
  TYPE(AMP4N31),TARGET :: PSI 
  CHARACTER(LEN=*) filename
  CHARACTER(LEN=*) ampname
  INTEGER(HID_T) writef_id
  INTEGER(HID_T) group_amp_id
  INTEGER(HID_T) dset_nn_id,dsp_nn_id,dsp_hyper_id
  INTEGER(HID_T) msp_nn_id,msp_hyper_id
  INTEGER(HSIZE_T) dims(5)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
  INTEGER(HSIZE_T) start(5)     ! storage for offset in output arrays
  INTEGER(HSIZE_T) block(5)     ! storage for blocksize in output arrays
  INTEGER(HSIZE_T) stride(5)    ! storage for stride in output arrays
  INTEGER(HSIZE_T) count(5)     ! storage for count in output arrays
  INTEGER ierr
  REAL(spreal),ALLOCATABLE  :: buf(:,:,:,:,:)
  INTEGER(lint) blocksize

  writeamptime=writeamptime-MPI_WTIME()
  
  ! first open the hdf file
  ! is always assumed that the file exits
  ! copy empty file if necessary
  
  CALL open_write_h5file(filename,commall,writef_id)
  
  ! then create the group <ampname>

  CALL h5gcreate_f(writef_id,trim(ampname),group_amp_id, ierr)

  ! file should contain all info to interprete the data 
  !    1. meshpoints for P12 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  CALL writep12p(group_amp_id)
   
  !    2. meshpoints for P3 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  CALL writep3p(group_amp_id)

  !    3. meshpoints for Q4 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  CALL writeq4p(group_amp_id)
  
  !    4. 4N31 channels qnalpha4N31(:,:) including number of channels 
  CALL write_4n31_channels(group_amp_id)
  
  !    5. physical constants from constant mod are also written 
  CALL writephysconst(group_amp_id) 
  
  !    6. 4N31 amplitude if it exists 
  
  hdftime=hdftime-MPI_WTIME()
  
  IF(allocated(PSI%amp)) THEN
  ! this needs to be written in parallel since the data 
  ! is distributed 
  ! make sure that that only one processor writes the data 
  ! in case the distribution is not complete 
  
  ! here data is distributed in in commp12 and commalpha 
  ! that writing processors should have myp3id=myq4id=myenerid=0
  ! put the hyperslabs for other processors to empty 
  
  !     need to create a data set in the file 
  !     with complete number amp-elements 
  !     outline shape of global set of data 
   dims(1)=2                  !  complex = 2 * real
   dims(2)=P12N          !  write needs dimensions of arrays 
   dims(3)=P3N
   dims(4)=Q4N
   dims(5)=alpha4N31cdepmax     !  here array size in data file  
  !     create corresponding data space            
   CALL h5screate_simple_f(5,dims,dsp_nn_id,ierr);

  !     describe the memory layout of the data
   dims(1)=2                  ! complex = 2 * real
   dims(2)=mynp12        ! write needs dimension of arrays
   dims(3)=mynp3
   dims(4)=mynq4
   dims(5)=mynalpha4N31         ! here size of array in memory 
  !     create corresponding mem space    
   CALL h5screate_simple_f(5,dims,msp_nn_id,ierr);
  
  ! select the hyperslap in the file  that corresponds to the local 
  ! 4N31 data 
   IF(myenerid.EQ.0) THEN
    start(1)=0
    start(2)=myp12-1
    start(3)=myp3-1
    start(4)=myq4-1
    start(5)=myalphaid 
    block(1)=2
    block(2)=mynp12
    block(3)=mynp3
    block(4)=mynq4
    block(5)=1
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=1
    stride(5)=npe_alpha
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=1
    count(5)=mynalpha4N31
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF
  
   CALL h5scopy_f(dsp_nn_id,dsp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(dsp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)

  ! select the hyperslap in memorythat corresponds to the local 
  ! 4N31 data 

   IF(myenerid.EQ.0) THEN
    start(1)=0
    start(2)=0             
    start(3)=0
    start(4)=0
    start(5)=0
    block(1)=2
    block(2)=mynp12
    block(3)=mynp3
    block(4)=mynq4
    block(5)=mynalpha4N31
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=1
    stride(5)=1
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=1
    count(5)=1
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF   
     
   CALL h5scopy_f(msp_nn_id,msp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(msp_hyper_id,H5S_SELECT_SET_F,&
      &                      start, count,ierr,stride,block)
       
      
!    create the dataset in file 
   CALL h5dcreate_f(group_amp_id,trim('4n31amp'),H5T_NATIVE_REAL,&
      &             dsp_nn_id,dset_nn_id,ierr)
    

      
! resolve compatiblity issue by copying to REAL buffer 
   ALLOCATE(buf(2,mynp12,mynp3,mynq4,mynalpha4N31))
   buf(1,:,:,:,:)=REAL(PSI%amp(:,:,:,:))
   buf(2,:,:,:,:)=AIMAG(PSI%amp(:,:,:,:))

#ifdef HDF_MPIO
      ! check for consistency with MPI
   blocksize=4_lint*2*mynp12*mynp3*mynq4*mynalpha4N31  
   IF(blocksize.GE.blocklimit) THEN
    WRITE(*,*) 'block too large to be handled with MPI-IO'
    WRITE(*,*) 'blocksize = ',blocksize
    WRITE(*,*) 'blocklimit = ',blocklimit
    CALL abort
   END IF
#endif
   
! write data of amplitude (use here collective communication)       
   CALL h5dwrite_f(dset_nn_id, H5T_NATIVE_REAL,& 
     &      buf,dims,ierr,&
     &      mem_space_id = msp_hyper_id,&
     &      file_space_id = dsp_hyper_id,& 
     &      xfer_prp = pcollectwrite_id)

   DEALLOCATE(buf)

! close data set   
   CALL h5dclose_f(dset_nn_id,ierr)
  
! do not need mem_space anymore      
   CALL h5sclose_f(msp_hyper_id,ierr)
   CALL h5sclose_f(msp_nn_id,ierr)
! and not data space anymore    
   CALL h5sclose_f(dsp_hyper_id,ierr)
   CALL h5sclose_f(dsp_nn_id,ierr) 
  END IF  ! 3N exists 
  
  hdftime=hdftime+MPI_WTIME()
  
! close the group for this amplitude 
   CALL h5gclose_f(group_amp_id,ierr) 
! finally close the file 
   CALL close_h5file(writef_id)  
  
  writeamptime=writeamptime+MPI_WTIME()
  
 END SUBROUTINE writeamp_4n31

 ! 4N31  amplitude
 ! open file <filename>
 ! read amplitude PSI from group ampname 
 ! of file <filename>
 ! file is then closed
 ! print the parameters and give back the amplitude
 
 SUBROUTINE readamp_4n31(PSI,filename,ampname)
  IMPLICIT NONE 
  TYPE(AMP4N31),TARGET :: PSI 
  CHARACTER(LEN=*) filename
  CHARACTER(LEN=*) ampname
  INTEGER(HID_T) readf_id
  INTEGER(HID_T) group_amp_id
  INTEGER(HID_T) dset_nn_id,dsp_nn_id,dsp_hyper_id
  INTEGER(HID_T) msp_nn_id,msp_hyper_id
  INTEGER(HSIZE_T) dims(5)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
  INTEGER(HSIZE_T) start(5)     ! storage for offset in output arrays
  INTEGER(HSIZE_T) block(5)     ! storage for blocksize in output arrays
  INTEGER(HSIZE_T) stride(5)    ! storage for stride in output arrays
  INTEGER(HSIZE_T) count(5)     ! storage for count in output arrays
  INTEGER ierr,p12n_read,p3n_read,q4n_read
  INTEGER alpha4Ncdepmax_read,alpha4Nmax_read
  INTEGER ip12,ip3,iq4,alpha
  ! buffer to read the complex data 
  REAL(spreal),ALLOCATABLE  :: buf(:,:,:,:,:),&
     &                         psi_tmp_p12(:,:,:,:,:), &
     &                         psi_inter_p12(:,:,:,:,:),&
     &                         psi_tmp_p3(:,:,:,:,:), &
     &                         psi_inter_p3(:,:,:,:,:),&
     &                         psi_tmp_q4(:,:,:,:,:), &
     &                         psi_inter_q4(:,:,:,:,:),&     
     &                         psi_tmp_alpha(:,:,:,:,:)
  REAL(dpreal),ALLOCATABLE  :: p12p_read(:),p12w_read(:),&
                       &       p3p_read(:),p3w_read(:),&
                       &       q4p_read(:),q4w_read(:)
  INTEGER,ALLOCATABLE  :: qnalpha4N_read(:,:)
  INTEGER mynp12_read,myp12_read,myendp12_read
  INTEGER mynp3_read,myp3_read,myendp3_read
  INTEGER mynq4_read,myq4_read,myendq4_read
  INTEGER mynalpha4N_read
  REAL(dpreal),ALLOCATABLE :: spl12(:,:),spl3(:,:),spl4(:,:)
  INTEGER,ALLOCATABLE :: indx12(:,:),indx3(:,:),indx4(:,:),alphap(:)
  LOGICAL amp_exists
  REAL(dpreal) :: fakt 
  INTEGER alphaloc,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari
  INTEGER(lint) blocksize

  readamptime=readamptime-MPI_WTIME()
  
  ! first open the hdf file
  ! is always assumed that the file exits
  ! copy empty file if necessary
  CALL open_read_h5file(filename,commall,readf_id)
  
  ! then create the group <ampname>
  
  CALL h5gopen_f(readf_id,trim(ampname),group_amp_id, ierr)

  ! file should contain all info to interprete the data 
  !    1. meshpoints for P12/P3/Q4 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  !       the meshpoints and weights are allocated within the routine
  CALL readp12p(group_amp_id,p12p_read,p12w_read,p12n_read)
  CALL readp3p(group_amp_id,p3p_read,p3w_read,p3n_read)
  CALL readq4p(group_amp_id,q4p_read,q4w_read,q4n_read)
  
  !    2. 4N31 channels qn_read(:,:) including number of channels
  !       the quantum numbers are allocated within the routine 
  CALL read_4n31_channels(group_amp_id,qnalpha4N_read,&
     &                   alpha4Nmax_read,alpha4Ncdepmax_read)
  
  !    3. physical constants from constant mod are also read 
  !       and printed. No values are used. 
  CALL readphysconst(group_amp_id) 

! prepare the spline elements to interpolate
! from p12p_read to p12p 
! and p3p_read to p3p  
! and q4p_read to q4p 
  
  ALLOCATE(spl12(4,mynp12),indx12(4,mynp12))
  ALLOCATE(spl3(4,mynp3),indx3(4,mynp3))
  ALLOCATE(spl4(4,mynq4),indx4(4,mynq4))
  
  CALL cubfast_dp(p12p_read,p12n_read,&
        &    p12p(myp12:myp12+mynp12-1),mynp12, &
        &    spl12,indx12)
  CALL cubfast_dp(p3p_read,p3n_read,&
        &    p3p(myp3:myp3+mynp3-1),mynp3, &
        &    spl3,indx3)
  CALL cubfast_dp(q4p_read,q4n_read,&
        &    q4p(myq4:myq4+mynq4-1),mynq4, &
        &    spl4,indx4)
        
  ! distribution of mesh and channels 
  CALL distr_block(p12n_read,npe_p12,myp12id,myp12_read,myendp12_read,mynp12_read)
  CALL distr_block(p3n_read,npe_p3,myp3id,myp3_read,myendp3_read,mynp3_read)
  CALL distr_block(q4n_read,npe_q4,myq4id,myq4_read,myendq4_read,mynq4_read)
  CALL distr_piece(alpha4Ncdepmax_read,npe_alpha,myalphaid,mynalpha4N_read)
  
  !    4. 4N amplitude if it exists 
  
  CALL h5lexists_f(group_amp_id, trim('4n31amp'),amp_exists,ierr)
  
  hdftime=hdftime-MPI_WTIME()
  
  IF(amp_exists) THEN
  ! first read the amplitude using the mesh stored 
  ! prepare buffer for real and imaginary part of the amplitude 
  
   ALLOCATE(buf(2,mynp12_read,mynp3_read,mynq4_read,mynalpha4N_read))

  
  ! here data is distributed in in commp12 and commalpha 
  ! the reading processors should have myenerid=0
  ! this is element 0 in commener
  ! put the hyperslabs for other processors to empty 

  !     need to open the data set in the file 
   CALL h5dopen_f(group_amp_id,trim('4n31amp'),dset_nn_id,ierr)
  
  ! get the data space descriptor 
   CALL h5dget_space_f(dset_nn_id, dsp_nn_id, ierr)
     
  !     describe the memory layout of the data
   dims(1)=2                  ! complex = 2 * real
   dims(2)=mynp12_read        ! read needs dimension of arrays
   dims(3)=mynp3_read
   dims(4)=mynq4_read
   dims(5)=mynalpha4N_read      ! here size of array in memory 
  !     create corresponding mem space    
   CALL h5screate_simple_f(5,dims,msp_nn_id,ierr);
  
  ! select the hyperslap in the file  that corresponds to the local 
  ! 4N data 
   IF(myenerid.EQ.0) THEN
    start(1)=0
    start(2)=myp12_read-1
    start(3)=myp3_read-1
    start(4)=myq4_read-1
    start(5)=myalphaid 
    block(1)=2
    block(2)=mynp12_read
    block(3)=mynp3_read
    block(4)=mynq4_read
    block(5)=1
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=1
    stride(5)=npe_alpha
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=1
    count(5)=mynalpha4N_read
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF
  
   CALL h5scopy_f(dsp_nn_id,dsp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(dsp_hyper_id,H5S_SELECT_SET_F,&
      &                      start, count,ierr,stride,block)

  ! select the hyperslap in memorythat corresponds to the local 
  ! 4N31 data 
   IF(myenerid.EQ.0) THEN
    start(1)=0
    start(2)=0             
    start(3)=0
    start(4)=0
    start(5)=0
    block(1)=2
    block(2)=mynp12_read
    block(3)=mynp3_read 
    block(4)=mynq4_read 
    block(5)=mynalpha4N_read
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=1
    stride(5)=1
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=1
    count(5)=1
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF   
     
   CALL h5scopy_f(msp_nn_id,msp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(msp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)
         
!     describe the memory layout of the data again 
   dims(1)=2                  ! complex = 2 * real
   dims(2)=mynp12_read        ! read needs dimension of arrays
   dims(3)=mynp3_read
   dims(4)=mynq4_read
   dims(5)=mynalpha4N_read      ! here size of array in memory  

#ifdef HDF_MPIO
      ! check for consistency with MPI
   blocksize=4_lint*2*mynp12_read*mynp3_read*mynq4_read*mynalpha4N_read  
   IF(blocksize.GE.blocklimit) THEN
    WRITE(*,*) 'block too large to be handled with MPI-IO'
    WRITE(*,*) 'blocksize = ',blocksize
    WRITE(*,*) 'blocklimit = ',blocklimit
    CALL abort
   END IF
#endif
     
! read data 

   CALL h5dread_f(dset_nn_id, H5T_NATIVE_REAL,& 
     &      buf,dims,ierr,&
     &      mem_space_id = msp_hyper_id,&
     &      file_space_id = dsp_hyper_id,& 
     &      xfer_prp = pcollectread_id)
     
     
! here all data is read into buffers on myq4enerid=0      
     
! close hdf handles 
     ! close data set   
   CALL h5dclose_f(dset_nn_id,ierr)
  
! do not need mem_space anymore      
   CALL h5sclose_f(msp_hyper_id,ierr)
   CALL h5sclose_f(msp_nn_id,ierr)
! and not data space anymore    
   CALL h5sclose_f(dsp_hyper_id,ierr)
   CALL h5sclose_f(dsp_nn_id,ierr) 

! need data on all PEs in commener 
   CALL MPI_BCAST(buf,2*mynp12_read*mynp3_read*mynq4_read*mynalpha4N_read,&
    &            MPI_REAL4,0,commener,ierr)
 
! set up the relation between current channels and the channels 
! read in 
   ALLOCATE(alphap(mynalpha4N31))
   IF(cdep4N) THEN  ! use all channels 
    CALL set_alpha_4n31(qnalpha4N_read,alpha4Ncdepmax_read,alphap,alpha4Ncdepmax_read,cdep4N)
   ELSE ! only use the charge independent ones 
    CALL set_alpha_4n31(qnalpha4N_read,alpha4Ncdepmax_read,alphap,alpha4Nmax_read,cdep4N)
   END IF
  
  !  devide buffer by proper l factor for interpolation
  
  DO alphaloc=1,mynalpha4N_read   
    alpha=myalphaid+1+npe_alpha*(alphaloc-1)
    l12=qnalpha4N_read(1,alpha)
    l3=qnalpha4N_read(5,alpha)
    l4=qnalpha4N_read(9,alpha)
    
    DO iq4=1,mynq4_read
     DO ip3=1,mynp3_read
      DO ip12=1,mynp12_read   
       fakt=max(p12p_read(myp12_read-1+ip12)**l12*p3p_read(myp3_read-1+ip3)**l3*q4p_read(myq4_read-1+iq4)**l4,tolerance)   
       buf(1:2,ip12,ip3,iq4,alphaloc)=buf(1:2,ip12,ip3,iq4,alphaloc)/fakt
      END DO ! ip12 
     END DO  ! ip3 
    END DO   ! iq4 
  END DO   ! alphaloc  
   
! continue with interpolation of the data 
! collect all p12 mesh points to the processors 
   ALLOCATE(psi_tmp_p12(2,p12n_read,mynp3_read,mynq4_read,mynalpha4N_read))
   CALL collect_part_spreal(buf,psi_tmp_p12,2,p12n_read,&
        &                   mynalpha4N_read*mynp3_read*mynq4_read,commp12)
   DEALLOCATE(buf)
! then interpolate  
   ALLOCATE(psi_inter_p12(2,mynp12,mynp3_read,mynq4_read,mynalpha4N_read))

! now interpolation of psi/p12**l12/p3**l3/q4**l4
     
   psi_inter_p12=0.0
   DO alpha=1,mynalpha4N_read
    DO iq4=1,mynq4_read
     DO ip3=1,mynp3_read
      DO ip12=1,mynp12     
       psi_inter_p12(1:2,ip12,ip3,iq4,alpha) = &
     &    spl12(1,ip12)*psi_tmp_p12(1:2,indx12(1,ip12),ip3,iq4,alpha) &
     &  + spl12(2,ip12)*psi_tmp_p12(1:2,indx12(2,ip12),ip3,iq4,alpha) &
     &  + spl12(3,ip12)*psi_tmp_p12(1:2,indx12(3,ip12),ip3,iq4,alpha) &
     &  + spl12(4,ip12)*psi_tmp_p12(1:2,indx12(4,ip12),ip3,iq4,alpha)
      END DO
     END DO
    END DO
   END DO 
   DEALLOCATE(psi_tmp_p12) 

! continue with interpolation of the p3 data 
! collect all p3 mesh points to the processors 
   ALLOCATE(psi_tmp_p3(2,mynp12,p3n_read,mynq4_read,mynalpha4N_read))
   CALL collect_part_spreal(psi_inter_p12,psi_tmp_p3,&
        &                   2*mynp12,p3n_read,&
        &                   mynq4_read*mynalpha4N_read,commp3)
   DEALLOCATE(psi_inter_p12)
! then interpolate  
   ALLOCATE(psi_inter_p3(2,mynp12,mynp3,mynq4_read,mynalpha4N_read))

! now interpolation of psi/p12**l12/p3**l3/q4**l4
   
   psi_inter_p3=0.0
   DO alpha=1,mynalpha4N_read
    DO iq4=1,mynq4_read
     DO ip3=1,mynp3
      psi_inter_p3(1:2,1:mynp12,ip3,iq4,alpha) = &
     &    spl3(1,ip3)*psi_tmp_p3(1:2,1:mynp12,indx3(1,ip3),iq4,alpha) &
     &  + spl3(2,ip3)*psi_tmp_p3(1:2,1:mynp12,indx3(2,ip3),iq4,alpha) &
     &  + spl3(3,ip3)*psi_tmp_p3(1:2,1:mynp12,indx3(3,ip3),iq4,alpha) &
     &  + spl3(4,ip3)*psi_tmp_p3(1:2,1:mynp12,indx3(4,ip3),iq4,alpha)
     END DO
    END DO
   END DO 
   DEALLOCATE(psi_tmp_p3) 
  
! continue with interpolation of the q4 data 
! collect all q4 mesh points to the processors 
   ALLOCATE(psi_tmp_q4(2,mynp12,mynp3,q4n_read,mynalpha4N_read))
   CALL collect_part_spreal(psi_inter_p3,psi_tmp_q4,&
        &                   2*mynp12*mynp3,q4n_read,&
        &                   mynalpha4N_read,commq4)
   DEALLOCATE(psi_inter_p3)
! then interpolate  
   ALLOCATE(psi_inter_q4(2,mynp12,mynp3,mynq4,mynalpha4N_read))
 
   psi_inter_q4=0.0
   DO alpha=1,mynalpha4N_read
    DO iq4=1,mynq4
     psi_inter_q4(1:2,1:mynp12,1:mynp3,iq4,alpha) = &
     &    spl4(1,iq4)*psi_tmp_q4(1:2,1:mynp12,1:mynp3,indx4(1,iq4),alpha) &
     &  + spl4(2,iq4)*psi_tmp_q4(1:2,1:mynp12,1:mynp3,indx4(2,iq4),alpha) &
     &  + spl4(3,iq4)*psi_tmp_q4(1:2,1:mynp12,1:mynp3,indx4(3,iq4),alpha) &
     &  + spl4(4,iq4)*psi_tmp_q4(1:2,1:mynp12,1:mynp3,indx4(4,iq4),alpha)
    END DO
   END DO 
   DEALLOCATE(psi_tmp_q4) 
   
   
   
! now collect all channels to reorder the channels 
! to current bookkeeping, drop channels that cannot be assigned 
! leave components zero that are not found in file   
! start with collecting all channels  
  
   ALLOCATE(psi_tmp_alpha(2,mynp12,mynp3,mynq4,alpha4Ncdepmax_read))
   CALL collect_piece_spreal(psi_inter_q4,psi_tmp_alpha, &
  &                  2*mynp12*mynp3*mynq4,alpha4Ncdepmax_read,1,commalpha)
  
   DEALLOCATE(psi_inter_q4)
  
   IF(.NOT. allocated(PSI%amp)) &
    &   ALLOCATE(PSI%amp(mynp12,mynp3,mynq4,mynalpha4N31))

   PSI%amp=0.0  
   DO alphaloc=1,mynalpha4N31     
    IF(alphap(alphaloc).NE.0) THEN
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)   
     call get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari) 
     DO iq4=1,mynq4
      DO ip3=1,mynp3
       DO ip12=1,mynp12
        fakt=p12p(myp12-1+ip12)**l12*p3p(myp3-1+ip3)**l3*q4p(myq4-1+iq4)**l4   
        PSI%amp(ip12,ip3,iq4,alphaloc) &
     &  = fakt*psi_tmp_alpha(1,ip12,ip3,iq4,alphap(alphaloc)) &
         + fakt*IU*psi_tmp_alpha(2,ip12,ip3,iq4,alphap(alphaloc))
       END DO
      END DO
     END DO
    END IF
   END DO
  
   DEALLOCATE(psi_tmp_alpha)
   DEALLOCATE(alphap)
  END IF  ! 4N exists 
  
  hdftime=hdftime+MPI_WTIME()   
  
  DEALLOCATE(spl12,indx12)
  DEALLOCATE(spl3,indx3)
  DEALLOCATE(spl4,indx4)
  DEALLOCATE(p12p_read,p12w_read)
  DEALLOCATE(p3p_read,p3w_read)
  DEALLOCATE(q4p_read,q4w_read)
  DEALLOCATE(qnalpha4N_read)
 
! close the group for this amplitude 
  CALL h5gclose_f(group_amp_id,ierr) 
! finally close the file 
  CALL close_h5file(readf_id)  

  readamptime=readamptime+MPI_WTIME()
   
 END SUBROUTINE readamp_4n31

 SUBROUTINE set_alpha_4n31(qn,alphamaxcdep_read,alphap,alphamax_read,cdep)
  IMPLICIT NONE
  INTEGER alphamax_read,alphamaxcdep_read
  INTEGER qn(16,alphamaxcdep_read),alphap(mynalpha4N31)
  LOGICAL cdep
  INTEGER alphaneu,alphaloc,alphaalt
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari
  INTEGER l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p
  
  alphap=0
  
  DO alphaloc=1,mynalpha4N31
   alphaneu=myalphaid+1+npe_alpha*(alphaloc-1)
   call get4N31qn(alphaneu,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari)
   DO alphaalt=1,alphamax_read
    l12p =qn(1,alphaalt)
    s12p =qn(2,alphaalt)
    j12p =qn(3,alphaalt)
    t12p =qn(4,alphaalt)
    l3p =qn(5,alphaalt)
    I3p =qn(6,alphaalt)
    j3p =qn(7,alphaalt)
    tau3p =qn(8,alphaalt)
    l4p =qn(9,alphaalt)
    I4p =qn(10,alphaalt)
    j4p =qn(11,alphaalt)
    tau4p =qn(12,alphaalt)
    mtau4p =qn(13,alphaalt)
    

    IF((l12p.eq.l12).AND.&
    &  (s12p.eq.s12).AND.&
    &  (j12p.eq.j12).AND.&
    &  (t12p.eq.t12).AND.&
    &  (l3p.eq.l3).AND.&
    &  (I3p.eq.I3).AND.&
    &  (j3p.eq.j3).AND.&
    &  (tau3p.eq.tau3).AND.&
    &  (l4p.eq.l4).AND.&
    &  (I4p.eq.I4).AND.&
    &  (j4p.eq.j4).AND.&
    &  (tau4p.eq.tau4).AND.&
    &  ((mtau4p.eq.mtau4).OR. (.NOT. cdep))) THEN
     alphap(alphaloc)=alphaalt
     EXIT
    END IF
    IF(alphaalt.EQ.alphamax_read &
   &  .and.myenerid.EQ.0  &
   &  .and.mymeshid.EQ.0) THEN
     WRITE(*,*) '4N31 Channel not found:',alphaneu
    ENDIF
   END DO
  END DO ! alphaloc
  
 END SUBROUTINE set_alpha_4n31
 
 ! 4N22 amplitude
 ! open file <filename>
 ! write amplitude PSI to file <filename>
 ! in the group ampname
 ! file is then closed 
 
 SUBROUTINE writeamp_4n22(PSI,filename,ampname)
  IMPLICIT NONE 
  TYPE(AMP4N22),TARGET :: PSI 
  CHARACTER(LEN=*) filename
  CHARACTER(LEN=*) ampname
  INTEGER(HID_T) writef_id
  INTEGER(HID_T) group_amp_id
  INTEGER(HID_T) dset_nn_id,dsp_nn_id,dsp_hyper_id
  INTEGER(HID_T) msp_nn_id,msp_hyper_id
  INTEGER(HSIZE_T) dims(5)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
  INTEGER(HSIZE_T) start(5)     ! storage for offset in output arrays
  INTEGER(HSIZE_T) block(5)     ! storage for blocksize in output arrays
  INTEGER(HSIZE_T) stride(5)    ! storage for stride in output arrays
  INTEGER(HSIZE_T) count(5)     ! storage for count in output arrays
  INTEGER ierr
  REAL(spreal),ALLOCATABLE  :: buf(:,:,:,:,:)
  INTEGER(lint) blocksize

  writeamptime=writeamptime-MPI_WTIME()
  
  ! first open the hdf file
  ! is always assumed that the file exits
  ! copy empty file if necessary
  
  CALL open_write_h5file(filename,commall,writef_id)
  
  ! then create the group <ampname>

  CALL h5gcreate_f(writef_id,trim(ampname),group_amp_id, ierr)

  ! file should contain all info to interprete the data 
  !    1. meshpoints for P12 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  
  CALL writep12p(group_amp_id)
   
  !    2. meshpoints for P34 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  CALL writep34p(group_amp_id)

  !    3. meshpoints for Q grid including number of mesh points
  !       and mesh parameters as far as they are known 
  CALL writeqp(group_amp_id)
  
  !    4. 4N22 channels qnalpha4N22(:,:) including number of channels 
  CALL write_4n22_channels(group_amp_id)
  
  !    5. physical constants from constant mod are also written 
  CALL writephysconst(group_amp_id) 
  
  !    6. 4N31 amplitude if it exists 
  
  hdftime=hdftime-MPI_WTIME()
  
  IF(allocated(PSI%amp)) THEN
  ! this needs to be written in parallel since the data 
  ! is distributed 
  ! make sure that that only one processor writes the data 
  ! in case the distribution is not complete 
  
  ! here data is distributed in in commp12 and commalpha 
  ! that writing processors should have myp3id=myq4id=myenerid=0
  ! put the hyperslabs for other processors to empty 
  
  !     need to create a data set in the file 
  !     with complete number amp-elements 
  !     outline shape of global set of data 
   dims(1)=2                  !  complex = 2 * real
   dims(2)=P12N          !  write needs dimensions of arrays 
   dims(3)=P34N
   dims(4)=QN
   dims(5)=beta4N22cdepmax     !  here array size in data file  
  !     create corresponding data space            
   CALL h5screate_simple_f(5,dims,dsp_nn_id,ierr);

  !     describe the memory layout of the data
   dims(1)=2                  ! complex = 2 * real
   dims(2)=mynp12        ! write needs dimension of arrays
   dims(3)=mynp34
   dims(4)=mynq
   dims(5)=mynbeta4N22         ! here size of array in memory 
  !     create corresponding mem space    
   CALL h5screate_simple_f(5,dims,msp_nn_id,ierr);
  
  ! select the hyperslap in the file  that corresponds to the local 
  ! 4N31 data 
   IF(myenerid.EQ.0) THEN
    start(1)=0
    start(2)=myp12-1
    start(3)=myp34-1
    start(4)=myq-1
    start(5)=myalphaid 
    block(1)=2
    block(2)=mynp12
    block(3)=mynp34
    block(4)=mynq
    block(5)=1
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=1
    stride(5)=npe_alpha
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=1
    count(5)=mynbeta4N22
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF
  
   CALL h5scopy_f(dsp_nn_id,dsp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(dsp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)

  ! select the hyperslap in memorythat corresponds to the local 
  ! 4N31 data 

   IF(myenerid.EQ.0) THEN
    start(1)=0
    start(2)=0             
    start(3)=0
    start(4)=0
    start(5)=0
    block(1)=2
    block(2)=mynp12
    block(3)=mynp34
    block(4)=mynq
    block(5)=mynbeta4N22
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=1
    stride(5)=1
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=1
    count(5)=1
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF   
     
   CALL h5scopy_f(msp_nn_id,msp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(msp_hyper_id,H5S_SELECT_SET_F,&
      &                      start, count,ierr,stride,block)
       
      
!    create the dataset in file 
   CALL h5dcreate_f(group_amp_id,trim('4n22amp'),H5T_NATIVE_REAL,&
      &             dsp_nn_id,dset_nn_id,ierr)
    

      
! resolve compatiblity issue by copying to REAL buffer 
   ALLOCATE(buf(2,mynp12,mynp34,mynq,mynbeta4N22))
   buf(1,:,:,:,:)=REAL(PSI%amp(:,:,:,:))
   buf(2,:,:,:,:)=AIMAG(PSI%amp(:,:,:,:))
#ifdef HDF_MPIO
      ! check for consistency with MPI
   blocksize=4_lint*2*mynp12*mynp34*mynq*mynbeta4N22  
   IF(blocksize.GE.blocklimit) THEN
    WRITE(*,*) 'block too large to be handled with MPI-IO'
    WRITE(*,*) 'blocksize = ',blocksize
    WRITE(*,*) 'blocklimit = ',blocklimit
    CALL abort
   END IF
#endif
         
! write data of amplitude (use here collective communication)       
   CALL h5dwrite_f(dset_nn_id, H5T_NATIVE_REAL,& 
     &      buf,dims,ierr,&
     &      mem_space_id = msp_hyper_id,&
     &      file_space_id = dsp_hyper_id,& 
     &      xfer_prp = pcollectwrite_id)

   DEALLOCATE(buf)

! close data set   
   CALL h5dclose_f(dset_nn_id,ierr)
  
! do not need mem_space anymore      
   CALL h5sclose_f(msp_hyper_id,ierr)
   CALL h5sclose_f(msp_nn_id,ierr)
! and not data space anymore    
   CALL h5sclose_f(dsp_hyper_id,ierr)
   CALL h5sclose_f(dsp_nn_id,ierr) 
  END IF  ! 3N exists 
  
  hdftime=hdftime+MPI_WTIME()
  
! close the group for this amplitude 
   CALL h5gclose_f(group_amp_id,ierr) 
! finally close the file 
   CALL close_h5file(writef_id)  
  
  writeamptime=writeamptime+MPI_WTIME()
  
 END SUBROUTINE writeamp_4n22

 ! 4N22  amplitude
 ! open file <filename>
 ! read amplitude PSI from group ampname 
 ! of file <filename>
 ! file is then closed
 ! print the parameters and give back the amplitude
 
 SUBROUTINE readamp_4n22(PSI,filename,ampname)
  IMPLICIT NONE 
  TYPE(AMP4N22),TARGET :: PSI 
  CHARACTER(LEN=*) filename
  CHARACTER(LEN=*) ampname
  INTEGER(HID_T) readf_id
  INTEGER(HID_T) group_amp_id
  INTEGER(HID_T) dset_nn_id,dsp_nn_id,dsp_hyper_id
  INTEGER(HID_T) msp_nn_id,msp_hyper_id
  INTEGER(HSIZE_T) dims(5)      ! storage for dimensions of output arrays (maximal 4 dim arrays)
  INTEGER(HSIZE_T) start(5)     ! storage for offset in output arrays
  INTEGER(HSIZE_T) block(5)     ! storage for blocksize in output arrays
  INTEGER(HSIZE_T) stride(5)    ! storage for stride in output arrays
  INTEGER(HSIZE_T) count(5)     ! storage for count in output arrays
  INTEGER ierr,p12n_read,p34n_read,qn_read
  INTEGER alpha4Ncdepmax_read,alpha4Nmax_read
  INTEGER ip12,ip34,iq,alpha
  ! buffer to read the complex data 
  REAL(spreal),ALLOCATABLE  :: buf(:,:,:,:,:),&
     &                         psi_tmp_p12(:,:,:,:,:), &
     &                         psi_inter_p12(:,:,:,:,:),&
     &                         psi_tmp_p34(:,:,:,:,:), &
     &                         psi_inter_p34(:,:,:,:,:),&
     &                         psi_tmp_q(:,:,:,:,:), &
     &                         psi_inter_q(:,:,:,:,:),&     
     &                         psi_tmp_alpha(:,:,:,:,:)
  REAL(dpreal),ALLOCATABLE  :: p12p_read(:),p12w_read(:),&
                       &       p34p_read(:),p34w_read(:),&
                       &       qp_read(:),qw_read(:)
  INTEGER,ALLOCATABLE  :: qnalpha4N_read(:,:)
  INTEGER mynp12_read,myp12_read,myendp12_read
  INTEGER mynp34_read,myp34_read,myendp34_read
  INTEGER mynq_read,myq_read,myendq_read
  INTEGER mynalpha4N_read
  REAL(dpreal),ALLOCATABLE :: spl12(:,:),spl34(:,:),spl(:,:)
  INTEGER,ALLOCATABLE :: indx12(:,:),indx34(:,:),indx(:,:),alphap(:)
  LOGICAL amp_exists
  REAL(dpreal) :: fakt 
  INTEGER alphaloc,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari
  INTEGER(lint) blocksize

  readamptime=readamptime-MPI_WTIME()
  
  ! first open the hdf file
  ! is always assumed that the file exits
  ! copy empty file if necessary
  CALL open_read_h5file(filename,commall,readf_id)
  
  ! then create the group <ampname>
  
  CALL h5gopen_f(readf_id,trim(ampname),group_amp_id, ierr)

  ! file should contain all info to interprete the data 
  !    1. meshpoints for P12/P3/Q4 grid including number of mesh points
  !       and mesh parameters as far as they are known 
  !       the meshpoints and weights are allocated within the routine
  CALL readp12p(group_amp_id,p12p_read,p12w_read,p12n_read)
  CALL readp34p(group_amp_id,p34p_read,p34w_read,p34n_read)
  CALL readqp(group_amp_id,qp_read,qw_read,qn_read)
  
  !    2. 4N31 channels qn_read(:,:) including number of channels
  !       the quantum numbers are allocated within the routine 
  CALL read_4n22_channels(group_amp_id,qnalpha4N_read,&
     &                   alpha4Nmax_read,alpha4Ncdepmax_read)
  
  !    3. physical constants from constant mod are also read 
  !       and printed. No values are used. 
  CALL readphysconst(group_amp_id) 

! prepare the spline elements to interpolate
! from p12p_read to p12p 
! and p34p_read to p34p  
! and qp_read to qp 
  
  ALLOCATE(spl12(4,mynp12),indx12(4,mynp12))
  ALLOCATE(spl34(4,mynp34),indx34(4,mynp34))
  ALLOCATE(spl(4,mynq),indx(4,mynq))
  
  CALL cubfast_dp(p12p_read,p12n_read,&
        &    p12p(myp12:myp12+mynp12-1),mynp12, &
        &    spl12,indx12)
  CALL cubfast_dp(p34p_read,p34n_read,&
        &    p34p(myp34:myp34+mynp34-1),mynp34, &
        &    spl34,indx34)
  CALL cubfast_dp(qp_read,qn_read,&
        &    qp(myq:myq+mynq-1),mynq, &
        &    spl,indx)
        
  ! distribution of mesh and channels 
  CALL distr_block(p12n_read,npe_p12,myp12id,myp12_read,myendp12_read,mynp12_read)
  CALL distr_block(p34n_read,npe_p3,myp3id,myp34_read,myendp34_read,mynp34_read)
  CALL distr_block(qn_read,npe_q4,myq4id,myq_read,myendq_read,mynq_read)
  CALL distr_piece(alpha4Ncdepmax_read,npe_alpha,myalphaid,mynalpha4N_read)
  
  !    4. 4N amplitude if it exists 
  
  CALL h5lexists_f(group_amp_id, trim('4n22amp'),amp_exists,ierr)
  
  hdftime=hdftime-MPI_WTIME()
  
  IF(amp_exists) THEN
  ! first read the amplitude using the mesh stored 
  ! prepare buffer for real and imaginary part of the amplitude 
  
   ALLOCATE(buf(2,mynp12_read,mynp34_read,mynq_read,mynalpha4N_read))

  
  ! here data is distributed in in commp12 and commalpha 
  ! the reading processors should have myenerid=0
  ! this is element 0 in commener
  ! put the hyperslabs for other processors to empty 

  !     need to open the data set in the file 
   CALL h5dopen_f(group_amp_id,trim('4n22amp'),dset_nn_id,ierr)
  
  ! get the data space descriptor 
   CALL h5dget_space_f(dset_nn_id, dsp_nn_id, ierr)
     
  !     describe the memory layout of the data
   dims(1)=2                  ! complex = 2 * real
   dims(2)=mynp12_read        ! read needs dimension of arrays
   dims(3)=mynp34_read
   dims(4)=mynq_read
   dims(5)=mynalpha4N_read      ! here size of array in memory 
  !     create corresponding mem space    
   CALL h5screate_simple_f(5,dims,msp_nn_id,ierr);
  
  ! select the hyperslap in the file  that corresponds to the local 
  ! 4N data 
   IF(myenerid.EQ.0) THEN
    start(1)=0
    start(2)=myp12_read-1
    start(3)=myp34_read-1
    start(4)=myq_read-1
    start(5)=myalphaid 
    block(1)=2
    block(2)=mynp12_read
    block(3)=mynp34_read
    block(4)=mynq_read
    block(5)=1
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=1
    stride(5)=npe_alpha
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=1
    count(5)=mynalpha4N_read
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF
  
   CALL h5scopy_f(dsp_nn_id,dsp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(dsp_hyper_id,H5S_SELECT_SET_F,&
      &                      start, count,ierr,stride,block)

  ! select the hyperslap in memorythat corresponds to the local 
  ! 4N31 data 
   IF(myenerid.EQ.0) THEN
    start(1)=0
    start(2)=0             
    start(3)=0
    start(4)=0
    start(5)=0
    block(1)=2
    block(2)=mynp12_read
    block(3)=mynp34_read 
    block(4)=mynq_read 
    block(5)=mynalpha4N_read
    stride(1)=1
    stride(2)=1
    stride(3)=1
    stride(4)=1
    stride(5)=1
    count(1)=1
    count(2)=1
    count(3)=1
    count(4)=1
    count(5)=1
   ELSE
    start=0
    block=1
    stride=1
    count=0
   END IF   
     
   CALL h5scopy_f(msp_nn_id,msp_hyper_id,ierr)
   CALL h5sselect_hyperslab_f(msp_hyper_id,H5S_SELECT_SET_F,&
     &                      start, count,ierr,stride,block)
         
!     describe the memory layout of the data again 
   dims(1)=2                  ! complex = 2 * real
   dims(2)=mynp12_read        ! read needs dimension of arrays
   dims(3)=mynp34_read
   dims(4)=mynq_read
   dims(5)=mynalpha4N_read      ! here size of array in memory  

#ifdef HDF_MPIO
      ! check for consistency with MPI
   blocksize=4_lint*2*mynp12_read*mynp34_read*mynq_read*mynalpha4N_read  
   IF(blocksize.GE.blocklimit) THEN
    WRITE(*,*) 'block too large to be handled with MPI-IO'
    WRITE(*,*) 'blocksize = ',blocksize
    WRITE(*,*) 'blocklimit = ',blocklimit
    CALL abort
   END IF
#endif
     
! read data 

   CALL h5dread_f(dset_nn_id, H5T_NATIVE_REAL,& 
     &      buf,dims,ierr,&
     &      mem_space_id = msp_hyper_id,&
     &      file_space_id = dsp_hyper_id,& 
     &      xfer_prp = pcollectread_id)
     
     
! here all data is read into buffers on myq4enerid=0      
     
! close hdf handles 
     ! close data set   
   CALL h5dclose_f(dset_nn_id,ierr)
  
! do not need mem_space anymore      
   CALL h5sclose_f(msp_hyper_id,ierr)
   CALL h5sclose_f(msp_nn_id,ierr)
! and not data space anymore    
   CALL h5sclose_f(dsp_hyper_id,ierr)
   CALL h5sclose_f(dsp_nn_id,ierr) 

! need data on all PEs in commener 
   CALL MPI_BCAST(buf,2*mynp12_read*mynp34_read*mynq_read*mynalpha4N_read,&
    &            MPI_REAL4,0,commener,ierr)
 
! set up the relation between current channels and the channels 
! read in 
   ALLOCATE(alphap(mynbeta4N22))
   IF(cdep4N) THEN  ! use all channels 
    CALL set_alpha_4n22(qnalpha4N_read,alpha4Ncdepmax_read,alphap,alpha4Ncdepmax_read,cdep4N)
   ELSE ! only use the charge independent ones 
    CALL set_alpha_4n22(qnalpha4N_read,alpha4Ncdepmax_read,alphap,alpha4Nmax_read,cdep4N)
   END IF
  
   !  devide buffer by proper l factor for interpolation
  
   DO alphaloc=1,mynalpha4N_read   
    alpha=myalphaid+1+npe_alpha*(alphaloc-1)
    l12=qnalpha4N_read(1,alpha)
    l34=qnalpha4N_read(5,alpha)
    lam=qnalpha4N_read(9,alpha)
    
    DO iq=1,mynq_read
     DO ip34=1,mynp34_read
      DO ip12=1,mynp12_read   
       fakt=max(p12p_read(myp12_read-1+ip12)**l12*p34p_read(myp34_read-1+ip34)**l34*qp_read(myq_read-1+iq)**lam,tolerance)   
       buf(1:2,ip12,ip34,iq,alphaloc)=buf(1:2,ip12,ip34,iq,alphaloc)/fakt
      END DO ! ip12 
     END DO  ! ip3 
    END DO   ! iq4 
   END DO   ! alphaloc  
   
   
! continue with interpolation of the data 
! collect all p12 mesh points to the processors 
   ALLOCATE(psi_tmp_p12(2,p12n_read,mynp34_read,mynq_read,mynalpha4N_read))
   CALL collect_part_spreal(buf,psi_tmp_p12,2,p12n_read,&
        &                   mynalpha4N_read*mynp34_read*mynq_read,commp12)
   DEALLOCATE(buf)
! then interpolate  
   ALLOCATE(psi_inter_p12(2,mynp12,mynp34_read,mynq_read,mynalpha4N_read))

! now interpolation of psi/p12**l12/p34**l34/q**lam
   
   psi_inter_p12=0.0
   DO alpha=1,mynalpha4N_read
    DO iq=1,mynq_read
     DO ip34=1,mynp34_read
      DO ip12=1,mynp12     
       psi_inter_p12(1:2,ip12,ip34,iq,alpha) = &
     &    spl12(1,ip12)*psi_tmp_p12(1:2,indx12(1,ip12),ip34,iq,alpha) &
     &  + spl12(2,ip12)*psi_tmp_p12(1:2,indx12(2,ip12),ip34,iq,alpha) &
     &  + spl12(3,ip12)*psi_tmp_p12(1:2,indx12(3,ip12),ip34,iq,alpha) &
     &  + spl12(4,ip12)*psi_tmp_p12(1:2,indx12(4,ip12),ip34,iq,alpha)
      END DO
     END DO
    END DO
   END DO 
   DEALLOCATE(psi_tmp_p12) 
   
! continue with interpolation of the p3 data 
! collect all p3 mesh points to the processors 
   ALLOCATE(psi_tmp_p34(2,mynp12,p34n_read,mynq_read,mynalpha4N_read))
   CALL collect_part_spreal(psi_inter_p12,psi_tmp_p34,&
        &                   2*mynp12,p34n_read,&
        &                   mynq_read*mynalpha4N_read,commp3)
   DEALLOCATE(psi_inter_p12)
! then interpolate  
   ALLOCATE(psi_inter_p34(2,mynp12,mynp34,mynq_read,mynalpha4N_read))

! now interpolation of psi/p12**l12/p34**l34/q**lam
   
   psi_inter_p34=0.0
   DO alpha=1,mynalpha4N_read
    DO iq=1,mynq_read
     DO ip34=1,mynp34
      psi_inter_p34(1:2,1:mynp12,ip34,iq,alpha) = &
     &    spl34(1,ip34)*psi_tmp_p34(1:2,1:mynp12,indx34(1,ip34),iq,alpha) &
     &  + spl34(2,ip34)*psi_tmp_p34(1:2,1:mynp12,indx34(2,ip34),iq,alpha) &
     &  + spl34(3,ip34)*psi_tmp_p34(1:2,1:mynp12,indx34(3,ip34),iq,alpha) &
     &  + spl34(4,ip34)*psi_tmp_p34(1:2,1:mynp12,indx34(4,ip34),iq,alpha)
     END DO
    END DO
   END DO 
   DEALLOCATE(psi_tmp_p34) 

! continue with interpolation of the q4 data 
! collect all q4 mesh points to the processors 
   ALLOCATE(psi_tmp_q(2,mynp12,mynp34,qn_read,mynalpha4N_read))
   CALL collect_part_spreal(psi_inter_p34,psi_tmp_q,&
        &                   2*mynp12*mynp34,qn_read,&
        &                   mynalpha4N_read,commq4)
   DEALLOCATE(psi_inter_p34)
! then interpolate  
   ALLOCATE(psi_inter_q(2,mynp12,mynp34,mynq,mynalpha4N_read))

! now interpolation of psi/p12**l12/p34**l34/q**lam
    
   psi_inter_q=0.0
   DO alpha=1,mynalpha4N_read
    DO iq=1,mynq
     psi_inter_q(1:2,1:mynp12,1:mynp34,iq,alpha) = &
     &    spl(1,iq)*psi_tmp_q(1:2,1:mynp12,1:mynp34,indx(1,iq),alpha) &
     &  + spl(2,iq)*psi_tmp_q(1:2,1:mynp12,1:mynp34,indx(2,iq),alpha) &
     &  + spl(3,iq)*psi_tmp_q(1:2,1:mynp12,1:mynp34,indx(3,iq),alpha) &
     &  + spl(4,iq)*psi_tmp_q(1:2,1:mynp12,1:mynp34,indx(4,iq),alpha)
    END DO
   END DO 
   DEALLOCATE(psi_tmp_q) 
   
   
   
! now collect all channels to reorder the channels 
! to current bookkeeping, drop channels that cannot be assigned 
! leave components zero that are not found in file   
! start with collecting all channels  
  
   ALLOCATE(psi_tmp_alpha(2,mynp12,mynp34,mynq,alpha4Ncdepmax_read))
   CALL collect_piece_spreal(psi_inter_q,psi_tmp_alpha, &
  &                  2*mynp12*mynp34*mynq,alpha4Ncdepmax_read,1,commalpha)
  
   DEALLOCATE(psi_inter_q)
  
   IF(.NOT. allocated(PSI%amp)) &
    &   ALLOCATE(PSI%amp(mynp12,mynp34,mynq,mynbeta4N22))

   PSI%amp=0.0  
   DO alphaloc=1,mynbeta4N22     
    IF(alphap(alphaloc).NE.0) THEN
     alpha=myalphaid+1+npe_alpha*(alphaloc-1)   
     CALL get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari)
     DO iq=1,mynq
      DO ip34=1,mynp34
       DO ip12=1,mynp12
        fakt=p12p(myp12-1+ip12)**l12*p34p(myp34-1+ip34)**l34*qp(myq-1+iq)**lam   
        PSI%amp(ip12,ip34,iq,alphaloc) &
     &  = fakt*psi_tmp_alpha(1,ip12,ip34,iq,alphap(alphaloc)) &
         + fakt*IU*psi_tmp_alpha(2,ip12,ip34,iq,alphap(alphaloc))
       END DO
      END DO
     END DO
    END IF
   END DO 
  
   DEALLOCATE(psi_tmp_alpha)
   DEALLOCATE(alphap)
  END IF  ! 4N exists 
  
  hdftime=hdftime+MPI_WTIME()   
  
  DEALLOCATE(spl12,indx12)
  DEALLOCATE(spl34,indx34)
  DEALLOCATE(spl,indx)
  DEALLOCATE(p12p_read,p12w_read)
  DEALLOCATE(p34p_read,p34w_read)
  DEALLOCATE(qp_read,qw_read)
  DEALLOCATE(qnalpha4N_read)
 
! close the group for this amplitude 
  CALL h5gclose_f(group_amp_id,ierr) 
! finally close the file 
  CALL close_h5file(readf_id)  

  readamptime=readamptime+MPI_WTIME()
   
 END SUBROUTINE readamp_4n22

 SUBROUTINE set_alpha_4n22(qn,alphamaxcdep_read,alphap,alphamax_read,cdep)
  IMPLICIT NONE
  INTEGER alphamax_read,alphamaxcdep_read
  INTEGER qn(16,alphamaxcdep_read),alphap(mynbeta4N22)
  LOGICAL cdep
  INTEGER alphaneu,alphaloc,alphaalt
  INTEGER l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari
  INTEGER l12p,s12p,j12p,t12p,l34p,s34p,j34p,t34p,lamp,Ip,j4p,tau4p,mtau4p
  
  alphap=0
  
  DO alphaloc=1,mynbeta4N22
   alphaneu=myalphaid+1+npe_alpha*(alphaloc-1)
   call get4N22qn(alphaneu,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari)
   DO alphaalt=1,alphamax_read
    l12p =qn(1,alphaalt)
    s12p =qn(2,alphaalt)
    j12p =qn(3,alphaalt)
    t12p =qn(4,alphaalt)
    l34p =qn(5,alphaalt)
    s34p =qn(6,alphaalt)
    j34p =qn(7,alphaalt)
    t34p =qn(8,alphaalt)
    lamp =qn(9,alphaalt)
    Ip =qn(10,alphaalt)
    j4p =qn(11,alphaalt)
    tau4p =qn(12,alphaalt)
    mtau4p =qn(13,alphaalt)
    

    IF((l12p.eq.l12).AND.&
    &  (s12p.eq.s12).AND.&
    &  (j12p.eq.j12).AND.&
    &  (t12p.eq.t12).AND.&
    &  (l34p.eq.l34).AND.&
    &  (s34p.eq.s34).AND.&
    &  (j34p.eq.j34).AND.&
    &  (t34p.eq.t34).AND.&
    &  (lamp.eq.lam).AND.&
    &  (Ip.eq.I).AND.&
    &  (j4p.eq.j4).AND.&
    &  (tau4p.eq.tau4).AND.&
    &  ((mtau4p.eq.mtau4).OR. (.NOT. cdep))) THEN
     alphap(alphaloc)=alphaalt
     EXIT
    END IF
    IF(alphaalt.EQ.alphamax_read &
   &  .and.myenerid.EQ.0  &
   &  .and.mymeshid.EQ.0 ) THEN
     WRITE(*,*) '4N22 Channel not found:',alphaneu
    ENDIF
   END DO
  END DO ! alphaloc

 END SUBROUTINE set_alpha_4n22
  
 ! module procedures for assigment statement involving r-spade 
 ! amplitudes 
 
 ! first part: assigment of amplitude to amplitude 
 
 SUBROUTINE rset_psipsinn(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMPNN),INTENT(IN)    :: PSI2
  TYPE (RAMPNN),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynalphann))
  END IF
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  
  PSI1%amp=PSI2%amp   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psipsinn
 
 SUBROUTINE rset_psipsi3n(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP3N),INTENT(IN)    :: PSI2
  TYPE (RAMP3N),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr3,mynalpha3n))
  END IF
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  
  PSI1%amp=PSI2%amp   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psipsi3n
 
 SUBROUTINE rset_psipsi4n31(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP4N31),INTENT(IN)    :: PSI2
  TYPE (RAMP4N31),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr3,mynr4,mynalpha4n31))
  END IF
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  
  PSI1%amp=PSI2%amp   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psipsi4n31
 
 SUBROUTINE rset_psipsi4n22(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP4N22),INTENT(IN)    :: PSI2
  TYPE (RAMP4N22),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr34,mynr,mynbeta4n22))
  END IF
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  
  PSI1%amp=PSI2%amp   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psipsi4n22
 
 ! assignment of spcmplx to amplitude 
 
 SUBROUTINE rset_psispcmplxnn(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(spreal),INTENT(IN) :: lambda
  TYPE (RAMPNN),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynalphann))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psispcmplxnn
 
 SUBROUTINE rset_psispcmplx3n(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(spreal),INTENT(IN)    :: lambda
  TYPE (RAMP3N),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr3,mynalpha3n))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psispcmplx3n
 
 SUBROUTINE rset_psispcmplx4n31(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(spreal),INTENT(IN)    :: lambda
  TYPE (RAMP4N31),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr3,mynr4,mynalpha4n31))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psispcmplx4n31
 
 SUBROUTINE rset_psispcmplx4n22(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(spreal),INTENT(IN)    :: lambda
  TYPE (RAMP4N22),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr34,mynr,mynbeta4n22))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psispcmplx4n22
 
 
 ! assignment of spreal to amplitude 
 
 SUBROUTINE rset_psisprealnn(PSI1,lambda)
  IMPLICIT NONE
  REAL(spreal),INTENT(IN) :: lambda
  TYPE (RAMPNN),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynalphann))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psisprealnn
 
 SUBROUTINE rset_psispreal3n(PSI1,lambda)
  IMPLICIT NONE
  REAL(spreal),INTENT(IN)    :: lambda
  TYPE (RAMP3N),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr3,mynalpha3n))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psispreal3n
 
 SUBROUTINE rset_psispreal4n31(PSI1,lambda)
  IMPLICIT NONE
  REAL(spreal),INTENT(IN)    :: lambda
  TYPE (RAMP4N31),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr3,mynr4,mynalpha4n31))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psispreal4n31
 
 SUBROUTINE rset_psispreal4n22(PSI1,lambda)
  IMPLICIT NONE
  REAL(spreal),INTENT(IN)    :: lambda
  TYPE (RAMP4N22),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr34,mynr,mynbeta4n22))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psispreal4n22
 
 
 ! assignment of dpcmplx to amplitude 
 
 SUBROUTINE rset_psidpcmplxnn(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(dpreal),INTENT(IN) :: lambda
  TYPE (RAMPNN),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynalphann))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psidpcmplxnn
 
 SUBROUTINE rset_psidpcmplx3n(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  TYPE (RAMP3N),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr3,mynalpha3n))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psidpcmplx3n
 
 SUBROUTINE rset_psidpcmplx4n31(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  TYPE (RAMP4N31),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr3,mynr4,mynalpha4n31))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psidpcmplx4n31
 
 SUBROUTINE rset_psidpcmplx4n22(PSI1,lambda)
  IMPLICIT NONE
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  TYPE (RAMP4N22),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr34,mynr,mynbeta4n22))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psidpcmplx4n22
 
 
 ! assignment of dpreal to amplitude 
 
 SUBROUTINE rset_psidprealnn(PSI1,lambda)
  IMPLICIT NONE
  REAL(dpreal),INTENT(IN) :: lambda
  TYPE (RAMPNN),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynalphann))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psidprealnn
 
 SUBROUTINE rset_psidpreal3n(PSI1,lambda)
  IMPLICIT NONE
  REAL(dpreal),INTENT(IN)    :: lambda
  TYPE (RAMP3N),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr3,mynalpha3n))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psidpreal3n
 
 SUBROUTINE rset_psidpreal4n31(PSI1,lambda)
  IMPLICIT NONE
  REAL(dpreal),INTENT(IN)    :: lambda
  TYPE (RAMP4N31),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr3,mynr4,mynalpha4n31))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psidpreal4n31
 
 SUBROUTINE rset_psidpreal4n22(PSI1,lambda)
  IMPLICIT NONE
  REAL(dpreal),INTENT(IN)    :: lambda
  TYPE (RAMP4N22),INTENT(INOUT) :: PSI1
  
  ! only allocate if not allocated 
  
  IF(.NOT.allocated(PSI1%amp)) THEN
   ALLOCATE(PSI1%amp(mynr12,mynr34,mynr,mynbeta4n22))
  END IF
  
  PSI1%amp=lambda   ! copy amplitude element by element 
  
 END SUBROUTINE rset_psidpreal4n22


 ! now define functions for scalar products of r-space amplitudes 
 ! the scalar product results in COMPLEX(spreal)
 ! communication is required, therefore scalarproduct has to be 
 ! call by all PE's in commamp 
 
 ! 1. case:  nn amplitude 
 
 FUNCTION rdotprod_nn(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMPNN),INTENT(IN)    :: PSI2
  TYPE (RAMPNN),INTENT(IN)    :: PSI1
  COMPLEX(dpreal) :: rdotprod_nn
  COMPLEX(dpreal) :: locprod
  INTEGER ir12,alpha,ierr
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rdotprod !!'
  END IF
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rdotprod !!'
  END IF
#endif
  
  ! do local summation with r12 weights
  
  locprod=0.0_dpreal    
  DO alpha=1,mynalphaNN
   DO ir12=1,mynr12
    locprod=locprod &
         + conjg(PSI1%amp(ir12,alpha))*PSI2%amp(ir12,alpha) &
         *r12weight(ir12)
   END DO
  END DO
  
  ! sum local results using MPI 
  CALL MPI_ALLREDUCE(locprod,rdotprod_nn,1,MPI_COMPLEX16,MPI_SUM,commampnn,ierr) ! here it is assumed that product is COMPLEX(dpreal)
  
 END FUNCTION rdotprod_nn
 
 ! 2. case:  3n amplitude 
 
 FUNCTION rdotprod_3n(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP3N),INTENT(IN)    :: PSI2
  TYPE (RAMP3N),INTENT(IN)    :: PSI1
  COMPLEX(dpreal) :: rdotprod_3n
  COMPLEX(dpreal) :: locprod
  INTEGER ir12,ir3,alpha,ierr
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rdotprod !!'
  END IF
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rdotprod !!'
  END IF
#endif
  
  ! do local summation with p12,p3 weights
  
  locprod=0.0_dpreal
  DO alpha=1,mynalpha3N
   DO ir3=1,mynr3
    DO ir12=1,mynr12
     locprod=locprod &
          + conjg(PSI1%amp(ir12,ir3,alpha))*PSI2%amp(ir12,ir3,alpha) &
          *r12r3weight(ir12,ir3)
    END DO
   END DO
  END DO
  
  ! sum local results using MPI 
  CALL MPI_ALLREDUCE(locprod,rdotprod_3n,1,MPI_COMPLEX16,MPI_SUM,commamp3n,ierr) ! here it is assumed that product is COMPLEX(dpreal)
  
 END FUNCTION rdotprod_3n
 
 ! 3. case:  4n31 amplitude 
 
 FUNCTION rdotprod_4n31(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP4N31),INTENT(IN)    :: PSI2
  TYPE (RAMP4N31),INTENT(IN)    :: PSI1
  COMPLEX(dpreal) :: rdotprod_4n31
  COMPLEX(dpreal) :: locprod
  INTEGER ir12,ir3,ir4,alpha,ierr
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rdotprod !!'
  END IF
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rdotprod !!'
  END IF
#endif
  
  ! do local summation with p12,p3,q4 weights
  
  locprod=0.0_dpreal
  DO alpha=1,mynalpha4N31
   DO ir4=1,mynr4
    DO ir3=1,mynr3
     DO ir12=1,mynr12
      locprod=locprod &
           + conjg(PSI1%amp(ir12,ir3,ir4,alpha))*PSI2%amp(ir12,ir3,ir4,alpha) &
           *r12r3r4weight(ir12,ir3,ir4)
     END DO
    END DO
   END DO
  END DO
  
  ! sum local results using MPI 
  CALL MPI_ALLREDUCE(locprod,rdotprod_4n31,1,MPI_COMPLEX16,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is COMPLEX(dpreal)
  
 END FUNCTION rdotprod_4n31
 
 ! 4. case:  4n22 amplitude 
 
 FUNCTION rdotprod_4n22(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP4N22),INTENT(IN)    :: PSI2
  TYPE (RAMP4N22),INTENT(IN)    :: PSI1
  COMPLEX(dpreal) :: rdotprod_4n22
  COMPLEX(dpreal) :: locprod
  INTEGER ir12,ir34,ir,beta,ierr
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rdotprod !!'
  END IF
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rdotprod !!'
  END IF
#endif
  
  ! do local summation with r12,r34,r weights
  
  locprod=0.0_dpreal
  DO beta=1,mynbeta4N22
   DO ir=1,mynr
    DO ir34=1,mynr34
     DO ir12=1,mynr12
      locprod=locprod &
           + conjg(PSI1%amp(ir12,ir34,ir,beta))*PSI2%amp(ir12,ir34,ir,beta) &
           *r12r34rweight(ir12,ir34,ir)
     END DO
    END DO
   END DO
  END DO
  
  ! sum local results using MPI 
  CALL MPI_ALLREDUCE(locprod,rdotprod_4n22,1,MPI_COMPLEX16,MPI_SUM,commamp4n,ierr) ! here it is assumed that product is COMPLEX(dpreal)
  
 END FUNCTION rdotprod_4n22
 
 ! now define functions for products of r-space amplitudes
 ! with scalars 
 ! the product results in corresponding type of amplitude
 ! communication is not required 
 ! 1. case:  nn amplitude 
 ! 1a. scalar = REAL(spreal)
 FUNCTION rprodsprealnn(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMPNN) :: rprodsprealnn
  TYPE (RAMPNN),INTENT(IN)    :: PSI1
  REAL(spreal),INTENT(IN)    :: lambda
  INTEGER ir12,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rprodsprealnn !!'
  END IF
#endif
  rprodsprealnn=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalphaNN   ! set the amplitude 
   DO ir12=1,mynr12
    rprodsprealnn%amp(ir12,alpha)=lambda*PSI1%amp(ir12,alpha)
   END DO
  END DO
  
 END FUNCTION rprodsprealnn
 
 ! 1b. scalar = REAL(dpreal)
 FUNCTION rproddprealnn(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMPNN) :: rproddprealnn
  TYPE (RAMPNN),INTENT(IN)    :: PSI1
  REAL(dpreal),INTENT(IN)    :: lambda
  INTEGER ir12,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rproddprealnn !!'
  END IF
#endif
  rproddprealnn=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalphaNN   ! set the amplitude 
   DO ir12=1,mynr12
    rproddprealnn%amp(ir12,alpha)=lambda*PSI1%amp(ir12,alpha)
   END DO
  END DO
  
 END FUNCTION rproddprealnn
 
 ! 1c. scalar = COMPLEX(spreal)
 FUNCTION rprodspcmplxnn(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMPNN) :: rprodspcmplxnn
  TYPE (RAMPNN),INTENT(IN)    :: PSI1
  COMPLEX(spreal),INTENT(IN)    :: lambda
  INTEGER ir12,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rprodspcmplxnn !!'
  END IF
#endif
  rprodspcmplxnn=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalphaNN   ! set the amplitude 
   DO ir12=1,mynr12
    rprodspcmplxnn%amp(ir12,alpha)=lambda*PSI1%amp(ir12,alpha)
   END DO
  END DO
  
 END FUNCTION rprodspcmplxnn
 
 ! 1d. scalar = COMPLEX(dpreal)
 FUNCTION rproddpcmplxnn(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMPNN) :: rproddpcmplxnn
  TYPE (RAMPNN),INTENT(IN)    :: PSI1
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  INTEGER ir12,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rproddpcmplxnn !!'
  END IF
#endif
  rproddpcmplxnn=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalphaNN   ! set the amplitude 
   DO ir12=1,mynr12
    rproddpcmplxnn%amp(ir12,alpha)=lambda*PSI1%amp(ir12,alpha)
   END DO
  END DO
  
 END FUNCTION rproddpcmplxnn
 
 
 ! 2. case:  3n amplitude 
 ! 2a. scalar = REAL(spreal)
 FUNCTION rprodspreal3n(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP3N) :: rprodspreal3n
  TYPE (RAMP3N),INTENT(IN)    :: PSI1
  REAL(spreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir3,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rprodspreal3n !!'
  END IF
#endif
  rprodspreal3n=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha3N   ! set the amplitude 
   DO ir3=1,mynr3
    DO ir12=1,mynr12
     rprodspreal3n%amp(ir12,ir3,alpha)=lambda*PSI1%amp(ir12,ir3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION rprodspreal3n
 
 ! 2b. scalar = REAL(dpreal)
 FUNCTION rproddpreal3n(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP3N) :: rproddpreal3n
  TYPE (RAMP3N),INTENT(IN)    :: PSI1
  REAL(dpreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir3,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddpreal3n !!'
  END IF
#endif
  rproddpreal3n=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha3N   ! set the amplitude 
   DO ir3=1,mynr3
    DO ir12=1,mynr12
     rproddpreal3n%amp(ir12,ir3,alpha)=lambda*PSI1%amp(ir12,ir3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION rproddpreal3n
 
 ! 2c. scalar = COMPLEX(spreal)
 FUNCTION rprodspcmplx3n(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP3N) :: rprodspcmplx3n
  TYPE (RAMP3N),INTENT(IN)    :: PSI1
  COMPLEX(spreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir3,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rprodspcmplx3n !!'
  END IF
#endif
  rprodspcmplx3n=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha3N   ! set the amplitude 
   DO ir3=1,mynr3
    DO ir12=1,mynr12
     rprodspcmplx3n%amp(ir12,ir3,alpha)=lambda*PSI1%amp(ir12,ir3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION rprodspcmplx3n
 
 ! 2d. scalar = COMPLEX(dpreal)
 FUNCTION rproddpcmplx3n(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP3N) :: rproddpcmplx3n
  TYPE (RAMP3N),INTENT(IN)    :: PSI1
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir3,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddpcmplx3n !!'
  END IF
#endif
  rproddpcmplx3n=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha3N   ! set the amplitude 
   DO ir3=1,mynr3
    DO ir12=1,mynr12
     rproddpcmplx3n%amp(ir12,ir3,alpha)=lambda*PSI1%amp(ir12,ir3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION rproddpcmplx3n
 
 ! 3. case:  4n31 amplitude 
 ! 3a. scalar = REAL(spreal)
 FUNCTION rprodspreal4n31(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP4N31) :: rprodspreal4n31
  TYPE (RAMP4N31),INTENT(IN)    :: PSI1
  REAL(spreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir3,ir4,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rprodspreal4n31 !!'
  END IF
#endif
  rprodspreal4n31=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha4N31   ! set the amplitude 
   DO ir4=1,mynr4
    DO ir3=1,mynr3
     DO ir12=1,mynr12
      rprodspreal4n31%amp(ir12,ir3,ir4,alpha)=lambda*PSI1%amp(ir12,ir3,ir4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rprodspreal4n31
 
 ! 3b. scalar = REAL(dpreal)
 FUNCTION rproddpreal4n31(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP4N31) :: rproddpreal4n31
  TYPE (RAMP4N31),INTENT(IN)    :: PSI1
  REAL(dpreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir3,ir4,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in proddpreal4n31 !!'
  END IF
#endif
  rproddpreal4n31=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha4N31   ! set the amplitude 
   DO ir4=1,mynr4
    DO ir3=1,mynr3
     DO ir12=1,mynr12
      rproddpreal4n31%amp(ir12,ir3,ir4,alpha)=lambda*PSI1%amp(ir12,ir3,ir4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rproddpreal4n31
 
 ! 3c. scalar = COMPLEX(spreal)
 FUNCTION rprodspcmplx4n31(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP4N31) :: rprodspcmplx4n31
  TYPE (RAMP4N31),INTENT(IN)    :: PSI1
  COMPLEX(spreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir3,ir4,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rprodspcmplx4n31 !!'
  END IF
#endif
  rprodspcmplx4n31=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha4N31   ! set the amplitude 
   DO ir4=1,mynr4
    DO ir3=1,mynr3
     DO ir12=1,mynr12
      rprodspcmplx4n31%amp(ir12,ir3,ir4,alpha)=lambda*PSI1%amp(ir12,ir3,ir4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rprodspcmplx4n31
 
 ! 3d. scalar = COMPLEX(dpreal)
 FUNCTION rproddpcmplx4n31(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP4N31) :: rproddpcmplx4n31
  TYPE (RAMP4N31),INTENT(IN)    :: PSI1
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir3,ir4,alpha
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rproddpcmplx4n31 !!'
  END IF
#endif
  rproddpcmplx4n31=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO alpha=1,mynalpha4N31   ! set the amplitude 
   DO ir4=1,mynr4
    DO ir3=1,mynr3
     DO ir12=1,mynr12
      rproddpcmplx4n31%amp(ir12,ir3,ir4,alpha)=lambda*PSI1%amp(ir12,ir3,ir4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rproddpcmplx4n31
 
 ! 4. case:  4n22 amplitude 
 ! 4a. scalar = REAL(spreal)
 FUNCTION rprodspreal4n22(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP4N22) :: rprodspreal4n22
  TYPE (RAMP4N22),INTENT(IN)    :: PSI1
  REAL(spreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir34,ir,beta
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rprodspreal4n22 !!'
  END IF
#endif
  rprodspreal4n22=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO beta=1,mynbeta4N22   ! set the amplitude 
   DO ir=1,mynr
    DO ir34=1,mynr34
     DO ir12=1,mynr12
      rprodspreal4n22%amp(ir12,ir34,ir,beta)=lambda*PSI1%amp(ir12,ir34,ir,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rprodspreal4n22
 
 ! 4b. scalar = REAL(dpreal)
 FUNCTION rproddpreal4n22(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP4N22) :: rproddpreal4n22
  TYPE (RAMP4N22),INTENT(IN)    :: PSI1
  REAL(dpreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir34,ir,beta
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rproddpreal4n22 !!'
  END IF
#endif
  rproddpreal4n22=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO beta=1,mynbeta4N22   ! set the amplitude 
   DO ir=1,mynr
    DO ir34=1,mynr34
     DO ir12=1,mynr12
      rproddpreal4n22%amp(ir12,ir34,ir,beta)=lambda*PSI1%amp(ir12,ir34,ir,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rproddpreal4n22
 
 ! 4c. scalar = COMPLEX(spreal)
 FUNCTION rprodspcmplx4n22(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP4N22) :: rprodspcmplx4n22
  TYPE (RAMP4N22),INTENT(IN)    :: PSI1
  COMPLEX(spreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir34,ir,beta
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rprodspcmplx4n22 !!'
  END IF
#endif
  rprodspcmplx4n22=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO beta=1,mynbeta4N22   ! set the amplitude 
   DO ir=1,mynr
    DO ir34=1,mynr34
     DO ir12=1,mynr12
      rprodspcmplx4n22%amp(ir12,ir34,ir,beta)=lambda*PSI1%amp(ir12,ir34,ir,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rprodspcmplx4n22
 
 ! 4d. scalar = COMPLEX(dpreal)
 FUNCTION rproddpcmplx4n22(lambda,PSI1)
  IMPLICIT NONE
  TYPE (RAMP4N22) :: rproddpcmplx4n22
  TYPE (RAMP4N22),INTENT(IN)    :: PSI1
  COMPLEX(dpreal),INTENT(IN)    :: lambda
  INTEGER ir12,ir34,ir,beta
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rproddpcmplx4n22 !!'
  END IF
#endif
  rproddpcmplx4n22=0.0_spreal ! allocate memory and predefine amplitude 
  
  DO beta=1,mynbeta4N22   ! set the amplitude 
   DO ir=1,mynr
    DO ir34=1,mynr34
     DO ir12=1,mynr12
      rproddpcmplx4n22%amp(ir12,ir34,ir,beta)=lambda*PSI1%amp(ir12,ir34,ir,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rproddpcmplx4n22
 
 ! now the set of procedure for the OPERATOR(+) for r-space amplitudes
 ! 1. case nn amplitude 
 
 FUNCTION rplus_psinn(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMPNN),INTENT(IN)    :: PSI2
  TYPE (RAMPNN),INTENT(IN)    :: PSI1
  TYPE (RAMPNN)  :: rplus_psinn
  INTEGER alpha,ir12
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rplus_psinn !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rplus_psinn !!'
  END IF
#endif
  
  rplus_psinn=0.0_spreal   ! allocate & predefine psinn 
  
  DO alpha=1,mynalphaNN
   DO ir12=1,mynr12
    rplus_psinn%amp(ir12,alpha)=PSI1%amp(ir12,alpha)+PSI2%amp(ir12,alpha)
   END DO
  END DO
  
 END FUNCTION rplus_psinn
 
 ! 2. case 3n amplitude 
 
 FUNCTION rplus_psi3n(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP3N),INTENT(IN)    :: PSI2
  TYPE (RAMP3N),INTENT(IN)    :: PSI1
  TYPE (RAMP3N)   :: rplus_psi3n
  INTEGER alpha,ir12,ir3
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rplus_psi3n !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rplus_psi3n !!'
  END IF
#endif
  
  rplus_psi3n=0.0_spreal   ! allocate & predefine psi3n 
  
  DO alpha=1,mynalpha3N
   DO ir3=1,mynr3
    DO ir12=1,mynr12
     rplus_psi3n%amp(ir12,ir3,alpha)=PSI1%amp(ir12,ir3,alpha)+PSI2%amp(ir12,ir3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION rplus_psi3n
 
 ! 3. case 4n31 amplitude 
 
 FUNCTION rplus_psi4n31(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP4N31),INTENT(IN)    :: PSI2
  TYPE (RAMP4N31),INTENT(IN)    :: PSI1
  TYPE (RAMP4N31)   :: rplus_psi4n31
  INTEGER alpha,ir12,ir3,ir4
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rplus_psi4n31 !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rplus_psi4n31 !!'
  END IF
#endif
  
  rplus_psi4n31=0.0_spreal   ! allocate & predefine psi4n31 
  
  DO alpha=1,mynalpha4N31
   DO ir4=1,mynr4
    DO ir3=1,mynr3
     DO ir12=1,mynr12
      rplus_psi4n31%amp(ir12,ir3,ir4,alpha)=PSI1%amp(ir12,ir3,ir4,alpha)&
           +PSI2%amp(ir12,ir3,ir4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rplus_psi4n31
 
 ! 4. case 4n22 amplitude 
 
 FUNCTION rplus_psi4n22(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP4N22),INTENT(IN)    :: PSI2
  TYPE (RAMP4N22),INTENT(IN)    :: PSI1
  TYPE (RAMP4N22)   :: rplus_psi4n22
  INTEGER beta,ir12,ir34,ir
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rplus_psi4n22 !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rplus_psi4n22 !!'
  END IF
#endif
  
  rplus_psi4n22=0.0_spreal   ! allocate & predefine psi4n22 
  
  DO beta=1,mynbeta4N22
   DO ir=1,mynr
    DO ir34=1,mynr34
     DO ir12=1,mynr12
      rplus_psi4n22%amp(ir12,ir34,ir,beta)=PSI1%amp(ir12,ir34,ir,beta)&
           +PSI2%amp(ir12,ir34,ir,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rplus_psi4n22
 
 ! now the set of procedure for the OPERATOR(-) for r-space amplitudes
 ! 1. case nn amplitude 
 
 FUNCTION rminus_psinn(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMPNN),INTENT(IN)    :: PSI2
  TYPE (RAMPNN),INTENT(IN)    :: PSI1
  TYPE (RAMPNN)  :: rminus_psinn
  INTEGER alpha,ir12
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rminus_psinn !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rminus_psinn !!'
  END IF
#endif
  
  rminus_psinn=0.0_spreal   ! allocate & predefine psinn 
  
  DO alpha=1,mynalphaNN
   DO ir12=1,mynr12
    rminus_psinn%amp(ir12,alpha)=PSI1%amp(ir12,alpha)-PSI2%amp(ir12,alpha)
   END DO
  END DO
  
 END FUNCTION rminus_psinn
 
 ! 2. case 3n amplitude 
 
 FUNCTION rminus_psi3n(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP3N),INTENT(IN)    :: PSI2
  TYPE (RAMP3N),INTENT(IN)    :: PSI1
  TYPE (RAMP3N)   :: rminus_psi3n
  INTEGER alpha,ir12,ir3
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rminus_psi3n !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rminus_psi3n !!'
  END IF
#endif
  
  rminus_psi3n=0.0_spreal   ! allocate & predefine psi3n 
  
  DO alpha=1,mynalpha3N
   DO ir3=1,mynr3
    DO ir12=1,mynr12
     rminus_psi3n%amp(ir12,ir3,alpha)=PSI1%amp(ir12,ir3,alpha)-PSI2%amp(ir12,ir3,alpha)
    END DO
   END DO
  END DO
  
 END FUNCTION rminus_psi3n
 
 ! 3. case 4n31 amplitude 
 
 FUNCTION rminus_psi4n31(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP4N31),INTENT(IN)    :: PSI2
  TYPE (RAMP4N31),INTENT(IN)    :: PSI1
  TYPE (RAMP4N31)   :: rminus_psi4n31
  INTEGER alpha,ir12,ir3,ir4
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rminus_psi4n31 !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rminus_psi4n31 !!'
  END IF
#endif
  
  rminus_psi4n31=0.0_spreal   ! allocate & predefine psi4n31 
  
  DO alpha=1,mynalpha4N31
   DO ir4=1,mynr4
    DO ir3=1,mynr3
     DO ir12=1,mynr12
      rminus_psi4n31%amp(ir12,ir3,ir4,alpha)=PSI1%amp(ir12,ir3,ir4,alpha)&
           -PSI2%amp(ir12,ir3,ir4,alpha)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rminus_psi4n31
 
 ! 4. case 4n22 amplitude 
 
 FUNCTION rminus_psi4n22(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP4N22),INTENT(IN)    :: PSI2
  TYPE (RAMP4N22),INTENT(IN)    :: PSI1
  TYPE (RAMP4N22)   :: rminus_psi4n22
  INTEGER beta,ir12,ir34,ir
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI1%amp)) THEN
   STOP 'PSI1 not allocated in rminus_psi4n22 !!'
  END IF
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in rminus_psi4n22 !!'
  END IF
#endif
  
  rminus_psi4n22=0.0_spreal   ! allocate & predefine psi4n22 
  
  DO beta=1,mynbeta4N22
   DO ir=1,mynr
    DO ir34=1,mynr34
     DO ir12=1,mynr12
      rminus_psi4n22%amp(ir12,ir34,ir,beta)=PSI1%amp(ir12,ir34,ir,beta)&
           -PSI2%amp(ir12,ir34,ir,beta)
     END DO
    END DO
   END DO
  END DO
  
 END FUNCTION rminus_psi4n22

 
!! preparation of splines for Fourier trafo 
 
 ! subroutine initializes the splines from original r12,r3,r4,r34 and r
 ! grids to the dense rgrid and from the p12,p3,q4,p34 and q
 ! to the dense pgrid
 ! required for Fourier transformation 

 SUBROUTINE initspline
  IMPLICIT NONE 
  INTEGER ir12 

  bessspltime=bessspltime-MPI_WTIME()
  
  ! allocate arrays 
  ! r-space to p-space 
  ALLOCATE(r12spl(rintn,4),r12indx(rintn,4))
  ALLOCATE(r3spl(rintn,4),r3indx(rintn,4))
  ALLOCATE(r4spl(rintn,4),r4indx(rintn,4))
  ALLOCATE(r34spl(rintn,4),r34indx(rintn,4))
  ALLOCATE(rspl(rintn,4),rindx(rintn,4))
  ! p-space to r-space 
  ALLOCATE(p12spl(pintn,4),p12indx(pintn,4))
  ALLOCATE(p3spl(pintn,4),p3indx(pintn,4))
  ALLOCATE(q4spl(pintn,4),q4indx(pintn,4))
  ALLOCATE(p34spl(pintn,4),p34indx(pintn,4))
  ALLOCATE(qspl(pintn,4),qindx(pintn,4))
  
  ! set splines for r-space 
  CALL cubherm(r12p,r12n,rintp,rintn,r12spl(:,:),r12indx(:,:))
  CALL cubherm(r3p,r3n,rintp,rintn,r3spl(:,:),r3indx(:,:))
  CALL cubherm(r4p,r4n,rintp,rintn,r4spl(:,:),r4indx(:,:))
  CALL cubherm(r34p,r34n,rintp,rintn,r34spl(:,:),r34indx(:,:))
  CALL cubherm(rp,rn,rintp,rintn,rspl(:,:),rindx(:,:))
               
  ! set splines for p-space 
  CALL cubherm(p12p,p12n,pintp,pintn,p12spl(:,:),p12indx(:,:))
  CALL cubherm(p3p,p3n,pintp,pintn,p3spl(:,:),p3indx(:,:))
  CALL cubherm(q4p,q4n,pintp,pintn,q4spl(:,:),q4indx(:,:))
  CALL cubherm(p34p,p34n,pintp,pintn,p34spl(:,:),p34indx(:,:))
  CALL cubherm(qp,qn,pintp,pintn,qspl(:,:),qindx(:,:))             
               
  bessspltime=bessspltime+MPI_WTIME()

 END SUBROUTINE initspline
 
!! prepartion of bessel functions 
 
 SUBROUTINE initbessel
  IMPLICIT NONE
  LOGICAL,SAVE :: first=.true.          
  
  REAL(dpreal),ALLOCATABLE ::  BF(:,:,:)
  REAL(spreal),ALLOCATABLE ::  bes_l(:,:,:)
  REAL(dpreal) :: faktp,fakts,faktr

  INTEGER is,ir12,ip12,ir12int,ip12int
  INTEGER ir3,ip3,ir3int,ip3int
  INTEGER ir4,iq4,ir4int,iq4int
  INTEGER ir34,ip34,ir34int,ip34int
  INTEGER ir,iq,irint,iqint
  
  INTEGER l12,l3,l4,l34,lam,ierr
  INTEGER irintloc,mynr12int,mynr3int,mynr4int
  INTEGER ipintloc,mynp12int,mynp3int,mynq4int
  
  IF(.NOT. first) return 
  first=.false. 
  
  CALL initspline 

  ! distribute the interpolation/integration points 
  CALL distr_piece(rintn,npe_p3*npe_q4*npe_ener*npe_alpha,myp3q4alphaenerid,mynr12int)
  CALL distr_piece(pintn,npe_p3*npe_q4*npe_ener*npe_alpha,myp3q4alphaenerid,mynp12int)
  CALL distr_piece(rintn,npe_p12*npe_q4*npe_ener*npe_alpha,myp12q4alphaenerid,mynr3int)
  CALL distr_piece(pintn,npe_p12*npe_q4*npe_ener*npe_alpha,myp12q4alphaenerid,mynp3int)
  CALL distr_piece(rintn,npe_p12*npe_p3*npe_ener*npe_alpha,myp12p3alphaenerid,mynr4int)
  CALL distr_piece(pintn,npe_p12*npe_p3*npe_ener*npe_alpha,myp12p3alphaenerid,mynq4int)
  
  ! prepare first besselfunction for r to p-space: 
  !     besslr2p12(mynp12,r12n,0:l12max)
  !     besslr2p3(mynp3,r3n,0:l3max)
  !     besslr2p4(mynq4,r4n,0:l4max)
  !     besslr2p34(mynp34,r34n,0:l12max)
  !     besslr2p(mynq,rn,0:lammax)
  
  ! start with 12 coordinate 
  besstime=besstime-MPI_WTIME()
  ! temporary array for bessel functions depending on dense grid dimension
  ALLOCATE(BF(mynp12,0:l12max,mynr12int))
  BF=0.0
   
  DO irintloc=1,mynr12int          ! parallization only for ir12ind (p12bess will be identical on all nodes)
   ir12int=myp3q4alphaenerid+1+(irintloc-1)*npe_p3*npe_q4*npe_alpha*npe_ener
   DO l12=0,min(l12max,16) 
    CALL SBF1(l12,rintp(ir12int),mynp12,BF(:,l12,irintloc),P12P(myp12))
   END DO
  END DO
  besstime=besstime+MPI_WTIME()
  bessinttime=bessinttime-MPI_WTIME()
 
  ALLOCATE(besslr2p12(mynp12,r12n,0:l12max))
  ALLOCATE(bes_l(mynp12,r12n,0:l12max))
  
  ! perform integeral over bessel functions 
  ! for given l,p,r over rint 

  bes_l=0.0

  DO irintloc=1,mynr12int       ! parallization only for ir12ind
   ir12int=myp3q4alphaenerid+1+(irintloc-1)*npe_p3*npe_q4*npe_alpha*npe_ener
   faktr=rintp(ir12int)**2*rintw(ir12int)*sqrt(2.0/pi)  ! phase needs to be added while applied   
   DO is=1,4
    ir12=r12indx(ir12int,is)
    fakts=faktr*r12spl(ir12int,is) 
    DO l12=0,l12max    
     DO ip12=1,mynp12
      bes_l(ip12,ir12,l12) &
       &  =bes_l(ip12,ir12,l12) + fakts*BF(ip12,l12,irintloc)
      
     END DO ! ip12
    END DO ! l12 
   END DO ! is
  END DO ! ir12int 
  
  bessinttime=bessinttime+MPI_WTIME()
  bessredtime=bessredtime-MPI_WTIME() 
  
  CALL MPI_ALLREDUCE(bes_l,besslr2p12,r12n*mynp12*(l12max+1),&
      &     MPI_REAL4,MPI_SUM,commp3q4alphaener,ierr)
      
  bessredtime=bessredtime+MPI_WTIME() 
      
  DEALLOCATE(BF,bes_l)
  
! then  3 coordinate 
  besstime=besstime-MPI_WTIME()
  ! temporary array for bessel functions depending on dense grid dimension
  ALLOCATE(BF(mynp3,0:l3max,mynr3int))
  BF=0.0
   
  DO irintloc=1,mynr3int          ! parallization only for ir3ind
   ir3int=myp12q4alphaenerid+1+(irintloc-1)*npe_p12*npe_q4*npe_alpha*npe_ener
   DO l3=0,min(l3max,16) 
    CALL SBF1(l3,rintp(ir3int),mynp3,BF(:,l3,irintloc),P3P(myp3))
   END DO
  END DO
  besstime=besstime+MPI_WTIME()
  bessinttime=bessinttime-MPI_WTIME()
 
  ALLOCATE(besslr2p3(mynp3,r3n,0:l3max))
  ALLOCATE(bes_l(mynp3,r3n,0:l3max))
  
  ! perform integeral over bessel functions 
  ! for given l,p,r over rint 

  bes_l=0.0

  DO irintloc=1,mynr3int       ! parallization only for ir3ind
   ir3int=myp12q4alphaenerid+1+(irintloc-1)*npe_p12*npe_q4*npe_alpha*npe_ener
   faktr=rintp(ir3int)**2*rintw(ir3int)*sqrt(2.0/pi)  ! phase needs to be added while applied   
   DO is=1,4
    ir3=r3indx(ir3int,is)
    fakts=faktr*r3spl(ir3int,is) 
    DO l3=0,l3max    
     DO ip3=1,mynp3
      bes_l(ip3,ir3,l3) &
       &  =bes_l(ip3,ir3,l3) + fakts*BF(ip3,l3,irintloc)
      
     END DO ! ip3
    END DO ! l3
   END DO ! is
  END DO ! ir3int 
  
  bessinttime=bessinttime+MPI_WTIME()
  bessredtime=bessredtime-MPI_WTIME() 
  
  CALL MPI_ALLREDUCE(bes_l,besslr2p3,r3n*mynp3*(l3max+1),&
      &     MPI_REAL4,MPI_SUM,commp12q4alphaener,ierr)
      
  bessredtime=bessredtime+MPI_WTIME() 
      
  DEALLOCATE(BF,bes_l)
  
  
  
! then  4 coordinate 
  besstime=besstime-MPI_WTIME()
  ! temporary array for bessel functions depending on dense grid dimension
  ALLOCATE(BF(mynq4,0:l4max,mynr4int))
  BF=0.0
   
  DO irintloc=1,mynr4int          ! parallization only for ir3ind
   ir4int=myp12p3alphaenerid+1+(irintloc-1)*npe_p12*npe_p3*npe_alpha*npe_ener
   DO l4=0,min(l4max,16) 
    CALL SBF1(l4,rintp(ir4int),mynq4,BF(:,l4,irintloc),Q4P(myq4))
   END DO
  END DO
  besstime=besstime+MPI_WTIME()
  bessinttime=bessinttime-MPI_WTIME()
 
  ALLOCATE(besslr2p4(mynq4,r4n,0:l4max))
  ALLOCATE(bes_l(mynq4,r4n,0:l4max))
  
  ! perform integeral over bessel functions 
  ! for given l,p,r over rint 

  bes_l=0.0

  DO irintloc=1,mynr4int       ! parallization only for ir3ind
   ir4int=myp12p3alphaenerid+1+(irintloc-1)*npe_p12*npe_p3*npe_alpha*npe_ener
   faktr=rintp(ir4int)**2*rintw(ir4int)*sqrt(2.0/pi)  ! phase needs to be added while applied   
   DO is=1,4
    ir4=r4indx(ir4int,is)
    fakts=faktr*r4spl(ir4int,is) 
    DO l4=0,l4max    
     DO iq4=1,mynq4
      bes_l(iq4,ir4,l4) &
       &  =bes_l(iq4,ir4,l4) + fakts*BF(iq4,l4,irintloc)
      
     END DO ! iq4
    END DO ! l4
   END DO ! is
  END DO ! ir4int 
  
  bessinttime=bessinttime+MPI_WTIME()
  bessredtime=bessredtime-MPI_WTIME() 
  
  CALL MPI_ALLREDUCE(bes_l,besslr2p4,r4n*mynq4*(l4max+1),&
      &     MPI_REAL4,MPI_SUM,commp12p3alphaener,ierr)
      
  bessredtime=bessredtime+MPI_WTIME() 
      
  DEALLOCATE(BF,bes_l)

! then  34 coordinate 
  besstime=besstime-MPI_WTIME()
  ! temporary array for bessel functions depending on dense grid dimension
  ALLOCATE(BF(mynp34,0:l12max,mynr3int))
  BF=0.0
   
  DO irintloc=1,mynr3int          ! parallization only for ir3ind
   ir34int=myp12q4alphaenerid+1+(irintloc-1)*npe_p12*npe_q4*npe_alpha*npe_ener
   DO l34=0,min(l12max,16) 
    CALL SBF1(l34,rintp(ir34int),mynp34,BF(:,l34,irintloc),P34P(myp34))
   END DO
  END DO
  besstime=besstime+MPI_WTIME()
  bessinttime=bessinttime-MPI_WTIME()
 
  ALLOCATE(besslr2p34(mynp34,r34n,0:l12max))
  ALLOCATE(bes_l(mynp34,r34n,0:l12max))
  
  ! perform integeral over bessel functions 
  ! for given l,p,r over rint 

  bes_l=0.0

  DO irintloc=1,mynr3int       ! parallization only for ir3ind
   ir34int=myp12q4alphaenerid+1+(irintloc-1)*npe_p12*npe_q4*npe_alpha*npe_ener
   faktr=rintp(ir34int)**2*rintw(ir34int)*sqrt(2.0/pi)  ! phase needs to be added while applied   
   DO is=1,4
    ir34=r34indx(ir34int,is)
    fakts=faktr*r34spl(ir34int,is) 
    DO l34=0,l12max    
     DO ip34=1,mynp34
      bes_l(ip34,ir34,l34) &
       &  =bes_l(ip34,ir34,l34) + fakts*BF(ip34,l34,irintloc)
      
     END DO ! ip34
    END DO ! l34
   END DO ! is
  END DO ! ir34int 
  
  bessinttime=bessinttime+MPI_WTIME()
  bessredtime=bessredtime-MPI_WTIME() 
  
  CALL MPI_ALLREDUCE(bes_l,besslr2p34,r34n*mynp34*(l12max+1),&
      &     MPI_REAL4,MPI_SUM,commp12q4alphaener,ierr)
      
  bessredtime=bessredtime+MPI_WTIME() 
      
  DEALLOCATE(BF,bes_l)
  
  
! then  2+2 coordinate 
  besstime=besstime-MPI_WTIME()
  ! temporary array for bessel functions depending on dense grid dimension
  ALLOCATE(BF(mynq,0:lammax,mynr4int))
  BF=0.0
   
  DO irintloc=1,mynr4int          ! parallization only for ir4ind
   ir4int=myp12p3alphaenerid+1+(irintloc-1)*npe_p12*npe_p3*npe_alpha*npe_ener
   DO lam=0,min(lammax,16) 
    CALL SBF1(lam,rintp(ir4int),mynq,BF(:,lam,irintloc),QP(myq))
   END DO
  END DO
  besstime=besstime+MPI_WTIME()
  bessinttime=bessinttime-MPI_WTIME()
 
  ALLOCATE(besslr2p(mynq,rn,0:lammax))
  ALLOCATE(bes_l(mynq,rn,0:lammax))
  
  ! perform integeral over bessel functions 
  ! for given l,p,r over rint 

  bes_l=0.0

  DO irintloc=1,mynr4int       ! parallization only for ir3ind
   irint=myp12p3alphaenerid+1+(irintloc-1)*npe_p12*npe_p3*npe_alpha*npe_ener
   faktr=rintp(irint)**2*rintw(irint)*sqrt(2.0/pi)  ! phase needs to be added while applied   
   DO is=1,4
    ir=rindx(irint,is)
    fakts=faktr*rspl(irint,is) 
    DO lam=0,lammax    
     DO iq=1,mynq
      bes_l(iq,ir,lam) &
       &  =bes_l(iq,ir,lam) + fakts*BF(iq,lam,irintloc)
      
     END DO ! iq
    END DO ! lam
   END DO ! is
  END DO ! irint 
  
  bessinttime=bessinttime+MPI_WTIME()
  bessredtime=bessredtime-MPI_WTIME() 
  
  CALL MPI_ALLREDUCE(bes_l,besslr2p,rn*mynq*(lammax+1),&
      &     MPI_REAL4,MPI_SUM,commp12p3alphaener,ierr)
      
  bessredtime=bessredtime+MPI_WTIME() 
      
  DEALLOCATE(BF,bes_l)

  DEALLOCATE(r12spl,r3spl,r4spl,r34spl,rspl)
  DEALLOCATE(r12indx,r3indx,r4indx,r34indx,rindx)
  
  
  ! prepare first besselfunction for p to r-space: 
  !     besslp2r12(mynr12,p12n,0:l12max)
  !     besslp2r3(mynr3,p3n,0:l3max)
  !     besslp2r4(mynr4,q4n,0:l4max)
  !     besslp2r34(mynr34,p34n,0:l12max)
  !     besslp2r(mynr,qn,0:lammax)
  
  ! start with 12 coordinate 
  besstime=besstime-MPI_WTIME()
  ! temporary array for bessel functions depending on dense grid dimension
  ALLOCATE(BF(mynr12,0:l12max,mynp12int))
  BF=0.0
   
  DO ipintloc=1,mynp12int          ! parallization only for ip12ind 
   ip12int=myp3q4alphaenerid+1+(ipintloc-1)*npe_p3*npe_q4*npe_alpha*npe_ener
   DO l12=0,min(l12max,16) 
    CALL SBF1(l12,pintp(ip12int),mynr12,BF(:,l12,ipintloc),R12P(myr12))
   END DO
  END DO
  besstime=besstime+MPI_WTIME()
  bessinttime=bessinttime-MPI_WTIME()
 
  ALLOCATE(besslp2r12(mynr12,p12n,0:l12max))
  ALLOCATE(bes_l(mynr12,p12n,0:l12max))
  
  ! perform integeral over bessel functions 
  ! for given l,p,r over pint 

  bes_l=0.0

  DO ipintloc=1,mynp12int       ! parallization only for ip12ind
   ip12int=myp3q4alphaenerid+1+(ipintloc-1)*npe_p3*npe_q4*npe_alpha*npe_ener
   faktp=pintp(ip12int)**2*pintw(ip12int)*sqrt(2.0/pi)  ! phase needs to be added while applied   
   DO is=1,4
    ip12=p12indx(ip12int,is)
    fakts=faktp*p12spl(ip12int,is) 
    DO l12=0,l12max    
     DO ir12=1,mynr12
      bes_l(ir12,ip12,l12) &
       &  =bes_l(ir12,ip12,l12) + fakts*BF(ir12,l12,ipintloc)
      
     END DO ! ip12
    END DO ! l12 
   END DO ! is
  END DO ! ip12int 
  
  bessinttime=bessinttime+MPI_WTIME()
  bessredtime=bessredtime-MPI_WTIME() 
  
  CALL MPI_ALLREDUCE(bes_l,besslp2r12,p12n*mynr12*(l12max+1),&
      &     MPI_REAL4,MPI_SUM,commp3q4alphaener,ierr)
      
  bessredtime=bessredtime+MPI_WTIME() 
      
  DEALLOCATE(BF,bes_l)
  
! then  3 coordinate 
  besstime=besstime-MPI_WTIME()
  ! temporary array for bessel functions depending on dense grid dimension
  ALLOCATE(BF(mynr3,0:l3max,mynp3int))
  BF=0.0
   
  DO ipintloc=1,mynr3int          ! parallization only for ip3ind
   ip3int=myp12q4alphaenerid+1+(ipintloc-1)*npe_p12*npe_q4*npe_alpha*npe_ener
   DO l3=0,min(l3max,16) 
    CALL SBF1(l3,pintp(ip3int),mynr3,BF(:,l3,ipintloc),R3P(myr3))
   END DO
  END DO
  besstime=besstime+MPI_WTIME()
  bessinttime=bessinttime-MPI_WTIME()
 
  ALLOCATE(besslp2r3(mynr3,p3n,0:l3max))
  ALLOCATE(bes_l(mynr3,p3n,0:l3max))
  
  ! perform integeral over bessel functions 
  ! for given l,p,r over pint 

  bes_l=0.0

  DO ipintloc=1,mynp3int       ! parallization only for ip3ind
   ip3int=myp12q4alphaenerid+1+(ipintloc-1)*npe_p12*npe_q4*npe_alpha*npe_ener
   faktp=pintp(ip3int)**2*pintw(ip3int)*sqrt(2.0/pi)  ! phase needs to be added while applied   
   DO is=1,4
    ip3=p3indx(ip3int,is)
    fakts=faktp*p3spl(ip3int,is) 
    DO l3=0,l3max    
     DO ir3=1,mynr3
      bes_l(ir3,ip3,l3) &
       &  =bes_l(ir3,ip3,l3) + fakts*BF(ir3,l3,ipintloc)
      
     END DO ! ip3
    END DO ! l3
   END DO ! is
  END DO ! ip3int 
  
  bessinttime=bessinttime+MPI_WTIME()
  bessredtime=bessredtime-MPI_WTIME() 
  
  CALL MPI_ALLREDUCE(bes_l,besslp2r3,p3n*mynr3*(l3max+1),&
      &     MPI_REAL4,MPI_SUM,commp12q4alphaener,ierr)
      
  bessredtime=bessredtime+MPI_WTIME() 
      
  DEALLOCATE(BF,bes_l)
  
  
  
! then  4 coordinate 
  besstime=besstime-MPI_WTIME()
  ! temporary array for bessel functions depending on dense grid dimension
  ALLOCATE(BF(mynr4,0:l4max,mynq4int))
  BF=0.0
   
  DO ipintloc=1,mynq4int          ! parallization only for ir3ind
   iq4int=myp12p3alphaenerid+1+(ipintloc-1)*npe_p12*npe_p3*npe_alpha*npe_ener
   DO l4=0,min(l4max,16) 
    CALL SBF1(l4,pintp(iq4int),mynr4,BF(:,l4,ipintloc),R4P(myr4))
   END DO
  END DO
  besstime=besstime+MPI_WTIME()
  bessinttime=bessinttime-MPI_WTIME()
 
  ALLOCATE(besslp2r4(mynr4,q4n,0:l4max))
  ALLOCATE(bes_l(mynr4,q4n,0:l4max))
  
  ! perform integeral over bessel functions 
  ! for given l,p,r over rint 

  bes_l=0.0

  DO ipintloc=1,mynq4int       ! parallization only for iq4ind
   iq4int=myp12p3alphaenerid+1+(ipintloc-1)*npe_p12*npe_p3*npe_alpha*npe_ener
   faktp=pintp(iq4int)**2*pintw(iq4int)*sqrt(2.0/pi)  ! phase needs to be added while applied   
   DO is=1,4
    iq4=q4indx(iq4int,is)
    fakts=faktp*q4spl(iq4int,is) 
    DO l4=0,l4max    
     DO ir4=1,mynr4
      bes_l(ir4,iq4,l4) &
       &  =bes_l(ir4,iq4,l4) + fakts*BF(ir4,l4,ipintloc)
      
     END DO ! iq4
    END DO ! l4
   END DO ! is
  END DO ! iq4int 
  
  bessinttime=bessinttime+MPI_WTIME()
  bessredtime=bessredtime-MPI_WTIME() 
  
  CALL MPI_ALLREDUCE(bes_l,besslp2r4,q4n*mynr4*(l4max+1),&
      &     MPI_REAL4,MPI_SUM,commp12p3alphaener,ierr)
      
  bessredtime=bessredtime+MPI_WTIME() 
      
  DEALLOCATE(BF,bes_l)

! then  34 coordinate 
  besstime=besstime-MPI_WTIME()
  ! temporary array for bessel functions depending on dense grid dimension
  ALLOCATE(BF(mynr34,0:l12max,mynp3int))
  BF=0.0
   
  DO ipintloc=1,mynp3int          ! parallization only for ip34ind
   ip34int=myp12q4alphaenerid+1+(ipintloc-1)*npe_p12*npe_q4*npe_alpha*npe_ener
   DO l34=0,min(l12max,16) 
    CALL SBF1(l34,pintp(ip34int),mynr34,BF(:,l34,ipintloc),R34P(myr34))
   END DO
  END DO
  besstime=besstime+MPI_WTIME()
  bessinttime=bessinttime-MPI_WTIME()
 
  ALLOCATE(besslp2r34(mynr34,p34n,0:l12max))
  ALLOCATE(bes_l(mynr34,p34n,0:l12max))
  
  ! perform integeral over bessel functions 
  ! for given l,p,r over rint 

  bes_l=0.0

  DO ipintloc=1,mynp3int       ! parallization only for irp34ind
   ip34int=myp12q4alphaenerid+1+(ipintloc-1)*npe_p12*npe_q4*npe_alpha*npe_ener
   faktp=pintp(ip34int)**2*pintw(ip34int)*sqrt(2.0/pi)  ! phase needs to be added while applied   
   DO is=1,4
    ip34=p34indx(ip34int,is)
    fakts=faktp*p34spl(ip34int,is) 
    DO l34=0,l12max    
     DO ir34=1,mynr34
      bes_l(ir34,ip34,l34) &
       &  =bes_l(ir34,ip34,l34) + fakts*BF(ir34,l34,ipintloc)
      
     END DO ! ip34
    END DO ! l34
   END DO ! is
  END DO ! ip34int 
  
  bessinttime=bessinttime+MPI_WTIME()
  bessredtime=bessredtime-MPI_WTIME() 
  
  CALL MPI_ALLREDUCE(bes_l,besslp2r34,p34n*mynr34*(l12max+1),&
      &     MPI_REAL4,MPI_SUM,commp12q4alphaener,ierr)
      
  bessredtime=bessredtime+MPI_WTIME() 
      
  DEALLOCATE(BF,bes_l)
  
  
! then  2+2 coordinate 
  besstime=besstime-MPI_WTIME()
  ! temporary array for bessel functions depending on dense grid dimension
  ALLOCATE(BF(mynr,0:lammax,mynq4int))
  BF=0.0
   
  DO ipintloc=1,mynq4int          ! parallization only for iqind
   iqint=myp12p3alphaenerid+1+(ipintloc-1)*npe_p12*npe_p3*npe_alpha*npe_ener
   DO lam=0,min(lammax,16) 
    CALL SBF1(lam,pintp(iqint),mynr,BF(:,lam,ipintloc),RP(myr))
   END DO
  END DO
  besstime=besstime+MPI_WTIME()
  bessinttime=bessinttime-MPI_WTIME()
 
  ALLOCATE(besslp2r(mynr,qn,0:lammax))
  ALLOCATE(bes_l(mynr,qn,0:lammax))
  
  ! perform integeral over bessel functions 
  ! for given l,p,r over pint 

  bes_l=0.0

  DO ipintloc=1,mynq4int       ! parallization only for iqind
   iqint=myp12p3alphaenerid+1+(ipintloc-1)*npe_p12*npe_p3*npe_alpha*npe_ener
   faktp=pintp(iqint)**2*pintw(iqint)*sqrt(2.0/pi)  ! phase needs to be added while applied   
   DO is=1,4
    iq=qindx(iqint,is)
    fakts=faktp*qspl(iqint,is) 
    DO lam=0,lammax    
     DO ir=1,mynr
      bes_l(ir,iq,lam) &
       &  =bes_l(ir,iq,lam) + fakts*BF(ir,lam,ipintloc)
      
     END DO ! iq
    END DO ! lam
   END DO ! is
  END DO ! iqint 
  
  bessinttime=bessinttime+MPI_WTIME()
  bessredtime=bessredtime-MPI_WTIME() 
  
  CALL MPI_ALLREDUCE(bes_l,besslp2r,qn*mynr*(lammax+1),&
      &     MPI_REAL4,MPI_SUM,commp12p3alphaener,ierr)
      
  bessredtime=bessredtime+MPI_WTIME() 
      
  DEALLOCATE(BF,bes_l)

  DEALLOCATE(p12spl,p3spl,q4spl,p34spl,qspl)
  DEALLOCATE(p12indx,p3indx,q4indx,p34indx,qindx)
 
     
 END SUBROUTINE initbessel 
 
 
!! now add Fourier trafos from p-space to r-space "p2r"
!! 1. NN amplitudes  
 
 SUBROUTINE fourier_p2r_psinn(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMPNN),INTENT(IN)     :: PSI2
  TYPE (RAMPNN),INTENT(INOUT) :: PSI1
  INTEGER ip12,alpha,alphaloc
  INTEGER l12,s12,j12,t12,mt12
  COMPLEX(dpreal) :: phase 
  COMPLEX(spreal),ALLOCATABLE :: psicollect12(:,:)
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  ! init besselfunctions if not done yet 
  CALL initbessel
  
  ! set zero and allocate if necessary 
  PSI1=0.0
    
  ALLOCATE(psicollect12(P12N,mynalphann))
  CALL collect_part_spcmplx(PSI2%amp,psicollect12,1,P12N,mynalphann,commp12)
  
  DO alphaloc=1,mynalphann
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL getnnqn(alpha,l12,s12,j12,t12,mt12)
   phase=iu**l12    ! should be correct for p- to r-space 
   DO ip12=1,P12N
    PSI1%amp(:,alphaloc)=PSI1%amp(:,alphaloc)+phase*psicollect12(ip12,alphaloc)*besslp2r12(:,ip12,l12)
   END DO ! ip12
  END DO  ! alpha
  
  DEALLOCATE(psicollect12)
  
 !! 2. 3N amplitudes  
 END SUBROUTINE fourier_p2r_psinn
 
 SUBROUTINE fourier_p2r_psi3n(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP3N),INTENT(IN)     :: PSI2
  TYPE (RAMP3N),INTENT(INOUT) :: PSI1
  INTEGER ip12,ip3,alpha,alphaloc,ir12
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  COMPLEX(dpreal) :: phase 
  COMPLEX(spreal),ALLOCATABLE :: psicollect12(:,:,:),psicollect3(:,:,:),psitrans12(:,:,:)
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  ! init besselfunctions if not done yet 
  CALL initbessel
  
  ! set zero and allocate if necessary 
  PSI1=0.0
    
  ALLOCATE(psicollect12(P12N,mynp3,mynalpha3N))
  CALL collect_part_spcmplx(PSI2%amp,psicollect12,1,P12N,mynp3*mynalpha3N,commp12)
  
  ALLOCATE(psitrans12(mynr12,mynp3,mynalpha3N))
  psitrans12=0.0
  
  DO alphaloc=1,mynalpha3N
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)  
   phase=iu**l12    ! should be correct for p- to r-space 
   DO ip3=1,mynp3
    DO ip12=1,P12N
     psitrans12(:,ip3,alphaloc)=psitrans12(:,ip3,alphaloc)+phase*psicollect12(ip12,ip3,alphaloc)*besslp2r12(:,ip12,l12)
    END DO ! ip12
   END DO  ! ip3 
  END DO  ! alpha

  DEALLOCATE(psicollect12)
  ALLOCATE(psicollect3(mynr12,P3N,mynalpha3N))
  CALL collect_part_spcmplx(psitrans12,psicollect3,mynr12,P3N,mynalpha3N,commp3)
  DEALLOCATE(psitrans12)
  
  
  DO alphaloc=1,mynalpha3N
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)  
   phase=iu**l3    ! should be correct for p- to r-space 
   DO ip3=1,P3N
    DO ir12=1,mynr12
     PSI1%amp(ir12,:,alphaloc)=PSI1%amp(ir12,:,alphaloc)+phase*psicollect3(ir12,ip3,alphaloc)*besslp2r3(:,ip3,l3)
    END DO ! ir12
   END DO  ! ip3 
  END DO  ! alpha
  
  DEALLOCATE(psicollect3)
  
 END SUBROUTINE fourier_p2r_psi3n

 !! 3. 4N31 amplitudes 
 SUBROUTINE fourier_p2r_psi4n31(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP4N31),INTENT(IN)     :: PSI2
  TYPE (RAMP4N31),INTENT(INOUT) :: PSI1
  INTEGER ip12,ip3,iq4,alpha,alphaloc,ir12,ir3
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari
  COMPLEX(dpreal) :: phase 
  COMPLEX(spreal),ALLOCATABLE :: psicollect12(:,:,:,:),psicollect3(:,:,:,:),psicollect4(:,:,:,:),&
         &     psitrans12(:,:,:,:),psitrans3(:,:,:,:)
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  ! init besselfunctions if not done yet 
  CALL initbessel
  
  ! set zero and allocate if necessary 
  PSI1=0.0
    
  ALLOCATE(psicollect12(P12N,mynp3,mynq4,mynalpha4N31))
  CALL collect_part_spcmplx(PSI2%amp,psicollect12,1,P12N,mynp3*mynq4*mynalpha4N31,commp12)
  
  ALLOCATE(psitrans12(mynr12,mynp3,mynq4,mynalpha4N31))
  psitrans12=0.0
  
  DO alphaloc=1,mynalpha4N31
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari)  
   phase=iu**l12    ! should be correct for p- to r-space 
   DO iq4=1,mynq4
    DO ip3=1,mynp3
     DO ip12=1,P12N
      psitrans12(:,ip3,iq4,alphaloc)=psitrans12(:,ip3,iq4,alphaloc)+phase*psicollect12(ip12,ip3,iq4,alphaloc)*besslp2r12(:,ip12,l12)
     END DO ! ip12
    END DO  ! ip3 
   END DO   ! iq4
  END DO  ! alpha

  DEALLOCATE(psicollect12)
  ALLOCATE(psicollect3(mynr12,P3N,mynq4,mynalpha4N31))
  CALL collect_part_spcmplx(psitrans12,psicollect3,mynr12,P3N,mynq4*mynalpha4N31,commp3)
  DEALLOCATE(psitrans12)
  
  ALLOCATE(psitrans3(mynr12,mynr3,mynq4,mynalpha4N31))
  psitrans3=0.0
  
  DO alphaloc=1,mynalpha4N31
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari)    
   phase=iu**l3    ! should be correct for p- to r-space 
   DO iq4=1,mynq4
    DO ip3=1,P3N
     DO ir12=1,mynr12
      psitrans3(ir12,:,iq4,alphaloc)=psitrans3(ir12,:,iq4,alphaloc)+phase*psicollect3(ir12,ip3,iq4,alphaloc)*besslp2r3(:,ip3,l3)
     END DO ! ir12
    END DO  ! ip3 
   END DO   ! iq4
  END DO  ! alpha
  DEALLOCATE(psicollect3)
  
  ALLOCATE(psicollect4(mynr12,mynr3,Q4N,mynalpha4N31))
  CALL collect_part_spcmplx(psitrans3,psicollect4,mynr12*mynr3,Q4N,mynalpha4N31,commq4)
  DEALLOCATE(psitrans3)
  
  
  DO alphaloc=1,mynalpha4N31
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari)    
   phase=iu**l4    ! should be correct for p- to r-space 
   DO iq4=1,Q4N
    DO ir3=1,mynr3
     DO ir12=1,mynr12
      PSI1%amp(ir12,ir3,:,alphaloc)=PSI1%amp(ir12,ir3,:,alphaloc)+phase*psicollect4(ir12,ir3,iq4,alphaloc)*besslp2r4(:,iq4,l4)
     END DO ! ip12
    END DO  ! ir3 
   END DO   ! iq4
  END DO  ! alpha
  
  DEALLOCATE(psicollect4)
  
 END SUBROUTINE fourier_p2r_psi4n31
 
 !! 4. 4N22 amplitudes 
 
 SUBROUTINE fourier_p2r_psi4n22(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (AMP4N22),INTENT(IN)    :: PSI2
  TYPE (RAMP4N22),INTENT(INOUT) :: PSI1
  INTEGER ip12,ip34,iq,alpha,alphaloc,ir12,ir34
  INTEGER l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari
  COMPLEX(dpreal) :: phase 
  COMPLEX(spreal),ALLOCATABLE :: psicollect12(:,:,:,:),psicollect34(:,:,:,:),psicollect(:,:,:,:),&
              &    psitrans12(:,:,:,:),psitrans34(:,:,:,:)
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  ! init besselfunctions if not done yet 
  CALL initbessel
  
  ! set zero and allocate if necessary 
  PSI1=0.0
    
  ALLOCATE(psicollect12(P12N,mynp34,mynq,mynbeta4N22))
  CALL collect_part_spcmplx(PSI2%amp,psicollect12,1,P12N,mynp34*mynq*mynbeta4N22,commp12)
  
  ALLOCATE(psitrans12(mynr12,mynp34,mynq,mynbeta4N22))
  psitrans12=0.0
  
  DO alphaloc=1,mynbeta4N22
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari) 
   phase=iu**l12    ! should be correct for p- to r-space 
   DO iq=1,mynq
    DO ip34=1,mynp34
     DO ip12=1,P12N
      psitrans12(:,ip34,iq,alphaloc)=psitrans12(:,ip34,iq,alphaloc)+phase*psicollect12(ip12,ip34,iq,alphaloc)*besslp2r12(:,ip12,l12)
     END DO ! ip12
    END DO  ! ip34 
   END DO   ! iq
  END DO  ! alpha

  DEALLOCATE(psicollect12)
  ALLOCATE(psicollect34(mynr12,P34N,mynq,mynbeta4N22))
  CALL collect_part_spcmplx(psitrans12,psicollect34,mynr12,P34N,mynq*mynbeta4N22,commp3)
  DEALLOCATE(psitrans12)
  
  ALLOCATE(psitrans34(mynr12,mynr34,mynq,mynbeta4N22))
  psitrans34=0.0
  
  DO alphaloc=1,mynbeta4N22
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari)   
   phase=iu**l34    ! should be correct for p- to r-space 
   DO iq=1,mynq
    DO ip34=1,P34N
     DO ir12=1,mynr12
      psitrans34(ir12,:,iq,alphaloc)=psitrans34(ir12,:,iq,alphaloc)+phase*psicollect34(ir12,ip34,iq,alphaloc)*besslp2r34(:,ip34,l34)
     END DO ! ir12
    END DO  ! ip34 
   END DO   ! iq
  END DO  ! alpha
  DEALLOCATE(psicollect34)
  
  ALLOCATE(psicollect(mynr12,mynr34,QN,mynbeta4N22))
  CALL collect_part_spcmplx(psitrans34,psicollect,mynr12*mynr34,QN,mynbeta4N22,commq4)
  DEALLOCATE(psitrans34)
  
  
  DO alphaloc=1,mynbeta4N22
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari)   
   phase=iu**lam    ! should be correct for p- to r-space 
   DO iq=1,QN
    DO ir34=1,mynr34
     DO ir12=1,mynr12
      PSI1%amp(ir12,ir34,:,alphaloc)=PSI1%amp(ir12,ir34,:,alphaloc)+phase*psicollect(ir12,ir34,iq,alphaloc)*besslp2r(:,iq,lam)
     END DO ! ip12
    END DO  ! ir34 
   END DO   ! iq
  END DO  ! alpha
  
  DEALLOCATE(psicollect)
  
 END SUBROUTINE fourier_p2r_psi4n22

 
 !! now add Fourier trafos from r-space to p-space "r2p"
 !! 1. NN amplitudes  
 
 SUBROUTINE fourier_r2p_psinn(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMPNN),INTENT(IN)     :: PSI2
  TYPE (AMPNN),INTENT(INOUT) :: PSI1
  INTEGER ir12,alpha,alphaloc
  INTEGER l12,s12,j12,t12,mt12
  COMPLEX(dpreal) :: phase 
  COMPLEX(spreal),ALLOCATABLE :: psicollect12(:,:)
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  ! init besselfunctions if not done yet 
  CALL initbessel
  
  ! set zero and allocate if necessary 
  PSI1=0.0
    
  ALLOCATE(psicollect12(R12N,mynalphann))
  CALL collect_part_spcmplx(PSI2%amp,psicollect12,1,R12N,mynalphann,commp12)
  
  DO alphaloc=1,mynalphann
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL getnnqn(alpha,l12,s12,j12,t12,mt12)
   phase=(-iu)**l12    ! should be correct for r- to p-space 
   DO ir12=1,R12N
    PSI1%amp(:,alphaloc)=PSI1%amp(:,alphaloc)+phase*psicollect12(ir12,alphaloc)*besslr2p12(:,ir12,l12)
   END DO ! ir12
  END DO  ! alpha
  
  DEALLOCATE(psicollect12)
  
 END SUBROUTINE fourier_r2p_psinn

  !! 2. 3N amplitudes 
 
 SUBROUTINE fourier_r2p_psi3n(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP3N),INTENT(IN)     :: PSI2
  TYPE (AMP3N),INTENT(INOUT) :: PSI1
  INTEGER ir12,ir3,alpha,alphaloc,ip12
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  COMPLEX(dpreal) :: phase 
  COMPLEX(spreal),ALLOCATABLE :: psicollect12(:,:,:),psicollect3(:,:,:),psitrans12(:,:,:)
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  ! init besselfunctions if not done yet 
  CALL initbessel
  
  ! set zero and allocate if necessary 
  PSI1=0.0
    
  ALLOCATE(psicollect12(R12N,mynr3,mynalpha3N))
  CALL collect_part_spcmplx(PSI2%amp,psicollect12,1,R12N,mynr3*mynalpha3N,commp12)
  
  ALLOCATE(psitrans12(mynp12,mynr3,mynalpha3N))
  psitrans12=0.0
  
  DO alphaloc=1,mynalpha3N
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)  
   phase=(-iu)**l12    ! should be correct for r- to p-space 
   DO ir3=1,mynr3
    DO ir12=1,R12N
     psitrans12(:,ir3,alphaloc)=psitrans12(:,ir3,alphaloc)+phase*psicollect12(ir12,ir3,alphaloc)*besslr2p12(:,ir12,l12)
    END DO ! ir12
   END DO  ! ir3 
  END DO  ! alpha

  DEALLOCATE(psicollect12)
  ALLOCATE(psicollect3(mynp12,R3N,mynalpha3N))
  CALL collect_part_spcmplx(psitrans12,psicollect3,mynp12,R3N,mynalpha3N,commp3)
  DEALLOCATE(psitrans12)
  
  
  DO alphaloc=1,mynalpha3N
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)  
   phase=(-iu)**l3    ! should be correct for r- to p-space 
   DO ir3=1,R3N
    DO ip12=1,mynp12
     PSI1%amp(ip12,:,alphaloc)=PSI1%amp(ip12,:,alphaloc)+phase*psicollect3(ip12,ir3,alphaloc)*besslr2p3(:,ir3,l3)
    END DO ! ip12
   END DO  ! ir3 
  END DO  ! alpha
  
  DEALLOCATE(psicollect3)
  
 END SUBROUTINE fourier_r2p_psi3n

 !! 3. 4N31 amplitudes 
 SUBROUTINE fourier_r2p_psi4n31(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP4N31),INTENT(IN)     :: PSI2
  TYPE (AMP4N31),INTENT(INOUT) :: PSI1
  INTEGER ir12,ir3,ir4,alpha,alphaloc,ip12,ip3
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari
  COMPLEX(dpreal) :: phase 
  COMPLEX(spreal),ALLOCATABLE :: psicollect12(:,:,:,:),psicollect3(:,:,:,:),psicollect4(:,:,:,:),&
         &       psitrans12(:,:,:,:),psitrans3(:,:,:,:)
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  ! init besselfunctions if not done yet 
  CALL initbessel
  
  ! set zero and allocate if necessary 
  PSI1=0.0
    
  ALLOCATE(psicollect12(R12N,mynr3,mynr4,mynalpha4N31))
  CALL collect_part_spcmplx(PSI2%amp,psicollect12,1,R12N,mynr3*mynr4*mynalpha4N31,commp12)
  
  ALLOCATE(psitrans12(mynp12,mynr3,mynr4,mynalpha4N31))
  psitrans12=0.0
  
  DO alphaloc=1,mynalpha4N31
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari)  
   phase=(-iu)**l12    ! should be correct for r- to p-space 
   DO ir4=1,mynr4
    DO ir3=1,mynr3
     DO ir12=1,R12N
      psitrans12(:,ir3,ir4,alphaloc)=psitrans12(:,ir3,ir4,alphaloc)+phase*psicollect12(ir12,ir3,ir4,alphaloc)*besslr2p12(:,ir12,l12)
     END DO ! ir12
    END DO  ! ir3 
   END DO   ! ir4
  END DO  ! alpha

  DEALLOCATE(psicollect12)
  ALLOCATE(psicollect3(mynp12,R3N,mynr4,mynalpha4N31))
  CALL collect_part_spcmplx(psitrans12,psicollect3,mynp12,R3N,mynr4*mynalpha4N31,commp3)
  DEALLOCATE(psitrans12)
  
  ALLOCATE(psitrans3(mynp12,mynp3,mynr4,mynalpha4N31))
  psitrans3=0.0
  
  DO alphaloc=1,mynalpha4N31
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari)    
   phase=(-iu)**l3    ! should be correct for r- to p-space 
   DO ir4=1,mynr4
    DO ir3=1,R3N
     DO ip12=1,mynp12
      psitrans3(ip12,:,ir4,alphaloc)=psitrans3(ip12,:,ir4,alphaloc)+phase*psicollect3(ip12,ir3,ir4,alphaloc)*besslr2p3(:,ir3,l3)
     END DO ! ip12
    END DO  ! ir3 
   END DO   ! ir4
  END DO  ! alpha
  DEALLOCATE(psicollect3)
  
  ALLOCATE(psicollect4(mynp12,mynp3,R4N,mynalpha4N31))
  CALL collect_part_spcmplx(psitrans3,psicollect4,mynp12*mynp3,R4N,mynalpha4N31,commq4)
  DEALLOCATE(psitrans3)
  
  
  DO alphaloc=1,mynalpha4N31
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari)    
   phase=(-iu)**l4    ! should be correct for r- to p-space 
   DO ir4=1,R4N
    DO ip3=1,mynp3
     DO ip12=1,mynp12
      PSI1%amp(ip12,ip3,:,alphaloc)=PSI1%amp(ip12,ip3,:,alphaloc)+phase*psicollect4(ip12,ip3,ir4,alphaloc)*besslr2p4(:,ir4,l4)
     END DO ! ip12
    END DO  ! ip3 
   END DO   ! ir4
  END DO  ! alpha
  
  DEALLOCATE(psicollect4)
  
 END SUBROUTINE fourier_r2p_psi4n31
 
 !! 4. 4N22 amplitudes 
 
 SUBROUTINE fourier_r2p_psi4n22(PSI1,PSI2)
  IMPLICIT NONE
  TYPE (RAMP4N22),INTENT(IN)     :: PSI2
  TYPE (AMP4N22),INTENT(INOUT) :: PSI1
  INTEGER ir12,ir34,ir,alpha,alphaloc,ip12,ip34
  INTEGER l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari
  COMPLEX(dpreal) :: phase 
  COMPLEX(spreal),ALLOCATABLE :: psicollect12(:,:,:,:),psicollect34(:,:,:,:),psicollect(:,:,:,:),&
           &   psitrans12(:,:,:,:),psitrans34(:,:,:,:)
  
#ifdef DEBUG
  IF(.NOT.allocated(PSI2%amp)) THEN
   STOP 'PSI2 not allocated in assigment !!'
  END IF
#endif
  ! init besselfunctions if not done yet 
  CALL initbessel
  
  ! set zero and allocate if necessary 
  PSI1=0.0
    
  ALLOCATE(psicollect12(R12N,mynr34,mynr,mynbeta4N22))
  CALL collect_part_spcmplx(PSI2%amp,psicollect12,1,R12N,mynr34*mynr*mynbeta4N22,commp12)
  
  ALLOCATE(psitrans12(mynp12,mynr34,mynr,mynbeta4N22))
  psitrans12=0.0
  
  DO alphaloc=1,mynbeta4N22
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari) 
   phase=(-iu)**l12    ! should be correct for r- to p-space 
   DO ir=1,mynr
    DO ir34=1,mynr34
     DO ir12=1,R12N
      psitrans12(:,ir34,ir,alphaloc)=psitrans12(:,ir34,ir,alphaloc)+phase*psicollect12(ir12,ir34,ir,alphaloc)*besslr2p12(:,ir12,l12)
     END DO ! ir12
    END DO  ! ir34 
   END DO   ! ir
  END DO  ! alpha

  DEALLOCATE(psicollect12)
  ALLOCATE(psicollect34(mynp12,R34N,mynr,mynbeta4N22))
  CALL collect_part_spcmplx(psitrans12,psicollect34,mynp12,R34N,mynr*mynbeta4N22,commp3)
  DEALLOCATE(psitrans12)
  
  ALLOCATE(psitrans34(mynp12,mynp34,mynr,mynbeta4N22))
  psitrans34=0.0
  
  DO alphaloc=1,mynbeta4N22
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari)   
   phase=(-iu)**l34    ! should be correct for r- to p-space 
   DO ir=1,mynr
    DO ir34=1,R34N
     DO ip12=1,mynp12
      psitrans34(ip12,:,ir,alphaloc)=psitrans34(ip12,:,ir,alphaloc)+phase*psicollect34(ip12,ir34,ir,alphaloc)*besslr2p34(:,ir34,l34)
     END DO ! ip12
    END DO  ! ir34 
   END DO   ! ir
  END DO  ! alpha
  DEALLOCATE(psicollect34)
  
  ALLOCATE(psicollect(mynp12,mynp34,RN,mynbeta4N22))
  CALL collect_part_spcmplx(psitrans34,psicollect,mynp12*mynp34,RN,mynbeta4N22,commq4)
  DEALLOCATE(psitrans34)
  
  
  DO alphaloc=1,mynbeta4N22
   alpha=myalphaid+1+(alphaloc-1)*npe_alpha
   CALL get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alpha12,alpha34,pari)   
   phase=(-iu)**lam    ! should be correct for r- to p-space 
   DO ir=1,RN
    DO ip34=1,mynp34
     DO ip12=1,mynp12
      PSI1%amp(ip12,ip34,:,alphaloc)=PSI1%amp(ip12,ip34,:,alphaloc)+phase*psicollect(ip12,ip34,ir,alphaloc)*besslr2p(:,ir,lam)
     END DO ! ip12
    END DO  ! ir34 
   END DO   ! iq
  END DO  ! alpha
  
  DEALLOCATE(psicollect)
  
 END SUBROUTINE fourier_r2p_psi4n22
 
 
 ! add routines to prepare simple test or starting amplitudes 
 
 SUBROUTINE testamp_nn(phiNN)
  IMPLICIT NONE
  TYPE(AMPNN) :: phiNN
  INTEGER localpha,ip12,ip3
  INTEGER alpha,l12,s12,j12,t12,mt12
  
  phiNN=0.0
  
  DO localpha=1,mynalphaNN                    ! all local 3N channels
   ! calculate global channel number 
   alpha=(localpha-1)*npe_alpha+1+myalphaid      
   ! and get quantum numbers 
   call getNNqn(alpha,l12,s12,j12,t12,mt12)

   DO ip12=1,mynp12                             ! all local p12 points 
    phiNN%amp(ip12,localpha) & 
 &    =P12P(myp12+ip12-1)**l12  & 
 &        *exp(-0.1*(P12P(myp12+ip12-1)**4)) &
 &        /(l12+1.)**2*(-1.0)**s12
   END DO ! ip12 

  END DO ! localpha
 END SUBROUTINE testamp_nn 
 
 
 SUBROUTINE testamp_3n(phi3N)
  IMPLICIT NONE
  TYPE(AMP3N) :: phi3N
  INTEGER localpha,ip12,ip3
  INTEGER alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alphaNN12,pari
  
  phi3N=0.0
  
  DO localpha=1,mynalpha3N                    ! all local 3N channels
   ! calculate global channel number 
   alpha=(localpha-1)*npe_alpha+1+myalphaid      
   ! and get quantum numbers 
   call get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alphaNN12,pari)

   DO ip3=1,mynp3                                ! all local p3 points 
    DO ip12=1,mynp12                             ! all local p12 points 
     phi3N%amp(ip12,ip3,localpha) & 
  &    =P12P(myp12+ip12-1)**l12*P3P(myp3+ip3-1)**l3  & 
  &        *exp(-0.1*(P12P(myp12+ip12-1)**4+P3P(myp3+ip3-1)**4)) &
  &        /((l12+1.)**2*(l3+1.)**2)*(-1.0)**(s12+(I3+1)/2)
    END DO ! ip12 
   END DO ! ip3
  END DO ! localpha
 END SUBROUTINE testamp_3n 
 
 SUBROUTINE testamp_4n31(phi4N31)
  IMPLICIT NONE
  TYPE(AMP4N31) :: phi4N31
  INTEGER localpha,ip12,ip3,iq4
  INTEGER alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alphaNN12,alpha3N,pari
  phi4N31=0.0
  
  DO localpha=1,mynalpha4N31                    ! all local 4N channels
   ! calculate global channel number 
   alpha=(localpha-1)*npe_alpha+1+myalphaid      
   ! and get quantum numbers 
   call get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3, &
  &                tau3,l4,I4,j4,tau4,mtau4,alphaNN12,alpha3N,pari)
  
!   call get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,  & 
!  &             tau3,mtau3,alphaNN12,pari)

   DO iq4=1,mynq4
    DO ip3=1,mynp3                                ! all local p3 points 
     DO ip12=1,mynp12                             ! all local p12 points 
!     sum=P12P(myp12+ip12-1)**2+0.75*P3P(myp3+ip3-1)**2
!     IF((l12.eq.0).AND.l3.eq.0) THEN
!       phifad%amp(ip12,ip3,localpha) = exp(-sum/1.0)*(-1.0)**s12
!     END IF
      phi4N31%amp(ip12,ip3,iq4,localpha) & 
   &    =P12P(myp12+ip12-1)**l12*P3P(myp3+ip3-1)**l3*Q4P(myq4+iq4-1)**l4  & 
   &        *exp(-0.1*(P12P(myp12+ip12-1)**4+P3P(myp3+ip3-1)**4+Q4P(myq4+iq4-1)**4)) &
   &        /((l12+1.)**2*(l3+1.)**2*(l4+1.)**2)*(-1.0)**(s12+(I3+1)/2)
     END DO ! ip12 
    END DO ! ip3
   END DO ! iq4 
  END DO ! localpha
 END SUBROUTINE testamp_4n31 
  
 SUBROUTINE testamp_4n22(phi4N22)
  IMPLICIT NONE
  TYPE(AMP4N22) :: phi4N22
  INTEGER localpha,ip12,ip34,iq
  INTEGER alpha,l12,s12,j12,t12,l34,s34,j34,t34,lam,I,j4,tau4,mtau4,alphaNN12,alphaNN34,pari
  
  phi4N22=0.0
  
  DO localpha=1,mynbeta4N22                    ! all local 4N channels
   ! calculate global channel number 
   alpha=(localpha-1)*npe_alpha+1+myalphaid      
   ! and get quantum numbers 
   call get4N22qn(alpha,l12,s12,j12,t12,l34,s34,j34,t34,  &
                  &   lam,I,j4,tau4,mtau4,alphaNN12,alphaNN34,pari)
   
!   call get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,  & 
!  &             tau3,mtau3,alphaNN12,pari)

   DO iq=1,mynq
    DO ip34=1,mynp34                                ! all local p3 points 
     DO ip12=1,mynp12                             ! all local p12 points 
!     sum=P12P(myp12+ip12-1)**2+0.75*P3P(myp3+ip3-1)**2
!     IF((l12.eq.0).AND.l3.eq.0) THEN
!       phifad%amp(ip12,ip3,localpha) = exp(-sum/1.0)*(-1.0)**s12
!     END IF
      phi4N22%amp(ip12,ip34,iq,localpha) & 
   &    =P12P(myp12+ip12-1)**l12*P34P(myp34+ip34-1)**l34*QP(myq+iq-1)**lam  & 
   &        *exp(-0.1*(P12P(myp12+ip12-1)**4+P34P(myp34+ip34-1)**4+QP(myq+iq-1)**4)) &
   &        /((l12+1.)**2*(l34+1.)**2*(lam+1.)**2)*(-1.0)**(s12+s34+I)
     END DO ! ip12 
    END DO ! ip34
   END DO ! iq
  END DO ! localpha
  
 END SUBROUTINE testamp_4n22
 
 !! print radii NN correlation and third particle 
 
 !! NN system (only sqrt(<r12**2>) exists, but this is printed and scaled 
    
 SUBROUTINE printrms_nn(PHI)
  IMPLICIT NONE 
  TYPE(AMPNN) PHI
  TYPE(RAMPNN) PHIR
  INTEGER ir12,ierr
  REAL(dpreal) :: rms12loc,rms12,norm(-1:1),normloc(-1:1),rmsloc(-1:1),rms(-1:1),fakt   ! rms(1) = proton, rms(-1)=neutron, rms(0)=nucleon 
  REAL(dpreal),ALLOCATABLE :: isofakt(:,:,:) 
  COMPLEX(spreal),ALLOCATABLE :: psiampNN(:,:)
  INTEGER alphap,alpha,localpha
  INTEGER l12,s12,j12,t12,mt12
  INTEGER l12p,s12p,j12p,t12p,mt12p
  
  PHIR=PHI  ! Fourier trafo to r-space 
  
  ! set isospin factors for 1/2*(1+/-tau_3(2)) projectors 
  ALLOCATE(isofakt(-1:1,alphaNNcdepmax,mynalphaNN))
  isofakt=0.0_dpreal 
  
  DO alphap=1,alphaNNcdepmax
   call getNNqn(alphap,l12p,s12p,j12p,t12p,mt12p)
   DO localpha=1,mynalphaNN
    alpha=myalphaid+1+npe_alpha*(localpha-1)    
    call getNNqn(alpha,l12,s12,j12,t12,mt12)
    IF(l12.EQ.l12p .AND. s12.EQ. s12p .AND. j12.EQ.j12p &
   &  .AND. mt12 .EQ. mt12p) THEN
     IF(t12.EQ.t12p) THEN     
      isofakt(-1,alphap,localpha)=0.5_dpreal
      isofakt( 0,alphap,localpha)=1.0_dpreal
      isofakt( 1,alphap,localpha)=0.5_dpreal
     END IF
     
     isofakt(-1,alphap,localpha)=isofakt(-1,alphap,localpha) &
     & - 0.5_dpreal*CG(2*t12,2,2*t12p,2*mt12,0,2*mt12p)  &
     &     *6*sqrt(2.0_dpreal*t12+1)*C9J(1,1,0,1,1,2,2*t12p,2*t12,2)
     isofakt(1,alphap,localpha)=isofakt(1,alphap,localpha) &
     & + 0.5_dpreal*CG(2*t12,2,2*t12p,2*mt12,0,2*mt12p)  &
     &     *6*sqrt(2.0_dpreal*t12+1)*C9J(1,1,0,1,1,2,2*t12p,2*t12,2)
    END IF
     
   END DO
  END DO
  
  
  ! collect all partial waves on one side 
  ALLOCATE(psiampNN(mynr12,alphaNNcdepmax))
  CALL collect_piece_spcmplx(PHIR%amp,psiampNN,mynr12,alphaNNcdepmax,1,commalpha)
          
  
  ! to r12 radius first
  
  rms12loc=0.0
  rmsloc=0.0
  normloc=0.0
  
  DO alphap=1,alphaNNcdepmax
   DO alpha=1,mynalphaNN
    IF(isofakt(-1,alphap,alpha).NE.0.0) THEN 
     DO ir12=1,mynr12 
      fakt=conjg(PHIR%amp(ir12,alpha))*psiampNN(ir12,alphap)*r12weight(ir12)
      rms12loc=rms12loc+ fakt*R12P(ir12+myr12-1)**2*isofakt(0,alphap,alpha)
      normloc(-1)=normloc(-1)+fakt*isofakt(-1,alphap,alpha)
      normloc( 0)=normloc( 0)+fakt*isofakt( 0,alphap,alpha)
      normloc( 1)=normloc( 1)+fakt*isofakt( 1,alphap,alpha)
      rmsloc(-1)=rmsloc(-1)+fakt*R12P(ir12+myr12-1)**2*isofakt(-1,alphap,alpha)
      rmsloc( 0)=rmsloc( 0)+fakt*R12P(ir12+myr12-1)**2*isofakt( 0,alphap,alpha)
      rmsloc( 1)=rmsloc( 1)+fakt*R12P(ir12+myr12-1)**2*isofakt( 1,alphap,alpha)
     END DO ! ir12 
    END IF
   END DO ! alpha
  END DO !  alphap
  
  DEALLOCATE(isofakt)
  
  CALL MPI_ALLREDUCE(rms12loc,rms12,1,MPI_REAL8,MPI_SUM,commampnn,ierr) 
  CALL MPI_ALLREDUCE(rmsloc,rms,3,MPI_REAL8,MPI_SUM,commampnn,ierr) 
  CALL MPI_ALLREDUCE(normloc,norm,3,MPI_REAL8,MPI_SUM,commampnn,ierr) 
  
  rms12=sqrt(rms12)
  rms=sqrt(rms)*1.0_dpreal/2.0_dpreal   ! assuming equal masses of nucleons 
  
  
  IF(master) THEN
   WRITE(*,*)
   WRITE(*,'(A,3E15.6)') 'NN RMS RADII in fm: (R12 and Rnucl,rspace norm) :',rms12,rms(0),norm(0)
   WRITE(*,'(A,2E15.6)') 'NN P-RMS RADIUS in fm (and norm): ',rms(1)/sqrt(norm(1)),norm(1)
   WRITE(*,'(A,2E15.6)') 'NN N-RMS RADIUS in fm (and norm): ',rms(-1)/sqrt(norm(-1)),norm(-1)
   WRITE(*,*) 
  END IF
  
  DEALLOCATE(PHIR%amp,psiampNN)
 END SUBROUTINE printrms_nn
 
 SUBROUTINE printrms_3n(PHI)
  IMPLICIT NONE 
  TYPE(AMP3N) PHI
  TYPE(RAMP3N) PHIR
  INTEGER ir12,ir3,ierr
  REAL(dpreal) :: rms12loc,rms12,norm(-1:1),normloc(-1:1),rmsloc(-1:1),rms(-1:1),fakt
  REAL(dpreal),ALLOCATABLE :: isofakt(:,:,:) 
  COMPLEX(spreal),ALLOCATABLE :: psiamp3N(:,:,:)
  INTEGER alphap,alpha,localpha
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari
  INTEGER l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,mtau3p,alpha2bp,parip
  
  PHIR=PHI  ! Fourier trafo to r-space 
  
  ! set isospin factors for 1/2*(1+/-tau_3(3)) projectors 
  ALLOCATE(isofakt(-1:1,alpha3Ncdepmax,mynalpha3N))
  isofakt=0.0_dpreal 
  
  DO alphap=1,alpha3Ncdepmax
   call get3Nqn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,mtau3p,alpha2bp,parip)
   DO localpha=1,mynalpha3N
    alpha=myalphaid+1+npe_alpha*(localpha-1)    
    call get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alpha2b,pari)
    IF(l12.EQ.l12p .AND. s12.EQ. s12p .AND. j12.EQ.j12p &
   &  .AND. l3.EQ.l3p .AND. I3.EQ.I3p .AND. j3.EQ.j3p & 
   &  .AND. mtau3 .EQ. mtau3p) THEN
     IF(t12.EQ.t12p .AND. tau3.EQ.tau3p) THEN     
      isofakt(-1,alphap,localpha)=0.5_dpreal
      isofakt( 0,alphap,localpha)=1.0_dpreal
      isofakt( 1,alphap,localpha)=0.5_dpreal
     END IF
     
     isofakt(-1,alphap,localpha)=isofakt(-1,alphap,localpha) &
     & - 0.5_dpreal*CG(tau3,2,tau3p,mtau3,0,mtau3p)  &
     &     *3*sqrt(2*(2.0_dpreal*t12p+1)*(tau3+1))*C9J(2*t12p,2*t12,0,1,1,2,tau3p,tau3,2)
     isofakt(1,alphap,localpha)=isofakt(1,alphap,localpha) &
     & + 0.5_dpreal*CG(tau3,2,tau3p,mtau3,0,mtau3p)  &
     &     *3*sqrt(2*(2.0_dpreal*t12p+1)*(tau3+1))*C9J(2*t12p,2*t12,0,1,1,2,tau3p,tau3,2)
    END IF     
   END DO
  END DO
   
  ! collect all partial waves on one side 
  ALLOCATE(psiamp3N(mynr12,mynr3,alpha3Ncdepmax))
  CALL collect_piece_spcmplx(PHIR%amp,psiamp3N,mynr12*mynr3,alpha3Ncdepmax,1,commalpha)
          
  ! to r12 radius first
  
  rms12loc=0.0
  rmsloc=0.0
  normloc=0.0
 
  DO alphap=1,alpha3Ncdepmax
   DO alpha=1,mynalpha3N
    IF(isofakt(-1,alphap,alpha).NE.0.0) THEN 
     DO ir3=1,mynr3
      DO ir12=1,mynr12 
       fakt=conjg(PHIR%amp(ir12,ir3,alpha))*psiamp3N(ir12,ir3,alphap)*r12r3weight(ir12,ir3)
       rms12loc=rms12loc+fakt*R12P(ir12+myr12-1)**2*isofakt( 0,alphap,alpha)
       normloc(-1)=normloc(-1)+fakt*isofakt(-1,alphap,alpha)
       normloc( 0)=normloc( 0)+fakt*isofakt( 0,alphap,alpha)
       normloc( 1)=normloc( 1)+fakt*isofakt( 1,alphap,alpha)
       rmsloc(-1)=rmsloc(-1)+fakt*R3P(ir3+myr3-1)**2*isofakt(-1,alphap,alpha)
       rmsloc( 0)=rmsloc( 0)+fakt*R3P(ir3+myr3-1)**2*isofakt( 0,alphap,alpha)
       rmsloc( 1)=rmsloc( 1)+fakt*R3P(ir3+myr3-1)**2*isofakt( 1,alphap,alpha)
      END DO ! ir12 
     END DO  ! ir3 
    END IF
   END DO ! alpha
  END DO !  alphap
  
  
  CALL MPI_ALLREDUCE(rms12loc,rms12,1,MPI_REAL8,MPI_SUM,commamp3n,ierr) 
  CALL MPI_ALLREDUCE(rmsloc,rms,3,MPI_REAL8,MPI_SUM,commamp3n,ierr) 
  CALL MPI_ALLREDUCE(normloc,norm,3,MPI_REAL8,MPI_SUM,commamp3n,ierr) 
  
  rms12=sqrt(rms12)
  rms=sqrt(rms)*2.0_dpreal/3.0_dpreal   ! assuming equal masses of nucleons 
  
  IF(master) THEN
   WRITE(*,*)
   WRITE(*,'(A,3E15.6)') '3N RMS RADII in fm: (R12 and Rnucl,rspace norm) :',rms12,rms(0),norm(0)
   WRITE(*,'(A,2E15.6)') '3N P-RMS RADIUS in fm (and norm): ',rms(1)/sqrt(norm(1)),norm(1)
   WRITE(*,'(A,2E15.6)') '3N N-RMS RADIUS in fm (and norm): ',rms(-1)/sqrt(norm(-1)),norm(-1)
   WRITE(*,*) 
  END IF
  
  DEALLOCATE(PHIR%amp,psiamp3N,isofakt)
 END SUBROUTINE printrms_3n
 
SUBROUTINE printrms_4n31(PHI)
  IMPLICIT NONE 
  TYPE(AMP4N31) PHI
  TYPE(RAMP4N31) PHIR
  INTEGER ir12,ir3,ir4,ierr
  REAL(dpreal) :: rms12loc,rms12,norm(-1:1),normloc(-1:1),rmsloc(-1:1),rms(-1:1),fakt
  REAL(dpreal),ALLOCATABLE :: isofakt(:,:,:) 
  COMPLEX(spreal),ALLOCATABLE :: psiamp4N31(:,:,:,:)
  INTEGER alphap,alpha,localpha
  INTEGER l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari
  INTEGER l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p,alpha2bp,alpha3bp,parip
  
  PHIR=PHI  ! Fourier trafo to r-space 
  
   ! set isospin factors for 1/2*(1+/-tau_3(3)) projectors 
  ALLOCATE(isofakt(-1:1,alpha4N31cdepmax,mynalpha4N31))
  isofakt=0.0_dpreal 
  
  DO alphap=1,alpha4N31cdepmax
   call get4N31qn(alphap,l12p,s12p,j12p,t12p,l3p,I3p,j3p,tau3p,l4p,I4p,j4p,tau4p,mtau4p,alpha2bp,alpha3bp,parip)
   DO localpha=1,mynalpha4N31 
    alpha=myalphaid+1+npe_alpha*(localpha-1)   
    call get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,tau3,l4,I4,j4,tau4,mtau4,alpha2b,alpha3b,pari) 
    IF(l12.EQ.l12p .AND. s12.EQ. s12p .AND. j12.EQ.j12p &
   &  .AND. l3.EQ.l3p .AND. I3.EQ.I3p .AND. j3.EQ.j3p & 
   &  .AND. l4.EQ.l4p .AND. I4.EQ.I4p .AND. j4.EQ.j4p &
   &  .AND. mtau4 .EQ. mtau4p) THEN
     IF(t12.EQ.t12p .AND. tau3.EQ.tau3p .AND. tau4.EQ.tau4p) THEN     
      isofakt(-1,alphap,localpha)=0.5_dpreal
      isofakt( 0,alphap,localpha)=1.0_dpreal
      isofakt( 1,alphap,localpha)=0.5_dpreal
     END IF
     
     isofakt(-1,alphap,localpha)=isofakt(-1,alphap,localpha) &
     & - 0.5_dpreal*CG(2*tau4,2,2*tau4p,2*mtau4,0,2*mtau4p)  &
     &     *3*sqrt(2*(2.0_dpreal*tau4+1)*(tau3p+1))*C9J(tau3p,tau3,0,1,1,2,2*tau4p,2*tau4,2)
     isofakt(1,alphap,localpha)=isofakt(1,alphap,localpha) &
     & + 0.5_dpreal*CG(2*tau4,2,2*tau4p,2*mtau4,0,2*mtau4p)  &
     &     *3*sqrt(2*(2.0_dpreal*tau4+1)*(tau3p+1))*C9J(tau3p,tau3,0,1,1,2,2*tau4p,2*tau4,2)
    END IF     
   END DO
  END DO
   
  ! collect all partial waves on one side 
  ALLOCATE(psiamp4N31(mynr12,mynr3,mynr4,alpha4N31cdepmax))
  CALL collect_piece_spcmplx(PHIR%amp,psiamp4N31,mynr12*mynr3*mynr4,alpha4N31cdepmax,1,commalpha)
 
  ! to r12 radius first
  
  rms12loc=0.0
  rmsloc=0.0
  normloc=0.0
  
  DO alphap=1,alpha4N31cdepmax
   DO alpha=1,mynalpha4N31
    IF(isofakt(-1,alphap,alpha).NE.0.0) THEN 
     DO ir4=1,mynr4
      DO ir3=1,mynr3
       DO ir12=1,mynr12 
        fakt=conjg(PHIR%amp(ir12,ir3,ir4,alpha))*psiamp4N31(ir12,ir3,ir4,alphap)*r12r3r4weight(ir12,ir3,ir4)
        rms12loc=rms12loc+ fakt*R12P(ir12+myr12-1)**2*isofakt( 0,alphap,alpha) 
        normloc(-1)=normloc(-1)+fakt*isofakt(-1,alphap,alpha)
        normloc( 0)=normloc( 0)+fakt*isofakt( 0,alphap,alpha)
        normloc( 1)=normloc( 1)+fakt*isofakt( 1,alphap,alpha)
        rmsloc(-1)=rmsloc(-1)+fakt*R4P(ir4+myr4-1)**2*isofakt(-1,alphap,alpha)
        rmsloc( 0)=rmsloc( 0)+fakt*R4P(ir4+myr4-1)**2*isofakt( 0,alphap,alpha)
        rmsloc( 1)=rmsloc( 1)+fakt*R4P(ir4+myr4-1)**2*isofakt( 1,alphap,alpha)
       END DO ! ir12 
      END DO  ! ir3 
     END DO   ! ir4 
    END IF
   END DO ! alpha
  END DO !  alphap
  
  
  CALL MPI_ALLREDUCE(rms12loc,rms12,1,MPI_REAL8,MPI_SUM,commamp4n,ierr) 
  CALL MPI_ALLREDUCE(rmsloc,rms,3,MPI_REAL8,MPI_SUM,commamp4n,ierr) 
  CALL MPI_ALLREDUCE(normloc,norm,3,MPI_REAL8,MPI_SUM,commamp4n,ierr) 
  
  rms12=sqrt(rms12)
  rms=sqrt(rms)*3.0_dpreal/4.0_dpreal   ! assuming equal masses of nucleons 
  
  IF(master) THEN
   WRITE(*,*)
   WRITE(*,'(A,3E15.6)') '4N31 RMS RADII in fm: (R12 and Rnucl,rspace norm) :',rms12,rms(0),norm(0)
   WRITE(*,'(A,2E15.6)') '4N31 P-RMS RADIUS in fm (and norm): ',rms(1)/sqrt(norm(1)),norm(1)
   WRITE(*,'(A,2E15.6)') '4N31 N-RMS RADIUS in fm (and norm): ',rms(-1)/sqrt(norm(-1)),norm(-1)
   WRITE(*,*) 
  END IF
  
  DEALLOCATE(PHIR%amp,psiamp4N31,isofakt)
 END SUBROUTINE printrms_4n31

!> subroutine prints RMS radii based on the Fourier transformed wave function in 2+2 coordinates.
!! @param[in] PHI 4N wave function in 2+2 coordinates
!! The radii are not devided by the norm. 
 SUBROUTINE printrms_4n22(PHI)
  IMPLICIT NONE 
  TYPE(AMP4N22) PHI
  TYPE(RAMP4N22) PHIR
  INTEGER ir12,ir34,ir,alpha,ierr
  REAL(dpreal) :: rms12loc,rms12,rms34loc,rms34,rmsRloc,rmsR,norm,normloc,fakt
  
  PHIR=PHI  ! Fourier trafo to r-space 
  
  ! to r12,r34 and R radius 
  
  rms12loc=0.0
  rms34loc=0.0
  rmsRloc=0.0
  normloc=0.0
  DO alpha=1,mynbeta4N22
   DO ir=1,mynr
    DO ir34=1,mynr34
     DO ir12=1,mynr12 
      fakt=abs(PHIR%amp(ir12,ir34,ir,alpha))**2*r12r34rweight(ir12,ir34,ir)
      rms12loc=rms12loc+ fakt*R12P(ir12+myr12-1)**2
      rms34loc=rms34loc+ fakt*R34P(ir34+myr34-1)**2
      rmsRloc=rmsRloc+ fakt*RP(ir+myr-1)**2
      normloc=normloc+fakt
     END DO ! ir12 
    END DO  ! ir3 
   END DO   ! ir4 
  END DO !  alpha
  
 
  CALL MPI_ALLREDUCE(rms12loc,rms12,1,MPI_REAL8,MPI_SUM,commamp4n,ierr) 
  CALL MPI_ALLREDUCE(rms34loc,rms34,1,MPI_REAL8,MPI_SUM,commamp4n,ierr) 
  CALL MPI_ALLREDUCE(rmsRloc,rmsR,1,MPI_REAL8,MPI_SUM,commamp4n,ierr) 
  CALL MPI_ALLREDUCE(normloc,norm,1,MPI_REAL8,MPI_SUM,commamp4n,ierr) 
  
  rms12=sqrt(rms12)
  rms34=sqrt(rms34)
  rmsR=sqrt(rmsR)
  
  IF(master) THEN
   WRITE(*,*)
   WRITE(*,'(A,2E15.6)') '4N22 RMS RADII in fm: (R12 and rspace norm) :',rms12,norm
   WRITE(*,'(A,2E15.6)') '4N22 RMS RADII in fm: (R34 and R)           :',rms34,rmsR
   WRITE(*,*) 
  END IF
  
  DEALLOCATE(PHIR%amp)
 END SUBROUTINE printrms_4n22
 
END MODULE amplitudes
