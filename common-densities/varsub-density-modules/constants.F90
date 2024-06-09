#include "fdefs.h"
module constants
 use precision
 implicit none


#ifdef SU3MASSES
 real(dpreal),parameter :: hbarc=197.3269718_dpreal         ! NIST Ref 02.12.2014              
 real(dpreal),parameter :: mprot=938.94298_dpreal/hbarc    ! NIST Ref 02.12.2014/PDG 2014
 real(dpreal),parameter :: mneu=938.94298_dpreal/hbarc     ! NIST Ref 02.12.2014/PDG 2014
 real(dpreal),parameter :: amunit=931.4940954_dpreal/hbarc  ! CODATA 2014 (Rev. Mod. Phys. 88,035009)
 real(dpreal),parameter :: mlam=938.94298_dpreal/hbarc       ! PDG 2014
 real(dpreal),parameter :: msigp=938.94298_dpreal/hbarc       ! PDG 2014
 real(dpreal),parameter :: msig0=938.94298_dpreal/hbarc      ! PDG 2014
 real(dpreal),parameter :: msigm=938.94298_dpreal/hbarc      ! PDG 2014
 real(dpreal),parameter :: mxi0=938.94298_dpreal/hbarc       ! PDG 2014
 real(dpreal),parameter :: mxim=938.94298_dpreal/hbarc       ! PDG 2014
 real(dpreal),parameter :: mxiave=(mxi0+mxim)/2.0_dpreal
 real(dpreal),parameter :: msigave=(msigp+msig0+msigm)/3.0_dpreal
 real(dpreal),parameter :: mnuc=2*mneu*mprot/(mneu+mprot)
 real(dpreal),parameter :: m_pip=138.0_dpreal/hbarc     ! PDG 2014
 real(dpreal),parameter :: m_pi0=138.0_dpreal/hbarc      ! PDG 2014
 real(dpreal),parameter :: m_pi=(2.0_dpreal*m_pip+m_pi0)/3.0_dpreal
 real(dpreal),parameter :: alpha_fein = 1.0_dpreal/137.035999074_dpreal ! NIST Ref 02.12.2014
 real(dpreal),parameter :: echarge2=alpha_fein*hbarc        ! e^2 = alpha * hbarc (no 4 pi e0) 
 real(dpreal),parameter :: m_eta=138.0_dpreal/hbarc       ! PDG 2014
 real(dpreal),parameter :: m_rho=775.26_dpreal/hbarc        ! PDG 2014
 real(dpreal),parameter :: m_omega=782.65_dpreal/hbarc      ! PDG 2014
 real(dpreal),parameter :: m_kaonp=138.0_dpreal/hbarc     ! PDG 2014
 real(dpreal),parameter :: m_kaon0=138.0_dpreal/hbarc     ! PDG 2014
 real(dpreal),parameter :: m_kaon=(m_kaonp+m_kaon0)/2.0_dpreal
 real(dpreal),parameter :: m_xim=938.94298_dpreal/hbarc
 real(dpreal),parameter :: m_xi0=938.94298_dpreal/hbarc
 real(dpreal),parameter :: m_xiave=(m_xim+m_xi0)/2.0_dpreal 
#else 
 real(dpreal),parameter :: hbarc=197.3269718_dpreal         ! NIST Ref 02.12.2014              
 real(dpreal),parameter :: mprot=938.272046_dpreal/hbarc    ! NIST Ref 02.12.2014/PDG 2014
 real(dpreal),parameter :: mneu=939.565379_dpreal/hbarc     ! NIST Ref 02.12.2014/PDG 2014
 real(dpreal),parameter :: amunit=931.4940954_dpreal/hbarc  ! CODATA 2014 (Rev. Mod. Phys. 88,035009)
 real(dpreal),parameter :: mlam=1115.683_dpreal/hbarc       ! PDG 2014
 real(dpreal),parameter :: msigp=1189.37_dpreal/hbarc       ! PDG 2014
 real(dpreal),parameter :: msig0=1192.642_dpreal/hbarc      ! PDG 2014
 real(dpreal),parameter :: msigm=1197.449_dpreal/hbarc      ! PDG 2014
 real(dpreal),parameter :: mxi0=1314.86_dpreal/hbarc       ! PDG 2014
 real(dpreal),parameter :: mxim=1321.71_dpreal/hbarc       ! PDG 2014
 real(dpreal),parameter :: mxiave=(mxi0+mxim)/2.0_dpreal
 real(dpreal),parameter :: msigave=(msigp+msig0+msigm)/3.0_dpreal
 real(dpreal),parameter :: mnuc=2*mneu*mprot/(mneu+mprot)
 real(dpreal),parameter :: m_pip=139.57018_dpreal/hbarc     ! PDG 2014
 real(dpreal),parameter :: m_pi0=134.9766_dpreal/hbarc      ! PDG 2014
 real(dpreal),parameter :: m_pi=(2.0_dpreal*m_pip+m_pi0)/3.0_dpreal
 real(dpreal),parameter :: alpha_fein = 1.0_dpreal/137.035999074_dpreal ! NIST Ref 02.12.2014
 real(dpreal),parameter :: echarge2=alpha_fein*hbarc        ! e^2 = alpha * hbarc (no 4 pi e0)
 real(dpreal),parameter :: m_eta=547.862_dpreal/hbarc       ! PDG 2014
 real(dpreal),parameter :: m_rho=775.26_dpreal/hbarc        ! PDG 2014
 real(dpreal),parameter :: m_omega=782.65_dpreal/hbarc      ! PDG 2014
 real(dpreal),parameter :: m_kaonp=493.677_dpreal/hbarc     ! PDG 2014
 real(dpreal),parameter :: m_kaon0=497.614_dpreal/hbarc     ! PDG 2014
 real(dpreal),parameter :: m_kaon=(m_kaonp+m_kaon0)/2.0_dpreal
 real(dpreal),parameter :: m_xim=mxim
 real(dpreal),parameter :: m_xi0=mxi0
 real(dpreal),parameter :: m_xiave=(m_xim+m_xi0)/2.0_dpreal
 ! real(dpreal),parameter ::mprot=938.94298/hbarc) ! gemittelte physikalische Massen 
 ! real(dpreal),parameter ::mneu=938.94298/hbarc) ! hbar**2/m=41.47
 ! real(dpreal),parameter ::mprot=938.92635/hbarc)     ! Miyagawa's Massen 
 ! real(dpreal),parameter ::mneu=938.92635/hbarc)     
 ! real(dpreal),parameter :: mlam=1115.60/hbarc)
#endif
 
CONTAINS

SUBROUTINE printconstants
 IMPLICIT NONE
 
#ifdef SU3MASSES
 WRITE(*,*) 'SU3 symmetric, non-physical masses are used!' 
#else
 WRITE(*,*) 'physical masses are used!'
#endif 
 
 WRITE(*,'(A,E15.6,5X,A)') 'hbarc  =   ',hbarc,'! NIST Ref 02.12.2014'              
 WRITE(*,'(A,E15.6,5X,A)') 'mprot  =   ',mprot*hbarc,'! NIST Ref 02.12.2014/PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'mneu   =   ',mneu*hbarc,'! NIST Ref 02.12.2014/PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'amunit =   ',amunit*hbarc,'! CODATA 2014 (Rev. Mod. Phys. 88,035009)'
 WRITE(*,'(A,E15.6,5X,A)') 'mlam   =   ',mlam*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'msigp  =   ',msigp*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'msig0  =   ',msig0*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'msigm  =   ',msigm*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'mxi0   =   ',mxi0*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'mxim   =   ',mxim*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'mxiave =   ',mxiave*hbarc,'' 
 WRITE(*,'(A,E15.6,5X,A)') 'msigave=   ',msigave*hbarc,'' 
 WRITE(*,'(A,E15.6,5X,A)') 'mnuc   =   ',mnuc*hbarc,'' 
 WRITE(*,'(A,E15.6,5X,A)') 'm_pip  =   ',m_pip*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'm_pi0  =   ',m_pi0*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'm_pi   =   ',m_pi*hbarc 
 WRITE(*,'(A,E15.6,5X,A)') '1/alpha=   ',1.0_dpreal/alpha_fein,'! NIST Ref 02.12.2014'
 WRITE(*,'(A,E15.6,5X,A)') 'echarge2=  ',echarge2,'! e^2 = alpha * hbarc (no 4 pi e0)'
 WRITE(*,'(A,E15.6,5X,A)') 'm_eta  =   ',m_eta*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'm_rho  =   ',m_rho*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'm_omega=   ',m_omega*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'm_kaonp=   ',m_kaonp*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'm_kaon0=   ',m_kaon0*hbarc,'! PDG 2014'
 WRITE(*,'(A,E15.6,5X,A)') 'm_kaon =   ',m_kaon*hbarc,'' 
 WRITE(*,'(A,E15.6,5X,A)') 'm_xim  =   ',m_xim*hbarc,'' 
 WRITE(*,'(A,E15.6,5X,A)') 'm_xi0  =   ',m_xi0*hbarc,'' 
 WRITE(*,'(A,E15.6,5X,A)') 'm_xiave=   ',m_xiave*hbarc,'' 
 
 
END SUBROUTINE printconstants
end module constants

