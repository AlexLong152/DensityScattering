c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018
c
c     This file based on deuteron onebody/polnsamp.f Version 7 hgrie Feb 2013 -- Nov 2014
c
c     contains subroutines
c             trafo1Nampxx: Initial photon polarization=x, final photon polarization=x
c             trafo1Nampyx: Initial photon polarization=y, final photon polarization=x
c             trafo1Nampxy: Initial photon polarization=x, final photon polarization=y
c             trafo1Nampyy: Initial photon polarization=y, final photon polarization=y
c      
c     NOTE:
c     Input are the amplitudes Ai in the basis of [review (2.1)],
c     i.e. where ε&k structures are unit vectors!
c     NB: the deuteron subroutines in onebody/polnsamp.f do NOT work in that basis,
c     but in a basis where A[2,4-6] carry an additional factor of ω²!
c     I have removed these instances here. 
c
c     NOTE: in the  deuteron subroutines in onebody/polnsamp.f,
c     variable hold() get an overall "-" sign  for trafo1NampXY and trafo1NampYX (final do-loop in these routines!)
c     In the routines below, I ELMINIATED this final do-loop and implemented the sign directly in the amp defs.
c     Comment in polnamp.f by Daniel: "The minus sign arises because the mathematica script was
c                run with a choice for epsx and epsxp which corresponds to a left-handed coordinate system"
c      
c     New output extensively checked against module (several incommensurate angles & energies) in
c     check-trafo1Namps-to-helicity.nb:
c      
c      spinstructure[{i1p_, i1_}(*individual nucleon spin combination*),
c         {\[Epsilon]vec1_(*in*), \[Epsilon]vec2_(*out*)}(*polarisation vectors contracted with \[Sigma][i]*),
c         {kvec1_(*in*),kvec2_(*out*)}(*momentum vectors contracted with \[Sigma][i] -- NORMALISED TO UNITY !!!*),
c         {A1_, A2_, A3_, A4_, A5_, A6_}(*fundamental 1N amplitudes*)]
c
c     taken from mathematica notebook DEUTERON/UNifIED/fortran/mathematica/make-onebody-amplitudes.v7.0.nb
c
c     which also takes the Ai of the SAME basis [review (2.1)] as input (1N amps in "hatted" basis).
c
c     This confirms multiple checks performed by Daniel, Deepshikha etc.
c     PRE-HISTORY from deuteron code (with factors ω² in A[2,4-6]):
c
c        Checked in nucleon code for arbitrary angles-7/98
c        Corrected missing factor of i in spin-dependent pieces-5/04
c        Corrected inconsistent notation for initial/final poln state-5/04
c
c     twoSmax/twoMz dependence: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: EXTENSIVE doCUMENTATION FOR trafo1Nampxx, minimal docu for the others.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c      
c     hgrie May 2018: more extensive description of changes in main.*.f
c                     rewritten such that magnetic quantum numbers are "2xMz" etc
c
c**********************************************************************
      subroutine trafo1Nampxx(hold,A1,A2,A3,A4,A5,A6,
     &     alpha,beta,alphap,betap,verbosity)
c
c**********************************************************************
c   Calculates the single-nucleon c.m. frame amplitude for given "invariant"
c     amplitudes A1-A6 and incoming photon angles (alpha,beta), outgoing
c     (alphap,betap). (For notation see BKM review).
c
c   Initial photon polarization=x, final photon polarization=x
c
c**********************************************************************
      implicit none
c
      include '../common-densities/constants.def'
c
c======================================================================
c
c  OUTPUT: amplitude in (twom1Np,twom1N) basis: 2*1N spin projection: integer  
c
      complex*16,intent(out) :: hold(-1:1,-1:1)
c
c     Holds the two-body gammaNN->gammaNN amplitude for
c         polarization y->x amplitude on a basis of single-nucleon spins.
c     nucleon 1 is struck
c      
c     hgrie May 2018: Using mathematica routine spinstructure[] above, checked also that
c     in hold(i,j), i= is OUTGOING nucleon spin, j= INCOMING nucleon spin (sign of (1,2) vs (2,1)!).
c     NB: The final sign-change-do-loop in onebody/polnsamp.f for
c         "Comptonxy" and "Conptonyx" also confirms (yet again) that ordering of indices
c
c======================================================================
c
c  INPUT:
c
      real*8,intent(in) :: A1,A2,A3,A4,A5,A6
      real*8,intent(in) :: alpha,beta,alphap,betap
      integer,intent(in) :: verbosity ! set verbosity level of stdout
c
c     A1-A6-invariant amplitudes for given omega and theta (c.m. frame)
c     alpha,beta-orientations theta and phi of incoming photon (typically both 0)
c     alphap,betap-orientations theta and phi of outgoing photon
c
c======================================================================
c
c  LOCAL: 
c
      real*8 ca,cap,sa,sap,cb,cbp,sb,sbp ! trig functions of photon angles
c
c**********************************************************************
c
      if (verbosity.eq.1000) continue
      hold=c0
      ca=dcos(alpha)
      cap=dcos(alphap)
      sa=dsin(alpha)
      sap=dsin(alphap)
      cb=dcos(beta)
      cbp=dcos(betap)
      sb=dsin(beta)
      sbp=dsin(betap)
      hold(1,1)=(ca*cap*cb*cbp + sa*sap + ca*cap*sb*sbp)*A1
     &     + (-cap*cb*cbp*sa +  ca*sap - cap*sa*sb*sbp)
     &     *(cap*sa - ca*cb*cbp*sap - ca*sap*sb*sbp)*A2
      hold(-1,-1)=hold(1,1)
c
c  This does not contribute to the spin-flip
c
      hold(1,-1)=c0
      hold(1,1)=hold(1,1) + ca*cap*dsin(beta - betap)*ci*A3
      hold(1,-1)=hold(1,-1) + ((-cap*sa*sbp + ca*sb*sap)
     &                        - ci*(cbp*cap*sa - ca*cb*sap))*ci*A3
      hold(-1,-1)=hold(-1,-1) - ca*cap*dsin(beta - betap)*ci*A3
      hold(1,1)=hold(1,1) + sa*dsin(beta - betap)*sap
     &     *(ca*cb*cbp*cap + ca*cap*sb*sbp + sa*sap)*ci*A4
      hold(1,-1)=hold(1,-1) + (ca*cb*cbp*cap + ca*cap*sb*sbp + sa*sap)
     &     *(-cap*sa*sb+ca*sbp*sap - ci*(cb*cap*sa - ca*cbp*sap))*ci*A4
      hold(-1,-1)=hold(-1,-1) - sa*dsin(beta - betap)*sap
     &     *(ca*cb*cbp*cap + ca*cap*sb*sbp + sa*sap)*ci*A4
      hold(1,1)=hold(1,1)
     &  - 1/8.0d0*(4.0d0 - dcos(2.0d0*(alpha - alphap) + beta - betap)
     &  - dcos(2.0d0*(alpha - alphap) - beta + betap)
     &  - 2.0d0*(dcos(2.0d0*(alpha - alphap))+dcos(2.0d0*(alpha + alphap)))
     &  + dcos(2.0d0*(alpha + alphap) + beta - betap)
     &  + dcos(2.0d0*(alpha + alphap) - beta + betap))*dsin(beta-betap)*ci*A5
      hold(1,-1)=hold(1,-1) + 1/8.0d0*(dcos((beta + betap)/2.0d0)
     &  - ci*dsin((beta + betap)/2.0d0))*
     &  (dcos((4.0d0*alpha + 3.0d0*(beta - betap) + 4.0d0*alphap)/2.0d0)
     &  - dcos((4.0d0*alpha - 3.0d0*(beta - betap) + 4.0d0*alphap)/2.0d0)
     &  - 3.0d0*dcos((4.0d0*alpha + beta - betap + 4.0d0*alphap)/2.0d0)
     &  + 3.0d0*dcos((4.0d0*alpha - beta + betap + 4.0d0*alphap)/2.0d0)
     &  - ci*(dsin((4.0d0*alpha + 3.0d0*(beta - betap) - 4.0d0*alphap)/2.0d0)
     &  + dsin((4.0d0*alpha - 3.0d0*(beta - betap) - 4.0d0*alphap)/2.0d0)
     &  + 3.0d0*dsin((4.0d0*alpha + beta - betap - 4.0d0*alphap)/2.0d0)
     &  + 3.0d0*dsin((4.0d0*alpha - beta + betap - 4.0d0*alphap)/2.0d0)))*ci*A5
      hold(-1,-1)=hold(-1,-1) + 1/8.0d0*(4.0d0 - dcos(2.0d0*alpha + beta
     &   - betap - 2.0d0*alphap) - dcos(2.0d0*alpha - beta
     &   + betap - 2.0d0*alphap) - 2.0d0*(dcos(2.0d0*(alpha - alphap))
     &   + dcos(2.0d0*(alpha + alphap))) + dcos(2.0d0*alpha + beta
     &   - betap + 2.0d0*alphap) + dcos(2.0d0*alpha - beta + betap
     &   + 2.0d0*alphap))*dsin(beta-betap)*ci*A5
      hold(1,-1)=hold(1,-1)
     &   + ((ci*cb + sb)*(-cb*cbp*cap*sa - cap*sa*sb*sbp + ca*sap)
     &   - (ci*cbp + sbp)*(cap*sa - ca*cb*cbp*sap - ca*sb*sbp*sap))*ci*A6
c
c     The operator which multiplies A6 only contributes to the spin-flip amplitude
c     Now take the complex conjugate to get hold(-1,1)
      hold(-1,1)=-(Real(hold(1,-1)) - ci*Imag(hold(1,-1)))

      if (verbosity.ge.3) then
         write(*,*) "   trafo1Nampxx"
         write(*,*) "      hold(1,1) = ",hold(1,1)
         write(*,*) "      hold(-1,1) = ",hold(-1,1)
         write(*,*) "      hold(1,-1) = ",hold(1,-1)
         write(*,*) "      hold(-1,-1) = ",hold(-1,-1)
      end if
      
      return
      end
c
c**********************************************************************
c**********************************************************************
c

      subroutine trafo1Nampyx(hold,A1,A2,A3,A4,A5,A6,
     &     alpha,beta,alphap,betap,verbosity)
c
c**********************************************************************
c   Calculates the single-nucleon c.m. frame amplitude for given
c     "invariant"  amplitudes A1-A6 and incoming photon angles
c     (alpha,beta), outgoing (alphap,betap). (For notation see BKM review.)
c
c   Initial photon polarization=y, final photon polarization=x
c
c**********************************************************************
c doCU SEE SUBROUTINE trafo1Nampxx
c**********************************************************************
c
      implicit none
c
      include '../common-densities/constants.def'
c
c======================================================================
c
c  OUTPUT:
c
      complex*16,intent(out) :: hold(-1:1,-1:1)
c
c======================================================================
c
c  INPUT:
c
      real*8,intent(in) :: A1,A2,A3,A4,A5,A6
      real*8,intent(in) :: alpha,beta,alphap,betap
      integer,intent(in) :: verbosity
      
c======================================================================
c
c  LOCAL:
c
      real*8 ca,cap,sa,sap,cb,cbp,sb,sbp
c
c**********************************************************************
c
      if (verbosity.eq.1000) continue
      hold=c0
      ca=dcos(alpha)
      cap=dcos(alphap)
      sa=dsin(alpha)
      sap=dsin(alphap)
      cb=dcos(beta)
      cbp=dcos(betap)
      sb=dsin(beta)
      sbp=dsin(betap)
c     hgrie May 2018: flipped overall signs in following expressions -- see note on top of file
      hold(1,1)=-cap*dsin(beta-betap)*A1
     & + sap*dsin(beta-betap)*(ca*sap - cap*cb*cbp*sa - cap*sa*sb*sbp)*A2
      hold(-1,-1)=hold(1,1)
      hold(1,-1)=c0
      hold(1,1)=hold(1,1) + cap*dcos(beta-betap)*ci*A3
      hold(1,-1)=hold(1,-1) + sap*(cb - ci*sb)*ci*A3
      hold(-1,-1)=hold(-1,-1) - cap*dcos(beta-betap)*ci*A3
      hold(1,1)=hold(1,1) - sa*dsin(2.0d0*alphap)*dsin(beta-betap)**2
     &     *1/2.0d0*ci*A4
      hold(1,-1)=hold(1,-1) - cap*dsin(beta-betap)*(ca*sap*sbp - cap*sa*sb
     &     -ci*(cap*cb*sa - ca*cbp*sap))*ci*A4
      hold(-1,-1)=hold(-1,-1) + sa*dsin(2.0d0*alphap)*dsin(beta-betap)**2
     &     *1/2.0d0*ci*A4
      hold(1,1)=hold(1,1) - (sa*dsin(2.0d0*alphap)*dsin(beta-betap)**2/2.0d0
     & + dcos(beta-betap)*sap
     &       *(ca*sap - cap*cb*cbp*sa - cap*sa*sb*sbp))*ci*A5
      hold(1,-1)=hold(1,-1) - (sap*dsin(beta-betap)*(ci*(ca*cap*cbp + cb*sa*sap)
     & + sa*sap*sb + ca*cap*sbp) + cap*(cap*cb*cbp*sa - ca*sap + cap*sa*sb*sbp)
     & *(cb - ci*sb))*ci*A5
      hold(-1,-1)=hold(-1,-1) + (sa*dsin(2.0d0*alphap)*dsin(beta-betap)**2/2.0d0
     & + dcos(beta-betap)*sap*
     &     (ca*sap - cap*cb*cbp*sa - cap*sa*sb*sbp))*ci*A5
      hold(1,1)=hold(1,1) -
     &     sa*(ca*sap - cap*cb*cbp*sa - cap*sa*sb*sbp)*ci*A6
      hold(1,-1)=hold(1,-1) - (sap*dsin(beta-betap)*(ci*cbp + sbp)
     &     + ca*(cb - ci*sb)*
     &     (cap*cb*cbp*sa - ca*sap + cap*sa*sb*sbp))*ci*A6
      hold(-1,-1)=hold(-1,-1)
     &     + sa*(ca*sap - cap*cb*cbp*sa - cap*sa*sb*sbp)*ci*A6
c Now take the complex conjugate to get hold(-1,1)
      hold(-1,1)=-(Real(hold(1,-1)) - ci*Imag(hold(1,-1)))
      
      if (verbosity.ge.3) then
         write(*,*) "   trafo1Nampyx"
         write(*,*) "      hold(1,1) = ",hold(1,1)
         write(*,*) "      hold(-1,1) = ",hold(-1,1)
         write(*,*) "      hold(1,-1) = ",hold(1,-1)
         write(*,*) "      hold(-1,-1) = ",hold(-1,-1)
      end if
      
      return
      end
c
c**********************************************************************
c**********************************************************************
c
      subroutine trafo1Nampxy(hold,A1,A2,A3,A4,A5,A6,
     &     alpha,beta,alphap,betap,verbosity)
c
c**********************************************************************
c   Calculates the single-nucleon c.m. frame amplitude for given
c     "invariant"  amplitudes A1-A6 and incoming photon angles
c     (alpha,beta), outgoing (alphap,betap). (For notation see BKM review.)
c
c   Initial photon polarization=x, final photon polarization=y
c
c**********************************************************************
c doCU SEE SUBROUTINE trafo1Nampxx
c**********************************************************************
c
c      implicit none
c
      include '../common-densities/constants.def'
c
c======================================================================
c
c  OUTPUT:
c
      complex*16,intent(out) :: hold(-1:1,-1:1)
      
c======================================================================
c
c  INPUT:
c
      real*8,intent(in) :: A1,A2,A3,A4,A5,A6
      real*8,intent(in) :: alpha,beta,alphap,betap
      integer,intent(in) :: verbosity
c
c======================================================================
c
c  LOCAL:
c
      real*8 ca,cap,sa,sap,cb,cbp,sb,sbp
      
c**********************************************************************
      if (verbosity.eq.1000) continue
      hold=c0
      ca=dcos(alpha)
      cap=dcos(alphap)
      sa=dsin(alpha)
      sap=dsin(alphap)
      cb=dcos(beta)
      cbp=dcos(betap)
      sb=dsin(beta)
      sbp=dsin(betap)
c     hgrie May 2018: flipped overall signs in following expressions -- see note on top of file
      hold(1,1)=-ca*dsin(betap-beta)*A1
     & - sa*dsin(betap-beta)*(cap*sa - ca*cb*cbp*sap - ca*sap*sb*sbp)*A2
      hold(-1,-1)=hold(1,1)
      hold(1,-1)=c0
      hold(1,1)=hold(1,1) - ca*dcos(beta-betap)*ci*A3
      hold(1,-1)=hold(1,-1) - sa*(cbp - ci*sbp)*ci*A3
      hold(-1,-1)=hold(-1,-1) + ca*dcos(beta-betap)*ci*A3
      hold(1,1)=hold(1,1) + sap*dsin(2.0d0*alpha)*dsin(beta-betap)**2
     &     *1/2.0d0*ci*A4
      hold(1,-1)=hold(1,-1) + ca*dsin(betap-beta)*(cap*sa*sb - ca*sap*sbp
     &     -ci*(ca*cbp*sap - cap*cb*sa))*ci*A4
      hold(-1,-1)=hold(-1,-1) - sap*dsin(2.0d0*alpha)*dsin(beta-betap)**2
     &     *1/2.0d0*ci*A4
      hold(1,1)=hold(1,1) + (sap*dsin(2.0d0*alpha)*dsin(beta-betap)**2/2.0d0
     & + dcos(beta-betap)
     &     *sa*(cap*sa - ca*cb*cbp*sap - ca*sap*sb*sbp))*ci*A5
      hold(1,-1)=hold(1,-1) - (-sa*dsin(betap-beta)*(ci*(ca*cap*cb + cbp*sa*sap)
     & + sa*sap*sbp + ca*cap*sb) + ca*(-ca*cb*cbp*sap + cap*sa - ca*sap*sb*sbp)
     & *(cbp - ci*sbp))*ci*A5
      hold(-1,-1)=hold(-1,-1) - (sap*dsin(2.0d0*alpha)*dsin(beta-betap)**2/2.0d0
     & + dcos(beta-betap)
     &     *sa*(cap*sa - ca*cb*cbp*sap - ca*sap*sb*sbp))*ci*A5
      hold(1,1)=hold(1,1) +
     &     sap*(cap*sa - ca*cb*cbp*sap - ca*sap*sb*sbp)*ci*A6
      hold(1,-1)=hold(1,-1) - (sa*dsin(beta-betap)*(ci*cb + sb)
     & + cap*(cbp-ci*sbp)
     &     *(-ca*cb*cbp*sap + cap*sa - ca*sap*sb*sbp))*ci*A6
      hold(-1,-1)=hold(-1,-1)
     &     - sap*(cap*sa - ca*cb*cbp*sap - ca*sap*sb*sbp)*ci*A6
c Now take the complex conjugate to get hold(-1,1)
      hold(-1,1)=-(Real(hold(1,-1)) - ci*Imag(hold(1,-1)))

      if (verbosity.ge.3) then
         write(*,*) "   trafo1Nampxy"
         write(*,*) "      hold(1,1) = ",hold(1,1)
         write(*,*) "      hold(-1,1) = ",hold(-1,1)
         write(*,*) "      hold(1,-1) = ",hold(1,-1)
         write(*,*) "      hold(-1,-1) = ",hold(-1,-1)
      end if
      
      return
      end
c
c**********************************************************************
c**********************************************************************
c
      subroutine trafo1Nampyy(hold,A1,A2,A3,A4,A5,A6,
     &     alpha,beta,alphap,betap,verbosity)
c
c**********************************************************************
c   Calculates the single-nucleon c.m. frame amplitude for given
c     "invariant"  amplitudes A1-A6 and incoming photon angles
c     (alpha,beta), outgoing (alphap,betap). (For notation see BKM review.)
c
c   Initial photon polarization=y, final photon polarization=y
c
c**********************************************************************
c doCU SEE SUBROUTINE trafo1Nampxx
c**********************************************************************
c
      implicit none
c
      include '../common-densities/constants.def'
c
c======================================================================
c
c  OUTPUT:
c
      complex*16,intent(out) :: hold(-1:1,-1:1)
      
c======================================================================
c
c  INPUT:
c
      real*8,intent(in) :: A1,A2,A3,A4,A5,A6
      real*8,intent(in) :: alpha,beta,alphap,betap
      integer,intent(in) :: verbosity
      
c======================================================================
c
c  LOCAL:
c
      real*8 ca,cap,sa,sap,cb,cbp,sb,sbp
      
c**********************************************************************
      if (verbosity.eq.1000) continue
      hold=c0
      ca=dcos(alpha)
      cap=dcos(alphap)
      sa=dsin(alpha)
      sap=dsin(alphap)
      cb=dcos(beta)
      cbp=dcos(betap)
      sb=dsin(beta)
      sbp=dsin(betap)
      hold(1,1)=dcos(beta-betap)*A1 - sa*sap*dsin(beta-betap)**2*A2
      hold(-1,-1)=hold(1,1)
      hold(1,-1)=c0
      hold(1,1)=hold(1,1) + dsin(beta-betap)*ci*A3
      hold(-1,-1)=hold(-1,-1) - dsin(beta-betap)*ci*A3
c
c  No contribution from A1->A3 to spin-flip amplitude
c
      hold(1,1)=hold(1,1) + sa*sap*dsin(2.0d0*(beta-betap))*1/2.0d0*ci*A4
      hold(1,-1)=hold(1,-1) + dcos(beta-betap)*(ca*sap*sbp - cap*sa*sb
     &     - ci*(cap*cb*sa - ca*cbp*sap))*ci*A4
      hold(-1,-1)=hold(-1,-1) - sa*sap*dsin(2.0d0*(beta-betap))*1/2.0d0*ci*A4
      hold(1,1)=hold(1,1) + sa*sap*dsin(2.0d0*(beta-betap))*ci*A5
      hold(1,-1)=hold(1,-1) - (ca*sap*(cbp - ci*sbp) + cap*sa*(cb - ci*sb))
     &     *dsin(beta-betap)*ci*A5
      hold(-1,-1)=hold(-1,-1) - sa*sap*dsin(2.0d0*(beta-betap))*ci*A5
      hold(1,1)=hold(1,1) - (dcos(2.0d0*alpha) + dcos(2.0d0*alphap) - 2.0d0)
     &     *dsin(beta-betap)*1/2.0d0*ci*A6
      hold(1,-1)=hold(1,-1) - dsin(beta-betap)*(dsin(2.0d0*alpha)*(cb - ci*sb)
     &      + dsin(2.0d0*alphap)*(cbp - ci*sbp))*1/2.0d0*ci*A6
      hold(-1,-1)=hold(-1,-1) + (dcos(2.0d0*alpha) + dcos(2.0d0*alphap) - 2.0d0)
     &     *dsin(beta-betap)*1/2.0d0*ci*A6
c     Now take the complex conjugate to get hold(-1,1)
      hold(-1,1)=-(Real(hold(1,-1)) - ci*Imag(hold(1,-1)))
      
      if (verbosity.ge.3) then
         write(*,*) "   trafo1Nampyy"
         write(*,*) "      hold(1,1) = ",hold(1,1)
         write(*,*) "      hold(-1,1) = ",hold(-1,1)
         write(*,*) "      hold(1,-1) = ",hold(1,-1)
         write(*,*) "      hold(-1,-1) = ",hold(-1,-1)
      end if
      
      return
      end
