ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              Setquad12 : Set up the quadratures for the radial, theta, & phi integrations of the (12) systems
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, identical to setquads.f of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     Set up the quadratures for the radial, theta, & phi integrations of the (12) systems:
c     use of LebedevLaikov or Gaussian integration for theta & phi separately,
c     for solid angle integral in (12) system 
c     combined wth*wphi*sin(th) (weight of angles theta and phi) into
c     one array angweight12(,) -- so the sum of all weights is 4\pi.
c      
c     twoSmax/twoMz dependence: none
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Setquad12(th12,Nth12,phi12,Nphi12,
     &     Nordth12,Nthbins12,Nordphi12,Nphibins12,
     &     AngularType12,angweight12,Nanggrid12,verbosity)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
      implicit none
      include '../common-densities/params.def'
      include '../common-densities/constants.def'
      include '../common-densities/calctype.def'
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES:
c     
c     Nth-total number of theta quadratures
c     Nphi12-total number of phi quadratures
c     angweight12: combined weight of angular integrations, _including_ sin^2(theta)
c     so sum of all weights is 4*Pi
c     
      real*8,intent(out)  :: th12(Nangmax)
      real*8,intent(out)  :: phi12(Nangmax)
      integer,intent(out) :: Nth12,Nphi12
      real*8,intent(out)  :: angweight12(Nangmax,Nangmax)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
c     
c     Nordp,Npbins-number of quadratures/bin and number of bins for
c     radial integration
c     Radialtype-determines nature of mapping from [0,infty] to [0,1]
c     Nordth,Nthbins-number of quadratures/bin and number of bins 
c     for theta integration
c     Nordphi,Nphibins-number of quadratures/bin and number of bins 
c     for phi integration
c     densitytype-used to determine whether we need to cut off the p
c     quadratures
c     hgrie 20 June 2014:
c     AngularType12-determines nature of angular integration in (12) system
c     Nanggrid12: number of points on solid angle grid of Lebedev-Laikov grid
c     
      integer,intent(in) :: Nordth12,Nthbins12,Nordphi12,Nphibins12

      integer,intent(in) :: AngularType12,Nanggrid12

      integer,intent(in) :: verbosity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
      
c     th12,wth-theta quadratures and weights, set up on [0,PI]
c     phi12,wphi-phi quadratures and weights, set up on [0,2 PI]

      integer ith,iphi
      real*8 wth12(Nangmax)
      real*8 wphi(Nangmax)
      
c     angwgth: local variable passed from LebedevLaikov routine 
      real*8 angwgth(Nangmax)
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Set up  two-dimensional angular integration mesh
c     CalculateIntegralI2 ((12)-integration) needs theta & phi mesh
c     Fill array with ZEROES
      angweight12=0.0E0
      if (AngularType12.eq.1) then ! Gauss-Legendre integration in angles
         call AnglePtsWts(Nordth12,Nthbins12,Nangmax,0.d0,PI,th12,wth12,Nth12,verbosity)
         call AnglePtsWts(Nordphi12,Nphibins12,Nangmax,0.d0,2.0d0*PI,phi12,wphi,Nphi12,verbosity)
         do ith=1,Nth12
            do iphi=1,Nphi12
               angweight12(ith,iphi) = wth12(ith) * wphi(iphi) * dsin(th12(ith))
            end do
         end do
      else if (AngularType12.eq.2) then ! LebedevLaikov integration in angles
         call LebedevLaikovWts(th12,phi12,angwgth,Nanggrid12,verbosity)
         do ith=1,Nanggrid12
            angweight12(ith,ith) = angwgth(ith)
         end do
c     be brutal and just re-define Nth12 to be Nanggrid12
         Nth12 = Nanggrid12
      end if   
c     
      if (verbosity.eq.1) then  ! sum should be 4\pi, for total solid angle
         write(*,*) "Sum of all angular weights (expect 4Pi) = ",SUM(angweight12)/(4*Pi),"*4Pi"
      endif
      if (verbosity.eq.1000) continue
c     end hgrie mod       
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
