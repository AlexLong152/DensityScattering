ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              Calculate2BIntegralI2 :  Σ_(msp,ms) ∫ dΩ12 ∫ dΩ12p of kernel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, near-identical to calculate2BI2.f of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Calculate2BIntegralI2(Int2B,Mnucl,
c     &     Int2Bx,Int2By,Int2Bpx,Int2Bpy, ! for STUMP, see below
     &     extQnumlimit,
     &     j12p,m12p,l12p,s12p,t12p,mt12p,j12,m12,
     &     l12,s12,t12,mt12,p12,p12p,th12,phi12,Nth12,Nphi12,
     &     thetacm,k,
     &     AngularType12,angweight12,calctype,verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      USE clebsch
      implicit none
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
      
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: m12p,m12,j12p,s12p,l12p,j12,s12,l12,Nth12,Nphi12
      integer,intent(in) :: t12p,t12,mt12p,mt12
      
      real*8,intent(in) :: thetacm,k,th12(Nangmax),phi12(Nangmax)
      real*8, intent(in) :: Mnucl
      integer,intent(in) :: AngularType12
      real*8,intent(in) :: angweight12(Nangmax,Nangmax)
      integer,intent(in) :: calctype,verbosity
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES:

      complex*16,intent(out) :: Int2B(1:extQnumlimit)
      
c      complex*16,intent(out) :: Int2Bx,Int2By,Int2Bpx,Int2Bpy ! for STUMP, see below
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
      
c     index limits of iphi, depending on AngularType12:
      integer imin,imax,jmin,jmax
      
      complex*16 Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1) ! was Compton2Bxx/xy/yx/yy
      
c      complex*16 Compton2Bx(0:1,-1:1,0:1,-1:1)  ! for STUMP, see below
c      complex*16 Compton2By(0:1,-1:1,0:1,-1:1)  ! for STUMP, see below
c      complex*16 Compton2Bpx(0:1,-1:1,0:1,-1:1) ! for STUMP, see below
c      complex*16 Compton2Bpy(0:1,-1:1,0:1,-1:1) ! for STUMP, see below
     
      integer extQnum           ! counter of combined external quantum numbers of in and out state
      
      integer ith,iphi,jth,jphi,msp,ms,ml12p,ml12
      complex*16 Yl12(-5:5),Yl12p(-5:5)
      complex*16 Int(1:extQnumlimit,-5:5,-5:5)
c      complex*16 Intx(-5:5,-5:5),Inty(-5:5,-5:5)   ! for STUMP, see below
c      complex*16 Intpx(-5:5,-5:5),Intpy(-5:5,-5:5) ! for STUMP, see below
      complex*16 Yl12pstar
      real*8 cgcp,cgc,p12x,p12y,p12z,p12px,p12py,p12pz,p12,p12p
      real*8 offset
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      if (verbosity.eq.1000) continue ! unused variable, kept for future use
c     
      if ((l12p .gt. 5) .or. (l12 .gt. 5)) then
         goto 100
      endif 
      Int2B=c0
      Kernel2B=c0
      
c      Int2Bx=c0      ! for STUMP, see below
c      Int2By=c0      ! for STUMP, see below
c      Int2Bpx=c0     ! for STUMP, see below
c      Int2Bpy=c0     ! for STUMP, see below
c      Compton2Bx=c0  ! for STUMP, see below
c      Compton2By=c0  ! for STUMP, see below
c      Compton2Bpx=c0 ! for STUMP, see below
c      Compton2Bpy=c0 ! for STUMP, see below
     
      call initclebsch                ! Initializing the factorial array
c     Loop  to sum over the spin projections of the (12) system: ms12 and ms12p (called ms and msp here). 
c     The value of ms and msp together with m12 & m12p determine ml12 and ml12p. 
      do msp=-s12p,s12p,1
         ml12p=m12p-msp
         do ms=-s12,s12,1
            ml12=m12-ms
c     Initializing to zero
            cgc=0.d0
            cgcp=0.d0
            Yl12=c0
            Yl12p=c0
            Yl12pstar=c0
            Int=c0
c            Intx=c0  ! for STUMP, see below
c            Inty=c0  ! for STUMP, see below
c            Intpx=c0  ! for STUMP, see below
c            Intpy=c0  ! for STUMP, see below
            if ((abs(ml12p) .le. l12p) .and. (abs(ml12) .le. l12)) then     
c     angle integral: θ of p12
               do ith=1,Nth12
c     c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c     hgrie 20 June 2014: pick theta& phi summation parameters following AngularType12
c     for LebedevLaikov, only sum over diagonal elements of angweight12 (all others are zero)
                  if (AngularType12.eq.1) then !Gaussian in theta and phi separately
                     imin=1
                     imax=Nphi12
                  else if (AngularType12.eq.2) then !LebedevLaikov
                     imin=ith
                     imax=ith
                  else
                     write(*,*) "*** ERROR: Something went wrong with imin/imax in Calculate2BI2. -- Exiting."
                     stop
                  end if    
c     angle integral: φ of p12
                  do iphi=imin,imax
                     offset=10.d0*sqrt(2.d0)*PI/180
c                    rotate by an irrational number of degrees to avoid singularities

c     Inputs are theta and phi, outputs are x,y,z
                    
                     call CalculatePVector(p12x,p12y,p12z,p12,
     &                    mod(th12(ith)+offset,PI),mod(phi12(iphi)+offset,2*PI),verbosity)
                     call getsphericalharmonics(Yl12,l12,mod(th12(ith)+offset,PI),mod(phi12(iphi)+offset,2.d0*PI))
c     angle integral: θprime of p12p
                     do jth=1,Nth12
c     c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c     hgrie 20 June 2014: pick theta& phi summation parameters following AngularType12
c     for LebedevLaikov, only sum over diagonal elements of angweight12 (all others are zero)
                        if (AngularType12.eq.1) then !Gaussian in theta and phi separately
                           jmin=1
                           jmax=Nphi12
                        else if (AngularType12.eq.2) then !LebedevLaikov
                           jmin=jth
                           jmax=jth
                        else
                           write(*,*) "*** ERROR: Something went wrong with imin/imax in Calculate2BI2. -- Exiting."
                           stop
                        end if    
c     angle integral: φprime of p12p
                        do jphi=jmin,jmax
                           call CalculatePVector(p12px,p12py,p12pz,p12p,
     &                          th12(jth),phi12(jphi),verbosity)
                           call getsphericalharmonics(Yl12p,l12p,th12(jth),phi12(jphi))
                           Yl12pstar=Real(Yl12p(ml12p))-ci*Imag(Yl12p(ml12p))
                           call Calc2Bspinisospintrans(Kernel2B,Mnucl,
     &                          extQnumlimit,
     &                          t12,mt12,t12p,mt12p,l12,
     &                          s12,l12p,s12p,thetacm,k,p12x,p12y,p12z,
     &                          p12px,p12py,p12pz,calctype,verbosity)
                           
                           do extQnum=1,extQnumlimit
                              Int(extQnum,ml12p,ml12) = Int(extQnum,ml12p,ml12) + Yl12pstar*Yl12(ml12)*
     &                             angweight12(ith,iphi)*angweight12(jth,jphi)*Kernel2B(extQnum,s12p,msp,s12,ms)
                           end do   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie Nov 2023: Following is a STUMP from the Compton code, used there only for OQ4 -- NOT YET IMPLEMENTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     I leave this here because maybe some of this can be recycled later for boost corrections or so?
c                           Intx(ml12p,ml12)=Intx(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
c     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
c     &                          Compton2Bx(s12p,msp,s12,ms)
c                           Inty(ml12p,ml12)=Inty(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
c     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
c     &                          Compton2By(s12p,msp,s12,ms)
c                           Intpx(ml12p,ml12)=Intpx(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
c     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
c     &                          Compton2Bpx(s12p,msp,s12,ms)
c                           Intpy(ml12p,ml12)=Intpy(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
c     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
c     &                          Compton2Bpy(s12p,msp,s12,ms)
c     END OF STUMP    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
                        end do  ! jphi
                     end do     ! jth
                  end do        ! iphi
               end do           ! ith
            end if
c     Clebsches for unprimed and primed quantum numbers
            cgc=CG(2*l12,2*s12,2*j12,2*ml12,2*ms,2*m12)
            cgcp=CG(2*l12p,2*s12p,2*j12p,2*ml12p,2*msp,2*m12p)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
c     hgrie Nov 2023: Following is a STUMP from the Compton code, used there only for OQ4 -- NOT YET IMPLEMENTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     I leave this here because maybe some of this can be recycled later for boost corrections or so?
c            Int2Bxx=Int2Bxx+Intxx(ml12p,ml12)*cgc*cgcp
c            Int2Bxy=Int2Bxy+Intxy(ml12p,ml12)*cgc*cgcp
c            Int2Byx=Int2Byx+Intyx(ml12p,ml12)*cgc*cgcp
c            Int2Byy=Int2Byy+Intyy(ml12p,ml12)*cgc*cgcp
c            Int2Bx=Int2Bx+Intx(ml12p,ml12)*cgc*cgcp
c            Int2By=Int2By+Inty(ml12p,ml12)*cgc*cgcp
c            Int2Bpx=Int2Bpx+Intpx(ml12p,ml12)*cgc*cgcp
c            Int2Bpy=Int2Bpy+Intpy(ml12p,ml12)*cgc*cgcp
c     END OF STUMP    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            do extQnum=1,extQnumlimit
               Int2B(extQnum) = Int2B(extQnum) + Int(extQnum,ml12p,ml12)*cgc*cgcp
            end do
         end do                 !ms12
      end do                    !ms12p
      if (verbosity.eq.1000) continue
 100  return
      end
