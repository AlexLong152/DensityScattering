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
      subroutine Calculate2BIntegralI2(Int2B,ppVecs,Mnucl,
     &     extQnumlimit,
     &     j12p,m12p,l12p,s12p,t12p,mt12p,j12,m12,
     &     l12,s12,t12,mt12,p12,p12p,th12,phi12,Nth12,Nphi12,
     &     thetacm,k,
     &     AngularType12,angweight12,calctype,numDiagrams,verbosity)
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
      
      integer,intent(in) :: extQnumlimit, numDiagrams
      integer,intent(in) :: m12p,m12,j12p,s12p,l12p,j12,s12,l12,Nth12,Nphi12
      integer,intent(in) :: t12p,t12,mt12p,mt12
      
      real*8,intent(in) :: thetacm,k,th12(Nangmax),phi12(Nangmax)
      real*8, intent(in) :: Mnucl
      integer,intent(in) :: AngularType12
      real*8,intent(in) :: angweight12(Nangmax,Nangmax)
      integer,intent(in) :: calctype,verbosity
      real*8 ppVecs(1:numDiagrams,1:3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES:

      complex*16,intent(out) :: Int2B(1:numDiagrams,1:extQnumlimit)
      
c      complex*16,intent(out) :: Int2Bx,Int2By,Int2Bpx,Int2Bpy ! for STUMP, see below
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
      
c     index limits of iphi, depending on AngularType12:
      integer imin,imax,jmin,jmax
      
      complex*16 Kernel2B(1:numDiagrams,1:extQnumlimit,0:1,-1:1,0:1,-1:1) ! was Compton2Bxx/xy/yx/yy
      
c      complex*16 Compton2Bx(0:1,-1:1,0:1,-1:1)  ! for STUMP, see below
c      complex*16 Compton2By(0:1,-1:1,0:1,-1:1)  ! for STUMP, see below
c      complex*16 Compton2Bpx(0:1,-1:1,0:1,-1:1) ! for STUMP, see below
c      complex*16 Compton2Bpy(0:1,-1:1,0:1,-1:1) ! for STUMP, see below
     
      integer extQnum,diagNum           ! counter of combined external quantum numbers of in and out state
      
      integer ith,iphi,jth,jphi,msp,ms,ml12p,ml12
      complex*16 Yl12(-5:5),Yl12p(-5:5)
      complex*16 Ylm(-5:5)!for testing
      complex*16 Int(1:numDiagrams,1:extQnumlimit,-5:5,-5:5)
c      complex*16 Intx(-5:5,-5:5),Inty(-5:5,-5:5)   ! for STUMP, see below
c      complex*16 Intpx(-5:5,-5:5),Intpy(-5:5,-5:5) ! for STUMP, see below
      complex*16 Yl12pstar
      real*8 cgcp,cgc,p12x,p12y,p12z,p12,p12p
      real*8 ux,uy,uz
      real*8 pVec(3), uVec(3)
      integer l,mPrint,lMax
      real*8 x,y,z
      real*8 vec(3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     TODO: try more values
c     lMax=2
c     x=2.d0
c     y=2.d0
c     z=2.d0
c     vec=(/x,y,z/)
c     do l=0,lMax
c         call getHarmonicCart(Ylm,l,vec,verbosity)
c         write(*,'(A,F5.3,A,F5.3,A,F5.3,A,I1)') " Ylm(x=",x,", y=",y,", z=",z,") l=",l

c         do mPrint=-l,l
c           write(*,'(A,I2,A,F7.3, " + ", F5.3, "i")')" m=",mPrint,":", real(Ylm(mPrint)), aimag(Ylm(mPrint))
c         end do
c         write(*,*) ""
c         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
c     end do
c     stop
      
      if (verbosity.eq.1000) continue ! unused variable, kept for future use
c     
      if ((l12p .gt. 5) .or. (l12 .gt. 5)) then
         goto 100
      endif 
      Int2B=c0
      Kernel2B=c0
      
c     write(*,*) "In calculate2BI2.f: extQnumlimit=",extQnumlimit 
c     write(*,*) "In calculate2BI2.f: numDiagrams=",numDiagrams 
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

c     Inputs are theta and phi, outputs are x,y,z
                    
                     call sphere2cart(p12x,p12y,p12z,p12,
     &                    th12(ith),phi12(iphi),verbosity)
                     pVec=(/p12x,p12y,p12z/)
                     call getsphericalharmonics(Yl12,l12,th12(ith),phi12(iphi))
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
                           call sphere2cart(ux,uy,uz,p12p,
     &                          th12(jth),phi12(jphi),verbosity)
                           uVec=(/ux,uy,uz/)!generic integration variable
c                          write(*,*) "diff of |uVec|-p12p should be zero, it is"
c                          write(*,*) sqrt(dot_product(uVec,uVec))-p12p

c                          call getsphericalharmonics(Yl12p,l12p,th12(jth),phi12(jphi))
c                          Yl12pstar=Real(Yl12p(ml12p))-ci*Imag(Yl12p(ml12p))
                        
                           call Calc2Bspinisospintrans(Kernel2B,ppVecs,Mnucl,
     &                          extQnumlimit,ml12,ml12p,
     &                          t12,mt12,t12p,mt12p,l12,
     &                          s12,l12p,s12p,thetacm,k,pVec,
     &                          uVec,calctype,numDiagrams,verbosity)

                          do diagNum=1,numDiagrams
                              call getHarmonicCart(Yl12p,l12p,ppVecs(diagNum,:),verbosity)
                              Yl12pstar=Real(Yl12p(ml12p))-ci*Imag(Yl12p(ml12p))

                               do extQnum=1,extQnumlimit
                                   Int(diagNum,extQnum,ml12p,ml12) = Int(diagNum,extQnum,ml12p,ml12)+Yl12(ml12)*Yl12pstar*
     &                                      angweight12(ith,iphi)*angweight12(jth,jphi)*Kernel2B(diagNum,extQnum,s12p,msp,s12,ms)
                                end do!extQnum   
                           end do!diagNum
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
            do diagNum=1,numDiagrams
            do extQnum=1,extQnumlimit
               Int2B(diagNum,extQnum) = Int2B(diagNum,extQnum) + Int(diagNum,extQnum,ml12p,ml12)*cgc*cgcp
            end do
            end do
         end do                 !ms12
      end do                    !ms12p
      if (verbosity.eq.1000) continue
 100  return
      end


      subroutine getHarmonicCart(Yl,l12,uVec,verbosity)
c     Gets the spherical harmonic Y for uVec which is given in cartiesian
c     useful for evaluating Y_l(\hat{uVec+offsetVec})
c     this is "trivial" but its easy to make a typo so just offload it to this subroutine

      implicit none
c     Inputs ---------------------------------------
      real*8, intent(in) :: uVec(3)
      integer, intent(in):: l12,verbosity

c     Outputs----------------------------------------
      complex*16, intent(out) ::  Yl(-5:5)

c     Internal --------------------------------------
      real*8 r,theta,phi

      call cart2Sphere(uVec(1),uVec(2),uVec(3),
     &       r,theta,phi,verbosity)
      
      call getsphericalharmonics(Yl,l12,theta,phi)
      if (any(Yl.ne.Yl))then
          write(*,*) "In calculate2BI2.f: uVec=",uVec 
          write(*,*) "Yl12pstar is NaN, exiting"
          write(*,*) "In calculate2BI2.f: r=",r 
          write(*,*) "In calculate2BI2.f: theta=",theta 
          write(*,*) "In calculate2BI2.f: phi=",phi 
          stop
c     else
c         write(*,*) "Yl12pstar wasnt nan"
      end if
      return 
      end subroutine
