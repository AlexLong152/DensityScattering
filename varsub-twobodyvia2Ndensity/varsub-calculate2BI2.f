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
c     subroutine Calculate2BIntegralI2(Int2B,Mnucl,
c    &     extQnumlimit,
c    &     j12p,m12p,l12p,s12p,t12p,mt12p,j12,m12,
c    &     l12,s12,t12,mt12,pAbs,uVecRs,th12,phi12,Nth12,Nphi12,NP12,
c    &     thetacm,Eprobe,
c    &     AngularType12,angweight12,calctype,numDiagrams,ip12p,twoSnucl,twoMzp,twoMz,verbosity)

      subroutine Calculate2BIntegralI2(Int2B,ppVecs,Mnucl,
     &    extQnumlimit,
     &    j12p,m12p,l12p,s12p,t12p,mt12p,j12,m12,
     &    l12,s12,t12,mt12, pAbs,uVecRs,th12,
     &    phi12,Nth12,Nphi12,NP12,thetacm,Eprobe,
     &    AngularType12,angweight12,calctype,numDiagrams,ip12p,twoSnucl,twoMzp,twoMz,verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      USE clebsch
      USE CompDens
      implicit none
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
      
      integer,intent(in) :: extQnumlimit, numDiagrams,ip12p,twoSnucl
      integer,intent(in) :: m12p,m12,j12p,s12p,l12p,j12,s12,l12,Nth12,Nphi12
      integer,intent(in) :: t12p,t12,mt12p,mt12
      
      real*8,intent(in) :: thetacm,th12(Nangmax),phi12(Nangmax)
      real*8, intent(in) :: Mnucl, Eprobe
      integer,intent(in) :: AngularType12
      real*8,intent(in) :: angweight12(Nangmax,Nangmax)
      integer,intent(in) :: calctype,verbosity
      real*8 ppVecs(1:numDiagrams,1:3),ppAbs
      real*8 tmpRho 
      integer twoMz,twoMzp
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
      integer NP12
      integer ith,iphi,jth,jphi,msp,ms,ml12p,ml12
      complex*16 Yl12(-5:5),Yl12p(-5:5)
      complex*16 Int(1:numDiagrams,1:extQnumlimit,-5:5,-5:5)
c      complex*16 Intx(-5:5,-5:5),Inty(-5:5,-5:5)   ! for STUMP, see below
c      complex*16 Intpx(-5:5,-5:5),Intpy(-5:5,-5:5) ! for STUMP, see below
      complex*16 Yl12pstar
      real*8 cgcp,cgc,p12x,p12y,p12z,pAbs,uVecRs(NP12)
      real*8 ux,uy,uz
      real*8 pVec(3), uVec(3)
      real*8,allocatable :: P12_MeV(:)
      integer alpha2N, alpha2Np, rindx

c     write(*,*) "numDiagrams=", numDiagrams 
c     write(*,*) "extQnumlimit=", extQnumlimit 
      if (.not.allocated(P12_MeV)) then
          allocate(P12_MeV, mold=P12P_density)
      end if
      if (verbosity.eq.1000) continue ! unused variable, kept for future use
c     
      if ((l12p .gt. 5) .or. (l12 .gt. 5)) then
         goto 100
      endif 
      Int2B=c0
      Kernel2B=c0
      
      alpha2N = get2Nchannum(l12,s12,j12,mt12,m12,twoMz)
      alpha2Np = get2Nchannum(l12p,s12p,j12p,mt12p,m12p,twoMzp)
      rindx=rhoindx(alpha2N,alpha2Np)

      call initclebsch                ! Initializing the factorial array
c     Loop  to sum over the spin projections of the (12) system: ms12 and ms12p (called ms and msp here). 
c     The value of ms and msp together with m12 & m12p determine ml12 and ml12p. 
      P12_MeV=P12P_density*HC
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
                    
                     call sphere2cart(p12x,p12y,p12z,pAbs,
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
                           call sphere2cart(ux,uy,uz,uVecRs(ip12p),
     &                          th12(jth),phi12(jphi),verbosity)!uVecR used to be p12p

                           uVec=(/ux,uy,uz/)!generic integration variable

                           ppVecs=0.d0
c                          write(*,*) "varsub-calculate2BI2.f:173 numDiagrams,calctype,verbosity=", numDiagrams,calctype,verbosity 
                           call Calc2Bspinisospintrans(Kernel2B,ppVecs,Mnucl,
     &                          extQnumlimit,ml12,ml12p,
     &                          t12,mt12,t12p,mt12p,l12,
     &                          s12,l12p,s12p,thetacm,Eprobe,pVec,
     &                          uVec,numDiagrams,calctype,verbosity)
c                         write(*,*) "ppVecs(1,:)=", ppVecs(1,:) 
c                         tmpRho=0.0000001
                          do diagNum=1,numDiagrams
c                             radVec=(/uVecR,th12(jth),phi12(jphi)/)
                              call getHarmonicCart(Yl12p,l12p,ppVecs(diagNum,:),verbosity)
                              Yl12pstar=Real(Yl12p(ml12p))-ci*Imag(Yl12p(ml12p))

                              ppAbs = sqrt(DOT_PRODUCT(ppVecs(diagNum,:),ppVecs(diagNum,:)))!in MeV
c                             write(*,*) "In varsub-calculate2BI2: ppAbs=",ppAbs 
c                             write(*,*) "In varsub-calculate2BI2: ppVecs(diagNum,:)=",ppVecs(diagNum,:) 
c                             stop
                              call interpolate(tmpRho, real(rhoDensity(:,:,rindx),8), P12_MeV, ppAbs, pAbs, size(P12_MeV))

                              do extQnum=1,extQnumlimit
                                   Int(diagNum,extQnum,ml12p,ml12) = Int(diagNum,extQnum,ml12p,ml12)+Yl12(ml12)*Yl12pstar*
     &                              angweight12(ith,iphi)*angweight12(jth,jphi)*Kernel2B(diagNum,extQnum,s12p,msp,s12,ms)*
     &                                  tmpRho
c                                   Would need ppAbs**2 here if it wasn't in varsub-finalstatesums.f already
                                end do!extQnum   
                           end do!diagNum
                        end do  ! jphi
                     end do     ! jth
                  end do        ! iphi
               end do           ! ith
            end if
c     Clebsches for unprimed and primed quantum numbers
            cgc=CG(2*l12,2*s12,2*j12,2*ml12,2*ms,2*m12)
            cgcp=CG(2*l12p,2*s12p,2*j12p,2*ml12p,2*msp,2*m12p)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            do diagNum=1,numDiagrams
            do extQnum=1,extQnumlimit
               Int2B(diagNum,extQnum) = Int2B(diagNum,extQnum) + Int(diagNum,extQnum,ml12p,ml12)*cgc*cgcp
            end do
            end do
         end do                 !ms12
      end do                    !ms12p

        
c     do extQnum=1,extQnumlimit
c     do diagNum=1,numDiagrams
c        ppAbs= sqrt(DOT_PRODUCT(ppVecs(diagNum,:),ppVecs(diagNum,:)))!in MeV
c        call interpolate(tmpRho, real(rhoDensity(:,:,rindx),8), P12_MeV, ppAbs, pAbs, size(P12_MeV))
c        Int2B(diagNum,extQnum)=Int2B(diagNum,extQnum)*tmpRho
c     end do !diagNum
c     end do !extQnum
      
      if (verbosity.eq.1000) continue
 100  return
      end subroutine


      subroutine interpolate(tmpRho,rho,P12_MeV, ppAbs, pAbs, Nlength)
c     Units of all momenta in MeV
c     Uses inlined binary search O(log N) to find bracketing indices
c     Assumes P12_MeV is sorted in ascending order
      implicit none
      real*8 tmpRho,ppAbs,pAbs,P12_MeV(Nlength)
      integer Nlength
      real*8 rho(Nlength,Nlength)
      real*8 fQ11, fQ12,fQ21,fQ22
      integer locp, loc2p, locpp, loc2pp, lo, hi, mid
      real*8 x1,x2,y1,y2,term1,term2

c     Inlined binary search for pAbs bracket
      if (pAbs .le. P12_MeV(1)) then
         locp = 1
         loc2p = 2
      else if (pAbs .ge. P12_MeV(Nlength)) then
         locp = Nlength - 1
         loc2p = Nlength
      else
         lo = 1
         hi = Nlength
         do while (hi - lo .gt. 1)
            mid = (lo + hi) / 2
            if (P12_MeV(mid) .le. pAbs) then
               lo = mid
            else
               hi = mid
            end if
         end do
         locp = lo
         loc2p = hi
      end if

c     Inlined binary search for ppAbs bracket
      if (ppAbs .le. P12_MeV(1)) then
         locpp = 1
         loc2pp = 2
      else if (ppAbs .ge. P12_MeV(Nlength)) then
         locpp = Nlength - 1
         loc2pp = Nlength
      else
         lo = 1
         hi = Nlength
         do while (hi - lo .gt. 1)
            mid = (lo + hi) / 2
            if (P12_MeV(mid) .le. ppAbs) then
               lo = mid
            else
               hi = mid
            end if
         end do
         locpp = lo
         loc2pp = hi
      end if

      x1=P12_MeV(locp)
      x2=P12_MeV(loc2p)

      y1=P12_MeV(locpp)
      y2=P12_MeV(loc2pp)

c https://en.wikipedia.org/wiki/Bilinear_interpolation
c fQ11 =f(x1, y1), fQ12 =f(x1, y2), fQ21 =f(x2, y1), fQ22 =f(x2, y2)

      fQ11=rho(locp,locpp)
      fQ12=rho(locp,loc2pp)
      fQ21=rho(loc2p,locpp)
      fQ22=rho(loc2p,loc2pp)

c     Inlined bilinear interpolation
      term1=((y2-ppAbs)/(y2-y1))
      term1=term1*(((x2-pAbs)/(x2-x1))*fQ11+(((pAbs-x1)/(x2-x1))*fQ21))

      term2=((ppAbs-y1)/(y2-y1))
      term2=term2*(
     &   ((x2-pAbs)/(x2-x1))*fQ12+
     &   ((pAbs-x1)/(x2-x1))*fQ22)
      tmpRho=term1+term2

      end subroutine

      subroutine getHarmonicCart(Yl,l12,ppVec,verbosity)
c     Gets the spherical harmonic Y for uVec which is given in cartiesian
c     useful for evaluating Y_l(\hat{uVec+offsetVec})
c     this is "trivial" but its easy to make a typo so just offload it to this subroutine

      implicit none
      include '../common-densities/constants.def'!can be removed once PI is not longer needed for debugging
c     Inputs ---------------------------------------
      real*8, intent(in) :: ppVec(3)
      integer, intent(in):: l12,verbosity
      

c     Outputs----------------------------------------
      complex*16, intent(out) ::  Yl(-5:5)

c     Internal --------------------------------------
      real*8 r,theta,phi

c   If there isn't a substitution then ppVecs(diagNum,:)==uVec so can check we get uVecR, th12(jth), phi12(jphi) in the
c   back transform which is what radVec and passFlag are for

c     To be deleted after debugging
c     real*8, intent(in) :: radVec(3)
c     Variables for debugging
c     complex*16 testYl(-5:5)
c     logical cond
c     logical passFlag
c     real*8 tmp(3)
c     integer ell

      call cart2Sphere(ppVec(1),ppVec(2),ppVec(3),
     &       r,theta,phi,verbosity)
      call getsphericalharmonics(Yl,l12,theta,phi)
    

c     tmp=abs(radVec-(/r,theta,phi/))/(2*PI)
c     cond=all(abs(tmp-nint(tmp)).le.1e-5)!nint gets closest whole number, dmod=double precision mod

c     if (cond) then
c         passFlag=.true.
c     else
c         write(*,*) "test array is", dmod(abs(radVec-(/r,theta,phi/)), 2*PI)
c         write(*,'(A)') "Test failed, as passed to spherical harmonic"
c         write(*,'(A,F12.5,",",F12.5",",F12.5,A)') "ppVec=",ppVec(1),ppVec(2),ppVec(3) ,"     : This is actually ppVecs(diagNum,:)"
c         write(*,'(A,F10.5,",",F10.5",",F10.5,A)') "r, theta, phi= ", radVec(1), radVec(2), radVec(3), " : radvec"
c         write(*,'(A,F10.5,",",F10.5",",F10.5,A)') "r, theta, phi= ", r, theta, phi, " : cart2sphere in getHarmonicCart"
c         write(*,'(A,F10.5,",",F10.5",",F10.5)') "sphere diffs = ", radVec(1)-r, radVec(2)-theta, radVec(3)-phi
c         call getsphericalharmonics(testYl,l12,radVec(2),radVec(3))
c         write(*,'(A)') "Spherical Harmonic diff is: "
c         do ell=-l12,l12
c             write(*,'(F7.5,SP,F7.5,"I, l=",I3)'), testYl(ell)-Yl(ell), ell
c         end do
c         passFlag=.false.
c     end if

c     if (any(Yl.ne.Yl))then
c         write(*,*) "In calculate2BI2.f: uVec=",uVec 
c         write(*,*) "Yl12pstar is NaN, exiting"
c         write(*,*) "In calculate2BI2.f: r=",r 
c         write(*,*) "In calculate2BI2.f: theta=",theta 
c         write(*,*) "In calculate2BI2.f: phi=",phi 
c         stop
c     end if
      return 
      end subroutine
