cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              densitybracket  : bracket |p| on the density grid and return linear-interpolation weights
c              getHarmonicCart : Y_l of the direction of a Cartesian vector
c
c     Calculate2BIntegralI2 (Σ_(msp,ms) ∫ dΩ12 ∫ dΩ12p of kernel) is now an internal subroutine of
c     twobodyfinalstatesumsvia2Ndensity in finalstatesums.twobodyvia2Ndensity.f -- see there.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, near-identical to calculate2BI2.f of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility 
c     v2.0 Nov 2026 Alexander Long: Calculate2BIntegralI2 moved into finalstatesums.twobodyvia2Ndensity.f and
c           restructured for speed (see there). subroutine interpolate() replaced by densitybracket():
c           the bilinear interpolation of rho is now split into p12 part (once per p12) and p12p part,
c           with the same arithmetic, so rho values are bit-identical to before.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine densitybracket(P12_MeV,Nlength,pAbs,loc,loc2,w1,w2)
c     Find loc, loc2=loc+1 with P12_MeV(loc) <= pAbs < P12_MeV(loc2) (binary search, O(log N)),
c     and the weights of linear interpolation in pAbs between these points:
c         f(pAbs) = w1*f(loc) + w2*f(loc2),   w1=(x2-pAbs)/(x2-x1), w2=(pAbs-x1)/(x2-x1)
c     Outside the grid, the first/last interval is used, i.e. linear extrapolation.
c     Bilinear interpolation of rho(p,pp) as in former subroutine interpolate(), with the same arithmetic:
c         rho = wpp1*(wp1*rho(locp,locpp)+wp2*rho(loc2p,locpp)) + wpp2*(wp1*rho(locp,loc2pp)+wp2*rho(loc2p,loc2pp))
c     Units of all momenta in MeV. Assumes P12_MeV is sorted in ascending order.
      implicit none
      integer,intent(in) :: Nlength
      real*8,intent(in)  :: P12_MeV(Nlength),pAbs
      integer,intent(out):: loc,loc2
      real*8,intent(out) :: w1,w2
      integer lo,hi,mid
      real*8 x1,x2

      if (pAbs .le. P12_MeV(1)) then
         loc = 1
         loc2 = 2
      else if (pAbs .ge. P12_MeV(Nlength)) then
         loc = Nlength - 1
         loc2 = Nlength
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
         loc = lo
         loc2 = hi
      end if
      x1=P12_MeV(loc)
      x2=P12_MeV(loc2)
      w1=(x2-pAbs)/(x2-x1)
      w2=(pAbs-x1)/(x2-x1)
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
      real*8 vec(3)
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
      vec=ppVec
      if(all(vec .eq. 0.d0)) then
        vec(3)=0.001
      end if
      call cart2Sphere(vec(1),vec(2),vec(3),
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
