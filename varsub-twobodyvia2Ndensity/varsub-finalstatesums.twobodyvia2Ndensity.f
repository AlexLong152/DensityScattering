cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              twobodyfinalstatesumsvia2Ndensity : perform
c                  Σ_(mt12,j12,s12,l12,m12) Σ_(mt12p=mt12,j12p,s12p,l12p,m12p,Mzp,Mz)
c                  ∫ dp12 p12² ∫ dp12p p12p²/(2π)³ (HC)³ * rho(p12,p12p) * KERNEL
c              with internal subroutine
c                  Calculate2BIntegralI2 : Σ_(msp,ms) ∫ dΩ12 ∫ dΩ12p of kernel*rho, for ALL channel pairs
c                                          at fixed (p12,p12p)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, based on finalstatesums.twobodyvia2Ndensity.f of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code.
c           No back-compatibility
c     v2.0 Nov 2026 Alexander Long: restructured for speed. Same result up to round-off, for ANY kernel.
c           Previously, main looped over (mt12,j12,s12,l12,m12,p12), this routine over (j12p,s12p,l12p,m12p,p12p,Mzp,Mz),
c           and Calculate2BIntegralI2 over (msp,ms) and then over the angles. So at every integration point
c           (p12,Ω12,p12p,Ω12p), the kernel was called again with IDENTICAL arguments for every (Mzp,Mz), every
c           (j12,j12p) and every (m12,ms),(m12p,msp) with the same ml12=m12-ms, ml12p=m12p-msp -- none of
c           (Mzp,Mz,j12,m12,j12p,m12p) is an argument of the kernel.
c           Now the loop order is inverted: integration points outside, channels inside. At each integration point,
c           the kernel is called once for each distinct set of its arguments [see "kappa" below]; Y_l12p(p12p-hat) and
c           the bracketing of |p12p| on the density grid are computed once per distinct p12p returned by the kernel;
c           the interpolation of rho in p12 is done once per p12 (it is fixed during the inner loops). Only the
c           products of these with the Clebsches are formed per channel pair and (Mzp,Mz). Channel pairs with no
c           density entry are skipped, as are kernel values which are exactly zero at a point.
c           The sum over all channels, formerly in main, is therefore now done here.
c           ASSUMPTION ON KERNEL -- the only one: Calc2Bspinisospintrans() is a function of its arguments, i.e. same
c           arguments give same Kernel2B and ppVecs. EVERY argument is treated as something the kernel may depend
c           on, including l12,ml12,l12p,ml12p which no kernel uses as of Nov 2026. So kernels of other processes and
c           of higher orders (e.g. with t12/=t12p pieces as in Compton at Odelta4) need no change here.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c
c     twoMz dependence: arrays and do-loops
c
c     Defined twobody ME to INCLUDE the factor 1/(2π)³, so that final amplitudes for onebody and twobody have SAME sizes.
c
c     Note structure of loops is almost the same as in one-body case, only real difference is in computation of I2, and
c     fact that p_{12}' integral now runs over full range [0,infty). Partly for this reason, the order of the p_{12}'
c     and p_3' loops has been interchanged, c.f. the one-body version of this routine.
c
c     hgrie 20 June 2014: modified for use of LebedevLaikov or Gaussian
c     integration for theta & phi separately, for solid angle integral in (12) system
c
c     Units: p12, p12p in fm⁻¹ in the mantle; converted to MeV (pAbs, uAbs) before the kernel is called.
c     Note P12MAG(ip12p) is the magnitude of the generic integration variable uVec, not of the physical momentum p12p.
c     The physical p12p (ppVecs) is returned by the kernel (variable substitution).
c
c     Bookkeeping used below:
c        "state"  : (12) quantum numbers |(l12 s12) j12 m12, t12> without mt12 -- same set for in- and out-state
c        "pair"   : (mt12, in-state, out-state) with at least one (Mzp,Mz) for which the density has an entry
c        "half"   : distinct (s12,l12,ml12) [t12 follows from l12,s12] -- what the kernel is told about one side
c        "kappa"  : distinct full sets of quantum-number arguments of the kernel: (mt12, in-half, out-half), i.e.
c                   (ml12,ml12p,t12,mt12,t12p,mt12p=mt12,l12,s12,l12p,s12p). The kernel is called once per kappa
c                   and integration point; the block Kernel2B(:,:,s12p,msp,s12,ms), all (msp,ms), is kept.
c                   Pair (in-state,out-state) with (ms,msp) uses the kappa with ml12=m12-ms, ml12p=m12p-msp,
c                   and of it the element (msp,ms) -- exactly the kernel call the former code made for it.
c        "slot"   : distinct physical p12p vector returned by the kernel at this integration point (ppVecs).
c                   Usually one, but each kappa and each diagram may come with its own.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine twobodyfinalstatesumsvia2Ndensity(
     &     Result,Mnucl,
     &     Anucl,twoSnucl,extQnumlimit,
     &     Eprobe,thetacm,
     &     P12MAG,AP12MAG,NP12,
     &     th12,phi12,Nth12,Nphi12,j12max,
     &     AngularType12,angweight12,calctype,symmetry,numDiagrams,verbosity)
c
      USE clebsch
      USE CompDens ! needs module CompDens.mod
c     CompDens makes available the following variables (along with some others not listed here)
c     rho, rhoDensity, P12P_density, P12W_density
      implicit NONE
c
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
      integer,intent(in) :: j12max
      real*8,intent(in)  :: Eprobe,thetacm
      real*8,intent(in)  :: P12MAG(Npmax),AP12MAG(Npmax)
      integer,intent(in) :: NP12
      real*8,intent(in)  :: th12(Nangmax),phi12(Nangmax)
      integer,intent(in) :: Nth12,Nphi12
      integer,intent(in) :: Anucl,twoSnucl
      integer,intent(in) :: extQnumlimit, numDiagrams

      integer,intent(in) :: AngularType12
      real*8,intent(in)  :: angweight12(Nangmax,Nangmax)
      integer,intent(in) :: calctype,symmetry,verbosity
      real*8, intent(in) :: Mnucl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES: contributions are ADDED to Result

      complex*16,intent(inout) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
c
c     states: quantum numbers; spin projections ms allowed by |m12-ms|<=l12 (Nvms of them, in vms); Clebsches
      integer Nst,ist,istp
      integer,allocatable :: stJ(:),stS(:),stL(:),stT(:),stM(:)
      integer,allocatable :: Nvms(:),vms(:,:)
      real*8,allocatable  :: cgcSt(:,:)        ! (ms,ist) = CG(l12 ml12, s12 ms | j12 m12)
c     halfs: quantum numbers; stHf(ms,ist) = half of state ist with spin projection ms (0 if ms not allowed)
      integer Nhf,ihf,ihfp
      integer,allocatable :: hfS(:),hfL(:),hfMl(:),hfT(:),stHf(:,:)
c     kappas: mt12, in-half kapHu, out-half kapHp; kapIdx(in-half,out-half,mt12) = kappa (0 if not needed)
      integer Nkap
      integer,allocatable :: kapMt(:),kapHu(:),kapHp(:),kapIdx(:,:,:)
c     pairs: in-state prU, out-state prP, kappa prKap(msp,ms,ipair); prNz valid (Mzp,Mz) with density slot prSlot
      integer Npair,maxpair,NMz,nz
      integer,allocatable :: prU(:),prP(:),prKap(:,:,:),prNz(:)
      integer,allocatable :: prMzp(:,:),prMz(:,:),prSlot(:,:)
c     kernel at one integration point: Kblk(idiag,extQnum,msp,ms,ikap) = Kernel2B(idiag,extQnum,s12p,msp,s12,ms);
c     kapNZ: not all zero; kapSlot(idiag,ikap): slot of its ppVecs
      complex*16,allocatable :: Kblk(:,:,:,:,:)
      logical,allocatable :: kapNZ(:)
      integer,allocatable :: kapSlot(:,:)
c     slots: ppVec; bracket and weights of |ppVec| on density grid;
c            YPst(msp,istp,slot) = CG * conjg( Y_l12p^(m12p-msp)(ppVec-hat) )
      integer maxslot
      real*8,allocatable :: slotPP(:,:),slotWy1(:),slotWy2(:)
      integer,allocatable :: slotLoc(:),slotLoc2(:)
      complex*16,allocatable :: YPst(:,:,:)
c     density: slices rhoDensity(:,:,rindx) used (rUsed), map rindx->slot (rSlot);
c              rhoP(:,slot): rho interpolated in p12 at the current p12, as function of p12p on the density grid
      integer Nr,ir,Nrho,Ndens,jdens
      integer,allocatable :: rUsed(:),rSlot(:)
      real*8,allocatable  :: rhoP(:,:),P12_MeV(:)
c     angular grid, flattened to one index for Gaussian (θ,φ) and LebedevLaikov alike
      integer Nang,ia,ith,iphi
      real*8,allocatable :: angTh(:),angPhi(:),angW(:)
c     Yu(ms,ist,ia) = CG * Y_l12^(m12-ms)(Ω12(ia)) -- independent of p12, so computed once
      complex*16,allocatable :: Yu(:,:,:)
      complex*16 Yl12(-5:5)
c     momenta: pVecs(:,ia) physical p12 (fixed p12), uVecs(:,ib) integration variable (fixed p12p), in MeV
      real*8,allocatable :: pVecs(:,:),uVecs(:,:)
      integer ip12,ip12p
      real*8 pAbs,uAbs,wx1,wx2
      integer locp,loc2p
      complex*16 Int2B(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16 fact
c
      integer mt12,j12,s12,l12,t12,m12,ms,msp,k,kp,lmax
      integer twoMz,twoMzp,alpha2N,alpha2Np,rindx

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     states of (12) subsystem
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      lmax=j12max+1
      if (lmax.gt.5) then
         write(*,*) "*** ERROR: Y_lm only implemented up to l=5, i.e. j12max<=4. -- Exiting."
         stop
      end if
      call initclebsch                ! Initializing the factorial array

      Nst=0
      do j12=0,j12max
         do s12=0,1
            do l12=abs(j12-s12),j12+s12
               Nst=Nst+2*j12+1
            end do
         end do
      end do
      allocate(stJ(Nst),stS(Nst),stL(Nst),stT(Nst),stM(Nst))
      allocate(Nvms(Nst),vms(3,Nst),cgcSt(-1:1,Nst),stHf(-1:1,Nst))
      allocate(hfS(3*Nst),hfL(3*Nst),hfMl(3*Nst),hfT(3*Nst))
      Nhf=0
      stHf=0
      ist=0
      do j12=0,j12max                      ! total ang mom (12); usually Jmax=1 for 1% convergence
         do s12=0,1                        ! spin (12)
            do l12=abs(j12-s12),j12+s12    ! angular mom. (12)
               t12=(1-(-1)**(l12+s12+1))/2 ! isospin (12)
               do m12=-j12,j12             ! spin projection (12)
                  ist=ist+1
                  stJ(ist)=j12
                  stS(ist)=s12
                  stL(ist)=l12
                  stT(ist)=t12
                  stM(ist)=m12
                  Nvms(ist)=0
                  cgcSt(:,ist)=0.d0
                  do ms=-s12,s12
                     if (abs(m12-ms).le.l12) then
                        Nvms(ist)=Nvms(ist)+1
                        vms(Nvms(ist),ist)=ms
                        cgcSt(ms,ist)=CG(2*l12,2*s12,2*j12,2*(m12-ms),2*ms,2*m12)
c     find half (s12,l12,ml12=m12-ms), or add new one
                        do ihf=1,Nhf
                           if ((hfS(ihf).eq.s12).and.(hfL(ihf).eq.l12).and.(hfMl(ihf).eq.m12-ms)) then
                              stHf(ms,ist)=ihf
                              exit
                           end if
                        end do
                        if (stHf(ms,ist).eq.0) then
                           Nhf=Nhf+1
                           hfS(Nhf)=s12
                           hfL(Nhf)=l12
                           hfMl(Nhf)=m12-ms
                           hfT(Nhf)=t12
                           stHf(ms,ist)=Nhf
                        end if
                     end if
                  end do
               end do
            end do
         end do
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     channel pairs, kappas, and density slices used
c     mt12: only charged pion exchange at OQ4 => (12) subsystem is (pn) -- but all mt12 kept,
c     kernel decides; pairs with zero kernel are skipped at run time.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      NMz=(twoSnucl+1)**2
      maxpair=3*Nst*Nst
      Nrho=size(rhoDensity,3)
      allocate(prU(maxpair),prP(maxpair),prKap(-1:1,-1:1,maxpair),prNz(maxpair))
      allocate(prMzp(NMz,maxpair),prMz(NMz,maxpair),prSlot(NMz,maxpair))
      allocate(kapMt(3*Nhf*Nhf),kapHu(3*Nhf*Nhf),kapHp(3*Nhf*Nhf),kapIdx(Nhf,Nhf,-1:1))
      allocate(rUsed(maxpair*NMz),rSlot(Nrho))
      kapIdx=0
      prKap=0
      rSlot=0
      Npair=0
      Nkap=0
      Nr=0
      do mt12=-1,1
         do ist=1,Nst
            do istp=1,Nst       ! mt12p=mt12
               if ((Nvms(ist).eq.0).or.(Nvms(istp).eq.0)) cycle
               nz=0
               do twoMzp=twoSnucl,-twoSnucl,-2
                  do twoMz=twoSnucl,-twoSnucl,-2
                     alpha2N = get2Nchannum(stL(ist),stS(ist),stJ(ist),mt12,stM(ist),twoMz)
                     alpha2Np = get2Nchannum(stL(istp),stS(istp),stJ(istp),mt12,stM(istp),twoMzp)
                     if (alpha2N.le.0 .or. alpha2Np.le.0) then
                        rindx = 0
                     else
                        rindx=rhoindx(alpha2N,alpha2Np)
                     end if
                     if (rindx.lt.1 .or. rindx.gt.Nrho) cycle ! not in density: contributes nothing
                     if (rSlot(rindx).eq.0) then
                        Nr=Nr+1
                        rUsed(Nr)=rindx
                        rSlot(rindx)=Nr
                     end if
                     nz=nz+1
                     prMzp(nz,Npair+1)=twoMzp
                     prMz(nz,Npair+1)=twoMz
                     prSlot(nz,Npair+1)=rSlot(rindx)
                  end do
               end do
               if (nz.eq.0) cycle
               Npair=Npair+1
               prU(Npair)=ist
               prP(Npair)=istp
               prNz(Npair)=nz
c     kappas needed by this pair: one per (ms,msp); add those not yet there
               do kp=1,Nvms(istp)
                  msp=vms(kp,istp)
                  ihfp=stHf(msp,istp)
                  do k=1,Nvms(ist)
                     ms=vms(k,ist)
                     ihf=stHf(ms,ist)
                     if (kapIdx(ihf,ihfp,mt12).eq.0) then
                        Nkap=Nkap+1
                        kapMt(Nkap)=mt12
                        kapHu(Nkap)=ihf
                        kapHp(Nkap)=ihfp
                        kapIdx(ihf,ihfp,mt12)=Nkap
                     end if
                     prKap(msp,ms,Npair)=kapIdx(ihf,ihfp,mt12)
                  end do
               end do
            end do
         end do
      end do
      if (verbosity.ge.1) then
         write(*,'(A,I6,A,I4,A,I6,A)') "    Convolution: ",Npair," channel pairs with density entries, ",
     &        Nkap," kernel calls per integration point, ",Nr," density slices."
      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     angular grid, flattened; Yu; density grid
c     hgrie 20 June 2014: pick theta& phi summation parameters following AngularType12
c     for LebedevLaikov, only sum over diagonal elements of angweight12 (all others are zero)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (AngularType12.eq.1) then !Gaussian in theta and phi separately
         Nang=Nth12*Nphi12
      else if (AngularType12.eq.2) then !LebedevLaikov
         Nang=Nth12
      else
         write(*,*) "*** ERROR: Something went wrong with AngularType12 in twobodyfinalstatesums. -- Exiting."
         stop
      end if
      allocate(angTh(Nang),angPhi(Nang),angW(Nang))
      ia=0
      do ith=1,Nth12
         if (AngularType12.eq.1) then
            do iphi=1,Nphi12
               ia=ia+1
               angTh(ia)=th12(ith)
               angPhi(ia)=phi12(iphi)
               angW(ia)=angweight12(ith,iphi)
            end do
         else
            ia=ia+1
            angTh(ia)=th12(ith)
            angPhi(ia)=phi12(ith)
            angW(ia)=angweight12(ith,ith)
         end if
      end do

      allocate(Yu(-1:1,Nst,Nang))
      Yu=c0
      do ia=1,Nang
         do ist=1,Nst
            call getsphericalharmonics(Yl12,stL(ist),angTh(ia),angPhi(ia))
            do k=1,Nvms(ist)
               ms=vms(k,ist)
               Yu(ms,ist,ia)=cgcSt(ms,ist)*Yl12(stM(ist)-ms)
            end do
         end do
      end do

c     P12P_density units in  fm^-1, so convert to MeV
      Ndens=size(P12P_density)
      if ((size(rhoDensity,1).ne.Ndens).or.(size(rhoDensity,2).ne.Ndens)) then
         write(*,*) "*** ERROR: rhoDensity not on (P12P_density x P12P_density) grid. -- Exiting."
         stop
      end if
      allocate(P12_MeV(Ndens),rhoP(Ndens,max(Nr,1)))
      P12_MeV=P12P_density*HC
      allocate(pVecs(3,Nang),uVecs(3,Nang))
      maxslot=max(numDiagrams*Nkap,1)
      allocate(Kblk(1:numDiagrams,1:extQnumlimit,-1:1,-1:1,max(Nkap,1)),kapNZ(max(Nkap,1)))
      allocate(kapSlot(numDiagrams,max(Nkap,1)))
      allocate(slotPP(3,maxslot),slotWy1(maxslot),slotWy2(maxslot),slotLoc(maxslot),slotLoc2(maxslot))
      allocate(YPst(-1:1,Nst,maxslot))

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     radial integrations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do ip12=1,NP12 ! mag of momentum (12) subsystem
         pAbs=P12MAG(ip12)*HC   !pAbs in units of MeV
c     rho is linear in p12 between grid points: interpolate all slices used to pAbs once
         call densitybracket(P12_MeV,Ndens,pAbs,locp,loc2p,wx1,wx2)
         do ir=1,Nr
            do jdens=1,Ndens
               rhoP(jdens,ir)=wx1*real(rhoDensity(locp,jdens,rUsed(ir)),8)
     &              +wx2*real(rhoDensity(loc2p,jdens,rUsed(ir)),8)
            end do
         end do
         do ia=1,Nang
            call sphere2cart(pVecs(1,ia),pVecs(2,ia),pVecs(3,ia),pAbs,angTh(ia),angPhi(ia),verbosity)
         end do

         do ip12p=1,NP12 ! mag of generic integration variable uVec
            uAbs=P12MAG(ip12p)*HC
            do ia=1,Nang
               call sphere2cart(uVecs(1,ia),uVecs(2,ia),uVecs(3,ia),uAbs,angTh(ia),angPhi(ia),verbosity)
            end do

            call Calculate2BIntegralI2(Int2B)

c   not including explcit (ppAbs/HC)**2 so lose two powers of HC
            fact=(Anucl*(Anucl-1)/2)*(P12MAG(ip12)**2)*AP12MAG(ip12)*AP12MAG(ip12p)*(P12MAG(ip12p)**2)/(2*Pi)**3*
     &           HC**3.d0
            Result = Result + fact*Int2B
         end do                 !ip12p
      end do                    !ip12

      deallocate(stJ,stS,stL,stT,stM,Nvms,vms,cgcSt,stHf,hfS,hfL,hfMl,hfT)
      deallocate(prU,prP,prKap,prNz,prMzp,prMz,prSlot)
      deallocate(kapMt,kapHu,kapHp,kapIdx,Kblk,kapNZ,kapSlot)
      deallocate(slotPP,slotWy1,slotWy2,slotLoc,slotLoc2,YPst)
      deallocate(rUsed,rSlot,rhoP,P12_MeV,angTh,angPhi,angW,Yu,pVecs,uVecs)

      if (symmetry.eq.1000) continue
      return

      contains

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Σ_(msp,ms) ∫ dΩ12 ∫ dΩ12p of Y_l12(p12-hat) Y*_l12p(p12p-hat) * kernel * rho, times Clebsches,
c     for all channel pairs and (Mzp,Mz), at fixed p12 (pVecs, rhoP) and p12p (uVecs).
c     Output Int2B(extQnum,Mzp,Mz), summed over diagrams and channel pairs.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Calculate2BIntegralI2(Int2B)

      complex*16,intent(out) :: Int2B(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)

      complex*16 Kernel2B(1:numDiagrams,1:extQnumlimit,0:1,-1:1,0:1,-1:1)
      real*8 ppVecs(1:numDiagrams,1:3),ppAbs
      complex*16 Yl12p(-5:5,0:5)
c     T(:,j): Σ_(msp,ms) of CG CGp Y Y* Kernel, for those (msp,ms) of a pair whose kappa has slot tSlot(j)
      complex*16 T(1:extQnumlimit,9),c
      integer tSlot(9),Nt,jt
      real*8 wab,rhoval,f
      integer Nslot,islot,jslot
      integer ia,ib,idiag,s,sp,ms,msp,k,kp,iz,extQnum,l
      integer ipair,ikap,ist,istp,ir,ihf,ihfp

      Int2B=c0
      do ia=1,Nang              ! angle integral: Ω of p12
         do ib=1,Nang           ! angle integral: Ω of p12p (in fact of generic integration variable uVec)
            wab=angW(ia)*angW(ib)
c     kernel for all kappas, each with exactly the arguments the former code used for it;
c     Y_l12p and density bracket for each distinct p12p returned
            Nslot=0
            do ikap=1,Nkap
               ihf=kapHu(ikap)
               ihfp=kapHp(ikap)
               s=hfS(ihf)
               sp=hfS(ihfp)
               ppVecs=0.d0
               call Calc2Bspinisospintrans(Kernel2B,ppVecs,Mnucl,
     &              extQnumlimit,hfMl(ihf),hfMl(ihfp),
     &              hfT(ihf),kapMt(ikap),hfT(ihfp),kapMt(ikap),hfL(ihf),
     &              s,hfL(ihfp),sp,thetacm,Eprobe,pVecs(:,ia),
     &              uVecs(:,ib),numDiagrams,calctype,verbosity)
               Kblk(:,:,:,:,ikap)=Kernel2B(:,:,sp,:,s,:)
               kapNZ(ikap)=any(Kblk(:,:,:,:,ikap).ne.c0)
               if (.not.kapNZ(ikap)) cycle ! contributes exactly zero at this integration point
               do idiag=1,numDiagrams
                  kapSlot(idiag,ikap)=0
                  do jslot=1,Nslot
                     if (all(ppVecs(idiag,:).eq.slotPP(:,jslot))) then
                        kapSlot(idiag,ikap)=jslot
                        exit
                     end if
                  end do
                  if (kapSlot(idiag,ikap).ne.0) cycle
                  Nslot=Nslot+1
                  kapSlot(idiag,ikap)=Nslot
                  slotPP(:,Nslot)=ppVecs(idiag,:)
                  do l=0,lmax
                     call getHarmonicCart(Yl12p(:,l),l,ppVecs(idiag,:),verbosity)
                  end do
                  ppAbs = sqrt(DOT_PRODUCT(ppVecs(idiag,:),ppVecs(idiag,:)))!in MeV
                  call densitybracket(P12_MeV,Ndens,ppAbs,slotLoc(Nslot),slotLoc2(Nslot),
     &                 slotWy1(Nslot),slotWy2(Nslot))
                  do istp=1,Nst
                     do k=1,Nvms(istp)
                        msp=vms(k,istp)
                        YPst(msp,istp,Nslot)=cgcSt(msp,istp)*conjg(Yl12p(stM(istp)-msp,stL(istp)))
                     end do
                  end do
               end do           !idiag
            end do              !ikap
c     sum over channel pairs
            do ipair=1,Npair
               ist=prU(ipair)
               istp=prP(ipair)
               do idiag=1,numDiagrams
c     T = Σ_(msp,ms) CG CGp Y_l12^ml12(p12-hat) Y*_l12p^ml12p(p12p-hat) Kernel(s12p,msp,s12,ms),
c     kept apart for each p12p (slot) since rho(p12,p12p) differs
                  Nt=0
                  do kp=1,Nvms(istp)
                     msp=vms(kp,istp)
                     do k=1,Nvms(ist)
                        ms=vms(k,ist)
                        ikap=prKap(msp,ms,ipair)
                        if (.not.kapNZ(ikap)) cycle
                        islot=kapSlot(idiag,ikap)
                        do jt=1,Nt
                           if (tSlot(jt).eq.islot) exit
                        end do
                        if (jt.gt.Nt) then
                           Nt=jt
                           tSlot(jt)=islot
                           T(:,jt)=c0
                        end if
                        c=Yu(ms,ist,ia)*YPst(msp,istp,islot)
                        do extQnum=1,extQnumlimit
                           T(extQnum,jt)=T(extQnum,jt)+c*Kblk(idiag,extQnum,msp,ms,ikap)
                        end do
                     end do
                  end do
c     times rho(p12,p12p), bilinear interpolation (in p12 already done in rhoP)
                  do jt=1,Nt
                     islot=tSlot(jt)
                     do iz=1,prNz(ipair)
                        ir=prSlot(iz,ipair)
                        rhoval=slotWy1(islot)*rhoP(slotLoc(islot),ir)+slotWy2(islot)*rhoP(slotLoc2(islot),ir)
                        f=wab*rhoval
c                                   Would need ppAbs**2 here if it wasn't in fact below already
                        do extQnum=1,extQnumlimit
                           Int2B(extQnum,prMzp(iz,ipair),prMz(iz,ipair)) =
     &                          Int2B(extQnum,prMzp(iz,ipair),prMz(iz,ipair)) + f*T(extQnum,jt)
                        end do
                     end do
                  end do
               end do           !idiag
            end do              !ipair
         end do                 !ib
      end do                    !ia

      if (any(Int2B.ne.Int2B)) then
         write(*,*) "Int array is nan somewhere, stopping"
         write(*,*) "   at P12MAG(ip12) =",P12MAG(ip12)," P12MAG(ip12p) =",P12MAG(ip12p)
         error stop
      end if
      end subroutine Calculate2BIntegralI2

      end subroutine twobodyfinalstatesumsvia2Ndensity
