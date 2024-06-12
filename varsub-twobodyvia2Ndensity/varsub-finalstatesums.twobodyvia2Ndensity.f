cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              twobodyfinalstatesumsvia2Ndensity : perform Σ_(mt12p=mt12,j12p,s12p,l12p,m12p,Mzp,Mz) ∫ dp12 p12² ∫ dp12p p12p²/(2π)³ (HC)³ * rho(p12,p12p) * KERNEL
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, based on finalstatesums.twobodyvia2Ndensity.f of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility 
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine twobodyfinalstatesumsvia2Ndensity(
     &     Result,Mnucl,
     &     Anucl,twoSnucl,extQnumlimit,j12,m12,l12,s12,t12,mt12,
     &     k,thetacm,
     &     ip12,p12,wp12,
     &     P12MAG,AP12MAG,NP12,
     &     th12,phi12,Nth12,Nphi12,j12max,
     &     AngularType12,angweight12,calctype,symmetry,numDiagrams,verbosity)
c     
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
      integer,intent(in) :: j12,m12,l12,s12,t12,mt12 ! are automatically integers, so do NOT multiply by 2, unlike for Mz=>twoMz
      integer,intent(in) :: j12max
      real*8,intent(in)  :: k,thetacm
      real*8,intent(in)  :: p12,wp12
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
c     OUTPUT VARIABLES:
      
      complex*16,intent(out) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
      
      integer alpha2N,alpha2Np,rindx
      integer extQnum           ! counter of combined external quantum numbers of in and out state
      
      integer mt12p,j12p,s12p,l12p,t12p,m12p ! quantum #s of (12) system -- integer already, so no factor 2
      integer ip12,ip12p
      integer diagNum
      complex*16 Int2B(1:numDiagrams,1:extQnumlimit)
      complex*16 fact
      
      real*8 ppVecs(1:numDiagrams,1:3)
c      complex*16 Int2Bx,Int2By,Int2Bpx,Int2Bpy ! for STUMP, see below
c      complex*16 Int3x, Int3y, Int3px, Int3py  ! for STUMP, see below
c      complex*16 factx, facty, factpx, factpy  ! for STUMP, see below
      
      integer twoMz,twoMzp
c     real*8, allocatable :: rhotmp(:,:,:)
c     real*4, allocatable :: rhotmp1(:,:,:)
c     real*8 getBaseArray
c      logical invokesymmetrytwoMzp  ! function to invoke symmetry/-ies of amplitudes, dependent on process
c      logical invokesymmetrytwoMz   ! function to invoke symmetry/-ies of amplitudes, dependent on process
c      logical invokesymmetryextQnum ! function to invoke symmetry/-ies of amplitudes, dependent on process
      real*8 tmpRho, ppAbs,pAbs
      real*8,allocatable :: diffsp(:), diffspp(:)
      real*8,allocatable :: P12_MeV(:)
      integer locp,loc2p,loc2pp, locpp!,loc
      real*8 x1,x2,y1,y2,fQ11,fQ12,fQ21,fQ22
c     real*8 switch1, switch2,switch

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mt12p=mt12                              ! isospin projection in (12), fixes charge of spectator(s)    
      do j12p=0,j12max                        ! total ang mom (12); usually j12max=1 for 1% convergence
         do s12p=0,1                          ! spin (12) subsystem
            do l12p=abs(j12p-s12p),j12p+s12p  ! orbital angular mom. (12)
               t12p=(1-(-1)**(l12p+s12p+1))/2 ! isospin (12) subsystem
               do m12p=-j12p,j12p             ! total ang mom projection of out-(12)
c                  
c     Angular-momentum sums are implemented exactly as in one-body version of this routine
c     
                  do ip12p=1,NP12 ! mag of momentum (12) subsystem
c
c                   write(*,*) "In finalstatesums.twobodyvia2Ndensity.f: P12MAG(ip12p)*HC=",P12MAG(ip12p)*HC 
                     call Calculate2BIntegralI2(Int2B,ppVecs,Mnucl,
     &                    extQnumlimit,
     &                    j12p,m12p,l12p,s12p,t12p,mt12p,
     &                    j12,m12,l12,s12,t12,mt12,
     &                    p12*HC,P12MAG(ip12p)*HC,th12,
     &                    phi12,Nth12,Nphi12,thetacm,k,
     &                    AngularType12,angweight12,calctype,numDiagrams,verbosity)
c     Note P12MAG(ip12p) is the generic integration variable, not physical momenta 

                     do twoMzp=twoSnucl,-twoSnucl,-2
c     check if ME to be skipped and reconstructed via a symmetry of the amplitudes
c                        if (invokesymmetrytwoMzp(symmetry,twoSnucl,twoMzp,verbosity)) exit
                        
                        do twoMz=twoSnucl,-twoSnucl,-2
c     check if ME to be skipped and reconstructed via a symmetry of the amplitudes
c                           if (invokesymmetrytwoMz(symmetry,twoSnucl,twoMzp,twoMz,verbosity)) exit
c                           
c     now call density
c                           write(*,*) l12,s12,j12,mt12,m12,twoMz
c                           write(*,*) l12p,s12p,j12p,mt12p,m12p,twoMzp
                           alpha2N = get2Nchannum(l12,s12,j12,mt12,m12,twoMz)
                           alpha2Np = get2Nchannum(l12p,s12p,j12p,mt12p,m12p,twoMzp)
                           rindx=rhoindx(alpha2N,alpha2Np)
c                           write(*,*) "rindx = ",rindx," ; alpha2N  = ",alpha2N," ; alpha2Np = ",alpha2Np
c                           write(*,*) "             ρ12 = ",rho(ip12,ip12,rindx)
c                           write(*,'(3(A,I5),A,E15.8,SP,E16.8," I")')
c     &                          "    ρ(rindx=",rindx,",ip12=",ip12,",ip12p=",ip12p,") = ",rho(ip12,ip12p,rindx)
                           
c     Next line: integration volume of p12 and p12p magnitudes and 2Ndensity ρ
c     hgrie May 2018/July 2020:
c     original factor 3.d0 in 3He-code is replaced by number of nucleon pairs inside nucleus -- see 3He-densities paper
c     hgrie Oct 2022: include here factor 1/(2π)³ so that final amplitudes for onebody and twobody have SAME sizes.
c     multiplication by HC**3.d0 transforms ME units from fm^-3 to MeV^3

c     TODO: check rhoDensity assignment doesn't fail in CompDens.F90 for large rho due to block reading etc. check line 1315
c     Note that P12MAG(ip12p) isn't actually physical momenta, but p12 is.
                            
c                          allocate(diffspp(size(P12P_density)))
                           if (.not.allocated(diffspp)) then
                                allocate(diffspp, mold=P12P_density)
                           end if 

                           if (.not.allocated(diffsp)) then
                               allocate(diffsp, mold=P12P_density)
                           end if 
                           if (.not.allocated(P12_MeV)) then
                               allocate(P12_MeV, mold=P12P_density)
                           end if
c                          allocate(diffsp(size(P12P_density)))
c                          allocate(P12_MeV(size(P12P_density)))

c P12P_density units in  fm^-1, so convert to MeV 
                           P12_MeV=P12P_density*HC
                           pAbs=p12*HC!pAbs in units of MeV now
c                          ppVecs=ppVecs        !ppVecs in units of MeV from 2Bkernel already

                           do diagNum=1,numDiagrams
                               ppAbs= sqrt(DOT_PRODUCT(ppVecs(diagNum,:),ppVecs(diagNum,:)))!in MeV

c   Find the closest values in the P12 array to ppAbs
                               diffspp=abs((P12_MeV-ppAbs))
                               locpp = minloc(diffspp,DIM=1) !closest value
                               loc2pp = minloc(diffspp,DIM=1,
     &                                      mask=.not.(diffspp.eq.diffspp(locpp)))! second closest value
c   Now do the same thing for p                                       
                               diffsp=abs((P12_MeV-pAbs))
                               locp = minloc(diffsp,DIM=1)
                               loc2p = minloc(diffsp,DIM=1,mask=.not.(diffsp.eq.diffsp(locp)))

                               x1=P12_MeV(locp)
                               x2=P12_MeV(loc2p)

                               y1=P12_MeV(locpp)
                               y2=P12_MeV(loc2pp)

c https://en.wikipedia.org/wiki/Bilinear_interpolation
c use the names Q to decrease the likelihood of a transcription error

                               fQ11=rhoDensity(locp,locpp,rindx)
                               fQ12=rhoDensity(locp,loc2pp,rindx)
                               fQ21=rhoDensity(loc2p,locpp,rindx)
                               fQ22=rhoDensity(loc2p,loc2pp,rindx)

                               call bilinear_interpolate(tmpRho, fQ11,fQ12,fQ21,fQ22,x1,x2,y1,y2,pAbs,ppAbs)
c                              if (tmpRho.ge.0.1) then
c                                  write(*,*) "In finalstatesums.twobodyvia2Ndensity.f"
c                                  write(*,*) "locp,locpp=",locp,loc2p,"-> x1,x2=",x1,x2
c                                  write(*,*) "loc2p,loc2pp=",locpp,loc2pp, "-> y1,y2=",y1,y2
c                                  write(*,*) "tmpRho, rho(ip12,ip12p,rindx)=",tmpRho, rho(ip12,ip12p,rindx)
c                                  write(*,*) "fQ21,fQ22", fQ21,fQ22
c                                  write(*,*) "fQ11,fQ12,=",fQ11,fQ12 
c                                  write(*,*) "x,y=",pAbs,ppAbs
c                                  stop
c                              end if

                               ppAbs=ppAbs/HC!TODO combine this into the HC**3 in fact
                               pAbs=pAbs/HC
                               fact=(Anucl*(Anucl-1)/2)*(pAbs**2)*wp12*(ppAbs**2)*AP12MAG(ip12p)/(2*Pi)**3*
     &                                 tmpRho*HC**3.d0
c                              fact=Anucl*(Anucl-1)/2*p12**2*wp12*P12MAG(ip12p)**2*AP12MAG(ip12p)/(2*Pi)**3*
c    &                                 rho(ip12,ip12p,rindx)*HC**3.d0
                           do extQnum=1,extQnumlimit
                              Result(extQnum,twoMzp,twoMz) = Result(extQnum,twoMzp,twoMz) + fact*Int2B(diagNum,extQnum)
                           end do! extQnum
                           end do !diagNum
                        end do    !twoMz
                     end do       !twoMzp
                  end do          !ip12p
               end do             !m12p
            end do                !l12p
         end do                   !s12p
      end do                      !j12p
      
      if (symmetry.eq.1000) continue
      if (verbosity.eq.1000) continue
      return
      end

      subroutine bilinear_interpolate(tmpRho, fQ11,fQ12,fQ21,fQ22,x1,x2,y1,y2,x,y)
c Implimentation of this: https://en.wikipedia.org/wiki/Bilinear_interpolation
c TODO: check if order matters 
      implicit none
      real*8 fQ11, fQ12, fQ21, fQ22
      real*8 x1,x2,y1,y2,x,y
      real*8 tmpRho 
      real*8 term1,term2 

c     real*8 tmpRho2
c     real*8 xVec(2),yVec(2)
c     real*8 pAbs, ppAbs
c     real*8 myMat(2,2), output(2,1)
      term1=((y2-y)/(y2-y1))
      term1=term1*(((x2-x)/(x2-x1))*fQ11+(((x-x1)/(x2-x1))*fQ21))

      term2=((y-y1)/(y2-y1))
      term2=term2*(
     &   ((x2-x)/(x2-x1))*fQ12+
     &   ((x-x1)/(x2-x1))*fQ22)
      tmpRho=term1+term2
      

c     Below this is in theory the same implimentation as above, but its implimented twice as a cross check

c     xVec=(/x2-x,x-x1/)
c     yVec=(/y2-y,y-y1/)
c     myMat=reshape((/fQ11,fQ12,fQ21,fQ22/),shape(myMat))
c     tmp(1)=DOT_PRODUCT(myMat(:,1),yVec)!matrix dot product with a vector
c     tmp(2)=DOT_PRODUCT(myMat(:,2),yVec)

c     tmpRho2=DOT_PRODUCT(xVec, tmp)/((x2-x1)*(y2-y1))

c     if (abs(tmpRho-tmpRho2).ge.1e-8) then
c         write(*,*) ""
c         write(*,*) "abs(tmpRho-tmpRho2).ge.1e-8 evaluated true"
c         write(*,*) "In finalstatesums.twobodyvia2Ndensity.f: tmpRho=",tmpRho 
c         write(*,*) "In finalstatesums.twobodyvia2Ndensity.f: tmpRho2=",tmpRho2
c         write(*,*) "fQ11,fQ12,fQ21,fQ22=",fQ11,fQ12,fQ21,fQ22 
c         write(*,*) "Values dont match, stopping"
c         stop
c     end if

      end subroutine
