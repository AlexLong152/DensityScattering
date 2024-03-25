cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              singlesigmasym : compute (σ1.A) + (1↔2)
c              doublesigmasym : compute (σ1.A)(σ2.B) + (1↔2)
c              singlesigmaasy : compute (σ1.A) - (1↔2)
c              doublesigmaasy : compute (σ1.A)(σ2.B) - (1↔2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, based on spintricks.f, spintricksasy.f of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     hgrie 17 Nov 2023 v2.0 fewbody-Compton: split the following subroutines into new file spinstructures.f and renamed two for more intuitive names:

c         singlesigma     => singlesigmasym
c         Calchold        => doublesigmasym
c         singlesigmaasy  => singlesigmaasy
c         Calcholdasy     => doublesigmaasy
c      
c     This way, spintricks*f only contains individual diagram
c     contributions and not these routines which are generally relevant for spin structures.
c
c     Following comments are relevant to these subroutines and copied from the old spintrics*f files:
c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie 25 Sep 2023: correction in subroutine singlesigma()
c     which computes (σ1+σ2).A and enters in Compton at e²δ⁴ (i.e. affects no publication!)
c     Alex Long found a missing sign in line 2420, which should read
c          hold(1,-1,1,-1)= - hold(1,1,1,1) ! minus sign added!!!
c     Rationale and confirmation by Daniel documented in 
c          documentation/corrections/spintrick-singlesigma.signerror-corrected.20230925.pdf
c     which are scans of manuscript Compton Densities Approach pp. 65a-f.
c     These notes show computation of the MEs of (σ1±σ2).A and claim to show that these are
c     now correctly implemented.
c     These notes therefore also clarify that there is NO mistake in spintrickasy.f's
c     corresponding subroutine singlesigmaasy(),
c     and also checked that the implementation in the deuteron's twobodypolnamp.f,
c     subroutine CalcSdotA(), is correct.
c     Fortunately, we have not yet actually RUN anything of the twobody code at e²δ⁴, so
c     the coding mistake does NOT affect published or produced results.
c     Remember that the e²δ⁴ twobody parts were coded by Arman in analogy to the Beame/... 2003
c     paper for the deuteron,e xtended to inlcude the spin-asymmetric piece (σ1-σ2).A .
c     However, we have as of Sep 2023 not yet CHECKED that/if Arman's code is a correct implementation.
c     So this is the first mistake we find in that part. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c====================================================================
c
c     Order of arguments is
c     hold(Sp,Msp,S,Ms)

      subroutine twosigmas(hold,A,Sp,S,verbosity)
c     Calculates part of spin structure <Sp,mp|σ1.A σ2.A|S,m>
      implicit none
      include '../common-densities/constants.def'
c     
c********************************************************************
c     OUTPUT VARIABLES:
c     
      real*8,intent(out) :: hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
      real*8,intent(in)  :: A(3)
      integer,intent(in) :: verbosity
      integer,intent(in) :: Sp,S
c     
c     Sp,S-final- and initial-state spin
c     
c********************************************************************
c     LOCAL VARIABLES:
c

      complex*16 holdtmp(0:1,-1:1,0:1,-1:1)
      real*8 Bx,By,Bz
      Bx=A(1)
      By=A(2)
      Bz=A(3)

      call doublesigmasym(holdtmp,A(1),A(2),A(3),Bx,By,Bz,Sp,S,verbosity)
      hold=0.5*real(holdtmp)
      end subroutine




      subroutine doublesigmasym(hold,Ax,Ay,Az,Bx,By,Bz,Sp,S,verbosity)
c     
c     Calculates symmetric part of spin structure σ1.A σ2.B: σ1.A σ2.B + σ1.B σ2.A
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     OUTPUT VARIABLES:
c     
      complex*16,intent(out) :: hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8,intent(in)  :: Ax,Ay,Az,Bx,By,Bz
      integer,intent(in) :: verbosity
      integer,intent(in) :: Sp,S
c     
c     Sp,S-final- and initial-state spin
c     
c********************************************************************
c     LOCAL VARIABLES:
c
      real*8 AdotB
      complex*16 Aplus,Aminus,Bplus,Bminus
c     
c     
c********************************************************************
c     
      hold=c0
      AdotB=Ax*Bx+Ay*By+Az*Bz
      
      if ((Sp .eq. 0) .and. (S .eq. 0)) then
         hold(0,0,0,0)=-2.d0*AdotB
      else if ((Sp .eq. 1) .and. (S .eq. 1)) then
         
         Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
         Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))
         Bplus=-(Bx+ci*By)/(dsqrt(2.d0))
         Bminus=(Bx-ci*By)/(dsqrt(2.d0))
         
         hold(1,1,1,1)=(2.d0*Az*Bz)
         hold(1,0,1,1)=-2.d0*(Az*Bplus+Aplus*Bz)
         hold(1,-1,1,1)=4.d0*Aplus*Bplus
         hold(1,1,1,0)=2.d0*(Aminus*Bz+Az*Bminus)
         hold(1,0,1,0)=-2.d0*(Aplus*Bminus+Aminus*Bplus+Az*Bz)
         hold(1,-1,1,0)=-hold(1,0,1,1)
         hold(1,1,1,-1)=4.d0*Aminus*Bminus
         hold(1,0,1,-1)=-hold(1,1,1,0)
         hold(1,-1,1,-1)=hold(1,1,1,1)
      end if
      if (verbosity.eq.1000) continue
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
cccc  new subroutine to calculate the matrix elements of A.σ 
      subroutine singlesigmasym(hold,Ax,Ay,Az,Sp,S,verbosity)
c     Calculates symmetric part of spin structure σ.A: (σ1+σ2).A -- checked in manu-script Compton Densities pp. 65aff
c********************************************************************

      implicit none
      include '../common-densities/constants.def'

c********************************************************************
c     OUTPUT VARIABLE:
      complex*16,intent(out) :: hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8,intent(in)  :: Ax,Ay,Az
      integer,intent(in) :: verbosity
      integer,intent(in) :: Sp,S
c     
c     Sp,S-final- and initial-state spin
c     
c********************************************************************
c     LOCAL VARIABLES:
c     
      complex*16 Aplus,Aminus
c********************************************************************
c     
      hold=c0 ! set all MEs to zero -- only Sp=S=1 MEs are nonzero
      
      if ((Sp .eq. 1) .and. (S .eq. 1)) then
         
         Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
         Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))
         
         hold(1, 1,1, 1) =  2.d0*Az
c     hgrie 25 Sep 2023: following contains added MINUS sign as described at top of file
         hold(1,-1,1,-1) = -hold(1,1,1,1)
         hold(1, 0,1, 1) = -2.d0*Aplus
         hold(1, 1,1, 0) =  2.d0*Aminus   
         hold(1,-1,1, 0) = -2.d0*Aplus !check: is -CONJG(hold(1,1,1,0))
         hold(1, 0,1,-1) =  2.d0*Aminus
      end if
      
      if (verbosity.eq.1000) continue
      end
c     

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine doublesigmaasy(hold,Ax,Ay,Az,Bx,By,Bz,Sp,S,verbosity)   
c     
c     Calculates anti-symmetric part of spin structure σ1.A σ2.B: σ1.A σ2.B - σ1.B σ2.A
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     OUTPUT VARIABLES:
      complex*16,intent(out) :: hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8,intent(in)  :: Ax,Ay,Az,Bx,By,Bz
      integer,intent(in) :: Sp,S
      integer,intent(in) :: verbosity
c
c     Sp,S-final- and initial-state spin
c     
c********************************************************************
c     LOCAL VARIABLES:
      complex*16 Aplus,Aminus,Bplus,Bminus
c     
c********************************************************************
c     
      hold=c0
      Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
      Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))
      Bplus=-(Bx+ci*By)/(dsqrt(2.d0))
      Bminus=(Bx-ci*By)/(dsqrt(2.d0))
      if ((Sp .eq. 0) .and. (S .eq. 1)) then
         hold(0,0,1,1)=2.d0*(-Aplus*Bz+Az*Bplus)
         hold(0,0,1,0)=-2.d0*(Aplus*Bminus-Aminus*Bplus)
         hold(0,0,1,-1)=2.d0*(Aminus*Bz-Az*Bminus)
      else if ((Sp .eq. 1) .and. (S .eq. 0)) then
         hold(1,1,0,0)=2.d0*(Aminus*Bz-Az*Bminus)
         hold(1,0,0,0)=2.d0*(Aplus*Bminus-Aminus*Bplus)
         hold(1,-1,0,0)=2.d0*(-Aplus*Bz+Az*Bplus)
      end if
      
      if (verbosity.eq.1000) continue
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine singlesigmaasy(hold,Ax,Ay,Az,Sp,S,verbosity)
c     
c     Calculates anti-symmetric part of spin structure σ.A: (σ1-σ2).A -- checked in manu-script Compton Densities pp. 65aff
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     OUTPUT VARIABLES:
c     
      complex*16,intent(out) :: hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8,intent(in)  :: Ax,Ay,Az
      integer,intent(in) :: verbosity
      integer,intent(in) :: Sp,S
c     
c********************************************************************
c     LOCAL VARIABLES:
c     
      complex*16 Aplus,Aminus
c     
c********************************************************************
c     
      Aplus=-(Ax+ci*Ay)/(dsqrt(2.d0))
      Aminus=(Ax-ci*Ay)/(dsqrt(2.d0))
      hold=c0
      
      if ((Sp .eq. 0) .and. (S .eq. 1)) then
         hold(0,0,1,1)=-2.d0*Aplus
         hold(0,0,1,0)=2.d0*Az
         hold(0,0,1,-1)=-2.d0*Aminus
      else if ((Sp .eq. 1) .and. (S .eq. 0)) then
         hold(1,1,0,0)=-2.d0*Aminus
         hold(1,0,0,0)=2.d0*Az
         hold(1,-1,0,0)=-2.d0*Aplus
      end if
      if (verbosity.eq.1000) continue
      end
c

      real*8 function ddelta(a,b)
      implicit none
      integer, intent(in) :: a,b
      if (a.eq.b) then
          ddelta=1.d0  
      else
          ddelta=0.d0
      end if
      end function ddelta
