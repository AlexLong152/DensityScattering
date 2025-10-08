c     hgrie Oct 2022: v2.0 fewbody-Compton
c     hgrie Aug 2020: for usesymmetry.and.Mzp=0.andMz=0, Resultyx and Resultyx not calculated
c             They must be zero by symmetry, see manu-script "Compton Densities Approach" p.53
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018: produce onebody amplitudes from 1Ndensities.
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c     Does output of multiple angles & energies to same output file work? May need adjusting output file name.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c
c     hgrie Oct 2022: eliminated mistake in computation of ResultXYXY(): these used to go over ALL twoMzp, twoMzp.
c     That is wrong when one uses usesymmetry() of the amplitudes. In that case, the New Version computes now 2*(2Snucl+1)²
c     independent cartesian ResultXYXY(). However, transcarttosphere() translates them into 2*(2Snucl+1)² +2 spherical amplitudes
c     using that Resultxy=Resultyx=0 for Mzp=Mz=0! 
c     
c     hgrie Aug/Sep 2020: rewrote makedensityfilename() to deal with extracting densities from a .gz or downloading from server
c     hgrie June 2018: renamed "parity" to "symmetry -- see notes in usesymmetry+*.f
c       
c     hgrie May 2018: version 1 based on traditional main.onebody.f
c      
c           All quantum numbers which start with "two" run over integers
c                  Examples:
c                     twoMz,twoMzp: magnetic quantum numbers of in/out target nucleus, times 2.
c                     twoSnucl: 2 x spin of target nucleus
c           For all such variables, run over 2xQM values.
c                  Examples:
c                     in do-loops: twoMz runs from +twoSnucl to -twoSnucl, in steps of -2
c                     in array: Resultxx(twoMzp,twoMz) runs over   (-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
c
c     This numbering agrees with Andreas' assignments of magnetic quantum numbers.
c     In ARRAYs, this means a slight waste of space since we leave many array entries unused.
c     Example S=1/2: 2x2=4 entries needed, but array is (-1:1,-1:1): 3x3=9  : 225% of needed size
c     Example S=1:   3x3=9 entries needed, but array is (-2:2,-2:2): 5x5=25 : 280% of needed size
c
c     Still, our arrays are for small nuclei -- we do not really waste a lot. Not a time/storage issue.
c
c     hgrie May 2018: outsourced symmetry+output into sub routine outputroutine(), identical for onebody and twobody
c      
c     Implemented symmetry for arbitrary nucleon spin:
c     Use Mzp>=0, and for Mzp=0, run only over Mz>=0
c     -- that's still 2 more than necessary since ME(+0->+0) = ME(-0->-0) and ME(+0->-0) = ME(-0->+0)
c     but it's good enough, saving lots of CPU time.
c     see manuscript "Compton Densities Approach" pp11-12
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

      PROGRAM onebodydensitymain

      USE CompDens              ! needs module CompDens.mod

      IMPLICIT NONE
c**********************************************************************
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
      include '../common-densities/constants.def'
c**********************************************************************
c     
c     program argument management
      integer narg              ! number of arguments
      character*500 inputfile   ! argument string is written to this
c     
c*********************************************************************
c     Input information:
c     
c     inUnitno-unit where input information is stored
c     outUnitno-unit where output is written
c     kgamma-photon momentum in cm frame in units of MeV
c     Nangles-number of angles at which calculation is to be done
c     thetaL-photon scattering angle in lab frame, in degrees
c     thetacm-photon scattering angle in cm frame, in radians
c     calctype-which calculation to do
c     
c     frame-which frame the results are to be presented for
c     1=lab. frame
c     2=c.m. frame
c     outfile-name of output file
c     

      integer inUnitno,outUnitno,extQnumlimit

      real*8 Egamma,kgamma,thetaL,thetacm,Elow,Ehigh,Einterval
      
      ! real*8 FT_sPlus, FT_sMinus, FL_sPlus, FL_sMinus
      real*8 E_prot, E_neut, L_prot, L_neut, K1N
      ! real*8 E1N        

      real*8 thetaLow,thetaHigh,thetaInterval
      integer calctype,frame,Nangles,Nenergy,ienergy,j ! number of energies/angles; index for energies/angles
      character*200 descriptors  ! additional descriptors of calculation

      character*3 nucleus ! name of nucleus to be considered, for output file name
      integer Anucl             ! target nucleus mass number
      integer twoSnucl          ! 2 x target nucleus spin
      real*8 Mnucl               ! mass of target nucleus
      
      character*500 outfile
      character*500 densityFileName,originaldensityFileName ! second for multiple energies or angles
      
c*********************************************************************
c     Quadrature variables: 
c     onebody knows only about 1N amplitude's Feynman parameter integration
c     
      integer Nx                ! grid size 
      real*8 xq(Nxmax),wx(Nxmax)! points & weights
      
c     
c----------------------------------------------------------------------
c     
c     Momentum variables:
c     
c     pp-magnitude of final-state relative three-momentum vector,
c     in IA=p + 1/2(k - k') as a vector sum. In IA the
c     kinematics are
c     
c                 \    /
c                  \  /   
c     k' + k/2 + p  \/        -k/2 + p
c     ---------------------------------------
c     
c     
c     
c     ---------------------------------------
c     - k/2 - p
c     
      real*8 k,kth,kphi,kp,kpth,kpphi,Qk,Qkth,Qkphi
      real*8 t,omega
c      
c**********************************************************************
c
c     2* projections of single-nucleon's isospin & spin, both in & out
c     NB: this is the nucleon struck by the photons
      integer twomt1N,twomt1Np,twom1N,twom1Np
      
      integer rindx
      integer i
c     projections of target nucleus' in-spin, out-spin
      integer twoMz,twoMzp

c     hgrie June 2014: added variable twoMzplimit; changed by flag "nosymmetry" in input file.
c     Value twoMzplimit = 0 calculates half of the amplitudes, the other amps then from symmetry
c     Value twoMzplimit = -twoSnucl calculates all amplitudes
      ! integer twoMzplimit
      
c     integer twoMzlimit ! for symmetry calculation: Mzp>=0 *and* for Mzp=0, only Mz>=0, else Mz between +Snucl and -Snucl
      real*8 frac

c     complex*16, allocatable :: SResultx(:,:),SResulty(:,:), SResultz(:,:) ! twoMz from -twoSnucl to
c     complex*16, allocatable :: VResultx(:,:),VResulty(:,:), VResultz(:,:) ! twoMz from -twoSnucl to
      complex*16, allocatable :: FSPlusV(:,:), FSMinusV(:,:)
c     complex*16 :: nucS(3)
c     That means arrays are less than 2^2=4 times bigger than need be, but that's ok since quite small anyway. 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     added by hgrie May 2018: arrays which hold the 1N amplitude in the basis (twomt1N,twom1Np,twom1N,L1N,ML1N)
      
c     NB: these are the quantum numbers of the nucleon struck by the photons
c         (in traditional approach, index "3" was used for spectator) -- see also docu to read1Ndensity()
c     oneNspinbasisXY, where
c     XY: photon helicities in Cartesian basis X: out; Y: in
c     index 1: twomt1N     conserved isospin of nucleon (1: proton, -1: neutron)
c     index 2: twom1Np     1N out-spin                  (1: up, -1: down)
c     index 3: twom1N      1N in-spin                   (1: up, -1: down)
c     index 4: L1N         angular momentum of 1N op.   (integer starts at 0, up to L1Nmax; "K" in Andreas' notes)
c     index 5: ML1N        mag. quantum of L1Nop        (-L1N to L1N: 2L1N+1 integers; "κ" in Andreas' notes)
      
c     At the moment, L1N & ML1N are meaningless (L=ML=0), but they are implemented here as stump already. 

      integer,parameter :: L1Nmax=0    
      integer L1N, ML1N
      complex*16, allocatable :: Result(:,:,:) ! extQnum from 1 to extQnumlimit; twoMzp from -twoSnucl to twoSnucl, stepsize 2; twoMz from -twoSnucl to twoSnucl, stepsize 2; rest blank.
      real*8 aveE0!,aveL0
      complex*16 :: tmpPlus(-1:1, -1:1)
      complex*16 tmpMinus(-1:1, -1:1)
      integer variedA           !BS: integer variedA to indicate which A is varied by calctype=VaryA
      logical cartesian         !hgrie Oct 2014: for output in Cartesian basis of photon polarisations
      
      integer verbosity         !verbosity index for stdout hgrie June 2014
      
      integer test              ! a generic integer for testing i/o
      logical testtf            ! a generic logical for testing i/o
c     if density file generated from a .gz, delete that temporary file after each energy/angle
c     if downloaded and .gz, also delete the download.
c     0: do not delete; 1: delete un-gz'd file; 2: delete downloaded and un-gz'd file 
      integer rmDensityFileLater  

      real*8 :: Iden(-1:1,-1:1)  ! (ms3p,ms3): sigma-0=unit matrix
      complex*16 :: sigmax(-1:1,-1:1)  ! (ms3p,ms3): sigma-x
      complex*16 :: sigmay(-1:1,-1:1) ! (ms3p,ms3): sigma-y
      complex*16 :: sigmaz(-1:1,-1:1)  ! (ms3p,ms3): sigma-z
      complex*16 :: Sigma(-1:1,-1:1)  ! (ms3p,ms3): sigma-z
      complex*16 :: SigmaVec(3,-1:1,-1:1)

      complex*16 :: SpinOnex(-2:2,-2:2)  ! (ms3p,ms3): sigma-x
      complex*16 :: SpinOney(-2:2,-2:2) ! (ms3p,ms3): sigma-y
      complex*16 :: SpinOnez(-2:2,-2:2)  ! (ms3p,ms3): sigma-z
      complex*16 :: SpinOneVec(3,-2:2,-2:2)
      complex*16 :: factorx,factory

      complex*16 :: SpinZerox(-0:0,-0:0)  ! (ms3p,ms3): sigma-x
      complex*16 :: SpinZeroy(-0:0,-0:0) ! (ms3p,ms3): sigma-y
      complex*16 :: SpinZeroz(-0:0,-0:0)  ! (ms3p,ms3): sigma-z
      complex*16 :: SpinZeroVec(3,-0:0,-0:0)

      real*8,parameter ::  munucleon(-1:1) = (/kappan,0.d0,kappap+1.d0/)! indices of entries: (-1,0,+1)!!!!
      real*8 :: kVec(3)  ! Momentum vector of k
      real*8 :: kHat(3)  ! normalized kVec
      integer lab, cm
      character(len=512) :: errmsg
      real*8 eps(3,3)
      integer ieps
      integer ii, jj
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end OF VARIABLE DECLARATIONS, BEGINNING OF CODING
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Reading the input file from command line
c**********************************************************************
c     get the number of arguments
      narg=command_argument_count()

c     if you have 1 argument, write it to inputfile, otherwise stop 
      if (narg.eq.1) then
         call get_command_argument(1, inputfile)
      else
         write(*,*) "*** ERROR: Pass one input file as argument!"
         stop
      end if
      lab=1
      cm=2 
      frame=2 !Assume CM frame
c     
c     
c**********************************************************************
c     Reading in data from the input file
c**********************************************************************
      inUnitno=13
      outUnitno=10
      open(unit=inUnitno, file=inputfile, status= 'OLD',iostat=test)

      if (test .ne. 0) stop "*** ERROR: Could not open input file!!! Aborting."

      call ReadinputCommon(Elow,Ehigh,Einterval,
     &     thetaLow,thetaHigh,thetaInterval,
     &     outfile,descriptors,densityFileName,inUnitno,
     &     nucleus,Anucl,twoSnucl,Mnucl,extQnumlimit,
     &     cartesian,verbosity)
      write(*,*) ""

      if (.not.(extQnumlimit.eq.3)) then
         write(*,*) "Assertion failed: (extQnumlimit.eq.3) evaluated False stopping"
         stop
      end if
      call ReadinputOnebody(inUnitno,calctype,variedA,descriptors,
c---- Variable to control Feynman quadrature settings------------------------
     &     Nx,
     &     verbosity)

      call ReadinputCommonComments(descriptors,inUnitno,verbosity)
      
      close(unit=inUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close input file!!! Aborting."
c
      open(unit=outUnitno, file=outfile, status='unknown', iostat=test,iomsg=errmsg)

      write(*,*) "Pion Photoproduction, onebody, at threshold"
      write(outUnitno,*) "Pion Photoproduction, onebody, at threshold"
      if (test .ne. 0) stop "*** ERROR: In main, could not open output file!!! Aborting."

      write(*,*) "densityFileName=", densityFileName 
      write(outUnitno,*) "densityFileName=", densityFileName 
      call makeoutputfilename(outfile,calctype,nucleus,descriptors,densityFileName,variedA,
     &     Elow,Ehigh,Einterval,thetaLow,thetaHigh,thetaInterval,verbosity)
c**********************************************************************
c     FINAL PRELIMINARIES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     hgrie June 2017: keep original filename: needed for replacements of energy & angle later 
      originaldensityFileName = densityFileName

      thetaL=0.d0
      thetacm=0.d0
c**********************************************************************
c     Setting up quadratures for the Feynman integrals
      call AnglePtsWts(Nx,1,Nxmax,0.d0,1.0d0,xq,wx,Nx,verbosity)
c**********************************************************************
c     BS: new Nangles algebra due to changed input form
      Nangles=int((thetaHigh-thetaLow)/thetaInterval)+1
      Nenergy=int((Ehigh-Elow)/Einterval)+1
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Nov 30 2022 Hard coded values from https://arxiv.org/abs/1103.3400v2 for 3He
c     3He values
c**********************************************************************
      open(unit=outUnitno, file=outfile,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not open output file!!! Aborting."
c**********************************************************************
c     Loop over Energies
c**********************************************************************

      do j=1,Nenergy

         Egamma=Elow+Einterval*(j-1)
         ienergy=int(Egamma)

          kVec = (/0.0,0.0, real(Egamma)/)
          kHat = kVec/SQRT(DOT_PRODUCT(kVec, kVec))
c         write(*,*) "kVec=", kVec
c         write(*,*) "kHat=", kHat
c**********************************************************************
c     Loop over angles
c**********************************************************************
         do i=1,Nangles
            if (frame.eq.lab) then
               kgamma=Egamma/sqrt(1.d0 + 2.0d0*Egamma/Mnucl)
c     BS: new theta algebra due to changed input
c     hgrie Sep 2014: if thetaL is ZERO degrees, actual calculated at 1 Degree  
               thetaL=(thetaLow+real(i-1)*thetaInterval)*Pi/180.d0
               if (thetaL.eq.0.0d0) then
                  thetaL=1.0d0*Pi/180.d0
                  write(*,*) "   Replaced input angle 0 deg with 1 deg."
               end if   

               frac=(Mnucl + (Mnucl + Egamma)*(dcos(thetaL) - 1.0d0))/
     &              (Mnucl + Egamma*(1.d0 - dcos(thetaL)))
               thetacm=dacos(frac)
            else
               thetacm=(thetaLow+real(i-1)*thetaInterval)*Pi/180.d0
               if (thetacm.eq.0.0d0) then
                  thetacm=1.0d0*Pi/180.d0
                  write(*,*) "   Replaced input angle 0 deg with 1 deg."
               end if   
               kgamma=Egamma      
            end if
c**********************************************************************
            call calcphotonmomenta(k,kth,kphi,t,kp,kpth,kpphi,omega,
     &           Qk,Qkth,Qkphi,kgamma,thetacm,verbosity)
            
c**********************************************************************
c      Initialise everything to 0, overwriting entries from previous ω/θ
c**********************************************************************


            allocate(FSPlusV(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(FSMinusV(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))

            FSPlusV=c0
            FSMinusV=c0
            Result=c0


            sigmax=0.d0
            sigmay=0.d0
            sigmaz=0.d0
c     
            sigmax(1,-1)=dcmplx(1.d0,0)
            sigmax(-1,1)=dcmplx(1.d0,0)
c
            sigmay(1,-1)=dcmplx(0, -1.d0)
            sigmay(-1,1)=dcmplx (0, 1.d0) 
c
            sigmaz(1,1)=dcmplx(1.d0,0)
            sigmaz(-1,-1)=dcmplx(-1.d0,0)

            Iden=0.d0
            Iden(-1,-1)=1
            Iden(1,1)=1

            SigmaVec(1,:,:)=sigmax
            SigmaVec(2,:,:)=sigmay
            SigmaVec(3,:,:)=sigmaz
            SigmaVec=SigmaVec/2.d0
            SpinOnex = c0
            SpinOney = c0
            SpinOnez = c0

c--- Raw integer pattern for spin-1 matrices (before normalization)
c--- Sx pattern: [[0,1,0],[1,0,1],[0,1,0]]
            SpinOnex(2,0) = dcmplx(1.d0,0.d0)
            SpinOnex(0,2) = dcmplx(1.d0,0.d0)
            SpinOnex(0,-2) = dcmplx(1.d0,0.d0)
            SpinOnex(-2,0) = dcmplx(1.d0,0.d0)

c--- Sy pattern: [[0,1,0],[-1,0,1],[0,-1,0]]
            SpinOney(2,0) = dcmplx(1.d0,0.d0)
            SpinOney(0,2) = dcmplx(-1.d0,0.d0)
            SpinOney(0,-2) = dcmplx(1.d0,0.d0)
            SpinOney(-2,0) = dcmplx(-1.d0,0.d0)

c--- Sz = diag(1,0,-1)
            SpinOnez(2,2) = dcmplx(1.d0,0.d0)
            SpinOnez(0,0) = dcmplx(0.d0,0.d0)
            SpinOnez(-2,-2) = dcmplx(-1.d0,0.d0)

c--- Apply normalization factors
            factorx = 1.d0 / dsqrt(2.d0)
            factory = 1.d0 / (ci * dsqrt(2.d0))

            SpinOnex = factorx * SpinOnex
            SpinOney = factory * SpinOney

c--- Assemble vector
            SpinOneVec(1,:,:) = SpinOnex
            SpinOneVec(2,:,:) = SpinOney
            SpinOneVec(3,:,:) = SpinOnez

            SpinZerox(:,:)=1.d0
            SpinZeroy(:,:)=1.d0
            SpinZeroz(:,:)=1.d0

            SpinZeroVec(1,:,:) = SpinZerox
            SpinZeroVec(2,:,:) = SpinZeroy
            SpinZeroVec(3,:,:) = SpinZeroz
c**********************************************************************
c     hgrie June 2017: create name of 1Ndensity file for given energy and angle, unpack it
c     define correct formats for energy and angle
c     hgrie May 2018: outsourced into subroutine common-densities/makedensityfilename.f
            densityFileName = originaldensityFileName
            call makedensityfilename(densityFileName,Egamma,thetacm,rmDensityFileLater,verbosity)
c**********************************************************************
c     hgrie May 2018: read 1N density
            call read1Ndensity(densityFileName,Anucl,twoSnucl,omega,thetacm,verbosity)

c           tmpPlus and tmpMinus combines spin and isospin part of diagrams
            tmpPlus = Iden+sigmaz
            tmpMinus = Iden-sigmaz
            eps = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
            do ieps=1,3
            FSMinusV=c0
            FSPlusV=c0
            do rindx=1,maxrho1bindex
                    CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
                    If (L1N.eq.0) then
                        Sigma=SigmaVec(ieps,:,:)
                        FSPlusV(twoMzp,twoMz)=FSPlusV(twoMzp,twoMz)+Anucl*rho1b(rindx)*Sigma(twom1Np,twom1N)
     &                          *tmpPlus(twomt1Np,twomt1N)
                        FSMinusV(twoMzp,twoMz)=FSMinusV(twoMzp,twoMz)+Anucl*rho1b(rindx)*Sigma(twom1Np,twom1N)
     &                          *tmpMinus(twomt1Np,twomt1N)
                    end if ! L1N
            end do              !rindx   

            E_prot = -1.16E-3
            E_neut = 2.13E-3
            L_prot = -1.35E-3
            L_neut = -2.41E-3

            K1N= ((Mnucleon+mpi)/(Mnucl+mpi))*(mnucl/Mnucleon)
            if (twoSnucl.eq.1) then
              FSPlusV=FSPlusV/Sigma
              FSMinusV=FSMinusV/Sigma
            else if (twoSnucl.eq.2) then
              FSPlusV=FSPlusV/SpinOneVec(i,:,:)
              FSMinusV=FSMinusV/SpinOneVec(i,:,:)
            else if (twoSnucl.eq.0) then
              FSPlusV=FSPlusV/SpinZeroVec(i,:,:)!Trivial in this case
              FSMinusV=FSMinusV/SpinZeroVec(i,:,:)
            end if
            Result(ieps,:,:)=(K1N/2.d0)*(E_prot*FSPlusV + E_neut*FSMinusV)*(10**3)

            do ii = -twoSnucl, twoSnucl
            do jj = -twoSnucl, twoSnucl
              if (Result(ieps,ii,jj).ne.Result(ieps,ii,jj) ) then
                Result(ieps,ii,jj)=c0 !If its NaN from dividing by zero, set it to zero
              end if
            end do  
            end do


            write(outUnitno,*) ""
            write(outUnitno,'(A)') "############################################"
            write(outUnitno,'(A,I1,",",I1,",",I1,A)') "eps=", int(eps(ieps,:)), " Result"
            write(*,*) ""
            write(*,'(A)') "############################################"
            write(*,'(A,I1,",",I1,",",I1,A)') "eps=", int(eps(ieps,:)), " Result"

            if (twoSnucl.eq.1) then
              call ResultWrite(FSPlusV,FSMinusV,Sigma,twoSnucl,outUnitno)
            else if (twoSnucl.eq.2) then
              call ResultWrite(FSPlusV,FSMinusV,SpinOneVec,twoSnucl,outUnitno)
            else if  (twoSnucl.eq.0) then
              call ResultWrite(FSPlusV,FSMinusV,SpinZeroVec,twoSnucl,outUnitno)
            end if
            end do!ieps

            if (twoSnucl.eq.1) then
              ! write(*,*) "Result(1,1,-1)=", Result(1,1,-1) 
              ! write(*,*) "Result(2,-1,1)=", Result(2,-1,1) 
              aveE0=(real(Result(1,1,-1))+real(Result(2,-1,1)))/2.d0
            else if (twoSnucl.eq.2) then
              aveE0=(real(Result(1,2,0))+aimag(Result(2,0,2)))/2.d0
            else if  (twoSnucl.eq.0) then
              aveE0=(real(Result(1,0,0))+aimag(Result(2,0,0)))/2.d0
              write(*,*) ""
            end if

            write(*,*) ""
            write(*,*) ""
            write(outUnitno,*) ""
            write(outUnitno,*) ""
            write(*,'(A)') "Result is in units of 10^-3/M_pi^+"
            write(*,'(A,F24.19)') "Average E_0+^1N=",aveE0
            write(outUnitno,'(A)') "Result is in units of 10^-3/M_pi^+"
            write(outUnitno,'(A,F24.19)') "Average E_0+^1N=",aveE0

            call secondOutput(Result,extQnumlimit,twoSnucl,outUnitno)

            write(*,*) ""
            write(outUnitno,*) ""

            write(*,'(A)') "Repeated output for automatted file reading"
            write(outUnitno,'(A)')"Repeated output for automatted file reading" 
            call outputroutine(outUnitno,twoSnucl,extQnumlimit,
     &           Result,verbosity)


            ! In units of m_pi^+
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     be a good boy and deallocate arrays. Compilers do that automatically for simple programs. Better safe than sorry.
            deallocate (FSMinusV,FSPlusV, STAT=test ) ! test becomes nonzero if this fails
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie Aug/Sep 2020: delete the local .dat file if one was generated from .gz
            if (rmDensityFileLater.gt.0) then
               call EXECUTE_COMMAND_LINE("rm "//densityFileName, WAIT=.True., EXITSTAT=test )
               if (test.ne.0) stop "*** ERROR: Could not remove .dat file created from .gz"
               write(*,*) "   Removed .dat file unzipped from .gz."
               if (rmDensityFileLater.ge.2) then
                  INQUIRE(FILE=TRIM(densityFileName)//".gz", EXIST=testtf)
                  if ( testtf ) then
                     call EXECUTE_COMMAND_LINE("rm "//TRIM(densityFileName)//".gz", WAIT=.True., EXITSTAT=test )
                     if (test.ne.0) stop "*** ERROR: Could not remove .dat file created from .gz"
                     write(*,*) "   Removed .dat.gz file downloaded."
                  end if
               end if
            end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         end do                 !Nangles
      end do                    !Nenergy
      close(outUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close output file!!!"
     
      write (*,*) '*** Wrote output to file: ',TRIM(outfile)
      
      ! write(*,*) "Remember to change the value of the energy to match the density "
      stop
      
c20   format(' ',A,I6,A,8I8,A,E24.15,SP,E25.15," I")
c30   format(' ',A,5I4,A,F20.13,SP,F21.13," I")
c40   format(A,2F18.13)
      
      end PROGRAM

      subroutine secondOutput(Result,extQnumlimit,twoSnucl,outUnitno)
      implicit none
      include '../common-densities/constants.def'
      integer twoSnucl, outUnitno,extQnumlimit
      complex*16, intent(in) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      integer i, j,k
      ! character (len=29) fmt

      do k=1,extQnumlimit
      do i=twoSnucl,-twoSnucl,-2
      do j=twoSnucl, -twoSnucl,-2
        write(*,*) Result(k,i,j)
        write(outUnitno,*) Result(k,i,j)
      end do
      end do
      end do 
      end subroutine



      subroutine ResultWrite(FPlus,FMinus,Sigma,twoSnucl,outUnitno)
      implicit none
      include '../common-densities/constants.def'
      integer twoSnucl, outUnitno
      complex*16, intent(in) :: FPlus(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16, intent(in) :: FMinus(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16, intent(in) :: Sigma(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)  ! (ms3p,ms3): sigma-z
      integer i, j
      character (len=29) fmt


      ! write(*,*) "twoSnucl=", twoSnucl 
      if (twoSnucl.ne.0) then!Actually want to see the zeros in this case
      fmt='(A,I0,A,I0,A,F8.4,SP,F8.4,A)'
        do i=-twoSnucl,twoSnucl
        do j=-twoSnucl,twoSnucl
          if ((Sigma(i,j).ne.c0).and.(FMinus(i,j).ne.c0)) then
          if (FMinus(i,j).eq.FMinus(i,j)) then
                  WRITE(*,fmt) 'FMinus(', i, ',', j, ') = ', FMinus(i,j),'i'!/Sigma(i,j),'j'
          end if
          end if        
        end do
        end do

        do i=-twoSnucl,twoSnucl
        do j=-twoSnucl,twoSnucl
          if ((Sigma(i,j).ne.c0).and.(FPlus(i,j).ne.c0)) then

          if (FPlus(i,j).eq.FPlus(i,j)) then
                  WRITE(*,fmt) 'FPlus(', i, ',', j, ') = ', FPlus(i,j),'i'!/Sigma(i,j),'j'
          end if
          end if        
        end do
        end do
        
        do i=-twoSnucl,twoSnucl
        do j=-twoSnucl,twoSnucl
          if ((Sigma(i,j).ne.c0).and.(FMinus(i,j).ne.c0)) then
          if (FMinus(i,j).eq.FMinus(i,j)) then
                  WRITE(outUnitno,fmt) 'FMinus(', i, ',', j, ') = ', FMinus(i,j),'i'!/Sigma(i,j),'j'
          end if 
          end if        
        end do
        end do

        do i=-twoSnucl,twoSnucl
        do j=-twoSnucl,twoSnucl
          if ((Sigma(i,j).ne.c0).and.(FPlus(i,j).ne.c0)) then
          if (FPlus(i,j).eq.FPlus(i,j)) then
                  WRITE(outUnitno,fmt) 'FPlus(', i, ',', j, ') = ', FPlus(i,j),'i'!/Sigma(i,j),'j'
          end if
          end if        
        end do
        end do
      else
      ! fmt='(A,I0,A,I0,A,F8.12,SP,F8.12,A)'
      fmt='(A,I0,A,I0,A,F8.4,SP,F8.4,A)'
        do i=-twoSnucl,twoSnucl
        do j=-twoSnucl,twoSnucl
                  WRITE(*,fmt) 'FMinus(', i, ',', j, ') = ', FMinus(i,j),'i'!/Sigma(i,j),'j'
        end do
        end do

        do i=-twoSnucl,twoSnucl
        do j=-twoSnucl,twoSnucl
                  WRITE(*,fmt) 'FPlus(', i, ',', j, ') = ', FPlus(i,j),'i'!/Sigma(i,j),'j'
        end do
        end do

        do i=-twoSnucl,twoSnucl
        do j=-twoSnucl,twoSnucl
                  WRITE(outUnitno,fmt) 'FMinus(', i, ',', j, ') = ', FMinus(i,j),'i'!/Sigma(i,j),'j'
        end do
        end do

        do i=-twoSnucl,twoSnucl
        do j=-twoSnucl,twoSnucl
                  WRITE(outUnitno,fmt) 'FPlus(', i, ',', j, ') = ', FPlus(i,j),'i'!/Sigma(i,j),'j'
        end do
        end do
      end if
      end subroutine
