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
      integer inUnitno,outUnitno,lab,cm

      real*8 Egamma,kgamma,thetaL,thetacm,Elow,Ehigh,Einterval
      
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
      
      integer i
      
c     projections of target nucleus' in-spin, out-spin
      integer twoMz,twoMzp

c     hgrie June 2014: added variable twoMzplimit; changed by flag "nosymmetry" in input file.
c     Value twoMzplimit = 0 calculates half of the amplitudes, the other amps then from symmetry
c     Value twoMzplimit = -twoSnucl calculates all amplitudes
      integer twoMzplimit
      
      integer twoMzlimit ! for symmetry calculation: Mzp>=0 *and* for Mzp=0, only Mz>=0, else Mz between +Snucl and -Snucl

      real*8 frac

      complex*16, allocatable :: Resultxx(:,:),Resultxy(:,:) ! twoMz from -twoSnucl to twoSnucl, stepsize 2; rest blank.
      complex*16, allocatable :: Resultyx(:,:),Resultyy(:,:) ! twoMz from -twoSnucl to twoSnucl, stepsize 2; rest blank.
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
      integer rindx

      complex*16,allocatable :: oneNspinbasisxx(:,:,:,:,:) ! (twomt1N,twom1Np,twom1N,L1N,ML1N)
      complex*16,allocatable :: oneNspinbasisxy(:,:,:,:,:) ! (twomt1N,twom1Np,twom1N,L1N,ML1N)
      complex*16,allocatable :: oneNspinbasisyx(:,:,:,:,:) ! (twomt1N,twom1Np,twom1N,L1N,ML1N)
      complex*16,allocatable :: oneNspinbasisyy(:,:,:,:,:) ! (twomt1N,twom1Np,twom1N,L1N,ML1N)
      
c     1N amps of proton or neutron outside fewbody loops: array with 1: proton (twomt1N=+1); -1: neutron (twomt1N=-1)
      real*8 A1(-1:1),A2(-1:1),A3(-1:1),A4(-1:1),A5(-1:1),A6(-1:1)
c     real*8 A1p,A2p,A3p,A4p,A5p,A6p
c     real*8 A1n,A2n,A3n,A4n,A5n,A6n
      
      integer variedA           !BS: integer variedA to indicate which A is varied by calctype=VaryA
      logical cartesian         !hgrie Oct 2014: for output in Cartesian basis of photon polarisations
      
      integer verbosity         !verbosity index for stdout hgrie June 2014
      
      integer test              ! a generic integer for testing i/o
      logical testtf            ! a generic logical for testing i/o
c     if density file generated from a .gz, delete that temporary file after each energy/angle
c     if downloaded and .gz, also delete the download.
c     0: do not delete; 1: delete un-gz'd file; 2: delete downloaded and un-gz'd file 
      integer rmDensityFileLater  

      real*8 dummy

c     for calculating magnetic-moment insertions on the way: (twoMzp,twoMzp,twomt1N)
      real*8,allocatable     :: insertion0(:,:,:),insertionz(:,:,:),insertiony(:,:,:),insertionx(:,:,:)
      real*8 :: sigma0(-1:1,-1:1)  ! (ms3p,ms3): sigma-0=unit matrix
      real*8 :: sigmax(-1:1,-1:1)  ! (ms3p,ms3): sigma-x
      real*8 :: isigmay(-1:1,-1:1) ! (ms3p,ms3): I times sigma-y !!!!
      real*8 :: sigmaz(-1:1,-1:1)  ! (ms3p,ms3): sigma-z
      real*8,parameter ::  munucleon(-1:1) = (/kappan,0.d0,kappap+1.d0/)! indices of entries: (-1,0,+1)!!!!
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end OF VARIABLE DECLARATIONS, BEGINNING OF CODING
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) 
      write(*,*) "================================================================================"
      write(*,*) "Onebody Contributions to Few-Nucleon Compton Scattering Calculated Via 1N-Density Matrix"
      write(*,*) "================================================================================"
      write(*,*) "   Version 1.0"
      write(*,*) "      D. Phillips/A. Nogga/hgrie starting August 2020   "
      write(*,*) "      based on 3He codes: D. Phillips/A. Nogga/hgrie starting May 2018"
      write(*,*) "                          D. Phillips/B. Strasberg/hgrie starting June 2014"
      write(*,*) "                          with A. Margaryan 2016-17, modifying codes by D. Shukla 2007/8"
      write(*,*)
c**********************************************************************
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
c     
c     
c**********************************************************************
c     Reading in data from the input file
c**********************************************************************
      inUnitno=13
      outUnitno=10
      open(unit=inUnitno, file=inputfile, status= 'OLD',iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not open input file!!! Aborting."

      call ReadinputCommon(Elow,Ehigh,Einterval,frame,
     &     thetaLow,thetaHigh,thetaInterval,
     &     outfile,descriptors,densityFileName,inUnitno,
     &     nucleus,Anucl,twoSnucl,Mnucl,
     &     twoMzplimit,cartesian,verbosity)
c      
      call ReadinputOnebody(inUnitno,calctype,variedA,descriptors,
c---- Variable to control Feynman quadrature settings------------------------
     &     Nx,
     &     verbosity)

      call ReadinputCommonComments(descriptors,inUnitno,verbosity)
      
      close(unit=inUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close input file!!! Aborting."
c
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
      
      write (*,*) "***************************** END OF INITIALISATION *****************************"
c**********************************************************************
      open(unit=outUnitno, file=outfile,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not open output file!!! Aborting."
c**********************************************************************
c     Loop over Energies
c**********************************************************************
      do j=1,Nenergy
         Egamma=Elow+Einterval*(j-1)
         ienergy=int(Egamma)
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
         write(*,*)
         write(*,*) "Incoming photon energy (rounded) = ", ienergy, " MeV"
         write(*,*) "Number of Angles Nangles =   ", Nangles
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
            if (frame.eq.cm) then
               write (outUnitno,*) "cm ","omega =",Egamma,"thetacm =",thetacm*180.0/Pi
            else if (frame.eq.lab) then
               write (outUnitno,*) "lab ","omega =",Egamma,"thetalab =",thetaL*180.0/Pi
            end if
            write(*,*)
            write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            write(*,*) 'Calculating amps: theta =',thetacm*180.0/Pi,'deg: angle #',i
c**********************************************************************
            call calcphotonmomenta(k,kth,kphi,t,kp,kpth,kpphi,omega,
     &           Qk,Qkth,Qkphi,kgamma,thetacm,verbosity)
            
c**********************************************************************
c      be a good boy and initialise everything to 0, overwriting entries from previous ω/θ
c**********************************************************************
            write(*,*) "   Allocating 1N operators: At present, only L1Nmax=0 implemented  (K=0 in Andreas' notes)."
            allocate (oneNspinbasisxx(-1:1,-1:1,-1:1,0:L1Nmax,-L1Nmax:L1Nmax))
            allocate (oneNspinbasisxy(-1:1,-1:1,-1:1,0:L1Nmax,-L1Nmax:L1Nmax))
            allocate (oneNspinbasisyx(-1:1,-1:1,-1:1,0:L1Nmax,-L1Nmax:L1Nmax))
            allocate (oneNspinbasisyy(-1:1,-1:1,-1:1,0:L1Nmax,-L1Nmax:L1Nmax))
            oneNspinbasisxx = c0
            oneNspinbasisxy = c0
            oneNspinbasisyx = c0
            oneNspinbasisyy = c0

            allocate(Resultxx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(Resultxy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(Resultyx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            allocate(Resultyy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            Resultxx=c0
            Resultxy=c0
            Resultyx=c0
            Resultyy=c0

c     for calculating electric FF on the way
            dummy=0.d0
c     for calculating magnetic moment insertions on the way
            allocate(insertion0(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
            allocate(insertionx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
            allocate(insertiony(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
            allocate(insertionz(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl,-1:1))
            insertion0 =0.0
            insertionx =0.0
            insertiony =0.0
            insertionz =0.0
            sigma0=0.d0
            sigmax=0.d0
            isigmay=0.d0
            sigmaz=0.d0
c     
            sigma0(1,1)=1.d0
            sigma0(-1,-1)=1.d0
c     
            sigmax(1,-1)=1.d0
            sigmax(-1,1)=1.d0
c            
            isigmay(1,-1)=1.d0  ! this is I σy, NOT σy !
            isigmay(-1,1)=-1.d0
c
            sigmaz(1,1)=1.d0
            sigmaz(-1,-1)=-1.d0     
c**********************************************************************
c     hgrie June 2017: create name of 1Ndensity file for given energy and angle, unpack it
c     define correct formats for energy and angle
c     hgrie May 2018: outsourced into subroutine common-densities/makedensityfilename.f
            densityFileName = originaldensityFileName
            call makedensityfilename(densityFileName,Egamma,thetacm,rmDensityFileLater,verbosity)
c**********************************************************************
c     hgrie May 2018: read 1N density
            call read1Ndensity(densityFileName,Anucl,twoSnucl,omega,thetacm,verbosity)
c            write(*,*) "J12MAX",j12max_rho
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
c**********************************************************************
c     hgrie June 2014:
c     now calculate A1 to A6 for the particular gamma-N energy and angle
c     *outside* the few-body loops: speeds up computation.
c     at present, the "print1Namps" option is hardwired to .True.:
c     print 1N amplitudes used to screen.
            
c     BS: Implemented A variations to ConstructAmps. calctype and
c     variedA determine which amplitude to vary
            write(*,*) "  Now creating one-nucleon helicity amplitudes."
            do L1N=0,L1Nmax
               do ML1N=-L1N,L1N
c     at present, the "print1Namps" option is hardwired to .True.:
c     print 1N amplitudes used to screen.
                  call construct1NAmps(A1(1),A2(1),A3(1),A4(1),A5(1),A6(1),
     &                 A1(-1),A2(-1),A3(-1),A4(-1),A5(-1),A6(-1),
     &                 thetacm,xq,wx,Nx,t,omega,calctype,.True.,
     &                 variedA,verbosity)
                  
c                  write(*,*) "   Proton:"
c                  write(*,*) "       A1p = ",A1(1)
c                  write(*,*) "       A2p = ",A2(1)
c                  write(*,*) "       A3p = ",A3(1)
c                  write(*,*) "       A4p = ",A4(1)
c                  write(*,*) "       A5p = ",A5(1)
c                  write(*,*) "       A6p = ",A6(1)
c                  write(*,*) "   Neutron:"
c                  write(*,*) "       A1n = ",A1(-1)
c                  write(*,*) "       A2n = ",A2(-1)
c                  write(*,*) "       A3n = ",A3(-1)
c                  write(*,*) "       A4n = ",A4(-1)
c                  write(*,*) "       A5n = ",A5(-1)
c                  write(*,*) "       A6n = ",A6(-1)
c                  
c     transform 1N amplitudes into spin/helicity basis -- one routine per photon-helicity pair
c              in trafo1NampAB: A: helicity in-photon; B: helicity out-photon
                  do twomt1N=1,-1,-2
                     call trafo1Nampxx(oneNspinbasisxx(twomt1N,:,:,L1N,ML1N),
     &                    A1(twomt1N),A2(twomt1N),A3(twomt1N),A4(twomt1N),A5(twomt1N),A6(twomt1N),
     &                    kth,kphi,kpth,kpphi,verbosity)
                     call trafo1Nampyx(oneNspinbasisyx(twomt1N,:,:,L1N,ML1N),
     &                    A1(twomt1N),A2(twomt1N),A3(twomt1N),A4(twomt1N),A5(twomt1N),A6(twomt1N),
     &                    kth,kphi,kpth,kpphi,verbosity)
                     call trafo1Nampxy(oneNspinbasisxy(twomt1N,:,:,L1N,ML1N),
     &                    A1(twomt1N),A2(twomt1N),A3(twomt1N),A4(twomt1N),A5(twomt1N),A6(twomt1N),
     &                    kth,kphi,kpth,kpphi,verbosity)
                     call trafo1Nampyy(oneNspinbasisyy(twomt1N,:,:,L1N,ML1N),
     &                    A1(twomt1N),A2(twomt1N),A3(twomt1N),A4(twomt1N),A5(twomt1N),A6(twomt1N),
     &                    kth,kphi,kpth,kpphi,verbosity)
                  end do        !twomt1N
               end do           !ML1N
            end do              !L1N 
cccccccccccccccccccccccc                 
            if (verbosity.ge.2) then
               write(*,*) "1N amplitudes in helicity basis:"
               do L1N=0,L1Nmax
                  do ML1N=-L1N,L1N
                     do twomt1N=1,-1,-2
                        write(*,*) "   (twomt1N,twom1Np,twom1N,L1N,ML1N); Photon helicities: XX"
                        do twom1Np=1,-1,-2
                           do twom1N=1,-1,-2
                              write(*,30) "  ",twomt1N,twom1Np,twom1N,L1N,ML1N," : ",
     &                             oneNspinbasisxx(twomt1N,twom1Np,twom1N,L1N,ML1N)
                           end do !twom1N
                        end do  !twom1Np
                        write(*,*) "   (twomt1N,twom1Np,twom1N,L1N,ML1N); Photon helicities: YX"
                        do twom1Np=1,-1,-2
                           do twom1N=1,-1,-2
                              write(*,30) "  ",twomt1N,twom1Np,twom1N,L1N,ML1N," : ",
     &                             oneNspinbasisyx(twomt1N,twom1Np,twom1N,L1N,ML1N)
                           end do !twom1N
                        end do  !twom1Np
                        write(*,*) "   (twomt1N,twom1Np,twom1N,L1N,ML1N); Photon helicities: XY"
                        do twom1Np=1,-1,-2
                           do twom1N=1,-1,-2
                              write(*,30) "  ",twomt1N,twom1Np,twom1N,L1N,ML1N," : ",
     &                             oneNspinbasisxy(twomt1N,twom1Np,twom1N,L1N,ML1N)
                           end do !twom1N
                        end do  !twom1Np
                        write(*,*) "   (twomt1N,twom1Np,twom1N,L1N,ML1N); Photon helicities: YY"
                        do twom1Np=1,-1,-2
                           do twom1N=1,-1,-2
                              write(*,30) "  ",twomt1N,twom1Np,twom1N,L1N,ML1N," : ",
     &                             oneNspinbasisyy(twomt1N,twom1Np,twom1N,L1N,ML1N)
                           end do !twom1N
                        end do  !twom1Np
                     end do     !twomt1N
                  end do        !ML1N
               end do           !L1N
            end if              !verbosity      
c**********************************************************************
c**********************************************************************
            write(*,*) "*********Now convoluting 1N helicity amplitudes with 1N density matrix.*********"
c     following is the explicit summation directly over all quantum numbers -- only in verbosity mode
            if (verbosity.ge.5) then
               write(*,*) "   rindx:    twoMzp, twoMz, twomt1Np, twomt1N, twom1Np, twom1N, L1N, ML1N:   ρ1(rindx)"
               do twoMzp=twoSnucl,twoMzplimit,-2
c         for Mzp=0, run only over Mz>=0 -- that's still 2 more than necessary, but good enough -- cured below
                  if ((twoMzp.eq.0).and.(twoMzplimit.eq.0)) then
                     twoMzlimit = 0
                  else
                     twoMzlimit = -twoSnucl
                  end if   
                  do twoMz=twoSnucl,twoMzlimit,-2
                     do L1N=0,L1Nmax
                        do ML1N=-L1N,L1N
                           do twomt1N=1,-1,-2
                              do twom1Np=1,-1,-2
                                 do twom1N=1,-1,-2
                                    twomt1Np=twomt1N
c     following is the empirical formula found from output of readdensity subroutine by Andreas
c     with that, the reuslts of summing over ridndx or over quantum numbers are identical
                                    rindx = 1+ (1+twom1N)/2 + 2*(1+twom1Np)/2 + 4*(1+twomt1N)/2 +
     &                                   8*(1+twoMz)/2 + 16*(1+twoMzp)/2 + 32*L1N + 64*(ML1N+L1N)
                                 
                                    if (verbosity.ge.5) then
                                       write(*,20) "  ",rindx,":    ",
     &                                      twoMzp,twoMz,twomt1Np,twomt1N,twom1Np,twom1N,L1N,ML1N,": ",rho1b(rindx)
                                    end if   
                                    Resultxx(twoMzp,twoMz) = Resultxx(twoMzp,twoMz) +
     &                                   oneNspinbasisxx(twomt1N,twom1Np,twom1N,L1N,ML1N)*rho1b(rindx)*Anucl
                                    Resultyy(twoMzp,twoMz) = Resultyy(twoMzp,twoMz) +
     &                                   oneNspinbasisyy(twomt1N,twom1Np,twom1N,L1N,ML1N)*rho1b(rindx)*Anucl
c hgrie Aug 2020: now cure: for Mzp=Mz=0, only calculate xx and yy, since xy and yx must be zero                                    
                                    if ((twoMzplimit.eq.0).and.(twoMzp.eq.0).and.(twoMz.eq.0)) then
                                       continue
                                    else
                                       Resultyx(twoMzp,twoMz) = Resultyx(twoMzp,twoMz) +
     &                                      oneNspinbasisyx(twomt1N,twom1Np,twom1N,L1N,ML1N)*rho1b(rindx)*Anucl
                                       Resultxy(twoMzp,twoMz) = Resultxy(twoMzp,twoMz) +
     &                                      oneNspinbasisxy(twomt1N,twom1Np,twom1N,L1N,ML1N)*rho1b(rindx)*Anucl
                                    end if
c end cure                  
                                 end do !twom1N
                              end do !twom1Np
                           end do !twomt1N
                        end do  !ML1N
                     end do     !L1N
                  end do        !twoMz
               end do           !twoMzp
ccccccccccccccccccccccccccccccc following to compare results of explicit (above) and implicit (below) summations
               do twoMzp=twoSnucl,twoMzplimit,-2
c         for Mzp=0, run only over Mz>=0 -- that's still 2 more than necessary, but good enough -- cured below
                  if ((twoMzp.eq.0).and.(twoMzplimit.eq.0)) then
                     twoMzlimit = 0
                  else
                     twoMzlimit = -twoSnucl
                  end if   
                  do twoMz=twoSnucl,twoMzlimit,-2
                     write (*,*) "Resultxx(twoMzp,twoMz): ",twoMzp,twoMz,Resultxx(twoMzp,twoMz)
                     if ((twoMzplimit.eq.0).and.(twoMzp.eq.0).and.(twoMz.eq.0)) then
                        continue
                     else
                        write (*,*) "Resultxy(twoMzp,twoMz): ",twoMzp,twoMz,Resultxy(twoMzp,twoMz)
                        write (*,*) "Resultyx(twoMzp,twoMz): ",twoMzp,twoMz,Resultyx(twoMzp,twoMz)
                     end if
                     write (*,*) "Resultyy(twoMzp,twoMz): ",twoMzp,twoMz,Resultyy(twoMzp,twoMz)
                  end do
               end do
               Resultxx=c0      ! must reset Result
               Resultxy=c0
               Resultyx=c0
               Resultyy=c0
            end if              ! verbosity   
ccccccccccccccccccccccccccccccc
c     Now the actual sumation over quantum numbers
            write(*,*) "Beep boop in main in onebodyvia oneN"
            stop
            if (verbosity.ge.3) then ! provide a header for table
               write(*,*) "   rindx:    twoMzp, twoMz, twomt1Np, twomt1N, twom1Np, twom1N, L1N, ML1N:   ρ1(rindx)"
            end if              ! verbosity
c            
            do rindx=1,maxrho1bindex
               CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
c     only part of sum when quantum numbers match: 1N isospin unchanged, matching L1N
c     and must translate FROM Andreas' definition of quantum numbers TO ours
               if ((twomt1N.eq.twomt1Np).and.(L1N.le.L1Nmax)) then
c     hgrie Oct 2022: run only over those MEs needed when usesymetry() is set; otherwise, run over all
                  if ((twoMzplimit.eq.0).and.(twoMzp.lt.0)) then
                     continue
                  else
                     if ((twoMzplimit.eq.0).and.(twoMzp.eq.0).and.(twoMz.lt.0)) then
                        continue
                     else
                        if (verbosity.ge.3) then
                           write(*,20) "  ",rindx,":",twoMzp,twoMz,twomt1Np,twomt1N,twom1Np,twom1N,L1N,ML1N,": ",rho1b(rindx)
c                          write(*,*) "    oneNspinbasisxx:",oneNspinbasisxx(twomt1N,twom1Np,twom1N,L1N,ML1N)
                        end if  ! verbosity
c     following gives results of insertions of σ[0-3]:
                        if (verbosity.ge.3) then
                           if (twom1N.eq.twom1Np) then
                              dummy = dummy + rho1b(rindx)/2.d0 ! calculate electric FF
                           end if
                           insertion0(twoMzp,twoMz,twomt1N)=
     &                          insertion0(twoMzp,twoMz,twomt1N)+rho1b(rindx)*sigma0(twom1Np,twom1N)*Anucl ! sigma0 insertion
                           insertionx(twoMzp,twoMz,twomt1N)=
     &                          insertionx(twoMzp,twoMz,twomt1N)+rho1b(rindx)*sigmax(twom1Np,twom1N)*Anucl ! sigmax insertion
                           insertiony(twoMzp,twoMz,twomt1N)=
     &                          insertiony(twoMzp,twoMz,twomt1N)+rho1b(rindx)*isigmay(twom1Np,twom1N)*Anucl ! I*sigmay insertion
                           insertionz(twoMzp,twoMz,twomt1N)=
     &                          insertionz(twoMzp,twoMz,twomt1N)+rho1b(rindx)*sigmaz(twom1Np,twom1N)*Anucl ! sigmaz insertion
                        end if  ! verbosity/σ[0-3] insertions
c
c     hgrie May 2018: factor Anucl because each of the nucleons can be struck
                        Resultxx(twoMzp,twoMz) = Resultxx(twoMzp,twoMz) +
     &                       oneNspinbasisxx(twomt1N,twom1Np,twom1N,L1N,ML1N)*rho1b(rindx)*Anucl
                        Resultyy(twoMzp,twoMz) = Resultyy(twoMzp,twoMz) +
     &                       oneNspinbasisyy(twomt1N,twom1Np,twom1N,L1N,ML1N)*rho1b(rindx)*Anucl
c     hgrie Aug 2020: now cure: if usesymmetry(): for Mzp=Mz=0, only calculate xx and yy, since xy and yx must be zero                               
                        if ((twoMzplimit.eq.0).and.(twoMzp.eq.0).and.(twoMz.eq.0)) then
                           continue
                        else
                           Resultyx(twoMzp,twoMz) = Resultyx(twoMzp,twoMz) +
     &                          oneNspinbasisyx(twomt1N,twom1Np,twom1N,L1N,ML1N)*rho1b(rindx)*Anucl
                           Resultxy(twoMzp,twoMz) = Resultxy(twoMzp,twoMz) +
     &                          oneNspinbasisxy(twomt1N,twom1Np,twom1N,L1N,ML1N)*rho1b(rindx)*Anucl
                        end if ! ((twoMzplimit.eq.0).and.(twoMzp.eq.0).and.(twoMz.eq.0))
c     end cure
                     end if ! ((twoMzplimit.eq.0).and.(twoMzp.eq.0).and.(twoMz.lt.0))
                  end if ! ((twoMzplimit.eq.0).and.(twoMzp.lt.0))
               end if ! ((twomt1N.eq.twomt1Np).and.(L1N.le.L1Nmax))
            end do              !rindx   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            if (verbosity.ge.3) then
               write(*,*) "  qval [fm⁻¹] = ", qval
               write(*,*) "  Electric Form Factor from 1Ndensity: ",dummy
               write(*,*) "                                  (twoMz=+1)   (twoMz=-1)  "
               write(*,*) "   ME: insert  σ0 into proton : "
               write(*,40) "                  (twoMzp=+1) ",insertion0(1,1,1),insertion0(1,-1,1)
               write(*,40) "                  (twoMzp=-1) ",insertion0(-1,1,1),insertion0(-1,-1,1)
               write(*,*) "   ME: insert  σx into proton : "
               write(*,40) "                  (twoMzp=+1) ",insertionx(1,1,1),insertionx(1,-1,1)
               write(*,40) "                  (twoMzp=-1) ",insertionx(-1,1,1),insertionx(-1,-1,1)
               write(*,*) "   ME: insert iσy into proton : "
               write(*,40) "                  (twoMzp=+1) ",insertiony(1,1,1),insertiony(1,-1,1)
               write(*,40) "                  (twoMzp=-1) ",insertiony(-1,1,1),insertiony(-1,-1,1)
               write(*,*) "   ME: insert  σz into proton : "
               write(*,40) "                  (twoMzp=+1) ",insertionz(1,1,1),insertionz(1,-1,1)
               write(*,40) "                  (twoMzp=-1) ",insertionz(-1,1,1),insertionz(-1,-1,1)
               write(*,*)
               write(*,*) "   ME: insert  σ0 into neutron : "
               write(*,40) "                  (twoMzp=+1) ",insertion0(1,1,-1),insertion0(1,-1,-1)
               write(*,40) "                  (twoMzp=-1) ",insertion0(-1,1,-1),insertion0(-1,-1,-1)
               write(*,*) "   ME: insert  σx into neutron : "
               write(*,40) "                  (twoMzp=+1) ",insertionx(1,1,-1),insertionx(1,-1,-1)
               write(*,40) "                  (twoMzp=-1) ",insertionx(-1,1,-1),insertionx(-1,-1,-1)
               write(*,*) "   ME: insert iσy into neutron : "
               write(*,40) "                  (twoMzp=+1) ",insertiony(1,1,-1),insertiony(1,-1,-1)
               write(*,40) "                  (twoMzp=-1) ",insertiony(-1,1,-1),insertiony(-1,-1,-1)
               write(*,*) "   ME: insert  σz into neutron : "
               write(*,40) "                  (twoMzp=+1) ",insertionz(1,1,-1),insertionz(1,-1,-1)
               write(*,40) "                  (twoMzp=-1) ",insertionz(-1,1,-1),insertionz(-1,-1,-1)
c
               write(*,*)  "      Derived from that: "
               write(*,*)  "      Mag. FF via  σx: ",insertionx(::2,::2,1)*munucleon(1)+insertionx(::2,::2,-1)*munucleon(-1)
               write(*,*)  "      Mag. FF via iσy: ",insertiony(::2,::2,1)*munucleon(1)+insertiony(::2,::2,-1)*munucleon(-1)
               write(*,*)  "      Mag. FF via  σz: ",insertionz(::2,::2,1)*munucleon(1)+insertionz(::2,::2,-1)*munucleon(-1)
            end if
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c CALCULATE FF
c            if (verbosity.ge.1) then
c               dummy=0.d0
c               do rindx=1,maxrho1bindex
c                  CALL get1Nqnnum(rindx,twom1N,twomt1N,twoMz,twom1Np,twomt1Np,twoMzp,L1N,ML1N)
cc     calculate FF
c                  if ((twomt1N.eq.twomt1Np).and.(L1N.eq.0).and.(ML1N.eq.0).and.(twom1N.eq.twom1Np)) then
c                     if (verbosity.ge.1) then
c                        write(*,20) "  ",rindx,":    ",twoMzp,twoMz,twomt1Np,twomt1N,twom1Np,twom1N,L1N,ML1N,": ",rho1b(rindx)
c                     end if
c                     dummy = dummy + rho1b(rindx)/2.d0
c                  end if
c               end do           !rindx
c               write(*,*) "FF from density file: ",dummy
c            end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
            if (verbosity.ge.3) then
               do twoMzp=twoSnucl,twoMzplimit,-2
c         for Mzp=0, run only over Mz>=0 -- that's still 2 more than necessary, but good enough -- see cure below
                  if ((twoMzp.eq.0).and.(twoMzplimit.eq.0)) then
                     twoMzlimit = 0
                  else
                     twoMzlimit = -twoSnucl
                  end if   
                  do twoMz=twoSnucl,twoMzlimit,-2
                     write (*,*) "Resultxx(twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Resultxx(twoMzp,twoMz)
c hgrie Aug 2020: now cure: for Mzp=Mz=0, only calculate xx and yy, since xy and yx must be zero               
                     if ((twoMzplimit.eq.0).and.(twoMzp.eq.0).and.(twoMz.eq.0)) then
                        continue
                     else
                        write (*,*) "Resultxy(twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Resultxy(twoMzp,twoMz)
                        write (*,*) "Resultyx(twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Resultyx(twoMzp,twoMz)
                     end if
                     write (*,*) "Resultyy(twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Resultyy(twoMzp,twoMz)
                  end do
               end do
            end if
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: symmetry and output
            call outputroutine(outUnitno,cartesian,twoSnucl,twoMzplimit,
     &           Resultxx,Resultxy,Resultyx,Resultyy,verbosity)
            
c     be a good boy and deallocate arrays. Compilers do that automatically for simple programs. Better safe than sorry.
            deallocate (Resultxx,Resultxy,Resultyx,Resultyy, STAT=test ) ! test becomes nonzero if this fails
            if (test .ne. 0) stop "*** ERROR: Arrays ResulyAB: Deallocation error. Abort."
            deallocate (oneNspinbasisxx,oneNspinbasisxy,oneNspinbasisyx,oneNspinbasisyy, STAT=test ) ! test nonzero if fails
            if (test .ne. 0) stop "*** ERROR: Arrays oneNspinbasisAB: Deallocation error. Abort."
            deallocate (insertion0,insertionx,insertiony,insertionz, STAT=test ) ! test becomes nonzero if this fails
            if (test .ne. 0) stop "*** ERROR: Arrays insertion: Deallocation error. Abort."
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         end do                 !Nangles
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
      end do                    !Nenergy
      close(outUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close output file!!!"
     
      write (*,*) '*** Wrote output to file: ',TRIM(outfile)
      
      stop
      
 20   format(' ',A,I6,A,8I8,A,E24.15,SP,E25.15," I")
 30   format(' ',A,5I4,A,F20.13,SP,F21.13," I")
 40   format(A,2F18.13)
      
      end PROGRAM
