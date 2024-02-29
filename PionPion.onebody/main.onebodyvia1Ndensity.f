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
c     real*8,allocatable     :: insertion0(:,:,:),insertionz(:,:,:),insertiony(:,:,:),insertionx(:,:,:)
      real*8 aproton(3), aneutron(3), prefactor, PionMass(3), output(3)
      real*8 aPlus, aMinus
      integer*8 Znucl, extQuantNum
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
      if (nucleus.eq."3He") then
          Znucl=2
      else if(nucleus.eq."3H") then
          Znucl=1
      else if(nucleus.eq."2H") then
          Znucl=1
      else if(nucleus.eq."4He") then
          Znucl=2
      else if(nucleus.eq."6Li") then
          Znucl=3
      end if

      call ReadinputOnebody(inUnitno,calctype,variedA,descriptors,
c---- Variable to control Feynman quadrature settings------------------------
     &     Nx,
     &     verbosity)

      call ReadinputCommonComments(descriptors,inUnitno,verbosity)
      
      close(unit=inUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close input file!!! Aborting."
c
      aPlus=1.5d0*(10**-3) ! in units of inverse pion mass
      aMinus=85.2d0*(10**-3)
      PionMass=(/139.57,134.977,139.57/)
      aproton=(/aPlus+aMinus, aPlus, aPlus-aMinus/)
      aneutron=(/aPlus-aMinus, aPlus, aPlus+aMinus/)
      
      prefactor=(1+(PionMass/938))/(1+PionMass/Mnucl)!can make 938 more accurate later
      output =prefactor*(Znucl*aproton + (Anucl-Znucl)*aneutron)

      call makeoutputfilename(outfile,calctype,nucleus,descriptors,densityFileName,variedA,
     &     Elow,Ehigh,Einterval,thetaLow,thetaHigh,thetaInterval,verbosity)
c**********************************************************************
c     FINAL PRELIMINARIES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     hgrie June 2017: keep original filename: needed for replacements of energy & angle later 
      originaldensityFileName = densityFileName
      

                  
c**********************************************************************
      open(unit=outUnitno, file=outfile,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not open output file!!! Aborting."
c**********************************************************************
c     Loop over Energies
c**********************************************************************
      close(outUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close output file!!!"
     
      write (*,*) '*** Wrote output to file: ',TRIM(outfile)
      
      stop
      
 20   format(' ',A,I6,A,8I8,A,E24.15,SP,E25.15," I")
 30   format(' ',A,5I4,A,F20.13,SP,F21.13," I")
 40   format(A,2F18.13)
      
      end PROGRAM
