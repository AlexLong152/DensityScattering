cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              read2Ndensity : read input specific for 2N density computations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, based on read2Ndensity.f of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     based on Andreas Nogga's "template" files common-densities/2Ndensity-module/testcompdens.F90 in May 2017/2018.
c                                           and common-densities/2Ndensity-module/CompDens.F90 in May 2017/2018.
c      
c     Purpose: read in a set of Compton densities.
c              Result is available as "rho(ip12,ip12p,rindx)"
c              Provides also routines
cc                  get2Nchannum(l12,s12,j12,mt12,m12,twoMz): alpha-index of given (12) quantum numbers
cc              and
cc                  rhoindx(alpha2N,alpha2Np): entry of rho-matrix for given alpha & alphap
c
c     rhoindx(alpha2N,alpha2Np) has units of fm³ [Andreas email 24 May 2018]
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine read2Ndensity(densityFileName,Anucl,twoSnucl,omega,theta,j12max,P12MAG,AP12MAG,NP12,verbosity)
      
      USE CompDens              ! this needs module CompDens.mod

c      USE precision
      USE pwbook
      USE meshpoints
c      USE mpi_const
c      USE constants 
      USE spline
      
      USE HDF5
      USE parallel
      USE amplitudes
      USE hdf_tool
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
c**********************************************************************
      include '../common-densities/constants.def'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT vARIABLES:
c     rho(ip12,ip12p,rindx) , defined as "PUBLIC" via module CompDen.
c     ip12 and ip12p are indices of the momentum array P12mag(ip12) of the INTEGRATION,
c                          NOT of the grid with which the densities were produced!
c                                     rindx is the quantum-number combination of incoming and outgoing
c                                          (12) subsystem (alpha2N, alpha2Np)
c             omval [in fm^-1], thetaval [in rad]:  defined as "PUBLIC" via module CompDen.
c                                     are energy and angle to which 2Ndensity file applies.      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
c      
      character*500,intent(in) :: densityFileName  ! abuse of language: for Compton with density, this is name of density file
      
      integer,intent(in) :: NP12, j12max

      real*8,intent(in) :: P12MAG(NP12)  ! Input array of p12 momentum magnitudes at which 2Ndensity is needed, having NP12 entries.
      real*8,intent(in) :: AP12MAG(NP12) ! Input array of p12 momentum weights at which 2Ndensity is needed, having NP12 entries.
                                ! Must be passed into this routine so that CompDens.F90 subroutines can calculate electric FF

      real*8,intent(in) :: omega, theta        ! energy (in MeV) and angle (in rad) of input file. Is input to subroutine
                                ! to check that they match what's available in 2Ndensity file. 
      real*8,parameter :: eps=1.d0 ! a "small parameter" to compare omega and theta to be numerically identical
      
      integer,intent(in) :: Anucl
      integer,intent(in) :: twoSnucl
      integer,intent(in) :: verbosity         ! verbosity of stdoutout
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
c      
c     for outputting 2Ndensity at particular values of momenta (NOT ALL!!) in mathematica-friendly format
c     at two momentum-pairs with indices (idxp12,idxp12+jumpidx) and (idxp12+jumpidx,idxp12) 
      character*500    :: mathoutputfilename
      integer          :: test        ! a generic integer for testing i/o
      integer          :: idxp12      ! index of first momentum ip12
      integer          :: jumpidx     ! by how much ip12 and ip12p should differ: |ip12-ip12p|=jumpidx
      integer          :: jump        ! dummy
      
      real*8,ALLOCATABLE :: omsets(:),thetasets(:) ! arrays of energies & angles in density file. Length: nsets
      integer nsets,iset ! number of  energy-angle cobinations covered by density file, counter

c      logical testtf            ! a generic logical for testing i/o
      
c     parallelisable code: master-slave relationship set in common-densities/2Ndensity-module/parallel.dat as "just the master",
c      i.e. no parallelisation. Maybe in future. But we do only 2 integrals; should be able to do that outright.
      LOGICAL,parameter :: master=.true.

c     The following variables are pre-defined as "PUBLIC" via modules in common-densities/2Ndensity-mudules: 
c     p12n => our NP12: number or momentum grid points 
c     p12p => our P12MAG: the array of grid points (not quite, little rounding relative to our P12MAG)
c     However, we will NOT use them UNLESS for cross-checking.
c     We rely on our won momentum grid, also for the FF calculation
      
c     following are definitions for possible cross checks below 
      integer rindx
      integer ip12,l12,s12,j12,t12,mt12,m12 ! automatic integers, so not mutiplied with 2, unlike twoMz 
      integer ip12p,l12p,s12p,j12p,t12p,mt12p,m12p ! automatic integers, so not mutiplied with 2, unlike twoMz
      integer twoMz,twoMzp  ! magnetic quantum numbers of in/out target nucleus, times 2. 
      integer alpha2N,alpha2Np    
      real*8 :: ffpairs,ffnn,ffnp,ffpp

      character*64  hashtag
      character*4   rowstring
      integer       rownumber
      character*500 uniquefilename
      character*500 dummy

      integer fileSize
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) "*********************** 2N DENSITY MATRIX PARAMETERS ***************************"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if densityFileName contains hashtag of 64 characters (marked by "hash="), then isolate it
            
      if (index(densityFileName,'hash=').ne.0.) then
         dummy = densityFileName(INDEX(densityFileName,'hash=')+5:)
         hashtag = dummy(:INDEX(dummy,'.'))
         if (len(hashtag).ne.64) then
            write(*,*) "ERROR: Presumed hash ",hashtag," is NOT 64 characters long. -- Exiting."
            stop
         end if
         if (VERIFY(hashtag,"12345677890abcdef").ne.0) then
            write(*,*) "ERROR: Presumed hash ",hashtag," has characters which are not compatible with a hexadecimal. -- Exiting."
            stop
         end if
      else
         hashtag = "X"
         write(*,*) "   Density filename does not contain hashtag. => Cannot provide uniquefilename."
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if densityFileName contains rownumber (marked by "row="), then isolate it: rownumber is interpreted as everything before ".h5"
         if (index(densityFileName,'row=').ne.0.) then
            dummy = densityFileName(INDEX(densityFileName,'row=')+4:)
            rowstring = dummy(:INDEX(dummy,'.h5')-1)
            if (VERIFY(rowstring," 1234567890").ne.0) then
               write(*,*) "ERROR: Presumed rownumber ",rowstring," not a number. --  Proceeding without it."
               rownumber = 0
               stop
            else
               read(rowstring,*) rownumber
               write(rowstring,'(I0.4)') rownumber
            end if
         else
            rownumber = 0
            write(*,*) "   Density filename does not contain row number."
         end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     print hashtag and rownumber information to stdout, if it exists
      if (rownumber.ne.0) write(*,*) "Density file has row number:   ",rownumber
      if (hashtag.ne."X") then
         write(*,*) "Density file has hashtag:      ",hashtag
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     find the unique filename from the hashtag via the python script and write to stdout
c     wait for up to 15 seconds for a reply from the script. If no reply, proceed without

c          call EXECUTE_COMMAND_LINE("../venv/bin/python3 ../common-densities/find_uniquefilename_from_hashtag.py "//'"'//
c      &        "'"//hashtag//" > "//hashtag//'"', WAIT=.True.,EXITSTAT=test)
c     The following doulc be used to implement a timeout after 15 seconds -- you will need to set WATI=.False. above.
c     However, the python process continues in the background and would have to be killed separately, even after fortran has ended.
c     Is something like this is needed, it may be better to implement it via a timeout in the python script,
c     see second answer in https://stackoverflow.com/questions/492519/timeout-on-a-function-call
c         do jump=1,15           ! wait for up to 15 seconds
c            call sleep(1)
            inquire(file=hashtag,size=fileSize)
c            if (fileSize.gt.0) exit
c         end do

c     if file size is big enough, then assume we deal with a legitimate python output
            if (fileSize.gt.200) then
c           identify last line of file as the uniquefilename    
               open(unit=99,file=hashtag)
               do
                  read(99,*,IOSTAT=jump) dummy
                  if (jump.eq.0) then
                     uniquefilename = dummy
                  else
                     exit
                  end if
               end do
               close(99)
c     if what is extracted as "uniquefilename" has maore than 200 characters, it's very likely the actual unique filename
            if (len(uniquefilename).ge.200) then
               write(*,*) "Density hashtag translates into unique filename: "
               write(*,*) "    ",uniquefilename
               call EXECUTE_COMMAND_LINE("rm -rf "//hashtag, WAIT=.True.,EXITSTAT=test)
            else
               write(*,*) "   Python script did not return uniquefilename. -- Proceeding without uniquefilename."
            end if
c         else
c            write(*,*) "   Python script did not return uniquefilename within 15 seconds. -- Proceeding without uniquefilename."
         end if
      end if
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     prepare MPI and distributions of processors; parameters are in parallel.dat
      CALL initparallel

c     cannot call following because defs of parameters in precision.mod clashes with our defs
c      if (verbosity.ge.1) CALL print_eps ! prints information to STDOUT

c printout statistics on the processors used; info on parallel module   
      if (verbosity.ge.1) CALL printparallel ! prints information to STDOUT
      
c initialize hdf
      CALL init_hdf 
  
c call preparation routine for bookkeeping
c this prepares the partial wave channels for the 
c NN,3N and 4N system (even in 3N runs)  
c constraints on channels are read in from book-para.dat  
      CALL preppwbook
  
c and print the tables with channel numbers 
      if ((verbosity.ge.1).and.(master)) CALL printpwbook ! prints information to STDOUT
  
c now generate an input file for meshinit that matches 
c the generated points 
  
      CALL writemeshpara(P12MAG,AP12MAG,NP12)

c call the meshpoint initialization routine
c prepares the arrays with grid points and 
c integration weights based on the definitions in meshpoints.dat
c 4N is prepared too. Should be set to one point per q4 to switch to 
c 3N mode   
      CALL meshinit

c call the meshpoint output routine to printout grids 
      if ((verbosity.ge.1).and.(master)) CALL meshprint ! prints information to STDOUT
      
c call the initialization of the amplitudes modules 
c this call mainly prepares the distribution of the 
c grid points over processors   
      CALL initamp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     open rho file for reading, read bookkeeping, dimensions etc.
c     write(*,*) ""
c     write(*,*) "in read2Ndensity.f next line is CALL openread_rho_hdffile(densityFileName)"
c     write(*,*) ""
      CALL openread_rho_hdffile(densityFileName)
c     write(*,*) ""
c     write(*,*) "in read2Ndensity.f next line is readsets"
c     write(*,*) ""
      CALL readsets(omsets,thetasets,nsets)
     
      write(*,*) 'Sets in file:'
      write(*,*) 'set number     energy [MeV]     scatt. angle [deg]'
      
      do iset=1,nsets
         write(*,'(A,I4,2F21.2)') '   ',iset,omsets(iset)*HC,thetasets(iset)*180/Pi
      end do
      
      write(*,*) "   Finished reading bookkeeping information. Starting to read data sets."
      
      write(*,*) "Finished reading dataset of 2Ndensity."
      write(*,*) "Computing Form Factors of NN pairs using densities grid, normalised to F_pairs(q²=0)=1 for A=3,4; =2S+1 else."
      write(*,*) "      q [fm^-1]      F_nn           F_np           F_pp           F_pairs"

c     find out how many energy-angle combinations in density file      
      do iset=1,nsets
         CALL readhdf_rhoset(Anucl,omsets(iset),thetasets(iset))!,rhoRead, p12Read) 
      end do

c     output pair form factors * number of pairs: cross check not using the densities grid but the integration grid and j12max
      ffpairs= 0.d0
      ffnn   = 0.d0
      ffnp   = 0.d0
      ffpp   = 0.d0

      write(*,*) "Computing Form Factors of NN Pairs using integration and j12max from input file, normalised to F_pairs(q²=0)=1."
      write(*,'(AI3A)') " Total # pairs:",Anucl*(Anucl-1)/2,"    nn             np             pp             pairs"
      do iset=1,nsets
      write(*,'(AI3AF10.5A)') " Set #",iset,": q = ",sqrt(2*omsets(iset)**2*(1-Cos(thetasets(iset))))," MeV"
         do j12=0,j12max                      ! total ang mom (12)
            do s12=0,1                        ! spin (12) subsystem
               do l12=abs(j12-s12),j12+s12    ! orbital angular mom. (12)
                  t12=(1-(-1)**(l12+s12+1))/2 ! isospin (12) subsystem
                  do m12=-j12,j12             ! total ang mom projection of in-(12)
                     do ip12=1,NP12           ! mag of momentum (12) subsystem
                        do mt12=-t12,t12      ! mag isospin (12)
                           do twoMz=-twoSnucl,twoSnucl,2 ! average over all nucleon spin projections ||
                              alpha2N = get2Nchannum(l12,s12,j12,mt12,m12,twoMz)
                              rindx   = rhoindx(alpha2N,alpha2N)
                              if (mt12.eq.1) then
                                 ffpp = ffpp + P12MAG(ip12)**2*AP12MAG(ip12)*rho(ip12,ip12,rindx)/(twoSnucl+1)
                              else if  (mt12.eq.0) then
                                 ffnp = ffnp + P12MAG(ip12)**2*AP12MAG(ip12)*rho(ip12,ip12,rindx)/(twoSnucl+1)
                              else if  (mt12.eq.-1) then
                                 ffnn = ffnn + P12MAG(ip12)**2*AP12MAG(ip12)*rho(ip12,ip12,rindx)/(twoSnucl+1)
                              end if
                              ffpairs = ffpairs + P12MAG(ip12)**2*AP12MAG(ip12)*rho(ip12,ip12,rindx)/(twoSnucl+1)
                           end do ! twoMz   
                        end do ! mt12
                     end do ! integration i12
                  end do ! m12
               end do ! l12
            end do ! s12
         end do ! j12
         write(*,'(A(4F15.6))') "  Form Factors:",ffnn,ffnp,ffpp,ffpairs
         write(*,'(A(4F15.6))') "  Eff. # pairs:",
     &               Anucl*(Anucl-1)/2*ffnn,Anucl*(Anucl-1)/2*ffnp,Anucl*(Anucl-1)/2*ffpp,Anucl*(Anucl-1)/2*ffpairs
      end do ! set
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     following could be used for verbosity and cross checks....
c     
c check alphas: output to screen (that's a LOT!!!!)     
      
c      write(*,*) "Alpha Conversion:"
c      write(*,'(A,A5,2X,4A5,2X,3A6)') '        ','alpha','l12','s12','j12','t12','mt12','m12','twoMz'    
c      do alpha2N=1,num2Nchan
c         CALL getalpha2N(alpha2N,l12,s12,j12,t12,mt12,m12,twoMz)
cc         write(*,'(A,2I5,2X,4I5,2X,3I6)') 'ALPHA2N: ',alpha2N,l12,s12,j12,t12,mt12,m12,twoMz
c      end do   
c      do alpha2Np=1,num2Nchan
c         CALL getalpha2N(alpha2Np,l12p,s12p,j12p,t12p,mt12p,m12p,twoMzp)
cc         write(*,'(A,2I5,2X,4I5,2X,3I6)') 'ALPHA2Np: ',alpha2Np,l12p,s12p,j12p,t12p,mt12p,m12p,twoMzp
cc      end do
c
cc     check MEs
c
c      write(*,*) "P12MAG = ",P12MAG
c      write(*,*) "p12p   = ",p12p
c      write(*,*) "p12n   = ",p12n
c      do alpha2N=1,num2Nchan
c         do alpha2Np=1,num2Nchan
c            rindx=rhoindx(alpha2N,alpha2Np)
c            write(*,*) rindx
c            if(rindx.NE.0) THEN ! if rindx = 0 no coupling of channels, or matrix element not available 
c               do ip12=1,P12N   ! dimension of momentum grid P12P
c                  do ip12p=1,P12N ! dimension of momentum grid P12P
c                     call getalpha2N(alpha2N,l12,s12,j12,t12,mt12,m12,twoMz)
c                     call getalpha2N(alpha2Np,l12p,s12p,j12p,t12p,mt12p,m12p,twoMzp)
c                     if (rho(ip12,ip12p,rindx).ne.0) then
c                        write(*,'(A,I5,2X,4I5,2X,3I6)') 'ALPHA2N:  ',alpha2N,l12,s12,j12,t12,mt12,m12,twoMz
c                        write(*,'(A,I5,2X,4I5,2X,3I6)') 'ALPHA2Np: ',alpha2Np,l12p,s12p,j12p,t12p,mt12p,m12p,twoMzp
c                        write (*,*) "       rindx = ",rindx," momenta: ",p12p(ip12)," ",p12p(ip12p),"  ;   ρ12= ",
c     &                       rho(ip12,ip12p,rindx)
c                     end if   
c                  end do
c               end do
c            end if              ! rindx
c         end do                 ! alpha 2Np 
c      end do                    ! alpha 2N
c
cc         write(*,*) rho(1,1,2)
cc         write(*,*) rho(1,1,533)
cc      write(*,*) "XXXXXXXXXXXXXXXXXXXX COUNTER ",6,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
cc         write(*,*) omval
cc         write(*,*) thetaval*180.0d0/Pi
c
      
      if (verbosity.eq.4) then
c     for large verbosity number, create file "2Ndensity-formath.XXX.m" of the density which can be read in by mathematica
         mathoutputfilename = TRIM('2Ndensity-for-mathematica.'//
     &        TRIM(densityFileName(INDEX(densityFileName,'/',back=.True.)+1:INDEX(densityFileName,'.',back=.True.)-1))//'.m')
         write(*,*) '   ******* Write 2Ndensities to mathematica-readable file: '
         open(unit=15, file=mathoutputfilename,iostat=test)
         if (test .ne. 0) stop "*** ERROR: Could not open output file!!! Aborting."
c     following makes sure we get nice output that can be pasted into mathematica
         write(*,*) '           number of points in density grid: ',SIZE(P12P_density)
         write(15,*) 'ToExpression[StringReplace["{'
         write(15,'(A)')
     &        " {{ ip,  i},{        pp,         p},index,{αp,  α},{lp,sp,jp,tp,mtp,mp,twoMp},{ l, s, j, t,mt,m,twoM},  ρ2N}"
c     uncomment following if you want range of momenta -- would be VERY sizable file!!!
c         do ip12=7,7     ! 1,P12N   ! unprimed momentum grid ****SET HERE WHAT YOU WANT***
c           do ip12p=4,4        ! 1,P12N   ! primed momentum grid   ****SET HERE WHAT YOU WANT***
c     ****SET NOW WHAT YOU WANT***: momentum-pairs with indices (idxp12,idxp12+jumpidx) and (idxp12+jumpidx,idxp12)
         idxp12  = 2  ! ****SET HERE WHAT YOU WANT***
         jumpidx = 1  ! ****SET HERE WHAT YOU WANT***
         write(*,*) '           Selected grid points:'
c     check if these grid points are actually inside the rage of allowed grid points
         if ((idxp12.gt.NP12).or.((idxp12-jumpidx).gt.NP12).or.(idxp12.lt.1).or.((idxp12-jumpidx).lt.1)) then
            write(*,*) " **** Pair of requested momenta not inside input file -- ABORT. ****"
            stop
         end if
c     now produce output
         do jump=0,jumpidx,jumpidx
            ip12p = idxp12 - jump
            ip12  = idxp12 - ( jumpidx - jump )
            write(*,*) "     indices of momentum pair: (",ip12p,",",ip12,")"
               do alpha2N=1,num2Nchan ! unprimed quantum #s
                  do alpha2Np=1,num2Nchan ! primed quantum #s
                     rindx=rhoindx(alpha2N,alpha2Np)
                     if(rindx.NE.0) THEN ! if rindx = 0 no coupling of channels, or matrix element not available
c     produce quantum nimbers for channel
                        call getalpha2N(alpha2N,l12,s12,j12,t12,mt12,m12,twoMz)
                        call getalpha2N(alpha2Np,l12p,s12p,j12p,t12p,mt12p,m12p,twoMzp)
                        if (rho(ip12,ip12p,rindx).ne.0) then
                           write(15,25)
     &                          ',{{',ip12p,ip12,'},{',p12p(ip12p),p12p(ip12),'},', ! {ip12p,ip12},{p12p,p12},
     &                          rindx,',{',alpha2Np,alpha2N,'},{', ! rindx,{alpha2np,alpha2n},
     &                          l12p,s12p,j12p,t12p,mt12p,m12p,twoMzp,'},{', ! primed quantum #s from alpha12np
     &                          l12 ,s12 ,j12 ,t12 ,mt12 ,m12 ,twoMz, '},', ! unprimed quantum #s from alpha12n
     &                          rho(ip12,ip12p,rindx),'}' ! value of density
                           
                        end if   
                     end if     ! rindx
                  end do        ! alpha 2Np 
               end do           ! alpha 2N
            end do              ! index-pairing
c            end do              ! ip12p
c         end do                 ! ip12
         write(15,*) '}","E" -> "*10^"]]'
         close(unit=15,iostat=test)
         if (test .ne. 0) stop "*** ERROR: Could not close output file of 2N densities!!! Aborting."
         write(*,*) '           file name: ',mathoutputfilename
         write(*,*) '           with format:'
         write(*,'(A)')
     &        " {{ ip,  i},{     p(ip),      p(i)},index,{αp,  α},{lp,sp,jp,tp,mtp,mp,twoMp},{ l, s, j, t,mt,m,twoM},  ρ2N}"
         write(*,*) '   ******* Successfully written 2Ndensities to mathematica-readable file.'
      end if                    ! verbosity = 4
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     perform checks that density file parameters match requested parameters -- exit if not.

         if (nsets.ne.1) then
            write(*,*) "*** ERROR: 2Ndensity file contains more than 1 energy/angle combination --- "
            write(*,*)             " --- THAT IS NOT YET IMPLEMENTED ---"
            write(*,*) "-- Exiting."
            stop
         end if   
      
         if (abs(omega-omval*hc).ge.eps) then
            write(*,*) "*** ERROR: Energy of 2Ndensity file does not match:"
            write(*,*)             "    omega = ",omega," MeV   <===>   omega2N = ",omval*hc," MeV"
            write(*,'(A,F15.6,A)') "    omega - omega(2Ndensity) = ",omega-omval*hc," MeV"
            write(*,*) "-- Exiting."
            stop
         else if (abs(theta-thetaval).ge.eps) then
            write(*,*) "*** ERROR: Angle of 2Ndensity file does not match:"
            write(*,*)             "    theta = ",theta*180.0d0/Pi," deg   <===>   theta2N = ",thetaval*180.0d0/Pi," deg"
            write(*,'(A,F15.6,A)') "    theta - theta(2Ndensity) = ",(theta-thetaval)*180.0d0/Pi," deg"
            write(*,*) "-- Exiting."
            stop
         else
            write(*,*)             "   Energy and angle of 2Ndensity file match request." 
         end if

         if (j12max_rho.lt.j12max) then
            write(*,*) "*** ERROR: j12max_rho of 2Ndensity file smaller than requested. You should reduce j12max!"
            write(*,'(A,I4,A,I4)')             "    j12max = ",j12max,"   <===>   j12max2N = ",j12max_rho
            write(*,*) "-- Exiting."
            stop
         else
            write(*,'(A,I4)')             "      j12max     = ",j12max
            write(*,'(A,I4)')             "      j12max2N   = ",j12max_rho
            write(*,'(A,I4)')             "    2Ndensity file covers requested j12max. Proceeding with j12max = ",j12max    
         end if   

c     hgrie Aug 2020: deallocate the PUBLIC() ALLOCATABLE() arrays in meshpoints
c         -- needs to be done here since PUBLIC does not span all the way into main.*.f         
         if (allocated(p12p)) deallocate(p12p)
         if (allocated(p12w)) deallocate(p12w)
         if (allocated(p3p)) deallocate(p3p)
         if (allocated(p3w)) deallocate(p3w)
         if (allocated(q4p)) deallocate(q4p)
         if (allocated(q4w)) deallocate(q4w)
         if (allocated(qp)) deallocate(qp)
         if (allocated(qw)) deallocate(qw)
         if (allocated(p34p)) deallocate(p34p)
         if (allocated(p34w)) deallocate(p34w)
         
         if (allocated(xintp)) deallocate(xintp)
         if (allocated(xintw)) deallocate(xintw)
         if (allocated(phiintp)) deallocate(phiintp)
         if (allocated(phiintw)) deallocate(phiintw)
         if (allocated(xintpolp)) deallocate(xintpolp)
         if (allocated(xintpolw)) deallocate(xintpolw)
         
         if (allocated(r12p)) deallocate(r12p)
         if (allocated(r12w)) deallocate(r12w)
         if (allocated(r3p)) deallocate(r3p)
         if (allocated(r3w)) deallocate(r3w)
         if (allocated(r4p)) deallocate(r4p)
         if (allocated(r4w)) deallocate(r4w)
         if (allocated(rp)) deallocate(rp)
         if (allocated(rw)) deallocate(rw)
         if (allocated(r34p)) deallocate(r34p)
         if (allocated(r34w)) deallocate(r34w)

c         
         if (verbosity.eq.1000) continue
         
      return
c The format is tweaked such that it synchronises with the "header",    || <= here and || <= here
      
 25   format(A,I3,",",I3,A,F10.6,",",F10.6,A,I4,A,I3,",",I3,A,6(I2,","),I6,A,6(I2,","),I3,A,E24.15,A)

      end subroutine


      subroutine read2Ndensity_nointerp(densityFileName,Anucl,twoSnucl,omega,theta,j12max,P12MAG,AP12MAG,NP12,verbosity)
      
      USE CompDens              ! this needs module CompDens.mod

c      USE precision
      USE pwbook
      USE meshpoints
c      USE mpi_const
c      USE constants 
      USE spline
      
      USE HDF5
      USE parallel
      USE amplitudes
      USE hdf_tool
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
c**********************************************************************
      include '../common-densities/constants.def'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT vARIABLES:
c     rho(ip12,ip12p,rindx) , defined as "PUBLIC" via module CompDen.
c     ip12 and ip12p are indices of the momentum array P12mag(ip12) of the INTEGRATION,
c                          NOT of the grid with which the densities were produced!
c                                     rindx is the quantum-number combination of incoming and outgoing
c                                          (12) subsystem (alpha2N, alpha2Np)
c             omval [in fm^-1], thetaval [in rad]:  defined as "PUBLIC" via module CompDen.
c                                     are energy and angle to which 2Ndensity file applies.      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
c      
      character*500,intent(in) :: densityFileName  ! abuse of language: for Compton with density, this is name of density file
      
      integer,intent(in) :: NP12, j12max

      real*8,intent(in) :: P12MAG(NP12)  ! Input array of p12 momentum magnitudes at which 2Ndensity is needed, having NP12 entries.
      real*8,intent(in) :: AP12MAG(NP12) ! Input array of p12 momentum weights at which 2Ndensity is needed, having NP12 entries.
                                ! Must be passed into this routine so that CompDens.F90 subroutines can calculate electric FF

      real*8,intent(in) :: omega, theta        ! energy (in MeV) and angle (in rad) of input file. Is input to subroutine
                                ! to check that they match what's available in 2Ndensity file. 
      real*8,parameter :: eps=1.d0 ! a "small parameter" to compare omega and theta to be numerically identical
      
      integer,intent(in) :: Anucl
      integer,intent(in) :: twoSnucl
      integer,intent(in) :: verbosity         ! verbosity of stdoutout
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
c      
c     for outputting 2Ndensity at particular values of momenta (NOT ALL!!) in mathematica-friendly format
c     at two momentum-pairs with indices (idxp12,idxp12+jumpidx) and (idxp12+jumpidx,idxp12) 
      character*500    :: mathoutputfilename
      integer          :: test        ! a generic integer for testing i/o
      integer          :: idxp12      ! index of first momentum ip12
      integer          :: jumpidx     ! by how much ip12 and ip12p should differ: |ip12-ip12p|=jumpidx
      integer          :: jump        ! dummy
      
      real*8,ALLOCATABLE :: omsets(:),thetasets(:) ! arrays of energies & angles in density file. Length: nsets
      integer nsets,iset ! number of  energy-angle cobinations covered by density file, counter

c      logical testtf            ! a generic logical for testing i/o
      
c     parallelisable code: master-slave relationship set in common-densities/2Ndensity-module/parallel.dat as "just the master",
c      i.e. no parallelisation. Maybe in future. But we do only 2 integrals; should be able to do that outright.
      LOGICAL,parameter :: master=.true.

c     The following variables are pre-defined as "PUBLIC" via modules in common-densities/2Ndensity-mudules: 
c     p12n => our NP12: number or momentum grid points 
c     p12p => our P12MAG: the array of grid points (not quite, little rounding relative to our P12MAG)
c     However, we will NOT use them UNLESS for cross-checking.
c     We rely on our won momentum grid, also for the FF calculation
      
c     following are definitions for possible cross checks below 
      integer rindx
      integer ip12,l12,s12,j12,t12,mt12,m12 ! automatic integers, so not mutiplied with 2, unlike twoMz 
      integer ip12p,l12p,s12p,j12p,t12p,mt12p,m12p ! automatic integers, so not mutiplied with 2, unlike twoMz
      integer twoMz,twoMzp  ! magnetic quantum numbers of in/out target nucleus, times 2. 
      integer alpha2N,alpha2Np    
      real*8 :: ffpairs,ffnn,ffnp,ffpp

      character*64  hashtag
      character*4   rowstring
      integer       rownumber
      character*500 uniquefilename
      character*500 dummy

      integer fileSize
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) "*********************** 2N DENSITY MATRIX PARAMETERS ***************************"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if densityFileName contains hashtag of 64 characters (marked by "hash="), then isolate it
      write(*,*) "densityFileName=", densityFileName 
      if (index(densityFileName,'hash=').ne.0.) then
         dummy = densityFileName(INDEX(densityFileName,'hash=')+5:)
         hashtag = dummy(:INDEX(dummy,'.'))
         if (len(hashtag).ne.64) then
            write(*,*) "ERROR: Presumed hash ",hashtag," is NOT 64 characters long. -- Exiting."
            stop
         end if
         if (VERIFY(hashtag,"12345677890abcdef").ne.0) then
            write(*,*) "ERROR: Presumed hash ",hashtag," has characters which are not compatible with a hexadecimal. -- Exiting."
            stop
         end if
      else
         hashtag = "X"
         write(*,*) "   Density filename does not contain hashtag. => Cannot provide uniquefilename."
      end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if densityFileName contains rownumber (marked by "row="), then isolate it: rownumber is interpreted as everything before ".h5"
         if (index(densityFileName,'row=').ne.0.) then
            dummy = densityFileName(INDEX(densityFileName,'row=')+4:)
            rowstring = dummy(:INDEX(dummy,'.h5')-1)
            if (VERIFY(rowstring," 1234567890").ne.0) then
               write(*,*) "ERROR: Presumed rownumber ",rowstring," not a number. --  Proceeding without it."
               rownumber = 0
               stop
            else
               read(rowstring,*) rownumber
               write(rowstring,'(I0.4)') rownumber
            end if
         else
            rownumber = 0
            write(*,*) "   Density filename does not contain row number."
         end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     print hashtag and rownumber information to stdout, if it exists
      if (rownumber.ne.0) write(*,*) "Density file has row number:   ",rownumber
      if (hashtag.ne."X") then
         write(*,*) "Density file has hashtag:      ",hashtag
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     find the unique filename from the hashtag via the python script and write to stdout
c     wait for up to 15 seconds for a reply from the script. If no reply, proceed without
c          call EXECUTE_COMMAND_LINE("../venv/bin/python3 ../common-densities/find_uniquefilename_from_hashtag.py "//
c      &        "'"//hashtag//" > "//hashtag, WAIT=.True.,EXITSTAT=test)
c     The following doulc be used to implement a timeout after 15 seconds -- you will need to set WATI=.False. above.
c     However, the python process continues in the background and would have to be killed separately, even after fortran has ended.
c     Is something like this is needed, it may be better to implement it via a timeout in the python script,
c     see second answer in https://stackoverflow.com/questions/492519/timeout-on-a-function-call
c         do jump=1,15           ! wait for up to 15 seconds
c            call sleep(1)
            inquire(file=hashtag,size=fileSize)
c            if (fileSize.gt.0) exit
c         end do

c     if file size is big enough, then assume we deal with a legitimate python output
            if (fileSize.gt.200) then
c           identify last line of file as the uniquefilename    
               open(unit=99,file=hashtag)
               do
                  read(99,*,IOSTAT=jump) dummy
                  if (jump.eq.0) then
                     uniquefilename = dummy
                  else
                     exit
                  end if
               end do
               close(99)
c     if what is extracted as "uniquefilename" has maore than 200 characters, it's very likely the actual unique filename
            if (len(uniquefilename).ge.200) then
               write(*,*) "Density hashtag translates into unique filename: "
               write(*,*) "    ",uniquefilename
               call EXECUTE_COMMAND_LINE("rm -rf "//hashtag, WAIT=.True.,EXITSTAT=test)
            else
               write(*,*) "   Python script did not return uniquefilename. -- Proceeding without uniquefilename."
            end if
c         else
c            write(*,*) "   Python script did not return uniquefilename within 15 seconds. -- Proceeding without uniquefilename."
         end if
      end if
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     prepare MPI and distributions of processors; parameters are in parallel.dat
      CALL initparallel

c     cannot call following because defs of parameters in precision.mod clashes with our defs
c      if (verbosity.ge.1) CALL print_eps ! prints information to STDOUT

c printout statistics on the processors used; info on parallel module   
      if (verbosity.ge.1) CALL printparallel ! prints information to STDOUT
      
c initialize hdf
      CALL init_hdf 
  
c call preparation routine for bookkeeping
c this prepares the partial wave channels for the 
c NN,3N and 4N system (even in 3N runs)  
c constraints on channels are read in from book-para.dat  
      CALL preppwbook
  
c and print the tables with channel numbers 
      if ((verbosity.ge.1).and.(master)) CALL printpwbook ! prints information to STDOUT
  
c now generate an input file for meshinit that matches 
c the generated points 
  
c     CALL writemeshpara(P12P_density,P12W_density,len(P12P_density))!Changed
      CALL writemeshpara(P12MAG,AP12MAG,NP12)

c call the meshpoint initialization routine
c prepares the arrays with grid points and 
c integration weights based on the definitions in meshpoints.dat
c 4N is prepared too. Should be set to one point per q4 to switch to 
c 3N mode   
      CALL meshinit

c call the meshpoint output routine to printout grids 
      if ((verbosity.ge.1).and.(master)) CALL meshprint ! prints information to STDOUT
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     open rho file for reading, read bookkeeping, dimensions etc.
c     write(*,*) ""
c     write(*,*) "in read2Ndensity.f next line is CALL openread_rho_hdffile(densityFileName)"
c     write(*,*) ""

      CALL openread_rho_hdffile(densityFileName)!openread_rho_hdffile defindes P12P_density
c     write(*,*) "In read2Ndensity.f after first call: P12P_density=",P12P_density
c     write(*,*) "In read2Ndensity.f after first call: size(P12P_density)=",size(P12P_density)
c     call execute_command_line("cat meshpoints.dat",wait=.True.)
      CALL writemeshpara(P12P_density,P12W_density,size(P12P_density))!Now overwrite the old meshpoints with the new ones
c     write(*,*) "meshpoints after new call"
c     write(*,*) ""
c     call execute_command_line("cat meshpoints.dat",wait=.True.)
      write(*,*) "In read2Ndensity.f: shape(rho)=",shape(rho) 
c     write(*,*) ""
c     write(*,*) "in read2Ndensity.f next line is readsets"
c     write(*,*) ""
      CALL readsets(omsets,thetasets,nsets)
     
c call the initialization of the amplitudes modules 
c this call mainly prepares the distribution of the 
c grid points over processors   
      CALL initamp

      write(*,*) 'Sets in file:'
      write(*,*) 'set number     energy [MeV]     scatt. angle [deg]'
      
      do iset=1,nsets
         write(*,'(A,I4,2F21.2)') '   ',iset,omsets(iset)*HC,thetasets(iset)*180/Pi
      end do
      
      write(*,*) "   Finished reading bookkeeping information. Starting to read data sets."
      
      write(*,*) "Finished reading dataset of 2Ndensity."
      write(*,*) "Computing Form Factors of NN pairs using densities grid, normalised to F_pairs(q²=0)=1 for A=3,4; =2S+1 else."
      write(*,*) "      q [fm^-1]      F_nn           F_np           F_pp           F_pairs"

c     find out how many energy-angle combinations in density file      
      do iset=1,nsets
         CALL readhdf_rhoset(Anucl,omsets(iset),thetasets(iset))!,rhoRead, p12Read) 
      end do

c     output pair form factors * number of pairs: cross check not using the densities grid but the integration grid and j12max
      ffpairs= 0.d0
      ffnn   = 0.d0
      ffnp   = 0.d0
      ffpp   = 0.d0

      write(*,*) "Computing Form Factors of NN Pairs using integration and j12max from input file, normalised to F_pairs(q²=0)=1."
      write(*,'(AI3A)') " Total # pairs:",Anucl*(Anucl-1)/2,"    nn             np             pp             pairs"
      do iset=1,nsets
      write(*,'(AI3AF10.5A)') " Set #",iset,": q = ",sqrt(2*omsets(iset)**2*(1-Cos(thetasets(iset))))," MeV"
         do j12=0,j12max                      ! total ang mom (12)
            do s12=0,1                        ! spin (12) subsystem
               do l12=abs(j12-s12),j12+s12    ! orbital angular mom. (12)
                  t12=(1-(-1)**(l12+s12+1))/2 ! isospin (12) subsystem
                  do m12=-j12,j12             ! total ang mom projection of in-(12)
                     do ip12=1,NP12           ! mag of momentum (12) subsystem
                        do mt12=-t12,t12      ! mag isospin (12)
                           do twoMz=-twoSnucl,twoSnucl,2 ! average over all nucleon spin projections ||
                              alpha2N = get2Nchannum(l12,s12,j12,mt12,m12,twoMz)
                              rindx   = rhoindx(alpha2N,alpha2N)
                              if (mt12.eq.1) then
                                 ffpp = ffpp + P12MAG(ip12)**2*AP12MAG(ip12)*rho(ip12,ip12,rindx)/(twoSnucl+1)
                              else if  (mt12.eq.0) then
                                 ffnp = ffnp + P12MAG(ip12)**2*AP12MAG(ip12)*rho(ip12,ip12,rindx)/(twoSnucl+1)
                              else if  (mt12.eq.-1) then
                                 ffnn = ffnn + P12MAG(ip12)**2*AP12MAG(ip12)*rho(ip12,ip12,rindx)/(twoSnucl+1)
                              end if
                              ffpairs = ffpairs + P12MAG(ip12)**2*AP12MAG(ip12)*rho(ip12,ip12,rindx)/(twoSnucl+1)
                           end do ! twoMz   
                        end do ! mt12
                     end do ! integration i12
                  end do ! m12
               end do ! l12
            end do ! s12
         end do ! j12
         write(*,'(A(4F15.6))') "  Form Factors:",ffnn,ffnp,ffpp,ffpairs
         write(*,'(A(4F15.6))') "  Eff. # pairs:",
     &               Anucl*(Anucl-1)/2*ffnn,Anucl*(Anucl-1)/2*ffnp,Anucl*(Anucl-1)/2*ffpp,Anucl*(Anucl-1)/2*ffpairs
      end do ! set
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
      
      if (verbosity.eq.4) then
c     for large verbosity number, create file "2Ndensity-formath.XXX.m" of the density which can be read in by mathematica
         mathoutputfilename = TRIM('2Ndensity-for-mathematica.'//
     &        TRIM(densityFileName(INDEX(densityFileName,'/',back=.True.)+1:INDEX(densityFileName,'.',back=.True.)-1))//'.m')
         write(*,*) '   ******* Write 2Ndensities to mathematica-readable file: '
         open(unit=15, file=mathoutputfilename,iostat=test)
         if (test .ne. 0) stop "*** ERROR: Could not open output file!!! Aborting."
c     following makes sure we get nice output that can be pasted into mathematica
         write(*,*) '           number of points in density grid: ',SIZE(P12P_density)
         write(15,*) 'ToExpression[StringReplace["{'
         write(15,'(A)')
     &        " {{ ip,  i},{        pp,         p},index,{αp,  α},{lp,sp,jp,tp,mtp,mp,twoMp},{ l, s, j, t,mt,m,twoM},  ρ2N}"
c     uncomment following if you want range of momenta -- would be VERY sizable file!!!
c         do ip12=7,7     ! 1,P12N   ! unprimed momentum grid ****SET HERE WHAT YOU WANT***
c           do ip12p=4,4        ! 1,P12N   ! primed momentum grid   ****SET HERE WHAT YOU WANT***
c     ****SET NOW WHAT YOU WANT***: momentum-pairs with indices (idxp12,idxp12+jumpidx) and (idxp12+jumpidx,idxp12)
         idxp12  = 2  ! ****SET HERE WHAT YOU WANT***
         jumpidx = 1  ! ****SET HERE WHAT YOU WANT***
         write(*,*) '           Selected grid points:'
c     check if these grid points are actually inside the rage of allowed grid points
         if ((idxp12.gt.NP12).or.((idxp12-jumpidx).gt.NP12).or.(idxp12.lt.1).or.((idxp12-jumpidx).lt.1)) then
            write(*,*) " **** Pair of requested momenta not inside input file -- ABORT. ****"
            stop
         end if
c     now produce output
         do jump=0,jumpidx,jumpidx
            ip12p = idxp12 - jump
            ip12  = idxp12 - ( jumpidx - jump )
            write(*,*) "     indices of momentum pair: (",ip12p,",",ip12,")"
               do alpha2N=1,num2Nchan ! unprimed quantum #s
                  do alpha2Np=1,num2Nchan ! primed quantum #s
                     rindx=rhoindx(alpha2N,alpha2Np)
                     if(rindx.NE.0) THEN ! if rindx = 0 no coupling of channels, or matrix element not available
c     produce quantum nimbers for channel
                        call getalpha2N(alpha2N,l12,s12,j12,t12,mt12,m12,twoMz)
                        call getalpha2N(alpha2Np,l12p,s12p,j12p,t12p,mt12p,m12p,twoMzp)
                        if (rho(ip12,ip12p,rindx).ne.0) then
                           write(15,25)
     &                          ',{{',ip12p,ip12,'},{',p12p(ip12p),p12p(ip12),'},', ! {ip12p,ip12},{p12p,p12},
     &                          rindx,',{',alpha2Np,alpha2N,'},{', ! rindx,{alpha2np,alpha2n},
     &                          l12p,s12p,j12p,t12p,mt12p,m12p,twoMzp,'},{', ! primed quantum #s from alpha12np
     &                          l12 ,s12 ,j12 ,t12 ,mt12 ,m12 ,twoMz, '},', ! unprimed quantum #s from alpha12n
     &                          rho(ip12,ip12p,rindx),'}' ! value of density
                           
                        end if   
                     end if     ! rindx
                  end do        ! alpha 2Np 
               end do           ! alpha 2N
            end do              ! index-pairing
c            end do              ! ip12p
c         end do                 ! ip12
         write(15,*) '}","E" -> "*10^"]]'
         close(unit=15,iostat=test)
         if (test .ne. 0) stop "*** ERROR: Could not close output file of 2N densities!!! Aborting."
         write(*,*) '           file name: ',mathoutputfilename
         write(*,*) '           with format:'
         write(*,'(A)')
     &        " {{ ip,  i},{     p(ip),      p(i)},index,{αp,  α},{lp,sp,jp,tp,mtp,mp,twoMp},{ l, s, j, t,mt,m,twoM},  ρ2N}"
         write(*,*) '   ******* Successfully written 2Ndensities to mathematica-readable file.'
      end if                    ! verbosity = 4
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     perform checks that density file parameters match requested parameters -- exit if not.

         if (nsets.ne.1) then
            write(*,*) "*** ERROR: 2Ndensity file contains more than 1 energy/angle combination --- "
            write(*,*)             " --- THAT IS NOT YET IMPLEMENTED ---"
            write(*,*) "-- Exiting."
            stop
         end if   
      
         if (abs(omega-omval*hc).ge.eps) then
            write(*,*) "*** ERROR: Energy of 2Ndensity file does not match:"
            write(*,*)             "    omega = ",omega," MeV   <===>   omega2N = ",omval*hc," MeV"
            write(*,'(A,F15.6,A)') "    omega - omega(2Ndensity) = ",omega-omval*hc," MeV"
            write(*,*) "-- Exiting."
            stop
         else if (abs(theta-thetaval).ge.eps) then
            write(*,*) "*** ERROR: Angle of 2Ndensity file does not match:"
            write(*,*)             "    theta = ",theta*180.0d0/Pi," deg   <===>   theta2N = ",thetaval*180.0d0/Pi," deg"
            write(*,'(A,F15.6,A)') "    theta - theta(2Ndensity) = ",(theta-thetaval)*180.0d0/Pi," deg"
            write(*,*) "-- Exiting."
            stop
         else
            write(*,*)             "   Energy and angle of 2Ndensity file match request." 
         end if

         if (j12max_rho.lt.j12max) then
            write(*,*) "*** ERROR: j12max_rho of 2Ndensity file smaller than requested. You should reduce j12max!"
            write(*,'(A,I4,A,I4)')             "    j12max = ",j12max,"   <===>   j12max2N = ",j12max_rho
            write(*,*) "-- Exiting."
            stop
         else
            write(*,'(A,I4)')             "      j12max     = ",j12max
            write(*,'(A,I4)')             "      j12max2N   = ",j12max_rho
            write(*,'(A,I4)')             "    2Ndensity file covers requested j12max. Proceeding with j12max = ",j12max    
         end if   

c     hgrie Aug 2020: deallocate the PUBLIC() ALLOCATABLE() arrays in meshpoints
c         -- needs to be done here since PUBLIC does not span all the way into main.*.f         
         if (allocated(p12p)) deallocate(p12p)
         if (allocated(p12w)) deallocate(p12w)
         if (allocated(p3p)) deallocate(p3p)
         if (allocated(p3w)) deallocate(p3w)
         if (allocated(q4p)) deallocate(q4p)
         if (allocated(q4w)) deallocate(q4w)
         if (allocated(qp)) deallocate(qp)
         if (allocated(qw)) deallocate(qw)
         if (allocated(p34p)) deallocate(p34p)
         if (allocated(p34w)) deallocate(p34w)
         
         if (allocated(xintp)) deallocate(xintp)
         if (allocated(xintw)) deallocate(xintw)
         if (allocated(phiintp)) deallocate(phiintp)
         if (allocated(phiintw)) deallocate(phiintw)
         if (allocated(xintpolp)) deallocate(xintpolp)
         if (allocated(xintpolw)) deallocate(xintpolw)
         
         if (allocated(r12p)) deallocate(r12p)
         if (allocated(r12w)) deallocate(r12w)
         if (allocated(r3p)) deallocate(r3p)
         if (allocated(r3w)) deallocate(r3w)
         if (allocated(r4p)) deallocate(r4p)
         if (allocated(r4w)) deallocate(r4w)
         if (allocated(rp)) deallocate(rp)
         if (allocated(rw)) deallocate(rw)
         if (allocated(r34p)) deallocate(r34p)
         if (allocated(r34w)) deallocate(r34w)

c         
         if (verbosity.eq.1000) continue
         
      return
c The format is tweaked such that it synchronises with the "header",    || <= here and || <= here
      
 25   format(A,I3,",",I3,A,F10.6,",",F10.6,A,I4,A,I3,",",I3,A,6(I2,","),I6,A,6(I2,","),I3,A,E24.15,A)

      end subroutine
