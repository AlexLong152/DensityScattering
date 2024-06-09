C     define macros for compiling optional code 
#include "fdefs.h" 

C     this module is setting up the bookkeeping
C     for the few body scattering code

      MODULE pwbook
       USE precision 
       USE mpi_const
       USE HDF5
       USE hdf_tool
       USE parallel 
       
       PRIVATE

C     module is not parallelized, but timing routines
C     from MPI are called
       
C     arrays and variables for partial wave channels of NN

       INTEGER,PUBLIC :: alphaNNmax,alphaNNcdepmax,l12max,
     $                   mt12min,mt12max  
       INTEGER :: j12max
       LOGICAL,PUBLIC :: evenNN,cdepNN
       INTEGER,ALLOCATABLE :: qnalphaNN(:,:)

C     arrays and variables for partial wave channels of 3N

       INTEGER,PUBLIC :: alpha3Nmax,alpha3Ncdepmax,alpha3NLScdepmax,l3max,
     $                   I3max,j3min,j3max,
     $                   tau3min,tau3max,pari3Nmin,pari3Nmax,
     $                   mtau3min,mtau3max
       LOGICAL,PUBLIC :: cdep3N
       INTEGER,ALLOCATABLE :: qnalpha3N(:,:),qnalpha3NLS(:,:)

C     arrays and variables for partial wave channels of 4N
C> @page 4N bookkeeping
C> @section beta4N22max
C> number of 2+2 channels without charge dependence
C> @section beta4N22cdepmax
C> number of 2+2 channels including charge dependence
C> @section alpha4N31max
C> number of 3+1 channels without charge dependence
C> @section alpha4N31cdepmax
C> number of 3+1 channels including charge dependence
C> @section l4max
C> maximum l4
C> @section I4max
C> twice the maximum I4 
C> @section j4min
C> minimum total J
C> @section j4max
C> maximum total J 
C> @section lammax
C> maximum lambda 
C> @section Imax
C> maximum I 
C> @section tau4min
C> minimum tau4 
C> @section tau4max
C> maximum tau4 
C> @section pari4Nmin
C> minimal parity +-1
C> @section pari4Nmax
C> maximal parity +-1
C> @section mtau4min
C> minimum third component of total isospin 
C> @section mtau4max
C> maximum third component of total isospin
C> @section lsummax31
C> maximum of the sum of l12+l3+l4 
C> @section lsummax22
C> maximum of the sum of l12+l34+lambda 
C> @section cdep4N
C> flag to allow for charge dependendent 4N states 
C> @section qnbeta4N22
C> qnbeta4N22(1:16,1:\ref beta4N22cdepmax) contains channel information <br>
C>               - qnbeta4N22(1,beta) = l12 
C>               - qnbeta4N22(2,beta) = s12
C>               - qnbeta4N22(3,beta) = j12
C>               - qnbeta4N22(4,beta) = t12
C>               - qnbeta4N22(5,beta) = l34
C>               - qnbeta4N22(6,beta) = s34
C>               - qnbeta4N22(7,beta) = j34
C>               - qnbeta4N22(8,beta) = t34
C>               - qnbeta4N22(9,beta) = lam
C>               - qnbeta4N22(10,beta)= I
C>               - qnbeta4N22(11,beta)= j4
C>               - qnbeta4N22(12,beta)= tau4
C>               - qnbeta4N22(13,beta)= mtau4
C>               - qnbeta4N22(14,beta)= alphaNN12
C>               - qnbeta4N22(15,beta)= alphaNN34
C>               - qnbeta4N22(16,beta)= pari
C> @section qnalpha4N31
C> qnalpha4N31(1:16,1:\ref alpha4N31cdepmax) contains channel information <br>
C>              - qnalpha4N31(1,alpha) = l12
C>              - qnalpha4N31(2,alpha) = s12
C>              - qnalpha4N31(3,alpha) = j12
C>              - qnalpha4N31(4,alpha) = t12
C>              - qnalpha4N31(5,alpha) = l3
C>              - qnalpha4N31(6,alpha) = twice I3
C>              - qnalpha4N31(7,alpha) = twice j3
C>              - qnalpha4N31(8,alpha) = twice tau3
C>              - qnalpha4N31(9,alpha) = l4
C>              - qnalpha4N31(10,alpha)= twice I4
C>              - qnalpha4N31(11,alpha)= j4
C>              - qnalpha4N31(12,alpha)= tau4
C>              - qnalpha4N31(13,alpha)= mtau4
C>              - qnalpha4N31(14,alpha)= alphaNN12
C>              - qnalpha4N31(15,alpha)= alpha3N
C>              - qnalpha4N31(16,alpha)= pari           
       INTEGER,PUBLIC :: beta4N22max,beta4N22cdepmax,
     $                   alpha4N31max,alpha4N31cdepmax,l4max,
     $                   I4max,j4min,j4max,lammax,Imax,
     $                   tau4min,tau4max,pari4Nmin,pari4Nmax,
     $                   mtau4min,mtau4max,lsummax31,lsummax22
       LOGICAL,PUBLIC :: cdep4N
       INTEGER,ALLOCATABLE :: qnbeta4N22(:,:),qnalpha4N31(:,:)

C     module contains two public subroutine 
C     preppwbook prepares the bookkeeping based on the parameters 
C         in the file book-para.dat
C     printpwbook prints out all channels for checking purposes 
C     getNNqn returns quantum numbers of the NN system for a specific channel  
C     get3Nqn returns quantum numbers of the 3N system for a specific channel  
C     note that prepbook opens and closes UNIT=75 for input parameters 
 
       PUBLIC preppwbook,printpwbook,getNNqn,get3Nqn,get3NLSqn,
     $        get4N31qn,get4N22qn,writennchannels,
     $        readnnchannels,write_3n_channels,read_3n_channels,
     $        write_4n31_channels,read_4n31_channels,
     $        write_4n22_channels,read_4n22_channels

C     the public variable pwbooktime measures time 
C     needed for calling preppwbook 
C     to give example of possible timing features 
C     for other modules 

       REAL(dpreal),PUBLIC :: pwbooktime=0.0

      CONTAINS

      SUBROUTINE printpwbook
       IMPLICIT NONE 

       INTEGER l12,s12,j12,t12,mt12
       INTEGER l3,I3,j3,tau3,mtau3,pari
       INTEGER l34,s34,j34,t34,l4,I4,lam,I,j4,tau4,mtau4
       INTEGER bl,bs
       INTEGER alpha,beta,alphaNN12,alphaNN34,alpha3N

       WRITE(*,*) 'Parameters of NN bookkeeping'
       WRITE(*,'(A,I5)') 'j12max        = ',j12max
       WRITE(*,'(A,I5)') 'l12max        = ',l12max
       WRITE(*,*)        'evenNN        = ',evenNN
       WRITE(*,*)        'cdepNN        = ',cdepNN
       WRITE(*,*)        'mt12 min-max  = ',mt12min,mt12max
       WRITE(*,*) 
       WRITE(*,'(A,I5)') 'alphaNNmax     = ',alphaNNmax
       WRITE(*,'(A,I5)') 'alphaNNcdepmax = ',alphaNNcdepmax
       WRITE(*,*) 


       WRITE(*,*) 
       WRITE(*,*) 'PW bookkeeping:'
       WRITE(*,*) '         alphaNN   l12 s12 j12 t12 mt12'
       WRITE(*,*) 

       DO alpha=1,alphaNNcdepmax
        call getNNqn(alpha,l12,s12,j12,t12,mt12)

        WRITE(*,'(A,I7,11I5)') 'ALPHANN: ',alpha,l12,s12,j12,t12,mt12
       END DO

       WRITE(*,*) 
       WRITE(*,*) 

       WRITE(*,*) 'Parameters of 3N bookkeeping'
       WRITE(*,'(A,I5)')  'l3max               = ',l3max
       WRITE(*,'(A,I5)')  'I3max               = ',I3max
       WRITE(*,'(A,2I5)') 'j3min,j3max         = ',j3min,j3max
       WRITE(*,'(A,2I5)') 'pari3Nmin,pari3Nmax = ',pari3Nmin,pari3Nmax
       WRITE(*,'(A,2I5)') 'tau3min,tau3max     = ',tau3min,tau3max
       WRITE(*,*)         'cdep3N              = ',cdep3N
       WRITE(*,*)         'mtau3 min-max       = ',mtau3min,mtau3max
       WRITE(*,*) 
       WRITE(*,'(A,I5)') 'alpha3Nmax     = ',alpha3Nmax
       WRITE(*,'(A,I5)') 'alpha3Ncdepmax = ',alpha3Ncdepmax
       WRITE(*,*) 


       WRITE(*,*) 
       WRITE(*,*)
     $  'PW bookkeeping:'
       WRITE(*,*)
     $  '         alpha3N   l12 s12 j12 t12 l3  I3  j3'
     $  // '  tau3  mtau3  alphaNN12  pari'
       WRITE(*,*) 

       DO alpha=1,alpha3Ncdepmax
        call get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,
     $                               tau3,mtau3,alphaNN12,pari)

        WRITE(*,'(A,I7,11I5)')
     $   'ALPHA3N: ',alpha,l12,s12,j12,t12,l3,I3,j3,
     $                                 tau3,mtau3,alphaNN12,pari
       END DO

       WRITE(*,*) 
       WRITE(*,*)
     $  'LS bookkeeping:'
       WRITE(*,'(A8,A7,10A5)') ' ','beta3N','l12','l3','L','s12','S',
     $          'j3','t12','tau3','mtau3','pari'
       WRITE(*,*) 

       DO alpha=1,alpha3NLScdepmax
        call get3NLSqn(alpha,l12,l3,bl,s12,bs,j3,t12,tau3,mtau3,pari)

        WRITE(*,'(A8,I7,10I5)')
     $   'BETA3N: ',alpha,l12,l3,bl,s12,bs,j3,t12,tau3,mtau3,pari
       END DO
       
       
       WRITE(*,*) 
       WRITE(*,*) 

       WRITE(*,*) 'Parameters of 4N 3+1 bookkeeping'
       WRITE(*,'(A,I5)')  'l4max               = ',l4max
       WRITE(*,'(A,I5)')  'I4max               = ',I4max
       WRITE(*,'(A,2I5)') 'j4min,j4max         = ',j4min,j4max
       WRITE(*,'(A,I5)')  'lsummax31           = ',lsummax31
       WRITE(*,'(A,2I5)') 'pari4Nmin,pari4Nmax = ',pari4Nmin,pari4Nmax
       WRITE(*,'(A,2I5)') 'tau4min,tau4max     = ',tau4min,tau4max
       WRITE(*,*)         'cdep4N              = ',cdep4N
       WRITE(*,*)         'mtau4 min-max       = ',mtau4min,mtau4max
       WRITE(*,*) 
       WRITE(*,'(A,I5)') 'alpha4N31max     = ',alpha4N31max
       WRITE(*,'(A,I5)') 'alpha4N31cdepmax = ',alpha4N31cdepmax
       WRITE(*,*) 


       WRITE(*,*) 
       WRITE(*,*)
     $  'PW bookkeeping:'
       WRITE(*,*)
     $  '         alpha4N   l12 s12 j12 t12 l3  I3  j3'
     $  // '  tau3  l4  I4  j4 tau4 mtau4 alphaNN12 alpha3N  pari'
       WRITE(*,*) 

       DO alpha=1,alpha4N31cdepmax
        call get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,
     $                 tau3,l4,I4,j4,tau4,mtau4,alphaNN12,alpha3N,pari)

        WRITE(*,'(A,I7,14I5,I7,I5)')
     $   'ALPHA4N31: ',alpha,l12,s12,j12,t12,l3,I3,j3,
     $               tau3,l4,I4,j4,tau4,mtau4,alphaNN12,alpha3N,pari
       END DO

       WRITE(*,*) 
       WRITE(*,*) 


       WRITE(*,*) 'Parameters of 4N 2+2 bookkeeping'
       WRITE(*,'(A,I5)')  'lammax              = ',lammax
       WRITE(*,'(A,I5)')  'Imax                = ',Imax
       WRITE(*,'(A,I5)')  'lsummax22           = ',lsummax22
       WRITE(*,*) 
       WRITE(*,'(A,I5)')  'beta4N22max          = ',beta4N22max
       WRITE(*,'(A,I5)')  'beta4N22cdepmax      = ',beta4N22cdepmax
       WRITE(*,*)  


       WRITE(*,*) 
       WRITE(*,*)
     $  'PW bookkeeping:'
       WRITE(*,*)
     $  '         beta4N   l12 s12 j12 t12 l34  s34  j34  t34'
     $  // '  lam  I  j4  tau4 mtau4 alphaNN12  alphaNN34  pari'
       WRITE(*,*) 

       DO beta=1,beta4N22cdepmax
        call get4N22qn(beta,l12,s12,j12,t12,l34,s34,j34,t34,
     $                 lam,I,j4,tau4,mtau4,alphaNN12,alphaNN34,pari)

        WRITE(*,'(A,I7,16I5)')
     $   'BETA4N22: ',beta,l12,s12,j12,t12,l34,s34,j34,t34,
     $                     lam,I,j4,tau4,mtau4,alphaNN12,alphaNN34,pari
       END DO

       WRITE(*,*) 
       WRITE(*,*) 

      END SUBROUTINE printpwbook


C     subroutine for obtaining the NN quantum numbers of NN channel alpha

      SUBROUTINE getNNqn(alpha,l12,s12,j12,t12,mt12)
       IMPLICIT NONE 
       INTEGER alpha,l12,s12,j12,t12,mt12

#ifdef DEBUG
C       check for consistency only then debugging 
       IF(alpha.LE.0) STOP 'alpha too small'
       IF(alpha.GT.alphaNNcdepmax) STOP 'alpha too large'
#endif 

       l12=qnalphaNN(1,alpha)
       s12=qnalphaNN(2,alpha)
       j12=qnalphaNN(3,alpha)
       t12=qnalphaNN(4,alpha)
       mt12=qnalphaNN(5,alpha)

       
      END SUBROUTINE getNNqn

      SUBROUTINE get3Nqn(alpha,l12,s12,j12,t12,l3,I3,j3,
     $                                     tau3,mtau3,alphaNN,pari)
       IMPLICIT NONE 
       INTEGER alpha,l12,s12,j12,t12,l3,I3,j3,tau3,mtau3,alphaNN,pari

#ifdef DEBUG
C       check for consistency only then debugging 
       IF(alpha.LE.0) STOP 'alpha too small'
       IF(alpha.GT.alpha3Ncdepmax) STOP 'alpha too large'
#endif 

       l12=qnalpha3N(1,alpha)
       s12=qnalpha3N(2,alpha)
       j12=qnalpha3N(3,alpha)
       t12=qnalpha3N(4,alpha)
       l3=qnalpha3N(5,alpha)
       I3=qnalpha3N(6,alpha)
       j3=qnalpha3N(7,alpha)
       tau3=qnalpha3N(8,alpha)
       mtau3=qnalpha3N(9,alpha)

C     redundent information might help to identify 
C     channel relations 

       alphaNN=qnalpha3N(10,alpha)
       pari=qnalpha3N(11,alpha)

      END SUBROUTINE get3Nqn

      SUBROUTINE get3NLSqn(alpha,l12,l3,bl,s12,bs,j3,
     $                           t12,tau3,mtau3,pari)
       IMPLICIT NONE 
       INTEGER alpha,l12,s12,t12,l3,j3,tau3,mtau3,pari,bl,bs

#ifdef DEBUG
C       check for consistency only then debugging 
       IF(alpha.LE.0) STOP 'alpha too small'
       IF(alpha.GT.alpha3NLScdepmax) STOP 'alpha too large'
#endif 

       l12=qnalpha3NLS(1,alpha)
       l3=qnalpha3NLS(2,alpha)
       bl=qnalpha3NLS(3,alpha)
       s12=qnalpha3NLS(4,alpha)
       bs=qnalpha3NLS(5,alpha)
       j3=qnalpha3NLS(6,alpha)
       t12=qnalpha3NLS(7,alpha)
       tau3=qnalpha3NLS(8,alpha)
       mtau3=qnalpha3NLS(9,alpha)

C     redundent information might help to identify 
C     channel relations 
       pari=qnalpha3NLS(10,alpha)

      END SUBROUTINE get3NLSqn

C> This subroutine returns the quantum numbers of 3+1 coordinates 
C> @param[in] alpha index of the channel 1, ..., \ref alpha4N31cdepmax
C> @param[out] l12 l12 
C> @param[out] s12 s12
C> @param[out] j12 j12 
C> @param[out] t12 t12
C> @param[out] l3 l3 
C> @param[out] I3 twice I3
C> @param[out] j3 twice j3
C> @param[out] tau3 twice tau3
C> @param[out] l4 l4 
C> @param[out] I4 twice I4 
C> @param[out] j4 total angular momentum  
C> @param[out] tau4 total isospin 
C> @param[out] mtau4 third component of total isospin
C> @param[out] alphaNN12 index of NN channel 
C> @param[out] alpha3N index of 3N channel 
C> @param[out] pari parity      
      SUBROUTINE get4N31qn(alpha,l12,s12,j12,t12,l3,I3,j3,
     $                 tau3,l4,I4,j4,tau4,mtau4,alphaNN12,alpha3N,pari)
       IMPLICIT NONE 
       INTEGER alpha,l12,s12,j12,t12,l3,I3,j3,tau3,
     $         l4,I4,j4,tau4,mtau4,alphaNN12,alpha3N,pari

#ifdef DEBUG
C       check for consistency only then debugging 
       IF(alpha.LE.0) STOP 'alpha too small'
       IF(alpha.GT.alpha4N31cdepmax) STOP 'alpha too large'
#endif 

       l12       =qnalpha4N31(1,alpha) 
       s12       =qnalpha4N31(2,alpha) 
       j12       =qnalpha4N31(3,alpha) 
       t12       =qnalpha4N31(4,alpha) 
       l3        =qnalpha4N31(5,alpha) 
       I3        =qnalpha4N31(6,alpha) 
       j3        =qnalpha4N31(7,alpha) 
       tau3      =qnalpha4N31(8,alpha) 
       l4        =qnalpha4N31(9,alpha) 
       I4        =qnalpha4N31(10,alpha)
       j4        =qnalpha4N31(11,alpha)
       tau4      =qnalpha4N31(12,alpha)
       mtau4     =qnalpha4N31(13,alpha)
       alphaNN12 =qnalpha4N31(14,alpha)
       alpha3N   =qnalpha4N31(15,alpha)
       pari      =qnalpha4N31(16,alpha)

      END SUBROUTINE get4N31qn

C> This subroutine returns the quantum numbers of 2+2 coordinates 
C> @param[in] beta index of the channel 1, ..., \ref beta4N22cdepmax
C> @param[out] l12 l12 
C> @param[out] s12 s12
C> @param[out] j12 j12 
C> @param[out] t12 t12
C> @param[out] l34 l34
C> @param[out] s34 s34 
C> @param[out] j34 j34 
C> @param[out] t34 t34 
C> @param[out] lam lambda 
C> @param[out] I I 
C> @param[out] j4 total angular momentum  
C> @param[out] tau4 total isospin 
C> @param[out] mtau4 third component of total isospin
C> @param[out] alphaNN12 index of NN channel 
C> @param[out] alphaNN34 index of NN channel 
C> @param[out] pari parity 
      SUBROUTINE get4N22qn(beta,l12,s12,j12,t12,l34,s34,j34,t34,
     $                 lam,I,j4,tau4,mtau4,alphaNN12,alphaNN34,pari)
       IMPLICIT NONE 
       INTEGER beta,l12,s12,j12,t12,l34,s34,j34,t34,
     $         lam,I,j4,tau4,mtau4,alphaNN12,alphaNN34,pari

#ifdef DEBUG
C       check for consistency only then debugging 
       IF(beta.LE.0) STOP 'beta too small'
       IF(beta.GT.beta4N22cdepmax) STOP 'beta too large'
#endif 

       l12      =qnbeta4N22(1,beta) 
       s12      =qnbeta4N22(2,beta) 
       j12      =qnbeta4N22(3,beta) 
       t12      =qnbeta4N22(4,beta) 
       l34      =qnbeta4N22(5,beta) 
       s34      =qnbeta4N22(6,beta) 
       j34      =qnbeta4N22(7,beta) 
       t34      =qnbeta4N22(8,beta) 
       lam      =qnbeta4N22(9,beta) 
       I        =qnbeta4N22(10,beta)
       j4       =qnbeta4N22(11,beta)
       tau4     =qnbeta4N22(12,beta)
       mtau4    =qnbeta4N22(13,beta)
       alphaNN12=qnbeta4N22(14,beta)
       alphaNN34=qnbeta4N22(15,beta)
       pari     =qnbeta4N22(16,beta)

      END SUBROUTINE get4N22qn

C subroutine preparing the pw channels 

      SUBROUTINE preppwbook
       IMPLICIT NONE

       pwbooktime=pwbooktime-MPI_WTIME()

C     opens the parameter file and reads first line with comments 

       OPEN(UNIT=75,FILE='book-para.dat',FORM='FORMATTED',STATUS='OLD')
       READ(75,*,END=100,ERR=100) 

C     prepare the NN bookkeeping 
       CALL prepalphaNN

C     prepare the 3N bookkeeping (jj and LS coupling)

       CALL prepalpha3N
       CALL prepalpha3NLS

C     prepare the 4N bookkeeping 

       CALL prepalpha4N

C     close the parameter file again 

       CLOSE(75)

C     timing of subroutine as an example 

       pwbooktime=pwbooktime+MPI_WTIME()
       
       RETURN 
C     regular end of subroutine 

 100   CONTINUE 
C     treatment of I/O errors 
       STOP 'problem with book-para.dat'
      END SUBROUTINE preppwbook

C     subroutine that prepares the NN channels according to parameters 
C     in UNIT=75 
C     UNIT should be opened at positioned correctly
C     when the subroutine is called
 
      SUBROUTINE prepalphaNN
       IMPLICIT NONE 
       INTEGER l12,s12,j12,t12,mt12,alpha,mt12cdep
       
       IF(allocated(qnalphaNN)) DEALLOCATE(qnalphaNN)

C     first read in parameters relevant for the NN bookkeeping 
       READ(75,*,END=100,ERR=100)   ! first line is comment
       READ(75,*,END=100,ERR=100) l12max,j12max,evenNN,cdepNN,
     $                            mt12min,mt12max      
C     parameters are printed in printpwbook
C     no printing is necessary 

C     count channels 

       alpha=0

C     calculate the mt12 that covers all relevant t12 
      IF(mt12min*mt12max>0) THEN   ! includes zero ? 
       mt12cdep=min(abs(mt12min),abs(mt12max))*(mt12max/abs(mt12max))
      ELSE
       mt12cdep=0
      END IF
      
C     count mt12=mt12cdep channels 
       mt12=mt12cdep
       DO j12=0,j12max
        DO s12=0,1

C     for given j12,s12 define channels
C      first l12=j12  for s12=0,1 (except for j12==0)
C      second l12=|j12-s12|,j12+s12 for s12=1  
C      in this way the singlet-triplet (coupled/uncoupled) is next to each other 
C      and the coupled triplet 
         
         IF(j12.NE.0.OR.s12.EQ.0) THEN
          l12=j12
          t12=(1+(-1)**(l12+s12))/2
          IF(abs(mt12).LE.t12) THEN
           IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $        .AND.l12.LE.l12max) THEN
            alpha=alpha+1
           END IF
          END IF
#ifdef CHUNPHYS          
          t12=(1+(-1)**(l12+s12+1))/2
          IF(abs(mt12).LE.t12) THEN
           IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $        .AND.l12.LE.l12max) THEN
            alpha=alpha+1
           END IF
          END IF
#endif

         END IF

         IF(s12.EQ.1) THEN
          DO l12=abs(j12-s12),j12+s12,2
           t12=(1+(-1)**(l12+s12))/2
           IF(abs(mt12).LE.t12) THEN
            IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $         .AND.l12.LE.l12max) THEN
             alpha=alpha+1
            END IF
           END IF
#ifdef CHUNPHYS          
           t12=(1+(-1)**(l12+s12+1))/2
           IF(abs(mt12).LE.t12) THEN
            IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $         .AND.l12.LE.l12max) THEN
             alpha=alpha+1
            END IF
           END IF
#endif
          END DO
         END IF

        END DO                  ! s12
       END DO                   ! j12 
          
       alphaNNmax=alpha

       IF(cdepNN) THEN
C     count mt12=mt12min,mt12max channels except mt12=0 
        DO mt12=mt12min,mt12max
         IF(mt12.EQ.mt12cdep) cycle     ! mt12=0 is done above
         DO j12=0,j12max
          DO s12=0,1
           
C     for given j12,s12 define channels
C     first l12=j12  for s12=0,1 (except for j12==0)
C     second l12=|j12-s12|,j12+s12 for s12=1  
C     in this way the singlet-triplet (coupled/uncoupled) is next to each other 
C     and the coupled triplet 

           IF(j12.NE.0.OR.s12.EQ.0) THEN
            l12=j12
            t12=(1+(-1)**(l12+s12))/2
            IF(abs(mt12).LE.t12) THEN
             IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN
              alpha=alpha+1
             END IF
            END IF
#ifdef CHUNPHYS          
            t12=(1+(-1)**(l12+s12+1))/2
            IF(abs(mt12).LE.t12) THEN
             IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN
              alpha=alpha+1
             END IF
            END IF
#endif
           END IF
           
           IF(s12.EQ.1) THEN
            DO l12=abs(j12-s12),j12+s12,2
             t12=(1+(-1)**(l12+s12))/2
             IF(abs(mt12).LE.t12) THEN
              IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $            .AND.l12.LE.l12max) THEN               
               alpha=alpha+1
              END IF
             END IF
#ifdef CHUNPHYS          
             t12=(1+(-1)**(l12+s12+1))/2
             IF(abs(mt12).LE.t12) THEN
              IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN               
               alpha=alpha+1
              END IF
             END IF
#endif
            END DO
           END IF
           
          END DO                ! s12
         END DO                 ! j12 
        END DO                  ! mt12 

       END IF                   ! cdepNN

       alphaNNcdepmax=alpha

C     now allocate memory for bookkeeping and redo the 
C     counting with storing 
       
       ALLOCATE(qnalphaNN(5,alphaNNcdepmax))

       alpha=0
C     count mt12=0 channels 
       mt12=mt12cdep
       DO j12=0,j12max
        DO s12=0,1

C     for given j12,s12 define channels
C      first l12=j12  for s12=0,1 (except for j12==0)
C      second l12=|j12-s12|,j12+s12 for s12=1  
C      in this way the singlet-triplet (coupled/uncoupled) is next to each other 
C      and the coupled triplet 
         
         IF(j12.NE.0.OR.s12.EQ.0) THEN
          l12=j12
          t12=(1+(-1)**(l12+s12))/2
          IF(abs(mt12).LE.t12) THEN
           IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN               
            alpha=alpha+1
            qnalphaNN(1,alpha)=l12
            qnalphaNN(2,alpha)=s12
            qnalphaNN(3,alpha)=j12
            qnalphaNN(4,alpha)=t12
            qnalphaNN(5,alpha)=mt12
           END IF
          END IF
#ifdef CHUNPHYS          
          t12=(1+(-1)**(l12+s12+1))/2
          IF(abs(mt12).LE.t12) THEN
           IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN               
            alpha=alpha+1
            qnalphaNN(1,alpha)=l12
            qnalphaNN(2,alpha)=s12
            qnalphaNN(3,alpha)=j12
            qnalphaNN(4,alpha)=t12
            qnalphaNN(5,alpha)=mt12
           END IF
          END IF
#endif
         END IF

         IF(s12.EQ.1) THEN
          DO l12=abs(j12-s12),j12+s12,2
           t12=(1+(-1)**(l12+s12))/2
           IF(abs(mt12).LE.t12) THEN
            IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN               
             alpha=alpha+1
             qnalphaNN(1,alpha)=l12
             qnalphaNN(2,alpha)=s12
             qnalphaNN(3,alpha)=j12
             qnalphaNN(4,alpha)=t12
             qnalphaNN(5,alpha)=mt12
            END IF
           END IF
#ifdef CHUNPHYS          
           t12=(1+(-1)**(l12+s12+1))/2
           IF(abs(mt12).LE.t12) THEN
            IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN               
             alpha=alpha+1
             qnalphaNN(1,alpha)=l12
             qnalphaNN(2,alpha)=s12
             qnalphaNN(3,alpha)=j12
             qnalphaNN(4,alpha)=t12
             qnalphaNN(5,alpha)=mt12
            END IF
           END IF
#endif
          END DO
         END IF

        END DO                  ! s12
       END DO                   ! j12 

#ifdef DEBUG 
C     check consistency of both runs 
       IF(alphaNNmax.NE.alpha) STOP 'inconsistent alphaNN'
#endif 

       IF(cdepNN) THEN
C     count mt12=mt12min,mt12max channels 
        DO mt12=mt12min,mt12max
         IF(mt12.EQ.mt12cdep) cycle
         DO j12=0,j12max
          DO s12=0,1
           
C     for given j12,s12 define channels
C     first l12=j12  for s12=0,1 (except for j12==0)
C     second l12=|j12-s12|,j12+s12 for s12=1  
C     in this way the singlet-triplet (coupled/uncoupled) is next to each other 
C     and the coupled triplet 
           
           IF(j12.NE.0.OR.s12.EQ.0) THEN
            l12=j12
            t12=(1+(-1)**(l12+s12))/2
            IF(abs(mt12).LE.t12) THEN
             IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN               
              alpha=alpha+1
              qnalphaNN(1,alpha)=l12
              qnalphaNN(2,alpha)=s12
              qnalphaNN(3,alpha)=j12
              qnalphaNN(4,alpha)=t12
              qnalphaNN(5,alpha)=mt12
             END IF
            END IF
#ifdef CHUNPHYS          
            t12=(1+(-1)**(l12+s12+1))/2
            IF(abs(mt12).LE.t12) THEN
             IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN               
              alpha=alpha+1
              qnalphaNN(1,alpha)=l12
              qnalphaNN(2,alpha)=s12
              qnalphaNN(3,alpha)=j12
              qnalphaNN(4,alpha)=t12
              qnalphaNN(5,alpha)=mt12
             END IF
            END IF
#endif
           END IF
           
           IF(s12.EQ.1) THEN
            DO l12=abs(j12-s12),j12+s12,2
             t12=(1+(-1)**(l12+s12))/2
             IF(abs(mt12).LE.t12) THEN
              IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN               
               alpha=alpha+1
               qnalphaNN(1,alpha)=l12
               qnalphaNN(2,alpha)=s12
               qnalphaNN(3,alpha)=j12
               qnalphaNN(4,alpha)=t12
               qnalphaNN(5,alpha)=mt12
              END IF
             END IF
#ifdef CHUNPHYS          
             t12=(1+(-1)**(l12+s12+1))/2
             IF(abs(mt12).LE.t12) THEN
              IF((.NOT. evenNN .OR. mod(l12,2).EQ.0)
     $           .AND.l12.LE.l12max) THEN               
               alpha=alpha+1
               qnalphaNN(1,alpha)=l12
               qnalphaNN(2,alpha)=s12
               qnalphaNN(3,alpha)=j12
               qnalphaNN(4,alpha)=t12
               qnalphaNN(5,alpha)=mt12
              END IF
             END IF
#endif
            END DO
           END IF
           
          END DO                ! s12
         END DO                 ! j12 
        END DO                  ! mt12 

       END IF                   ! cdepNN

#ifdef DEBUG
       IF(alpha.NE.alphaNNcdepmax) STOP 'inconsistent alphaNN'
#endif 

       RETURN
C     regular end of subroutine 
 100   CONTINUE 
C     treatment of I/O errors 
       STOP 'problem with book-para.dat'

      END SUBROUTINE prepalphaNN

C     subroutine that prepares the 3N channels according to parameters 
C     in UNIT=75 
C     UNIT should be opened at positioned correctly
C     when the subroutine is called
 
      SUBROUTINE prepalpha3N
       IMPLICIT NONE 
       INTEGER l12,s12,j12,t12,mt12,alpha,mtau3cdep
       INTEGER l3,I3,j3,tau3,mtau3,pari,alphaNN

       IF(allocated(qnalpha3N)) DEALLOCATE(qnalpha3N)

C     first read in parameters relevant for the 3N bookkeeping 
       READ(75,*,END=100,ERR=100)   ! first line is comment
       READ(75,*,END=100,ERR=100) l3max,I3max,j3min,j3max,
     $                            pari3Nmin,pari3Nmax     
       READ(75,*,END=100,ERR=100) tau3min,tau3max,cdep3N,
     $                            mtau3min,mtau3max
C     parameters are printed in printpwbook
C     no printing is necessary 

C     count channels 

       alpha=0

C     calculate the mtau3 that covers all relevant tau3 
      IF(mtau3min*mtau3max>0) THEN   ! includes zero ? 
       mtau3cdep=min(abs(mtau3min),abs(mtau3max))*mtau3max/abs(mtau3max)
      ELSE
       mtau3cdep=1
      END IF
 
C     count mtau3=mtau3cdep channels first 
       mtau3=mtau3cdep
       DO tau3=max(tau3min,abs(mtau3)),tau3max,2 
        DO j3=j3min,j3max,2
         DO pari=pari3Nmin,pari3Nmax,2          

C     for given total quantum numbers run spectator quantum numbers 
          
          DO I3=1,I3max,2 
           DO l3=abs(I3-1)/2,min((I3+1)/2,l3max)
            
C go through all NN channels and check whether they fit 
            DO alphaNN=1,alphaNNmax
             CALL getNNqn(alphaNN,l12,s12,j12,t12,mt12)
             
             IF(     pari.EQ.(-1)**(l12+l3)
     $        .AND.  j12.GE.abs(I3-j3)/2
     $        .AND.  j12.LE.(I3+j3)/2
     $        .AND.  t12.GE.abs(tau3-1)/2
     $        .AND.  t12.LE.(tau3+1)/2     ) THEN 
              
              alpha=alpha+1
             END IF
            END DO              ! alphaNN

           END DO               ! l3 
          END DO                ! I3
          
         END DO                 ! pari
        END DO                  ! j3 
       END DO                   ! tau3 
         
          
       alpha3Nmax=alpha

       IF(cdep3N) THEN
C     count mtau3=mtau3min,mtau3max,2 channels 
        DO mtau3=mtau3min,mtau3max,2 
         IF(mtau3.EQ.mtau3cdep) cycle    ! skip the mtau3=mtau3cdep part since it is done           
         DO tau3=max(tau3min,abs(mtau3)),tau3max,2 
          DO j3=j3min,j3max,2
           DO pari=pari3Nmin,pari3Nmax,2          
            
C     for given total quantum numbers run spectator quantum numbers 
          
            DO I3=1,I3max,2 
             DO l3=abs(I3-1)/2,min((I3+1)/2,l3max)
            
C go through all NN channels and check whether they fit 
              DO alphaNN=1,alphaNNmax
               CALL getNNqn(alphaNN,l12,s12,j12,t12,mt12)
             
               IF(     pari.EQ.(-1)**(l12+l3)
     $          .AND.  j12.GE.abs(I3-j3)/2
     $          .AND.  j12.LE.(I3+j3)/2
     $          .AND.  t12.GE.abs(tau3-1)/2
     $          .AND.  t12.LE.(tau3+1)/2     ) THEN 
                
                alpha=alpha+1
               END IF
              END DO            ! alphaNN

             END DO             ! l3 
            END DO              ! I3
          
           END DO               ! pari
          END DO                ! j3 
         END DO                 ! tau3 
        END DO                  ! mtau3

       END IF                   ! cdep3N

       alpha3Ncdepmax=alpha

C     now allocate memory for bookkeeping and redo the 
C     counting with storing 
       
       ALLOCATE(qnalpha3N(11,alpha3Ncdepmax))

C     count channels 

       alpha=0

C     count mtau3=mtau3cdep channels first 
       mtau3=mtau3cdep
       DO tau3=max(tau3min,abs(mtau3)),tau3max,2 
        DO j3=j3min,j3max,2
         DO pari=pari3Nmin,pari3Nmax,2          

C     for given total quantum numbers run spectator quantum numbers 
          
          DO I3=1,I3max,2 
           DO l3=abs(I3-1)/2,min((I3+1)/2,l3max)
            
C go through all NN channels and check whether they fit 
            DO alphaNN=1,alphaNNmax
             CALL getNNqn(alphaNN,l12,s12,j12,t12,mt12)
             
             IF(     pari.EQ.(-1)**(l12+l3)
     $        .AND.  j12.GE.abs(I3-j3)/2
     $        .AND.  j12.LE.(I3+j3)/2
     $        .AND.  t12.GE.abs(tau3-1)/2
     $        .AND.  t12.LE.(tau3+1)/2     ) THEN 
              
              alpha=alpha+1
              qnalpha3N(1,alpha)=l12
              qnalpha3N(2,alpha)=s12
              qnalpha3N(3,alpha)=j12
              qnalpha3N(4,alpha)=t12
              qnalpha3N(5,alpha)=l3
              qnalpha3N(6,alpha)=I3
              qnalpha3N(7,alpha)=j3
              qnalpha3N(8,alpha)=tau3
              qnalpha3N(9,alpha)=mtau3
              qnalpha3N(10,alpha)=alphaNN
              qnalpha3N(11,alpha)=pari
             END IF
            END DO              ! alphaNN

           END DO               ! l3 
          END DO                ! I3
          
         END DO                 ! pari
        END DO                  ! j3 
       END DO                   ! tau3 
         
          
#ifdef DEBUG 
C     check consistency of both runs 
       IF(alpha3Nmax.NE.alpha) STOP 'inconsistent alpha3N'
#endif 

       IF(cdep3N) THEN
C     count mtau3=mtau3min,mtau3max,2 channels 
        DO mtau3=mtau3min,mtau3max,2 
         IF(mtau3.EQ.mtau3cdep) cycle    ! skip the mtau3=1 part since it is done  
         DO tau3=max(tau3min,abs(mtau3)),tau3max,2 
          DO j3=j3min,j3max,2
           DO pari=pari3Nmin,pari3Nmax,2          
            
C     for given total quantum numbers run spectator quantum numbers 
          
            DO I3=1,I3max,2 
             DO l3=abs(I3-1)/2,min((I3+1)/2,l3max)
            
C go through all NN channels and check whether they fit 
              DO alphaNN=1,alphaNNmax
               CALL getNNqn(alphaNN,l12,s12,j12,t12,mt12)
             
               IF(     pari.EQ.(-1)**(l12+l3)
     $          .AND.  j12.GE.abs(I3-j3)/2
     $          .AND.  j12.LE.(I3+j3)/2
     $          .AND.  t12.GE.abs(tau3-1)/2
     $          .AND.  t12.LE.(tau3+1)/2     ) THEN 
                
                alpha=alpha+1
                qnalpha3N(1,alpha)=l12
                qnalpha3N(2,alpha)=s12
                qnalpha3N(3,alpha)=j12
                qnalpha3N(4,alpha)=t12
                qnalpha3N(5,alpha)=l3
                qnalpha3N(6,alpha)=I3
                qnalpha3N(7,alpha)=j3
                qnalpha3N(8,alpha)=tau3
                qnalpha3N(9,alpha)=mtau3
                qnalpha3N(10,alpha)=alphaNN
                qnalpha3N(11,alpha)=pari
               END IF
              END DO            ! alphaNN

             END DO             ! l3 
            END DO              ! I3
          
           END DO               ! pari
          END DO                ! j3 
         END DO                 ! tau3 
        END DO                  ! mtau3

       END IF                   ! cdep3N

#ifdef DEBUG
       IF(alpha.NE.alpha3Ncdepmax) STOP 'inconsistent alpha3N'
#endif 

       RETURN
C     regular end of subroutine 
 100   CONTINUE 
C     treatment of I/O errors 
       STOP 'problem with book-para.dat'

      END SUBROUTINE prepalpha3N

      
      SUBROUTINE prepalpha3NLS
       IMPLICIT NONE 
       INTEGER l12,s12,j12,t12,mt12,alpha
       INTEGER l3,bl,bs,j3,tau3,mtau3,pari,mtau3cdep

       IF(allocated(qnalpha3NLS)) DEALLOCATE(qnalpha3NLS)

!     parmeters are read in in prepalpha3NLS
!     no reading is necessary here 
!     parameters are printed in printpwbook
!     no printing is necessary 

!     calculate the mtau3 that covers all relevant tau3 
      IF(mtau3min*mtau3max>0) THEN   ! includes zero ? 
       mtau3cdep=min(abs(mtau3min),abs(mtau3max))*mtau3max/abs(mtau3max)
      ELSE
       mtau3cdep=1
      END IF
      
!     count channels 

       alpha=0
       IF(cdep3N) THEN
        DO mtau3=mtau3min,mtau3max,2 
         DO tau3=max(tau3min,abs(mtau3)),tau3max,2 
          DO j3=j3min,j3max,2
           DO pari=pari3Nmin,pari3Nmax,2          

 !     for given total quantum numbers run spin quantum numbers 
            DO s12=0,1
             DO bs=abs(2*s12-1),2*s12+1,2            
 !     now run orbital quantum numbers             
              DO l12=0,l12max 
               DO l3=0,l3max
                IF(pari.NE.(-1)**(l12+l3)) cycle ! check parity 
                DO bl=max(abs(l12-l3),abs(j3-bs)/2),
     $               min(l12+l3,(j3+bs)/2) 
 !      change isospin 
                 DO t12=0,1 
 !     check whether the channels are physical 
                  IF(t12.LT.abs(tau3-1)/2) cycle 
                  IF(t12.GT.(tau3+1)/2) cycle  
#ifndef CHUNPHYS 
                  IF((-1)**(l12+s12+t12).GT.0) cycle
#endif                
                  alpha=alpha+1
                 END DO    ! t12
                END DO     ! bl 
               END DO      ! l3 
              END DO       ! l12
             END DO        ! bs 
            END DO         ! s12          
           END DO          ! pari
          END DO           ! j3 
         END DO            ! tau3
        END DO             ! mtau3
       ELSE ! cdep3N 
        mtau3=mtau3cdep 
        DO tau3=max(tau3min,abs(mtau3)),tau3max,2 
         DO j3=j3min,j3max,2
          DO pari=pari3Nmin,pari3Nmax,2          

!     for given total quantum numbers run spin quantum numbers 
           DO s12=0,1
            DO bs=abs(2*s12-1),2*s12+1,2            
!     now run orbital quantum numbers             
             DO l12=0,l12max 
              DO l3=0,l3max
               IF(pari.NE.(-1)**(l12+l3)) cycle ! check parity 
               DO bl=max(abs(l12-l3),abs(j3-bs)/2),
     $               min(l12+l3,(j3+bs)/2) 
!      change isospin 
                DO t12=0,1 
!     check whether the channels are physical 
                 IF(t12.LT.abs(tau3-1)/2) cycle 
                 IF(t12.GT.(tau3+1)/2) cycle  
#ifndef CHUNPHYS 
                 IF((-1)**(l12+s12+t12).GT.0) cycle
#endif                
                 alpha=alpha+1
                END DO    ! t12
               END DO     ! bl 
              END DO      ! l3 
             END DO       ! l12
            END DO        ! bs 
           END DO         ! s12          
          END DO          ! pari
         END DO           ! j3 
        END DO            ! tau3
       END IF             ! cdep3N
         
       alpha3NLScdepmax=alpha

!     now allocate memory for bookkeeping and redo the 
!     counting with storing 
       
       ALLOCATE(qnalpha3NLS(10,alpha3NLScdepmax))
       
!     count channels up

       alpha=0

       IF(cdep3N) THEN
        DO mtau3=mtau3min,mtau3max,2 
         DO tau3=max(tau3min,abs(mtau3)),tau3max,2 
          DO j3=j3min,j3max,2
           DO pari=pari3Nmin,pari3Nmax,2          

 !     for given total quantum numbers run spin quantum numbers 
            DO s12=0,1
             DO bs=abs(2*s12-1),2*s12+1,2            
 !     now run orbital quantum numbers             
              DO l12=0,l12max 
               DO l3=0,l3max
                IF(pari.NE.(-1)**(l12+l3)) cycle ! check parity 
                DO bl=max(abs(l12-l3),abs(j3-bs)/2),
     $               min(l12+l3,(j3+bs)/2)
 !      change isospin 
                 DO t12=0,1 
 !     check whether the channels are physical 
                  IF(t12.LT.abs(tau3-1)/2) cycle 
                  IF(t12.GT.(tau3+1)/2) cycle  
#ifndef CHUNPHYS 
                  IF((-1)**(l12+s12+t12).GT.0) cycle
#endif                
                  alpha=alpha+1
                  qnalpha3NLS(1,alpha)=l12
                  qnalpha3NLS(2,alpha)=l3
                  qnalpha3NLS(3,alpha)=bl
                  qnalpha3NLS(4,alpha)=s12
                  qnalpha3NLS(5,alpha)=bs
                  qnalpha3NLS(6,alpha)=j3
                  qnalpha3NLS(7,alpha)=t12
                  qnalpha3NLS(8,alpha)=tau3
                  qnalpha3NLS(9,alpha)=mtau3
                  qnalpha3NLS(10,alpha)=pari
                 END DO    ! t12
                END DO     ! bl 
               END DO      ! l3 
              END DO       ! l12
             END DO        ! bs 
            END DO         ! s12          
           END DO          ! pari
          END DO           ! j3 
         END DO            ! tau3
        END DO             ! mtau3
       ELSE ! cdep3N
        mtau3=mtau3cdep
        DO tau3=max(tau3min,abs(mtau3)),tau3max,2 
         DO j3=j3min,j3max,2
          DO pari=pari3Nmin,pari3Nmax,2          

!     for given total quantum numbers run spin quantum numbers 
           DO s12=0,1
            DO bs=abs(2*s12-1),2*s12+1,2            
!     now run orbital quantum numbers             
             DO l12=0,l12max 
              DO l3=0,l3max
               IF(pari.NE.(-1)**(l12+l3)) cycle ! check parity 
               DO bl=max(abs(l12-l3),abs(j3-bs)/2),
     $               min(l12+l3,(j3+bs)/2)
!      change isospin 
                DO t12=0,1 
!     check whether the channels are physical 
                 IF(t12.LT.abs(tau3-1)/2) cycle 
                 IF(t12.GT.(tau3+1)/2) cycle  
#ifndef CHUNPHYS 
                 IF((-1)**(l12+s12+t12).GT.0) cycle
#endif                
                 alpha=alpha+1
                 qnalpha3NLS(1,alpha)=l12
                 qnalpha3NLS(2,alpha)=l3
                 qnalpha3NLS(3,alpha)=bl
                 qnalpha3NLS(4,alpha)=s12
                 qnalpha3NLS(5,alpha)=bs
                 qnalpha3NLS(6,alpha)=j3
                 qnalpha3NLS(7,alpha)=t12
                 qnalpha3NLS(8,alpha)=tau3
                 qnalpha3NLS(9,alpha)=mtau3
                 qnalpha3NLS(10,alpha)=pari
                END DO    ! t12
               END DO     ! bl 
              END DO      ! l3 
             END DO       ! l12
            END DO        ! bs 
           END DO         ! s12          
          END DO          ! pari
         END DO           ! j3 
        END DO            ! tau3
       END IF       ! cdep3N
#ifdef DEBUG
       IF(alpha.NE.alpha3NLScdepmax) STOP 'inconsistent alpha3NLS'
#endif 

      END SUBROUTINE prepalpha3NLS
      
C     subroutine that prepares the 4N channels according to parameters 
C     in UNIT=75 
C     UNIT should be opened at positioned correctly
C     when the subroutine is called
 
      SUBROUTINE prepalpha4N
       IMPLICIT NONE 
       INTEGER l12,s12,j12,t12,mt12,alpha,mtau4cdep
       INTEGER l3,I3,j3,tau3,mtau3,pari3N,alphaNN12
       INTEGER l4,I4,j4,tau4,mtau4,pari,alpha3N
       INTEGER lam,I,beta,l34,s34,j34,t34,mt34,alphaNN34

       IF(allocated(qnbeta4N22)) DEALLOCATE(qnbeta4N22)
       IF(allocated(qnalpha4N31)) DEALLOCATE(qnalpha4N31)

C     first read in parameters relevant for the 4N bookkeeping 
       READ(75,*,END=100,ERR=100)   ! first line is comment
       READ(75,*,END=100,ERR=100) l4max,I4max,j4min,j4max,
     $                            pari4Nmin,pari4Nmax     
       READ(75,*,END=100,ERR=100) lammax,Imax     
       READ(75,*,END=100,ERR=100) tau4min,tau4max,cdep4N,
     $                            mtau4min,mtau4max
       READ(75,*,END=100,ERR=100) lsummax31,lsummax22
C     parameters are printed in printpwbook
C     no printing is necessary 


C  3+1 bookkeeping 

C     count channels 

       alpha=0

C     calculate the mtau4 that covers all relevant tau4 
      IF(mtau4min*mtau4max>0) THEN   ! includes zero ? 
       mtau4cdep=min(abs(mtau4min),abs(mtau4max))*mtau4max/abs(mtau4max)
      ELSE
       mtau4cdep=0
      END IF
 
C     count mtau4=mtau4cdep channels first 
       mtau4=mtau4cdep
       DO tau4=max(tau4min,abs(mtau4)),tau4max
        DO j4=j4min,j4max
         DO pari=pari4Nmin,pari4Nmax,2          

C     for given total quantum numbers run spectator quantum numbers 
          
          DO I4=1,I4max,2 
           DO l4=abs(I4-1)/2,min((I4+1)/2,l4max)
            
C go through all 3N channels and check whether they fit 
            DO alpha3N=1,alpha3Nmax
             call get3Nqn(alpha3N,l12,s12,j12,t12,l3,I3,j3,
     $                           tau3,mtau3,alphaNN12,pari3N)
            IF(l12+l3+l4 .GT.lsummax31) CYCLE

             IF(     pari.EQ.(-1)**(l12+l3+l4)
     $        .AND.  j3.GE.abs(I4-2*j4)
     $        .AND.  j3.LE.(I4+2*j4)
     $        .AND.  tau3.GE.abs(2*tau4-1)
     $        .AND.  tau3.LE.(2*tau4+1)     ) THEN 
              
              alpha=alpha+1
             END IF
            END DO              ! alpha3N

           END DO               ! l4 
          END DO                ! I4
          
         END DO                 ! pari
        END DO                  ! j4 
       END DO                   ! tau4 
         
          
       alpha4N31max=alpha

       IF(cdep4N) THEN
C     count mtau4=mtau4min,mtau4max channels 
        DO mtau4=mtau4min,mtau4max 
         IF(mtau4.EQ.mtau4cdep) cycle    ! skip the mtau4=mtau4cdep part since it is done           
         DO tau4=max(tau4min,abs(mtau4)),tau4max 
          DO j4=j4min,j4max
           DO pari=pari4Nmin,pari4Nmax,2          
            
C     for given total quantum numbers run spectator quantum numbers 
          
            DO I4=1,I4max,2 
             DO l4=abs(I4-1)/2,min((I4+1)/2,l4max)
            
C go through all 3N channels and check whether they fit 
              DO alpha3N=1,alpha3Nmax
               call get3Nqn(alpha3N,l12,s12,j12,t12,l3,I3,j3,
     $                           tau3,mtau3,alphaNN12,pari3N)
               IF(l12+l3+l4 .GT.lsummax31) CYCLE
               
               IF(     pari.EQ.(-1)**(l12+l3+l4)
     $          .AND.  j3.GE.abs(I4-2*j4)
     $          .AND.  j3.LE.(I4+2*j4)
     $          .AND.  tau3.GE.abs(2*tau4-1)
     $          .AND.  tau3.LE.(2*tau4+1)     ) THEN 
              
                alpha=alpha+1
               END IF
              END DO            ! alpha3N

             END DO             ! l4
            END DO              ! I4
          
           END DO               ! pari
          END DO                ! j4 
         END DO                 ! tau4
        END DO                  ! mtau4

       END IF                   ! cdep3N

       alpha4N31cdepmax=alpha

C     now allocate memory for bookkeeping and redo the 
C     counting with storing 
       
       ALLOCATE(qnalpha4N31(16,alpha4N31cdepmax))

C     start with  channel 0

       alpha=0

C      mtau4=mtau4cdep channels first 
       mtau4=mtau4cdep
       DO tau4=max(tau4min,abs(mtau4)),tau4max
        DO j4=j4min,j4max
         DO pari=pari4Nmin,pari4Nmax,2          

C     for given total quantum numbers run spectator quantum numbers 
          
          DO I4=1,I4max,2 
           DO l4=abs(I4-1)/2,min((I4+1)/2,l4max)
            
C go through all 3N channels and check whether they fit 
            DO alpha3N=1,alpha3Nmax
             call get3Nqn(alpha3N,l12,s12,j12,t12,l3,I3,j3,
     $                           tau3,mtau3,alphaNN12,pari3N)
             IF(l12+l3+l4 .GT.lsummax31) CYCLE

             IF(     pari.EQ.(-1)**(l12+l3+l4)
     $        .AND.  j3.GE.abs(I4-2*j4)
     $        .AND.  j3.LE.(I4+2*j4)
     $        .AND.  tau3.GE.abs(2*tau4-1)
     $        .AND.  tau3.LE.(2*tau4+1)     ) THEN 
              
              alpha=alpha+1

              qnalpha4N31(1,alpha)=l12
              qnalpha4N31(2,alpha)=s12
              qnalpha4N31(3,alpha)=j12
              qnalpha4N31(4,alpha)=t12
              qnalpha4N31(5,alpha)=l3
              qnalpha4N31(6,alpha)=I3
              qnalpha4N31(7,alpha)=j3
              qnalpha4N31(8,alpha)=tau3
              qnalpha4N31(9,alpha)=l4
              qnalpha4N31(10,alpha)=I4
              qnalpha4N31(11,alpha)=j4
              qnalpha4N31(12,alpha)=tau4
              qnalpha4N31(13,alpha)=mtau4
              qnalpha4N31(14,alpha)=alphaNN12
              qnalpha4N31(15,alpha)=alpha3N
              qnalpha4N31(16,alpha)=pari

             END IF
            END DO              ! alpha3N

           END DO               ! l4 
          END DO                ! I4
          
         END DO                 ! pari
        END DO                  ! j4 
       END DO                   ! tau4 
         
          
#ifdef DEBUG
       IF(alpha.NE.alpha4N31max) STOP 'inconsistent alpha4N'
#endif 

       IF(cdep4N) THEN
C         mtau4=mtau4min,mtau4max channels without mtau=0 
        DO mtau4=mtau4min,mtau4max 
         IF(mtau4.EQ.mtau4cdep) cycle    ! skip the mtau4=mtau4cdep part since it is done           
         DO tau4=max(tau4min,abs(mtau4)),tau4max 
          DO j4=j4min,j4max
           DO pari=pari4Nmin,pari4Nmax,2          
            
C     for given total quantum numbers run spectator quantum numbers 
          
            DO I4=1,I4max,2 
             DO l4=abs(I4-1)/2,min((I4+1)/2,l4max)
            
C go through all 3N channels and check whether they fit 
              DO alpha3N=1,alpha3Nmax
               call get3Nqn(alpha3N,l12,s12,j12,t12,l3,I3,j3,
     $                           tau3,mtau3,alphaNN12,pari3N)
               IF(l12+l3+l4 .GT.lsummax31) CYCLE

               IF(     pari.EQ.(-1)**(l12+l3+l4)
     $          .AND.  j3.GE.abs(I4-2*j4)
     $          .AND.  j3.LE.(I4+2*j4)
     $          .AND.  tau3.GE.abs(2*tau4-1)
     $          .AND.  tau3.LE.(2*tau4+1)     ) THEN 
              
                alpha=alpha+1
                
                qnalpha4N31(1,alpha)=l12
                qnalpha4N31(2,alpha)=s12
                qnalpha4N31(3,alpha)=j12
                qnalpha4N31(4,alpha)=t12
                qnalpha4N31(5,alpha)=l3
                qnalpha4N31(6,alpha)=I3
                qnalpha4N31(7,alpha)=j3
                qnalpha4N31(8,alpha)=tau3
                qnalpha4N31(9,alpha)=l4
                qnalpha4N31(10,alpha)=I4
                qnalpha4N31(11,alpha)=j4
                qnalpha4N31(12,alpha)=tau4
                qnalpha4N31(13,alpha)=mtau4
                qnalpha4N31(14,alpha)=alphaNN12
                qnalpha4N31(15,alpha)=alpha3N
                qnalpha4N31(16,alpha)=pari
                
               END IF
              END DO            ! alpha3N

             END DO             ! l4
            END DO              ! I4
          
           END DO               ! pari
          END DO                ! j4 
         END DO                 ! tau4
        END DO                  ! mtau4

       END IF                   ! cdep4N

#ifdef DEBUG
       IF(alpha.NE.alpha4N31cdepmax) STOP 'inconsistent alpha4N'
#endif 
C     3+1 bookkeeping done 

C     2+2 bookkeeping 
C     count channels 
       beta=0

C     count mtau4=mtau4cdep channels first 
       mtau4=mtau4cdep
       DO tau4=max(tau4min,abs(mtau4)),tau4max
        DO j4=j4min,j4max
         DO pari=pari4Nmin,pari4Nmax,2          

C     for given total quantum numbers run 34 quantum numbers 
          DO alphaNN34=1,alphaNNmax
           CALL getNNqn(alphaNN34,l34,s34,j34,t34,mt34)

C     and 12 quantum numbers 
           DO alphaNN12=1,alphaNNmax
            CALL getNNqn(alphaNN12,l12,s12,j12,t12,mt12)
            
C     run inter cluster quantum numbers 
            DO I=abs(j4-j34),min(Imax,j4+j34) 
             DO lam=abs(I-j12),min(I+j12,lammax)
              IF(lam+l12+l34 .GT.lsummax22) CYCLE
C     check whether channels fit 
              IF(     pari.EQ.(-1)**(l12+l34+lam)
     $         .AND.  tau4.GE.abs(t12-t34)
     $         .AND.  tau4.LE.(t12+t34)     ) THEN 
               
               beta=beta+1
              END IF
             END DO             ! lam
            END DO              ! I
            
           END DO               ! alphaNN12
          END DO                ! alphaNN34
         END DO                 ! pari
        END DO                  ! j4 
       END DO                   ! tau4 
         
          
       beta4N22max=beta

       IF(cdep4N) THEN
C     count mtau4=mtau4min,mtau4max channels 
        DO mtau4=mtau4min,mtau4max 
         IF(mtau4.EQ.mtau4cdep) cycle    ! skip the mtau4=mtau4cdep part since it is done           

         DO tau4=max(tau4min,abs(mtau4)),tau4max
          DO j4=j4min,j4max
           DO pari=pari4Nmin,pari4Nmax,2          

C     for given total quantum numbers run 34 quantum numbers 
            DO alphaNN34=1,alphaNNmax
             CALL getNNqn(alphaNN34,l34,s34,j34,t34,mt34)

C     and 12 quantum numbers 
             DO alphaNN12=1,alphaNNmax
              CALL getNNqn(alphaNN12,l12,s12,j12,t12,mt12)
            
C     run inter cluster quantum numbers 
              DO I=abs(j4-j34),min(Imax,j4+j34) 
               DO lam=abs(I-j12),min(I+j12,lammax)
                IF(lam+l12+l34 .GT.lsummax22) CYCLE
C     check whether channels fit 
                IF(     pari.EQ.(-1)**(l12+l34+lam)
     $           .AND.  tau4.GE.abs(t12-t34)
     $           .AND.  tau4.LE.(t12+t34)     ) THEN 
               
                 beta=beta+1
                END IF
               END DO           ! lam
              END DO            ! I
            
             END DO             ! alphaNN12
            END DO              ! alphaNN34
           END DO               ! pari
          END DO                ! j4 
         END DO                 ! tau4 
         
       
        END DO                  ! mtau4
       END IF                   ! cdep4N

       beta4N22cdepmax=beta

C     now allocate memory for bookkeeping and redo the 
C     counting with storing 
       
       ALLOCATE(qnbeta4N22(16,beta4N22cdepmax))

C     start with  channel 0

       beta=0

C     count mtau4=mtau4cdep channels first 
       mtau4=mtau4cdep
       DO tau4=max(tau4min,abs(mtau4)),tau4max
        DO j4=j4min,j4max
         DO pari=pari4Nmin,pari4Nmax,2          

C     for given total quantum numbers run 34 quantum numbers 
          DO alphaNN34=1,alphaNNmax
           CALL getNNqn(alphaNN34,l34,s34,j34,t34,mt34)

C     and 12 quantum numbers 
           DO alphaNN12=1,alphaNNmax
            CALL getNNqn(alphaNN12,l12,s12,j12,t12,mt12)
            
C     run inter cluster quantum numbers 
            DO I=abs(j4-j34),min(Imax,j4+j34) 
             DO lam=abs(I-j12),min(I+j12,lammax)
              IF(lam+l12+l34 .GT.lsummax22) CYCLE
C     check whether channels fit 
              IF(     pari.EQ.(-1)**(l12+l34+lam)
     $         .AND.  tau4.GE.abs(t12-t34)
     $         .AND.  tau4.LE.(t12+t34)     ) THEN 
               
               beta=beta+1

               qnbeta4N22(1,beta)=l12
               qnbeta4N22(2,beta)=s12
               qnbeta4N22(3,beta)=j12
               qnbeta4N22(4,beta)=t12
               qnbeta4N22(5,beta)=l34
               qnbeta4N22(6,beta)=s34
               qnbeta4N22(7,beta)=j34
               qnbeta4N22(8,beta)=t34
               qnbeta4N22(9,beta)=lam
               qnbeta4N22(10,beta)=I
               qnbeta4N22(11,beta)=j4
               qnbeta4N22(12,beta)=tau4
               qnbeta4N22(13,beta)=mtau4
               qnbeta4N22(14,beta)=alphaNN12
               qnbeta4N22(15,beta)=alphaNN34
               qnbeta4N22(16,beta)=pari

              END IF
             END DO             ! lam
            END DO              ! I
            
           END DO               ! alphaNN12
          END DO                ! alphaNN34
         END DO                 ! pari
        END DO                  ! j4 
       END DO                   ! tau4 
         
          
#ifdef DEBUG
       IF(beta.NE.beta4N22max) STOP 'inconsistent beta4N'
#endif 

       IF(cdep4N) THEN
C     count mtau4=mtau4min,mtau4max channels 
        DO mtau4=mtau4min,mtau4max 
         IF(mtau4.EQ.mtau4cdep) cycle    ! skip the mtau4=mtau4cdep part since it is done           

         DO tau4=max(tau4min,abs(mtau4)),tau4max
          DO j4=j4min,j4max
           DO pari=pari4Nmin,pari4Nmax,2          

C     for given total quantum numbers run 34 quantum numbers 
            DO alphaNN34=1,alphaNNmax
             CALL getNNqn(alphaNN34,l34,s34,j34,t34,mt34)

C     and 12 quantum numbers 
             DO alphaNN12=1,alphaNNmax
              CALL getNNqn(alphaNN12,l12,s12,j12,t12,mt12)
            
C     run inter cluster quantum numbers 
              DO I=abs(j4-j34),min(Imax,j4+j34) 
               DO lam=abs(I-j12),min(I+j12,lammax)
                IF(lam+l12+l34 .GT.lsummax22) CYCLE
C     check whether channels fit 
                IF(     pari.EQ.(-1)**(l12+l34+lam)
     $           .AND.  tau4.GE.abs(t12-t34)
     $           .AND.  tau4.LE.(t12+t34)     ) THEN 
               
                 beta=beta+1

                 qnbeta4N22(1,beta)=l12
                 qnbeta4N22(2,beta)=s12
                 qnbeta4N22(3,beta)=j12
                 qnbeta4N22(4,beta)=t12
                 qnbeta4N22(5,beta)=l34
                 qnbeta4N22(6,beta)=s34
                 qnbeta4N22(7,beta)=j34
                 qnbeta4N22(8,beta)=t34
                 qnbeta4N22(9,beta)=lam
                 qnbeta4N22(10,beta)=I
                 qnbeta4N22(11,beta)=j4
                 qnbeta4N22(12,beta)=tau4
                 qnbeta4N22(13,beta)=mtau4
                 qnbeta4N22(14,beta)=alphaNN12
                 qnbeta4N22(15,beta)=alphaNN34
                 qnbeta4N22(16,beta)=pari
                 
                END IF
               END DO           ! lam
              END DO            ! I
            
             END DO             ! alphaNN12
            END DO              ! alphaNN34
           END DO               ! pari
          END DO                ! j4 
         END DO                 ! tau4 
         
       
        END DO                  ! mtau4
       END IF                   ! cdep4N

#ifdef DEBUG
       IF(beta.NE.beta4N22cdepmax) STOP 'inconsistent beta4N'
#endif 



       RETURN
C     regular end of subroutine 
 100   CONTINUE 
C     treatment of I/O errors 
       STOP 'problem with book-para.dat'

      END SUBROUTINE prepalpha4N

! now I/O routines for channels follow 
! the routines write the considered channel set to 
! to the group given as argument 
! file/group need to be opened beforehand 
   
! writing first       
      
! 1. write NN channels       
! the channels and the parameters of the bookkeeping are written:
!   alphaNNmax,alphaNNcdepmax,l12max,j12max,
!   mt12min,mt12max,evenNN,cdepNN
      
      
      SUBROUTINE writennchannels(group_amp_id)
       IMPLICIT NONE
       INTEGER(HID_T) group_amp_id
       
       CALL write_scalar_int(group_amp_id,'alphaNNmax',alphaNNmax,master)
       CALL write_scalar_int(group_amp_id,'alphaNNcdepmax',alphaNNcdepmax,master)
       CALL write_scalar_int(group_amp_id,'l12NNmax',l12max,master)
       CALL write_scalar_int(group_amp_id,'j12NNmax',j12max,master)
       CALL write_scalar_int(group_amp_id,'mt12NNmin',mt12min,master)
       CALL write_scalar_int(group_amp_id,'mt12NNmax',mt12max,master)
       
       IF(evenNN) THEN        
        CALL write_scalar_int(group_amp_id,'evenNN',1,master)
       ELSE
        CALL write_scalar_int(group_amp_id,'evenNN',0,master)
       ENDIF
       
       IF(cdepNN) THEN        
        CALL write_scalar_int(group_amp_id,'cdepNN',1,master)
       ELSE
        CALL write_scalar_int(group_amp_id,'cdepNN',0,master)
       ENDIF
       
       CALL write_2D_int(group_amp_id,'qnalphaNN',5,alphaNNcdepmax,qnalphaNN,master) 
       
      END SUBROUTINE writennchannels

! 2. write 3N channels       
! the channels and the parameters of the bookkeeping are written:
!   alpha3Nmax,alpha3Ncdepmax,l3max,I3max,j3min,j3max,
!       tau3min,tau3max,mtau3min,mtau3max,cdep3N,
!       pari3Nmin,pari3Nmax 
!   NN channels are written, too      
      
      
      SUBROUTINE write_3n_channels(group_amp_id)
       IMPLICIT NONE
       INTEGER(HID_T) group_amp_id
       
       CALL writennchannels(group_amp_id)
       
       CALL write_scalar_int(group_amp_id,'alpha3Nmax',alpha3Nmax,master)
       CALL write_scalar_int(group_amp_id,'alpha3Ncdepmax',alpha3Ncdepmax,master)
       CALL write_scalar_int(group_amp_id,'l33Nmax',l3max,master)
       CALL write_scalar_int(group_amp_id,'I33Nmax',I3max,master)
       CALL write_scalar_int(group_amp_id,'j33Nmin',j3min,master)
       CALL write_scalar_int(group_amp_id,'j33Nmax',j3max,master)
       CALL write_scalar_int(group_amp_id,'tau33Nmin',tau3min,master)
       CALL write_scalar_int(group_amp_id,'tau33Nmax',tau3max,master)
       CALL write_scalar_int(group_amp_id,'mtau33Nmin',mtau3min,master)
       CALL write_scalar_int(group_amp_id,'mtau33Nmax',mtau3max,master)
       CALL write_scalar_int(group_amp_id,'pari3Nmin',pari3Nmin,master)
       CALL write_scalar_int(group_amp_id,'pari3Nmax',pari3Nmax,master)
       
       IF(cdep3N) THEN        
        CALL write_scalar_int(group_amp_id,'cdep3N',1,master)
       ELSE
        CALL write_scalar_int(group_amp_id,'cdep3N',0,master)
       ENDIF
       
       CALL write_2D_int(group_amp_id,'qnalpha3N',11,alpha3Ncdepmax,qnalpha3N,master) 
       
      END SUBROUTINE write_3n_channels

! 3. write 4N31 channels       
! the channels and the parameters of the bookkeeping are written:
!   alpha4N31max,alpha4N31cdepmax,l4max,I4max,j4min,j4max,
!       tau4min,tau4max,mtau4min,mtau4max,cdep4N,
!       pari4Nmin,pari4Nmax 
!    3N channels are written, too       
      
      
      SUBROUTINE write_4n31_channels(group_amp_id)
       IMPLICIT NONE
       INTEGER(HID_T) group_amp_id
       
       CALL write_3n_channels(group_amp_id)
       
       CALL write_scalar_int(group_amp_id,'alpha4N31max',alpha4N31max,master)
       CALL write_scalar_int(group_amp_id,'alpha4N31cdepmax',alpha4N31cdepmax,master)
       CALL write_scalar_int(group_amp_id,'l44N31max',l4max,master)
       CALL write_scalar_int(group_amp_id,'I44N31Nmax',I4max,master)
       CALL write_scalar_int(group_amp_id,'j44N31min',j4min,master)
       CALL write_scalar_int(group_amp_id,'j44N31max',j4max,master)
       CALL write_scalar_int(group_amp_id,'lsummax4N31',lsummax31,master)
       CALL write_scalar_int(group_amp_id,'tau44N31Nmin',tau4min,master)
       CALL write_scalar_int(group_amp_id,'tau44N31Nmax',tau4max,master)
       CALL write_scalar_int(group_amp_id,'mtau44N31Nmin',mtau4min,master)
       CALL write_scalar_int(group_amp_id,'mtau44N31max',mtau4max,master)
       CALL write_scalar_int(group_amp_id,'pari4N31min',pari4Nmin,master)
       CALL write_scalar_int(group_amp_id,'pari4N31max',pari4Nmax,master)
       
       IF(cdep4N) THEN        
        CALL write_scalar_int(group_amp_id,'cdep4N',1,master)
       ELSE
        CALL write_scalar_int(group_amp_id,'cdep4N',0,master)
       ENDIF
       
       CALL write_2D_int(group_amp_id,'qnalpha4N31',16,alpha4N31cdepmax,qnalpha4N31,master) 
       
      END SUBROUTINE write_4n31_channels      

! 3. write 4N22 channels       
! the channels and the parameters of the bookkeeping are written:
!   beta4N22max,beta4N22cdepmax,lammax,Imax     
      
      SUBROUTINE write_4n22_channels(group_amp_id)
       IMPLICIT NONE
       INTEGER(HID_T) group_amp_id
       
       CALL writennchannels(group_amp_id)
       
       CALL write_scalar_int(group_amp_id,'beta4N22max',beta4N22max,master)
       CALL write_scalar_int(group_amp_id,'beta4N22cdepmax',beta4N22cdepmax,master)
       CALL write_scalar_int(group_amp_id,'lam4N22max',lammax,master)
       CALL write_scalar_int(group_amp_id,'I4N22max',Imax,master)
       CALL write_scalar_int(group_amp_id,'j44N22min',j4min,master)
       CALL write_scalar_int(group_amp_id,'j44N22max',j4max,master)
       CALL write_scalar_int(group_amp_id,'lsummax4N22',lsummax22,master)
       CALL write_scalar_int(group_amp_id,'tau44N22min',tau4min,master)
       CALL write_scalar_int(group_amp_id,'tau44N22max',tau4max,master)
       CALL write_scalar_int(group_amp_id,'mtau44N22min',mtau4min,master)
       CALL write_scalar_int(group_amp_id,'mtau44N22max',mtau4max,master)
       CALL write_scalar_int(group_amp_id,'pari4N22min',pari4Nmin,master)
       CALL write_scalar_int(group_amp_id,'pari4N22max',pari4Nmax,master)
       
       IF(cdep4N) THEN        
        CALL write_scalar_int(group_amp_id,'cdep4N',1,master)
       ELSE
        CALL write_scalar_int(group_amp_id,'cdep4N',0,master)
       ENDIF
       
       CALL write_2D_int(group_amp_id,'qnbeta4N22',16,beta4N22cdepmax,qnbeta4N22,master) 
       
      END SUBROUTINE write_4n22_channels      
      
! 1. read NN channels       
! the channels and the parameters of the bookkeeping are read in.
! The parameters are printed: 
!   alphaNNmax,alphaNNcdepmax,l12NNmax,j12NNmax,
!          t12NNmin,t12NNmax,mt12NNmin,mt12NNmax,evenNN,cdepNN
!  and channel array qnalphaNN_read is alloced 
!    and used to store the channels 
!  the dimension is given back via  alphaNNmax_read and 
!   alphaNNcdepmax_read       
      
      SUBROUTINE readnnchannels(group_amp_id,qnalphaNN_read,
     $                           alphaNNmax_read,alphaNNcdepmax_read)
       IMPLICIT NONE
       INTEGER(HID_T) group_amp_id
       INTEGER,ALLOCATABLE :: qnalphaNN_read(:,:) 
       INTEGER alphaNNmax_read,alphaNNcdepmax_read 
       INTEGER logint,alpha
       INTEGER j12NNmax_read,l12NNmax_read,
     $         mt12NNmin_read,mt12NNmax_read  
       LOGICAL cdepNN_read,evenNN_read 
!  read in parameters        
       CALL read_scalar_int(group_amp_id,'alphaNNmax',
     $   alphaNNmax_read,commall)
       CALL read_scalar_int(group_amp_id,'alphaNNcdepmax',
     $   alphaNNcdepmax_read,commall)
       CALL read_scalar_int(group_amp_id,'l12NNmax',
     $   l12NNmax_read,commall)
       CALL read_scalar_int(group_amp_id,'j12NNmax',
     $   j12NNmax_read,commall)
       CALL read_scalar_int(group_amp_id,'mt12NNmin',
     $   mt12NNmin_read,commall)
       CALL read_scalar_int(group_amp_id,'mt12NNmax',
     $   mt12NNmax_read,commall)
              
       CALL read_scalar_int(group_amp_id,'evenNN',logint,commall)
       evenNN_read=logint.eq.1 
       CALL read_scalar_int(group_amp_id,'cdepNN',logint,commall)
       cdepNN_read=logint.eq.1 
       
! allocate array for quantum numbers        
       ALLOCATE(qnalphaNN_read(5,alphaNNcdepmax_read))
 
! read quantum numbers        
       CALL read_2D_int(group_amp_id,'qnalphaNN',5,alphaNNcdepmax_read,
     $    qnalphaNN_read,commall) 

! print parameters and quantum numbers to check file content        
       IF(master) THEN
        WRITE(*,*) 'Read in parameters of NN bookkeeping'
        WRITE(*,'(A,I5)') 'j12NNmax        = ',j12NNmax_read
        WRITE(*,'(A,I5)') 'l12NNmax        = ',l12NNmax_read
        WRITE(*,*)        'evenNN        = ',evenNN_read
        WRITE(*,*)        'cdepNN        = ',cdepNN_read
        WRITE(*,*)        'mt12NN min-max  = ',
     $                     mt12NNmin_read,mt12NNmax_read
        WRITE(*,*) 
        WRITE(*,'(A,I5)') 'alphaNNmax     = ',alphaNNmax_read
        WRITE(*,'(A,I5)') 'alphaNNcdepmax = ',alphaNNcdepmax_read
        WRITE(*,*)  


        WRITE(*,*) 
        WRITE(*,*) 'PW NN bookkeeping:'
        WRITE(*,'(A10,A7,5A5)') 
     $      '#','alphaNN','l12','s12','j12','t12','mt12'
        WRITE(*,*) 

        DO alpha=1,alphaNNcdepmax_read
         WRITE(*,'(A,I7,5I5)') 'ALNNREAD: ',
     $        alpha,qnalphaNN_read(1,alpha),
     $              qnalphaNN_read(2,alpha),
     $              qnalphaNN_read(3,alpha),
     $              qnalphaNN_read(4,alpha),
     $              qnalphaNN_read(5,alpha)
        END DO       
       END IF
       
      END SUBROUTINE readnnchannels

! 2. read 3N channels       
! the channels and the parameters of the bookkeeping are read in.
! The parameters are printed: 
!  alpha3Nmax,alpha3Ncdepmax,l33Nmax,I33Nmax,j33Nmin,j33Nmax,
!       tau33Nmin,tau33Nmax,mtau33Nmin,mtau33Nmax,cdep3N,
!       pari3Nmin,pari3Nmax  
!  and channel array qnalpha3N_read is alloced 
!    and used to store the channels 
!  the dimension is given back via  alpha3Nmax_read and 
!   alpha3Ncdepmax_read       
      
      SUBROUTINE read_3n_channels(group_amp_id,qnalpha3N_read,
     $                           alpha3Nmax_read,alpha3Ncdepmax_read)
       IMPLICIT NONE
       INTEGER(HID_T) group_amp_id
       INTEGER,ALLOCATABLE :: qnalpha3N_read(:,:) 
       INTEGER logint,alpha
       INTEGER alpha3Nmax_read,alpha3Ncdepmax_read,
     $         l33Nmax_read,I33Nmax_read,
     $         j33Nmin_read,j33Nmax_read,
     $         tau33Nmin_read,tau33Nmax_read,
     $         mtau33Nmin_read,mtau33Nmax_read,
     $         pari3Nmin_read,pari3Nmax_read   
       LOGICAL cdep3N_read
!  read in parameters        
       CALL read_scalar_int(group_amp_id,'alpha3Nmax',
     $   alpha3Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'alpha3Ncdepmax',
     $   alpha3Ncdepmax_read,commall)
       CALL read_scalar_int(group_amp_id,'l33Nmax',
     $   l33Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'I33Nmax',
     $   I33Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'j33Nmin',
     $   j33Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'j33Nmax',
     $   j33Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'tau33Nmin',
     $   tau33Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'tau33Nmax',
     $   tau33Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'mtau33Nmin',
     $   mtau33Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'mtau33Nmax',
     $   mtau33Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'pari3Nmin',
     $   pari3Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'pari3Nmax',
     $   pari3Nmax_read,commall)
     
     
     
     
       CALL read_scalar_int(group_amp_id,'cdep3N',logint,commall)
       cdep3N_read=logint.eq.1 
       
       
! allocate array for quantum numbers        
       ALLOCATE(qnalpha3N_read(11,alpha3Ncdepmax_read))
 
! read quantum numbers        
       CALL read_2D_int(group_amp_id,'qnalpha3N',11,
     $    alpha3Ncdepmax_read,
     $    qnalpha3N_read,commall) 

! print parameters and quantum numbers to check file content        
       IF(master) THEN
        WRITE(*,*) 'Read in parameters of 3N bookkeeping'
        WRITE(*,*) 'Parameters of 3N bookkeeping'
        WRITE(*,'(A,I5)')  'l33Nmax              = ',l33Nmax_read
        WRITE(*,'(A,I5)')  'I33Nmax              = ',I33Nmax_read
        WRITE(*,'(A,2I5)') 'j33Nmin,j33Nmax     = ',
     $                      j33Nmin_read,j33Nmax_read
        WRITE(*,'(A,2I5)') 'pari3Nmin,pari3Nmax = ',
     $                      pari3Nmin_read,pari3Nmax_read
        WRITE(*,'(A,2I5)') 'tau33Nmin,tau33Nmax = ',
     $                      tau33Nmin_read,tau33Nmax_read
        WRITE(*,*)         'cdep3N               = ',cdep3N_read
        WRITE(*,*)         'mtau33N min-max      = ',
     $                     mtau33Nmin_read,mtau33Nmax_read
        WRITE(*,*) 
        WRITE(*,'(A,I5)') 'alpha3Nmax     = ',alpha3Nmax_read
        WRITE(*,'(A,I5)') 'alpha3Ncdepmax = ',alpha3Ncdepmax_read
        WRITE(*,*)  


        WRITE(*,*) 
        WRITE(*,*) 'PW YN bookkeeping:'
        WRITE(*,'(A10,A10,4A5,2X,5A5,2X,2A5)') 
     $      '#','alpha3N','l12','s12','j12','t12',
     $      'l3','I3','j3','tau3','mtau3','alNN','pari'
        WRITE(*,*) 

        DO alpha=1,alpha3Ncdepmax_read
         WRITE(*,'(A,I10,4I5,2X,5I5,2X,2I5)') 'AL3NREAD: ',
     $        alpha,qnalpha3N_read(1,alpha),
     $              qnalpha3N_read(2,alpha),
     $              qnalpha3N_read(3,alpha),
     $              qnalpha3N_read(4,alpha),
     $              qnalpha3N_read(5,alpha),
     $              qnalpha3N_read(6,alpha),
     $              qnalpha3N_read(7,alpha),
     $              qnalpha3N_read(8,alpha),
     $              qnalpha3N_read(9,alpha),
     $              qnalpha3N_read(10,alpha),
     $              qnalpha3N_read(11,alpha)

        END DO       
       END IF
       
      END SUBROUTINE read_3n_channels

! 3. read 4N31 channels       
! the channels and the parameters of the bookkeeping are read in.
! The parameters are printed: 
!  alpha4N31max,alpha4N31cdepmax,l44Nmax,I44Nmax,j44Nmin,j44Nmax,
!       tau44Nmin,tau44Nmax,mtau44Nmin,mtau44Nmax,cdep4N,
!       pari4Nmin,pari4Nmax  
!  and channel array qnalpha4N_read is alloced 
!    and used to store the channels 
!  the dimension is given back via  alpha4Nmax_read and 
!   alpha4Ncdepmax_read       
      
      SUBROUTINE read_4n31_channels(group_amp_id,qnalpha4N_read,
     $                           alpha4Nmax_read,alpha4Ncdepmax_read)
       IMPLICIT NONE
       INTEGER(HID_T) group_amp_id
       INTEGER,ALLOCATABLE :: qnalpha4N_read(:,:) 
       INTEGER logint,alpha
       INTEGER alpha4Nmax_read,alpha4Ncdepmax_read,
     $         l44Nmax_read,I44Nmax_read,
     $         j44Nmin_read,j44Nmax_read,
     $         tau44Nmin_read,tau44Nmax_read,
     $         mtau44Nmin_read,mtau44Nmax_read,
     $         pari4Nmin_read,pari4Nmax_read,lsummax31_read   
       LOGICAL cdep4N_read
!  read in parameters        
       CALL read_scalar_int(group_amp_id,'alpha4N31max',
     $   alpha4Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'alpha4N31cdepmax',
     $   alpha4Ncdepmax_read,commall)
       CALL read_scalar_int(group_amp_id,'l44N31max',
     $   l44Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'I44N31Nmax',
     $   I44Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'j44N31min',
     $   j44Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'j44N31max',
     $   j44Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'lsummax4N31',
     $   lsummax31_read,commall)
       CALL read_scalar_int(group_amp_id,'tau44N31Nmin',
     $   tau44Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'tau44N31Nmax',
     $   tau44Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'mtau44N31Nmin',
     $   mtau44Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'mtau44N31max',
     $   mtau44Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'pari4N31min',
     $   pari4Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'pari4N31max',
     $   pari4Nmax_read,commall)
     
     
     
     
       CALL read_scalar_int(group_amp_id,'cdep4N',logint,commall)
       cdep4N_read=logint.eq.1 
       
       
! allocate array for quantum numbers        
       ALLOCATE(qnalpha4N_read(16,alpha4Ncdepmax_read))
 
! read quantum numbers        
       CALL read_2D_int(group_amp_id,'qnalpha4N31',16,
     $    alpha4Ncdepmax_read,
     $    qnalpha4N_read,commall) 

! print parameters and quantum numbers to check file content        
       IF(master) THEN
        WRITE(*,*) 'Read in parameters of 4N31 bookkeeping'
        WRITE(*,*) 'Parameters of 4N31 bookkeeping'
        WRITE(*,'(A,I5)')  'l44Nmax              = ',l44Nmax_read
        WRITE(*,'(A,I5)')  'I44Nmax              = ',I44Nmax_read
        WRITE(*,'(A,2I5)') 'j44Nmin,j44Nmax      = ',
     $                      j44Nmin_read,j44Nmax_read
        WRITE(*,'(A,2I5)') 'lsummax31            = ',lsummax31_read     
        WRITE(*,'(A,2I5)') 'pari4Nmin,pari4Nmax = ',
     $                      pari4Nmin_read,pari4Nmax_read
        WRITE(*,'(A,2I5)') 'tau44Nmin,tau44Nmax = ',
     $                      tau44Nmin_read,tau44Nmax_read
        WRITE(*,*)         'cdep4N               = ',cdep4N_read
        WRITE(*,*)         'mtau44N min-max      = ',
     $                     mtau44Nmin_read,mtau44Nmax_read
        WRITE(*,*) 
        WRITE(*,'(A,I5)') 'alpha4N31max     = ',alpha4Nmax_read
        WRITE(*,'(A,I5)') 'alpha4N31cdepmax = ',alpha4Ncdepmax_read
        WRITE(*,*)  


        WRITE(*,*) 
        WRITE(*,*) 'PW 4N31 bookkeeping:'
        WRITE(*,'(A10,A10,4A5,2X,4A5,2X,5A5,2X,3A5)') 
     $      '#','alpha4N31','l12','s12','j12','t12',
     $      'l3','I3','j3','tau3','l4','I4','j4','tau4','mtau4',
     $      'alNN','al3N','pari'
        WRITE(*,*) 

        DO alpha=1,alpha4Ncdepmax_read
         WRITE(*,'(A,I10,4I5,2X,4I5,2X,5I5,2X,3I5)') 'AL4N31READ: ',
     $        alpha,qnalpha4N_read(1,alpha),
     $              qnalpha4N_read(2,alpha),
     $              qnalpha4N_read(3,alpha),
     $              qnalpha4N_read(4,alpha),
     $              qnalpha4N_read(5,alpha),
     $              qnalpha4N_read(6,alpha),
     $              qnalpha4N_read(7,alpha),
     $              qnalpha4N_read(8,alpha),
     $              qnalpha4N_read(9,alpha),
     $              qnalpha4N_read(10,alpha),
     $              qnalpha4N_read(11,alpha),
     $              qnalpha4N_read(12,alpha),
     $              qnalpha4N_read(13,alpha),
     $              qnalpha4N_read(14,alpha),
     $              qnalpha4N_read(15,alpha),
     $              qnalpha4N_read(16,alpha)

        END DO       
       END IF
       
      END SUBROUTINE read_4n31_channels

! 4. read 4N22 channels       
! the channels and the parameters of the bookkeeping are read in.
! The parameters are printed: 
!  alpha4N22max,alpha4N22cdepmax,lam4Nmax,I4Nmax
!  and channel array qnalpha4N_read is alloced 
!    and used to store the channels 
!  the dimension is given back via  alpha4Nmax_read and 
!   alpha4Ncdepmax_read       
      
      SUBROUTINE read_4n22_channels(group_amp_id,qnalpha4N_read,
     $                           alpha4Nmax_read,alpha4Ncdepmax_read)
       IMPLICIT NONE
       INTEGER(HID_T) group_amp_id
       INTEGER,ALLOCATABLE :: qnalpha4N_read(:,:) 
       INTEGER alpha
       INTEGER alpha4Nmax_read,alpha4Ncdepmax_read,
     $         lam4Nmax_read,I4Nmax_read,
     $         j44Nmin_read,j44Nmax_read,
     $         tau44Nmin_read,tau44Nmax_read,
     $         mtau44Nmin_read,mtau44Nmax_read,
     $         pari4Nmin_read,pari4Nmax_read,lsummax22_read 
       LOGICAL cdep4N_read
       INTEGER logint
     
!  read in parameters        
       CALL read_scalar_int(group_amp_id,'beta4N22max',
     $   alpha4Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'beta4N22cdepmax',
     $   alpha4Ncdepmax_read,commall)
       CALL read_scalar_int(group_amp_id,'lam4N22max',
     $   lam4Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'I4N22max',
     $   I4Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'j44N22min',
     $   j44Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'j44N22max',
     $   j44Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'lsummax4N22',
     $   lsummax22_read,commall)
       CALL read_scalar_int(group_amp_id,'tau44N22min',
     $   tau44Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'tau44N22max',
     $   tau44Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'mtau44N22min',
     $   mtau44Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'mtau44N22max',
     $   mtau44Nmax_read,commall)
       CALL read_scalar_int(group_amp_id,'pari4N22min',
     $   pari4Nmin_read,commall)
       CALL read_scalar_int(group_amp_id,'pari4N22max',
     $   pari4Nmax_read,commall)
         
       CALL read_scalar_int(group_amp_id,'cdep4N',logint,commall)
       cdep4N_read=logint.eq.1 
       
! allocate array for quantum numbers        
       ALLOCATE(qnalpha4N_read(16,alpha4Ncdepmax_read))
 
! read quantum numbers        
       CALL read_2D_int(group_amp_id,'qnbeta4N22',16,
     $    alpha4Ncdepmax_read,
     $    qnalpha4N_read,commall) 

! print parameters and quantum numbers to check file content        
       IF(master) THEN
        WRITE(*,*) 'Read in parameters of 4N22 bookkeeping'
        WRITE(*,*) 'Parameters of 4N22 bookkeeping'
        WRITE(*,'(A,I5)')  'lam4Nmax              = ',lam4Nmax_read
        WRITE(*,'(A,I5)')  'I44Nmax              = ',I4Nmax_read
        WRITE(*,'(A,2I5)') 'j44Nmin,j44Nmax     = ',
     $                      j44Nmin_read,j44Nmax_read
        WRITE(*,'(A,2I5)') 'lsummax22            = ',lsummax22_read      
        WRITE(*,'(A,2I5)') 'pari4Nmin,pari4Nmax = ',
     $                      pari4Nmin_read,pari4Nmax_read
        WRITE(*,'(A,2I5)') 'tau44Nmin,tau44Nmax = ',
     $                      tau44Nmin_read,tau44Nmax_read
        WRITE(*,*)         'cdep4N               = ',cdep4N_read
        WRITE(*,*)         'mtau44N min-max      = ',
     $                     mtau44Nmin_read,mtau44Nmax_read
        WRITE(*,*) 
        
        WRITE(*,'(A,I5)') 'beta4N22max     = ',alpha4Nmax_read
        WRITE(*,'(A,I5)') 'beta4N22cdepmax = ',alpha4Ncdepmax_read
        WRITE(*,*)  


        WRITE(*,*) 
        WRITE(*,*) 'PW 4N22 bookkeeping:'
        WRITE(*,'(A10,A10,4A5,2X,4A5,2X,5A5,2X,3A5)') 
     $      '#','alpha4N22','l12','s12','j12','t12',
     $      'l34','s34','j34','t34','lam','I','j4','tau4','mtau4',
     $      'al12','al34','pari'
        WRITE(*,*) 

        DO alpha=1,alpha4Ncdepmax_read
         WRITE(*,'(A,I10,4I5,2X,4I5,2X,5I5,2X,3I5)') 'AL4N22READ: ',
     $        alpha,qnalpha4N_read(1,alpha),
     $              qnalpha4N_read(2,alpha),
     $              qnalpha4N_read(3,alpha),
     $              qnalpha4N_read(4,alpha),
     $              qnalpha4N_read(5,alpha),
     $              qnalpha4N_read(6,alpha),
     $              qnalpha4N_read(7,alpha),
     $              qnalpha4N_read(8,alpha),
     $              qnalpha4N_read(9,alpha),
     $              qnalpha4N_read(10,alpha),
     $              qnalpha4N_read(11,alpha),
     $              qnalpha4N_read(12,alpha),
     $              qnalpha4N_read(13,alpha),
     $              qnalpha4N_read(14,alpha),
     $              qnalpha4N_read(15,alpha),
     $              qnalpha4N_read(16,alpha)

        END DO       
       END IF
       
      END SUBROUTINE read_4n22_channels

      END MODULE pwbook

