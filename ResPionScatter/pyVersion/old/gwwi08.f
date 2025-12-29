C GETS F AND G/SIN(THETA) AS FUNCTION OF E(MEV) AND C(COS THETA)        
C AMPLITUDES (F, G/SIN) ARE IN FERMIS                                   
C Modified 10/05/01 to do type 1 (CMKM form)
C WI08 : R.L.Workman et al., PRC 86, 035202 (2012).
      DIMENSION TR(4,8),TI(4,8),F(2),G(2),TTL(18),FASE(6,4)
      DATA FASE/4H    ,4HP11 ,4HD13 ,4HF15 ,4HG17 ,4HH19 
     C,4HS11 ,4HP13 ,4HD15 ,4HF17 ,4HG19 ,4HH111
     C,4H    ,4HP31 ,4HD33 ,4HF35 ,4HG37 ,4HH39 
     C,4HS31 ,4HP33 ,4HD35 ,4HF37 ,4HG39 ,4HH311 /
      DATA WI,WT,HBC/139.65,938.256,197.32/
      OPEN(7,FILE='SAID.TMP',STATUS='UNKNOWN')
      CALL PNAMPV(50.0,-1.0,0,TR,TI,F,G,TTL)                            
      WRITE(*,105) TTL                                                  
      WRITE(7,105) TTL                                                  
  105 FORMAT(//,' Partial waves are shown for L=0-5 and for N=1-4.'
     C,/' (I,J)N=(1/2,-1/2)1, (1/2,+1/2)2, (3/2,-1/2)3, (3/2,+1/2)4 '
     C,/' Amplitudes (f and g) are in Fm; DSG(mb)=10*(f**2+g**2)'
     C,//' Reaction is labeled as 1(pi+p), 2(pi-p), 3(cxs), 0(HADRONIC)'
     C,//' Solution is:',//,' ',18A4)
      IR=0
      DO 20 I=1,1
c     E=400+100*I
      E=100
      CALL PNFIXD(E,0,TR,TI,TTL)
      WRITE(*,102) E,IR                                                 
      WRITE(7,102) E,IR                                                 
  102 FORMAT(//,'  Tlab=',F8.2,'  React=',I2)
      DO 21 N=1,4                                                        
      WRITE(*,107) (FASE(J,N),J=1,6)
      WRITE(7,107) (FASE(J,N),J=1,6)
  107 FORMAT(/8X,6(3X,A4,3X))
      WRITE(*,101) N,(TR(N,J),J=1,6),(TI(N,J),J=1,6)                   
      WRITE(7,101) N,(TR(N,J),J=1,6),(TI(N,J),J=1,6)                    
  101 FORMAT(' Tr(',I1,')',6(2X,F8.4),/,' Ti   ',6(2X,F8.4))            
      DO 22 J=1,6
   22 CALL DOFT(TR(N,J),TI(N,J),0.0)
      WRITE(*,108) (TR(N,J),J=1,6),(TI(N,J),J=1,6)
      WRITE(7,108) (TR(N,J),J=1,6),(TI(N,J),J=1,6)
  108 FORMAT(' Delta',6(2X,F8.4),/' 1-e2 ',6(2X,F8.4))
   21 CONTINUE                                                          
   20 CONTINUE
      IF(E.GT.500.0) GO TO 99
      DO 4 II=1,7                                                       
      IR=II-1                                                           
      IF(II.GT.4) IR=4-II                                               
c     E=500.0                                                           
      E=100
      DO 5 I=1,1                                                        
      A=0.0                                                             
      DO 6 K=1,7                                                        
      CALL PNAMPV(E,A,IR,TR,TI,F,G,TTL)                                 
      IF(K.GT.1) GO TO 7                                                
      WRITE(*,102) E,IR                                                 
      WRITE(7,102) E,IR                                                 
      DO 1 N=1,4                                                        
      WRITE(*,107) (FASE(J,N),J=1,6)
      WRITE(7,107) (FASE(J,N),J=1,6)
      WRITE(*,101) N,(TR(N,J),J=1,6),(TI(N,J),J=1,6)                   
      WRITE(7,101) N,(TR(N,J),J=1,6),(TI(N,J),J=1,6)                   
    1 CONTINUE                                                          
      IF(IR.NE.0) WRITE(*,104)                                          
      IF(IR.NE.0) WRITE(7,104)                                          
  104 FORMAT('  Theta     DSG           f               g/sin(theta) ')
    7 IF(IR.EQ.0) GO TO 6                                               
      S=SIN(0.0174532*A)
      DSG=10.0*(F(1)**2+F(2)**2+S**2*(G(1)**2+G(2)**2))
      WRITE(*,106) A,DSG,F,G                                            
      WRITE(7,106) A,DSG,F,G                                            
c 106 FORMAT(F7.2,F8.3,2F8.4,4X,2F8.4)                                  
  106 FORMAT(F7.1,2X,F8.5,2X,2F9.6,'j,',4X,2F9.6,'j,')                                  
    6 A=A+30.0                                                         
    5 E=E+100.0                                                         
    4 CONTINUE                                                          
   99 STOP                                                              
      END                                                               
C **************************************************************        
      SUBROUTINE PNAMPV(E,A,IRX,TR,TI,F,G,TTL)                          
C GET PI-N AMPLITUDES FOR A VPI SOLUTION 11/8/94  ARNDT                 
C NOTE, IF A=0 BOTH COULOMB AND PHASE-ROTATIONS ARE SUPPRESSED          
C IF IR=0(FOR PI-N) THEN F,G=0 (ONLY T-MTX IS CALCULATED)               
C IR=0(Hadronic), 1(Pi+P), 2(Pi-P), 3(CXS); If IR < 0 kill Coulomb IR=-I
      DIMENSION TR(4,8),TI(4,8),F(2),G(2),TTL(18),SPH(9)                
      DATA PI/3.1415927/                                                
      DATA WI,WT,LMX/139.65,938.256,8/                                  
      DATA IRM,EM/-1,-1.0/                                              
      SAVE
      F(1)=0.                                                           
      F(2)=0.                                                           
      G(1)=0.                                                           
      G(2)=0.                                                           
      IF(IRX.NE.IRM) EM=-1.0                                            
      IF(E.EQ.EM) GO TO 15                                              
      IRM=IRX                                                           
      IR=IRX                                                            
      IF(IR.LT.0) IR=-IR                                                
      EM=E                                                              
      CALL PNFIXD(E,IR,TR,TI,TTL)                                       
      IF(A.LT.0.0.OR.IR.EQ.0) GO TO 99                                  
      S0=(WI+WT)**2                                                     
      W=SQRT(S0+2.*WT*E)                                                
      QPQ=(W**2-1074.7**2)*(W**2-804.61**2)                             
      QPQ=SQRT(QPQ/(W**2-S0)/(W**2-(WT-WI)**2))                         
      XKMEV=SQRT(E*(E+2.*WI)/((1.+WI/WT)**2+2.*E/WT))                   
      IF(IR.EQ.3) XKMEV=XKMEV*SQRT(QPQ)                                 
      XKM=197.32/XKMEV                                                  
      SPH(1)=0.                                                         
      QR=0.0                                                            
      IF(IR.EQ.1) QR=1.0                                                
      IF(IR.GT.1) QR=-1.0                                               
      IF(IRX.LT.0) QR=0.0                                               
      PZR=SQRT(XKMEV**2+WT**2)                                          
      QZR=SQRT(XKMEV**2+WI**2)                                          
      W=PZR+QZR                                                         
      Z=0.007297348*(QZR*PZR+XKMEV**2)/XKMEV/W                          
      ZZ=0.                                                             
      DO 17 L=2,8                                                       
      XL=L-1                                                            
      ZZ=ZZ+ATAN(Z/XL)                                                  
   17 SPH(L)=SIN(ZZ*SGLSCL(E,L-1,WI))                                   
   15 C=.0174532*A                                                      
      C=COS(C)                                                          
      S=SQRT(1.-C**2)                                                   
      IF(IR.EQ.3) QR=0.0                                                
      IF(QR.NE.0.0) CALL COUL(QR,E,WI,WT,A,F,F(2),G)                    
      XL=0.                                                             
      PL=1.                                                             
      PL1=0.                                                            
      DO 31 L=1,LMX                                                     
      SP=SPH(L)                                                         
      IF(IR.GT.1) SP=-SP                                                
      IF(IR.NE.3) SP=2.0*SP*SQRT(1.0-SP**2)                             
      IF(A.EQ.0.0) SP=0.0                                               
      IF(IRX.LT.0) SP=0.0                                               
      CP=SQRT(1.0-SP**2)                                                
      DO 32 N=1,4                                                       
      IF(IR.EQ.1.AND.N.LT.3) GO TO 32                                   
      TRX=TR(N,L)                                                       
      TIX=TI(N,L)                                                       
      IF(TRX.EQ.0.0.AND.TIX.EQ.0.0) GO TO 32                            
      ZR=TRX                                                            
      ZI=TIX                                                            
      TRX=ZR*CP-ZI*SP                                                   
      TIX=ZR*SP+ZI*CP                                                   
      F1=PL*XKM                                                         
      F2=PL1*XKM                                                        
      IF(N.EQ.1.OR.N.EQ.3) F2=-F2                                       
      IF(N.EQ.1.OR.N.EQ.3) F1=F1*XL                                     
      IF(N.EQ.2.OR.N.EQ.4) F1=F1*(XL+1.)                                
      IF(IR.EQ.1) GO TO 35                                              
      FS=1./3.                                                          
      IF(N.LT.3) FS=2./3.                                               
      IF(IR.EQ.3) FS=SQRT(2.)/3.                                        
      F1=F1*FS                                                          
      F2=F2*FS                                                          
      IF(IR.EQ.3.AND.N.LT.3) F1=-F1                                     
      IF(IR.EQ.3.AND.N.LT.3) F2=-F2                                     
   35 F(1)=F(1)+F1*TRX                                                  
      F(2)=F(2)+F1*TIX                                                  
      G(1)=G(1)+F2*TRX                                                  
      G(2)=G(2)+F2*TIX                                                  
   32 CONTINUE                                                          
      XL=XL+1.                                                          
      IF(L.GT.1) GO TO 36                                               
      PLM=1.                                                            
      PL=C                                                              
      PL1M=0.                                                           
      PL1=1.0                                                           
      GO TO 31                                                          
   36 PX=((2.*XL-1.)*C*PL-(XL-1)*PLM)/XL                                
      PLM=PL                                                            
      PL=PX                                                             
      PX=((2.*XL-1.)*C*PL1-XL*PL1M)/(XL-1.)                             
      PL1M=PL1                                                          
      PL1=PX                                                            
   31 CONTINUE                                                          
C     IF(ABS(C).GT.1.0) GO TO 99                                        
C     S=SQRT(1.0-C**2)                                                  
C     G(1)=G(1)*S                                                       
C     G(2)=G(2)*S                                                       
   99 RETURN                                                            
      END                                                               
C ******************************************************                
      SUBROUTINE PNFIXD(EX,IRR,TRZ,TIZ,NTL)                                
C Get SAID partial waves. Parameters are in DATA statements
C E is Tlab(MeV), IR=0(PiN),1(Pi+P),2(Pi-P),3(Cxs)
C T(N,L) is PW for l=L-1, and N=1(I=1/2,J-), 2(I=1/2,J+),3(I=3/2,J-)
C and 4(I=3/2,J+) eg (N,L)=(2,1) for S11, (4,1) for S31, (4,2) for P33
C (1,4) for F15 .......
C NTL is set on 1st call and is a "title" for the SAID solution encoded
      DIMENSION PP(460),NNTL(18),NNTX(5),PP1(70),PP2(70),PP3(70),PP4(70)
     C,PP5(70),PP6(19)
      DIMENSION TR(4,8),TI(4,8),NFM(4,8),P(30,4,8),TRZ(4,8),TIZ(4,8)
      DIMENSION NTL(18),W1(3),W2(3),DW2(3),V(8,3),VI(8,3),BPL(8)
      CHARACTER HTL*52
      EQUIVALENCE (HTL,NNTL),(PP,PP1),(PP(71),PP2),(PP(141),PP3)
     C,(PP(211),PP4),(PP(281),PP5),(PP(351),PP6)
      DATA NNTX/4H    ,4HPiN ,4HAnal,4Hysis,4H    /
      DATA HTL/'WI08 766276 57255/31876 P+=27190/13344 P-=22702/1196'/
      DATA IMX/369/
      DATA PP1/ 0.22100E+11, 0.36137E+00, -0.18963E+00, 0.00000E+00,
     C 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.17971E+02,-0.28969E+03,
     C 0.48324E+03, 0.16901E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,  
     C 0.96111E+01,-0.61530E+01, 0.00000E+00, 0.68147E+00, 0.72515E+01,
     C 0.00000E+00,-0.26534E+01, 0.18496E+02,-0.21453E+02, 0.47744E+02,
     C 0.73100E+11, 0.55000E+03, 0.24100E+11,-0.23838E+00,-0.42402E+01,
     C 0.43713E+01, 0.00000E+00, 0.00000E+00,-0.57780E+01, 0.14861E+02,
     C 0.00000E+00,-0.10495E+03,-0.19812E+01, 0.00000E+00,-0.16797E+02,
     C 0.00000E+00, 0.00000E+00, 0.16032E+01, 0.21200E+11,-0.14085E+01,
     C 0.13476E+02,-0.13900E+02, 0.00000E+00, 0.00000E+00, 0.53149E+01,
     C-0.24608E+01, 0.16797E+01, 0.00000E+00, 0.00000E+00, 0.42388E+01,
     C-0.24081E+01, 0.00000E+00,-0.36067E+01, 0.22200E+11,-0.47760E+00,
     C 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,-0.10382E+01,
     C 0.40793E+01,-0.64920E+01, 0.12585E+02, 0.34411E+01, 0.00000E+00, 
     C 0.00000E+00/ 
      DATA PP2/ 0.00000E+00,-0.58911E+01, 0.23200E+11,-0.85056E+00,
     C-0.59896E+00,-0.10914E+02, 0.16135E+02, 0.00000E+00,-0.65440E+01,
     C 0.11327E+02, 0.00000E+00, 0.39072E+01, 0.24200E+11, 0.20756E+01,
     C-0.35398E+01, 0.31662E+00, 0.00000E+00, 0.13802E+04, 0.22156E+01,
     C-0.10990E+01, 0.20011E+01, 0.00000E+00, 0.95673E+00,-0.12871E+00,
     C 0.00000E+00, 0.00000E+00, 0.16710E+01, 0.21300E+11, 0.29236E+00,
     C 0.61349E+01,-0.18106E+02, 0.00000E+00, 0.00000E+00, 0.80522E+01,
     C 0.00000E+00,-0.88144E+02, 0.17681E+03,-0.10658E+02,-0.12869E+02,
     C 0.00000E+00, 0.00000E+00,-0.40816E+02, 0.00000E+00, 0.65383E+00,
     C 0.22300E+11, 0.75747E+00,-0.12938E+01, 0.19522E+01, 0.00000E+00,  
     C 0.00000E+00, 0.00000E+00, 0.69812E+01,-0.42417E+02, 0.62856E+02,
     C-0.71467E+00, 0.00000E+00,-0.66570E+01, 0.00000E+00, 0.28884E+01,
     C 0.23300E+11, 0.23273E+00, 0.00000E+00,-0.29315E+01, 0.00000E+00,
     C 0.00000E+00, 0.85223E+00, 0.56647E+01,-0.61346E+01, 0.00000E+00, 
     C 0.23588E+01/
      DATA PP3/ 0.00000E+00,-0.94402E+01, 0.24300E+11,-0.37143E+00,
     C 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.50997E+00,
     C 0.00000E+00,-0.11167E+02, 0.13155E+02,-0.23223E+01, 0.98186E+00,
     C 0.21400E+11, 0.27615E+00, 0.19689E+01, 0.00000E+00, 0.00000E+00, 
     C 0.00000E+00, 0.44751E+00, 0.10608E+02,-0.82973E+02, 0.11923E+03,
     C 0.37889E+01,-0.18430E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00, 0.11985E+00,-0.13046E+03,-0.18221E+04, 0.22400E+11,
     C-0.72618E-01, 0.17841E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00, 0.77748E+00, 0.00000E+00, 0.12110E+01, 0.23400E+11,
     C 0.38234E-01,-0.52850E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C 0.10844E+01, 0.98635E+00,-0.12576E+02, 0.16757E+02, 0.19866E+01,
     C 0.24400E+11, 0.53591E+00, 0.00000E+00, 0.49010E+00, 0.00000E+00,
     C 0.00000E+00/
      DATA PP4/-0.59911E+00, 0.40633E+01,-0.23872E+02, 0.27209E+02,
     C 0.33391E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,-0.10827E+02,
     C 0.21500E+11, 0.15685E+00, 0.11355E+01,-0.39781E+00, 0.00000E+00,
     C 0.00000E+00, 0.15236E+01, 0.00000E+00,-0.65434E+01, 0.61106E+01,
     C 0.22500E+11, 0.22613E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00, 0.39141E+00, 0.45996E+00,-0.83784E+01, 0.76916E+01,
     C 0.23500E+11, 0.12716E+00,-0.20234E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00, 0.14470E+00, 0.70760E+00, 0.29952E+00, 0.24500E+11,
     C-0.68867E-01, 0.17605E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00, 0.86708E+00, 0.00000E+00, 0.16417E+01, 0.12169E+01,
     C 0.21600E+11, 0.17329E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00, 0.00000E+00,-0.54227E+00,-0.29442E+01, 0.00000E+00,
     C 0.54009E+01,-0.29036E+01, 0.00000E+00, 0.00000E+00, 0.38925E+01,
     C 0.22600E+11, 0.86222E-01,-0.14068E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00/
      DATA PP5/ 0.00000E+00, 0.45945E+00,-0.19433E+02, 0.17842E+02,
     C 0.23600E+11, 0.10598E+00,-0.38295E+00, 0.00000E+00, 0.00000E+00, 
     C 0.00000E+00, 0.71183E+00, 0.24600E+11, 0.23690E+00, 0.00000E+00,
     C 0.00000E+00, 0.00000E+00, 0.00000E+00,-0.88933E+00, 0.10984E+01,
     C-0.70039E+01, 0.62563E+01, 0.84874E+00, 0.21700E+11, 0.18615E+00,
     C 0.35318E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C 0.86938E+00, 0.00000E+00, 0.10973E+01, 0.22700E+11, 0.18664E+00,
     C-0.24847E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C-0.44347E+00, 0.00000E+00, 0.25134E+01, 0.23700E+11, 0.13887E+00,
     C 0.88100E+00,-0.12226E+01, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     C 0.10069E+01, 0.55355E+01,-0.48144E+01, 0.24700E+11, 0.10555E+00,
     C 0.00000E+00, 0.14526E+00, 0.00000E+00, 0.00000E+00, 0.10303E+01,
     C-0.11688E+01, 0.21800E+11, 0.23129E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.64200E+00, 0.22800E+11,
     C 0.18599E+00/
      DATA PP6/ 0.39711E+00,-0.72595E+00, 0.00000E+00, 0.00000E+00,
     C 0.00000E+00, 0.58167E+00, 0.23800E+11, 0.17148E+00,-0.73124E+00,
     C 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.11635E+01, 0.16493E+01,
     C 0.00000E+00,-0.12551E+02, 0.24800E+11, 0.22062E+00,-0.95906E-01/
      DATA W1/139.65,139.65,938.256/
      DATA W2/938.256,1212.0,547.3/
      DATA DW2/1.0,102.0,0.01/
      DATA NNL,NCH,WI,WT/8,3,139.65,938.256/
      DATA IRM,EM,S11L,NSTRT/-1,0.0,0.0,0/
      DATA WSUB,ETH/150.0,5.0/
      SAVE
      IF(NSTRT.EQ.1) GO TO 1
      DO 51 L=1,8
      DO 51 N=1,4
      NFM(N,L)=0
      TR(N,L)=0.0
      TI(N,L)=0.0
      DO 51 J=1,30
   51 P(J,N,L)=0.0
      NSTRT=1
      DO 54 N=1,18
      NTL(N)=NNTL(N)
      IF(N.GT.13) NTL(N)=NNTX(N-13)
   54 CONTINUE
      I=1
   52 Z=PP(I)/1.E8
      NL=Z+0.1
      NF=NL/100
      NL=NL-100*NF
      N=NL/10
      L=NL-10*N
C     IF(N.EQ.3.AND.L.EQ.1) NF=0
      NFM(N,L)=NF
      J=0
   53 I=I+1
      IF(I.GT.IMX) GO TO 1
      IF(PP(I).GT.1.E8) GO TO 52
      J=J+1
      P(J,N,L)=PP(I)
      IF(J.EQ.1.AND.N.EQ.3.AND.L.EQ.1) WSUB=PP(I)
      GO TO 53
    1 E=EX
      IF(E.LT.ETH) E=ETH
      IF(EM.NE.E) IRM=-27
      IF(TR(2,1).NE.S11L) IRM=-27
      IF(IRR.EQ.IRM) GO TO 98
      IRM=IRR
      IR=IRR
      IF(IR.GT.3) IR=IR-3
      EM=E
      TLB=E
      DO 55 L=1,8
      DO 55 N=1,4
      TR(N,L)=0.0
   55 TI(N,L)=0.0
      IF(IR.LT.0.OR.IR.GT.3) GO TO 99
C SMALLEST ENERGIES SEEM TO BREED TROUBLE
C     IF(TLB.LT.0.2) GO TO 99
      W=SQRT((WI+WT)**2+2.*TLB*WT)
      QPQ=(W**2-1074.7**2)*(W**2-804.6**2)
      QPQ=SQRT(QPQ/(W**2-(WI+WT)**2)/(W**2-(WT-WI)**2))
      DO 7 N=1,NCH
      WSU=W1(N)+W2(N)
      WC=WI+WT+140.0
      IF(N.EQ.1) WC=WSU-WSUB
      GU=-DW2(N)/2.
      IF(N.EQ.1.) GU=0.
      WIM=0.0
   13 CALL CMFN(W,WIM,WC,WSU,GU,NNL,V(1,N),VI(1,N))
      IF(V(1,N).NE.0.0) GO TO 7
      WIM=WIM+1.0
      GO TO 13
    7 CONTINUE
      ETA=IR
      IF(ETA.GT.1.0) ETA=-1.0
      XKMEV=SQRT(TLB*(TLB+2.0*WI)/((1.0+WI/WT)**2+2.0*TLB/WT))
      IF(IR.EQ.3) XKMEV=XKMEV*SQRT(QPQ)
      XKM=197.32/XKMEV
      PZR=SQRT(XKMEV**2+W2(1)**2)
      QZR=SQRT(XKMEV**2+W1(1)**2)
      ETA=ETA*.007297348*(QZR*PZR+XKMEV**2)/XKMEV/(PZR+QZR)
C  PUT COULOMB BARRIER FACTORS INTO VI(L,5) 8/26/82 ARNDT
      IF(ETA.GT.100.) ETA=100.
      Z=2.*3.1415927*ETA
      BL=1.
      IF(ETA.NE.0.) BL=Z/(EXP(Z)-1.)
      Z=0.
      ZZ=SQRT(QPQ)
      DO 8 L=1,NNL
      BPL(L)=BL
      Z=Z+1.
      IF(IR.NE.3) GO TO 8
      BPL(L)=SQRT(BL)*ZZ
      ZZ=ZZ*QPQ
    8 BL=BL*(1.+(ETA/Z)**2)
      DO 9 LL=1,NNL
      DO 9 NN=1,4
      NNF=NFM(NN,LL)
      TRX=0.
      TIX=0.
      NL=10*NN+LL-1
      IF(NN.EQ.1.AND.LL.EQ.1) GO TO 9
      IF(NN.EQ.3.AND.LL.EQ.1) GO TO 9
      IF(NNF.LT.3) CALL TMCM(TLB,NL,NNF,IRR,P(1,NN,LL),V,VI,TRX,TIX)
C  ENCODE C-M K-MTX FIT FOR FORM 4 9/23/81 ARNDT
      IF(NNF.LT.3) GO TO 14
      IF(NNF.LT.3.OR.NNF.GT.6) GO TO 10
      BL=BPL(LL)
      IF(BL.EQ.0.0) BL=1.0
      CER=V(LL,1)
      CEI=VI(LL,1)
      LI=LL
      IF(NN.EQ.1.OR.NN.EQ.3) LI=LL-2
      IF(LI.LT.1) LI=LI+2
      IF(LL.EQ.1) LI=3
      NIL=2
      IF(NIL.GT.NCH) NIL=NCH
      CIR=V(LI,NIL)
      CII=VI(LI,NIL)
      WTH=W1(1)+W2(1)
      WPITH=WTH+140.
      WCM=SQRT(WTH**2+2.*W2(1)*TLB)
      Z=(WCM-WPITH)/1000.
      ZZ=1.
      ZE=0.0
      WKP=P(5,NN,LL)
      DRL=1.0
C MASS-SPLIT K-MTX POLE PIECE FOR P33 1/95 ARNDT
      IF(WKP.EQ.0.0) GO TO 34
      DWK=P(18,NN,LL)/2.0
      IF(NL.NE.41) DWK=0.0
      IF(BL.GT.1.0001) WKP=WKP+DWK
      IF(BL.LT.0.9999) WKP=WKP-DWK
      ZZ=WKP-WTH
      DRL=WKP-WCM
      DGK=P(17,NN,LL)/2.0
      IF(DWK.EQ.0.0) DGK=0.0
      IF(DGK.EQ.0.0) GO TO 34
      IF(BL.GT.1.0001) ZE=DGK
      IF(BL.LT.0.9999) ZE=-DGK
   34 CONTINUE
      IF(NNF.GT.2) Z=(WCM-WTH)/1000.0
      LIP=LI+2
      IF(LIP.GT.8) LIP=LIP-2
      DO 12 J=1,4
      ZE=ZE+P(J,NN,LL)*ZZ
   12 ZZ=ZZ*Z
      ZEI=0.
      DIM=0.
      DO 31 J=1,3
      IF(J.NE.3) GO TO 33
      CIR=V(LIP,2)
      CII=VI(LIP,2)
      IF(NNF.EQ.3.OR.NNF.EQ.5) GO TO 33
      CIR=V(LL,3)
      CII=VI(LL,3)
   33 CONTINUE
C     IF(CII.LT.0.0) CII=0.0
      K=2+4*J
      Z0=Z*P(K,NN,LL)+Z**2*P(K+1,NN,LL)
      IF(Z0.EQ.0.0) GO TO 31
      ZZ=P(K+2,NN,LL)+Z*P(K+3,NN,LL)
      IF(NNF.GT.4) ZZ=ZZ*Z
      DIR=1.0-CIR*ZZ
      DII=-CII*ZZ
      Z2R=Z0**2*CIR
      Z2I=CII*Z0**2
      ZZ=Z2R
      Z2R=ZZ*DRL-Z2I*DIM
      Z2I=ZZ*DIM+Z2I*DRL
      ZZ=ZE
      ZE=ZE*DIR-ZEI*DII+Z2R
      ZEI=ZZ*DII+ZEI*DIR+Z2I
      ZZ=DRL
      DRL=DRL*DIR-DIM*DII
      DIM=ZZ*DII+DIM*DIR
   31 CONTINUE
      DRL=DRL-CER*ZE+CEI*ZEI
      DIM=DIM-CER*ZEI-CEI*ZE
      D2=DRL**2+DIM**2
      Z=CEI/D2
      TRX=Z*(ZE*DRL+ZEI*DIM)
      TIX=Z*(ZEI*DRL-ZE*DIM)
   14 IF(IRR.GT.3) BL=1.0
      IF(BL.EQ.1.0) GO TO 11
      IF(BL.GT.1.0.AND.NL.EQ.41) CALL ETA33(TLB,TRX,TIX)
      CALL PWCC(TLB,TRX,TIX,NL,BL,NFM(3,1),IRR)
   11 IF(TIX.GE.1.0) TRX=0.0
      IF(TIX.GT.1.0) TIX=1.0
      IF(NNF.NE.1) CALL ADDRES(E,TRX,TIX,NL,P(19,NN,LL))
      IF(TIX.GE.TRX**2+TIX**2) GO TO 10
      D2=1.0+TRX**2+TIX**2-2.0*TIX
      IF(D2.LT.1.E-20) WRITE(7,224) TLB,WCM,TRX,TIX
  224 FORMAT(' TLB, WCM=',2F8.2,' TR,TI=',2F9.5)
      IF(D2.LT.1.E-20) GO TO 10
      Z=TRX/D2
      TRX=Z/(1.+Z**2)
      TIX=Z*TRX
   10 TR(NN,LL)=TRX
    9 TI(NN,LL)=TIX
      S11L=TR(2,1)
      IF(WI.GT.150.0) GO TO 98
      IF(NFM(3,1).NE.4) GO TO 98
C add in f13 corrections for S11, P13 ONLY for TROMBERG
      L=1
      CALL TROMF13(TLB,L,IR,TR(2,L),TI(2,L),TR(4,L),TI(4,L))
      L=2
      CALL TROMF13(TLB,L,IR,TR(2,L),TI(2,L),TR(4,L),TI(4,L))
      S11L=TR(2,1)
   98 Z=1.0
      IF(EX.LT.ETH) Z=EX/ETH
      ZZ=SQRT(Z)
      DO 97 L=1,8
      IF(L.GT.1) ZZ=ZZ*Z
      DO 97 N=1,4
      TRZ(N,L)=TR(N,L)
      TIZ(N,L)=TI(N,L)
      IF(Z.EQ.1.0) GO TO 97
      TRZ(N,L)=ZZ*TR(N,L)
      TIZ(N,L)=TRZ(N,L)**2
   97 CONTINUE
   99 RETURN
      END
C *********************************************************
      SUBROUTINE TMCM(TLB,NL,NF,IRR,P,V,VI,TRX,TIX)
C Coupled-Channel CM-K-mtx for FORM=1 9/3/01 RAA
      COMMON/PGLOB/PG(20,4,8),NFG(4,8)
      DIMENSION V(8,3),VI(8,3),P(30)
      DIMENSION TRL(10),TIM(10),ARC(10),AIC(10),RR(10),RI(10)
     C,CCR(4),CCI(4),RH(4),VRR(9),VRI(9)
      DATA WI,WT/139.65,938.256/
      DATA WCMM/0.0/
      SAVE
      NN=NL/10
      LL=NL-10*NN+1
      IRX=IRR
      IF(IRR.GT.3) IRX=IRR-3
C ??? don't know WHAT this is
      WSE=WI+WT
      WRL=SQRT(WSE**2+2.0*WT*TLB)
C get CM functions for rhoN (NF=2)
      IF(WRL.NE.WCMM) CALL CMFN(WRL,-1.0,WSE+WI,1709.3,-75.0,8,VRR,VRI)
      WCMM=WRL
C do 4x4 K-matrix to channels pipi, pid, pieta(or pid+)
C P= K11(4), WkP, K12(2), K22(2), K13(2), K23(2), K33(2)
C K14(2), K24(2), K34(2), K44(2), dWk, dGk, addres(4)
      LI=LL
      IF(NN.EQ.1.OR.NN.EQ.3) LI=LL-2
      IF(LI.LT.1) LI=LI+2
      IF(LL.EQ.1) LI=3
      LIP=LI+2
      ZZR=(WRL-WSE)/1000.0
      Z2R=ZZR**2
      WKP=P(5)
      IF(IRX.EQ.2.OR.IRX.EQ.3) WKP=WKP+P(25)
      WKPR=1.0
      IF(WKP.NE.0.0) WKPR=(WKP-WRL)/1000.0
      NCX=4
      IF(NN.GT.2) NCX=3
      JMX=NCX*(NCX+1)
      JMX=JMX/2
      DO 2 J=2,JMX
      K=2*J+2
      ZR=P(K)*ZZR+P(K+1)*Z2R
      IF(J.LT.7.OR.J.GT.9) GO TO 2
      ZR=P(K)+(WRL/1000.0-1.4)*P(K+1)
    2 ARC(J)=ZR*WKPR
      ZR=1.0
      GE=0.0
      IF(WKP.EQ.0.0) GO TO 9
      ZR=(WKP-WSE)/1000.0
      IF(IRX.EQ.2.OR.IRX.EQ.3) GE=P(24)/1000.0
    9 DO 3 K=1,4
      GE=GE+P(K)*ZR
    3 ZR=ZR*ZZR
      ARC(1)=GE
      IF(ARC(7).EQ.0.0) NCX=3
      IF(ARC(4).EQ.0.0.AND.NCX.EQ.3) NCX=2
      IF(ARC(2).EQ.0.0.AND.NCX.EQ.2) NCX=1
      JD=1
      CRX=V(LL,1)
      CIX=VI(LL,1)
      DO 4 J=1,NCX
      IF(J.EQ.2) CRX=V(LI,2)
      IF(J.EQ.2) CIX=VI(LI,2)
      IF(J.EQ.3) CRX=V(LIP,2)
      IF(J.EQ.3) CIX=VI(LIP,2)
      IF(J.EQ.3.AND.NF.EQ.2) CRX=VRR(LL)
      IF(J.EQ.3.AND.NF.EQ.2) CIX=VRI(LL)
      IF(J.EQ.4) CRX=V(LL,3)
      IF(J.EQ.4) CIX=VI(LL,3)
      CCR(J)=CRX
      CCI(J)=CIX
      DO 5 K=1,J
      TRL(JD)=-ARC(JD)
      TIM(JD)=0.0
    5 JD=JD+1
      D2=CRX**2+CIX**2
      TRL(JD-1)=WKPR*CRX/D2-ARC(JD-1)
      TIM(JD-1)=-WKPR*CIX/D2
    4 CONTINUE
      CALL CMSINV(TRL,TIM,RR,RI,NCX)
      JK=1
      DO 6 J=1,NCX
      JD=J*(J-1)
      JD=JD/2
      RH(J)=CCI(J)
      IF(RH(J).LT.0.0) RH(J)=0.0
      RH(J)=SQRT(RH(J))
      DO 6 K=1,J
      KD=K*(K-1)
      KD=KD/2
      ZR=0.0
      ZI=0.0
      DO 7 M=1,NCX
      MD=M*(M-1)
      MD=MD/2
      JM=JD+M
      IF(M.GT.J) JM=MD+J
      MK=KD+M
      IF(M.GT.K) MK=MD+K
      ZR=ZR+ARC(JM)*RR(MK)
      ZI=ZI+ARC(JM)*RI(MK)
    7 CONTINUE
      DCR=CCR(K)
      DCI=CCI(K)
      D2=DCR**2+DCI**2
      TRX=RH(J)*RH(K)*(DCR*ZR+DCI*ZI)/D2
      TIX=RH(J)*RH(K)*(DCR*ZI-DCI*ZR)/D2
      IF(IRR.LT.5) GO TO 8
C T(pi-eta) = 0 if NO piN coupling
      IF(NCX.EQ.3) TRX=0.0
      IF(NCX.EQ.3) TIX=0.0
      IF(JK.EQ.7) GO TO 99
    6 JK=JK+1
    8 KP1=26
      CALL ADDRES(TLB,TRX,TIX,NL,P(KP1))
C     CALL ADDRESK(NFG(NN,LL),WRL,TRX,TIX,NL,PG(1,NN,LL))
      IF(IRX.GT.1.AND.NL.EQ.41) CALL ETA33(TLB,TRX,TIX)
C KILL BARRIER FACTOR IF HADRONIC IS NEEDED
   99 RETURN
      END
C *********************************************************
      SUBROUTINE CMSINV(AR,AI,AIR,AII,N)
      DIMENSION AR(10),AI(10),AIR(10),AII(10),WR(60),WI(60)
      SAVE
C  INVERT COMPLEX-SYMMETRIC MATRIX
C  MATRICES ARE STORED A11,A12,A22,A13,...  N=ORDER OF MATRIX
C  G=INV OF DIAGONAL ELEMENT (M)  M=SINGULAR ORDER  W=WORKING SPACE
      M=0
      JD=0
      GR=0.0
      GI=0.0
    1 M=M+1
      JD=JD+M
      GR=AR(JD)
      GI=AI(JD)
      IF(M.EQ.1) GO TO 2
      MM=M-1
      CALL MCMCV(AIR,AII,AR(JD-MM),AI(JD-MM),WR,WI,MM)
      JJ=JD-M
      DO 3 K=1,MM
      JJ=JJ+1
      GR=GR-AR(JJ)*WR(K)+AI(JJ)*WI(K)
    3 GI=GI-AR(JJ)*WI(K)-AI(JJ)*WR(K)
    2 D2=GR**2+GI**2
      IF(D2.LT.1.E-9) D2=1.E-9
C Note!! This is to take care of SINGULAR K-mtx  9/12/01 RAA
      GR=GR/D2
      GI=-GI/D2
      AIR(JD)=GR
      AII(JD)=GI
      IF(M.EQ.1) GO TO 5
      J0=1
      DO 4 J=1,MM
      ZR=GR*WR(J)-GI*WI(J)
      ZI=GR*WI(J)+GI*WR(J)
      KK=JD-M+J
      AIR(KK)=-ZR
      AII(KK)=-ZI
      DO 4 K=1,J
      ZZR=ZR*WR(K)-ZI*WI(K)
      ZZI=ZR*WI(K)+ZI*WR(K)
      AIR(J0)=AIR(J0)+ZZR
      AII(J0)=AII(J0)+ZZI
    4 J0=J0+1
    5 IF(M.LT.N) GO TO 1
      RETURN
      END
C *********************************
      SUBROUTINE MCMCV(AR,AI,VR,VI,PR,PI,N)
C MULTIPLY COMPLEX MATRIX(A) ON COMPLEX VECTOR(V) TO GET PRODUCT(P)
      DIMENSION AR(20),AI(20),VR(20),VI(20),PR(20),PI(20)
      J0=0
      DO 1 J=1,N
      ZR=0.0
      ZI=0.0
      DO 2 I=1,J
      IJ=J0+I
      ZR=ZR+AR(IJ)*VR(I)-AI(IJ)*VI(I)
    2 ZI=ZI+AR(IJ)*VI(I)+AI(IJ)*VR(I)
      IF(J.GE.N) GO TO 3
      JP=J+1
      DO 4 I=JP,N
      IJ=I*(I-1)
      IJ=IJ/2+J
      ZR=ZR+AR(IJ)*VR(I)-AI(IJ)*VI(I)
    4 ZI=ZI+AR(IJ)*VI(I)+AI(IJ)*VR(I)
    3 PR(J)=ZR
      PI(J)=ZI
    1 J0=J0+J
      RETURN
      END
C *******************************************
      SUBROUTINE PWCC(TLB,TRX,TIX,NL,BL,NCC,IRZ)
      DATA WI,PI/139.65,3.1415927/
C NCC=5(NO CC), 4(Nordita S,P waves Tl<500), 6(Barrier+"h")
C NCC=7(Nordita+Gibbs to 985, Barrier factors above)
C otherwise use "Barrier" multiplication of K(Hadronic)
      IF(TLB.LT.0.5) GO TO 99
      IF(BL.EQ.1.) GO TO 99
      IF(NCC.EQ.5) GO TO 99
      IF(NCC.LT.6) GO TO 2
      IF(NCC.EQ.7) GO TO 2
C Try adding "H" correction to Eff-Rng Charge-corrections 5/8/00
      T2=TRX**2+TIX**2
      D2=1.0+T2-2.0*TIX
      ZHR=TRX/D2
      ZHI=(TIX-T2)/D2
      B=SQRT(TLB*(TLB+2.0*WI))/(TLB+WI)
      E=1.0/B/137.06
      IF(BL.GT.1.0) E=-E
      h=0.0
      XX=0.0
      Z=1.0
      DO 3 J=1,10
      XX=XX+1.0/Z/(1.0+(E/Z)**2)
    3 Z=Z+1.0
      Z=2.0*PI*E
      IF(Z.GT.80.0) Z=80.0
      C2=Z/(EXP(Z)-1.0)
      H=2.0*E*(E**2*XX-0.57721-ALOG(ABS(E)))*BL/C2
      IF(NCC.EQ.7) H=0.0
      ZHR=ZHR
      ZHI=ZHI
      DR=1.0-ZHR*H
      DI=-ZHI*H
      D2=DR**2+DI**2
      ZCR=BL*(ZHR*DR+ZHI*DI)/D2
      ZCI=BL*(ZHI*DR-ZHR*DI)/D2
      Z2=ZCR**2+ZCI**2
      D=1.0+Z2+2.0*ZCI
      TRX=ZCR/D
      TIX=(ZCI+Z2)/D
      GO TO 99
    2 CONTINUE
      IF(NCC.NE.4.AND.NCC.NE.2.AND.NCC.NE.7) GO TO 1
C Nordita corrections to S,P waves for Tl<550, otherwise Barrier
      TCT=535.0
      IF(NCC.EQ.7) TCT=985.0
      IF(TLB.GT.TCT) GO TO 1
      NN=NL/10
      LL=NL-10*NN+1
      IF(LL.GT.2.OR.NL.EQ.11) GO TO 1
      IF(LL.GT.NCC/2) GO TO 1
C Nordita for S-waves(NNC=2), or S+P-waves(NCC=4)
C Use Barrier factors for all but S-waves MP 5/16/00 Nuts!!
      DD=0.0174532*DTROMB(NL,IRZ,TLB)
      S=SIN(DD)
      C=COS(DD)
      TCR=S*C
      TCI=S**2
      SR=1.0-2.0*TIX
      SI=2.0*TRX
      TRX=TRX+SR*TCR-SI*TCI
      TIX=TIX+SR*TCI+SI*TCR 
      GO TO 99
    1 BLX=BL
      IF(NCC.EQ.7) BLX=BLGIBBS(TLB,NL,BL)
      DR=1.-TIX*(1.-BLX)
      DI=TRX*(1.-BLX)
      D2=(DR**2+DI**2)/BLX
      Z=TRX
      TRX=(Z*DR+TIX*DI)/D2
      TIX=(TIX*DR-Z*DI)/D2
   99 RETURN
      END
C ***************************************
      FUNCTION DTROMB(NL,IS,T)
C     DO QUADRATIC TABLE LOOKUP OF TROMBERG PHASES 
C     K=1(S31),2(P31),3(P33)               Pi+P
C     4(S11),5(S31),6(P31),7(P13),8(P33)   Pi-P/CXS
C     n.b.!! corrections adjusted  so  corr =  del_nuc - del_had
C     (Tromborg had -1/3, -2/3 factors for I=3,1  pi- corrections)
C     TI = 10(25)535 MEV
C     modified June 16/00 by M.M. Pavan
c              Aug  21/01 by MMP  include P31- (Helv.Phys.Acta51,584,1978) 
c              Oct  09/01 by MMP  extend to 985 MeV, blend with Gibbs barrier
c                                 and reduce 'Coulomb only' part of Nordita
c                                 P33 15% for missing 'sigma' exchange 
c                                 (see Bugg NPB58 ('73) p397 Table 2)
      DIMENSION F(320)
      DATA F/0.110, 0.093, 0.091, 0.100, 0.100, 0.110, 0.121, 0.120 
     +     , 0.131, 0.130, 0.130, 0.130, 0.130, 0.132, 0.136, 0.139
     +     , 0.141, 0.142, 0.142, 0.142, 0.141, 0.140, 0.139, 0.139
     +     , 0.133, 0.125, 0.115, 0.109, 0.109, 0.127, 0.148, 0.156
     +     , 0.151, 0.149, 0.147, 0.148, 0.143, 0.145, 0.141, 0.143
     +     , 0.010, 0.012, 0.024, 0.040, 0.049, 0.068, 0.074, 0.090
     +     , 0.092, 0.105, 0.116, 0.127, 0.137, 0.146, 0.155, 0.164
     +     , 0.172, 0.180, 0.189, 0.197, 0.206, 0.216, 0.226, 0.237
     +     , 0.243, 0.249, 0.249, 0.257, 0.259, 0.260, 0.261, 0.260 
     +     , 0.264, 0.263, 0.261, 0.258, 0.258, 0.257, 0.254, 0.255
     +     ,-0.033,-0.105,-0.240,-0.439,-0.745,-1.092,-1.233,-0.949
     +     ,-0.522,-0.200, 0.006, 0.131, 0.197, 0.243, 0.261, 0.268
     +     , 0.275, 0.272, 0.263, 0.251, 0.238, 0.229, 0.220, 0.210
     +     , 0.199, 0.188, 0.175, 0.163, 0.150, 0.137, 0.123, 0.110
     +     , 0.098, 0.085, 0.074, 0.063, 0.053, 0.044, 0.036, 0.029
     +     , 0.238, 0.177, 0.129, 0.096, 0.069, 0.042, 0.016, 0.000
     +     ,-0.010,-0.024,-0.031,-0.038,-0.044,-0.049,-0.053,-0.052
     +     ,-0.048,-0.038,-0.022,-0.010, 0.008, 0.031, 0.059, 0.081
     +     , 0.080, 0.018,-0.034,-0.021, 0.009, 0.040, 0.075, 0.098
     +     , 0.119, 0.108, 0.061, 0.006,-0.041,-0.069,-0.082,-0.086
     +     ,-0.206,-0.145,-0.111,-0.098,-0.090,-0.084,-0.082,-0.080
     +     ,-0.076,-0.073,-0.071,-0.071,-0.072,-0.073,-0.075,-0.078
     +     ,-0.081,-0.085,-0.089,-0.093,-0.097,-0.100,-0.104,-0.107
     +     ,-0.103,-0.099,-0.089,-0.081,-0.086,-0.107,-0.135,-0.136
     +     ,-0.142,-0.142,-0.136,-0.140,-0.137,-0.136,-0.139,-0.142
     +     ,-0.022,-0.051,-0.076,-0.096,-0.112,-0.126,-0.135,-0.144
     +     ,-0.152,-0.161,-0.171,-0.180,-0.191,-0.199,-0.204,-0.209
     +     ,-0.214,-0.220,-0.226,-0.232,-0.237,-0.241,-0.244,-0.249
     +     ,-0.253,-0.258,-0.258,-0.264,-0.266,-0.266,-0.271,-0.268
     +     ,-0.264,-0.265,-0.264,-0.261,-0.262,-0.263,-0.255,-0.256
     +     ,-0.007,-0.021,-0.038,-0.057,-0.072,-0.092,-0.103,-0.114
     +     ,-0.129,-0.139,-0.151,-0.161,-0.176,-0.195,-0.218,-0.237
     +     ,-0.245,-0.244,-0.235,-0.221,-0.204,-0.187,-0.172,-0.160
     +     ,-0.146,-0.131,-0.115,-0.100,-0.084,-0.070,-0.058,-0.049
     +     ,-0.044,-0.041,-0.042,-0.044,-0.048,-0.053,-0.057,-0.062
     +     , 0.136, 0.336, 0.533, 0.726, 0.904, 0.939, 0.638, 0.134
     +     ,-0.213,-0.370,-0.417,-0.414,-0.384,-0.353,-0.319,-0.284
     +     ,-0.254,-0.227,-0.203,-0.182,-0.162,-0.145,-0.130,-0.117
     +     ,-0.106,-0.096,-0.087,-0.080,-0.074,-0.069,-0.064,-0.060
     +     ,-0.056,-0.053,-0.050,-0.046,-0.043,-0.039,-0.034,-0.029/
      SAVE
C-----Initialize
      DTROMB = 0.0
      K = 0
C-----Ignore if PW not covered by Tromborg
      IF(IS.LT.1.OR.IS.GT.3) GO TO 99
      IF(T.GT.1000.0) GO TO 99
C-----Select PW 
      IF(NL.EQ.40) K=1
      IF(NL.EQ.31.AND.IS.EQ.1) K=2
      IF(NL.EQ.41) K=3
      IF(NL.EQ.20) K=4
      IF(NL.EQ.40.AND.IS.GT.1) K=5
      IF(NL.EQ.31.AND.IS.GT.1) K=6
      IF(NL.EQ.21.AND.IS.GT.1) K=7
      IF(NL.EQ.41.AND.IS.GT.1) K=8
C-----Ignore if PW not covered by Tromborg
      IF(K.EQ.0) GO TO 99
C-----Find PW and energy TI near energy T
C-----n.b. need 3 points for quadratic interp.
      I = (T-10.0)/25.0
      I = I+2
      IF(I.GT.39) I=39
      TI = 25*I-15
      I  = 40*(K-1)+I
c-----Interpolate to energy T using nearest 3 PWs
      F0 = F(I)
      FM = F(I-1)-F0
      FP = F(I+1)-F0
C     WRITE(*,222) T,I,K,TI,F0,FM,FP
C 222 FORMAT(' T=',F5.1,' I,K,TI=',2I3,F6.1,' F0,FM,FP=',3F7.3)
      ZM = -25.0
      ZP =  25.0
      Z  =  T-TI
      D  =  ZM*ZP*(ZP-ZM)
      A  =  (ZP**2*FM-ZM**2*FP)/D
      B  =  (ZM*FP-ZP*FM)/D
      DTROMB = F0+Z*(A+Z*B)
   99 RETURN
      END
C ***************************************************
      FUNCTION BLGIBBS(TLB,NL,BL)
      DIMENSION P(17)
      DATA P/14.319,-6.877,1.611,6.033,-.4876,3.51,0.1133,2.368
     C,0.2207,1.756,0.208,1.41,0.169,1.183,0.1297,0.9958,0.1147/
      N=NL/10
      L=NL-10*N
      BLGIBBS=BL
      IF(L.LT.0.OR.L.GT.7) GO TO 99
      Z=TLB/1000.0
      K=1
      IF(L.GT.0) K=2+2*L
      ZZ=P(K)+Z*P(K+1)
      IF(L.EQ.0.0) ZZ=ZZ+Z**2*P(3)
      ZZ=ZZ*SQRT(Z)/1000.0    
      IF(BL.GT.1.0) ZZ=-ZZ
      BLGIBBS=BL+ZZ
   99 RETURN
      END
C **************************************************************
      SUBROUTINE TROMF13(TL,L,MM,T1R,T1I,T3R,T3I)
C add corrections to S11 and P13 for Tromberg's f13  5/16/01 RAA
      DIMENSION P(21),I1(6),NI(6)
      DATA P/0.039,0.366,0.4165,0.0,5.112,-44.511,66.85,70.59
     C,-0.1514,-0.211,0.475,9.98,55.499,-54.3
     C,0.0,0.0,-16.82,13482.0,-55386.0,27.6,-36.73/
      DATA I1/1,4,9,12,15,20/
      DATA NI/3,5,3,3,5,2/
      IF(TL.GT.550.0) GO TO 99
      IF(L.GT.2) GO TO 99
      IF(MM.LT.2) GO TO 99
      ID=1
      IF(L.GT.1) ID=2 
      IF(L.GT.1.AND.TL.GT.250.0) ID=3
      NID=NI(ID)
      ID=I1(ID)
      IE=4
      IF(L.GT.1) IE=5 
      IF(L.GT.1.AND.TL.GT.180.0) IE=6
      NIE=NI(IE)
      IE=I1(IE)
      Z=TL/1000.0
      ZZ=1.0
      D13=0.0
      DO 1 I=1,NID
      D13=D13+P(ID+I-1)*ZZ
    1 ZZ=ZZ*Z
      ZZ=1.0
      E13=0.0
      DO 2 I=1,NIE
      E13=E13+P(IE+I-1)*ZZ
    2 ZZ=ZZ*Z
      E13=E13/10000.0
      D13=0.0174532*D13
      FCT=-1.333333
      IF(MM.EQ.3) FCT=-FCT/2.0
      ZR=E13*FCT
      ZI=D13*FCT
      D1=ATAN(T1I/(T1R+1.0E-8))
      D3=ATAN(T3I/(T3R+1.0E-12))
      IF(D3.LT.0.0) D3=D3+3.1415927
      ZZR=COS(D1+D3)
      ZZI=SIN(D1+D3)
      T1R=T1R+ZR*ZZR-ZI*ZZI
      T1I=T1I+ZR*ZZI+ZI*ZZR
   99 RETURN
      END
C ***********************************************************
      SUBROUTINE ADDRES(E,TR,TI,NL,P)                                   
C  add resonance "bump" to partial-wave 8/94 Arndt                      
      DIMENSION P(30)                                                   
      DATA WPI,WN/139.65,938.256/                                       
      SAVE                                                              
      IF(P(2).GT.0.0) GO TO 10                                          
      DRL=1.0                                                           
      DIM=0.0                                                           
      WR=SQRT((WPI+WN)**2+2.0*WN*E)                                     
      WI=0.0                                                            
      CALL ADDRES2(WR,WI,DRL,DIM,TR,TI,NL,P)                            
      GO TO 99                                                          
   10 GT=P(2)                                                           
      GE=P(1)*GT                                                        
      IF(GE.EQ.0.0) GO TO 99                                            
      GI=GT-GE                                                          
      IF(GI.LE.0.0) GI=0.0                                              
      ER=P(3)                                                           
      IF(ER.EQ.0.0) GO TO 99                                            
C     Z=SQRT(2.0*E**2/(E**2+ER**2))                                     
      Z=2.0*E/(ER+E)                                                    
      RE=SQRT(Z)                                                        
      N=NL/10                                                           
      LE=NL-10*N                                                        
      IF(LE.LT.0.OR.LE.GT.8) LE=1                                       
      IF(LE.GT.0) RE=RE*Z**LE                                           
      RI=(E-150.0)/(ER-150.0)                                           
      IF(RI.LT.0.0) RI=0.0                                              
      RI=RE*RI**3                                                       
      GE=GE*RE                                                          
      GI=GI*RI                                                          
      GT=GE+GI                                                          
      Z=ER-E                                                            
      D=Z**2+GT**2                                                      
      ZR=GE*Z/D                                                         
      ZI=GE*GT/D                                                        
      IF(P(4).EQ.0.0) GO TO 1                                           
      SR=1.0-2.0*ZI                                                     
      SI=2.0*ZR                                                         
      Z2=ZR**2+ZI**2                                                    
      ZE=P(4)*Z2                                                        
      D2=1.0+ZE**2                                                      
      ZR=ZR+(SR*ZE-SI*ZE**2)/D2                                         
      ZI=ZI+(SR*ZE**2+SI*ZE)/D2                                         
    1 SR=1.0-2.0*TI                                                     
      SI=2.0*TR                                                         
      TR=TR+SR*ZR-SI*ZI                                                 
      TI=TI+SR*ZI+SI*ZR                                                 
   99 RETURN                                                            
      END                                                               
C ******************************************************                
      SUBROUTINE ADDRES2(WR,WI,DRL,DIM,TR,TI,NL,P)                      
C ADD RESONANCE FOR COMPLEX(WR,WI) ENERGY 10/18/94 ARNDT                
      DIMENSION P(4)                                                    
      DATA WE,WIN/1078.0,1218.0/                                        
      SAVE                                                              
      GT=-P(2)                                                          
      IF(GT.LT.20.0) GO TO 99                                           
      WRES=ABS(P(3))                                                    
      GE=P(1)*GT                                                        
      GI=GT-GE                                                          
      IF(GI.LT.0.0) GI=0.0                                              
      GE=GT-GI                                                          
      N=NL/10                                                           
      LE=NL-10*N                                                        
      D=1.0/(WR**2+WI**2)                                               
      ZR=1.0-D*WE*WR                                                    
      ZIR=1.0-D*WIN*WR                                                  
      Z=1.0-WE/WRES                                                     
      ZZ=1.0-WIN/WRES                                                   
      ZR=ZR/Z                                                           
      ZIR=ZIR/ZZ                                                        
      IF(WR.LT.WIN) ZIR=0.0                                             
      IF(WI.EQ.0.0) GO TO 1                                             
      ZI=D*WE*WI/Z                                                      
      ZII=D*WIN*WI/ZZ                                                   
      RER=ZR                                                            
      REI=ZI                                                            
      CALL SQZ(RER,REI)                                                 
      RIR=ZIR**2-ZII**2                                                 
      RII=ZIR*ZII*2.0                                                   
      IF(LE.LT.1) GO TO 2                                               
      DO 3 L=1,LE                                                       
      Z=RER                                                             
      RER=Z*ZR-REI*ZI                                                   
    3 REI=Z*ZI+REI*ZR                                                   
      GO TO 2                                                           
    1 REI=0.0                                                           
      IF(ZR.LE.0.0) ZR=0.0                                              
      RER=SQRT(ZR)                                                      
      IF(LE.GT.0) RER=RER*ZR**LE                                        
      RIR=ZIR**2                                                        
      RII=0.0                                                           
    2 GER=GE*RER                                                        
      GEI=GE*REI                                                        
      IF(P(3).GT.0.0) GO TO 8                                           
      ZR=(WR-WE)/(WRES-WE)                                              
      ZI=WI/(WRES-WE)                                                   
      Z=GER                                                             
      GER=Z*ZR-GI*ZI                                                    
      GEI=ZR*GEI+Z*ZI                                                   
    8 GIR=GI*RIR                                                        
      GII=GI*RII                                                        
      GTR=GER+GIR                                                       
      GTI=GEI+GII                                                       
      DR=WRES-WR+GTI                                                    
      DI=-WI-GTR                                                        
      Z=DRL                                                             
      DRL=Z*DR-DI*DIM                                                   
      DIM=Z*DI+DR*DIM                                                   
      D2=DR**2+DI**2                                                    
      TRR=(GER*DR+GEI*DI)/D2                                            
      TRI=(GEI*DR-GER*DI)/D2                                            
      IF(P(4).EQ.0.0) GO TO 5                                           
      SR=1.0-2.0*TRI                                                    
      SI=2.0*TRR                                                        
      ZB=P(4)*GER**2/((WRES-WR)**2+GTR**2)                              
      TBR=ZB/(1.0+ZB**2)                                                
      TBI=TBR*ZB                                                        
      TRR=TRR+SR*TBR-SI*TBI                                             
      TRI=TRI+SR*TBI+SI*TBR                                             
    5 SR=1.0-2.0*TI                                                     
      SI=2.0*TR                                                         
      TR=TR+SR*TRR-SI*TRI                                               
      TI=TI+SR*TRI+SI*TRR                                               
   99 RETURN                                                            
      END                                                               
C ****************************************************                  
      SUBROUTINE ETA33(X,TR,TI)                                         
      DIMENSION P(4)                                                    
      DATA P/70.71,160.1,221.0,0.0307/                                  
C CORRECTS PI-P,CXS P33 FOR N-G CROSS SECTION X=TLAB 11/91 ARNDT        
      BW=P(1)**2/((X-P(2))**2+P(1)**2)                                  
      Z=X**2/(X**2+P(3)**2)                                             
      YS=P(4)*BW*Z                                                      
      ETA=1.0-YS                                                        
      TR=ETA*TR                                                         
      TI=ETA*TI+(1.0-ETA)/2.0                                           
      RETURN                                                            
      END                                                               
C *********************************************************             
      SUBROUTINE CMFN(WR,WI,WZ,WTR,WTI,LMX,CR,CI)                       
C  CHEW-MANDELSTAM FUNCTIONS 7/17/80 ARNDT                              
C  INT(0,1) OF X**(L+1/2)/PI/(X-Z)                                      
C  Z=(W-WT)/(W-WZ)                                                      
      DIMENSION CR(20),CI(20)                                           
      DATA PI/3.1415927/                                                
      DATA ZC/2./                                                       
      IF(LMX.GT.10) LMX=10                                              
      DO 10 L=1,LMX                                                     
      CR(L)=0.                                                          
   10 CI(L)=0.                                                          
      DR=WR-ABS(WZ)                                                     
      D2=DR**2+WI**2                                                    
      IF(D2.LT.1.) GO TO 99                                             
      ZR=((WR-WTR)*DR+WI*(WI-WTI))/D2                                   
      ZI=(DR*(WI-WTI)-WI*(WR-WTR))/D2                                   
      AR=ZR                                                             
      AI=ZI                                                             
      CALL SQZ(AR,AI)                                                   
      IF(WZ.LT.0..AND.ZI.LT.0.) AR=-AR                                  
      IF(WZ.LT.0..AND.ZI.LT.0.) AI=-AI                                  
      Z2=ZR**2+ZI**2                                                    
      IF(Z2.LT.ZC**2) GO TO 11                                          
C  USE POWER SERIES FOR Z GTO INF                                       
      RR=ZR/Z2                                                          
      RI=-ZI/Z2                                                         
      L=0                                                               
   12 TL=2*L                                                            
      TL=(TL+3.)/2.                                                     
      SR=-RR/PI/TL                                                      
      SI=-RI/PI/TL                                                      
      TR=SR                                                             
      TI=SI                                                             
      RT=SQRT(TR**2+TI**2)                                              
      DO 13 N=1,20                                                      
      R=TL/(TL+1.)                                                      
      Z=TR                                                              
      TR=R*(RR*Z-RI*TI)                                                 
      TI=R*(RI*Z+RR*TI)                                                 
      SR=SR+TR                                                          
      SI=SI+TI                                                          
      IF(R*RT.LT.1.E-6) GO TO 14                                        
   13 TL=TL+1.                                                          
   14 L=L+1                                                             
      CR(L)=SR                                                          
      CI(L)=SI                                                          
      IF(L.GE.LMX) GO TO 99                                             
      GO TO 12                                                          
   11 A2=AR**2+AI**2                                                    
      D=1.+A2+2.*AR                                                     
      ZZR=(1.-A2)/D                                                     
      ZZI=-2.*AI/D                                                      
      CALL ALG(ZZR,ZZI)                                                 
      BR=2./PI-AI+(AR*ZZR-AI*ZZI)/PI                                    
      BI=AR+(AR*ZZI+AI*ZZR)/PI                                          
      ZL=.5                                                             
      L=0                                                               
    1 L=L+1                                                             
      CR(L)=BR                                                          
      CI(L)=BI                                                          
      IF(L.GE.LMX) GO TO 99                                             
      ZL=ZL+1.                                                          
      ZZ=BR                                                             
      BR=ZR*BR-ZI*BI+1./PI/ZL                                           
      BI=ZI*ZZ+ZR*BI                                                    
      GO TO 1                                                           
   99 RETURN                                                            
      END                                                               
C **********************************************************            
      SUBROUTINE ALG(ZR,ZI)                                             
C  TAKE NATURAL LOG OF Z  BRANCH CUT AT Z=0 TAKEN TO LEFT               
      DATA PI/3.1415927/                                                
      IF(ZR.EQ.0.) ZR=1.E-10                                            
      ZM=ZR**2+ZI**2                                                    
      PHI=ATAN(ZI/ZR)                                                   
      IF(ZR.GT.0.) GO TO 1                                              
      IF(ZI.GE.0.) PHI=PHI+PI                                           
      IF(ZI.LT.0.) PHI=PHI-PI                                           
    1 ZR=ALOG(ZM)/2.                                                    
      ZI=PHI                                                            
      RETURN                                                            
      END                                                               
C **************************************************************        
      SUBROUTINE SQZ(ZR,ZI)                                             
C  SQRT(Z)  BRANCH CUT TAKEN TO LEFT OF Z=0                             
      DATA PI/3.1415927/                                                
      ZM=SQRT(SQRT(ZR**2+ZI**2))                                        
      IF(ZR.EQ.0.) ZR=1.E-10                                            
      PHI=ATAN(ZI/ZR)                                                   
      IF(ZR.GT.0.) GO TO 1                                              
      IF(ZI.LT.0.) PHI=PHI-PI                                           
      IF(ZI.GT.0.) PHI=PHI+PI                                           
    1 PHI=PHI/2.                                                        
      ZR=ZM*COS(PHI)                                                    
      ZI=ZM*SIN(PHI)                                                    
      RETURN                                                            
      END                                                               
C **************************************************************        
      SUBROUTINE COUL(QR,TL,WI,WT,A,FR,FI,GR)                           
C  COULOMB AMPLITUDES FOR PI-N FROM HOHLER HANDBOOK 1979                
C  QR=REL CHARGE  TL=TLAB  A=THETA(CM)  WI=BEAM MASS  WT=TARGET MASS    
C  FR,FI=NON-FLIP AMPLITUDES  GR=SPIN FLIP AMPLITUDE                    
C  USES NATURAL UNITS (SCALED BY CHARGED PION MASS)                     
      DIMENSION TF(5),F1F(9),F2F(9),FPF(4)                              
      DATA TF/91.,130.,31.468,70.,140./                                 
      DATA F1F/.147,-.16934,.92,-.352,-.068,.47737,.04497,18.2,27.523/  
      DATA F2F/-2.2,.90454,-.072,-.008,.02,2.6676,.48086,13.758,30.943/ 
      DATA FPF/1.0538,35.093,.0538,76.8/                                
      DATA WPI/139.57/                                                  
      SAVE                                                              
      FR=0.                                                             
      FI=0.                                                             
      GR=0.                                                             
      IF(A.LT.1.E-3) GO TO 99                                           
      IF(TL.EQ.0.) GO TO 99                                             
      S=(WI+WT)**2+2.*WT*ABS(TL)                                        
      Q2=(S-(WI+WT)**2)*(S-(WI-WT)**2)/4./S/WPI**2                      
      Z=COS(.0174532*A)                                                 
      W=SQRT(S)/WPI                                                     
      WM=WT/WPI                                                         
      T=-2.*Q2*(1.-Z)                                                   
      E=SQRT(Q2+WM**2)+WM                                               
    1 FP=FPF(1)/(1.-T/FPF(2))-FPF(3)/(1.-T/FPF(4))                      
      F1=(F1F(6)+F1F(7)/(1.-T/F1F(8))**2)/(1.-T/F1F(9))                 
      F2=(F2F(6)+F2F(7)/(1.-T/F2F(8)))/(1.-T/F2F(9))                    
      DO 2 I=1,5                                                        
      ZZ=1.-T/TF(I)                                                     
      F1=F1+F1F(I)/ZZ                                                   
    2 F2=F2+F2F(I)/ZZ                                                   
      IF(WI.GT.150.) FP=1.                                              
      Q=SQRT(Q2)                                                        
      ALF=QR/137.036                                                    
      G=ALF*(TL+WI)/SQRT(TL*(TL+2.*WI))                                 
      FR=2.*Q*G/T+ALF/2./W*(W+WM)/E                                     
      FR=FR*F1+(W-WM+T/4./E)*F2*ALF/2./W/WM                             
      FR=FR*FP                                                          
      GR=(W+WM)/E*F1+(W+T/4./E)*F2/WM                                   
      GR=GR*ALF*FP/2./W*SQRT((1.+Z)/(1.-Z))                             
      PHI=-G*ALOG((1.-Z)/2.)                                            
      IF(WI.GT.150.) GO TO 4                                            
      IF(Q2.GT.4.) GO TO 3                                              
C THE INTEGRAL OVER FP, F1 IS FITTED TO AN ACCURACY OF .01% FOR 0<Q2<4  
C AND .6% FOR 4<Q2<25 (IN PION MASS UNITS) JOHN FORD 8/9/85             
      ZZ=Q2/25.                                                         
      PHI=PHI-G*ZZ*7.72829*EXP(ZZ*(-2.8305+ZZ*(4.79697-ZZ*3.1306)))     
      GO TO 4                                                           
    3 ZZ=ALOG(Q2)                                                       
      PHI=PHI-G*Q2*0.270088*EXP(-ZZ*(0.0313754+0.0895484*ZZ))           
    4 ZZ=197.32/WPI                                                     
      FR=FR*ZZ                                                          
      GR=GR*ZZ                                                          
      FI=FR*SIN(PHI)                                                    
      FR=FR*COS(PHI)                                                    
   99 RETURN                                                            
      END                                                               
C *************************************************************         
      FUNCTION SGLSCL(E,NL,WI)                                          
C Scale Coulomb rotation phases for Pi-N interaction 1/90 Arndt         
C This multiplies Point-Coulomb rotation to give Sigma's of             
C Tromberg et al; PRD15, 725(1977)                                      
      DIMENSION FF(6,7),DUM(42),FK(2,7),DK(14)                          
      EQUIVALENCE (DUM,FF(1,1)),(FK(1,1),DK)                            
      DATA EKX,EKY/200.0,1800.0/                                        
      DATA DK/0.54371, -0.00915,  0.27221, -0.00208,  0.22038, -0.00532,
     C  0.17635, -0.00485,  0.16644, -0.00632,  0.14591, -0.00563,      
     C  0.14045, -0.00602/                                              
      DATA EX,EY/200.0,1500.0/                                          
      DATA DUM/1.00130, 0.16663, 0.34645, 0.02116, 0.60983, 0.03425     
     C,  0.36022, 0.03279, 0.24906, 0.01096, 0.29734, 0.02064           
     C,  0.25217, 0.01239, 0.20328, 0.00508, 0.23804, 0.00746           
     C,  0.20849, 0.00462, 0.17804, 0.00173, 0.19015, 0.00399           
     C,  0.18365, 0.00097, 0.16188,-0.00025, 0.17895,-0.00025           
     C,  0.16718,-0.00096, 0.15047,-0.00147, 0.15690,-0.00069           
     C,  0.15529,-0.00207, 0.13828,-0.00122, 0.15097,-0.00207/          
      SAVE                                                              
      N=NL/10                                                           
      L=NL-10*N                                                         
      Z=1.0                                                             
      IF(L.EQ.0) GO TO 99                                               
      IF(L.GT.7) L=7                                                    
      J=5                                                               
      IF(N.EQ.1.OR.N.EQ.3) J=1                                          
      IF(N.EQ.2.OR.N.EQ.4) J=3                                          
      EB=E/EX                                                           
      A=FF(J,L)                                                         
      B=FF(J+1,L)                                                       
      IF(WI.LT.200.0) GO TO 1                                           
      A=FK(1,L)                                                         
      B=FK(2,L)                                                         
    1 Z=1.0/(1.0+EB*A+EB*(EB-1.0)*B)                                    
   99 SGLSCL=Z                                                          
      RETURN                                                            
      END                                                               
C ************************************************************
      SUBROUTINE DOFT(TR,TI,DLM)
C  GET PHASE AND 1-ETA**2 FROM TR,TI  10/30/79 ARNDT
      SR=1.-2.*TI
      SI=2.*TR
      CALL GTPHS(SR,SI,DLM,TR)
      TI=1.-SR**2-SI**2
      RETURN
      END
C **************************************************************
      SUBROUTINE GTPHS(SR,SI,DM,Z)
      SAVE
      Z=0.
      IF(SR.EQ.0.) GO TO 99
      Z=ATAN(SI/SR)*90./3.1415927
      SGN=SI/ABS(SI+1.E-8)
      IF(SR.GE.0.) GO TO 1
      Z=Z+90.*SGN
    1 IF(Z.EQ.DM) GO TO 99
      S=-180.
      IF(Z.LT.DM) S=-S
      Z=Z-DM
    3 ZP=Z+S
      IF(ZP**2.GT.Z**2) GO TO 2
      Z=ZP
      GO TO 3
    2 Z=Z+DM
   99 RETURN
      END
