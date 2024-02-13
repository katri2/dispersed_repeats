       program SCAN
        include 'mpif.h'
        common/myyy/myid,numprocs,I1000
c     s puti
c     sokratil dinamo - virav v  trube!!!
c   umensil razmer v dinamo
       CHARACTER*70 NAME(10)
       INTEGER SEQQ2(10000)
       INTEGER*1 SEQ5(200000000)
       COMMON /CAP/IPER,KAPPA
       REAL COR(4,4,5000),COR1(4,5000)
       INTEGER SEQ1(10000)
       COMMON/KORELL/COR,COR1
       COMMON/SLY/ISEED
       REAL SPE(100000)
       COMMON/PUT2/S22,IM,JM,I0,J0
      CHARACTER*1 TOT
      character*3 WWW1,WW1
      CHARACTER*4 WWW2,WW2
      character*7 WWW3,WW3
c----------------------------------
      INTEGER V1(10000),V2(10000)
      INTEGER VV1(10000),VV2(10000)
      COMMON/VIRAV/V1,V2,K10,SINF,M
      INTEGER M(20,800)
C-------------------------------- 
      COMMON/LET/AAA
      CHARACTER*1 AAA(0:30),AAA1(0:1000)
      COMMON/LET1/AAA1
C      
      COMMON/SAGSAG/ISAG
      REAL ZZ(4000000)
      INTEGER G1(4000000),G2(4000000),JJ(4000000)
c
      CHARACTER*3 CHR1,CHR2
      CHARACTER*6 FASTA,CHR3
      CHARACTER*13 CHR4
      CHARACTER*2 ZET,PPP
      COMMON/LEV/SLEV
C
      INTEGER MAT(4,4,5000)
      COMMON/MATMAT/MAT
      real MC(4,4,5000),M1(4,5000)
      character*5 FFF
C--------------------------------------------------------------
      REAL DD(0:1000) 
      COMMON/STOLB/DD 
C--------------------------------------------------------------
      integer myid,numprocs,ierr,rc
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD,myid,ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD,numprocs,ierr )
C                                          call LOADSETTINGS
C--------------------------------- 
C
c          OPEN(17,FILE='epdat_600.fasta',STATUS='old')
c          OPEN(18,FILE='scan_loc1.txt',STATUS='replace')
c          OPEN(19,FILE='600.txN',STATUS='old')
c       
c          OPEN(18,FILE='chr1_400_35rt1.txt',STATUS='replace')
       OPEN(17,FILE='genome.fasta',STATUS='old')
c       OPEN(17,FILE='ecoli_k12.fasta',STATUS='old')
c      OPEN(17,FILE='ecoli_5845.fasta',STATUS='old')
       OPEN(20,FILE='nom.txt',STATUS='old')
       OPEN(22,FILE='iseed.txt',STATUS='old')
C          OPEN(19,FILE='600.txN',STATUS='old')
c          OPEN(20,FILE='seq.txt',STATUS='new')
C---------------------------------------------
      read(22,*) ISEED
      PRINT 784, ISEED
784   FORMAT('DO_ISEED=',I40)
C                  np=4
          numprocs=4
c          md=0
          md=myid 
c---------------------------------------------
c          ISEED=4111
          id=0
          TOT='.'
      read (20,11) NOM,IPER
      PRINT 11, NOM,IPER
11    FORMAT(4X,I6,7X,I6)
      CLOSE (20)        
c
         write(WWW1,430) md+501
         write(WW1,430)  id+200
         WW2=TOT//WW1     
         WW3=WWW1//WW2
         OPEN(18,FILE=WW3,STATUS='replace')
c----------- prokrutka v konez
550          READ(18,551,END=552) FFF
551          FORMAT(5X,A5)
             go to 550
552        CONTINUE
           BACKSPACE 18
c
c**********************
c DLINA MATRIZI  !!!!!!!!!!!
       KOKO=1
C       IPER=400
c**********************
       KON=0
       KAPPA=4
       I100=0
c                          ISEED=5628
       NNN=0
       K7=0
c                          SMAX=6.0
       KOLS=0.0
       MAT=0       
C       ISAG=10       
       ISAG=IPER
C
c     
      call ALF
C         CALL READCOR
C         KAPPA=4
c  
30     CALL EEE_FASTA(IDL1,NAME,SEQ5,KON)
       IF(KON.NE.0) GO TO 250
c             call INVER(SEQ5,IDL1)
              WRITE(18,238) NAME(1),NAME(2),IDL1
c             WRITE(18,233)(SEQ5(I),I=1,IDL1)
233         FORMAT(60I1)
C  SLUTCHAINOST !!!!!!!!!!!!!!
c           L1=1
c           L2=IDL1
c           call SLY_POSL1(SEQ5,L1,L2)
c           call SLY_TRI(SEQ5,L2)
c            
C           WRITE(18,233)(SEQ5(I),I=L1,L2)
C233        FORMAT(60I1)
C           STOP           
C           
300    CONTINUE
C  ЭТО УРОВЕНЬ ФОНА ДЛЯ ПОИСКА МАКСИМУМОВ. НИЖЕ ЕГО ДАННЫХ НЕТ.
       SLEV=3.0
C  ЭТО УРОВЕНЬ ЗАПИСИ. нИЖЕ ЕГО В 500 ФАЙЛЫ НЕТ ЗАПИСИ.
       OPEN(25,FILE='szap.txt',STATUS='old')
       read(25,499)SZAP
499    format(f9.2)
       PRINT 498, SZAP      
498    FORMAT('SZAP=',F9.2)
       CLOSE(25)
C       SZAP=4.0  
C            IF(KOKO.EQ.1)SLEV=3.0
       I100=0
       K7=0
       NNN1=0
       ZZZ1=0.0
       KOL1=0
       KOL2=0
       KOL3=0.0
C----------------------------------------------------
C      VXODNIE MATRIZI ZAPOLNENEI******
       ZET='Z='
       MAT=0
c           write(WWW1,430) md+601
c           write(WW1,430)  id+100
430        format(I3)
c           WW2=TOT//WW1     
c           WW3=WWW1//WW2
           OPEN(19,FILE='mat.txt',STATUS='old')
c              READ(19,18) IPER
c              print 18, IPER
c18            FORMAT(5X,I6)
c
           DO 140 I=1,4
           DO 141 J=1,4
             read(19,325)I1,J1
325          FORMAT(9x,I6,4X,I1)
             backspace 19
             read(19,326)(MAT(I1,J1,K),K=1,IPER)
326          FORMAT(22X,1000I8)
141         CONTINUE
140         CONTINUE               
           close (19) 
c
C   COZDANIE KORRELAZIONNIX MATRIZ!!!!!!!!!!!!
      CALL NORMA_C(COR)
      CALL MATRI(MC,M1)
      KAPPA=4 
C      
        DO 840 K=1,IPER
        DO 810 I=1,KAPPA
        DO 820 J=1,KAPPA
           COR(I,J,K)=MC(I,J,K)
820     CONTINUE
810     CONTINUE
840     CONTINUE
C
        DO 850 K=1,IPER
        DO 860 I=1,KAPPA
          COR1(I,K)=M1(I,K)
860     CONTINUE
850     CONTINUE
C---------------------------------
c         write(WWW1,430) md+501
c         write(WW1,430)  id+100
c         WW2=TOT//WW1     
c         WW3=WWW1//WW2
c         OPEN(18,FILE=WW3,STATUS='old')
C-----------------------------------------------------------
               WRITE(18,226)
226            FORMAT('START-----------------------')
               DO 225 I=1,4
                WRITE(18,277)I,(COR1(I,K),K=1,IPER)
277            FORMAT('SINGL=',2X,'I=',2X,i8,2x,'Z=',1000F8.2)
225          CONTINUE
                   DO 172 I=1,4
                   DO 172 J=1,4
172                write(18,171) I,J,(COR(I,J,K),K=1,IPER)
171                FORMAT('COR=',2x,'I=',I6,2X,'J=',I1,2X,1000F8.2)
                   DO 272 I=1,4
                   DO 272 J=1,4
272                write(18,271) I,J,(MAT(I,J,K),K=1,IPER)
271                FORMAT('MCOR=',2x,'I=',I6,2X,'J=',I1,2X,1000I8)
C       MAT=0       
C   
238    FORMAT('>',2X,2A70,2X,I8)       
c       stop
C       
       IDL=IPER
c       
       KK=(1.0*IDL1)/(1.0*numprocs)
       I1=KK*myid+1
       I2=KK*(myid+1)
       WRITE(18,149)I1,I2,KK,IDL1
149    FORMAT('I1=',I12,2X,'I2=',I12,2X,'KK=',I12,2X,'IDL1=',I12)
       CALL SEQUENCE(IDL,SEQ1)
C            300    CONTINUE
C
       DD=0.0
       DO 20 J=I1,I2,ISAG
          IF((I2-J).LT.IDL) GO TO 20
C          
              NNN=0
              DO 21 I=1,IDL
                 IF(SEQ5(I-1+J).EQ.5) NNN=NNN+1
21            SEQQ2(I)=SEQ5(I-1+J)
              IF(NNN.GT.0) GO TO 20
                 K7=K7+1
C!!!!!!!!!!!!!!!!!!!!!!!!!!              
                 IF(K7.GT.1000) GO TO 100
C              print 888, (SEQQ2(I),I=1,IDL)
C888           FORMAT(2X,'SEQQ2=',60I1)              
C                                     IF(NNN.GT.1000) KON=1
C            PRINT 566,K7,J,IDL1
C                                     WRITE(18, 566) NNN
566        FORMAT(2X,'SLY, K7=',I8,2X,'J=',I12,2X,'IDL=',I12)
C                                     IF(KON.EQ.1) GO TO 100
c dlya dispersii i srednego!!           
           L1=1
           L2=IDL
           call SLY_POSL(SEQQ2,L1,L2)
C                                     call UNICOD
c                                     CALL SEQUENCE(IDL,SEQ1)
       CALL DINAMO(SEQ1,SEQQ2,IDL)
c                                     IF(SMAX.LT.S22) SMAX=S22
c                                     10     CONTINUE
       I100=I100+1
C       WRITE(18, 311) J,I100,S22
C311    FORMAT('J=',I10,2X,'I100=',I6,2X,'S22=',F8.1) 
       SPE(I100)=S22
20    CONTINUE
c--------------------------------------------------------------------------
100    CONTINUE
C       IF(I100.LT.10) GO TO 250    
C
       CALL SRR(SPE,I100,SRED,DISP)
       KON=0
       I5=0
       WRITE(18,179) SRED,SQRT(DISP)
       PRINT 179, SRED,SQRT(DISP)
179    FORMAT(5X,'SRED=',F8.2,2X,'SIG=',F8.2)       
C
c--------------------------------------------------------------------------
       NNN=0
C         K7=0
       CALL POINTS(I1,I2,IDL,ISAG,SEQ5,SEQQ2,SRED,DISP,J,SEQ1,
     8                   ZZ,G1,G2,JJ,K7,L2)
       IF(K7.GT.0) GO TO 453    
            SLEV=-10.0
            CALL POINTS(I1,I2,IDL,ISAG,SEQ5,SEQQ2,SRED,DISP,J,SEQ1,
     8                   ZZ,G1,G2,JJ,K7,L2)
            SLEV=3.0
453    CONTINUE

       CALL GRADIENT(ZZ,G1,G2,JJ,K7,K8)
C
C    *** NULI V MATRIZE ***  
       CALL NULL
C
       KOKO=0
       IF(K8.LE.0) GO TO 200
       DO 120 J5=1,K8
            J=JJ(J5)
            IF(J.LT.1) GO TO 120
            IF(J.GT.IDL1-IDL) GO TO 120
c                                    120     CALL EEE_FASTA(IDL,NAME,SEQQ2,KON)       
          IF((I2-J).LT.IDL) GO TO 120
              NNN=0
              DO 121 I=1,IDL
                 IF(SEQ5(I-1+J).NE.5) GO TO 122
123                    T1=GSU2R(ISEED)
                       IF(T1.EQ.1.0) GO TO 123    
                    IF(SEQ5(I-1+J).EQ.5) NNN=NNN+1
                    IF(SEQ5(I-1+J).EQ.5) SEQQ2(I)=INT(T1*4.0)+1
                    GO TO 121
122              SEQQ2(I)=SEQ5(I-1+J)
121           CONTINUE                 
c              IF(NNN.GT.50) GO TO 120
              IF(NNN.GT.0) GO TO 120
c
C                                    PRINT 567,NNN
C                                     567         FORMAT(2X,'REAL, NNN=',I6)           
C                                    IF(KON.EQ.1) GO TO 200
C                                    call UNICOD
c                      CALL SEQUENCE(IDL,SEQ1)
C                      DO 10 I=1,10000,100
C                      ISEED=I
       DD=0.0
       CALL DINAMO(SEQ1,SEQQ2,IDL)
       Z=(S22-SRED)/(SQRT(DISP))
C           IF(Z.LT.SLEV) GO TO 120
        IF(Z.GE.7.0) KOL3=KOL3+1
        IF(Z.GE.4.0) KOL2=KOL2+1
C
        IF(Z.LT.SZAP) GO TO 160
        KOKO=KOKO+1  
c                                 IF(Z.le.3.0) KOL1=KOL1+1
c                                 IF(Z.GT.SMAX) KOLS=KOLS+1
       PRINT 734,Z,J,S22,K10,ZZ(J5)
734    FORMAT('Z='F8.2,2X,'J=',I10,2X,'S22=',F8.2,2X,'K10=',I8,2X,
     8        'ZZ(J)=',F8.2)
C
      WRITE(18,733) Z,J,S22,I0,I0-1+J,IM,IM-1+J,K10,J0,JM,ZZ(J5),KOKO
c         WRITE(18,733) Z,J,S22,I0,IDL1-I0+2-J,IM,IDL1-IM+2-J,K10,J0,JM
733    FORMAT('Z='F18.10,2X,'J=',I10,2X,'S22=',F8.2,2X,'I0=',I3,2X,
     8 'J0=',I10,2X,'IM=',I3,2X,'JM=',I10,2X,'K10=',I8,5X,
     8 'J0=',I6,1X,'JM=',I6,2X,'ZZ(J)=',F8.2,2X,I8)
C
      write(18,144) (V2(I),I=1,K10) 
C                       write(18,145) (AAA(V1(I)),I=1,K10)
      write(18,146) (V1(I),I=1,K10)
144   FORMAT(5X,'V1=',1000I3)       
145   FORMAT(5X,'V2=',1000(2X,A1))       
146   FORMAT(5X,'V2=',1000(2X,I1))       
      KOL1=KOL1+1 
C
160   CONTINUE
C
C *** RATCHET MATRIZII ***
         NNN1=NNN1+1
C             IF((NNN1.GT.10).AND.(Z.GT.SLEV)) GO TO 501   
C             IF(NNN1.GT.10) GO TO 110   
c ***  level for matriza!!!!!
C           IF(Z.GT.4.0) GO TO 501   
c  ne menshe 10 slytchaev dlya matrizi !!!!
C           IF(NNN1.GT.10) GO TO 110
C
C501      continue
C             DO 500 I=2,K10
C                IF(V1(I-1).EQ.0) GO TO 500
C                IF(V1(I).EQ.0) GO TO 500
C                IF(V2(I).EQ.0) GO TO 500
C                MAT(V1(I-1),V1(I),V2(I))=MAT(V1(I-1),V1(I),V2(I))+1
C500         CONTINUE
C            ZZZ1=ZZZ1+Z
C            KOL1=KOL1+1
            CALL INF_PER(MAT,SINF)
C     
C     
110    CONTINUE
120   CONTINUE
200   CONTINUE
C
      IF(K8.EQ.0) PRINT 711
711   FORMAT('K8=0!!!!!!!!!!!2_1') 
      IF(K8.GT.0) GO TO 883
             DO 882 K=1,IPER
             DO 881 I=1,KAPPA
             DO 880 J=1,KAPPA
                MAT(I,J,K)=10
880          CONTINUE
881          CONTINUE
882          CONTINUE
883   CONTINUE 
C      CALL NORMA_C(COR)
C      CALL MATRI(MC,M1)
C      KAPPA=4 
C      
C        DO 800 K=1,IPER
C        DO 810 I=1,KAPPA
C        DO 820 J=1,KAPPA
C           COR(I,J,K)=MC(I,J,K)
C820     CONTINUE
C810     CONTINUE
C800     CONTINUE
C
C        DO 850 K=1,IPER
C        DO 860 I=1,KAPPA
C          COR1(I,K)=M1(I,K)
C860     CONTINUE
C850     CONTINUE
C
        WRITE(18,177)KOL1,KOL2,KOL3, NOM,SINF
        PRINT 177,KOL1, KOL2,KOL3,NOM,SINF
177     FORMAT('KOL1_DLY MAT=',I6,2X,'KOL2_>5.0=',I6,2X,'KOL3_>7.0='
     8 ,I6,2X,'NOM=',I6,2X,'SINF=',F8.1)
C           GO TO 30 
C       KOKO=KOKO+1
C       IF(KOKO.GT.100) STOP
C       GO TO 300 
C       
250   CONTINUE       
C---------------------------
      REWIND 22        
      WRITE(22,*) ISEED
      PRINT 783,  ISEED
      CLOSE (22)
783   FORMAT('AFTER_ISEED=',I10)
      close (17)
      close (18)    
      call MPI_FINALIZE(rc)
c---------------------------
       end
c===================================================================
      SUBROUTINE EEE_FASTA(IDL,NAME,SEQ1,KON)
C
      INTEGER*1 SEQ1(200000000)
      CHARACTER*1 PROM(1000),PU,ID,PP
      CHARACTER*70 NAME(10)
C----
      CHARACTER*1 AAA(0:30)
      COMMON/LET/AAA
      COMMON /CAP/IPER,KAPPA
C-----
      ID='>'
      KON=0
      PU=' '
C-----------------------
      SEQ1=99
10    READ(17,1,END=100) PP
C      print 22, PP
22       format(5x,a1)
1        FORMAT(A1)
      IF(PP.NE.ID) GO TO 10
      BACKSPACE 17
      READ(17,2,END=100) NAME(1),NAME(2)
c      WRITE(18,2) NAME(1)
C      print 2, NAME(1)
2     FORMAT(1X,2A70)
C
      I1=1
50    READ(17,3,END=41) (PROM(I),I=1,200)
C      print3,(PROM(I),I=1,80)
3     FORMAT(200A1)
      IF(PROM(1).EQ.ID) GO TO 40
      IKK=1
      DO 30 I=1,200
C                                              SEQ1(I1+I-1)=1
           DO 31 J=1,KAPPA+1
              IF(PROM(I).EQ.AAA(J)) SEQ1(I1+IKK-1)=J
              IF(PROM(I).EQ.AAA(J)) IKK=IKK+1
C                                              IF(PROM(I).EQ.PU) SEQ1(I1+I-1)=1
31         CONTINUE
30    CONTINUE
      I1=I1+IKK-1
      GO TO 50
C----------------
      go to 40
41       CONTINUE
      IDL=I1-1
C                          PRINT 771, KON           
C                          771      FORMAT(2X,'KON========',I8)
C                          DO 160 K=1,10000000
C                          IF(SEQ1(K).EQ.99) IDL=K-1
C                          IF(SEQ1(K).EQ.99) GO TO 161
C                          160        CONTINUE
C                          161        CONTINUE
      GO TO 100
C           
40    BACKSPACE 17
      IDL=I1-1
C                         DO 60 K=1,10000000
C                         IF(SEQ1(K).EQ.99) IDL=K-1
C                         IF(SEQ1(K).EQ.99) GO TO 61
C                         60    CONTINUE
C                         61    CONTINUE
      GO TO 200
100   CONTINUE
         KON=1
200   CONTINUE
c        WRITE(18,122) IDL,KON,NAME(1)
122   FORMAT(2X,'IDL=',2I8,2X,A40)
C        write(18,155) (SEQ1(I),I=1,IDL)
155   FORMAT(60I1)
C        PRINT 771, KON  
771   format('kon=',i8)      
      RETURN
      END
C====================================================
      SUBROUTINE ALF
      COMMON/LET/AAA
      CHARACTER*1 AAA(0:30),AAA1(0:1000)
      COMMON/LET1/AAA1
c
C
      AAA(0)= '*'
      AAA(1)= 'A'   
      AAA(2)= 'T'   
      AAA(3)= 'C'   
      AAA(4)= 'G'   
      AAA(5)= 'N'   
C      
C      AAA(1)= 'K'   
C      AAA(2)= 'N'   
C      AAA(3)= 'I'   
C      AAA(4)= 'M'   
c      AAA(5)= 'T'   
      AAA(6)= 'R'   
      AAA(7)= 'S'   
      AAA(8)= 'L'   
      AAA(9)= 'Y'   
      AAA(10)='F'   
      AAA(11)='C'   
      AAA(12)='W'   
      AAA(13)='P'   
      AAA(14)='H'   
      AAA(15)='Q'   
      AAA(16)='V'   
      AAA(17)='A'   
      AAA(18)='D'   
      AAA(19)='E'   
      AAA(20)='G'   
      AAA(21)='*'
      AAA(22)='.'
c
       AAA1(0)='*'   
       AAA1(1)='1'   
       AAA1(2)='2'   
       AAA1(3)='3'   
       AAA1(4)='4'   
       AAA1(5)='5'
       AAA1(6)='6'
       AAA1(7)='7'
       AAA1(8)='8'
       AAA1(9)='9'
       AAA1(10)='a'
       AAA1(11)='b'
       AAA1(12)='c'
       AAA1(13)='d'
       AAA1(14)='e'
       AAA1(15)='f'
       AAA1(15)='g'
       AAA1(16)='h'   
       AAA1(17)='i'   
       AAA1(18)='j'   
       AAA1(19)='k'   
       AAA1(20)='l'   
       AAA1(21)='m'   
       AAA1(22)='n'   
       AAA1(23)='o'   
       AAA1(24)='p'   
       AAA1(25)='q'   
       AAA1(26)='r'   
       AAA1(27)='s'   
       AAA1(28)='t'   
       AAA1(29)='y'
       AAA1(30)='v'
       AAA1(31)='w'
       AAA1(32)='x'
       AAA1(33)='y'
       AAA1(34)='x'
       AAA1(35)='z'
       AAA1(36)='A'
       AAA1(37)='B'
       AAA1(38)='C'
       AAA1(39)='D'
       AAA1(40)='E'
       AAA1(41)='F'
       AAA1(42)='G'
       AAA1(43)='H'
       AAA1(44)='I'
       AAA1(45)='J'
       AAA1(46)='K'
       AAA1(47)='L'
       AAA1(48)='M'
       AAA1(49)='N'
       AAA1(50)='O'
       AAA1(51)='P'   
       AAA1(52)='Q'   
       AAA1(53)='R'   
       AAA1(54)='S'   
       AAA1(55)='T'   
       AAA1(56)='Y'
       AAA1(57)='V'
       AAA1(58)='W'
       AAA1(59)='X'
       AAA1(60)='Y'
       AAA1(61)='X'
       AAA1(61)='Z'
       AAA1(62)='E'
       AAA1(63)='F'
       AAA1(64)='G'
       AAA1(65)='H'
       AAA1(66)='I'
       AAA1(67)='J'
       AAA1(68)='K'
       AAA1(69)='L'
       AAA1(70)='M'
       AAA1(71)='N'
       AAA1(72)='O'
       AAA1(73)='P'   
       AAA1(74)='Q'   
       AAA1(75)='R'   
       AAA1(76)='S'   
       AAA1(77)='T'   
       AAA1(78)='Y'
       AAA1(79)='V'
       AAA1(80)='W'
       AAA1(81)='X'
       AAA1(82)='Y'
       AAA1(83)='X'
       AAA1(84)='Z'
       AAA1(85)='E'
       AAA1(86)='F'
       AAA1(87)='G'
       AAA1(88)='H'
       AAA1(89)='I'
       AAA1(90)='J'
       AAA1(91)='K'
       AAA1(92)='L'
       AAA1(93)='M'
       AAA1(94)='N'
       AAA1(95)='O'
       AAA1(96)='P'   
       AAA1(97)='Q'   
       AAA1(98)='R'   
       AAA1(99)='S'   
       AAA1(100)='T'   
C
      RETURN
      END
C===========================================================================
      SUBROUTINE READCOR
      REAL COR(4,4,5000),PROM(4),COR1(4,5000)
      real MC(4,4,5000),M1(4,5000)
C
      COMMON/KORELL/COR,COR1
      COMMON /CAP/IPER,KAPPA
      COMMON/SLY/ISEED
C
C     *** ГЕНЕРАЦИЯ МАТРИЦЫ ***
         DO 30 K=1,IPER
         DO 40 J=1,KAPPA
         DO 10 I=1,KAPPA
            T1=GSU2R(ISEED)
            T1=T1-0.5
            COR(I,J,K)=T1
10       continue
40    CONTINUE 
30    CONTINUE 
C-----------------------------------------
c      КОРРЕКЦИЯ МАТРИЦЫ
        CALL MATRI(MC,M1)
C
        DO 800 K=1,IPER
        DO 810 I=1,4
        DO 820 J=1,4
           COR(I,J,K)=MC(I,J,K)
820     CONTINUE
810     CONTINUE
800     CONTINUE
C
        DO 850 K=1,IPER
        DO 860 I=1,4
          COR1(I,K)=M1(I,K)
860     CONTINUE
850     CONTINUE  
C        
C
C  SOZDANIE PROSTOI MATRIXI
c      DO 120 K=1,IPER
c      DO 20 J=1,4
c        S1=0.0
c        DO 21 I=1,4
c          S1=S1+COR(I,J,K)
c21      continue
c      COR1(J,K)=S1/4.0    
20    continue
120   continue 
      RETURN
      END
C===========================================================================      
      SUBROUTINE DINAMO(SEQ1,SEQ2,L2)
      INTEGER SEQ2(10000)
      INTEGER SEQ1(10000)
      REAL COR(4,4,5000),COR1(4,5000)
      COMMON/KORELL/COR,COR1
C      
      COMMON/SXOD10000/D
      REAL D(0:900,0:900)
      integer DI(0:900,0:900), DJ(0:900,0:900)
      REAL W(10)
C
      COMMON/LET/AAA
      CHARACTER*1 AAA(0:30)
C      
      COMMON/PUT2/S22,IM,JM,I0,J0
      COMMON /CAP/IPER, KAPPA
      COMMON/SAGSAG/ISAG
      REAL DD(0:1000) 
      COMMON/STOLB/DD
C   
C                 WRITE(18,175) L2
175   FORMAT('L2=',I8)
C                 write(18,191) (SEQ1(I),I=1,L2)
191   FORMAT(2X,'SEQ1=',100I2)
C                 write(18,192) (SEQ2(I),I=1,L2)
192   FORMAT(2X,'SEQQ2=',100I2)
      S22=0.0 
      SMAX=-1000000.0
      S2=0.0
c      W(1)=-10000.0
      W(1)=0.0
      V=35.0
      IS=25      
C       V=25000000.0
C  s1 - NUMBER OF ZEROS      
      S1=0.0
C      L1=1
C      L2=L2
C      call SLY_POSL(SEQ2,L1,L2)
C
c
C                     DO 50 J=1,2*L2
C          DO 50 J=1,L2
C             D(J,0)=-V*J
C             D(0,J)=-V*J
C50        CONTINUE             
C
c I -ПОСЛЕДОВАТЕЛЬНОСТЬ НУКЛЕОТИДОВ 
C J - ПОСЛЕДОВАТЕЛЬНОСТЬ НОМЕРОВ (1,2,3,...,L2)
c          D(0,0)=0.0
c          DO 55 I=1,L2
c             D(I,0)=0.0
c55        D(0,I)=0.0
           DO 50 I=0,L2 
50         D(I,0)=0.0

           DO 51 J=0,L2 
51         D(0,J)=DD(J)
C
C            DO 260 J=1,L2
C               I=J-IS-1
C               IF(I.LT.1)  go to 261
C                  D(I,J)=-100000.0
c                  
C261            I=J+IS+1
C               IF(I.GT.L2) GO TO 260
C                 D(I,J)=-100000.0
C260       continue             

C  
          DO 100 J=1,L2
          DO 200 I=1,L2
C           DO 200 I=J-IS,J+IS
C             IF(I.LT.1)  GO TO 200
C             IF(I.GT.L2) GO TO 200
          
             V1=V
             V2=V
             IF ((DI(I-1,J).EQ.I-2).AND.(DJ(I-1,J).EQ.J)) V1=V/4.0
             IF ((DI(I,J-1).EQ.I).AND.(DJ(I,J-1).EQ.J-2)) V2=V/4.0
c
233          FORMAT(2X,'V1=',F8.2,2X,'I=',I6,2X,'J=',I6)
234          FORMAT(2X,'V2=',F8.2,2X,'I=',I6,2X,'J=',I6)
            
             W(2)=D(I-1,J)-V1
             W(3)=D(I,J-1)-V2
             
C       WRITE(18,833) SEQ2(I),SEQ1(J),I,J
C833    FORMAT('SEQ2=',I4,2X,'SEQ1=',I4,2X,'I=',I4,2X,'J=',I4)
       if((I.LT.2).OR.(J.LT.2)) GO TO 933 
       if((DI(I-1,J-1).EQ.I-1).AND.(DJ(I-1,J-1).EQ.J-1)) GO TO 933 
          CALL SOSED(DI,DJ,I,J,KOR,KORJ)
          IF(KOR.LE.0) KOR=I-1
C              if(KORJ.NE.J-1) GO TO 933
          W(4)=D(I-1,J-1)+COR(SEQ2(KOR),SEQ2(I),SEQ1(J))
          GO TO 944 
933    W(4)=D(I-1,J-1)+COR1(SEQ2(I),SEQ1(J))
944    CONTINUE
C                                                                W(4)=D(I-1,J-1)+COR1(SEQ2(I),SEQ1(J))
       CALL MAXX(W,D(I,J),JJ)
       If(S22.GT.D(I,J)) GO TO 555
           IM=I
           JM=J
           S22=D(I,J)
          
C            IF(KOR.LT.1) GO TO 555
C               if(I.eq.J) WRITE(18,722) I,J,KOR,SEQ2(I),SEQ2(KOR),
C     8           SEQ1(J),COR(SEQ2(KOR),SEQ2(I),SEQ1(J)),D(I,J)
C722            FORMAT('I=',I4,2X,'J=',I4,2X,'KOR=',I4,2X,'SEQ2(I)=',
C     8 I4,2X,'SEQ2(KOR)=',I4,2X,'SEQ1(J)=',I4,2X,
C     8 'COR(SEQ2(KOR,SEQ2(I),SEQ1(J))=',2F8.2)
C             
555         GO TO (11,12,13,14) JJ
11          DI(I,J)=I
            DJ(I,J)=J
            GO TO 200
12          DI(I,J)=I-1
            DJ(I,J)=J
            GO TO 200
13          DI(I,J)=I
            DJ(I,J)=J-1
            GO TO 200
14          DI(I,J)=I-1
            DJ(I,J)=J-1
             
C                                  D(I,J)=AMAX1(S1,D(I-1,J)-V, D(I,J-1)-V, D(I-1,J-1)+
C                                  8                MX2(SEQ1(I),SEQ2(J)))
200       CONTINUE          
100       CONTINUE
C
C          write(18,333) (AAA(seq2(I)),I=1,L2)
C333       FORMAT(2X,1000A7)  
C          DO 291 J=1,L2
C291       write(18,292)AAA(SEQ1(J)),(D(I,J),I=1,L2)
C292       FORMAT(A1,2X,1000F7.1)
C
C
c          S22=0.0
c          DO 800 I=1,20
c             IF(S22.GE.D(L2-I,L2)) GO TO 810
c                S22=D(L2-I,L2)
c                IM=L2-I
c                JM=L2
c810       CONTINUE
c             IF(S22.GE.D(L2,L2-I)) GO TO 800
c                S22=D(L2,L2-I)
c                IM=L2
c                JM=L2-I
c800       CONTINUE                
C            S22=(S22-SRED)/(SQRT(DISP))
C          IM=L2
C          JM=L2
C          
C          
c
C          WRITE(18,169) S22,IM,JM
      call PUTI(SEQ1,SEQ2,DI,DJ,IM,JM,I0,J0,L2)
c          print 169, S22,IM,JM
c      write(18, 169) S22,IM,JM
          DO 53 J=0,L2 
53        DD(J)=D(L2,J)
169       FORMAT('S22=',F8.2,2X,'IM=',I6,2X,'JM=',I6)
c         stop
c
C       
          RETURN
          END              

C==================================================================
      SUBROUTINE SEQUENCE(IDL1,SEQ1)
        INTEGER SEQ1(10000)
        COMMON/SLY/ISEED
        COMMON /CAP/IPER,KAPPA
C
        KOK=0
        DO 10 I=1,IDL1 
           KOK=KOK+1
           IF(KOK.GT.IPER) KOK=1
10      SEQ1(I)=KOK
        RETURN
        END
C=================================================================
      SUBROUTINE MAXX(W,SMAX,IMAX)
      REAL W(10)
      SMAX=-100000.0
      DO 10 I=1,4
         IF(SMAX.GE.W(I)) GO TO 10
            SMAX=W(I)
            IMAX=I
10    CONTINUE
      RETURN
      END
C==================================================================
      subroutine SOSED(DI,DJ,I,J,KORI,KORJ)
      integer DI(0:900,0:900),DJ(0:900,0:900)
C  I -INDEX IZUTCHAEMOI POSLED
C  J - INDEX PERIODA    
       ISIG=-100
c      PRINT 711,I,J
711   FORMAT(2X,'I=',I6,2X,'J=',I6)
      I1=I-1
      J1=J-1
c   движение по диагонали
10    IF((DI(I1,J1).NE.I1-1).OR.(DJ(I1,J1).NE.J1-1)) GO TO 20
         KORI=I1
         KORJ=J1 
         ISIG=50
         GO TO 100
C-------------------------------
c   движение по горизонтали       
20    CONTINUE         
      I1=I-1
      J1=J-1
      IF((DI(I1,J1).NE.I1).OR. (DJ(I1,J1).NE.J1-1)) GO TO 30
21    IF((DI(I1,J1).Eq.I1).AND.(DJ(I1,J1).Eq.J1-1)) GO TO 22
             KORI=I1
             KORJ=J1 
             ISIG=50
             GO TO 100
22        J1=J1-1
          go to 21         
C---------
c   движение по вертикали
30    CONTINUE 
      I1=I-1
      J1=J-1
      if((DI(I1,J1).NE.I1-1).OR. (DJ(I1,J1).NE.J1))GO TO 100
31    if((DI(I1,J1).EQ.I1-1).AND.(DJ(I1,J1).Eq.J1))GO TO 32
              KORI=I1
              KORJ=J1
              ISIG=50
              GO TO 100
32         I1=I1-1
           GO TO 31
C      
100   CONTINUE
      IF(ISIG.NE.-100) GO TO 888
         I1=I-1
         J1=J-1
c              KORI=I1
c              KORJ=J1
         PRINT 777
777      FORMAT(5X,'ISIG=-100,  STOP!!!!!!!!!!!!!!!!!!!!!')
         WRITE(18,144) I1,J1, DI(I1,J1),DJ(I1,J1)
144      FORMAT(2X,'I1=',I6,2X,'J1=',I6,2X,'DI=',I6,2X,'DJ=',I6)
         IF(ISIG.EQ.-100) STOP
C
888   CONTINUE
C     WRITE(18, 778) I,J,KOR
778   FORMAT('I=',I6,2X,'J=',I6,2X,'KOR=',I6)
c      
      RETURN
      END
c====================================================================      
       SUBROUTINE UNICOD 
       REAL COR(4,4,5000),COR1(4,5000)
       INTEGER SEQ1(10000)
       COMMON/KORELL/COR,COR1
       COMMON /CAP/IPER,KAPPA
C      
C      J- PREDIDYSII SYMBOL
C      I-izythaemii (tekyshii) symbol 
       DO 20 K=1,IPER
       DO 20 J=1,KAPPA 
          SUM=0.0
          DO 10 I=1,KAPPA
10        SUM=SUM+COR(J,I,K)
C          
20     COR1(J,K)=SUM/(1.0*KAPPA) 

       DO 200 J=1,IPER
200    WRITE(18,155)J,(COR1(I,J),I=1,KAPPA)
155    format('J=',I6,2X,100F8.2)
       RETURN    
       END 
C=====================================================================      
      SUBROUTINE SLY_POSL(SEQ1,L1,L2)
      INTEGER SEQ1(10000),S1
C                          CHARACTER*1 SEQ1(10),S1
      COMMON/SLY/ISEED
C
      DO 10 I=L1,L2-1
20       T1=GSU2R(ISEED)
C                          20       CALL RANDOM(T1)
C                          20       T1=RAND(II5)  
         IF(T1.EQ.1.0) GO TO 20
C
         NOM=INT((L2-I)*T1)+1
         S1=SEQ1(NOM+I)
         SEQ1(NOM+I)=SEQ1(I)
         SEQ1(I)=S1
10    CONTINUE
100   RETURN
      END
C=========================================================================
      SUBROUTINE SLY_POSL1(SEQ1,L1,L2)
      INTEGER*1 SEQ1(200000000),S1
      INTEGER SEQ5(200000000)
      COMMON/SLY/ISEED
C
C  POMNIM 5 B MENYAEM NA SLY
      N5=0 
      DO 300 I=1,IDL
         IF(SEQ1(I).NE.5)GO TO 300
               N5=N5+1       
               SEQ5(N5)=I
301            T1=GSU2R(ISEED)
               IF(T1.GE.1.0) GO TO 301
               SEQ1(I)=INT(T1*4.0)+1
300   CONTINUE
C
C                          CHARACTER*1 SEQ1(10),S1
C
      DO 10 I=L1,L2-1
20       T1=GSU2R(ISEED)
C                          20       CALL RANDOM(T1)
C                          20       T1=RAND(II5)  
         IF(T1.EQ.1.0) GO TO 20
C
         NOM=INT((L2-I)*T1)+1
         IF(SEQ1(NOM+I).EQ.5) GO TO 10
         IF(SEQ1(I).EQ.5) GO TO 10
         S1=SEQ1(NOM+I)
         SEQ1(NOM+I)=SEQ1(I)
         SEQ1(I)=S1
10    CONTINUE
C
C     VSTAVLYAEM 5   
      DO 50 I=1,N5
         SEQ1(SEQ5(I))=5
50    CONTINUE
C
100   RETURN
      END
C=========================================================================
      REAL FUNCTION GSU2R(ISEED)
      INTEGER ISEED,D2P32M
      DOUBLE PRECISION Z,D2P31M,D2PN31,DMOD
      DATA  D2PN31/4.656612873077393D-10/,D2P31M/
     12147483647.D0/,D2P32M/16807/
      Z=DFLOAT(ISEED)
      Z=DMOD(D2P32M*Z,D2P31M)
      GSU2R=Z*D2PN31
      ISEED=Z
      RETURN
      END
C=================================================================
      SUBROUTINE SRR(SPE,I100,SRED,DISP)
      REAL SPE(100000)
      IF(I100.LE.1) I100=10
C
      SRED=0.0
      DO 300 I=1,I100
300   SRED=SRED+SPE(I)*1.0
      SRED=SRED/(I100*1.0)
      IF(SRED.LE.0.0) SRED=100.0
C
      DISP=0.0
      DO 310 I=1,I100
310   DISP=DISP+(1.0*SPE(I)-SRED)**2
      DISP=DISP/(1.0*I100-1.0)
      IF(DISP.LT.0.1) DISP=1.0
      RETURN
      END
C===============================================================
      SUBROUTINE PUTI(SEQ2,SEQ1,DI,DJ,IM,JM,I0,J0,L2)
c SEQ1=SEQ2
c SEQ2=SEQ1
C pomenyal mestami!!!!!
      REAL D(0:900,0:900)
      integer DI(0:900,0:900),DJ(0:900,0:900)
c
      INTEGER SEQ2(10000)
      INTEGER SEQ1(10000)
c      
      INTEGER VIRAV1(10000),V1(10000)
      integer VIRAV2(10000),V2(10000)
      COMMON/VIRAV/V1,V2,K,SINF,M
      INTEGER M(20,800)
      COMMON /CAP/IPER, KAPPA
      COMMON/SXOD10000/D
C      
      COMMON/LET/AAA
      CHARACTER*1 AAA(0:30)      
c      write(18,244) (SEQ1(I),I=1,L2) 
c      write(18,245) (SEQ2(I),I=1,L2)
244   FORMAT('SEQ1=',1000I3)       
245   FORMAT('SEQ2=',1000I3)       
C      INTEGER MCOR(4,4,200)
C      COMMON/CORREL/MCOR0
C      
C-----      
      I=IM
      J=JM
c поменялись местами SEQ1 и SEQ2, поэтому поменял местами координаты точки максимума
c      I=JM
c      J=IM
      K=0
c                         WRITE(18,149) D(I,J),I,J, AAA(SEQ1(I)),AAA(SEQ2(J))
C                         PRINT 511, I,J
C                         511   FORMAT(2X,'I=',I6,2X,'J=',I6)
c      VIRAV1(K)=SEQ1(I)
c      VIRAV2(K)=SEQ2(J)
20    CONTINUE      
c                K2=J-I
c                CALL KORDI(K1,K2,I,J)
         I1=DI(I,J)
         J1=DJ(I,J)
C         WRITE(18,149) D(I1,J1),I,J, AAA(SEQ1(I)),AAA(SEQ2(J))
C                WRITE(18,149) D(I1,J1),I1,J1
c         149      FORMAT(2X,'D=',F8.2,'I=',I6,2X,'J=',I6,2X,2A2)
c         IF((I1.NE.I).OR.(J1.NE.J)) GO TO 10
10       CONTINUE        
         K=K+1 
         IF(((I-I1).EQ.1).AND.((J-J1).EQ.1)) VIRAV1(K)=SEQ1(I)
         IF(((I-I1).EQ.1).AND.((J-J1).EQ.1)) VIRAV2(K)=SEQ2(J)
C         
         IF(((I-I1).EQ.1).AND.((J-J1).EQ.0)) VIRAV1(K)=SEQ1(I)
         IF(((I-I1).EQ.1).AND.((J-J1).EQ.0)) VIRAV2(K)=0
C
         IF(((I-I1).EQ.0).AND.((J-J1).EQ.1)) VIRAV1(K)=0
         IF(((I-I1).EQ.0).AND.((J-J1).EQ.1)) VIRAV2(K)=SEQ2(J)
c         
c         IF((I1.NE.0).AND.(J1.NE.0)) GO TO 11
c                I0=I
c                J0=J
C
C         GO TO 100
C11       CONTINUE
C         
C
         IF(D(I1,J1).GT.0.0)GO TO 15
                I0=I
                J0=J    
                GO TO 100
15       I=I1
         J=J1
      GO TO 20
C      
100   CONTINUE
      DO 50 I=1,K
         V1(I)=VIRAV1(K-I+1)
         V2(I)=VIRAV2(K-I+1)
50    CONTINUE
c      write(18,144) (V1(I),I=1,K) 
c      write(18,145) (V2(I),I=1,K)
144   FORMAT('V1=',1000I3)       
145   FORMAT('V2=',1000I3)       
      DO 60 I=1,KAPPA
      DO 60 J=1,IPER
60    M(I,J)=0
c      write(18,146) I0,J0,IM,JM
146   FORMAT('I0=',I4,2X,'J0=',I4,2X,
     8            'IM=',I4,2X,'JM=',I4) 
      DO 70 I=1,K
         IF((V1(I).EQ.0).OR.(V2(I).EQ.0)) GO TO 70
         M(V1(I),V2(I))=M(V1(I),V2(I))+1
70    CONTINUE
C      call MAT_COR
c      CALL INF_PER(M,SINF)
C         STOP
C      
200   CONTINUE
      RETURN
      END
C===================================================================
       subroutine POINTS(I1,I2,IDL,ISAG,SEQ5,SEQQ2,SRED,DISP,J,SEQ1,
     8                   ZZ,G1,G2,JJ,K7,L2)
       COMMON/SLY/ISEED
       INTEGER*1 SEQ5(200000000)
       INTEGER SEQQ2(10000)
       COMMON/PUT2/S22,IM,JM,I0,J0
       REAL ZZ(4000000)
       INTEGER G1(4000000),G2(4000000),JJ(4000000)
       INTEGER SEQ1(10000)
       COMMON/LEV/SLEV
C
        REAL DD(0:1000) 
        COMMON/STOLB/DD 
C
       NNN=0
       K7=0
       DD=0.0
       DO 120 J=I1,I2,ISAG
c                                    120     CALL EEE_FASTA(IDL,NAME,SEQQ2,KON)       
          IF((I2-J).LT.IDL) GO TO 120
              NNN=0
              DO 121 I=1,IDL
                 IF(SEQ5(I-1+J).NE.5) GO TO 122
123                    T1=GSU2R(ISEED)
                       IF(T1.EQ.1.0) GO TO 123    
                    IF(SEQ5(I-1+J).EQ.5) NNN=NNN+1
                    IF(SEQ5(I-1+J).EQ.5) SEQQ2(I)=INT(T1*4.0)+1
                    GO TO 121
122              SEQQ2(I)=SEQ5(I-1+J)
121           CONTINUE                 
c              IF(NNN.GT.50) GO TO 120
              IF(NNN.GT.0) GO TO 120
c
C     *** slytchainost!!!!! ***
c           L1=1
c           L2=IDL
c           call SLY_POSL(SEQQ2,L1,L2)
C           
C                                    PRINT 567,NNN
567         FORMAT(2X,'REAL, NNN=',I6)           
C                                    IF(KON.EQ.1) GO TO 200
C                                    call UNICOD
c       CALL SEQUENCE(IDL,SEQ1)
C                      DO 10 I=1,10000,100
C                      ISEED=I
       CALL DINAMO(SEQ1,SEQQ2,IDL)
       Z=(S22-SRED)/(SQRT(DISP))
       IF(Z.LT.SLEV) GO TO 120
          IF((I0.EQ.IM).OR.(J0.EQ.JM)) GO TO 120
            CALL KOOR(I0,IM,J0,JM,N0,NM,IDL)
            K7=K7+1
            ZZ(K7)=Z
            G1(K7)=N0-1+J
            G2(K7)=NM-1+J
C            JJ(K7)=J
             JJ(K7)=J+N0-1           
            K8=100*INT((1.0*K7+0.1)/100.0)      
            IF(K8.EQ.K7) PRINT 344, K7,J,myid
344         FORMAT(5X,'I7=',I12,2X,'J=',I12,2x,i4)            
C           PRINT 734,Z,J,S22,K10
C           734    FORMAT('Z='F8.2,2X,'J=',I10,2X,'S22=',F8.2,2X,'K10=',I8)
C
C           WRITE(18,733) Z,J,S22,I0,I0-1+J,IM,IM-1+J,K10,J0,JM
C           733    FORMAT('Z='F8.2,2X,'J=',I10,2X,'S22=',F8.2,2X,'I0=',I3,2X,
C           8 'J0=',I10,2X,'IM=',I3,2X,'JM=',I10,2X,'K10=',I8,5X,
C           8 'J0=',I6,1X,'JM=',I6)
C
C           write(18,144) (V2(I),I=1,K10) 
C           write(18,145) (AAA(V1(I)),I=1,K10)
C           144   FORMAT(5X,'V1=',1000I3)       
C           145   FORMAT(5X,'V2=',1000(2X,A1))       
C     
C     
C     
120   CONTINUE
c      WRITE(18,311) K7
311   FORMAT('K7=',I8)      
c      WRITE(18,224)(JJ(I),I=1,K7)
224   FORMAT('JJ_DO=',10I12)
c      
      RETURN
      END
C===========================================================================  
      SUBROUTINE KOOR(I0,IM,J0,JM,N0,NM,IDL)
      REAL A,B,C,D
      A=1.0*(IM-I0)
      B=1.0*(JM-J0)
      D=1.0*(IDL-J0)
      C=(A*D)/B
      NM=I0+INT(C)
      N0=(NM-IDL)         
C      PRINT*, I0,IM,J0,JM,N0,NM,IDL
      RETURN
      END
C==============================================================================           
      SUBROUTINE GRADIENT(ZZ,G1,G2,JJ,K7,K8)
      REAL ZZ(4000000),PROM(4000000)
      INTEGER G1(4000000),G2(4000000),JJ(4000000)
      INTEGER IP(4000000),KK(4000000)
      K8=0
      if(K7.EQ.0) GO TO 500
C      
      DO 10 I=1,K7
10    IP(I)=I
      CALL AVZ6R(ZZ,K7,IP)
c      write(18,125)(ZZ(I),I=1,K7)
c125   FORMAT('ZZ=',10F8.1)      
c      write(18,126)(IP(I),I=1,K7)
c126   FORMAT('IP=',10I8)
C
C ---- PEREVOROT
      DO 11 I=1,K7
        KK(I)=IP(K7-I+1)
11    PROM(I)=ZZ(K7-I+1)  
      DO 12 I=1,K7
      IP(I)=KK(I)
12    ZZ(I)=PROM(I)  
      
c      write(18,125)(ZZ(I),I=1,K7)
c      write(18,126)(IP(I),I=1,K7)
C      
C           
      DO 20 I=1,K7
20    KK(I)=JJ(IP(I))

      DO 21 I=1,K7
21    JJ(I)=KK(I)
C
c      write(18,127)(JJ(I),I=1,K7)
c127   FORMAT('JJ_0=',10I8)      

      DO 22 I=1,K7
22    KK(I)=G1(IP(I))
      DO 23 I=1,K7
23    G1(I)=KK(I)

      DO 24 I=1,K7
24    KK(I)=G2(IP(I))
      DO 25 I=1,K7
25    G2(I)=KK(I)
C   
      K8=0
      DO 50 I=1,K7
c
         IF(I.NE.1) GO TO 51
           K8=K8+1
           KK(K8)=JJ(I)
           GO TO 50
51       CONTINUE
c
         ISI=0
         DO 52 J=1,I-1
C                         G1(J)----G2(J)
C                  G1(I)-----------------G2(I)
            IF((G1(J).GE.G1(I)).AND.(G1(J).LE.G2(I))) ISI=100
            IF((G2(J).GE.G1(I)).AND.(G2(J).LE.G2(I))) ISI=100
C            
C                  G1(J)-----------------G2(J)
C                        G1(I)-----G2(I)
            IF((G1(I).GE.G1(J)).AND.(G1(I).LE.G2(J))) ISI=100
            IF((G2(I).GE.G1(J)).AND.(G2(I).LE.G2(J))) ISI=100
52       CONTINUE
C
         IF(ISI.EQ.100) GO TO 50
            K8=K8+1
            KK(K8)=JJ(I)
c
50    CONTINUE
       
      DO 60 I=1,K8
60    JJ(I)=KK(I)
c      write(18,233)(JJ(I),I=1,K8)
c233   FORMAT('JJ_1=',10I12)
500   CONTINUE
      RETURN
      END
C=========================================================================
      SUBROUTINE AVZ6R(A,N,IP)
      DIMENSION A(4000000),IP(4000000),IU(21),IL(21)
      INTEGER N,IP,IU,IL,LA,I,M,J,K,IJ,IT,L,ITT
      REAL A,R1,R2,R3,R4,T,TT,R
      DATA  R1/3.75E-01/,R2/5.898437E-01/,R3/3.90625E-02/,
     1R4/2.1875E-01/
      LA=N
      M=1
      I=1
      J=LA
      R=R1
    1 IF(I.EQ.J) GO TO 9
      IF(R.GT.R2) GO TO 2
      R=R+R3
      GO TO 3
    2 R=R-R4
    3 K=I
      IJ=I+(J-I)*IFIX(R)
      T=A(IJ)
      IT=IP(IJ)
      IF(A(I).LE.T) GO TO 4
      A(IJ)=A(I)
      A(I)=T
      T=A(IJ)
      IP(IJ)=IP(I)
      IP(I)=IT
      IT=IP(IJ)
    4 L=J
      IF(A(J).GE.T) GO TO 6
      A(IJ)=A(J)
      A(J)=T
      T=A(IJ)
      IP(IJ)=IP(J)
      IP(J)=IT
      IT=IP(IJ)
      IF(A(I).LE.T) GO TO 6
      A(IJ)=A(I)
      A(I)=T
      T=A(IJ)
      IP(IJ)=IP(I)
      IP(I)=IT
      IT=IP(IJ)
      GO TO 6
    5 TT=A(L)
      A(L)=A(K)
      A(K)=TT
      ITT=IP(L)
      IP(L)=IP(K)
      IP(K)=ITT
    6 L=L-1
      IF(A(L).GT.T) GO TO 6
    7 K=K+1
      IF(A(K).LT.T) GO TO 7
      IF(K.LE.L) GO TO 5
      IF(L-I.LE.J-K) GO TO 8
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 10
    8 IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 10
    9 M=M-1
      IF(M.EQ.0) GO TO 13
      I=IL(M)
      J=IU(M)
   10 IF(J-I.GE.1) GO TO 3
      IF(I.EQ.1) GO TO 1
      I=I-1
   11 I=I+1
      IF(I.EQ.J) GO TO 9
      T=A(I+1)
      IT=IP(I+1)
      IF(A(I).LE.T) GO TO 11
      K=I
   12 A(K+1)=A(K)
      IP(K+1)=IP(K)
      K=K-1
      IF(T.LT.A(K)) GO TO 12
      A(K+1)=T
      IP(K+1)=IT
      GO TO 11
   13 RETURN
      END
C==========================================================================
      SUBROUTINE INVER(SEQ5,IDL1)
      INTEGER*1 SEQ5(200000000),I1,I2
      INTEGER*1 COM(10)
      COM(1)=2 
      COM(2)=1 
      COM(3)=4 
      COM(4)=3 
      COM(5)=5 
      COM(6)=6 
      COM(7)=7 
      COM(8)=8 
      COM(9)=9 
C
      NOM=INT((1.0*IDL1+0.1)/2.0)
C
      DO 10 I=1,NOM
         I1=SEQ5(I)
         I2=SEQ5(IDL1-I+1)
         SEQ5(I)=I2
         SEQ5(IDL1-I+1)=I1
10    CONTINUE
C
c      DO 20 J=1,IDL1
c        I1=SEQ5(J)    
c20    SEQ5(J)=COM(I1)
C 
      RETURN
      END
C==============================================
        SUBROUTINE MATRI(MC,M1)
        REAL VER1(1000),VER2(20), PQ 
        COMMON /CAP/IPER,KAPPA
C        
        REAL COR(4,4,5000),COR1(4,5000)
        COMMON/KORELL/COR,COR1
C        
        REAL MX1(2,20,5000),MX2(20,5000),MX3(20,5000)
        REAL AA,BB,A,B,C,D,MM,S1,S2,S3
        INTEGER M(20,5000)
        real MC(4,4,5000),M1(4,5000)
C--------        
c        OPEN(19,FILE='277.txN',STATUS='old')
c        OPEN(18,FILE='277.res',STATUS='REPLACE')
c***************************************
c              IPER=177
c***************************************
c--------        
        KAPPA=16
C               CALL READCOR
c  racthet R 
c        CALL SXI(R,COR)
c        print 733, R
c733     FORMAT(2X,'R=',F10.2)
c   zadanie R ecli net, to soxranyaetsya R isxodnoe
        R=(45000.0/214.0)*IPER
        WRITE(18,133) R
133     FORMAT(2X,'R_ZADANNYI=',F10.2)                
C        
c        DO 21 I=1,4
c        DO 21 J=1,4
c21      WRITE(18,156)I,(COR(I,J,K),K=1,IPER)
c156     FORMAT('I=',I4,2X,'COR=',1000F6.1) 
        
        CALL BTRANS (COR,MX3)
C        DO 20 I=1,KAPPA
C20      WRITE(18,155)I,(MX3(I,J),J=1,IPER)
C155     FORMAT('I=',I4,2X,'MX3=',1000F6.1) 
C           real MX2(20,800),MC(4,4,800)
        KAPPA=16
        DO 10 I=1,IPER
10      VER1(I)=1.0/(IPER*1.0)
        DO 11 I=1,KAPPA
11      VER2(I)=1.0/(KAPPA*1.0)
C        
         call AB(VER1,VER2,PQ)
         
c
C        RR1=0.0
C        DO 80 I=1,KAPPA
C        DO 80 J=1,IPER
C80      RR1=RR1+MX3(I,J)**2.0
C        PRINT 777, RR1
C777     FORMAT(5X'RR1=',F10.2)        

         KK=1
         DO 810 I=1,KAPPA
         DO 810 J=1,IPER
810      MX1(KK,I,J)=MX3(I,J)
C  
         CALL MAT_PREOB(KK,VER1,VER2,SS,MX1,MX2,PQ,M,R)
c *** isxodnya matriza ***         
c         DO 12 I=1,KAPPA
c12       write(18,1)I,(MX1(KK,I,J),J=1,IPER)
c1        FORMAT('I=',I6,2X,'MX1=',600F8.2)
C *** posle preobrazovanya **** 
c         DO 13 I=1,KAPPA
c13       write(18,3)I,(MX2(I,J),J=1,IPER)
c3        FORMAT('I=',I6,2X,'MX2=',600F8.2)
C
         CALL DTRANS (MX2,MC,M1)
C ***  posle preobrazovanya edinitchaya matriza  v formate 4 na 5000 posle preob*** 
c            DO 25 K=1,IPER
c               WRITE(18,177)K,(M1(I,K),I=1,4)
c177            FORMAT('SINGL=',2X,'K=',2X,i8,2x,'Z=',100F8.2)
25          CONTINUE
C ***  posle preobrazovanya edinitchaya matriza  v formate 4 na 5000 posle preob*** 
c          DO 50 I=1,IPER
c                   DO 172 J=1,4
c172                write(18,171) I,J,(MC(J,K,I),K=1,4)
c171                FORMAT('COR=',2x,'I=',I6,2X,'J=',I1,2X,4F8.2)
50        CONTINUE    
        CALL SXI(R,MC)
        print 733,R
        WRITE(18,733) R
733     FORMAT(2X,'R_POSLE=',F10.2)
c            STOP
C 
         END         
C===================================================================
        SUBROUTINE MAT_PREOB(KI,VER1,VER2,SS,MX1,MX2,PQ,M,R)
C         DOUBLE PRECISION VER1(20),VER2(20),PQ 
        COMMON /CAP/IPER,KAPPA
        REAL VER1(1000),VER2(20), PQ
        REAL MX1(2,20,5000),MX2(20,5000)
        REAL K0,K,PI,S5
        REAL AA,BB,A,B,C,D,MM,S1,S2,S3
        INTEGER M(20,800)
C        DOUBLE PRECISION K0,K,PI,S5
C        DOUBLE PRECISION AA,BB,A,B,C,D,MM,S1,S2,S3
c        
C        PQ=0
C        DO 12 I=1,KAPPA
C        DO 12 J=1,IPER
C12      PQ=PQ+(VER1(J)*VER2(I))**2.0
C
c        K=-0.54
c
c      *****Kd*******
        K=-0.5
c      ***************
        K0=0.0
        DO 10 I=1,KAPPA
        DO 10 J=1,IPER
           MM=MX1(KI,I,J)
10      K0=K0+VER1(J)*VER2(I)*MM
        AA=(K-K0)/PQ
        BB=K/PQ
        S1=(BB**2.0)*PQ
C               WRITE(18,166) K0,AA,BB
166     FORMAT(2X,'KO=',F12.6,2X,'AA=',F12.6,2X,'BB=',F12.6)
C               STOP
C        
        S2=0.0
        S5=0.0
        DO 20 I=1,KAPPA
        DO 20 J=1,IPER
           PI=VER1(J)*VER2(I)
C           S5=S5+PI*(MX1(KI,I,J)+PI*AA)
            S5=S5+PI*BB*PI
c           WRITE(18,165) PI*BB, MX1(KI,I,J)+PI*AA, 
c     8                          MX1(KI,I,J)+PI*AA-PI*BB
165        FORMAT(2X,'Z=',F12.6,2X,'Y=',F12.6,2X,'R=',F12.6)
C
20      S2=S2+PI*(MX1(KI,I,J)+PI*(AA-BB))
C               WRITE(18,134) S5
134     FORMAT(2X,'S5=',F20.12)        
        S2=S2*2.0*BB
c               STOP
C        
        S3=0.0
        DO 30 I=1,KAPPA
        DO 30 J=1,IPER
           PI=VER1(J)*VER2(I)
30      S3=S3+(MX1(KI,I,J)+PI*(AA-BB))**2.0
        GO TO 255
C       IF(S3.GT.0.0) GO TO 255
              WRITE(18,256) PI,AA,BB,KI
256           FORMAT('PI=',F8.2,2X,'AA=',F8.2,2X,
     8               'BB=',F8.2,2X,'KI=',I6)
C
              WRITE(18,257) (VER1(J),J=1,IPER)
257           FORMAT('VER1=',2X,1000F8.4)
C 
              WRITE(18,258) (VER2(I),J=1,KAPPA)
258           FORMAT('VER2=',2X,1000F8.4)
C
              DO 259 I=1,KAPPA
259           WRITE(18,260)I,(MX1(KI,I,J),J=1,IPER)
260           FORMAT('MX1=',2X,I6,1X,1000F8.2)
              DO 261 I=1,KAPPA
261           WRITE(18,262)I,(M(I,J),J=1,IPER)
262           FORMAT('M=',2X,I6,1X,1000I6)
255     CONTINUE
        IF(S3.LE.0.0) GO TO 700
C       AX**2+BX+C=0
        A=S3
        B=S2
C                                     SXI=SQRT(R)
c                                     SXI=1050.0
c                                     SXI=800.0 
c        IF(IPER.EQ.2) C=-(SXI/(7.0**0.52))*(IPER**0.52)+S1
c        IF(IPER.EQ.3) C=-(SXI/(7.0**0.6))*(IPER**0.6)+S1
c        IF(IPER.EQ.4) C=-(SXI/(7.0**0.61))*(IPER**0.61)+S1
c        IF(IPER.EQ.5) C=-(SXI/(7.0**0.59))*(IPER**0.59)+S1
c        IF(IPER.EQ.6) C=-(SXI/(7.0**0.56))*(IPER**0.56)+S1
c        IF(IPER.GE.7) C=-(SXI/(7.0**0.55))*(IPER**0.55)+S1
c
c       C RAVNO S1 minys zelaemii R kvadrat (tut prosto R)
c
        C=-R+S1
c              C=-(800.0/(7.0**0.4))*(IPER**0.4)+S1
        D=B**2.0-4.0*A*C
        X1=(-B+(D**0.5))/(2.0*A)
        X2=(-B-(D**0.5))/(2.0*A)
C        
C               WRITE(18,155) A,B,C,D,PQ
155     FORMAT(2X,'A=',F12.6,2X,'B=',F12.6,2X,'C=',F12.6,2X,
     8        'D=',F18.6,2X,'PQ=',F12.6)
        I1=0     
        DO 60 I=1,KAPPA
        DO 60 J=1,IPER
           PI=VER1(J)*VER2(I)
           I1=I1+1
60      MX2(I,J)=PI*BB+X1*(MX1(KI,I,J)+PI*(AA-BB))
C
        SS=0.0
        DO 70 I=1,KAPPA
        DO 70 J=1,IPER
70      SS=SS+(VER1(J)*VER2(I))*MX2(I,J)
C
        SS1=0.0
        DO 80 I=1,KAPPA
        DO 80 J=1,IPER
80      SS1=SS1+MX2(I,J)**2.0
C
C              WRITE(18,156) KI,X1,X2
        WRITE(18,157) SS,SS1
156     FORMAT(2X,'KI=',I6,2X,'X1=',F14.8,2X,'X2=',F14.8)
157     FORMAT(2X,'SUMA1=',F25.2,2X,'SUMA2=',F25.2)
        GO TO 710 
700         DO 701 I=1,KAPPA
            DO 701 J=1,IPER
701         MX2(I,J)=0.0
C
710     CONTINUE
        RETURN
        END
C==================================================================
        SUBROUTINE AB(VER1,VER2,PQ)
        COMMON /CAP/IPER,KAPPA
        REAL VER1(1000),VER2(20),PQ
C          DOUBLE PRECISION VER1(20),VER2(20),PQ
C        
        PQ=0.0
        DO 10 I=1,KAPPA
        DO 10 J=1,IPER
10      PQ=PQ+(VER1(J)*VER2(I))**2.0
C           WRITE(18,155) PQ
155     FORMAT('PQ=',F20.12)
C
        RETURN
        END
C===================================================================
      subroutine DTRANS (MX2,MC,M1)
      real MX2(20,800),MC(4,4,800),M1(4,800)
      COMMON /CAP/IPER, KAPPA
C      
      DO 10 I=1,KAPPA
      DO 10 J=1,IPER
         S=(1.0*I-0.1)/(4.0)
         I2=INT(S)+1
         I1=I-(I2-1)*4
10    MC(I1,I2,J)=MX2(I,J)
C  SOZDANIE PROSTOI MATRIXI
      DO 20 K=1,IPER
      DO 20 J=1,4
        S1=0.0
        DO 21 I=1,4
21         S1=S1+MC(I,J,K)
20    M1(J,K)=S1/4.0    
c
      RETURN
      END
C=======================================================================           
      subroutine BTRANS (MC,MX2)
      real MX2(20,5000),MC(4,4,5000)
      COMMON /CAP/IPER, KAPPA
C      
      DO 10 I1=1,4
      DO 10 I2=1,4
      DO 10 J=1,IPER
         I=I1+4*(I2-1)
10    MX2(I,J)=MC(I1,I2,J)
      RETURN
      END
C======================================================================
      SUBROUTINE SXI(R,COR)
      REAL COR(4,4,5000)
C           COMMON/KORELL/COR,COR1
      COMMON /CAP/IPER,KAPPA
      R=0.0
                   DO 175 K=1,IPER
                   DO 175 J=1,4
                   DO 175 I=1,4          
175   R=R+COR(I,J,K)**2        
      RETURN
      END
C============================================================================
      SUBROUTINE NORMA_C(COR)
      COMMON /CAP/IPER, KAPPA
c         COMMON/MMM/MAT,X,Y,LONG
      INTEGER MAT(4,4,5000)
      REAL X(4),Y(4),LONG,XX
      REAL COR(4,4,5000),P(4)
      COMMON/MATMAT/MAT
cC            WRITE(18,157) (X(I),I=1,IPER)
C            157   FORMAT(2X,'X=',1000I6)
C            WRITE(18,158) (Y(I),I=1,4)
C            158   FORMAT(2X,'Y=',1000I6)
C            WRITE(18,159) LONG
C            159   FORMAT('LONG=',I8)
C
C      DO 200 K=1,IPER
C      DO 201 I=1,4
C         WRITE(18,202)K,I,(MAT(I,J,K),J=1,4)
C202      FORMAT('K=',I8,2X,'I=',I8,2X,'MAT=',4I8)
C201   CONTINUE
C200   CONTINUE
C      STOP
C
          X=0.0
          Y=0.0
          DO 40 K=1,IPER
             DO 141 I=1,4
             DO 41 J=1,4
                X(J)=X(J)+MAT(I,J,K)
                Y(I)=Y(I)+MAT(I,J,K)
41           CONTINUE
141           CONTINUE
40        CONTINUE
C
          SUM1=X(1)+X(2)+X(3)+X(4)
          SUM2=Y(1)+Y(2)+Y(3)+Y(4)
C          PRINT 777,SUM1
C777       FORMAT(2X,'SUM1=',F8.1)
C          
          DO 50 J=1,4
            X(J)=X(J)/SUM1
50        continue
C
          DO 51 I=1,4
           Y(I)=Y(I)/SUM2
51        continue
c                     
      DO 60 K=1,IPER
            XX=0.0
            DO 153 I=1,4
            DO 53 J=1,4
C  XX -obshee kol par            
              XX=XX+MAT(I,J,K)
53          continue
153         continue
C
             DO 152 I=1,4    
             DO 52 J=1,4
                PP=Y(I)*X(J)
                  IF(PP.LE.0.0) COR(I,J,K)=0.0
                  IF(PP.LE.0.0) GO TO 52
                  IF(XX.LE.0.0) COR(I,J,K)=0.0
                  IF(XX.LE.0.0) GO TO 52
                  IF(PP.GE.1.0) COR(I,J,K)=0.0
                  IF(PP.GE.1.0) GO TO 52
                SR=XX*PP
                SIG=SQRT(SR*(1.0-PP))
C                   IF(SIG.LT.1.0) SIG=1.0
                COR(I,J,K)=(MAT(I,J,K)-SR)/SIG
C                   COR(I,J,K)=((MAT(I,J,K)-SR)**2)/SR
C            
C                   WRITE(19,178) SR,SIG,P,KOL,PAR(J,K),COR(J,K,I),J,K,I
C178                FORMAT('SR=',F8.4,2X,'SIG=',F8.4,2X,'P=',F8.4,2X,
C     8             'KOL=',I6,2X,'PAR=',F8.4,2X,'COR=',F8.2,2X,
C     8             'J=',I4,2X,'K=',I4,2X,'I=',I6)
52           CONTINUE
152          continue
60    CONTINUE
c                   DO 172 I=1,4
c                   DO 172 J=1,4
c172                write(18,171) I,J,(COR(I,J,K),K=1,IPER)
c171                FORMAT('COR=',2x,'I=',I6,2X,'J=',I1,2X,100F8.2)
c                   DO 272 I=1,4
c                   DO 272 J=1,4
c272                write(18,271) I,J,(MAT(I,J,K),K=1,IPER)
c271                FORMAT('MCOR=',2x,'I=',I6,2X,'J=',I1,2X,100I8)
C       STOP
       RETURN
       END    
C====================================================================================
      SUBROUTINE NULL
      INTEGER MAT(4,4,5000)
      COMMON/MATMAT/MAT
      COMMON /CAP/IPER, KAPPA
      DO 10 K=1,IPER
      DO 20 I=1,4
      DO 30 J=1,4
         MAT(I,J,K)=0
30    CONTINUE
20    CONTINUE
10    CONTINUE
      RETURN
      END
C====================================================================================  
      SUBROUTINE INF_PER(M2,SINF)
      REAL X(20),Y(5000)
      REAL LOGAR(20010)
      COMMON/EE/I55,LOGAR
      COMMON /CAP/IPER, KAPPA
      REAL M(4,4,5000),M1(20,5000)
      INTEGER M2(4,4,5000)
C------
c      IF(I55.GT.10) GO TO 500
c               I55=10000
c               LOGAR(1)=0
c               DO 600 I=1,20000
c600            LOGAR(I+1)=I*1.0*ALOG(1.0*I)
c500   CONTINUE
      DO 24 K=1,IPER
      DO 25 I=1,KAPPA
      DO 26 J=1,KAPPA  
        M(I,J,K)=M2(I,J,K)*1.0
C        PRINT 777, M(I,J,K)
C777     FORMAT('M=',F8.1)
26    CONTINUE
25    CONTINUE
24    CONTINUE 
C
      call BTRANS (M,M1)
C-------
C      DO 224 K=1,IPER
C      DO 226 J=1,KAPPA  
C         PRINT 777, M1(J,K)
C777     FORMAT('M=',F8.1)
C226    CONTINUE
C224    CONTINUE 

      KAPPA1=16
C
      DO 10 I=1,KAPPA1
10    X(I)=0.0
      DO 20 I=1,IPER
20    Y(I)=0.0
      SLONG=0.0
      DO 30 I=1,KAPPA1
      DO 30 J=1,IPER
         SLONG=SLONG+M1(I,J)
         X(I)=X(I)+M1(I,J)
30    Y(J)=Y(J)+M1(I,J)
C      PRINT 888,SLONG
C888   FORMAT('SLONG=',F8.1)  
C------------------------------------------------
      IF(SLONG.LE.0.0) GO TO 1000
      S1=0.0
           DO 100 I=1,KAPPA1
              IF(X(I).LE.0.0) GO TO 100
              S1=S1+ALOG(1.0*X(I))
100        CONTINUE
      S2=0.0
           DO 200 J=1,IPER
              IF(Y(J).LE.0.0) GO TO 200
              S2=S2+ALOG(1.0*Y(J))
200        CONTINUE
      S3=0.0
           DO 300 I=1,KAPPA1
           DO 300 J=1,IPER
               IF(M1(I,J).LE.0.0) GO TO 300
               S3=S3+ALOG(1.0*M1(I,J))
300        CONTINUE
C                         WRITE(18,22) S1,S2,S3,LONG
22         FORMAT('S=',3F10.4,I6)
      SINF=2.0*S3-2.0*S1-2.0*S2+2.0*SLONG*ALOG(1.0*SLONG)+0.001
      SINF=SQRT(2.0*SINF)-SQRT(2.0*(KAPPA1-1)*(IPER-1.0)-1.0)
1000  CONTINUE
      RETURN
      END
C======================================================
c      SUBROUTINE SLY_TRI(SEQ1,IDL)
c      COMMON/SLY/ISEED
c      INTEGER*1 SEQ1(5000000)
c      INTEGER*1 SEQ2(83000000)
c      INTEGER SEQ5(10000000)
c
c      L=3*INT((1.0*IDL+0.1)/3.0)
c      IF(L.NE.IDL) PRINT 566
c566   FORMAT('L NE.IDL!!!!!!!!!!!!!!!!!')
c      IDL=L
C
C  POMNIM 5 B MENYAEM NA SLY
c      N5=0 
c      DO 300 I=1,IDL
c         IF(SEQ1(I).NE.5)GO TO 300
c               N5=N5+1       
c               SEQ5(N5)=I
c301            T1=GSU2R(ISEED)
c               IF(T1.GE.1.0) GO TO 301
c               SEQ1(I)=INT(T1*4.0)+1
c300   CONTINUE
C
C
c      N=0
c      DO 10 I=1,IDL-2,3
c         K=SEQ1(I)+(SEQ1(I+1)-1)*4+(SEQ1(I+2)-1)*16
c         N=N+1
c         SEQ2(N)=K
c10    CONTINUE     
C         
C
c         L1=1
c         CALL SLY_POSL5(SEQ2,L1,N)
C
c      DO 20 I=1,N
c         S=(1.0*SEQ2(I)-0.1)/(16.0)
c         N3=INT(S)+1
c         I1=SEQ2(I)-(N3-1)*16
C
c         S=(1.0*I1-0.1)/(4.0)
c         N2=INT(S)+1
c         N1=I1-(N2-1)*4
c         SEQ1(3*I)=N3
c         SEQ1(3*I-1)=N2 
c         SEQ1(3*I-2)=N1
c20    CONTINUE
C     VSTAVLYAEM 5   
c      DO 50 I=1,N5
c         SEQ1(SEQ5(I))=5
c50    CONTINUE
C 
c       RETURN
c       END
C==========================================================================
c      SUBROUTINE SLY_POSL5(SEQ1,L1,L2)
c      INTEGER*1 SEQ1(83000000),S1
c      COMMON/SLY/ISEED
C
c      DO 10 I=L1,L2-1
c20       T1=GSU2R(ISEED)
c         IF(T1.EQ.1.0) GO TO 20
C
c         NOM=INT((L2-I)*T1)+1
c         S1=SEQ1(NOM+I)
c         SEQ1(NOM+I)=SEQ1(I)
c         SEQ1(I)=S1
c10    CONTINUE
c100   RETURN
c      END
C===================================================================
      SUBROUTINE ADD_INV(SEQ5,IDL1)
      INTEGER*1 SEQ5(120000000)
      DO 10 I=1,IDL1
         J=IDL1-I+1
         SEQ5(I+IDL1)=SEQ5(J)
10    CONTINUE
      IDL1=IDL1*2
      RETURN
      END
C===================================================================

 

