       PROGRAM vstavka
       character*1 NNN,BO
       CHARACTER*1 SEQ(100000000)    
       CHARACTER*3 ZZ,PPP
       character*100 NAM
C
       OPEN(25,FILE='uni1.txt',STATUS='old')
       OPEN(17,FILE='genome.fasta',STATUS='old')
       ZZ='ZZ='
       NNN='N'
       BO='>'
CZZ=    5.41  NA=    825334  KO=    826026  NIT= +  NOM=     212
C         78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99100101102103104105106107108109110111112113114115116117118119120121122123124125126127128129130131132133134135136137138139140141142  *  *  *143144145146147148149150  *  *  *151152153154155156157158159160161162163164165166167168169170171172173174175176  *  *  *  *177178179180181182183184185186187188189190191192193194  *  *195196197  *198199200201202203204205206  *  *207208209210211212213214215216217218219220221222223224225226227228229230231232233234235236237238239240241242243244245  *246247248249250251252253254255256257258259260261262263264265266267268269270271272273274275276277278279280281282283284285286287288289  *  *  *  *  *  *290291292293294295296297298299300301  *  *302303304305306307308309310311312313314315316317318319320321322323324325326327328329330331332333  *  *  *334335336337338339340341342343344345346347348  *  *  *349350351352353354355356357358359360361362363364365366367368369370371372373374375376377378379380381382383384385386387388389390391392393394395396397398399  *400401402403404405406407408409410411412413414415416417418419420421422423424425426427428429430431432433434435436437438439440441442443444445446447448449450451452  *453454455456457458459460461462463464465466467468469470471472473474475476477478479480481482483484485486487488489490491492493494495496497498499500501502503504505506507508509510511512513514515516517518519520521522523524525526527528529530531532533534535536537538539540541542543544545  *546547548549550551552553  *554555556557558559560561562  *  *563564565566567568569570571572573574  *  *575576577578579580581582583584585586587588589590591592593594595596597598  *599600601602603604605606607608609610611612613614615616617618619620621622623624625626627628629630631632  *633634635636637638639640641642643644645646647648649650651652653654655656657658659660661662663664665666667668669670671672673674675  *  *676677678679680681682683684685686687688689690691692693694695696697698699700701702703704705706707708709710711712713714715716717718719720721722723724725726727728729730731732733734735
C          C  G  C  G  C  T  C  A  A  A  A  A  T  C  C  G  A  A  G  A  A  G  T  A  C  A  C  C  A  T  C  A  A  C  G  G  G  A  A  G  A  A  C  T  A  T  C  T  G  C  T  G  G  T  G  G  A  G  T  T  C  G  A  C  G  A  C  T  T  C  C  A  C  A  T  T  C  C  G  C  A  G  A  G  C  A  T  T  A  C  C  G  A  C  A  C  T  T  T  C  T  A  T  G  A  A  T  T  G  A  C  T  C  T  G  G  C  A  G  G  A  G  C  G  C  G  C  C  C  G  A  T  C  C  T  T  A  C  G  C  A  T  C  C  T  G  A  G  C  G  C  A  A  C  C  T  G  C  T  G  C  T  G  C  A  G  A  A  A  G  A  G  C  C  C  G  A  G  C  G  C  A  T  A  G  C  G  G  A  A  T  G  G  C  T  G  G  A  A  C  A  C  G  G  C  T  G  C  C  T  C  G  T  G  C  A  G  A  T  C  A  C  G  A  G  C  A  A  C  T  C  G  C  T  C  A  C  C  G  G  C  C  G  C  T  T  C  G  G  C  A  A  G  A  C  G  G  C  G  C  A  G  A  G  C  A  T  G  G  G  C  T  T  C  C  G  G  C  T  G  C  T  C  G  A  C  A  A  A  G  G  C  T  G  G  G  T  G  C  A  T  T  T  T  C  T  G  G  C  G  A  C  C  G  A  T  G  C  G  C  A  C  A  A  T  A  C  C  G  A  A  T  C  G  C  G  T  C  C  A  C  C  A  G  T  C  A  T  G  A  G  C  G  A  G  G  C  A  C  G  G  C  G  G  C  T  G  G  T  G  G  C  C  G  A  A  C  G  G  T  A  C  G  G  C  G  C  G  G  A  A  G  T  G  G  C  C  G  A  G  C  G  G  C  T  G  T  G  T  G  T  G  A  C  C  A  A  T  C  C  C  C  T  G  G  C  T  G  T  G  T  T  T  G  A  A  G  G  C  C  G  T  C  C  G  G  T  C  C  C  G  A  A  G  G  C  T  C  C  C  C  G  T  T  C  G  C  G  G  T  A  T  G  A  C  G  A  G  G  C  G  G  A  C  G  A  C  G  A  A  C  C  C  G  G  A  C  C  G  A  G  A  C  A  C  G  G  T  G  G  A  T  T  T  T  T  C  T  C  A  C  G  C  C  T  A  C  T  T  G  G  C  C  G  C  A  A  A  T  A  G  C  T  C  A  T  A  A  C  A  G  T  C  G  C  T  G  C  T  A  C  A  A  A  C  C  G  A  G  T  C  C  G  C  C  *  G  T  C  T  A  C  G  G  C  T  A  G  A  A  T  T  T  G  T  C  C  G  G  T  G  A  *  *  *  *  *  *  T  G  A  A  A  T  G  C  G  G  T  C  C  G  G  C  G  G  C  G  A  A  A  A  A  G  A  G  C  A  C  G  G  C  G  G  C  C  G  C  T  A  C  A  T  C  C  T  C  A  G  C  C  T  T  G  C  C  A  T  T  G  C  G  G  C  G  C  A  T  T  G  G  C  G  T  C  T  T  T  T  C  G  G  C  A  A  A  G  T  G  C
       CALL SEQUENCE(SEQ,IDL,NAM)
10     READ(25,1,END=100) PPP
1      FORMAT(A3)
       IF(PPP.NE.ZZ) GO TO 10
C
         BACKSPACE 25
         READ(25,2) I1,I2
C               PRINT 2, I1,I2
2        FORMAT(16X,I10,5X,I10)
          J2=MAX0(I1,I2)
          J1=MIN0(I1,I2)
          PRINT 2, J1,J2
          DO 20 I=J1,J2
             SEQ(I)=NNN
20        CONTINUE
          GO TO 10    
C------------------------------
100   continue
      REWIND 17 
      write(17,3) BO,NAM
3     FORMAT(A1,A100)
      WRITE(17,6)(SEQ(I),I=1,IDL)
6     FORMAT(60A1)
      WRITE(17,3) BO
C     
      END
C=====================================================
      SUBROUTINE SEQUENCE(SEQ,IDL,NAM) 
      CHARACTER*1 SEQ(100000000)    
      CHARACTER*100 NAM
      CHARACTER*1 PROM(100),BO,PU,PP
      BO='>'
      PU=' '
      NNN=0
C
10     READ(17,1,END=100) PP
1      FORMAT(A1)
       IF(PP.NE.BO) GO TO 10
         BACKSPACE 17
         READ(17,2,END=100) NAM
2        FORMAT(1X,A100)
C
50       READ(17,4) (PROM(I),I=1,100)
4        FORMAT(100A1)
         IF(PROM(1).EQ.BO) GO TO 100 
C
         DO 20 I=1,100
            I1=101-I
            IF(PROM(I1).NE.PU) K=I1
            IF(PROM(I1).NE.PU) GO TO 21
20       CONTINUE
21       CONTINUE

           DO 30 I=1,K
             SEQ(I+NNN)=PROM(I)
30         CONTINUE
           NNN=NNN+K
           GO TO 50
100   CONTINUE
      IDL=NNN
      RETURN
      END
C========================================      











       
