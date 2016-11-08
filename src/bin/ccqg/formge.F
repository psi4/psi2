      SUBROUTINE FORMGE(GEMAT,TAU,T2,T1,Z2,Z1,NO,NV,NO2,NV2,NTRO,NTRV,
     .                  IOFF,NTR,Q1,Q2,Q3,VALO,JOUT,OPTION,
     .                  ORBSYM,FLOV,ZLX,ITR,ITV,NT3,NIRRED,NTAU,NOV,
     .                  T3OFF,T3ADD,W3,Y3,NDIMT2,TOFF,TADD,
     .                  EOFF,EADD,NSLOV,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      COMMON/TAPES/ITAP81,ITAP82,ITAP83,ITAP99,itap76,i76,itap74
      CHARACTER*4 OPTION
      REAL*8 VAL,VAL1,VAL2,VAL3,VAL4,VAL5,VAL6,VAL7,VAL8,
     .       FACJK,FACIJ,FACIK
      REAL*8 GEMAT(NSLOV),Z2(NDIMT2),T2(NDIMT2),VALO(NO),
     .       T1(NO,NV),Z1(NO,NV),TAU(NDIMT2),Q1(NV),Q2(NV),Q3(NV),
     .       W3(NTAU),Y3(NTAU)
      INTEGER IOFF(NTR),ORBSYM(NO+NV),FLOV(NIRRED,4),ZLX(NV,NV),
     .        ITR(NO),ITV(NOV),NT3(NIRRED),
     .        T3OFF(NO,NO,2,NIRRED),T3ADD(NV,NV,NIRRED),
     .        TOFF(NO,NO,2),TADD(NV,NV),EOFF(NO*NV),EADD(NO*NO),
     .        FZO(NO),FZV(NV)
C
      CALL ZERO(GEMAT,NSLOV)
c
      if(option.eq.'CCT ')then
c        write(6,*)' reading eden',nslov,'  long'
         call wreadw(itap76,gemat,intowp(nslov),i76,i76)
c        write(6,'(8f10.6)')gemat
      endif
C
      IF(OPTION.EQ.'SDT1')THEN
      I82=1
      DO 16 A=1,NV
         IF(FZV(A).EQ.1)GO TO 16
         ASYM=ORBSYM(A+NO)
         DO 26 I=1,NO
            IF(FZO(I).EQ.1)GO TO 26
            ISYM=ORBSYM(I)
            IASYM=IEOR(ASYM,ISYM)+1
            DIM=NT3(IASYM)
            CALL WREADW(ITAP82,Y3,INTOWP(DIM),I82,I82)
            AI=ITR(A)+I
            DO 36 J=1,NO
               IF(FZO(J).EQ.1)GO TO 36
               JSYM=ORBSYM(J)
               KSYM=IEOR(IASYM-1,JSYM)
               FK=FLOV(KSYM+1,1)
               LK=FLOV(KSYM+1,2)
               DO 46 K=FK,LK
                  JK=ITR(MAX0(J,K))+MIN0(J,K)
                  FACJK=1.0D+00
                  IF(J.EQ.K) FACJK=2.0D+00
                  AIJK=EOFF(AI)+EADD(JK)
                  VAL1=0.0D+00
                  DO 56 M=1,NO
                     IF(FZO(M).EQ.1)GO TO 56
                     MSYM=ORBSYM(M)
                     MJSYM=IEOR(MSYM,JSYM)
                     MJ=ITR(M)+J
                     JM=ITR(J)+M
                     DO 178 C=1,NV
                        IF(FZV(C).EQ.1)GO TO 178
                        CSYM=ORBSYM(C+NO)
                        BSYM=IEOR(CSYM,MJSYM)+1
                        FB=FLOV(BSYM,3)-NO
                        LB=FLOV(BSYM,4)-NO
                        DO 173 B=FB,LB
                           ZLCB=ZLX(C,B)
                           MJCB=TOFF(M,J,ZLCB)+TADD(C,B)
                           KMBC=T3OFF(M,K,ZLCB,IASYM)+T3ADD(B,C,IASYM)
                           VAL1=VAL1+Y3(KMBC)*T2(MJCB)
  173                   CONTINUE
  178                CONTINUE
   56             CONTINUE
                  GEMAT(AIJK)=GEMAT(AIJK)+VAL1*FACJK
   46          CONTINUE
   36       CONTINUE
   26    CONTINUE
   16 CONTINUE
C
      I81=1
      DO 190 C=1,NV
      IF(FZV(C).EQ.1)GO TO 190
      CSYM=ORBSYM(C+NO)
      DO 185 K=1,NO
      IF(FZO(K).EQ.1)GO TO 185
      KSYM=ORBSYM(K)
      KCSYM=IEOR(KSYM,CSYM)+1
      DIM=NT3(KCSYM)
      CALL WREADW(ITAP81,W3,INTOWP(DIM),I81,I81)
      DO 180 B=1,NV
      IF(FZV(B).EQ.1)GO TO 180
      BSYM=ORBSYM(B+NO)
      CBSYM=IEOR(CSYM,BSYM)
      ZLCB=ZLX(C,B)
      CB=IOFF(MAX0(C,B))+MIN0(C,B)
      DO 175 I=1,NO
      IF(FZO(I).EQ.1)GO TO 175
      ISYM=ORBSYM(I)
      KISYM=IEOR(KSYM,ISYM)
      FACIK=1.0D0
      IF(I.EQ.K)FACIK=2.0D0
      DO 170 A=1,NV
      IF(FZV(A).EQ.1)GO TO 170
      AI=ITR(A)+I
      AK=ITR(A)+K
      ASYM=ORBSYM(A+NO)
      ACBSYM=IEOR(ASYM,CBSYM)
      MSYM=IEOR(ACBSYM,KISYM)
      FM=FLOV(MSYM+1,1)
      LM=FLOV(MSYM+1,2)
      DO 165 M=FM,LM
      JSYM=IEOR(MSYM,CBSYM)
      FJ=FLOV(JSYM+1,1)
      LJ=FLOV(JSYM+1,2)
      DO 160 J=FJ,LJ
      AJ=ITR(A)+J
      JK=ITR(MAX0(J,K))+MIN0(J,K)
      JI=ITR(MAX0(I,J))+MIN0(I,J)
      AIJK=EOFF(AI)+EADD(JK)
      AKJI=EOFF(AK)+EADD(JI)
      FACJK=1.0D0
      IF(J.EQ.K)FACJK=2.0D0
      FACIJ=1.0D0
      IF(J.EQ.I)FACIJ=2.0D0
      JMCB=TOFF(J,M,ZLCB)+TADD(C,B)
      ZLBA=ZLX(B,A)
      MIBA=T3OFF(M,I,ZLBA,KCSYM)+T3ADD(B,A,KCSYM)
      MIAB=T3OFF(I,M,ZLBA,KCSYM)+T3ADD(A,B,KCSYM)
      GEMAT(AIJK)=GEMAT(AIJK)-(W3(MIBA)+W3(MIBA)-W3(MIAB))
     .                       *Z2(JMCB)*FACJK
      GEMAT(AKJI)=GEMAT(AKJI)+ W3(MIBA)*Z2(JMCB)*FACIJ
  160 CONTINUE
  165 CONTINUE
  170 CONTINUE
  175 CONTINUE
  180 CONTINUE
  185 CONTINUE
  190 CONTINUE
C
      ENDIF
C
      DO 10 A=1,NV
         IF(FZV(A).EQ.1)GO TO 10
         ASYM=ORBSYM(A+NO)
         DO 20 I=1,NO
            IF(FZO(I).EQ.1)GO TO 20
            ISYM=ORBSYM(I)
            AISYM=IEOR(ASYM,ISYM)
            AI=ITR(A)+I
            DO 30 J=1,NO
               IF(FZO(J).EQ.1)GO TO 30
               JSYM=ORBSYM(J)
               KSYM=IEOR(AISYM,JSYM)
               FK=FLOV(KSYM+1,1)
               LK=FLOV(KSYM+1,2)
               FB=FLOV(KSYM+1,3)-NO
               LB=FLOV(KSYM+1,4)-NO
               DO 40 K=FK,LK
                  JK=ITR(MAX0(J,K))+MIN0(J,K)
                  FACJK=1.0D+00
                  IF(J.EQ.K) FACJK=2.0D+00
                  AIJK=EOFF(AI)+EADD(JK)
                  VAL1=0.0D+00
                  VAL2=0.0D+00
                  DO 50 B=FB,LB
                     ZLAB=ZLX(A,B)
                     IJAB=TOFF(I,J,ZLAB)+TADD(A,B)
                     JIAB=TOFF(J,I,ZLAB)+TADD(A,B)
                     VAL1=VAL1+(TAU(IJAB)+TAU(IJAB)-TAU(JIAB))*Z1(K,B)
                     VAL2=VAL2+T1(K,B)*Z2(IJAB)
   50             CONTINUE
                  GEMAT(AIJK)=GEMAT(AIJK)-(VAL1+VAL2)*FACJK
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
         FI=FLOV(ASYM+1,1)
         TI=FLOV(ASYM+1,2)
         DO 23 I=FI,TI
            AI=ITR(A)+I
            VAL1=0.0D+00
            VAL2=0.0D+00
            DO 82 L=1,NO
               IF(FZO(L).EQ.1)GO TO 82
               LSYM=ORBSYM(L)
               FD=FLOV(LSYM+1,3)-NO
               TD=FLOV(LSYM+1,4)-NO
               DO 83 D=FD,TD
                  ZLAD=ZLX(A,D)
                  ILAD=TOFF(I,L,ZLAD)+TADD(A,D)
                  LIAD=TOFF(L,I,ZLAD)+TADD(A,D)
                  VAL1=VAL1+T2(ILAD)*Z1(L,D)
                  VAL2=VAL2+TAU(LIAD)*Z1(L,D)
   83          CONTINUE
   82       CONTINUE
            DO 85 J=1,NO
               IF(FZO(J).EQ.1)GO TO 85
               JK=ITR(J)+J
               AIJK=EOFF(AI)+EADD(JK)
               GEMAT(AIJK)=GEMAT(AIJK)+VAL1+VAL1+VAL1+VAL1-VAL2-VAL2
     .                     +T1(I,A)+T1(I,A)+Z1(I,A)+Z1(I,A)
     .                     +VAL1+VAL1+VAL1+VAL1-VAL2-VAL2
     .                     +T1(I,A)+T1(I,A)+Z1(I,A)+Z1(I,A)
     .                     +T1(I,A)+T1(I,A)+T1(I,A)+T1(I,A)
   85       CONTINUE
   23    CONTINUE
C
         FJ=FLOV(ASYM+1,1)
         TJ=FLOV(ASYM+1,2)
         DO 310 J=FJ,TJ
            VAL1=0.0D+00
            VAL2=0.0D+00
            DO 311 L=1,NO
               IF(FZO(L).EQ.1)GO TO 311
               LSYM=ORBSYM(L)
               FD=FLOV(LSYM+1,3)-NO
               TD=FLOV(LSYM+1,4)-NO
               DO 312 D=FD,TD
                  ZLAD=ZLX(A,D)
                  JLAD=TOFF(J,L,ZLAD)+TADD(A,D)
                  LJAD=TOFF(L,J,ZLAD)+TADD(A,D)
                  VAL1=VAL1+T2(JLAD)*Z1(L,D)
                  VAL2=VAL2+TAU(LJAD)*Z1(L,D)
  312          CONTINUE
  311       CONTINUE
            DO 314 I=1,NO
               IF(FZO(I).EQ.1)GO TO 314
               AI=ITR(A)+I
               JK=ITR(MAX0(J,I))+MIN0(J,I)
               AIJK=EOFF(AI)+EADD(JK)
               FACJK=1.0D+00
               IF (I.EQ.J) FACJK=2.0D+00
               GEMAT(AIJK)=GEMAT(AIJK)-(VAL1+VAL1-VAL2
     .                  +T1(J,A)+T1(J,A)+Z1(J,A))*FACJK
  314       CONTINUE
  310    CONTINUE
   10 CONTINUE
C
C
      DO 11 J=1,NO
         IF(FZO(J).EQ.1)GO TO 11
         JSYM=ORBSYM(J)
         FK=FLOV(JSYM+1,1)
         TK=FLOV(JSYM+1,2)
         DO 21 K=FK,TK
            JK=ITR(MAX0(J,K))+MIN0(J,K)
            FACJK=1.0D+00
            IF (J.EQ.K) FACJK=2.0D+00
            VAL1=0.0D+00
            DO 31 L=1,NO
               IF(FZO(L).EQ.1)GO TO 31
               LSYM=ORBSYM(L)
               JLSYM=IEOR(JSYM,LSYM)
               DO 41 B=1,NV
                  IF(FZV(B).EQ.1)GO TO 41
                  BSYM=ORBSYM(B+NO)
                  CSYM=IEOR(JLSYM,BSYM)+1
                  FC=FLOV(CSYM,3)-NO
                  TC=FLOV(CSYM,4)-NO
                  DO 52 C=FC,TC
                     ZLBC=ZLX(B,C)
                     KLBC=TOFF(K,L,ZLBC)+TADD(B,C)
                     JLBC=TOFF(J,L,ZLBC)+TADD(B,C)
                     VAL1=VAL1+T2(KLBC)*Z2(JLBC)
   52             CONTINUE
   41          CONTINUE
   31       CONTINUE
            DO 32 I=1,NO
               IF(FZO(I).EQ.1)GO TO 32
               ISYM=ORBSYM(I)
               FA=FLOV(ISYM+1,3)-NO
               TA=FLOV(ISYM+1,4)-NO
               JI=ITR(MAX0(J,I))+MIN0(J,I)
               FACIJ=1.0D+00
               IF (I.EQ.J) FACIJ=2.0D+00
               DO 42 A=FA,TA
                  AI=ITR(A)+I
                  AK=ITR(A)+K
                  AIJK=EOFF(AI)+EADD(JK)
                  AKJI=EOFF(AK)+EADD(JI)
                  GEMAT(AIJK)=GEMAT(AIJK)-2.0D+00*VAL1*T1(I,A)*FACJK
                  GEMAT(AKJI)=GEMAT(AKJI)+VAL1*T1(I,A)*FACIJ
   42          CONTINUE
   32       CONTINUE
   21    CONTINUE
   11 CONTINUE
C
C
      DO 60 I=1,NO
         IF(FZO(I).EQ.1)GO TO 60
         ISYM=ORBSYM(I)
         DO 61 J=1,NO
            IF(FZO(J).EQ.1)GO TO 61
            JSYM=ORBSYM(J)
            IJSYM=IEOR(ISYM,JSYM)
            JI=ITR(MAX0(J,I))+MIN0(J,I)
            FACIJ=1.0D+00
            IF (I.EQ.J) FACIJ=2.0D+00
            DO 62 A=1,NV
               IF(FZV(A).EQ.1)GO TO 62
               ASYM=ORBSYM(A+NO)
               IASYM=IEOR(ISYM,ASYM)
               BSYM=IEOR(IJSYM,ASYM)+1
               FB=FLOV(BSYM,3)-NO
               TB=FLOV(BSYM,4)-NO
               AI=(A-1)*NO+I
               DO 63 B=FB,TB
                  Q1(B)=0.0D+00
                  Q2(B)=0.0D+00
                  Q3(B)=0.0D+00
                  DO 64 L=1,NO
                     IF(FZO(L).EQ.1)GO TO 64
                     LSYM=ORBSYM(L)
                     CSYM=IEOR(IASYM,LSYM)+1
                     FC=FLOV(CSYM,3)-NO
                     TC=FLOV(CSYM,4)-NO
                     DO 65 C=FC,TC
                        ZLAC=ZLX(A,C)
                        ZLBC=ZLX(B,C)
                        ILAC=TOFF(I,L,ZLAC)+TADD(A,C)
                        LIAC=TOFF(L,I,ZLAC)+TADD(A,C)
                        JLBC=TOFF(J,L,ZLBC)+TADD(B,C)
                        LJBC=TOFF(L,J,ZLBC)+TADD(B,C)
                        Q1(B)=Q1(B)+T2(ILAC)*Z2(JLBC)
                        Q2(B)=Q2(B)+T2(LIAC)*Z2(JLBC)
                        Q3(B)=Q3(B)+T2(LIAC)*Z2(LJBC)
   65                CONTINUE
   64             CONTINUE
   63          CONTINUE
C
               KSYM=IEOR(IJSYM,ASYM)+1
               FK=FLOV(KSYM,1)
               TK=FLOV(KSYM,2)
               FB=FLOV(KSYM,3)-NO
               TB=FLOV(KSYM,4)-NO
               DO 73 K=FK,TK
                  FACJK=1.0D+00
                  IF (J.EQ.K) FACJK=2.0D+00
                  AK=ITR(A)+K
                  JK=ITR(MAX0(J,K))+MIN0(J,K)
                  AIJK=EOFF(AI)+EADD(JK)
                  AKJI=EOFF(AK)+EADD(JI)
                  VAL1=0.0D+00
                  VAL2=0.0D+00
                  VAL3=0.0D+00
                  DO 74 B=FB,TB
                     VAL1=VAL1+Q1(B)*T1(K,B)
                     VAL2=VAL2+Q2(B)*T1(K,B)
                     VAL3=VAL3+Q3(B)*T1(K,B)
   74             CONTINUE
                  GEMAT(AIJK)=GEMAT(AIJK)-(VAL1+VAL1-VAL2)*FACJK
                  GEMAT(AKJI)=GEMAT(AKJI)+(VAL1+VAL3)*FACIJ
   73          CONTINUE
   62       CONTINUE
C
            DO 75 K=1,NO
               IF(FZO(K).EQ.1)GO TO 75
               KSYM=ORBSYM(K)
               JK=ITR(MAX0(J,K))+MIN0(J,K)
               FACJK=1.0D+00
               IF (J.EQ.K) FACJK=2.0D+00
               LSYM=IEOR(IJSYM,KSYM)+1
               FL=FLOV(LSYM,1)
               TL=FLOV(LSYM,2)
               DO 76 L=FL,TL
                  VALO(L)=0.0D+00
                  DO 77 B=1,NV
                     IF(FZV(B).EQ.1)GO TO 77
                     BSYM=ORBSYM(B+NO)
                     CSYM=IEOR(IJSYM,BSYM)+1
                     FC=FLOV(CSYM,3)-NO
                     TC=FLOV(CSYM,4)-NO
                     DO 78 C=FC,TC
                        ZLBC=ZLX(B,C)
                        LKBC=TOFF(L,K,ZLBC)+TADD(B,C)
                        IJBC=TOFF(I,J,ZLBC)+TADD(B,C)
                        VALO(L)=VALO(L)+TAU(IJBC)*Z2(LKBC)
   78                CONTINUE
   77             CONTINUE
   76          CONTINUE
               FA=FLOV(LSYM,3)-NO
               TA=FLOV(LSYM,4)-NO
               DO 80 A=FA,TA
                  AI=ITR(A)+I
                  AIJK=EOFF(AI)+EADD(JK)
                  VAL4=0.0D+00
                  DO 81 L=FL,TL
                     VAL4=VAL4+VALO(L)*T1(L,A)
   81             CONTINUE
                  GEMAT(AIJK)=GEMAT(AIJK)+VAL4*FACJK
   80          CONTINUE
   75       CONTINUE
   61    CONTINUE
   60 CONTINUE
C
C
      DO 87 I=1,NO
         IF(FZO(I).EQ.1)GO TO 87
         ISYM=ORBSYM(I)
         FU=FLOV(ISYM+1,1)
         TU=FLOV(ISYM+1,2)
         DO 88 U=FU,TU
            VALO(U)=0.0D+00
            DO 89 L=1,NO
               IF(FZO(L).EQ.1)GO TO 89
               LSYM=ORBSYM(L)
               ILSYM=IEOR(ISYM,LSYM)
               DO 90 BE=1,NV
                  IF(FZV(BE).EQ.1)GO TO 90
                  BESYM=ORBSYM(BE+NO)
                  DSYM=IEOR(ILSYM,BESYM)+1
                  FD=FLOV(DSYM,3)-NO
                  TD=FLOV(DSYM,4)-NO
                  DO 91 D=FD,TD
                     ZLBED=ZLX(BE,D)
                     ILBED=TOFF(I,L,ZLBED)+TADD(BE,D)
                     ULBED=TOFF(U,L,ZLBED)+TADD(BE,D)
                     VALO(U)=VALO(U)+T2(ILBED)*Z2(ULBED)
   91             CONTINUE
   90          CONTINUE
   89       CONTINUE
   88    CONTINUE
         FA=FLOV(ISYM+1,3)-NO
         TA=FLOV(ISYM+1,4)-NO
         DO 93 A=FA,TA
            AI=ITR(A)+I
            VAL1=0.0D+00
            DO 94 U=FU,TU
               VAL1=VAL1+VALO(U)*T1(U,A)
   94       CONTINUE
            DO 95 J=1,NO
               IF(FZO(J).EQ.1)GO TO 95
               JK=ITR(J)+J
               IJ=ITR(MAX0(I,J))+MIN0(I,J)
               FACIJ=1.0D+00
               IF (I.EQ.J) FACIJ=2.0D+00
               AK=ITR(A)+J
               AIJK=EOFF(AI)+EADD(JK)
               AKIJ=EOFF(AK)+EADD(IJ)
               GEMAT(AIJK)=GEMAT(AIJK)-VAL1-VAL1-VAL1-VAL1
               GEMAT(AKIJ)=GEMAT(AKIJ)+VAL1*FACIJ
   95       CONTINUE
   93    CONTINUE
   87 CONTINUE
C
C
      DO 96 A=1,NV
         IF(FZV(A).EQ.1)GO TO 96
         ASYM=ORBSYM(A+NO)
         FBE=FLOV(ASYM+1,3)-NO
         LBE=FLOV(ASYM+1,4)-NO
         DO 97 BE=FBE,LBE
            Q1(BE)=0.0D+00
            DO 98 U=1,NO
               IF(FZO(U).EQ.1)GO TO 98
               USYM=ORBSYM(U)
               DO 99 L=1,NO
                  IF(FZO(L).EQ.1)GO TO 99
                  LSYM=ORBSYM(L)
                  ULSYM=IEOR(USYM,LSYM)
                  DSYM=IEOR(ASYM,ULSYM)+1
                  FD=FLOV(DSYM,3)-NO
                  TD=FLOV(DSYM,4)-NO
                  DO 110 D=FD,TD
                     ZLAD=ZLX(A,D)
                     ZLBED=ZLX(BE,D)
                     ULAD=TOFF(U,L,ZLAD )+TADD(A,D)
                     ULBED=TOFF(U,L,ZLBED)+TADD(BE,D)
                     Q1(BE)=Q1(BE)+T2(ULAD)*Z2(ULBED)
  110             CONTINUE
   99          CONTINUE
   98       CONTINUE
   97    CONTINUE
C
         FI=FLOV(ASYM+1,1)
         TI=FLOV(ASYM+1,2)
         DO 119 I=FI,TI
            AI=ITR(A)+I
            VAL1=0.0D+00
            DO 120 BE=FBE,LBE
               VAL1=VAL1+Q1(BE)*T1(I,BE)
  120       CONTINUE
            DO 121 J=1,NO
               IF(FZO(J).EQ.1)GO TO 121
               JK=ITR(J)+J
               AK=ITR(A)+J
               IJ=ITR(MAX0(I,J))+MIN0(I,J)
               FACIJ=1.0D+00
               IF (I.EQ.J) FACIJ=2.0D+00
               AIJK=EOFF(AI)+EADD(JK)
               AKIJ=EOFF(AK)+EADD(IJ)
               GEMAT(AIJK)=GEMAT(AIJK)-VAL1-VAL1-VAL1-VAL1
               GEMAT(AKIJ)=GEMAT(AKIJ)+VAL1*FACIJ
  121       CONTINUE
  119    CONTINUE
   96 CONTINUE
C
C     CALL WRTEM (GEMAT,ITR,ORBSYM,FLOV,NIRRED,NO,NV,EOFF,EADD)
C
      RETURN
      END
