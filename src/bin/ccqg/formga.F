      SUBROUTINE FORMGA(GAMAT,TAU,T2,T1,Z2,Z1,NO,NV,NO2,NV2,NTRO,NTRV,
     .                  IOFF,NTR,JOUT,ZLX,FLOV,ORBSYM,NIRRED,
     .                  Q1,TZ,OPTION,MOO,NDIMT2,TOFF,TADD,
     .                  AOFF,AADD,ITR,NSGOO,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 OPTION
      REAL*8 VAL,VAL2,FAC,FAC1,FACIJ,FACKL,FACTR
      REAL*8 GAMAT(NSGOO),T1(NO,NV),Z1(NO,NV),
     .       Q1(NO,NO),TZ(NO,NO),MOO(NO,NO),
     .       T2(NDIMT2),Z2(NDIMT2),TAU(NDIMT2)
      INTEGER IOFF(NTR),ZLX(NV,NV),FLOV(NIRRED,4),ORBSYM(NO+NV),
     .        TOFF(NO,NO,2),TADD(NV,NV),AOFF(NO2),AADD(NO2),ITR(NO),
     .        FZO(NO),FZV(NV)
C
      CALL ZERO(GAMAT,NSGOO)
      CALL ZERO(Q1,NO*NO)
      CALL ZERO(TZ,NO*NO)
C
      DO 10 I=1,NO
         IF(FZO(I).EQ.1)GO TO 10
         ISYM=ORBSYM(I)
         DO 20 J=1,NO
            IF(FZO(J).EQ.1)GO TO 20
            JSYM=ORBSYM(J)
            IJSYM=IEOR(ISYM,JSYM)
            FACIJ=1.0D+00
            IF (I.EQ.J) FACIJ=2.0D+00
            IJ=ITR(MAX0(I,J))+MIN0(I,J)
            DO 30 K=1,NO
               IF(FZO(K).EQ.1)GO TO 30
               KSYM=ORBSYM(K)
               IKSYM=IEOR(ISYM,KSYM)
               LSYM=IEOR(KSYM,IJSYM)+1
               FL=FLOV(LSYM,1)
               LL=FLOV(LSYM,2)
               DO 40 L=FL,LL
                  FACKL=1.0D+00
                  IF (K.EQ.L) FACKL=2.0D+00
                  KL=ITR(MAX0(K,L))+MIN0(K,L)
                  IF (KL.GT.IJ) GOTO 40
                  IJKL=AOFF(IJ)+AADD(KL)
                  VAL=0.0D+00
                  DO 50 A=1,NV
                     IF(FZV(A).EQ.1)GO TO 50
                     ASYM=ORBSYM(A+NO)
                     BSYM=IEOR(IKSYM,ASYM)+1
                     FB=FLOV(BSYM,3)-NO
                     LB=FLOV(BSYM,4)-NO
                     DO 60 B=FB,LB
                        ZL=ZLX(A,B)
                        IKAB=TOFF(I,K,ZL)+TADD(A,B)
                        JLAB=TOFF(J,L,ZL)+TADD(A,B)
                        VAL=VAL+TAU(IKAB)*Z2(JLAB)
   60                CONTINUE
   50             CONTINUE
                  GAMAT(IJKL)=GAMAT(IJKL)+VAL*FACIJ*FACKL
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      DO 330 I=1,NO
         IF(FZO(I).EQ.1)GO TO 330
         ISYM=ORBSYM(I)
         FJ=FLOV(ISYM+1,1)
         LJ=FLOV(ISYM+1,2)
         IF (LJ.GT.I) LJ=I
         DO 320 J=FJ,LJ
            VAL=0.0D+00
         IF(OPTION.EQ.'SDT1'.or.option.eq.'CCT ') VAL=MOO(I,J)+MOO(J,I)
            VAL2=0.0D+00
            DO 31 A=1,NV
               IF(FZV(A).EQ.1)GO TO 31
               ASYM=ORBSYM(A+NO)
               IASYM=IEOR(ISYM,ASYM)
               VAL=VAL+T1(I,A)*Z1(J,A)
               VAL=VAL+T1(J,A)*Z1(I,A)
               DO 62 M=1,NO
                  IF(FZO(M).EQ.1)GO TO 62
                  MSYM=ORBSYM(M)
                  CSYM=IEOR(MSYM,IASYM)+1
                  FC=FLOV(CSYM,3)-NO
                  LC=FLOV(CSYM,4)-NO
                  DO 70 C=FC,LC
                     ZL=ZLX(A,C)
                     IMAC=TOFF(I,M,ZL)+TADD(A,C)
                     JMAC=TOFF(J,M,ZL)+TADD(A,C)
                     VAL2=VAL2+T2(IMAC)*Z2(JMAC)
                     VAL2=VAL2+T2(JMAC)*Z2(IMAC)
   70             CONTINUE
   62          CONTINUE
   31       CONTINUE
            TZ(I,J)=VAL
            TZ(J,I)=VAL
            Q1(I,J)=VAL2
            Q1(J,I)=VAL2
  320    CONTINUE
  330 CONTINUE
C
      DO 11 I=1,NO
         IF(FZO(I).EQ.1)GO TO 11
         ISYM=ORBSYM(I)
         FJ=FLOV(ISYM+1,1)
         LJ=FLOV(ISYM+1,2)
         IF (LJ.GT.I) LJ=I
         DO 21 J=FJ,LJ
            VAL=4.D0*(TZ(I,J)+Q1(I,J))
            IJ=ITR(I)+J
            DO 41 K=1,I
               IF(FZO(K).EQ.1)GO TO 41
               KSYM=ORBSYM(K)
               FL=FLOV(KSYM+1,1)
               LL=FLOV(KSYM+1,2)
               K1=K
               IF(K.EQ.I)K1=J
               IF (LL.GT.K1) LL=K1
               DO 42 L=FL,LL
                  IF(K.NE.L)GO TO 42
                  KL=ITR(K)+L
                  IJKL=AOFF(IJ)+AADD(KL)
                  GAMAT(IJKL)=GAMAT(IJKL)-VAL
   42          CONTINUE
   41       CONTINUE
   21    CONTINUE
   11 CONTINUE
C
      DO 911 I=1,NO
         IF(FZO(I).EQ.1)GO TO 911
         IJ=ITR(I)+I
         DO 941 K=1,I
            IF(FZO(K).EQ.1)GO TO 941
            KSYM=ORBSYM(K)
            FL=FLOV(KSYM+1,1)
            LL=FLOV(KSYM+1,2)
            K1=K
            IF(K.EQ.I)K1=I
            IF (LL.GT.K1) LL=K1
            DO 942 L=FL,LL
               KL=ITR(K)+L
               IJKL=AOFF(IJ)+AADD(KL)
               GAMAT(IJKL)=GAMAT(IJKL)-4.D0*(TZ(K,L)+Q1(K,L))
  942       CONTINUE
  941    CONTINUE
  911 CONTINUE
C
      DO 811 I=1,NO
         IF(FZO(I).EQ.1)GO TO 811
         ISYM=ORBSYM(I)
         FK=FLOV(ISYM+1,1)
         LK=FLOV(ISYM+1,2)
         IF(LK.GT.I)LK=I
         DO 821 J=1,I
            IF(FZO(J).EQ.1)GO TO 821
            IJ=ITR(I)+J
            DO 841 K=FK,LK
               VAL=TZ(I,K)+Q1(I,K)
               K1=K
               IF(K.EQ.I)K1=J
               DO 842 L=1,K1
                  IF(FZO(L).EQ.1)GO TO 842
                  IF(L.NE.J)GO TO 842
                  KL=ITR(K)+L
                  IJKL=AOFF(IJ)+AADD(KL)
                  GAMAT(IJKL)=GAMAT(IJKL)+VAL
  842          CONTINUE
  841       CONTINUE
  821    CONTINUE
  811 CONTINUE
C
      DO 711 I=1,NO
         IF(FZO(I).EQ.1)GO TO 711
         DO 721 J=1,I
            IF(FZO(J).EQ.1)GO TO 721
            JSYM=ORBSYM(J)
            FK=FLOV(JSYM+1,1)
            LK=FLOV(JSYM+1,2)
            IF(LK.GT.I)LK=I
            IJ=ITR(I)+J
            DO 741 K=FK,LK
               VAL=TZ(J,K)+Q1(J,K)
               K1=K
               IF(K.EQ.I)K1=J
               DO 742 L=1,K1
                  IF(FZO(L).EQ.1)GO TO 742
                  IF(L.NE.I)GO TO 742
                  KL=ITR(K)+L
                  IJKL=AOFF(IJ)+AADD(KL)
                  GAMAT(IJKL)=GAMAT(IJKL)+VAL
  742          CONTINUE
  741       CONTINUE
  721    CONTINUE
  711 CONTINUE
C
      DO 611 I=1,NO
         IF(FZO(I).EQ.1)GO TO 611
         ISYM=ORBSYM(I)
         DO 621 J=1,I
            IF(FZO(J).EQ.1)GO TO 621
            IJ=ITR(I)+J
            DO 641 K=1,I
               IF(FZO(K).EQ.1)GO TO 641
               IF(K.NE.J)GO TO 641
               FL=FLOV(ISYM+1,1)
               LL=FLOV(ISYM+1,2)
               IF(LL.GT.K)LL=K
               IF(K.EQ.I)LL=J
               DO 642 L=FL,LL
                  KL=ITR(K)+L
                  IJKL=AOFF(IJ)+AADD(KL)
                  GAMAT(IJKL)=GAMAT(IJKL)+TZ(I,L)+Q1(I,L)
  642          CONTINUE
  641       CONTINUE
  621    CONTINUE
  611 CONTINUE
C
      DO 511 I=1,NO
         IF(FZO(I).EQ.1)GO TO 511
         DO 521 J=1,I
            IF(FZO(J).EQ.1)GO TO 521
            JSYM=ORBSYM(J)
            IJ=ITR(I)+J
            DO 541 K=1,I
               IF(FZO(K).EQ.1)GO TO 541
               IF(K.NE.I)GO TO 541
               FL=FLOV(JSYM+1,1)
               LL=FLOV(JSYM+1,2)
               IF(LL.GT.J)LL=J
               DO 542 L=FL,LL
                  IF(FZO(L).EQ.1)GO TO 542
                  KL=ITR(K)+L
                  IJKL=AOFF(IJ)+AADD(KL)
                  GAMAT(IJKL)=GAMAT(IJKL)+TZ(J,L)+Q1(J,L)
  542          CONTINUE
  541       CONTINUE
  521    CONTINUE
  511 CONTINUE
C
    1 FORMAT (' GAMAT(',I3,',',I3,',',I3,',',I3,')=',F20.12)
      RETURN
      END
