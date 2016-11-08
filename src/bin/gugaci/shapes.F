      SUBROUTINE shapes(INT,C,S,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 sqrt2,sqrt3,sqt1p5,asqrt2,hsqrt3
      INTEGER ARR,TR1,TR2,ASM,AOS,OS,WTW,WTX,WTY,WAB,SS,SSI,SSJ,SYMORB
      INTEGER BMAX,ORBFRM
      REAL*8 INT(NMAX),C(NWKS),S(NWKS)
      DIMENSION KADD(SYMORB),LADD(SYMORB),IJADD(NUMIJ),WTW(ORBFRM,NSYM)
      DIMENSION WTX(ORBFRM,NSYM),WTY(ORBFRM),WAB(ORBFRM),SS(NORBS)
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /INTS/   NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
C  UNIVERSAL IDENTITY OF THE OBJECTS IN THESE COMMON
C     COMMON /ALL/ ARR,VAL1,VAL2,VAL3,ITR1,ITR2,IA,JA,ITYPE,ISEGT
C    *,LVFRM1,NLWKI,NLWKJ,IMAX,IMIN
C     COMMON /SYMM/ JSM,JFSYM,IFSYM,MAXSYM(8),MINSYM(8),ISMOFF(8)
C    #,             NUMSYM(8)
      COMMON /SYMM/ ASM,JS,IS,MX(8),MN(8),OS(8),NUMSYM(8)
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP99,ITAP94
     *,             ITAPE3,ITAP95,ITAP96
      COMMON /ALL/VAL1,VAL2,VAL3,ARR,TR1,TR2,IA,JA,M,ISEG,N,N1,N2
     *,           IMAX,IMIN
      COMMON /MINMAX/ IMING,IMAXG,JMING,JMAXG
      common /sqrts/ sqrt2,sqrt3,sqt1p5,asqrt2,hsqrt3
C
      REAL*8 Z,VAL1,VAL2,VAL3
cets111090      SQRT2=DSQRT(2.0D+00)
cets111090      SQRT3=DSQRT(3.0D+00)
cets111090      SQT1P5=DSQRT(1.5D+00)
cets110990      RETURN
C
C
cets110990      ENTRY SHAPES
C
C
      AOS=OS(ASM+1)
      GO TO (999,999,304,999,999,310,307,308,303,306,309,302,301,314
     1,314,314,317,318,319),M
C
      WRITE(ITAPE6,*) ' FUNNY M'
      CALL MABORT
C
C     ----- YZ NO 1 ISEG=22 -----
C***********************************************************************
  301 IF (ISEG.EQ.16) GO TO 3011
      IF (ISEG.EQ.18) GO TO 3012
      IF (ISEG.NE.22) GO TO 999
      IF (JS.NE.0) RETURN
      JJ=JA
      M1=MX(IS+1)
      M2=MN(IS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      II=IA+WTY(M2)
      CJJ=C(JJ)
      SJJ=0.0
      DO 152 I=M2,M1
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
         S(II)=S(II)+Z*CJJ
         SJJ=SJJ+Z*C(II)
  152 II=II+1
      S(JJ)=S(JJ)+SJJ
      RETURN
C
C    YZ   N0 1  ISEG=16
C
 3011 JJ=JA
      IF (JS.NE.0) RETURN
      M1=MX(IS+1)
      M2=MN(IS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      II=IA+WTY(M2)
      CJJ=C(JJ)
      SJJ=0.0
      DO 151 I=M2,M1
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*(INT(LAD+3)+INT(LAD+1))
         S(II)=S(II)+Z*CJJ
         SJJ=SJJ+Z*C(II)
  151 II=II+1
      S(JJ)=S(JJ)+SJJ
      RETURN
C
C  YZ  NO 1  ISEG=18
C
 3012 IF (JS.NE.0) RETURN
      M1=MX(IS+1)
      M2=MN(IS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      II=IA+WTY(M2)
      JJ=JA
      CJJ=C(JJ)
      SJJ=0.0
      DO 150 I=M2,M1
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*INT(LAD+3)
         S(II)=S(II)+Z*CJJ
         SJJ=SJJ+Z*C(II)
  150 II=II+1
      S(JJ)=S(JJ)+SJJ
      RETURN
C
C   XZ  NO 2 ISEG=13
C***********************************************************************
  302 IF (ISEG.NE.13) GO TO 999
      IF (JS.NE.0) RETURN
      JJ=JA
      CJJ=C(JJ)
      SJJ=0.0
      IF (N.LT.2) RETURN
      DO 250 I=2,N
         SSI=IEOR(SS(I),IS)
         IF (SSI.GT.SS(I)) GO TO 250
         M1=MX(SSI+1)
         M2=MN(SSI+1)
         IF (M2.GT.N) GO TO 250
         IF (SSI.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 250
         II=IA+WTX(I,IS+1)+WTY(M2)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         DO 251 J=M2,M1
            LAD1=LAD+LADD(J+LOS)
            Z=VAL1*(INT(LAD1+1)-INT(LAD1+3))
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  251    II=II+1
  250 CONTINUE
      S(JJ)=S(JJ)+SJJ
      RETURN
C
C   XY   NO 3 ISEG=18
C***********************************************************************
  303 CONTINUE
      IF (N.LT.2) RETURN
      IF (ISEG.EQ.16) GO TO 3031
      IF (ISEG.EQ.3) GO TO 3033
      IF (ISEG.EQ.22) GO TO 3034
      IF (ISEG.NE.18) GO TO 999
      SSI=IEOR(JS,IS)
      IF (SSI.GT.JS) GO TO 6601
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6601
      IF (M1.GT.N) M1=N
      M3=MX(SSI+1)
      M4=MN(SSI+1)
      IF (M4.GT.N) GO TO 6601
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6601
      JJ=JA+WTY(M2)
      DO 365 I=M2,M1
         IF (SSI.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 365
         II=IA+WTX(I,IS+1)+WTY(M4)
         CJJ=C(JJ)
         SJJ=0.0
         DO 366 J=M4,M3
            LAD=ARR+LADD(J+AOS)
            Z=-VAL1*INT(LAD+3)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  366    II=II+1
         S(JJ)=S(JJ)+SJJ
  365 JJ=JJ+1
C
C  NO 150
C
 6601 SSI=IEOR(JS,IS)
      IF (JS.GT.SSI) RETURN
      M1=MX(SSI+1)
      M2=MN(SSI+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) RETURN
      DO 5352 I=M2,M1
         IF (JS.EQ.SSI) M3=I-1
         IF (M3.LT.M4) GO TO 5352
         II=IA+WTX(I,IS+1)+WTY(M4)
         JJ=JA+WTY(M4)
         Z=VAL1*INT(ARR+LADD(I+AOS)+3)
         DO 5353 J=M4,M3
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
            II=II+1
 5353    JJ=JJ+1
 5352 CONTINUE
      RETURN
C
C  XY  NO 3 ISEG=16
C
 3031 SSI=IEOR(JS,IS)
      IF (SSI.GT.JS) GO TO 6602
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6602
      IF (M1.GT.N) M1=N
      M3=MX(SSI+1)
      M4=MN(SSI+1)
      IF (M4.GT.N) GO TO 6602
      IF (M2.LT.2) M2=2
      JJ=JA+WTY(M2)
      DO 363 I=M2,M1
         CJJ=C(JJ)
         SJJ=0.0
         IF (SSI.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 363
         II=IA+WTX(I,IS+1)+WTY(M4)
         DO 364 J=M4,M3
            LAD=ARR+LADD(J+AOS)
            Z=-VAL1*(INT(LAD+3)+INT(LAD+1))
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  364    II=II+1
         S(JJ)=S(JJ)+SJJ
  363 JJ=JJ+1
C
C  NO 150
C
 6602 SSI=IEOR(JS,IS)
      IF (JS.GT.SSI) RETURN
      M1=MX(SSI+1)
      M2=MN(SSI+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) RETURN
      DO 6352 I=M2,M1
         IF (JS.EQ.SSI) M3=I-1
         IF (M3.LT.M4) GO TO 6352
         II=IA+WTX(I,IS+1)+WTY(M4)
         JJ=JA+WTY(M4)
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*(INT(LAD+3)+INT(LAD+1))
         DO 6353 J=M4,M3
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
            II=II+1
 6353    JJ=JJ+1
 6352 CONTINUE
      RETURN
C
C  XY  NO 3 ISEG=22
C
 3034 SSI=IEOR(JS,IS)
      IF (SSI.GT.JS) GO TO 6603
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6603
      IF (M1.GT.N) M1=N
      M3=MX(SSI+1)
      M4=MN(SSI+1)
      IF (M4.GT.N) GO TO 6603
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6603
      JJ=JA+WTY(M2)
      DO 380 I=M2,M1
         IF (SSI.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 380
         II=IA+WTX(I,IS+1)+WTY(M4)
         CJJ=C(JJ)
         SJJ=0.0
         DO 381 J=M4,M3
            LAD=ARR+LADD(J+AOS)
            Z=-VAL1*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  381    II=II+1
         S(JJ)=S(JJ)+SJJ
  380 JJ=JJ+1
C
C NO 150
C
 6603 SSI=IEOR(JS,IS)
      IF (JS.GT.SSI) RETURN
      M1=MX(SSI+1)
      M2=MN(SSI+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) RETURN
      DO 7352 I=M2,M1
         IF (JS.EQ.SSI) M3=I-1
         IF (M3.LT.M4) GO TO 7352
         II=IA+WTX(I,IS+1)+WTY(M4)
         JJ=JA+WTY(M4)
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
         DO 7353 J=M4,M3
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
            II=II+1
 7353    JJ=JJ+1
 7352 CONTINUE
      RETURN
C
C  XY  NO 3  ISEG=3
C
 3033 SSI=IEOR(JS,IS)
      IF (SSI.GT.JS) GO TO 6604
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6604
      IF (M1.GT.N) M1=N
      M3=MX(SSI+1)
      M4=MN(SSI+1)
      IF (M4.GT.N) GO TO 6604
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6604
      JJ=JA+WTY(M2)
      DO 350 I=M2,M1
         IF (SSI.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 350
         II=IA+WTX(I,IS+1)+WTY(M4)
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)+KADD(I+KOS)
         CJJ=C(JJ)
         SJJ=0.0
         DO 351 J=M4,M3
            LAD1=LAD+LADD(J+AOS)
            Z=-VAL1*(INT(LAD1+2)-INT(LAD1+1))
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  351    II=II+1
         S(JJ)=S(JJ)+SJJ
  350 JJ=JJ+1
C
C  NO 4
C
 6604 SSI=IEOR(JS,IS)
      IF (JS.GT.SSI) GO TO 6605
      M1=MX(SSI+1)
      M2=MN(SSI+1)
      IF (M2.GT.N) GO TO 6605
      IF (M1.GT.N) M1=N
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) GO TO 6605
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6605
      DO 352 I=M2,M1
         IF (JS.EQ.SSI) M3=I-1
         IF (M3.LT.M4) GO TO 352
         JJ=JA+WTY(M4)
         II=IA+WTX(I,IS+1)+WTY(M4)
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)
         DO 353 J=M4,M3
            LSM=IEOR(KSM,SS(J))
            LOS=OS(LSM+1)
            LAD1=LAD+KADD(J+KOS)+LADD(J+LOS)
            Z=-VAL1*(INT(LAD1+1)-INT(LAD1+2))
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
            II=II+1
  353    JJ=JJ+1
  352 CONTINUE
C
C    NO   5
C
 6605 DO 354 I=3,N
         SSI=IEOR(SS(I),IS)
         IF (JS.GT.SSI) GO TO 354
         IF (SSI.GT.SS(I)) GO TO 354
         M1=MX(SSI+1)
         M2=MN(SSI+1)
         IF (M2.GT.N) GO TO 354
         IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
         IF (M2.LT.2) M2=2
         M3=MX(JS+1)
         M4=MN(JS+1)
         IF (M4.GT.N) GO TO 354
         II=IA+WTX(I,IS+1)+WTY(M2)
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)
         DO 355 J=M2,M1
            IF (JS.EQ.SSI) M3=J-1
            IF (M3.LT.M4) GO TO 355
            JJ=JA+WTY(M4)
            LSM=IEOR(KSM,SS(J))
            LOS=OS(LSM+1)
            LAD1=LAD+KADD(J+KOS)
            CII=C(II)
            SII=0.0
            DO 356 K=M4,M3
               LAD2=LAD1+LADD(K+LOS)
               Z=-VAL1*(INT(LAD2+1)-INT(LAD2+2))
               S(JJ)=S(JJ)+Z*CII
               SII=SII+Z*C(JJ)
  356       JJ=JJ+1
            S(II)=S(II)+SII
  355    II=II+1
  354 CONTINUE
C
C    NO   6
C
      DO 357 I=3,N
         SSI=IEOR(SS(I),IS)
         IF (JS.GT.SS(I)) GO TO 357
         IF (SSI.GT.JS) GO TO 357
         M1=MX(JS+1)
         M2=MN(JS+1)
         IF (M2.GT.N) GO TO 357
         IF (JS.EQ.SS(I).OR.M1.GT.N) M1=I-1
         IF (M2.LT.2) M2=2
         M3=MX(SSI+1)
         M4=MN(SSI+1)
         IF (M4.GT.N) GO TO 357
         JJ=JA+WTY(M2)
         IIA=IA+WTX(I,IS+1)
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)
         DO 358 J=M2,M1
            IF (SSI.EQ.JS) M3=J-1
            IF (M3.LT.M4) GO TO 358
            II=IIA+WTY(M4)
            LSM=IEOR(KSM,SS(J))
            LOS=OS(LSM+1)
            LAD1=LAD+KADD(J+KOS)
            CJJ=C(JJ)
            SJJ=0.0
            DO 359 K=M4,M3
               LAD2=LAD1+LADD(K+LOS)
               Z=-VAL1*(INT(LAD2+3)-INT(LAD2+2))
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
  359       II=II+1
            S(JJ)=S(JJ)+SJJ
  358    JJ=JJ+1
  357 CONTINUE
C
C   NO 7
C
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.3) M2=3
      IF (M1.LT.M2) RETURN
      JJ=JA+WTY(M2)
      DO 360 I=M2,M1
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)
         CJJ=C(JJ)
         SJJ=0.0
         DO 361 J=2,I-1
            SSI=IEOR(IS,SS(J))
            IF (SSI.GT.SS(J))GO TO 361
            M3=MX(SSI+1)
            M4=MN(SSI+1)
            IF (M4.GT.N) GO TO 361
            IF (SSI.EQ.SS(J)) M3=J-1
            IF (M3.LT.M4) GO TO 361
            II=IA+WTX(J,IS+1)+WTY(M4)
            LSM=IEOR(KSM,SS(J))
            LOS=OS(LSM+1)
            LAD1=LAD+KADD(J+KOS)
            DO 362 K=M4,M3
               LAD2=LAD1+LADD(K+LOS)
               Z=VAL1*(INT(LAD2+1)-INT(LAD2+3))
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
  362       II=II+1
  361    CONTINUE
         S(JJ)=S(JJ)+SJJ
         JJ=JJ+1
  360 CONTINUE
      RETURN
C
C   YY   NO 9  ISEG=4
C***********************************************************************
  304 IF (ISEG.EQ.5) GO TO 3043
      IF (ISEG.EQ.6) GO TO 3041
      IF (ISEG.EQ.8) GO TO 3042
      IF (ISEG.NE.4) GO TO 999
      VAL4=-VAL1/SQRT2
      IF (IS.GT.JS) GO TO 6606
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6606
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      M3=MX(IS+1)
      M4=MN(IS+1)
      IF (M4.GT.N) GO TO 6606
      JJ=JA+WTY(M2)
      DO 450 I=M2,M1
         IF (IS.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 450
         II=IA+WTY(M4)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         CJJ=C(JJ)
         SJJ=0.0
         DO 451 J=M4,M3
            LAD1=LAD+LADD(J+LOS)+1
            Z=VAL4*(INT(LAD1)-2.*INT(LAD1+1))
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  451    II=II+1
         S(JJ)=S(JJ)+SJJ
  450 JJ=JJ+1
C
C   NO 10
C
 6606 IF(IA.EQ.JA)RETURN
      IF (IS.NE.JS) RETURN
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      II=IA+WTY(M2)
      JJ=JA+WTY(M2)
      DO 460 I=M2,M1
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
         Z=VAL4*(INT(LAD+1)-2.*INT(LAD+2))
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
         II=II+1
  460 JJ=JJ+1
      RETURN
C
C YY  NO 9 ISEG=8
C
 3042 VAL4=VAL1*SQT1P5
      IF (IS.GT.JS) GO TO 6607
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6607
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6607
      M3=MX(IS+1)
      M4=MN(IS+1)
      IF (M4.GT.N) GO TO 6607
      JJ=JA+WTY(M2)
      DO 420 I=M2,M1
         IF (IS.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 420
         II=IA+WTY(M4)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         CJJ=C(JJ)
         SJJ=0.0
         DO 421 J=M4,M3
            LAD1=LAD+LADD(J+LOS)+TR1
            Z=VAL4*INT(LAD1)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  421    II=II+1
         S(JJ)=S(JJ)+SJJ
  420 JJ=JJ+1
C
C  NO 10
C
 6607 IF(IA.EQ.JA) RETURN
      IF (IS.NE.JS) GO TO 6598
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6598
      IF (M1.GT.N) M1=N
      II=IA+WTY(M2)
      JJ=JA+WTY(M2)
      DO 430 I=M2,M1
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
         Z=VAL4*INT(LAD+1)
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
         II=II+1
  430 JJ=JJ+1
C
C  NO 8
C
 6598 IF (JS.GT.IS) RETURN
      M1=MX(IS+1)
      M2=MN(IS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) RETURN
      II=IA+WTY(M2)
      DO 440 I=M2,M1
         IF (IS.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 440
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         JJ=JA+WTY(M4)
         DO 441 J=M4,M3
            LAD1=LAD+LADD(J+LOS)
            Z=VAL4*INT(LAD1+1)
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
  441    JJ=JJ+1
  440 II=II+1
      RETURN
C
C  YY  NO 10  ISEG=5
C
 3043 CONTINUE
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      II=IA+WTY(M2)
      JJ=JA+WTY(M2)
      XX=-VAL1/SQRT2+VAL3*SQT1P5
      IF (ABS(XX).LT.1.0D-06) GO TO 490
      VAL5=SQRT2*VAL1
      IF(IA.EQ.JA) GO TO 491
      IF (IS.NE.JS) GO TO 491
      DO 493 I=M2,M1
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
         Z=XX*INT(LAD+1)+VAL5*INT(LAD+2)
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
         II=II+1
  493 JJ=JJ+1
      GO TO 491
  490 VAL5=SQRT2*VAL1
      IF(IA.EQ.JA) GO TO 492
      IF (IS.NE.JS) GO TO 492
      DO 494 I=M2,M1
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
         Z=VAL5*INT(LAD+2)
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
         II=II+1
  494 JJ=JJ+1
      GO TO 492
C
C   NO   9
C
  491 IF (IS.GT.JS) RETURN
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      M3=MX(IS+1)
      M4=MN(IS+1)
      IF (M4.GT.N) RETURN
      JJ=JA+WTY(M2)
      DO 495 I=M2,M1
         IF (IS.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 495
         II=IA+WTY(M4)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         CJJ=C(JJ)
         SJJ=0.0
         DO 496 J=M4,M3
            LAD1=LAD+LADD(J+LOS)
            Z=XX*INT(LAD1+1)+VAL5*INT(LAD1+2)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  496    II=II+1
         S(JJ)=S(JJ)+SJJ
  495 JJ=JJ+1
      RETURN
C
  492 IF (IS.GT.JS) RETURN
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      M3=MX(IS+1)
      M4=MN(IS+1)
      IF (M4.GT.N) RETURN
      JJ=JA+WTY(M2)
      DO 497 I=M2,M1
         IF (IS.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 497
         II=IA+WTY(M4)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         CJJ=C(JJ)
         SJJ=0.0
         DO 498 J=M4,M3
            LAD1=LAD+LADD(J+LOS)+2
            Z=VAL5*INT(LAD1)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  498    II=II+1
         S(JJ)=S(JJ)+SJJ
  497 JJ=JJ+1
      RETURN
C
C   YY   NO 8  ISEG=6
C
 3041 XX=-VAL1/SQRT2+VAL3*SQT1P5
      IF (ABS(XX).LT.1.0D-06) GO TO 480
      VAL5=VAL1*SQRT2
      IF (JS.GT.IS) GO TO 481
      M1=MX(IS+1)
      M2=MN(IS+1)
      IF (M2.GT.N) GO TO 481
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 481
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) GO TO 481
      II=IA+WTY(M2)
      DO 472 I=M2,M1
         IF (IS.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 472
         JJ=JA+WTY(M4)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         CII=C(II)
         SII=0.0
         DO 473 J=M4,M3
            LAD1=LAD+LADD(J+LOS)
            Z=XX*INT(LAD1+1)+VAL5*INT(LAD1+2)
            SII=SII+Z*C(JJ)
            S(JJ)=S(JJ)+Z*CII
  473    JJ=JJ+1
         S(II)=S(II)+SII
  472 II=II+1
      GO TO 481
C
  480 VAL5=SQRT2*VAL1
      IF (JS.GT.IS) GO TO 483
      M1=MX(IS+1)
      M2=MN(IS+1)
      IF (M2.GT.N) GO TO 483
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 483
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) GO TO 483
      II=IA+WTY(M2)
      DO 474 I=M2,M1
         IF (IS.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 474
         JJ=JA+WTY(M4)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         CII=C(II)
         SII=0.0
         DO 475 J=M4,M3
            LAD1=LAD+LADD(J+LOS)
            Z=VAL5*INT(LAD1+2)
            SII=SII+Z*C(JJ)
            S(JJ)=S(JJ)+Z*CII
  475    JJ=JJ+1
         S(II)=S(II)+SII
  474 II=II+1
      GO TO 483
C
C    NO  10
C
  481 IF (IS.NE.JS) GO TO 482
      IF(IA.EQ.JA) GO TO 482
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 482
      IF (M1.GT.N) M1=N
      II=IA+WTY(M2)
      JJ=JA+WTY(M2)
      DO 470 I=M2,M1
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
         Z=XX*INT(LAD+1)+VAL5*INT(LAD+2)
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
         II=II+1
  470 JJ=JJ+1
      GO TO 482
C
  483 IF (IS.NE.JS) GO TO 484
      IF (IA.EQ.JA) GO TO 484
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 484
      IF (M1.GT.N) M1=N
      II=IA+WTY(M2)
      JJ=JA+WTY(M2)
      DO 471 I=M2,M1
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
         Z=VAL5*INT(LAD+2)
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
         II=II+1
  471 JJ=JJ+1
      GO TO 484
C
C  NO  9
C
  482 IF (IS.GT.JS) RETURN
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      M3=MX(IS+1)
      M4=MN(IS+1)
      IF (M4.GT.N) RETURN
      JJ=JA+WTY(M2)
      DO 476 I=M2,M1
         IF (IS.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 476
         II=IA+WTY(M4)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         CJJ=C(JJ)
         SJJ=0.0
         DO 477 J=M4,M3
            LAD1=LAD+LADD(J+LOS)
            Z=XX*INT(LAD1+3)+VAL5*INT(LAD1+2)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  477    II=II+1
         S(JJ)=S(JJ)+SJJ
  476 JJ=JJ+1
      RETURN
C
  484 IF (IS.GT.JS) RETURN
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      M3=MX(IS+1)
      M4=MN(IS+1)
      IF (M4.GT.N) RETURN
      JJ=JA+WTY(M2)
      DO 478 I=M2,M1
         IF (IS.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 478
         II=IA+WTY(M4)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         CJJ=C(JJ)
         SJJ=0.0
         DO 479 J=M4,M3
            LAD1=LAD+LADD(J+LOS)
            Z=VAL5*INT(LAD1+2)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  479    II=II+1
         S(JJ)=S(JJ)+SJJ
  478 JJ=JJ+1
      RETURN
C
C   WY  NO 27  ISEG=2
C***********************************************************************
  310 IF (ISEG.EQ.15) GO TO 3101
      IF (ISEG.EQ.17) GO TO 3102
      IF (ISEG.EQ.21) GO TO 3103
      IF (ISEG.NE.2) GO TO 999
C
C   NO 32
C
C
      VAL4=VAL1*SQRT2
      IF (IS.NE.0) GO TO 6587
      DO 1001 I=2,N
         IF (JS.GT.SS(I)) GO TO 1001
         II=IA+WAB(I)
         M1=MX(JS+1)
         M2=MN(JS+1)
         IF (M2.GT.N) GO TO 1001
         IF (JS.EQ.SS(I).OR.M1.GT.N) M1=I-1
         JJ=JA+WTY(M2)
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)+KADD(I+KOS)
         CII=C(II)
         SII=0.0
         DO 1002 J=M2,M1
            LAD1=LAD+LADD(J+AOS)
            Z=INT(LAD1+1)*VAL4
            SII=SII+Z*C(JJ)
            S(JJ)=S(JJ)+Z*CII
 1002    JJ=JJ+1
         S(II)=S(II)+SII
 1001 CONTINUE
C
C   NO  28
C
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6599
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6599
      JJ=JA+WTY(M2)
      DO 1003 I=M2,M1
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)
         CJJ=C(JJ)
         SJJ=0.0
         DO 1004 J=1,I-1
            II=IA+WAB(J)
            LSM=IEOR(KSM,SS(J))
            LOS=OS(LSM+1)
            LAD1=LAD+KADD(J+KOS)+LADD(J+LOS)
            Z=VAL4*INT(LAD1+1)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
 1004    CONTINUE
         S(JJ)=S(JJ)+SJJ
 1003 JJ=JJ+1
C
C   NO 30
C
 6587 M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6599
      IF (M1.GT.N) M1=N
      IF (M2.LT.3) M2=3
      IF (M1.LT.M2) GO TO 6599
      JJ=JA+WTY(M2)
      DO 1005 I=M2,M1
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)
         CJJ=C(JJ)
         SJJ=0.0
         DO 1006 J=2,I-1
            SSI=IEOR(IS,SS(J))
            IF (SSI.GT.SS(J))GO TO 1006
            M3=MX(SSI+1)
            M4=MN(SSI+1)
            IF (M4.GT.N) GO TO 1006
            IF (SSI.EQ.SS(J)) M3=J-1
            IF (M3.LT.M4) GO TO 1006
            II=IA+WTW(J,IS+1)+WTY(M4)
            LSM=IEOR(KSM,SS(J))
            LOS=OS(LSM+1)
            LAD1=LAD+KADD(J+KOS)
            DO 1007 K=M4,M3
               LAD2=LAD1+LADD(K+LOS)
               Z=VAL1*(INT(LAD2+1)+INT(LAD2+3))
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
 1007       II=II+1
 1006    CONTINUE
         S(JJ)=S(JJ)+SJJ
         JJ=JJ+1
 1005 CONTINUE
C
C   NO 32
C
 6599 IF(N.LT.3) GO TO 1018
      DO 1008 I=3,N
         SSI=IEOR(SS(I),IS)
         IF (JS.GT.SSI) GO TO 1008
         IF (SSI.GT.SS(I)) GO TO 1008
         M1=MX(SSI+1)
         M2=MN(SSI+1)
         IF (M2.GT.N) GO TO 1008
         IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
         IF (M2.LT.2) M2=2
         IF (M1.LT.M2) GO TO 1008
         M3=MX(JS+1)
         M4=MN(JS+1)
         IF (M4.GT.N) GO TO 1008
         II=IA+WTW(I,IS+1)+WTY(M2)
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)
         DO 1009 J=M2,M1
            IF (JS.EQ.SSI) M3=J-1
            IF (M3.LT.M4) GO TO 1009
            JJ=JA+WTY(M4)
            LSM=IEOR(KSM,SS(J))
            LOS=OS(LSM+1)
            LAD1=LAD+KADD(J+KOS)
            CII=C(II)
            SII=0.0
            DO 1010 K=M4,M3
               LAD2=LAD1+LADD(K+LOS)
               Z=VAL1*(INT(LAD2+1)+INT(LAD2+2))
               SII=SII+Z*C(JJ)
               S(JJ)=S(JJ)+Z*CII
 1010       JJ=JJ+1
            S(II)=S(II)+SII
 1009    II=II+1
 1008 CONTINUE
 1018 CONTINUE
C
C   NO 33
C
      SSI=IEOR(JS,IS)
      IF (SSI.GT.JS)GO TO  6608
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6608
      IF (M1.GT.N) M1=N
      M3=MX(SSI+1)
      M4=MN(SSI+1)
      IF (M4.GT.N) GO TO 6608
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6608
      JJ=JA+WTY(M2)
      DO 1011 I=M2,M1
         IF (SSI.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 1011
         II=IA+WTW(I,IS+1)+WTY(M4)
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)+KADD(I+KOS)
         CJJ=C(JJ)
         SJJ=0.0
         DO 1012 J=M4,M3
            LAD1=LAD+LADD(J+AOS)
            Z=VAL1*(INT(LAD1+2)+INT(LAD1+1))
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
 1012    II=II+1
         S(JJ)=S(JJ)+SJJ
 1011 JJ=JJ+1
C
C   NO 36
C
 6608 SSI=IEOR(JS,IS)
      IF (JS.GT.SSI) GO TO 6609
      M1=MX(SSI+1)
      M2=MN(SSI+1)
      IF (M2.GT.N) GO TO 6609
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) GO TO 6609
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6609
      DO 1013 I=M2,M1
         IF (JS.EQ.SSI) M3=I-1
         IF (M3.LT.M4) GO TO 1013
         JJ=JA+WTY(M4)
         II=IA+WTW(I,IS+1)+WTY(M4)
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)
         DO 1014 J=M4,M3
            LSM=IEOR(KSM,SS(J))
            LOS=OS(LSM+1)
            LAD1=LAD+KADD(J+KOS)+LADD(J+LOS)
            Z=VAL1*(INT(LAD1+1)+INT(LAD1+2))
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
            II=II+1
 1014    JJ=JJ+1
 1013 CONTINUE
C
C    NO 38
C
 6609 DO 1015 I=3,N
         SSI=IEOR(SS(I),IS)
         IF (JS.GT.SS(I)) GO TO 1015
         IF (SSI.GT.JS) GO TO 1015
         M1=MX(JS+1)
         M2=MN(JS+1)
         IF (M2.GT.N) GO TO 1015
         IF (JS.EQ.SS(I).OR.M1.GT.N) M1=I-1
         IF (M2.LT.2) M2=2
         M3=MX(SSI+1)
         M4=MN(SSI+1)
         IF (M4.GT.N) GO TO 1015
         JJ=JA+WTY(M2)
         IIA=IA+WTW(I,IS+1)
         KSM=IEOR(ASM,SS(I))
         KOS=OS(KSM+1)
         LAD=IJADD(ARR+I)
         DO 1016 J=M2,M1
            IF (SSI.EQ.JS) M3=J-1
            IF (M3.LT.M4) GO TO 1016
            II=IIA+WTY(M4)
            LSM=IEOR(KSM,SS(J))
            LOS=OS(LSM+1)
            LAD1=LAD+KADD(J+KOS)
            CJJ=C(JJ)
            SJJ=0.0
            DO 1017 K=M4,M3
               LAD2=LAD1+LADD(K+LOS)
               Z=VAL1*(INT(LAD2+3)+INT(LAD2+2))
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
 1017       II=II+1
            S(JJ)=S(JJ)+SJJ
 1016    JJ=JJ+1
 1015 CONTINUE
      RETURN
C
C   WY   NO 26   ISEG=21
C
 3103 VAL4=VAL1*SQRT2
      IF (IS.NE.0) GO TO 6610
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6610
      IF (M1.GT.N) M1=N
      JJ=JA+WTY(M2)
      DO 1020 I=M2,M1
         II=IA+WAB(I)
         LAD=ARR+LADD(I+AOS)
         Z=VAL4*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
 1020 JJ=JJ+1
C
C   NO 33
C
 6610 SSI=IEOR(JS,IS)
      IF (SSI.GT.JS) GO TO 6611
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6611
      IF (M1.GT.N) M1=N
      M3=MX(SSI+1)
      M4=MN(SSI+1)
      IF (M4.GT.N) GO TO 6611
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6611
      JJ=JA+WTY(M2)
      DO 1021 I=M2,M1
         IF (SSI.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 1021
         II=IA+WTW(I,IS+1)+WTY(M4)
         CJJ=C(JJ)
         SJJ=0.0
         DO 1022 J=M4,M3
            LAD=ARR+LADD(J+AOS)
            Z=VAL1*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
 1022    II=II+1
         S(JJ)=S(JJ)+SJJ
 1021 JJ=JJ+1
C
C NO 149
C
 6611 SSI=IEOR(JS,IS)
      IF (JS.GT.SSI) RETURN
      M1=MX(SSI+1)
      M2=MN(SSI+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) RETURN
      DO 6013 I=M2,M1
         IF (JS.EQ.SSI) M3=I-1
         IF (M3.LT.M4) GO TO 6013
         II=IA+WTW(I,IS+1)+WTY(M4)
         JJ=JA+WTY(M4)
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
         DO 6014 J=M4,M3
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
            II=II+1
 6014    JJ=JJ+1
 6013 CONTINUE
      RETURN
C
C    WY  NO 26   ISEG=15
C
 3101 VAL4=VAL1*SQRT2
      IF (IS.NE.0) GO TO 6612
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6612
      IF (M1.GT.N) M1=N
      JJ=JA+WTY(M2)
      DO 1030 I=M2,M1
         II=IA+WAB(I)
         LAD=ARR+LADD(I+AOS)
         Z=VAL4*(INT(LAD+1)+INT(LAD+2)+INT(LAD+3))
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
 1030 JJ=JJ+1
C
C    NO 33
C
 6612 SSI=IEOR(JS,IS)
      IF (SSI.GT.JS) GO TO 6613
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6613
      M3=MX(SSI+1)
      M4=MN(SSI+1)
      IF (M4.GT.N) GO TO 6613
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6613
      JJ=JA+WTY(M2)
      DO 1031 I=M2,M1
         IF (SSI.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 1031
         II=IA+WTW(I,IS+1)+WTY(M4)
         CJJ=C(JJ)
         SJJ=0.0
         DO 1032 J=M4,M3
            LAD=ARR+LADD(J+AOS)
            Z=VAL1*(INT(LAD+3)+INT(LAD+1))
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
 1032    II=II+1
         S(JJ)=S(JJ)+SJJ
 1031 JJ=JJ+1
C
C NO 149
C
 6613 SSI=IEOR(JS,IS)
      IF (JS.GT.SSI) RETURN
      M1=MX(SSI+1)
      M2=MN(SSI+1)
      IF (M2.GT.N) RETURN
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      DO 7013 I=M2,M1
         IF (JS.EQ.SSI) M3=I-1
         IF (M3.LT.M4) GO TO 7013
         II=IA+WTW(I,IS+1)+WTY(M4)
         JJ=JA+WTY(M4)
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*(INT(LAD+3)+INT(LAD+1))
         DO 7014 J=M4,M3
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
            II=II+1
 7014    JJ=JJ+1
 7013 CONTINUE
      RETURN
C
C   WY  NO 26 ISEG=17
C
 3102 VAL4=VAL1*SQRT2
      IF (IS.NE.0) GO TO 6614
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6614
      IF (M1.GT.N) M1=N
      JJ=JA+WTY(M2)
      DO 1023 I=M2,M1
         II=IA+WAB(I)
         LAD=ARR+LADD(I+AOS)
         Z=VAL4*(INT(LAD+3)+INT(LAD+2))
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
 1023 JJ=JJ+1
C
C   NO  33
C
 6614 SSI=IEOR(JS,IS)
      IF (SSI.GT.JS) GO TO 6615
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6615
      M3=MX(SSI+1)
      M4=MN(SSI+1)
      IF (M4.GT.N) GO TO 6615
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6615
      JJ=JA+WTY(M2)
      DO 1024 I=M2,M1
         IF (SSI.EQ.JS) M3=I-1
         IF (M3.LT.M4) GOTO 1024
         II=IA+WTW(I,IS+1)+WTY(M4)
         CJJ=C(JJ)
         SJJ=0.0
         DO 1025 J=M4,M3
            LAD=ARR+LADD(J+AOS)
            Z=VAL1*INT(LAD+3)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
 1025    II=II+1
         S(JJ)=S(JJ)+SJJ
 1024 JJ=JJ+1
C
C NO 149
C
 6615 SSI=IEOR(JS,IS)
      IF (JS.GT.SSI) RETURN
      M1=MX(SSI+1)
      M2=MN(SSI+1)
      IF (M2.GT.N) RETURN
      M3=MX(JS+1)
      M4=MN(JS+1)
      IF (M4.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      DO 5013 I=M2,M1
         IF (JS.EQ.SSI) M3=I-1
         IF (M3.LT.M4) GO TO 5013
         II=IA+WTW(I,IS+1)+WTY(M4)
         JJ=JA+WTY(M4)
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*INT(LAD+3)
         DO 5014 J=M4,M3
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
            II=II+1
 5014    JJ=JJ+1
 5013 CONTINUE
      RETURN
C
C   ZY   NO 17   ISEG=20
C***********************************************************************
  306 IF (ISEG.NE.20) GO TO 999
      IF (IS.NE.0) RETURN
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      JJ=JA+WTY(M2)
      II=IA
      CII=C(II)
      SII=0.0
      DO 601 I=M2,M1
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
         SII=SII+Z*C(JJ)
         S(JJ)=S(JJ)+Z*CII
  601 JJ=JJ+1
      S(II)=S(II)+SII
      RETURN
C
C     YX   NO  18   ISEG=20
C***********************************************************************
  307 IF (ISEG.NE.20) GO TO 999
      SSJ=IEOR(IS,JS)
      IF (SSJ.GT.IS) GO TO 6630
      M1=MX(IS+1)
      M2=MN(IS+1)
      IF (M2.GT.N) GO TO 6630
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6630
      II=IA+WTY(M2)
      M3=MX(SSJ+1)
      M4=MN(SSJ+1)
      IF (M4.GT.N) GO TO 6630
      DO 701 I=M2,M1
         IF (SSJ.EQ.IS) M3=I-1
         IF (M3.LT.M4) GO TO 701
         JJ=JA+WTX(I,JS+1)+WTY(M4)
         CII=C(II)
         SII=0.0
         DO 702 J=M4,M3
            LAD=ARR+LADD(J+AOS)
            Z=-VAL1*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
            SII=SII+Z*C(JJ)
            S(JJ)=S(JJ)+Z*CII
  702    JJ=JJ+1
         S(II)=S(II)+SII
  701 II=II+1
C
C NO 136
C
 6630 SSJ=IEOR(IS,JS)
      IF (IS.GT.SSJ) RETURN
      M1=MX(SSJ+1)
      M2=MN(SSJ+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      M3=MX(IS+1)
      M4=MN(IS+1)
      IF (M4.GT.N) RETURN
      DO 9080 I=M2,M1
         IF (IS.EQ.SSJ) M3=I-1
         IF (M3.LT.M4) GO TO 9080
         II=IA+WTY(M4)
         JJ=JA+WTX(I,JS+1)+WTY(M4)
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
         DO 9081 J=M4,M3
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
            II=II+1
 9081    JJ=JJ+1
 9080 CONTINUE
      RETURN
C
C    WZ   NO 23  ISEG=11
C***********************************************************************
  309 IF (ISEG.EQ.10) GO TO 3091
      IF (ISEG.EQ.14) RETURN
      IF (ISEG.NE.11) GO TO 999
      VAL4=VAL1*SQRT2
      IF (IS.NE.0.OR.JS.NE.0) GO TO 6631
      JJ=JA
      CJJ=C(JJ)
      SJJ=0.0
      DO 901 I=1,N
         II=IA+WAB(I)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)+LADD(I+LOS)
         Z=VAL4*INT(LAD+1)
         S(II)=S(II)+Z*CJJ
         SJJ=SJJ+Z*C(II)
  901 CONTINUE
      S(JJ)=S(JJ)+SJJ
C
C     NO 25
C
 6631 IF (JS.NE.0) RETURN
      JJ=JA
      CJJ=C(JJ)
      SJJ=0.0
      DO 902 I=2,N
         SSI=IEOR(SS(I),IS)
         IF (SSI.GT.SS(I)) GO TO 902
         M1=MX(SSI+1)
         M2=MN(SSI+1)
         IF (M2.GT.N) GO TO 902
         IF (SSI.EQ.SS(I)) M1=I-1
         II=IA+WTW(I,IS+1)+WTY(M2)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         DO 903 J=M2,M1
            LAD1=LAD+LADD(J+LOS)
            Z=VAL1*(INT(LAD1+1)+INT(LAD1+3))
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  903    II=II+1
  902 CONTINUE
      S(JJ)=S(JJ)+SJJ
      RETURN
C
C   YW   NO 19   ISEG=19
C***********************************************************************
  308 IF (ISEG.NE.19) GO TO 999
      VAL4=VAL1*SQRT2
      IF (JS.NE.0) GO TO 6632
      M1=MX(IS+1)
      M2=MN(IS+1)
      IF (M2.GT.N) GO TO 6632
      IF (M1.GT.N) M1=N
      II=IA+WTY(M2)
      DO 801 I=M2,M1
         JJ=JA+WAB(I)
         LAD=ARR+LADD(I+AOS)
         Z=VAL4*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
         II=II+1
  801 CONTINUE
C
C      NO 19A
C
 6632 SSJ=IEOR(IS,JS)
      IF (SSJ.GT.IS) GO TO 6633
      M1=MX(IS+1)
      M2=MN(IS+1)
      IF (M2.GT.N) GO TO 6633
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) GO TO 6633
      II=IA+WTY(M2)
      M3=MX(SSJ+1)
      M4=MN(SSJ+1)
      IF (M4.GT.N) GO TO 6633
      DO 802 I=M2,M1
         IF (SSJ.EQ.IS) M3=I-1
         IF (M3.LT.M4) GO TO 802
         JJ=JA+WTW(I,JS+1)+WTY(M4)
         CII=C(II)
         SII=0.0
         DO 803 J=M4,M3
            LAD=ARR+LADD(J+AOS)
            Z=VAL1*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
            SII=SII+Z*C(JJ)
            S(JJ)=S(JJ)+Z*CII
  803    JJ=JJ+1
         S(II)=S(II)+SII
  802 II=II+1
C
C NO 137
C
 6633 SSJ=IEOR(IS,JS)
      IF (IS.GT.SSJ) RETURN
      M1=MX(SSJ+1)
      M2=MN(SSJ+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      IF (M1.LT.M2) RETURN
      M3=MX(IS+1)
      M4=MN(IS+1)
      IF (M4.GT.N) RETURN
      DO 9082 I=M2,M1
         IF (IS.EQ.SSJ) M3=I-1
         IF (M3.LT.M4) GO TO 9082
         II=IA+WTY(M4)
         JJ=JA+WTW(I,JS+1)+WTY(M4)
         LAD=ARR+LADD(I+AOS)
         Z=VAL1*(INT(LAD+TR1)+VAL2*INT(LAD+TR2))
         DO 9083 J=M4,M3
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
            II=II+1
 9083    JJ=JJ+1
 9082 CONTINUE
      RETURN
C
C   WZ   NO 23 ISEG=10
C
 3091 IF (IS.NE.0.OR.JS.NE.0) GO TO 6634
      JJ=JA
      CJJ=C(JJ)
      SJJ=0.0
      DO 904 I=1,N
         II=IA+WAB(I)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)+LADD(I+LOS)
         Z=VAL1*INT(LAD+1)
         S(II)=S(II)+Z*CJJ
         SJJ=SJJ+Z*C(II)
  904 CONTINUE
      S(JJ)=S(JJ)+SJJ
C
C      NO  25
C
 6634 VAL4=SQRT2*VAL1
      IF (JS.NE.0) RETURN
      JJ=JA
      CJJ=C(JJ)
      SJJ=0.0
      DO 905 I=2,N
         SSI=IEOR(SS(I),IS)
         IF (SSI.GT.SS(I)) GO TO 905
         M1=MX(SSI+1)
         M2=MN(SSI+1)
         IF (M2.GT.N) GO TO 905
         IF (SSI.EQ.SS(I)) M1=I-1
         II=IA+WTW(I,IS+1)+WTY(M2)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+KADD(I+AOS)
         DO 906 J=M2,M1
            LAD1=LAD+LADD(J+LOS)
            Z=VAL4*INT(LAD1+1)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
  906    II=II+1
  905 CONTINUE
      S(JJ)=S(JJ)+SJJ
      RETURN
C
C    YY    NO  9    FOUR EXTERNALS
C***********************************************************************
C
  999 WRITE(ITAPE6,998) M,ISEG
  998 FORMAT (2I4)
      RETURN
C
C     CLOSED INTERNAL LOOP
C***********************************************************************
  318 Z=VAL1*(INT(ARR+TR1)+VAL2*INT(ARR+TR2))
      II=IA
      JJ=JA
      DO 4000 I=1,N1
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
         II=II+1
 4000 JJ=JJ+1
      RETURN
C
C***********************************************************************
  314 PRINT 315
  315 FORMAT(' INCORRECT ENTRY TO THIS SUBROUTINE. THIS IS FOR 4X.')
      CALL MABORT
C
C     ----- INTERNAL CASE FOR TRACKS OF (3,2,1) -----
C***********************************************************************
  319 CONTINUE
      Z=VAL1*(INT(ARR+1)+INT(ARR+2)+INT(ARR+3))
      II=IA
      JJ=JA
      DO 3999 I=1,N1
         S(JJ)=S(JJ)+Z*C(II)
         S(II)=S(II)+Z*C(JJ)
         II=II+1
 3999 JJ=JJ+1
      RETURN
C
C***********************************************************************
C
  317 Z=VAL1*(INT(ARR+TR1)+VAL2*INT(ARR+TR2))
C   ZZ   ARRIVAL,FOR COMPLETENESS
      WRITE(ITAPE6,*) ' ZZ ARRIVAL IN EXTERNAL'
      CALL MABORT
      END
