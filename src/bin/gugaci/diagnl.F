cets110990      SUBROUTINE INITDG(INT,D,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
      SUBROUTINE diagnl(INT,D,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ARR,TR1,TR2,OS,ASM,WTW,WTX,WTY,WAB,SS,SYMORB,AOS
      INTEGER BMAX,ORBFRM,SSJ
      REAL*8 INT(NMAX),D(NWKS)
      real*8 sqrt2,sqrt3,sqt1p5,asqrt2,hsqrt3
      DIMENSION IJADD(NUMIJ),KADD(SYMORB),LADD(SYMORB),WTW(ORBFRM,NSYM)
      DIMENSION WTX(ORBFRM,NSYM),WTY(ORBFRM),WAB(ORBFRM),SS(NORBS)
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /INTS/   NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
C  UNIVERSAL DESIGNATION FOR THESE COMMONS LOCAL VARIANT FOLLOWS
C     COMMON /ALL/ ARR,VAL1,VAL2,VAL3,ITR1,ITR2,IA,JA,ITYPE,ISEGT
C    *,LVFRM1,NLWKI,NLWKJ,IMAX,IMIN
C     COMMON /SYMM/ JSM,JFSYM,IFSYM,MAXSYM(8),MINSYM(8),ISMOFF(8)
C    #,             NUMSYM(8)
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP99,ITAP94
     *,             ITAPE3,ITAP95,ITAP96
      COMMON /ALL/VAL1,VAL2,VAL3,ARR,TR1,TR2,IA,JA,M,ISEG,N,N1,N2
     *,           IMAX,IMIN
      COMMON /SYMM/ ASM,JS,IS,MX(8),MN(8),OS(8),NUMSYM(8)
C:
      COMMON /MINMAX/ IMING,IMAXG,JMING,JMAXG
      common /sqrts/ sqrt2,sqrt3,sqt1p5,asqrt2,hsqrt3
C:
C
cets111090      SQRT2=DSQRT(2.0D+00)
cets111090      SQRT3=DSQRT(3.0D+00)
cets111090      SQT1P5=DSQRT(1.5D+00)
cets110990      DO 5 I=1,NSYM
cets110990         MN(I)=1000
cets110990    5 MX(I)=0
cets110990      MN(1)=1
cets110990      ISM=SS(1)
cets110990      DO 20 I=2,N
cets110990         IF(SS(I).EQ.ISM) GO TO 20
cets110990         MX(ISM)=I-1
cets110990         ISM=SS(I)
cets110990         MN(ISM)=I
cets110990   20 CONTINUE
cets110990      MX(ISM)=N
cets110990      WRITE(ITAPE3,29)
cets110990   29 FORMAT(/,1X,'IN SUBROUTINE INITDG ')
cets110990      WRITE(ITAPE3,30) (MN(J),J=1,NSYM)
cets110990   30 FORMAT(1H ,'MINSYM =',8I9)
cets110990      WRITE(ITAPE3,35) (MX(J),J=1,NSYM)
cets110990   35 FORMAT(1H ,'MAXSYM =',8I9)
cets110990      RETURN
C
cets110990      ENTRY DIAGNL
C
      AOS=OS(ASM+1)
      IA=JA
      GO TO (314,305,304,315,316,317,318,321),M
      WRITE(ITAPE6,10) M
   10 FORMAT(1X,'UNKNOWN TYPE VALUE--',I5)
      CALL MABORT
C    YY   NO   10   ISEG=4
  304 IF (ISEG.EQ.5) GO TO 3043
      IF (ISEG.EQ.6) GO TO 3043
      IF (ISEG.EQ.8) GO TO 3042
      IF (ISEG.NE.4) GO TO 999
      VAL4=-VAL1/SQRT2
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.EQ.1000) RETURN
      IF (M1.GT.N) M1=N
      IF (M1.LT.M2) RETURN
      II=IA+WTY(M2)
      JJ=II
      DO 460 I=M2,M1
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
         Z=VAL4*(INT(LAD+1)-2.*INT(LAD+2))
         D(JJ)=D(JJ)+Z
  460 JJ=JJ+1
      RETURN
C  YY   NO 10 ISEG=8
 3042 VAL4=VAL1*SQT1P5
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.EQ.1000) RETURN
      IF (M1.GT.N) M1=N
      IF (M1.LT.M2) RETURN
      JJ=JA+WTY(M2)
      DO 430 I=M2,M1
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
         Z=VAL4*INT(LAD+1)
         D(JJ)=D(JJ)+Z
  430 JJ=JJ+1
      RETURN
C  YY  NO 10  ISEG=5,6
 3043 CONTINUE
      XX=-VAL1/SQRT2+VAL3*SQT1P5
      IF (ABS(XX).LT.1.0D-06) GO TO 490
      VAL5=SQRT2*VAL1
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.EQ.1000) RETURN
      IF (M1.GT.N) M1=N
      IF (M1.LT.M2) RETURN
      JJ=JA+WTY(M2)
      DO 493 I=M2,M1
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
         Z=XX*INT(LAD+1)+VAL5*INT(LAD+2)
         D(JJ)=D(JJ)+Z
  493 JJ=JJ+1
      RETURN
  490 VAL5=SQRT2*VAL1
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.EQ.1000) RETURN
      IF (M1.GT.N) M1=N
      IF (M1.GT.M2) RETURN
      JJ=JA+WTY(M2)
      DO 494 I=M2,M1
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
         Z=VAL5*INT(LAD+2)
         D(JJ)=D(JJ)+Z
  494 JJ=JJ+1
      RETURN
C   WW    NO 68  ISEG=4,5,6
  314 IF (ISEG.EQ.8) RETURN
      IF (ISEG.EQ.6) GO TO 299
      IF (ISEG.EQ.4.OR.ISEG.EQ.5) GO TO 299
      GO TO 999
  299 IF (JS.NE.0) GO TO 6600
      VAL4=-VAL1*SQRT2
      DO 1401 I=1,N
         JJ=JA+WAB(I)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD=KADD(I+AOS)+LADD(I+LOS)+ARR
         Z=VAL4*(INT(LAD+1)-2.*INT(LAD+2))
         D(JJ)=D(JJ)+Z
 1401 CONTINUE
C   NO 81
 6600 VAL4=-VAL1/SQRT2
      DO 1406 I=2,N
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 1406
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 1406
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 1406
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
         Z1=VAL4*(INT(LAD1+1)-2.0D+00*INT(LAD1+2))
         JJ=JA+WTW(I,1+JS)+WTY(M2)
         DO 1407 J=M2,M1
            LSM=IEOR(ASM,SS(J))
            LOS=OS(LSM+1)
            LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
            Z=VAL4*(INT(LAD+1)-2.*INT(LAD+2))+Z1
            D(JJ)=D(JJ)+Z
 1407    JJ=JJ+1
 1406 CONTINUE
      RETURN
C    XX   NO   11   ISEG=4
 3051 CONTINUE
      VAL4=-VAL1/SQRT2
      DO 501 I=2,N
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 501
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 501
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 501
         JJ=JA+WTX(I,1+JS)+WTY(M2)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
         Z1=VAL4*(INT(LAD1+1)-2.0D+00*INT(LAD1+2))
         DO 502 J=M2,M1
            LSM=IEOR(ASM,SS(J))
            LOS=OS(LSM+1)
            LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
            Z=VAL4*(INT(LAD+1)-2.*INT(LAD+2))+Z1
            D(JJ)=D(JJ)+Z
  502    JJ=JJ+1
  501 CONTINUE
      RETURN
  305 IF (ISEG.EQ.4) GO TO 3051
      IF (ISEG.EQ.5) GO TO 3052
      IF (ISEG.EQ.8) GO TO 3054
      IF (ISEG.NE.6) GO TO 999
C XX   NO  11   ISEG=6,5
 3052 XX=-VAL1/SQRT2+VAL3
      VAL5=SQRT2*VAL1
      IF (ABS(XX).LT.1.0D-06) GO TO 550
      DO 560 I=2,N
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 560
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 560
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 560
         JJ=JA+WTX(I,1+JS)+WTY(M2)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
         Z1=XX*INT(LAD1+1)+VAL5*INT(LAD1+2)
         DO 561 J=M2,M1
            LSM=IEOR(ASM,SS(J))
            LOS=OS(LSM+1)
            LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
            Z=XX*INT(LAD+1)+VAL5*INT(LAD+2)+Z1
            D(JJ)=D(JJ)+Z
  561    JJ=JJ+1
  560 CONTINUE
      RETURN
  550 DO 562 I=2,N
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 562
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 562
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 562
         JJ=JA+WTX(I,1+JS)+WTY(M2)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
         Z1=VAL5*INT(LAD1+2)
         DO 563 J=M2,M1
            LSM=IEOR(ASM,SS(J))
            LOS=OS(LSM+1)
            LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
            Z=VAL5*INT(LAD+2)+Z1
            D(JJ)=D(JJ)+Z
  563    JJ=JJ+1
  562 CONTINUE
      RETURN
C   XX   NO 11  FOUR EXTERNALS
  316 JJA=JA
      M5=IMAX
      M6=IMIN
      IF (M6.LT.2) M6=2
      DO 2002 I=M6,M5
         IF (I.GT.IMAXG.OR.I.LT.IMING.OR.I.GT.JMAXG.OR.I.LT.JMING)GOTO
     *   2002
         LAD=IJADD(I*(I+1)/2)
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2002
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 2002
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 2002
         JJ=JA+WTX(I,1+JS)+WTY(M2)
         KOS=OS(1)
         LOS=OS(SS(I)+1)
         LAD2=LAD+KADD(I+KOS)+LADD(I+LOS)
         Z1=INT(LAD2+2)
         DO 2003 J=M2,M1
            LOS=OS(SS(J)+1)
C:    LAD1=LAD+KADD(J+KOS)+LADD(J+LOS)
            D(JJ)=D(JJ)+   Z1
 2003    JJ=JJ+1
 2002 CONTINUE
      JJA=JA
      M5=IMAX
      M6=IMIN
      IF (M6.LT.2) M6=2
      DO 8002 I=M6,M5
         IF (I.GT.IMAXG.OR.I.LT.IMING.OR.I.GT.JMAXG.OR.I.LT.JMING)GOTO
     *   8002
         LAD=IJADD(I*(I+1)/2)
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 8002
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 8002
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 8002
         JJ=JA+WTX(I,1+JS)+WTY(M2)
         KOS=OS(1)
         LOS=OS(SS(I)+1)
C:    LAD2=LAD+KADD(I+KOS)+LADD(I+LOS)
C:    Z1=INT(LAD2+2)
         DO 8003 J=M2,M1
            LOS=OS(SS(J)+1)
            LAD1=LAD+KADD(J+KOS)+LADD(J+LOS)
            D(JJ)=D(JJ)+(-INT(LAD1+1)+INT(LAD1+2))
 8003    JJ=JJ+1
 8002 CONTINUE
      DO 2082 I=2,N
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I))GO TO 2082
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 2082
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 2082
         M5=IMAX
         M6=IMIN
         IF (M5.GT.M1) M5=M1
         IF (M6.LT.M2) M6=M2
         JJ=JA+WTX(I,1+JS)+WTY(M6)
         KOS=OS(1)
         DO 2083 J=M6,M5
            IF (J.GT.IMAXG.OR.J.LT.IMING.OR.J.GT.JMAXG.OR.J.LT.JMING)
     *      GOTO 2083
            LOS=OS(SS(J)+1)
            LAD=IJADD(J*(J+1)/2)+KADD(J+KOS)+LADD(J+LOS)
            Z=INT(LAD+2)
            D(JJ)=D(JJ)+Z
 2083    JJ=JJ+1
 2082 CONTINUE
      RETURN
C   WW   NO   68  FOUR EXTERNALS
  317 CONTINUE
      M5=IMAX
      M6=IMIN
      IF(IS.NE.0) GO TO 6601
      DO 2030 I=M6,M5
         IF (I.GT.IMAXG.OR.I.LT.IMING.OR.I.GT.JMAXG.OR.I.LT.JMING)GOTO
     *   2030
         JJ=JA+WAB(I)
         KOS=OS(1)
         LOS=OS(SS(I)+1)
         LAD=IJADD(I*(I+1)/2)+KADD(I+KOS)+LADD(I+LOS)
         Z=INT(LAD+2)+INT(LAD+2)+INT(LAD+1)
         D(JJ)=D(JJ)+Z
 2030 CONTINUE
C   NO  81
 6601 M6=IMIN
      IF (M6.LT.2) M6=2
      DO 2046 I=M6,M5
         IF (I.GT.IMAXG.OR.I.LT.IMING.OR.I.GT.JMAXG.OR.I.LT.JMING)GOTO
     *   2046
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2046
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 2046
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 2046
         JJ=JA+WTW(I,1+JS)+WTY(M2)
         LAD=IJADD(I*(I+1)/2)
         KOS=OS(1)
         LOS=OS(SS(I)+1)
         LAD2=LAD+KADD(I+KOS)+LADD(I+LOS)
         Z1=INT(LAD2+2)
         DO 2047 J=M2,M1
C:    LOS=OS(SS(J)+1)
C:    LAD1=LAD+KADD(J+KOS)+LADD(J+LOS)
C:    Z=INT(LAD1+1)+INT(LAD1+2)+Z1
            D(JJ)=D(JJ)+Z1
 2047    JJ=JJ+1
 2046 CONTINUE
      M6=IMIN
      IF (M6.LT.2) M6=2
      DO 8046 I=M6,M5
         IF (I.GT.IMAXG.OR.I.LT.IMING.OR.I.GT.JMAXG.OR.I.LT.JMING)GOTO
     *   8046
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 8046
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 8046
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 8046
         JJ=JA+WTW(I,1+JS)+WTY(M2)
         LAD=IJADD(I*(I+1)/2)
         KOS=OS(1)
         LOS=OS(SS(I)+1)
C:    LAD2=LAD+KADD(I+KOS)+LADD(I+LOS)
C:    Z1=INT(LAD2+2)
         DO 8047 J=M2,M1
            LOS=OS(SS(J)+1)
            LAD1=LAD+KADD(J+KOS)+LADD(J+LOS)
            Z=INT(LAD1+1)+INT(LAD1+2)
            D(JJ)=D(JJ)+Z
 8047    JJ=JJ+1
 8046 CONTINUE
      DO 2084 I=2,N
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2084
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 2084
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 2084
         M5=IMAX
         M6=IMIN
         IF (M5.GT.M1) M5=M1
         IF (M6.LT.M2) M6=M2
         JJ=JA+WTW(I,1+JS)+WTY(M6)
         KOS=OS(1)
         DO 2085 J=M6,M5
            IF (J.GT.IMAXG.OR.J.LT.IMING.OR.J.GT.JMAXG.OR.J.LT.JMING)
     *      GOTO 2085
            LOS=OS(SS(J)+1)
            LAD=IJADD(J*(J+1)/2)+KADD(J+KOS)+LADD(J+LOS)
            Z=INT(LAD+2)
            D(JJ)=D(JJ)+Z
 2085    JJ=JJ+1
 2084 CONTINUE
      RETURN
C   YY   NO 10 FOUR EXTERNALS
  315 CONTINUE
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.EQ.1000) RETURN
      IF (M1.GT.N) M1=N
      IF (M1.LT.M2) RETURN
      M5=IMAX
      M6=IMIN
      IF (M5.GT.M1) M5=M1
      IF (M6.LT.M2) M6=M2
      JJ=JA+WTY(M6)
      DO 2080 I=M6,M5
         IF (I.GT.IMAXG.OR.I.LT.IMING.OR.I.GT.JMAXG.OR.I.LT.JMING)GOTO
     *   2080
         KOS=OS(1)
         LOS=OS(SS(I)+1)
         LAD=IJADD(I*(I+1)/2)+KADD(I+KOS)+LADD(I+LOS)+2
         D(JJ)=D(JJ)+INT(LAD)
 2080 JJ=JJ+1
      RETURN
C   XX    NO  11  ISEG=8
 3054 CONTINUE
      DO 510 I=2,N
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 510
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.EQ.1000) GO TO 510
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 510
         JJ=JA+WTX(I,1+JS)+WTY(M2)
         LSM=IEOR(ASM,SS(I))
         LOS=OS(LSM+1)
         LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)+1
         Z1=VAL1*INT(LAD1)
         DO 511 J=M2,M1
            LSM=IEOR(ASM,SS(J))
            LOS=OS(LSM+1)
            LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
            Z=VAL1*INT(LAD+1)+Z1
            D(JJ)=D(JJ)+Z
  511    JJ=JJ+1
  510 CONTINUE
      RETURN
C   CLOSED INTERNAL LOOP
  318 Z=VAL1*(INT(ARR+TR1)+VAL2*INT(ARR+TR2))
      JJ=JA
      DO 4000 I=1,N1
         D(JJ)=D(JJ)+Z
 4000 JJ=JJ+1
      RETURN
C ZZ   ARRIVAL IS AN ERROR
  321 CONTINUE
  999 WRITE(ITAPE6,998) M,ISEG
  998 FORMAT (2I4)
      RETURN
      END
