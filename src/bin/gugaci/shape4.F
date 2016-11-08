      SUBROUTINE shape4(INT,C,S,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
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
cets110990      ENTRY SHAPE4
C
      AOS=OS(ASM+1)
      GO TO (301,301,301,301,301,301,301,301,301,301,301,301,301,314
     1,315,316,317,318,319),M
C
      WRITE(ITAPE6,*) ' FUNNY M'
CALL MABORT
C
  301 PRINT 302
  302 FORMAT(' THIS ENTRY IS NOT POSSIBLE')
      CALL MABORT
C
C    YY    NO  9    FOUR EXTERNALS
C***********************************************************************
  314 IF (IS.GT.JS) RETURN
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      M3=MX(IS+1)
      M4=MN(IS+1)
      IF (M4.GT.N) RETURN
      M5=IMAX
      M6=IMIN
      IF (M5.GT.M1) M5=M1
      IF (M6.LT.M2) M6=M2
      IF (M5.LT.M6) RETURN
      JJ=JA+WTY(M6)
      DO 2000 I=M6,M5
C
         IF (I.GT.JMAXG.OR.I.LT.JMING) GO TO 2000
C
         IF (IS.EQ.JS) M3=I-1
         IF (M3.LT.M4) GO TO 2000
         KOS=OS(1)
         LOS=OS(SS(I)+1)
         LAD=IJADD(I*(I+1)/2)+KADD(I+KOS)+3
         II=IA+WTY(M4)
         CJJ=C(JJ)
         SJJ=0.0
         DO 2001 J=M4,M3
            LAD1=LAD+LADD(J+LOS)
            Z=INT(LAD1)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
 2001    II=II+1
         S(JJ)=S(JJ)+SJJ
 2000 JJ=JJ+1
      RETURN
C
C   XX   NO 13 & 161 FOUR EXTERNALS
C***********************************************************************
  315 M5=IMAX
      M6=IMIN
      IF (M6.LT.3) M6=3
      IF (M5.LT.M6) GO TO 2003
      DO 2004 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         SSI=IEOR(SS(I),IS)
         IF (SSJ.GT.SS(I)) GO TO 2004
         IF (SSI.GT.SSJ) GO TO 2004
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.GT.N) GO TO 2004
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M2.LT.2) M2=2
         IF (M1.LT.M2) GO TO 2004
         M3=MX(SSI+1)
         M4=MN(SSI+1)
         IF (M4.GT.N) GO TO 2004
         JJ=JA+WTX(I,JS+1)+WTY(M2)
         IIA=IA+WTX(I,IS+1)
         KOS=OS(1)
         LAD=IJADD(I*(I+1)/2)
         DO 2005 J=M2,M1
            IF (SSI.EQ.SSJ) M3=J-1
            IF (M3.LT.M4) GO TO 2005
            LOS=OS(SS(J)+1)
            LAD1=LAD+KADD(J+KOS)
            II=IIA+WTY(M4)
            CJJ=C(JJ)
            SJJ=0.0
            DO 2006 K=M4,M3
               LAD2=LAD1+LADD(K+LOS)
               Z=-INT(LAD2+1)+INT(LAD2+2)
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
 2006       II=II+1
            S(JJ)=S(JJ)+SJJ
 2005    JJ=JJ+1
 2004 CONTINUE
 2003 CONTINUE
C
C     NO 161
      IF (N.LT.3) GO TO 8084
      DO 7084 I=3,N
         SSJ=IEOR(SS(I),JS)
         SSI=IEOR(SS(I),IS)
         IF (SSJ.GT.SS(I)) GO TO 7084
         IF (SSI.GT.SSJ) GO TO 7084
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.GT.N) GO TO 7084
         IF (M2.LT.2) M2=2
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 7084
         M5=IMAX
         M6=IMIN
         IF (M5.GT.M1) M5=M1
         IF (M6.LT.M2) M6=M2
         M3=MX(SSI+1)
         M4=MN(SSI+1)
         IF (M4.GT.N) GO TO 7084
         JJ=JA+WTX(I,JS+1)+WTY(M6)
         IIA=IA+WTX(I,IS+1)
         KOS=OS(1)
         DO 7085 J=M6,M5
            IF (SSI.EQ.SSJ) M3=J-1
            IF (M3.LT.M4) GO TO 7085
            LOS=OS(SS(J)+1)
            LAD3=IJADD(J*(J+1)/2)+KADD(J+KOS)+3
            II=IIA+WTY(M4)
            CJJ=C(JJ)
            SJJ=0.0
            DO 7086 K=M4,M3
               Z=INT(LAD3+LADD(K+LOS))
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
 7086       II=II+1
            S(JJ)=S(JJ)+SJJ
 7085    JJ=JJ+1
 7084 CONTINUE
 8084 CONTINUE
C   NO  14
C
      M5=IMAX
      M6=IMIN
      IF (M6.LT.3) M6=3
      IF (M5.LT.M6) GO TO 8007
      DO 2007 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         SSI=IEOR(IS,SSJ)
         IF (SSJ.GT.SSI) GO TO 2007
         IF (SSI.GT.SS(I)) GO TO 2007
         M1=MX(SSI+1)
         M2=MN(SSI+1)
         IF (M2.GT.N) GO TO 2007
         IF (M2.LT.2) M2=2
         IF (SSI.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 2007
         M3=MX(SSJ+1)
         M4=MN(SSJ+1)
         IF (M4.GT.N) GO TO 2007
         JJA=JA+WTX(I,JS+1)
         KOS=OS(1)
         LOS1=OS(SS(I)+1)
         LAD3=IJADD(I*(I+1)/2)+KADD(I+KOS)+3
         DO 2008 J=M2,M1
            IF (SSI.EQ.SSJ) M3=J-1
            IF (M3.LT.M4) GO TO 2008
            Z1=INT(LAD3+LADD(J+LOS1))
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD1=IJADD(I*(I-1)/2+J)
            II=IA+WTX(J,IS+1)+WTY(M4)
            JJ=JJA+WTY(M4)
            DO 2009 K=M4,M3
               LSM=IEOR(KSM,SS(K))
               LOS=OS(LSM+1)
               LAD2=LAD1+KADD(K+KOS)+LADD(K+LOS)
               Z=-INT(LAD2+1)+INT(LAD2+2)+Z1
               S(JJ)=S(JJ)+Z*C(II)
               S(II)=S(II)+Z*C(JJ)
               II=II+1
 2009       JJ=JJ+1
 2008    CONTINUE
 2007 CONTINUE
 8007 CONTINUE
C
C    NO 16   BOTH
C
      IF (M5.LT.M6) GO TO 8010
      DO 2010 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         SSI=IEOR(IS,SSJ)
         IF (SSI.GT.SSJ) GO TO 2010
         IF (SSJ.GT.SS(I)) GO TO 2010
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.GT.N) GO TO 2010
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M2.LT.2) M2=2
         M3=MX(SSI+1)
         M4=MN(SSI+1)
         IF (M4.GT.N) GO TO 2010
         JJ=JA+WTX(I,JS+1)+WTY(M2)
         KOS=OS(1)
         LOS=OS(SS(I)+1)
         LAD2=IJADD(I*(I+1)/2)+KADD(I+KOS)
         DO 2011 J=M2,M1
            IF (SSI.EQ.SSJ) M3=J-1
            IF (M3.LT.M4) GO TO 2011
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD3=IJADD(I*(I-1)/2+J)+KADD(J+KOS)
            II=IA+WTX(J,IS+1)+WTY(M4)
            CJJ=C(JJ)
            SJJ=0.0
            DO 2012 K=M4,M3
               LAD4=LAD2+LADD(K+LOS)
               LAD5=LAD3+LADD(K+LOS)
               Z=-INT(LAD4+3)-INT(LAD5+2)+INT(LAD5+1)
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
 2012       II=II+1
            S(JJ)=S(JJ)+SJJ
 2011    JJ=JJ+1
 2010 CONTINUE
 8010 CONTINUE
C
C    NO  20
C
      IF (M5.LT.4) GO TO 8013
      M6=IMIN
      IF (M6.LT.4) M6=4
      DO 2013 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2013
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.GT.N) GO TO 2013
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M2.LT.3) M2=3
         IF (M1.LT.M2) GO TO 2013
         JJ=JA+WTX(I,JS+1)+WTY(M2)
         DO 2014 J=M2,M1
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD1=IJADD(I*(I-1)/2+J)
            CJJ=C(JJ)
            SJJ=0.0
            DO 2015 K=2,J-1
               SSI=IEOR(SS(K),IS)
               IF (SSI.GT.SS(K)) GO TO 2015
               M3=MX(SSI+1)
               M4=MN(SSI+1)
               IF (M4.GT.N) GO TO 2015
               IF (SS(K).EQ.SSI) M3=K-1
               IF (M3.LT.M4) GO TO 2015
               II=IA+WTX(K,IS+1)+WTY(M4)
               LSM=IEOR(KSM,SS(K))
               LOS=OS(LSM+1)
               LAD2=LAD1+KADD(K+KOS)
               DO 2016 L=M4,M3
                  LAD3=LAD2+LADD(L+LOS)
                  Z=INT(LAD3+1)-INT(LAD3+3)
                  S(II)=S(II)+Z*CJJ
                  SJJ=SJJ+Z*C(II)
 2016          II=II+1
 2015       CONTINUE
            S(JJ)=S(JJ)+SJJ
 2014    JJ=JJ+1
 2013 CONTINUE
 8013 CONTINUE
C
C      NO  21
C
      IF(M5.LT.M6) GO TO 8017
      DO 2017 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2017
         JJA=JA+WTX(I,JS+1)
         DO 2018 J=3,I-1
            SSI=IEOR(SS(J),IS)
            IF (SSI.GT.SS(J)) GO TO 2018
            IF (SSJ.GT.SSI) GO TO 2018
            M1=MX(SSI+1)
            M2=MN(SSI+1)
            IF (M2.GT.N) GO TO 2018
            IF (SSI.EQ.SS(J)) M1=J-1
            IF (M2.LT.2) M2=2
            IF (M1.LT.M2) GOTO 2018
            II=IA+WTX(J,IS+1)+WTY(M2)
            M3=MX(SSJ+1)
            M4=MN(SSJ+1)
            IF (M4.GT.N) GO TO 2018
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD1=IJADD(I*(I-1)/2+J)
            DO 2019 K=M2,M1
               IF (SSJ.EQ.SSI) M3=K-1
               IF (M3.LT.M4) GO TO 2019
               JJ=JJA+WTY(M4)
               LAD2=LAD1+KADD(K+KOS)
               LSM=IEOR(KSM,SS(K))
               LOS=OS(LSM+1)
               CII=C(II)
               SII=0.0
               DO 2020 L=M4,M3
                  LAD3=LAD2+LADD(L+LOS)
                  Z=-INT(LAD3+1)+INT(LAD3+2)
                  SII=SII+Z*C(JJ)
                  S(JJ)=S(JJ)+Z*CII
 2020          JJ=JJ+1
               S(II)=S(II)+SII
 2019       II=II+1
 2018    CONTINUE
 2017 CONTINUE
 8017 CONTINUE
C
C   NO  22
C
      IF(M5.LT.M6) RETURN
      DO 2021 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2021
         JJA=JA+WTX(I,JS+1)
         DO 2022 J=3,I-1
            SSI=IEOR(SS(J),IS)
            IF (SSJ.GT.SS(J))GO TO 2022
            IF (SSI.GT.SSJ) GO TO 2022
            M1=MX(SSJ+1)
            M2=MN(SSJ+1)
            IF (M2.GT.N) GO TO 2022
            IF (SSJ.EQ.SS(J)) M1=J-1
            IF (M2.LT.2) M2=2
            IF (M1.LT.M2) GO TO 2022
            JJ=JJA+WTY(M2)
            IIA=IA+WTX(J,IS+1)
            M3=MX(SSI+1)
            M4=MN(SSI+1)
            IF (M4.GT.N) GO TO 2022
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD1=IJADD(I*(I-1)/2+J)
            DO 2023 K=M2,M1
               IF (SSI.EQ.SSJ) M3=K-1
               IF (M3.LT.M4) GO TO 2023
               II=IIA+WTY(M4)
               LSM=IEOR(KSM,SS(K))
               LOS=OS(LSM+1)
               LAD2=LAD1+KADD(K+KOS)
               CJJ=C(JJ)
               SJJ=0.0
               DO 2024 L=M4,M3
                  LAD3=LAD2+LADD(L+LOS)
                  Z=-INT(LAD3+3)+INT(LAD3+2)
                  S(II)=S(II)+Z*CJJ
                  SJJ=SJJ+Z*C(II)
 2024          II=II+1
               S(JJ)=S(JJ)+SJJ
 2023       JJ=JJ+1
 2022    CONTINUE
 2021 CONTINUE
      RETURN
C
C   WW   NO 69A FOUR EXTERNALS
C***********************************************************************
  316 M5=IMAX
      M6=IMIN
      IF (M6.LT.2) M6=2
      IF (M5.LT.M6) GO TO 6586
      IF (IS.NE.0.OR.JS.NE.0) GO TO 6586
      DO 2031 I=M6,M5
         JJ=JA+WAB(I)
         LAD=IJADD(I*(I+1)/2)
         CJJ=C(JJ)
         KOS=OS(1)
         SJJ=0.0
         DO 2032 J=1,I-1
            LOS=OS(SS(J)+1)
            LAD1=KADD(J+KOS)+LADD(J+LOS)+LAD
            II=IA+WAB(J)
            Z=INT(LAD1+1)
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
 2032    CONTINUE
         S(JJ)=S(JJ)+SJJ
 2031 CONTINUE
C
C    NO   72A
C
 6586 IF (JS.NE.0) GO TO  6635
      IF (M5.LT.M6) GO TO 6635
      DO 2033 I=M6,M5
         SSI=IEOR(SS(I),IS)
         IF (SSI.GT.SS(I)) GO TO 2033
         M1=MX(SSI+1)
         M2=MN(SSI+1)
         IF (M2.GT.N) GO TO 2033
         IF (SSI.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GOTO 2033
         JJ=JA+WAB(I)
         II=IA+WTW(I,IS+1)+WTY(M2)
         KOS=OS(1)
         LOS=OS(SS(I)+1)
         LAD=IJADD(I*(I+1)/2)+KADD(I+KOS)
         CJJ=C(JJ)
         SJJ=0.0
         DO 2034 J=M2,M1
            LAD1=LAD+LADD(J+LOS)
            Z=SQRT2*(INT(LAD1+3)+INT(LAD1+1))
            S(II)=S(II)+Z*CJJ
            SJJ=SJJ+Z*C(II)
 2034    II=II+1
         S(JJ)=S(JJ)+SJJ
 2033 CONTINUE
C
C    NO  73
C
 6635 M6=IMIN
      IF (M6.LT.3) M6=3
      IF (M5.LT.M6) GO TO 6636
      IF (IS.NE.0) GO TO 6636
      DO 2035 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2035
         JJA=JA+WTW(I,JS+1)
         M3=MX(SSJ+1)
         M4=MN(SSJ+1)
         IF (M4.GT.N) GO TO 2035
         DO 2036 J=2,I-1
            IF (SSJ.GT.SS(J)) GO TO 2036
            IF (SSJ.EQ.SS(J)) M3=J-1
            IF (M3.LT.M4) GO TO 2036
            JJ=JJA+WTY(M4)
            II=IA+WAB(J)
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LOS=OS(SS(I)+1)
            LAD1=IJADD(I*(I-1)/2+J)+KADD(J+KOS)
            CII=C(II)
            SII=0.0
            DO 2037 K=M4,M3
               LAD2=LAD1+LADD(K+LOS)
               Z=SQRT2*INT(LAD2+1)
               SII=SII+Z*C(JJ)
               S(JJ)=S(JJ)+Z*CII
 2037       JJ=JJ+1
            S(II)=S(II)+SII
 2036    CONTINUE
 2035 CONTINUE
C
C   NO  74A
C
 6636 IF (JS.NE.0) GO TO 6637
      IF (M5.LT.M6) GO TO 6637
      DO 2038 I=M6,M5
         JJ=JA+WAB(I)
         LAD=IJADD(I*(I+1)/2)
         CJJ=C(JJ)
         KOS=OS(1)
         SJJ=0.0
         DO 2039 J=2,I-1
            SSI=IEOR(SS(J),IS)
            IF (SSI.GT.SS(J)) GO TO 2039
            M1=MX(SSI+1)
            M2=MN(SSI+1)
            IF (M2.GT.N) GO TO 2039
            IF (SSI.EQ.SS(J)) M1=J-1
            IF (M1.LT.M2) GO TO 2039
            II=IA+WTW(J,IS+1)+WTY(M2)
            LOS=OS(SS(J)+1)
            LAD1=LAD+KADD(J+KOS)
            DO 2040 K=M2,M1
               LAD2=LAD1+LADD(K+LOS)
               Z=SQRT2*INT(LAD2+1)
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
 2040       II=II+1
 2039    CONTINUE
         S(JJ)=S(JJ)+SJJ
 2038 CONTINUE
C
C   NO  76
C
 6637 IF (IS.NE.0) GO TO 6638
      IF (M5.LT.M6) GO TO 8041
      DO 2041 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2041
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.GT.N) GO TO 2041
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M2.LT.2) M2=2
         IF (M1.LT.M2) GO TO 2041
         JJ=JA+WTW(I,JS+1)+WTY(M2)
         DO 2042 J=M2,M1
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD1=IJADD(I*(I-1)/2+J)
            CJJ=C(JJ)
            SJJ=0.0
            DO 2043 K=1,J-1
               II=IA+WAB(K)
               LSM=IEOR(KSM,SS(K))
               LOS=OS(LSM+1)
               LAD2=LAD1+KADD(K+KOS)+LADD(K+LOS)
               Z=SQRT2*INT(LAD2+1)
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
 2043       CONTINUE
            S(JJ)=S(JJ)+SJJ
            JJ=JJ+1
 2042    CONTINUE
 2041 CONTINUE
 8041 CONTINUE
C
C    NO  78
C
      M6=IMIN
      IF (M6.LT.2) M6=2
      IF (M5.LT.M6) GO TO 6638
      DO 2044 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2044
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.GT.N) GO TO 2044
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M1.LT.M2) GO TO 2044
         JJ=JA+WTW(I,JS+1)+WTY(M2)
         KOS=OS(1)
         LOS=OS(SS(I)+1)
         LAD=IJADD(I*(I+1)/2)+KADD(I+KOS)
         DO 2045 J=M2,M1
            LAD1=LAD+LADD(J+LOS)
            II=IA+WAB(J)
            Z=SQRT2*(INT(LAD1+3)+INT(LAD1+2))
            S(JJ)=S(JJ)+Z*C(II)
            S(II)=S(II)+Z*C(JJ)
 2045    JJ=JJ+1
 2044 CONTINUE
C
C    NO   84A
C
 6638 M6=IMIN
      IF (M6.LT.3) M6=3
      IF (M5.LT.M6) GOTO 8054
      DO 2048 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         SSI=IEOR(IS,SSJ)
         IF (SSJ.GT.SSI) GO TO 2048
         IF (SSI.GT.SS(I)) GO TO 2048
         M1=MX(SSI+1)
         M2=MN(SSI+1)
         IF (M2.GT.N) GO TO 2048
         M3=MX(SSJ+1)
         M4=MN(SSJ+1)
         IF (M4.GT.N) GO TO 2048
         IF (SSI.EQ.SS(I)) M1=I-1
         IF (M2.LT.2) M2=2
         IF (M1.LT.M2) GO TO 2048
         JJA=JA+WTW(I,JS+1)
         KOS1=OS(1)
         LOS1=OS(SS(I)+1)
         LAD3=IJADD(I*(I+1)/2)+KADD(I+KOS1)+3
         DO 2049 J=M2,M1
            IF (SSI.EQ.SSJ) M3=J-1
            IF (M3.LT.M4) GO TO 2049
            Z1=INT(LAD3+LADD(J+LOS1))
            II=IA+WTW(J,IS+1)+WTY(M4)
            JJ=JJA+WTY(M4)
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD1=IJADD(I*(I-1)/2+J)
            DO 2050 K=M4,M3
               LSM=IEOR(KSM,SS(K))
               LOS=OS(LSM+1)
               LAD2=LAD1+KADD(K+KOS)+LADD(K+LOS)
               Z=INT(LAD2+1)+INT(LAD2+2)+Z1
               S(JJ)=S(JJ)+Z*C(II)
               S(II)=S(II)+Z*C(JJ)
               II=II+1
 2050       JJ=JJ+1
 2049    CONTINUE
 2048 CONTINUE
C
C   NO  88A
C
      DO 2051 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         SSI=IEOR(IS,SSJ)
         IF (SSI.GT.SSJ) GO TO 2051
         IF (SSJ.GT.SS(I)) GO TO 2051
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.GT.N) GO TO 2051
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M2.LT.2) M2=2
         IF (M1.LT.M2) GO TO 2051
         M3=MX(SSI+1)
         M4=MN(SSI+1)
         IF (M4.GT.N) GO TO 2051
         JJ=JA+WTW(I,JS+1)+WTY(M2)
         KOS1=OS(1)
         LOS1=OS(SS(I)+1)
         LAD1=IJADD(I*(I+1)/2)+KADD(I+KOS1)
         DO 2052 J=M2,M1
            IF (SSI.EQ.SSJ) M3=J-1
            IF (M3.LT.M4) GO TO 2052
            II=IA+WTW(J,IS+1)+WTY(M4)
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD2=IJADD(I*(I-1)/2+J)+KADD(J+KOS)
            CJJ=C(JJ)
            SJJ=0.0
            DO 2053 K=M4,M3
               Z=INT(LAD1+LADD(K+LOS1)+3)+INT(LAD2+LADD(K+LOS1)+2)
               Z=Z+INT(LAD2+LADD(K+LOS1)+1)
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
 2053       II=II+1
            S(JJ)=S(JJ)+SJJ
 2052    JJ=JJ+1
 2051 CONTINUE
C
C   NO   89 & 160
C
      DO 2054 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         SSI=IEOR(SS(I),IS)
         IF (SSJ.GT.SS(I)) GO TO 2054
         IF (SSI.GT.SSJ) GO TO 2054
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.GT.N) GO TO 2054
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M2.LT.2) M2=2
         IF (M1.LT.M2) GO TO 2054
         M3=MX(SSI+1)
         M4=MN(SSI+1)
         IF (M4.GT.N) GO TO 2054
         JJ=JA+WTW(I,JS+1)+WTY(M2)
         IIA=IA+WTW(I,IS+1)
         KOS=OS(1)
         LAD=IJADD(I*(I+1)/2)
         DO 2055 J=M2,M1
            IF (SSI.EQ.SSJ) M3=J-1
            IF (M3.LT.M4) GO TO 2055
            II=IIA+WTY(M4)
            LAD1=LAD+KADD(J+KOS)
            LOS=OS(SS(J)+1)
            CJJ=C(JJ)
            SJJ=0.0
            DO 2056 K=M4,M3
               LAD2=LAD1+LADD(K+LOS)
               Z=INT(LAD2+1)+INT(LAD2+2)
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
 2056       II=II+1
            S(JJ)=S(JJ)+SJJ
 2055    JJ=JJ+1
 2054 CONTINUE
 8054 CONTINUE
C
C  NO 160
      IF (N.LT.3) GO TO 7080
      DO 7087 I=3,N
         SSJ=IEOR(SS(I),JS)
         SSI=IEOR(SS(I),IS)
         IF (SSJ.GT.SS(I)) GO TO 7087
         IF (SSI.GT.SSJ) GO TO 7087
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.GT.N) GO TO 7087
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M2.LT.2) M2=2
         M5=IMAX
         M6=IMIN
         IF (M5.GT.M1) M5=M1
         IF (M6.LT.M2) M6=M2
         M3=MX(SSI+1)
         M4=MN(SSI+1)
         IF (M4.GT.N) GO TO 7087
         JJ=JA+WTW(I,JS+1)+WTY(M6)
         IIA=IA+WTW(I,IS+1)
         KOS=OS(1)
         DO 7088 J=M6,M5
            LOS=OS(SS(J)+1)
            LAD3=IJADD(J*(J+1)/2)+KADD(J+KOS)+3
            IF (SSI.EQ.SSJ) M3=J-1
            IF (M3.LT.M4) GO TO 7088
            II=IIA+WTY(M4)
            CJJ=C(JJ)
            SJJ=0.0
            DO 7089 K=M4,M3
               Z=INT(LAD3+LADD(K+LOS))
               S(II)=S(II)+Z*CJJ
               SJJ=SJJ+Z*C(II)
 7089       II=II+1
            S(JJ)=S(JJ)+SJJ
 7088    JJ=JJ+1
 7087 CONTINUE
 7080 CONTINUE
C     NO    99A
C
      M5=IMAX
      M6=IMIN
      IF (M6.LT.4) M6=4
      IF (M5.LT.M6) RETURN
      DO 2057 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2057
         JJA=JA+WTW(I,JS+1)
         DO 2058 J=3,I-1
            SSI=IEOR(SS(J),IS)
            IF(SSJ.GT.SS(J)) GO TO 2058
            IF (SSI.GT.SSJ) GO TO 2058
            M1=MX(SSJ+1)
            M2=MN(SSJ+1)
            IF (M2.GT.N) GO TO 2058
            IF (SSJ.EQ.SS(J)) M1=J-1
            IF (M2.LT.2) M2=2
            JJ=JJA+WTY(M2)
            IIA=IA+WTW(J,IS+1)
            M3=MX(SSI+1)
            M4=MN(SSI+1)
            IF (M4.GT.N) GO TO 2058
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD1=IJADD(I*(I-1)/2+J)
            DO 2059 K=M2,M1
               IF (SSI.EQ.SSJ) M3=K-1
               IF (M3.LT.M4) GO TO 2059
               II=IIA+WTY(M4)
               LSM=IEOR(KSM,SS(K))
               LOS=OS(LSM+1)
               LAD2=LAD1+KADD(K+KOS)
               CJJ=C(JJ)
               SJJ=0.0
               DO 2060 L=M4,M3
                  LAD3=LAD2+LADD(L+LOS)
                  Z=INT(LAD3+3)+INT(LAD3+2)
                  S(II)=S(II)+Z*CJJ
                  SJJ=SJJ+Z*C(II)
 2060          II=II+1
               S(JJ)=S(JJ)+SJJ
 2059       JJ=JJ+1
 2058    CONTINUE
 2057 CONTINUE
C
C     NO   100A
C
      DO 2061 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2061
         JJA=JA+WTW(I,JS+1)
         DO 2062 J=3,I-1
            SSI=IEOR(SS(J),IS)
            IF (SSI.GT.SS(J)) GO TO 2062
            IF (SSJ.GT.SSI) GO TO 2062
            M1=MX(SSI+1)
            M2=MN(SSI+1)
            IF (M2.GT.N) GO TO 2062
            IF (SSI.EQ.SS(J)) M1=J-1
            IF (M2.LT.2) M2=2
            M3=MX(SSJ+1)
            M4=MN(SSJ+1)
            IF (M4.GT.N) GO TO 2062
            II=IA+WTW(J,IS+1)+WTY(M2)
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD1=IJADD(I*(I-1)/2+J)
            DO 2063 K=M2,M1
               IF (SSJ.EQ.SSI) M3=K-1
               IF (M3.LT.M4) GOTO 2063
               JJ=JJA+WTY(M4)
               LSM=IEOR(KSM,SS(K))
               LOS=OS(LSM+1)
               LAD2=LAD1+KADD(K+KOS)
               CII=C(II)
               SII=0.0
               DO 2064 L=M4,M3
                  LAD3=LAD2+LADD(L+LOS)
                  Z=INT(LAD3+1)+INT(LAD3+2)
                  SII=SII+Z*C(JJ)
                  S(JJ)=S(JJ)+Z*CII
 2064          JJ=JJ+1
               S(II)=S(II)+SII
 2063       II=II+1
 2062    CONTINUE
 2061 CONTINUE
C
C   NO   101A
C
      DO 2065 I=M6,M5
         SSJ=IEOR(SS(I),JS)
         IF (SSJ.GT.SS(I)) GO TO 2065
         M1=MX(SSJ+1)
         M2=MN(SSJ+1)
         IF (M2.GT.N) GO TO 2065
         IF (SSJ.EQ.SS(I)) M1=I-1
         IF (M2.LT.3) M2=3
         IF (M1.LT.M2) GOTO 2065
         JJ=JA+WTW(I,JS+1)+WTY(M2)
         DO 2066 J=M2,M1
            KSM=IEOR(SS(I),SS(J))
            KOS=OS(KSM+1)
            LAD1=IJADD(I*(I-1)/2+J)
            CJJ=C(JJ)
            SJJ=0.0
            DO 2067 K=2,J-1
               SSI=IEOR(SS(K),IS)
               IF (SSI.GT.SS(K)) GO TO 2067
               M3=MX(SSI+1)
               M4=MN(SSI+1)
               IF (M4.GT.N) GO TO 2067
               IF (SS(K).EQ.SSI) M3=K-1
               IF (M3.LT.M4) GOTO 2067
               II=IA+WTW(K,IS+1)+WTY(M4)
               LSM=IEOR(KSM,SS(K))
               LOS=OS(LSM+1)
               LAD2=LAD1+KADD(K+KOS)
               DO 2068 L=M4,M3
                  LAD3=LAD2+LADD(L+LOS)
                  Z=INT(LAD3+1)+INT(LAD3+3)
                  S(II)=S(II)+Z*CJJ
                  SJJ=SJJ+Z*C(II)
 2068          II=II+1
 2067       CONTINUE
            S(JJ)=S(JJ)+SJJ
            JJ=JJ+1
 2066    CONTINUE
 2065 CONTINUE
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
