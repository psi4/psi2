      SUBROUTINE INIT4X(INT,C,S,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
      IMPLICIT REAL*8 (A-H,O-Z)
C     EXTENDED DUMMY INT,C,S
      REAL*8 Z,VAL1,VAL2,VAL3
      integer nmax,nwks,numij,nsym,norbs
      INTEGER ARR,TR1,TR2,ASM,AOS,OS,WTW,WTX,WTY,WAB,SS,SSI,SSJ,SYMORB
      INTEGER BMAX,ORBFRM
      REAL*8 INT(NMAX),C(NWKS),S(NWKS)
      DIMENSION KADD(SYMORB),LADD(SYMORB),IJADD(NUMIJ),WTW(ORBFRM,NSYM)
      DIMENSION WTX(ORBFRM,NSYM),WTY(ORBFRM),WAB(ORBFRM),SS(NORBS)
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,NROW4O,NWKSOC,NLEVOC,NORBOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /INTS/ NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
C  UNIVERSAL IDENTITY OF THE OBJECTS IN THESE COMMON
C     COMMON /ALL/ ARR,VAL1,VAL2,VAL3,ITR1,ITR2,IA,JA,ITYPE,ISEGT
C    *,LVFRM1,NLWKI,NLWKJ,IMAX,IMIN
C     COMMON /SYMM/ JSM,JFSYM,IFSYM,MAXSYM(8),MINSYM(8),ISMOFF(8)
C    #,             NUMSYM(8)
      COMMON /SYMM/ ASM,JS,IS,MX(8),MN(8),OS(8),NUMSYM(8)
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP03,ITAP04
     *,             ITAPE3,ITAP05,ITAP06,ITAP56
      COMMON /ALL/VAL1,VAL2,VAL3,ARR,TR1,TR2,IA,JA,M,ISEG,N,N1,N2
     *,           IMAX,IMIN
      COMMON /MINMAX/ IMING,IMAXG,JMING,JMAXG
      SQRT2=DSQRT(2.0D+00)
      SQRT3=DSQRT(3.0D+00)
      SQT1P5=DSQRT(1.5D+00)
CMC   WRITE(6,*) 'CI VECTOR'
CMC   DO IMIKE=1,55
CMC        WRITE(6,9010) IMIKE,C(IMIKE)
CMC   END DO
CMC10 FORMAT(I5,E20.12)
      RETURN
C
      ENTRY SHAPE4(INT,C,S,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
CMC   WRITE(6,*) ' ENTRY SHAPE4'
C
      AOS=OS(ASM+1)
      GO TO (301,301,301,301,301,301,301,301,301,301,301,301,301,314
C     GO TO ( YZ, XZ, XY, YY, XX, ZY, YX, YW, WZ, WY, WX, XW, WW,---YY
C
     1,315,316,301,301,301),M
C       XX, WW, ZZ,---4 INTERNALS),M
C
      STOP
C
  301 PRINT 302
  302 FORMAT(' THIS ENTRY IS NOT POSSIBLE')
      STOP
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
           IF(I.GT.IMAXG.OR.I.LT.IMING.OR.I.GT.JMAXG.OR.I.LT.JMING)GO TO
     *      2000
           IF (IS.EQ.JS) M3=I-1
           IF (M3.LT.M4) GO TO 2000
           KOS=OS(1)
           LOS=OS(SS(I)+1)
           LAD=IJADD(I*(I+1)/2)+KADD(I+KOS)+3
           II=IA+WTY(M4)
           CJJ=C(JJ)*2.0D+00
           DO 2001 J=M4,M3
                LAD1=LAD+LADD(J+LOS)
                INT(LAD1)=INT(LAD1)+CJJ*C(II)
 2001      II=II+1
 2000 JJ=JJ+1
      RETURN
C
C   XX   NO 13 & 161 FOUR EXTERNALS
C***********************************************************************
  315 M5=IMAX
      M6=IMIN
      IF (M6.LT.3) M6=3
clj IBM Fortran vers. 2 objects (level 4 error) the following branch:
clj   IF (M5.LT.M6) GO TO 2004
      IF (M5.LT.M6) GO TO 3209
      DO 2004 I=M6,M5
           IF(I.GT.JMAXG.OR.I.LT.JMING) GO TO 2004
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
                CJJ=C(JJ)*2.0D+00
                DO 2006 K=M4,M3
                     LAD2=LAD1+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD2+1)=INT(LAD2+1)-Z
                     INT(LAD2+2)=INT(LAD2+2)+Z
 2006           II=II+1
 2005      JJ=JJ+1
 2004 CONTINUE
 3209 continue
C
C     NO 161
clj   IF (N.LT.3) GO TO 7084
      IF (N.LT.3) GO TO 3998
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
                IF(J.GT.IMAXG.OR.J.LT.IMING.OR.J.GT.JMAXG.OR.J.LT.JMING)
     *          GO TO 7085
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 7085
                LOS=OS(SS(J)+1)
                LAD3=IJADD(J*(J+1)/2)+KADD(J+KOS)+3
                II=IIA+WTY(M4)
                CJJ=C(JJ)*2.0D+00
                DO 7086 K=M4,M3
                     INT(LAD3+LADD(K+LOS))=INT(LAD3+LADD(K+LOS))+CJJ*C(
     *               II)
 7086           II=II+1
 7085      JJ=JJ+1
 7084 CONTINUE
 3998 continue
C
C   NO  14
C
      M5=IMAX
      M6=IMIN
      IF (M6.LT.3) M6=3
clj   IF (M5.LT.M6) GO TO 2007
      IF (M5.LT.M6) GO TO 3997
      DO 2007 I=M6,M5
           IF(I.GT.IMAXG.OR.I.LT.IMING) GO TO 2007
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
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2008
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 2008
                KSM=IEOR(SS(I),SS(J))
                KOS=OS(KSM+1)
                LAD1=IJADD(I*(I-1)/2+J)
                II=IA+WTX(J,IS+1)+WTY(M4)
                JJ=JJA+WTY(M4)
                DO 2009 K=M4,M3
                     LSM=IEOR(KSM,SS(K))
                     LOS=OS(LSM+1)
                     LAD2=LAD1+KADD(K+KOS)+LADD(K+LOS)
                     Z=C(II)*C(JJ)*2.0D+00
                     INT(LAD2+1)=INT(LAD2+1)-Z
                     INT(LAD2+2)=INT(LAD2+2)+Z
                     II=II+1
 2009           JJ=JJ+1
 2008      CONTINUE
           IF(I.GT.JMAXG.OR.I.LT.JMING) GO TO 2007
           M3=MX(SSJ+1)
           DO 2083 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 2083
                II=IA+WTX(J,IS+1)+WTY(M4)
                JJ=JJA+WTY(M4)
                Z1=0.0D+00
                DO 2082 K=M4,M3
                     Z1=Z1+C(II)*C(JJ)
                     II=II+1
 2082           JJ=JJ+1
                INT(LAD3+LADD(J+LOS1))=INT(LAD3+LADD(J+LOS1))+Z1+Z1
 2083      CONTINUE
 2007 CONTINUE
 3997 continue
C
C    NO 16   BOTH
C
clj   IF (M5.LT.M6) GO TO 2010
      IF (M5.LT.M6) GO TO 3996
      DO 2010 I=M6,M5
           IF(I.GT.IMAXG.OR.I.LT.IMING) GO TO 2010
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
           LAD2=IJADD(I*(I+1)/2)+KADD(I+KOS)+3
           DO 2011 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 2011
                CJJ=C(JJ)*2.0D+00
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2071
                KSM=IEOR(SS(I),SS(J))
                KOS=OS(KSM+1)
                LAD3=IJADD(I*(I-1)/2+J)+KADD(J+KOS)
                II=IA+WTX(J,IS+1)+WTY(M4)
                DO 2012 K=M4,M3
                     LAD5=LAD3+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD5+2)=INT(LAD5+2)-Z
                     INT(LAD5+1)=INT(LAD5+1)+Z
 2012           II=II+1
 2071           IF(I.GT.JMAXG.OR.I.LT.JMING) GO TO 2011
                II=IA+WTX(J,IS+1)+WTY(M4)
                DO 2072 K=M4,M3
                     INT(LAD2+LADD(K+LOS))=INT(LAD2+LADD(K+LOS))-CJJ*C(
     *               II)
 2072           II=II+1
 2011      JJ=JJ+1
 2010 CONTINUE
 3996 continue
C
C    NO  20
C
clj   IF (M5.LT.4) GO TO 2013
      IF (M5.LT.4) GO TO 3995
      M6=IMIN
      IF (M6.LT.4) M6=4
      DO 2013 I=M6,M5
           IF(I.GT.IMAXG.OR.I.LT.IMING) GO TO 2013
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
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2014
                KSM=IEOR(SS(I),SS(J))
                KOS=OS(KSM+1)
                LAD1=IJADD(I*(I-1)/2+J)
                CJJ=C(JJ)*2.0D+00
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
                          Z=CJJ*C(II)
                          INT(LAD3+1)=INT(LAD3+1)+Z
                          INT(LAD3+3)=INT(LAD3+3)-Z
 2016                II=II+1
 2015           CONTINUE
 2014      JJ=JJ+1
 2013 CONTINUE
 3995 continue
C
C      NO  21
C
clj   IF(M5.LT.M6) GO TO 2017
      IF(M5.LT.M6) GO TO 3994
      DO 2017 I=M6,M5
           IF(I.GT.IMAXG.OR.I.LT.IMING) GO TO 2017
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 2017
           JJA=JA+WTX(I,JS+1)
           DO 2018 J=3,I-1
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2018
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
                     CII=C(II)*2.0D+00
                     DO 2020 L=M4,M3
                          LAD3=LAD2+LADD(L+LOS)
                          Z=CII*C(JJ)
                          INT(LAD3+1)=INT(LAD3+1)-Z
                          INT(LAD3+2)=INT(LAD3+2)+Z
 2020                JJ=JJ+1
 2019           II=II+1
 2018      CONTINUE
 2017 CONTINUE
 3994 continue
C
C   NO  22
C
      IF(M5.LT.M6) RETURN
      DO 2021 I=M6,M5
           IF(I.GT.IMAXG.OR.I.LT.IMING) GO TO 2021
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 2021
           JJA=JA+WTX(I,JS+1)
           DO 2022 J=3,I-1
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2022
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
                     CJJ=C(JJ)*2.0D+00
                     DO 2024 L=M4,M3
                          LAD3=LAD2+LADD(L+LOS)
                          Z=CJJ*C(II)
                          INT(LAD3+3)=INT(LAD3+3)-Z
                          INT(LAD3+2)=INT(LAD3+2)+Z
 2024                II=II+1
 2023           JJ=JJ+1
 2022      CONTINUE
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
           IF(I.GT.JMAXG.OR.I.LT.JMING)GO TO 2031
           LAD=IJADD(I*(I+1)/2)+1
           CJJ=C(JA+WAB(I))+C(JA+WAB(I))
           KOS=OS(1)
           DO 2032 J=1,I-1
                LAD1=KADD(J+KOS)+LADD(J+OS(SS(J)+1))+LAD
                INT(LAD1)=INT(LAD1)+C(IA+WAB(J))*CJJ
 2032      CONTINUE
 2031 CONTINUE
C
C    NO   72A
C
 6586 IF (JS.NE.0) GO TO  6635
      IF (M5.LT.M6) GO TO 6635
      DO 2033 I=M6,M5
           IF(I.GT.JMAXG.OR.I.LT.JMING) GO TO 2033
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
           CJJ=C(JJ)*2.0D+00*SQRT2
           DO 2034 J=M2,M1
                LAD1=LAD+LADD(J+LOS)
                Z=CJJ*C(II)
                INT(LAD1+3)=INT(LAD1+3)+Z
                INT(LAD1+1)=INT(LAD1+1)+Z
 2034      II=II+1
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
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2036
                IF (SSJ.GT.SS(J)) GO TO 2036
                IF (SSJ.EQ.SS(J)) M3=J-1
                IF (M3.LT.M4) GO TO 2036
                JJ=JJA+WTY(M4)
                II=IA+WAB(J)
                KSM=IEOR(SS(I),SS(J))
                KOS=OS(KSM+1)
                LOS=OS(SS(I)+1)
                LAD1=IJADD(I*(I-1)/2+J)+KADD(J+KOS)
                CII=C(II)*SQRT2*2.0D+00
                DO 2037 K=M4,M3
                     LAD2=LAD1+LADD(K+LOS)
                     INT(LAD2+1)=INT(LAD2+1)+CII*C(JJ)
 2037           JJ=JJ+1
 2036      CONTINUE
 2035 CONTINUE
C
C   NO  74A
C
 6636 IF (JS.NE.0) GO TO 6637
      IF (M5.LT.M6) GO TO 6637
      DO 2038 I=M6,M5
           IF(I.GT.JMAXG.OR.I.LT.JMING) GO TO 2038
           JJ=JA+WAB(I)
           LAD=IJADD(I*(I+1)/2)
           CJJ=C(JJ)*2.0D+00*SQRT2
           KOS=OS(1)
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
                     INT(LAD2+1)=INT(LAD2+1)+CJJ*C(II)
 2040           II=II+1
 2039      CONTINUE
 2038 CONTINUE
C
C   NO  76
C
 6637 IF (IS.NE.0) GO TO 6638
clj   IF (M5.LT.M6) GO TO 2041
      IF (M5.LT.M6) GO TO 3993
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
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2042
                KSM=IEOR(SS(I),SS(J))
                KOS=OS(KSM+1)
                LAD1=IJADD(I*(I-1)/2+J)
                CJJ=C(JJ)*2.0D+00*SQRT2
                DO 2043 K=1,J-1
                     II=IA+WAB(K)
                     LSM=IEOR(KSM,SS(K))
                     LOS=OS(LSM+1)
                     LAD2=LAD1+KADD(K+KOS)+LADD(K+LOS)
                     INT(LAD2+1)=INT(LAD2+1)+C(II)*CJJ
 2043           CONTINUE
 2042      JJ=JJ+1
 2041 CONTINUE
 3993 continue
C
C    NO  78
C
      M6=IMIN
      IF (M6.LT.2) M6=2
      IF (M5.LT.M6) GO TO 6638
      DO 2044 I=M6,M5
           IF(I.GT.JMAXG.OR.I.LT.JMING) GO TO 2044
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
                Z=C(II)*C(JJ)*SQRT2*2.0D+00
                INT(LAD1+3)=INT(LAD1+3)+Z
                INT(LAD1+2)=INT(LAD1+2)+Z
 2045      JJ=JJ+1
 2044 CONTINUE
C
C    NO   84A
C
 6638 M6=IMIN
      IF (M6.LT.3) M6=3
clj   IF (M5.LT.M6) GOTO 2054
      IF (M5.LT.M6) GOTO 3992
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
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2090
                II=IA+WTW(J,IS+1)+WTY(M4)
                JJ=JJA+WTY(M4)
                KSM=IEOR(SS(I),SS(J))
                KOS=OS(KSM+1)
                LAD1=IJADD(I*(I-1)/2+J)
                DO 2050 K=M4,M3
                     LSM=IEOR(KSM,SS(K))
                     LOS=OS(LSM+1)
                     LAD2=LAD1+KADD(K+KOS)+LADD(K+LOS)
                     Z=C(II)*C(JJ)*2.0D+00
                     INT(LAD2+1)=INT(LAD2+1)+Z
                     INT(LAD2+2)=INT(LAD2+2)+Z
                     II=II+1
 2050           JJ=JJ+1
 2090           IF(I.GT.JMAXG.OR.I.LT.JMING) GO TO 2049
                II=IA+WTW(J,IS+1)+WTY(M4)
                JJ=JJA+WTY(M4)
                Z1=0.0D+00
                DO 2084 K=M4,M3
                     Z1=Z1+C(II)*C(JJ)
                     II=II+1
 2084           JJ=JJ+1
                INT(LAD3+LADD(J+LOS1))=INT(LAD3+LADD(J+LOS1))+Z1+Z1
 2049      CONTINUE
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
           LAD1=IJADD(I*(I+1)/2)+KADD(I+KOS1)+3
           DO 2052 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 2052
                CJJ=C(JJ)*2.0D+00
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2085
                II=IA+WTW(J,IS+1)+WTY(M4)
                KSM=IEOR(SS(I),SS(J))
                KOS=OS(KSM+1)
                LAD2=IJADD(I*(I-1)/2+J)+KADD(J+KOS)
                DO 2053 K=M4,M3
                     LADX=LADD(K+LOS1)
                     Z=CJJ*C(II)
C     INT(LAD1+LADX+3)=INT(LAD1+LADX+3)+Z
                     INT(LAD2+LADX+2)=INT(LAD2+LADX+2)+Z
                     INT(LAD2+LADX+1)=INT(LAD2+LADX+1)+Z
 2053           II=II+1
 2085           IF(I.GT.JMAXG.OR.I.LT.JMING) GO TO 2052
                II=IA+WTW(J,IS+1)+WTY(M4)
                DO 2086 K=M4,M3
                     INT(LAD1+LADD(K+LOS1))=INT(LAD1+LADD(K+LOS1))+CJJ*C
     *               (II)
 2086           II=II+1
 2052      JJ=JJ+1
 2051 CONTINUE
C
C   NO   89 & 160
C
      DO 2054 I=M6,M5
           IF(I.GT.JMAXG.OR.I.LT.JMING)GO TO 2054
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
                CJJ=C(JJ)*2.0D+00
                DO 2056 K=M4,M3
                     LAD2=LAD1+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD2+1)=INT(LAD2+1)+Z
                     INT(LAD2+2)=INT(LAD2+2)+Z
 2056           II=II+1
 2055      JJ=JJ+1
 2054 CONTINUE
 3992 continue
C
C  NO 160
clj   IF (N.LT.3) GO TO 7087
      IF (N.LT.3) GO TO 3991
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
                IF(J.GT.IMAXG.OR.J.LT.IMING.OR.J.GT.JMAXG.OR.J.LT.JMING)
     *          GO TO 7088
                LOS=OS(SS(J)+1)
                LAD3=IJADD(J*(J+1)/2)+KADD(J+KOS)+3
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 7088
                II=IIA+WTY(M4)
                CJJ=C(JJ)*2.0D+00
                DO 7089 K=M4,M3
C     LADX=LAD3+LADD(K+LOS)
                     INT(LAD3+LADD(K+LOS))=INT(LAD3+LADD(K+LOS))+CJJ*C(
     *               II)
 7089           II=II+1
 7088      JJ=JJ+1
 7087 CONTINUE
 3991 continue
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
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2058
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
                     CJJ=C(JJ)*2.0D+00
                     DO 2060 L=M4,M3
                          LAD3=LAD2+LADD(L+LOS)
                          Z=CJJ*C(II)
                          INT(LAD3+3)=INT(LAD3+3)+Z
                          INT(LAD3+2)=INT(LAD3+2)+Z
 2060                II=II+1
 2059           JJ=JJ+1
 2058      CONTINUE
 2057 CONTINUE
C
C     NO   100A
C
      DO 2061 I=M6,M5
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 2061
           JJA=JA+WTW(I,JS+1)
           DO 2062 J=3,I-1
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2062
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
                     CII=C(II)*2.0D+00
                     DO 2064 L=M4,M3
                          LAD3=LAD2+LADD(L+LOS)
                          Z=CII*C(JJ)
                          INT(LAD3+1)=INT(LAD3+1)+Z
                          INT(LAD3+2)=INT(LAD3+2)+Z
 2064                JJ=JJ+1
 2063           II=II+1
 2062      CONTINUE
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
                IF(J.GT.JMAXG.OR.J.LT.JMING) GO TO 2066
                KSM=IEOR(SS(I),SS(J))
                KOS=OS(KSM+1)
                LAD1=IJADD(I*(I-1)/2+J)
                CJJ=C(JJ)*2.0D+00
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
                          Z=CJJ*C(II)
                          INT(LAD3+1)=INT(LAD3+1)+Z
                          INT(LAD3+3)=INT(LAD3+3)+Z
 2068                II=II+1
 2067           CONTINUE
 2066      JJ=JJ+1
 2065 CONTINUE
      RETURN
C
C***********************************************************************
c3-24-88  the below must be debug and is now commented out with 'cxc'
cxc999 WRITE (6,998) M,ISEG
cxc998 FORMAT (2I4)
cxc    RETURN
C
C     CLOSED INTERNAL LOOP
C***********************************************************************
cxc318 PRINT *,'  YOU SHOULD NOT HAVE M=18 IN SHAPE4!'
cxc    CALL DUMP
cxc318 Z=VAL1*(INT(ARR+TR1)+VAL2*INT(ARR+TR2))
cxc    II=IA
cxc    JJ=JA
cxc    Z=0.0D+00
cxc    DO 4000 I=1,N1
cxc         Z=Z+C(II)*C(JJ)
C GF  S(JJ)=S(JJ)+Z*C(II)
C GF  S(II)=S(II)+Z*C(JJ)
cxc         II=II+1
cxc4000 JJ=JJ+1
cxc    INT(ARR+TR1)=INT(ARR+TR1)+Z*VAL1*2.0D+00
cxc    INT(ARR+TR2)=INT(ARR+TR2)+Z*VAL1*VAL2*2.0D+00
cxc    RETURN
C
C     ----- INTERNAL CASE FOR TRACKS OF (3,2,1) -----
C***********************************************************************
cxc319 PRINT *,'  YOU SHOULD NOT HAVE M=19 IN SHAPE4!'
cxc      CALL DUMP
C GF  Z=VAL1*(INT(ARR+1)+INT(ARR+2)+INT(ARR+3))
cxc      Z=0.0D+00
cxc      II=IA
cxc      JJ=JA
cxc      DO 3999 I=1,N1
cxc           Z=Z+C(II)*C(JJ)
C GF  S(JJ)=S(JJ)+Z*C(II)
C GF  S(II)=S(II)+Z*C(JJ)
cxc           II=II+1
cxc 3999 JJ=JJ+1
cxc      Z=Z*VAL1*2.0D+00
cxc      INT(ARR+1)=INT(ARR+1)+Z
cxc      INT(ARR+2)=INT(ARR+2)+Z
cxc      INT(ARR+3)=INT(ARR+3)+Z
cxc      RETURN
C
C***********************************************************************
C
cxc 317 Z=VAL1*(INT(ARR+TR1)+VAL2*INT(ARR+TR2))
C   ZZ   ARRIVAL,FOR COMPLETENESS
cxc      STOP ' ZZ ARRIVAL IN EXTERNAL'
cxcC     RETURN
c3-24-88  the above must be debug and is now commented out with 'cxc'
      END
