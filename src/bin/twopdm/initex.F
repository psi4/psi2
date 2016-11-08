      SUBROUTINE INITEX(INT,C,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
      IMPLICIT REAL*8 (A-H,O-Z)
C      EXTENDED DUMMY INT,C
      INTEGER ARR,TR1,TR2,ASM,AOS,OS,WTW,WTX,WTY,WAB,SS,SSI,SSJ,SYMORB
      INTEGER BMAX,ORBFRM
      integer nmax,nwksmx,numij,nsym,norbs
      REAL*8 INT(NMAX),C(NWKSMX)
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
     *,             ITAPE3,ITAP05,ITAP06
      COMMON /ALL/VAL1,VAL2,VAL3,ARR,TR1,TR2,IA,JA,M,ISEG,N,N1,N2
     *,           IMAX,IMIN
      COMMON /MINMAX/ IMING,IMAXG,JMING,JMAXG
      REAL*8 Z,VAL1,VAL2,VAL3
      SQRT2=DSQRT(2.0D+00)
      SQRT3=DSQRT(3.0D+00)
      SQT1P5=DSQRT(1.5D+00)
      RETURN
C
C
      ENTRY SHAPES(INT,C,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
CMC   WRITE(ITAPE3,*) 'JUST ENTERED ENTRY SHAPES'
CMC   WRITE(ITAPE3,19010) M,ISEG
CMC10 FORMAT(/' M AND ISEG AS WE ENTER SHAPES',2I4)
C
C
      AOS=OS(ASM+1)
cets  GO TO (313,305,304,312,311,310,307,308,303,306,309,302,301,314
C GF  GO TO (999,999,304,999,999,310,307,308,303,306,309,302,301,314
C     GO TO (301,302,303,304,305,306,307,308,309,310,311,312,313,314
C     GO TO ( YZ, XZ, XY, YY, XX, ZY, YX, YW, WZ, WY, WX, XW, WW,---YY
C
cets 1,314,314,317,318,319),M
C       XX, WW, ZZ,---4 INTERNALS),M
      if(m.le.14)
     #GO TO (313,305,304,312,311,310,307,308,303,306,309,302,301,314),m
      goto (314,314,317,318,319),(m-14)
C
      PRINT *,'  FUNNY M:',M
      CALL DUMP
C
C   WX   NO 39  ISEG=7
C***********************************************************************
  311 IF (ISEG.NE.7) GO TO 999
      VAL4=SQT1P5*VAL1*2.0D+00
      IF (IS.NE.0) GO TO 6616
      DO 1201 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 1201
           II=IA+WAB(I)
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1201
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           JJ=JA+WTY(M2)+WTX(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+1
           CII=C(II)*VAL4
           DO 1202 J=M2,M1
                INT(LAD+LADD(J+LOS))=INT(LAD+LADD(J+LOS))-CII*C(JJ)
 1202      JJ=JJ+1
 1201 CONTINUE
C
C NO 129
C
 6616 VAL4=SQRT3*0.5D+00*VAL1*2.0D+00
      DO 5070 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SSI,JS)
           IF (SSJ.GT.SS(I)) GO TO 5070
           IF (SSI.GT.SSJ) GO TO 5070
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 5070
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 5070
           IIA=IA+WTW(I,IS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 5071 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 5071
                JJ=JA+WTX(J,JS+1)+WTY(M4)
                II=IIA+WTY(M4)
                LAD1=LAD+LADD(J+LOS)
                Z=0.0D+00
                DO 5072 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 5072           JJ=JJ+1
                INT(LAD1+1)=INT(LAD1+1)+VAL4*Z
 5071      CONTINUE
 5070 CONTINUE
C
C NO 132
C
      DO 5073 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSJ.GT.SSI) GO TO 5073
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 5073
           IF (SSI.GT.SS(I)) GO TO 5073
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 5073
           JJA=JA+WTX(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+TR1
           DO 5074 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 5074
                II=IA+WTW(J,IS+1)+WTY(M4)
                JJ=JJA+WTY(M4)
                Z=0.0D+00
                DO 5075 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 5075           JJ=JJ+1
                INT(LAD+LADD(J+LOS))=INT(LAD+LADD(J+LOS))+VAL4*Z
 5074      CONTINUE
 5073 CONTINUE
C
C   NO 40
C
      VAL4=SQT1P5*VAL1*2.0D+00
      IF (IS.NE.0) GO TO 6619
      DO 1203 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 1203
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1203
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+TR1
           DO 1204 J=M2,M1
                LAD1=LAD+LADD(J+LOS)
                INT(LAD1)=INT(LAD1)+C(IA+WAB(J))*C(JJ)*VAL4
                JJ=JJ+1
 1204      CONTINUE
 1203 CONTINUE
C
C   NO 45 AND 145
C
 6619 VAL4=-VAL1*SQRT3*0.5D+00*2.0D+00
      IF (IS.NE.JS) GO TO 6620
      DO 1205 I=2,N
           SSI=IEOR(SS(I),IS)
           IF (SSI.GT.SS(I)) GO TO 1205
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1205
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           II=IA+WTW(I,IS+1)+WTY(M2)
           JJ=JA+WTX(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
           Z1=0.0D+00
           DO 1206 J=M2,M1
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
                Z=VAL4*C(II)*C(JJ)
C****************************************
      IF(LAD.EQ.136)WRITE(ITAPE3,*)' Z,VAL4=',Z,VAL4
                Z1=Z1+Z
                INT(LAD+1)=INT(LAD+1)+Z
                II=II+1
 1206      JJ=JJ+1
           INT(LAD1+1)=INT(LAD1+1)-Z1
C****************************************
      IF(LAD1.EQ.136)WRITE(ITAPE3,*)' Z1,VAL4=',Z1,VAL4
 1205 CONTINUE
C
C    NO 49
C
 6620 VAL4=VAL1*SQRT3*0.5D+00*2.0D+00
      DO 1207 I=3,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 1207
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1207
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SSJ) GO TO 1207
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 1207
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+TR1
           DO 1208 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1208
                II=IA+WTW(J,IS+1)+WTY(M4)
                CJJ=C(JJ)*VAL4
                DO 1209 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+CJJ*C(II)
 1209           II=II+1
 1208      JJ=JJ+1
 1207 CONTINUE
C
C   NO 52
C
      VAL4=-VAL1*SQRT3*0.5D+00*2.0D+00
      DO 1210 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SS(I),JS)
           IF (SSI.GT.SSJ) GO TO 1210
           IF (SSJ.GT.SS(I)) GO TO 1210
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1210
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 1210
           IIA=IA+WTW(I,IS+1)
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           DO 1211 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1211
                II=IIA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+TR1
                CJJ=C(JJ)*VAL4
                DO 1212 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+CJJ*C(II)
 1212           II=II+1
 1211      JJ=JJ+1
 1210 CONTINUE
C
C    NO 53
C
      DO 1213 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SSI) GO TO 1213
           IF (SSI.GT.SS(I)) GO TO 1213
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1213
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 1213
           II=IA+WTW(I,IS+1)+WTY(M2)
           JJA=JA+WTX(I,JS+1)
           DO 1214 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1214
                JJ=JJA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+1
                CII=C(II)*VAL4
                DO 1215 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+C(JJ)*CII
 1215           JJ=JJ+1
 1214      II=II+1
 1213 CONTINUE
C
C NO 170
C
      VAL4=-VAL1*SQRT3*0.5D+00*2.0D+00
      DO 8090 I=3,N
           SSI=IEOR(SS(I),IS)
           IF (SSI.GT.SS(I)) GO TO 8090
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 8090
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           SSJ=IEOR(JS,SSI)
           IF (SSJ.GT.SSI) GO TO 8090
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 8090
           II=IA+WTW(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+1
           DO 8091 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 8091
                JJ=JA+WTX(J,JS+1)+WTY(M4)
                CII=C(II)*VAL4
                DO 8092 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+CII*C(JJ)
 8092           JJ=JJ+1
 8091      II=II+1
 8090 CONTINUE
      RETURN
C
C   XW   NO 39A   ISEG=9
C***********************************************************************
  312 IF (ISEG.NE.9) GO TO 999
      VAL4=-VAL1*SQT1P5*2.0D+00
      IF (JS.NE.0) GO TO 6621
      DO 1301 I=2,N
           SSI=IEOR(SS(I),IS)
           IF (SSI.GT.SS(I)) GO TO 1301
           JJ=JA+WAB(I)
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1301
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           II=IA+WTY(M2)+WTX(I,IS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+TR1
           CJJ=C(JJ)*VAL4
           DO 1302 J=M2,M1
                INT(LAD+LADD(J+LOS))=INT(LAD+LADD(J+LOS))+CJJ*C(II)
 1302      II=II+1
 1301 CONTINUE
C
C NO 130
C
 6621 VAL4=VAL1*SQRT3*0.5D+00*2.0D+00
      DO 9060 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(JS,SSI)
           IF (SSI.GT.SSJ) GO TO 9060
           IF (SSJ.GT.SS(I)) GO TO 9060
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 9060
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 9060
           IIA=IA+WTX(I,IS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+1
           DO 9061 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9061
                II=IIA+WTY(M4)
                JJ=JA+WTW(J,JS+1)+WTY(M4)
                Z=0.0D+00
                DO 9062 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9062           JJ=JJ+1
                INT(LAD+LADD(J+LOS))=INT(LAD+LADD(J+LOS))+VAL4*Z
 9061      CONTINUE
 9060 CONTINUE
C
C NO 131
C
      DO 9063 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SS(I)) GO TO 9063
           IF (SSJ.GT.SSI) GO TO 9063
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 9063
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 9063
           JJA=JA+WTW(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+TR1
           DO 9064 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9064
                II=IA+WTX(J,IS+1)+WTY(M4)
                JJ=JJA+WTY(M4)
                Z=0.0D+00
                DO 9065 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9065           JJ=JJ+1
                INT(LAD+LADD(J+LOS))=INT(LAD+LADD(J+LOS))+VAL4*Z
 9064      CONTINUE
 9063 CONTINUE
C
C   NO 40A
C
      VAL4=SQT1P5*VAL1*2.0D+00
      IF (JS.NE.0) GO TO 6622
      DO 1303 I=2,N
           SSI=IEOR(SS(I),IS)
           IF (SSI.GT.SS(I)) GO TO 1303
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1303
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           II=IA+WTX(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+1
           DO 1304 J=M2,M1
                INT(LAD+LADD(J+LOS))=INT(LAD+LADD(J+LOS))+VAL4*C(II)*C(
     *          JA+WAB(J))
 1304      II=II+1
 1303 CONTINUE
C
C   NO 45A AND 145A
C
 6622 VAL4=-VAL1*SQRT3*0.5D+00*2.0D+00
      IF (IS.NE.JS) GO TO 6623
      DO 1305 I=2,N
           SSI=IEOR(SS(I),IS)
           IF (SSI.GT.SS(I)) GO TO 1305
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1305
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
           Z1=0.0D+00
           II=IA+WTX(I,IS+1)+WTY(M2)
           JJ=JA+WTW(I,JS+1)+WTY(M2)
           DO 1306 J=M2,M1
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
                Z=VAL4*C(II)*C(JJ)
C****************************************
      IF(LAD.EQ.136)WRITE(ITAPE3,*)' Z,VAL4=',Z,VAL4
                Z1=Z1+Z
C****************************************
      IF(LAD1.EQ.136)WRITE(ITAPE3,*)' Z1,VAL4=',Z1,VAL4
                INT(LAD+1)=INT(LAD+1)+Z
                JJ=JJ+1
 1306      II=II+1
           INT(LAD1+1)=INT(LAD1+1)-Z1
 1305 CONTINUE
C
C   NO 49A
C
 6623 VAL4=VAL1*SQRT3*0.5D+00*2.0D+00
      DO 1307 I=3,N
           SSI=IEOR(SS(I),IS)
           IF (SSI.GT.SS(I)) GO TO 1307
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1307
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           SSJ=IEOR(JS,SSI)
           IF (SSJ.GT.SSI) GO TO 1307
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 1307
           II=IA+WTX(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+1
           DO 1308 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1308
                JJ=JA+WTW(J,JS+1)+WTY(M4)
                CII=C(II)*VAL4
                DO 1309 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+CII*C(JJ)
 1309           JJ=JJ+1
 1308      II=II+1
 1307 CONTINUE
C
C   NO 52A
C
      VAL4=-VAL1*SQRT3*0.5D+00*2.0D+00
      DO 1310 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SSI) GO TO 1310
           IF (SSI.GT.SS(I)) GO TO 1310
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1310
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 1310
           II=IA+WTX(I,IS+1)+WTY(M2)
           JJA=JA+WTW(I,JS+1)
           DO 1311 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1311
                JJ=JJA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+1
                CII=C(II)*VAL4
                DO 1312 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+C(JJ)*CII
 1312           JJ=JJ+1
 1311      II=II+1
 1310 CONTINUE
C
C   NO 53A
C
      DO 1313 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SS(I),JS)
           IF (SSI.GT.SSJ) GO TO 1313
           IF (SSJ.GT.SS(I)) GO TO 1313
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1313
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 1313
           IIA=IA+WTX(I,IS+1)
           JJ=JA+WTW(I,JS+1)+WTY(M2)
           DO 1314 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1314
                II=IIA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+TR1
                CJJ=C(JJ)*VAL4
                DO 1315 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+CJJ*C(II)
 1315           II=II+1
 1314      JJ=JJ+1
 1313 CONTINUE
C
C NO 171
C
      VAL4=-VAL1*SQRT3*0.5D+00*2.0D+00
      DO 8094 I=3,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 8094
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 8094
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SSJ) GO TO 8094
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 8094
           JJ=JA+WTW(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+TR1
           DO 8095 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 8095
                II=IA+WTX(J,IS+1)+WTY(M4)
                CJJ=C(JJ)*VAL4
                DO 8096 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+CJJ*C(II)
 8096           II=II+1
 8095      JJ=JJ+1
 8094 CONTINUE
      RETURN
C
C   WW    NO 68  ISEG=4,5
C***********************************************************************
  313 IF (ISEG.EQ.8) RETURN
      IF (ISEG.EQ.6) GO TO 3143
      IF (ISEG.EQ.4.OR.ISEG.EQ.5) GO TO 299
      GO TO 999
  299 IF(IA.EQ.JA) GO TO 8501
      VAL4=-VAL1*SQRT2*2.0D+00
      IF (IS.NE.JS) GO TO 8501
      IF (IS.NE.0) GO TO 8501
      DO 1401 I=1,N
           LOS=OS(IEOR(ASM,SS(I))+1)
           LAD=KADD(I+AOS)+LADD(I+LOS)+ARR
           Z=VAL4*C(IA+WAB(I))*C(JA+WAB(I))
C****************************************
      IF(LAD.EQ.136) WRITE(ITAPE3,*) '  Z,VAL4 = ',Z,VAL4,INT(LAD+1)
           INT(LAD+1)=INT(LAD+1)+Z
           INT(LAD+2)=INT(LAD+2)-Z-Z
 1401 CONTINUE
C
C   NO 72A
C
 8501 IF (JS.NE.0) GO TO 6624
      DO 1402 I=2,N
           SSI=IEOR(SS(I),IS)
           IF (SSI.GT.SS(I)) GO TO 1402
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1402
           IF (SSI.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M1.LT.M2) GO TO 1402
           JJ=JA+WAB(I)
           II=IA+WTW(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           CJJ=C(JJ)*VAL1*2.0D+00
           DO 1403 J=M2,M1
                LAD1=LAD+LADD(J+LOS)
                Z=C(II)*CJJ
                INT(LAD1+1)=INT(LAD1+1)-Z
                INT(LAD1+2)=INT(LAD1+2)+Z+Z
 1403      II=II+1
 1402 CONTINUE
C
C  NO 78
C
 6624 IF (IS.NE.0) GO TO 6625
      VAL4=VAL1+VAL1
      DO 1404 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 1404
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1404
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M1.LT.M2) GO TO 1404
           JJ=JA+WTW(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 1405 J=M2,M1
                II=IA+WAB(J)
                LAD1=LAD+LADD(J+LOS)
                Z=VAL4*C(IA+WAB(J))*C(JJ)
                INT(LAD1+1)=INT(LAD1+1)-Z
                INT(LAD1+2)=INT(LAD1+2)+Z+Z
 1405      JJ=JJ+1
 1404 CONTINUE
C
C   NO 81 AND 148
C
 6625 VAL4=-2.0D+00*VAL1/SQRT2
      IF (IS.NE.JS) GO TO 8502
      IF(IA.EQ.JA) GO TO 8502
      DO 1406 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 1406
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1406
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M1.LT.M2) GO TO 1406
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
           Z1=0.0D+00
           II=IA+WTW(I,IS+1)+WTY(M2)
           JJ=JA+WTW(I,IS+1)+WTY(M2)
           DO 1407 J=M2,M1
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
                Z=VAL4*C(II)*C(JJ)
C****************************************
      IF(LAD.EQ.136) WRITE(ITAPE3,*) '  Z,VAL4 = ',Z,VAL4,INT(LAD+1)
                Z1=Z1+Z
                INT(LAD+1)=INT(LAD+1)+Z
                INT(LAD+2)=INT(LAD+2)-Z-Z
                II=II+1
 1407      JJ=JJ+1
           INT(LAD1+1)=INT(LAD1+1)+Z1
           INT(LAD1+2)=INT(LAD1+2)-Z1-Z1
 1406 CONTINUE
C
C     NO  88A
C
 8502 DO 1408 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SSJ) GO TO 1408
           IF (SSJ.GT.SS(I)) GO TO 1408
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1408
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 1408
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 1408
           JJ=JA+WTW(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 1409 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1409
                II=IA+WTW(J,IS+1)+WTY(M4)
                CJJ=C(JJ)*VAL4
                DO 1410 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD1+1)=INT(LAD1+1)+Z
                     INT(LAD1+2)=INT(LAD1+2)-Z-Z
 1410           II=II+1
 1409      JJ=JJ+1
 1408 CONTINUE
C
C    NO 89
C
      DO 1411 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(SS(I),IS)
           IF (SSJ.GT.SS(I)) GO TO 1411
           IF (SSI.GT.SSJ) GO TO 1411
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1411
           IF (SSJ.EQ.SS(I).OR.M1.GT.N) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 1411
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 1411
           JJ=JA+WTW(I,JS+1)+WTY(M2)
           IIA=IA+WTW(I,IS+1)
           DO 1412 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1412
                II=IIA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)
                CJJ=C(JJ)*VAL4
                DO 1413 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD1+1)=INT(LAD1+1)+Z
                     INT(LAD1+2)=INT(LAD1+2)-Z-Z
 1413           II=II+1
 1412      JJ=JJ+1
 1411 CONTINUE
C
C NO 133
C
      VAL4=-2.0D+00*VAL1/SQRT2
      DO 9000 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSJ.GT.SSI) GO TO 9000
           IF (SSI.GT.SS(I)) GO TO 9000
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 9000
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 9000
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 9000
           JJA=JA+WTW(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 9001 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9001
                JJ=JJA+WTY(M4)
                II=IA+WTW(J,IS+1)+WTY(M4)
                LAD1=LAD+LADD(J+LOS)
                Z=0.0D+00
                DO 9002 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     JJ=JJ+1
 9002           II=II+1
                Z=Z*VAL4
                INT(LAD1+1)=INT(LAD1+1)+Z
                INT(LAD1+2)=INT(LAD1+2)-Z-Z
 9001      CONTINUE
 9000 CONTINUE
      RETURN
C
C   WW  NO 68  ISEG=6
C
 3143 VAL4=-VAL1*SQRT2*2.0D+00
      IF(IA.EQ.JA) GO TO 8503
      IF (IS.NE.JS) GO TO 8503
      IF (IS.NE.0) GO TO 8503
      DO 1420 I=1,N
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=KADD(I+AOS)+LADD(I+LOS)+ARR
           Z=C(IA+WAB(I))*C(JA+WAB(I))*VAL4
           INT(LAD+1)=INT(LAD+1)+Z
           INT(LAD+2)=INT(LAD+2)-Z-Z
 1420 CONTINUE
C
C   NO 72
C
 8503 VAL4=VAL1+VAL1
      IF (IS.NE.0) GO TO 6626
      DO 1421 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 1421
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1421
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M1.LT.M2) GO TO 1421
           JJ=JA+WTW(I,JS+1)+WTY(M2)
           II=IA+WAB(I)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           CII=C(II)*VAL4
           DO 1422 J=M2,M1
                LAD1=LAD+LADD(J+LOS)
                Z=C(JJ)*CII
                INT(LAD1+1)=INT(LAD1+1)-Z
                INT(LAD1+2)=INT(LAD1+2)+Z+Z
 1422      JJ=JJ+1
 1421 CONTINUE
C
C   NO 72A
C
 6626 IF (JS.NE.0) GO TO 6627
      DO 1423 I=2,N
           SSI=IEOR(SS(I),IS)
           IF (SSI.GT.SS(I)) GO TO 1423
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1423
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M1.LT.M2) GO TO 1423
           II=IA+WTW(I,IS+1)+WTY(M2)
           JJ=JA+WAB(I)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           CJJ=C(JJ)*VAL4
           DO 1424 J=M2,M1
                LAD1=LAD+LADD(J+LOS)
                Z=C(II)*CJJ
                INT(LAD1+3)=INT(LAD1+3)-Z
                INT(LAD1+2)=INT(LAD1+2)+Z+Z
 1424      II=II+1
 1423 CONTINUE
C
C     NO 78A
C
      DO 1425 I=2,N
           SSI=IEOR(SS(I),IS)
           IF (SSI.GT.SS(I)) GO TO 1425
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1425
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M1.LT.M2) GO TO 1425
           II=IA+WTW(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 1426 J=M2,M1
                LAD1=LAD+LADD(J+LOS)
                Z=VAL4*C(II)*C(JA+WAB(J))
                INT(LAD1+1)=INT(LAD1+1)-Z
                INT(LAD1+2)=INT(LAD1+2)+Z+Z
 1426      II=II+1
 1425 CONTINUE
C
C    NO 78
C
 6627 IF (IS.NE.0) GO TO 6628
      DO 1427 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 1427
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1427
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M1.LT.M2) GO TO 1427
           JJ=JA+WTW(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 1428 J=M2,M1
                LAD1=LAD+LADD(J+LOS)
                Z=VAL4*C(IA+WAB(J))*C(JJ)
                INT(LAD1+3)=INT(LAD1+3)-Z
                INT(LAD1+2)=INT(LAD1+2)+Z+Z
 1428      JJ=JJ+1
 1427 CONTINUE
C
C    NO 81  AND 148
C
 6628 VAL4=-2.0D+00*VAL1/SQRT2
      IF (IS.NE.JS) GO TO 8504
      IF(IA.EQ.JA) GO TO 8504
      DO 1429 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 1429
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1429
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M1.LT.M2) GO TO 1429
           II=IA+WTW(I,IS+1)+WTY(M2)
           JJ=JA+WTW(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
           Z1=0.0D+00
           DO 1430 J=M2,M1
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
                Z=VAL4*C(II)*C(JJ)
                Z1=Z1+Z
                INT(LAD+1)=INT(LAD+1)+Z
                INT(LAD+2)=INT(LAD+2)-Z-Z
                II=II+1
 1430      JJ=JJ+1
           INT(LAD1+1)=INT(LAD1+1)+Z1
           INT(LAD1+2)=INT(LAD1+2)-Z1-Z1
 1429 CONTINUE
C
C   NO 88
C
 8504 DO 1431 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(JS,SSI)
           IF (SSJ.GT.SSI) GO TO 1431
           IF (SSI.GT.SS(I)) GO TO 1431
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1431
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 1431
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 1431
           II=IA+WTW(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 1432 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1432
                JJ=JA+WTW(J,JS+1)+WTY(M4)
                CII=C(II)*VAL4
                DO 1433 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=C(JJ)*CII
                     INT(LAD1+1)=INT(LAD1+1)+Z
                     INT(LAD1+2)=INT(LAD1+2)-Z-Z
 1433           JJ=JJ+1
 1432      II=II+1
 1431 CONTINUE
C
C    NO  88A
C
      DO 1435 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SSJ) GO TO 1435
           IF (SSJ.GT.SS(I)) GO TO 1435
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1435
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 1435
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 1435
           JJ=JA+WTW(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 1436 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1436
                II=IA+WTW(J,IS+1)+WTY(M4)
                CJJ=C(JJ)*VAL4
                DO 1437 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD1+3)=INT(LAD1+3)+Z
                     INT(LAD1+2)=INT(LAD1+2)-Z-Z
 1437           II=II+1
 1436      JJ=JJ+1
 1435 CONTINUE
C
C      NO  89
C
      DO 1438 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(SS(I),IS)
           IF (SSJ.GT.SS(I)) GO TO 1438
           IF (SSI.GT.SSJ) GO TO 1438
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 1438
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 1438
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 1438
           JJ=JA+WTW(I,JS+1)+WTY(M2)
           IIA=IA+WTW(I,IS+1)
           DO 1439 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1439
                II=IIA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)
                CJJ=C(JJ)*VAL4
                DO 1440 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD1+3)=INT(LAD1+3)+Z
                     INT(LAD1+2)=INT(LAD1+2)-Z-Z
 1440           II=II+1
 1439      JJ=JJ+1
 1438 CONTINUE
C
C   NO 89A
C
      DO 1441 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SSI) GO TO 1441
           IF (SSI.GT.SS(I)) GO TO 1441
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 1441
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 1441
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 1441
           II=IA+WTW(I,IS+1)+WTY(M2)
           JJA=JA+WTW(I,JS+1)
           DO 1442 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 1442
                JJ=JJA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)
                CII=C(II)*VAL4
                DO 1443 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=C(JJ)*CII
                     INT(LAD1+1)=INT(LAD1+1)+Z
                     INT(LAD1+2)=INT(LAD1+2)-Z-Z
 1443           JJ=JJ+1
 1442      II=II+1
 1441 CONTINUE
C
C NO 133
C
      VAL4=-2.0D+00*VAL1/SQRT2
      DO 9004 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSJ.GT.SSI) GO TO 9004
           IF (SSI.GT.SS(I)) GO TO 9004
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 9004
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 9004
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 9004
           JJA=JA+WTW(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 9005 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9005
                JJ=JJA+WTY(M4)
                II=IA+WTW(J,IS+1)+WTY(M4)
                LAD1=LAD+LADD(J+LOS)
                Z=0.0D+00
                DO 9006 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     JJ=JJ+1
 9006           II=II+1
                Z=VAL4*Z
                INT(LAD1+3)=INT(LAD1+3)+Z
                INT(LAD1+2)=INT(LAD1+2)-Z-Z
 9005      CONTINUE
 9004 CONTINUE
C
C NO 134
C
      IF(IA.EQ.JA) RETURN
      VAL4=-2.0D+00*VAL1/SQRT2
      DO 9007 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SSI,JS)
           IF (SSI.GT.SSJ) GO TO 9007
           IF (SSJ.GT.SS(I)) GO TO 9007
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 9007
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 9007
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 9007
           IIA=IA+WTW(I,IS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 9008 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9008
                JJ=JA+WTW(J,JS+1)+WTY(M4)
                II=IIA+WTY(M4)
                LAD1=LAD+LADD(J+LOS)
                Z=0.0D+00
                DO 9009 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9009           JJ=JJ+1
                Z=Z*VAL4
                INT(LAD1+1)=INT(LAD1+1)+Z
                INT(LAD1+2)=INT(LAD1+2)-Z-Z
 9008      CONTINUE
 9007 CONTINUE
      RETURN
  305 CALL XXDM(INT,C,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
      RETURN
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
      CJJ=C(JJ)*VAL1*2.0D+00
      DO 152 I=M2,M1
           LAD=ARR+LADD(I+AOS)
           Z=CJJ*C(II)
           INT(LAD+TR1)=INT(LAD+TR1)+Z
           INT(LAD+TR2)=INT(LAD+TR2)+VAL2*Z
  152 II=II+1
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
      CJJ=C(JJ)*VAL1*2.0D+00
      DO 151 I=M2,M1
           LAD=ARR+LADD(I+AOS)
           Z=CJJ*C(II)
           INT(LAD+3)=INT(LAD+3)+Z
           INT(LAD+1)=INT(LAD+1)+Z
  151 II=II+1
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
      CJJ=C(JJ)*VAL1*2.0D+00
      DO 150 I=M2,M1
           INT(ARR+LADD(I+AOS)+3)=INT(ARR+LADD(I+AOS)+3)+CJJ*C(II)
  150 II=II+1
      RETURN
C
C   XZ  NO 2 ISEG=13
C***********************************************************************
  302 IF (ISEG.NE.13) GO TO 999
      IF (JS.NE.0) RETURN
      JJ=JA
      CJJ=C(JJ)*VAL1*2.0D+00
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
                Z=CJJ*C(II)
                INT(LAD1+1)=INT(LAD1+1)+Z
                INT(LAD1+3)=INT(LAD1+3)-Z
  251      II=II+1
  250 CONTINUE
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
      VAL4=VAL1+VAL1
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
           CJJ=C(JJ)*VAL4
           DO 366 J=M4,M3
                INT(ARR+LADD(J+AOS)+3)=INT(ARR+LADD(J+AOS)+3)-CJJ*C(II)
  366      II=II+1
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
           Z=0.0D+00
           DO 5353 J=M4,M3
                Z=Z+C(II)*C(JJ)
                II=II+1
 5353      JJ=JJ+1
           LADX=ARR+LADD(I+AOS)+3
           INT(LADX)=INT(LADX)+Z*VAL4
 5352 CONTINUE
      RETURN
C
C  XY  NO 3 ISEG=16
C
 3031 SSI=IEOR(JS,IS)
      VAL4=VAL1+VAL1
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
           CJJ=C(JJ)*VAL4
           IF (SSI.EQ.JS) M3=I-1
           IF (M3.LT.M4) GO TO 363
           II=IA+WTX(I,IS+1)+WTY(M4)
           DO 364 J=M4,M3
                LAD=ARR+LADD(J+AOS)
                Z=C(II)*CJJ
                INT(LAD+3)=INT(LAD+3)-Z
                INT(LAD+1)=INT(LAD+1)-Z
  364      II=II+1
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
           Z=0.0D+00
           DO 6353 J=M4,M3
                Z=Z+C(II)*C(JJ)
                II=II+1
 6353      JJ=JJ+1
           INT(LAD+3)=INT(LAD+3)+VAL4*Z
           INT(LAD+1)=INT(LAD+1)+VAL4*Z
 6352 CONTINUE
      RETURN
C
C  XY  NO 3 ISEG=22
C
 3034 SSI=IEOR(JS,IS)
      VAL4=VAL1+VAL1
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
           CJJ=C(JJ)*VAL4
           DO 381 J=M4,M3
                LAD=ARR+LADD(J+AOS)
                Z=CJJ*C(II)
                INT(LAD+TR1)=INT(LAD+TR1)-Z
                INT(LAD+TR2)=INT(LAD+TR2)-VAL2*Z
  381      II=II+1
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
           Z=0.0D+00
           DO 7353 J=M4,M3
                Z=Z+C(II)*C(JJ)
                II=II+1
 7353      JJ=JJ+1
           Z=Z*VAL4
           INT(LAD+TR1)=INT(LAD+TR1)+Z
           INT(LAD+TR2)=INT(LAD+TR2)+Z*VAL2
 7352 CONTINUE
      RETURN
C
C  XY  NO 3  ISEG=3
C
 3033 SSI=IEOR(JS,IS)
      VAL4=VAL1+VAL1
      IF (SSI.GT.JS) GO TO 6604
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6604
      IF (M1.GT.N) M1=N
      M3=MX(SSI+1)
      M4=MN(SSI+1)
      IF (M4.GT.N) GO TO 6604
      IF (M2.LT.2) M2=2
      M2=MAX(M2,JMING)
      M1=MIN(M1,JMAXG)
      IF (M1.LT.M2) GO TO 6604
      JJ=JA+WTY(M2)
      DO 350 I=M2,M1
           IF (SSI.EQ.JS) M3=I-1
           IF (M3.LT.M4) GO TO 350
           II=IA+WTX(I,IS+1)+WTY(M4)
           KSM=IEOR(ASM,SS(I))
           KOS=OS(KSM+1)
           LAD=IJADD(ARR+I)+KADD(I+KOS)
           CJJ=C(JJ)*VAL4
           DO 351 J=M4,M3
                LAD1=LAD+LADD(J+AOS)
                Z=CJJ*C(II)
                INT(LAD1+2)=INT(LAD1+2)-Z
                INT(LAD1+1)=INT(LAD1+1)+Z
  351      II=II+1
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
      M2=MAX(M2,JMING)
      M1=MIN(M1,JMAXG)
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
                Z=VAL4*C(II)*C(JJ)
                INT(LAD1+1)=INT(LAD1+1)-Z
                INT(LAD1+2)=INT(LAD1+2)+Z
                II=II+1
  353      JJ=JJ+1
  352 CONTINUE
C
C    NO   5
C
 6605 JMIN=MAX(3,JMING)
      JMAX=MIN(N,JMAXG)
      DO 354 I=JMIN,JMAX
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
                CII=C(II)*VAL4
                DO 356 K=M4,M3
                     LAD2=LAD1+LADD(K+LOS)
                     Z=C(JJ)*CII
                     INT(LAD2+1)=INT(LAD2+1)-Z
                     INT(LAD2+2)=INT(LAD2+2)+Z
  356           JJ=JJ+1
  355      II=II+1
  354 CONTINUE
C
C    NO   6
C
      DO 357 I=JMIN,JMAX
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
                CJJ=C(JJ)*VAL4
                DO 359 K=M4,M3
                     LAD2=LAD1+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD2+3)=INT(LAD2+3)-Z
                     INT(LAD2+2)=INT(LAD2+2)+Z
  359           II=II+1
  358      JJ=JJ+1
  357 CONTINUE
C
C   NO 7
C
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) RETURN
      IF (M1.GT.N) M1=N
      IF (M2.LT.3) M2=3
      M1=MIN(M1,JMAXG)
      M2=MAX(M2,JMING)
      IF (M1.LT.M2) RETURN
      JJ=JA+WTY(M2)
      DO 360 I=M2,M1
           KSM=IEOR(ASM,SS(I))
           KOS=OS(KSM+1)
           LAD=IJADD(ARR+I)
           CJJ=C(JJ)*VAL4
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
                     Z=CJJ*C(II)
                     INT(LAD2+1)=INT(LAD2+1)+Z
                     INT(LAD2+3)=INT(LAD2+3)-Z
  362           II=II+1
  361      CONTINUE
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
      VAL4=-2.0D+00*VAL1/SQRT2
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
           CJJ=C(JJ)*VAL4
           DO 451 J=M4,M3
                LAD1=LAD+LADD(J+LOS)+1
                Z=C(II)*CJJ
                INT(LAD1)=INT(LAD1)+Z
                INT(LAD1+1)=INT(LAD1+1)-Z-Z
  451      II=II+1
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
           Z=VAL4*C(II)*C(JJ)
           INT(LAD+1)=INT(LAD+1)+Z
           INT(LAD+2)=INT(LAD+2)-Z-Z
           II=II+1
  460 JJ=JJ+1
      RETURN
C
C YY  NO 9 ISEG=8
C
 3042 VAL4=VAL1*SQT1P5*2.0D+00
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
           CJJ=C(JJ)*VAL4
           DO 421 J=M4,M3
                LAD1=LAD+LADD(J+LOS)+TR1
                INT(LAD1)=INT(LAD1)+CJJ*C(II)
  421      II=II+1
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
           INT(LAD+1)=INT(LAD+1)+VAL4*C(II)*C(JJ)
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
                INT(LAD1+1)=INT(LAD1+1)+VAL4*C(II)*C(JJ)
  441      JJ=JJ+1
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
      IF (DABS(XX).LT.1.0D-9) GO TO 490
      XX=XX+XX
      VAL5=SQRT2*VAL1*2.0D+00
      IF(IA.EQ.JA) GO TO 491
      IF (IS.NE.JS) GO TO 491
      DO 493 I=M2,M1
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
           Z=C(II)*C(JJ)
           INT(LAD+1)=INT(LAD+1)+XX*Z
           INT(LAD+2)=INT(LAD+2)+VAL5*Z
           II=II+1
  493 JJ=JJ+1
      GO TO 491
  490 VAL5=2.0D+00*SQRT2*VAL1
      IF(IA.EQ.JA) GO TO 492
      IF (IS.NE.JS) GO TO 492
      DO 494 I=M2,M1
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+LADD(I+LOS)+KADD(I+AOS)
           INT(LAD+2)=INT(LAD+2)+VAL5*C(II)*C(JJ)
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
           DO 496 J=M4,M3
                LAD1=LAD+LADD(J+LOS)
                Z=C(II)*CJJ
                INT(LAD1+1)=INT(LAD1+1)+XX*Z
                INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
  496      II=II+1
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
           CJJ=C(JJ)*VAL5
           DO 498 J=M4,M3
                LAD1=LAD+LADD(J+LOS)+2
                INT(LAD1)=INT(LAD1)+C(II)*CJJ
  498      II=II+1
  497 JJ=JJ+1
      RETURN
C
C   YY   NO 8  ISEG=6
C
 3041 XX=-VAL1/SQRT2+VAL3*SQT1P5
      IF (DABS(XX).LT.1.0D-9) GO TO 480
      XX=XX+XX
      VAL5=VAL1*SQRT2*2.0D+00
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
           DO 473 J=M4,M3
                LAD1=LAD+LADD(J+LOS)
                Z=C(JJ)*CII
                INT(LAD1+1)=INT(LAD1+1)+XX*Z
                INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
  473      JJ=JJ+1
  472 II=II+1
      GO TO 481
C
  480 VAL5=SQRT2*VAL1*2.0D+00
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
           CII=C(II)*VAL5
           DO 475 J=M4,M3
                LAD1=LAD+LADD(J+LOS)
                INT(LAD1+2)=INT(LAD1+2)+CII*C(JJ)
  475      JJ=JJ+1
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
           Z=C(II)*C(JJ)
           INT(LAD+1)=INT(LAD+1)+XX*Z
           INT(LAD+2)=INT(LAD+2)+VAL5*Z
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
           INT(LAD+2)=INT(LAD+2)+VAL5*C(II)*C(JJ)
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
           DO 477 J=M4,M3
                LAD1=LAD+LADD(J+LOS)
                Z=CJJ*C(II)
                INT(LAD1+3)=INT(LAD1+3)+XX*Z
                INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
  477      II=II+1
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
           CJJ=C(JJ)*VAL5
           DO 479 J=M4,M3
                LAD1=LAD+LADD(J+LOS)
                INT(LAD1+2)=INT(LAD1+2)+CJJ*C(II)
  479      II=II+1
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
      VAL4=VAL1*SQRT2*2.0D+00
      IF (IS.NE.0) GO TO 6587
      JMIN=MAX(2,JMING)
      JMAX=MIN(N,JMAXG)
      DO 1001 I=JMIN,JMAX
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
           CII=C(II)*VAL4
           DO 1002 J=M2,M1
                LAD1=LAD+LADD(J+AOS)
                INT(LAD1+1)=INT(LAD1+1)+C(JJ)*CII
 1002      JJ=JJ+1
 1001 CONTINUE
C
C   NO  28
C
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6599
      IF (M1.GT.N) M1=N
      IF (M2.LT.2) M2=2
      M1=MIN(M1,JMAXG)
      M2=MAX(M2,JMING)
      IF (M1.LT.M2) GO TO 6599
      JJ=JA+WTY(M2)
      DO 1003 I=M2,M1
           KSM=IEOR(ASM,SS(I))
           KOS=OS(KSM+1)
           LAD=IJADD(ARR+I)
           CJJ=C(JJ)*VAL4
           DO 1004 J=1,I-1
                II=IA+WAB(J)
                LSM=IEOR(KSM,SS(J))
                LOS=OS(LSM+1)
                LAD1=LAD+KADD(J+KOS)+LADD(J+LOS)
                INT(LAD1+1)=INT(LAD1+1)+CJJ*C(II)
C****************************************
      IF(LAD1.EQ.136)WRITE(ITAPE3,*)' CJJ,VAL4=',
     +                        CJJ*C(II),VAL4,INT(LAD1+1)
 1004      CONTINUE
 1003 JJ=JJ+1
C
C   NO 30
C
 6587 M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6599
      IF (M1.GT.N) M1=N
      IF (M2.LT.3) M2=3
      M1=MIN(M1,JMAXG)
      M2=MAX(M2,JMING)
      IF (M1.LT.M2) GO TO 6599
      JJ=JA+WTY(M2)
      DO 1005 I=M2,M1
           KSM=IEOR(ASM,SS(I))
           KOS=OS(KSM+1)
           LAD=IJADD(ARR+I)
           CJJ=C(JJ)*VAL1*2.0D+00
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
                     Z=CJJ*C(II)
                     INT(LAD2+1)=INT(LAD2+1)+Z
                     INT(LAD2+3)=INT(LAD2+3)+Z
 1007           II=II+1
 1006      CONTINUE
           JJ=JJ+1
 1005 CONTINUE
C
C   NO 32
C
clj The next statement branched to the inside of a do loop:
c6599 IF(N.LT.3) GO TO 1008
 6599 IF(N.LT.3) GO TO 1999
      JMIN=MAX(3,JMING)
      JMAX=MIN(N,JMAXG)
      DO 1008 I=JMIN,JMAX
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
                CII=C(II)*VAL1*2.0D+00
                DO 1010 K=M4,M3
                     LAD2=LAD1+LADD(K+LOS)
                     Z=C(JJ)*CII
                     INT(LAD2+1)=INT(LAD2+1)+Z
                     INT(LAD2+2)=INT(LAD2+2)+Z
 1010           JJ=JJ+1
 1009      II=II+1
 1008 CONTINUE
 1999 continue
C
C   NO 33
C
      SSI=IEOR(JS,IS)
      IF (SSI.GT.JS)GO TO  6608
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6608
      IF (M1.GT.N) M1=N
      M1=MIN(M1,JMAXG)
      M3=MX(SSI+1)
      M4=MN(SSI+1)
      IF (M4.GT.N) GO TO 6608
      IF (M2.LT.2) M2=2
      M2=MAX(M2,JMING)
      IF (M1.LT.M2) GO TO 6608
      JJ=JA+WTY(M2)
      DO 1011 I=M2,M1
           IF (SSI.EQ.JS) M3=I-1
           IF (M3.LT.M4) GO TO 1011
           II=IA+WTW(I,IS+1)+WTY(M4)
           KSM=IEOR(ASM,SS(I))
           KOS=OS(KSM+1)
           LAD=IJADD(ARR+I)+KADD(I+KOS)
           CJJ=C(JJ)*VAL1*2.0D+00
           DO 1012 J=M4,M3
                LAD1=LAD+LADD(J+AOS)
                Z=CJJ*C(II)
                INT(LAD1+2)=INT(LAD1+2)+Z
                INT(LAD1+1)=INT(LAD1+1)+Z
 1012      II=II+1
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
      M1=MIN(M1,JMAXG)
      M2=MAX(M2,JMING)
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
                Z=VAL1*C(II)*C(JJ)*2.0D+00
                INT(LAD1+1)=INT(LAD1+1)+Z
                INT(LAD1+2)=INT(LAD1+2)+Z
                II=II+1
 1014      JJ=JJ+1
 1013 CONTINUE
C
C    NO 38
C
 6609 JMIN=MAX(2,JMING)
      JMAX=MIN(N,JMAXG)
      DO 1015 I=JMIN,JMAX
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
                CJJ=C(JJ)*VAL1*2.0D+00
                DO 1017 K=M4,M3
                     LAD2=LAD1+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD2+3)=INT(LAD2+3)+Z
                     INT(LAD2+2)=INT(LAD2+2)+Z
 1017           II=II+1
 1016      JJ=JJ+1
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
           Z=VAL4*C(II)*C(JJ)*2.0D+00
C****************************************
      IF(LAD.EQ.136)WRITE(ITAPE3,*)' Z,VAL4=',Z,VAL4,INT(LAD+1)
           INT(LAD+TR1)=INT(LAD+TR1)+Z
           INT(LAD+TR2)=INT(LAD+TR2)+VAL2*Z
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
           CJJ=C(JJ)*VAL1*2.0D+00
           DO 1022 J=M4,M3
                LAD=ARR+LADD(J+AOS)
                Z=CJJ*C(II)
                INT(LAD+TR1)=INT(LAD+TR1)+Z
                INT(LAD+TR2)=INT(LAD+TR2)+VAL2*Z
 1022      II=II+1
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
           Z=0.0D+00
           DO 6014 J=M4,M3
                Z=Z+C(II)*C(JJ)
                II=II+1
 6014      JJ=JJ+1
           Z=Z*VAL1*2.0D+00
           INT(LAD+TR1)=INT(LAD+TR1)+Z
           INT(LAD+TR2)=INT(LAD+TR2)+VAL2*Z
 6013 CONTINUE
      RETURN
C
C    WY  NO 26   ISEG=15
C
 3101 VAL4=VAL1*SQRT2*2.0D+00
      IF (IS.NE.0) GO TO 6612
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6612
      IF (M1.GT.N) M1=N
      JJ=JA+WTY(M2)
      DO 1030 I=M2,M1
           LAD=ARR+LADD(I+AOS)
           Z=VAL4*C(IA+WAB(I))*C(JJ)
           INT(LAD+1)=INT(LAD+1)+Z
           INT(LAD+2)=INT(LAD+2)+Z
           INT(LAD+3)=INT(LAD+3)+Z
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
           CJJ=C(JJ)*VAL1*2.0D+00
           DO 1032 J=M4,M3
                LAD=ARR+LADD(J+AOS)
                Z=CJJ*C(II)
                INT(LAD+3)=INT(LAD+3)+Z
                INT(LAD+1)=INT(LAD+1)+Z
 1032      II=II+1
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
           Z=0.0D+00
           DO 7014 J=M4,M3
                Z=Z+C(II)*C(JJ)
                II=II+1
 7014      JJ=JJ+1
           Z=Z*VAL1
           INT(LAD+3)=INT(LAD+3)+Z+Z
           INT(LAD+1)=INT(LAD+1)+Z+Z
 7013 CONTINUE
      RETURN
C
C   WY  NO 26 ISEG=17
C
 3102 VAL4=VAL1*SQRT2*2.0D+00
      IF (IS.NE.0) GO TO 6614
      M1=MX(JS+1)
      M2=MN(JS+1)
      IF (M2.GT.N) GO TO 6614
      IF (M1.GT.N) M1=N
      JJ=JA+WTY(M2)
      DO 1023 I=M2,M1
           LAD=ARR+LADD(I+AOS)
           Z=VAL4*C(JJ)*C(IA+WAB(I))
           INT(LAD+3)=INT(LAD+3)+Z
           INT(LAD+2)=INT(LAD+2)+Z
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
           CJJ=C(JJ)*VAL1*2.0D+00
           DO 1025 J=M4,M3
                INT(ARR+LADD(J+AOS)+3)=INT(ARR+LADD(J+AOS)+3)+CJJ*C(II)
 1025      II=II+1
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
           Z=0.0D+00
           DO 5014 J=M4,M3
                Z=Z+C(II)*C(JJ)
                II=II+1
 5014      JJ=JJ+1
           INT(ARR+LADD(I+AOS)+3)=INT(ARR+LADD(I+AOS)+3)+VAL1*Z*2.0D+00
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
      CII=C(II)*VAL1*2.0D+00
      DO 601 I=M2,M1
           LAD=ARR+LADD(I+AOS)
           Z=C(JJ)*CII
           INT(LAD+TR1)=INT(LAD+TR1)+Z
           INT(LAD+TR2)=INT(LAD+TR2)+VAL2*Z
  601 JJ=JJ+1
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
           CII=C(II)*VAL1*2.0D+00
           DO 702 J=M4,M3
                LAD=ARR+LADD(J+AOS)
                Z=C(JJ)*CII
                INT(LAD+TR1)=INT(LAD+TR1)-Z
                INT(LAD+TR2)=INT(LAD+TR2)-Z*VAL2
  702      JJ=JJ+1
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
           Z=0.0D+00
           DO 9081 J=M4,M3
                Z=Z+C(II)*C(JJ)
                II=II+1
 9081      JJ=JJ+1
           Z=2.0D+00*VAL1*Z
           INT(LAD+TR1)=INT(LAD+TR1)+Z
           INT(LAD+TR2)=INT(LAD+TR2)+Z*VAL2
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
      CJJ=C(JJ)*VAL4*2.0D+00
      DO 901 I=1,N
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+LADD(I+LOS)
           INT(LAD+1)=INT(LAD+1)+CJJ*C(IA+WAB(I))
C****************************************
      IF(LAD.EQ.136)WRITE(ITAPE3,*)' CJJ,VAL4=',CJJ*C(IA+WAB(I)),VAL4
  901 CONTINUE
C
C     NO 25
C
 6631 IF (JS.NE.0) RETURN
      JJ=JA
      CJJ=C(JJ)*VAL1*2.0D+00
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
                Z=CJJ*C(II)
                INT(LAD1+1)=INT(LAD1+1)+Z
                INT(LAD1+3)=INT(LAD1+3)+Z
  903      II=II+1
  902 CONTINUE
      RETURN
C
C   YW   NO 19   ISEG=19
C***********************************************************************
  308 IF (ISEG.NE.19) GO TO 999
      VAL4=VAL1*SQRT2*2.0D+00
      IF (JS.NE.0) GO TO 6632
      M1=MX(IS+1)
      M2=MN(IS+1)
      IF (M2.GT.N) GO TO 6632
      IF (M1.GT.N) M1=N
      II=IA+WTY(M2)
      DO 801 I=M2,M1
           JJ=JA+WAB(I)
           LAD=ARR+LADD(I+AOS)
           Z=VAL4*C(II)*C(JJ)
           INT(LAD+TR1)=INT(LAD+TR1)+Z
           INT(LAD+TR2)=INT(LAD+TR2)+VAL2*Z
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
           CII=C(II)*VAL1*2.0D+00
           DO 803 J=M4,M3
                LAD=ARR+LADD(J+AOS)
                Z=CII*C(JJ)
                INT(LAD+TR1)=INT(LAD+TR1)+Z
                INT(LAD+TR2)=INT(LAD+TR2)+VAL2*Z
  803      JJ=JJ+1
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
           Z=0.0D+00
           DO 9083 J=M4,M3
                Z=Z+C(II)*C(JJ)
                II=II+1
 9083      JJ=JJ+1
           Z=2.0D+00*VAL1*Z
           INT(LAD+TR1)=INT(LAD+TR1)+Z
           INT(LAD+TR2)=INT(LAD+TR2)+VAL2*Z
 9082 CONTINUE
      RETURN
C
C   WZ   NO 23 ISEG=10
C
 3091 IF (IS.NE.0.OR.JS.NE.0) GO TO 6634
      JJ=JA
      CJJ=C(JJ)*VAL1*2.0D+00
      DO 904 I=1,N
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+LADD(I+LOS)+1
           INT(LAD)=INT(LAD)+CJJ*C(IA+WAB(I))
  904 CONTINUE
C
C      NO  25
C
 6634 VAL4=SQRT2*VAL1
      IF (JS.NE.0) RETURN
      JJ=JA
      CJJ=C(JJ)*VAL4*2.0D+00
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
           LAD=ARR+KADD(I+AOS)+1
           DO 906 J=M2,M1
                INT(LAD+LADD(J+LOS))=INT(LAD+LADD(J+LOS))+CJJ*C(II)
  906      II=II+1
  905 CONTINUE
      RETURN
C
C    YY    NO  9    FOUR EXTERNALS
C***********************************************************************
C
C?K   WRITE(ITAPE3,*) ' PRINTING ISEG AND M IN INITEX '
  999 WRITE(ITAPE6,998) M,ISEG
  998 FORMAT (/' PRINTING M AND ISEG IN SHAPES',2I4)
      RETURN
C
C     CLOSED INTERNAL LOOP
C***********************************************************************
  318 CONTINUE
      Z=0.0D+00
      II=IA
      JJ=JA
      DO 4000 I=1,N1
           Z=Z+C(II)*C(JJ)
           II=II+1
 4000 JJ=JJ+1
      INT(ARR+TR1)=INT(ARR+TR1)+VAL1*Z*2.0D+00
      INT(ARR+TR2)=INT(ARR+TR2)+VAL1*VAL2*Z*2.0D+00
      RETURN
C
C***********************************************************************
  314 PRINT 315
  315 FORMAT(' INCORRECT ENTRY TO THIS SUBROUTINE. THIS IS FOR 4X.')
      STOP
C
C     ----- INTERNAL CASE FOR TRACKS OF (3,2,1) -----
C***********************************************************************
  319 CONTINUE
      Z=0.0D+00
      II=IA
      JJ=JA
      DO 3999 I=1,N1
           Z=Z+C(II)*C(JJ)
           II=II+1
 3999 JJ=JJ+1
      INT(ARR+1)=INT(ARR+1)+VAL1*Z*2.0D+00
      INT(ARR+2)=INT(ARR+2)+VAL1*Z*2.0D+00
      INT(ARR+3)=INT(ARR+3)+VAL1*Z*2.0D+00
      RETURN
C
C***********************************************************************
C
  317 Z=VAL1*(INT(ARR+TR1)+VAL2*INT(ARR+TR2))
C   ZZ   ARRIVAL,FOR COMPLETENESS
      STOP ' ZZ ARRIVAL IN EXTERNAL'
C     RETURN
      END
