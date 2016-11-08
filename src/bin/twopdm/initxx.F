      SUBROUTINE INITXX(INT,C,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
      IMPLICIT REAL*8 (A-H,O-Z)
C     EXTENDED DUMMY INT,C
      INTEGER ARR,TR1,TR2,ASM,AOS,OS,WTW,WTX,WTY,WAB,SS,SSI,SSJ,SYMORB
      INTEGER BMAX,ORBFRM
      integer nmax,nwksmx,numij,norbs,nsym
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
      ENTRY XXDM(INT,C,IJADD,KADD,LADD,WTW,WTX,WTY,WAB,SS)
      AOS=OS(ASM+1)
      GO TO 305
C
C    XX   NO 11  ISEG=8
C*********************************************************************
 3054 VAL4=VAL1+VAL1
C*********************************************************************
C     WRITE(ITAPE3,*) ' IN XXDM, IS,JS,IA,JA= ',IS,JS,IA,JA
C*********************************************************************
      IF (IS.NE.JS) GO TO 6629
      IF(IA.EQ.JA) GO TO 8505
      DO 510 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 510
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 510
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M1.LT.M2) GO TO 510
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
           Z1=0.0D+00
           JJ=JA+WTX(I,IS+1)+WTY(M2)
           II=IA+WTX(I,IS+1)+WTY(M2)
           DO 511 J=M2,M1
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
  511      II=II+1
           INT(LAD1+1)=INT(LAD1+1)+Z1
  510 CONTINUE
C
C    NO  12
C
 6629 DO 512 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SSI) GO TO 512
           IF (SSI.GT.SS(I)) GO TO 512
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 512
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 512
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 512
           II=IA+WTX(I,IS+1)+WTY(M2)
           JJA=JA+WTX(I,JS+1)
           DO 513 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 513
                JJ=JJA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+1
                CII=C(II)*VAL4
                DO 514 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+C(JJ)*CII
  514           JJ=JJ+1
  513      II=II+1
  512 CONTINUE
C
C    NO 15
C
      DO 518 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(JS,SSI)
           IF (SSJ.GT.SSI) GO TO 518
           IF (SSI.GT.SS(I)) GO TO 518
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 518
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 518
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 518
           II=IA+WTX(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+1
           DO 519 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 519
                JJ=JA+WTX(J,JS+1)+WTY(M4)
                CII=C(II)*VAL4
                DO 520 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))-C(JJ)*CII
  520           JJ=JJ+1
  519      II=II+1
  518 CONTINUE
C
C   NO 13
C
 8505 DO 515 I=3,N
C***********************************************************************
C     WRITE(ITAPE3,*) ' HERE WE ARE AT 8505 WITH I = ',I
C***********************************************************************
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(SS(I),IS)
           IF (SSJ.GT.SS(I)) GO TO 515
           IF (SSI.GT.SSJ) GO TO 515
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 515
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 515
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 515
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           IIA=IA+WTX(I,IS+1)
           DO 516 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 516
                II=IIA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+TR1
                CJJ=C(JJ)*VAL4
                DO 517 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+CJJ*C(II)
  517           II=II+1
  516      JJ=JJ+1
  515 CONTINUE
C
C    NO 16
C
      DO 521 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SSJ) GO TO 521
           IF (SSJ.GT.SS(I)) GO TO 521
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 521
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 521
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 521
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+TR1
           DO 522 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 522
                II=IA+WTX(J,IS+1)+WTY(M4)
                CJJ=C(JJ)*VAL4
                DO 523 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))-C(II)*CJJ
  523           II=II+1
  522      JJ=JJ+1
  521 CONTINUE
C
C NO 120
C
      DO 9030 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSJ.GT.SSI) GO TO 9030
           IF (SSI.GT.SS(I)) GO TO 9030
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 9030
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 9030
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 9030
           JJA=JA+WTX(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+TR1
           DO 9031 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9031
                JJ=JJA+WTY(M4)
                II=IA+WTX(J,IS+1)+WTY(M4)
                Z=0.0D+00
                DO 9032 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9032           JJ=JJ+1
                INT(LAD+LADD(J+LOS))=INT(LAD+LADD(J+LOS))+VAL4*Z
 9031      CONTINUE
 9030 CONTINUE
C
C  NO 121
C
C**********************************************************************
C     WRITE(ITAPE3,*) ' IN XXDM, IA,JA = ',IA,JA
C**********************************************************************
      IF (IA.EQ.JA) RETURN
      DO 9050 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SSI,JS)
           IF (SSI.GT.SSJ) GO TO 9050
           IF (SSJ.GT.SS(I)) GO TO 9050
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 9050
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 9050
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GOTO 9050
           IIA=IA+WTX(I,IS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+1
           DO 9051 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9051
                II=IIA+WTY(M4)
                JJ=JA+WTX(J,JS+1)+WTY(M4)
                Z=0.0D+00
                DO 9052 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9052           JJ=JJ+1
                INT(LAD+LADD(J+LOS))=INT(LAD+LADD(J+LOS))+VAL4*Z
 9051      CONTINUE
 9050 CONTINUE
      RETURN
C
C    XX   NO 11   ISEG=5
C
 3052 XX=-VAL1/SQRT2+VAL3
      IF (DABS(XX).LT.1.0D-9) GO TO 526
      XX=XX+XX
      VAL5=SQRT2*VAL1*2.0D+00
      IF(IA.EQ.JA) GO TO 524
      IF (IS.NE.JS) GO TO 524
      DO 529 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 529
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 529
           IF (SSJ.EQ.SS(I)) M1=I-1
           II=IA+WTX(I,IS+1)+WTY(M2)
           JJ=JA+WTX(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
           Z1=0.0D+00
           DO 530 J=M2,M1
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
                Z=C(II)*C(JJ)
C****************************************
      IF(LAD.EQ.136)WRITE(ITAPE3,*)' Z,VAL4=',Z,VAL4
                Z1=Z1+Z
C****************************************
      IF(LAD1.EQ.136)WRITE(ITAPE3,*)' Z1,VAL4=',Z1,VAL4
                INT(LAD+1)=INT(LAD+1)+Z*XX
                INT(LAD+2)=INT(LAD+2)+VAL5*Z
                JJ=JJ+1
  530      II=II+1
           INT(LAD1+1)=INT(LAD1+1)+Z1*XX
           INT(LAD1+2)=INT(LAD1+2)+VAL5*Z1
  529 CONTINUE
      GO TO 524
C
  526 VAL5=SQRT2*VAL1*2.0D+00
      IF(IA.EQ.JA) GO TO 527
      IF (IS.NE.JS) GO TO 527
      DO 531 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 531
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 531
           IF (SSJ.EQ.SS(I)) M1=I-1
           II=IA+WTX(I,IS+1)+WTY(M2)
           JJ=JA+WTX(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
           Z1=0.0D+00
           DO 532 J=M2,M1
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
                Z=VAL5*C(II)*C(JJ)
C****************************************
      IF(LAD.EQ.136)WRITE(ITAPE3,*)' Z,VAL4=',Z,VAL4
                Z1=Z1+Z
C****************************************
      IF(LAD1.EQ.136)WRITE(ITAPE3,*)' Z1,VAL4=',Z1,VAL4
                INT(LAD+2)=INT(LAD+2)+Z
                JJ=JJ+1
  532      II=II+1
           INT(LAD1+2)=INT(LAD1+2)+Z1
  531 CONTINUE
      GO TO 527
C
C      NO  16
C
  524 DO 533 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SSJ) GO TO 533
           IF (SSJ.GT.SS(I)) GO TO 533
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 533
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 533
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 533
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=KADD(I+AOS)+ARR
           DO 534 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 534
                II=IA+WTX(J,IS+1)+WTY(M4)
                CJJ=C(JJ)
                DO 535 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=C(II)*CJJ
                     INT(LAD1+1)=INT(LAD1+1)-Z*XX
                     INT(LAD1+2)=INT(LAD1+2)-VAL5*Z
  535           II=II+1
  534      JJ=JJ+1
  533 CONTINUE
      GO TO 525
C
  527 DO 536 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SSJ) GO TO 536
           IF (SSJ.GT.SS(I)) GO TO 536
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 536
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 536
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 536
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+2
           DO 537 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 537
                II=IA+WTX(J,IS+1)+WTY(M4)
                CJJ=C(JJ)*VAL5
                DO 538 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))-CJJ*C(II)
  538           II=II+1
  537      JJ=JJ+1
  536 CONTINUE
      GO TO 528
C
C    NO 13
C
  525 DO 539 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(SS(I),IS)
           IF (SSJ.GT.SS(I)) GO TO 539
           IF (SSI.GT.SSJ) GO TO 539
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 539
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 539
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 539
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           IIA=IA+WTX(I,IS+1)
           DO 540 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 540
                II=IIA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)
                CJJ=C(JJ)
                DO 541 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD1+1)=INT(LAD1+1)+Z*XX
                     INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
  541           II=II+1
  540      JJ=JJ+1
  539 CONTINUE
C
C NO 120
C
      DO 9010 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSJ.GT.SSI) GO TO 9010
           IF (SSI.GT.SS(I)) GO TO 9010
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 9010
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 9010
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 9010
           JJA=JA+WTX(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 9011 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9011
                LAD1=LAD+LADD(J+LOS)
                Z=0.0D+00
                JJ=JJA+WTY(M4)
                II=IA+WTX(J,IS+1)+WTY(M4)
                DO 9012 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9012           JJ=JJ+1
                INT(LAD1+1)=INT(LAD1+1)+Z*XX
                INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
 9011      CONTINUE
 9010 CONTINUE
      RETURN
C
  528 DO 542 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(SS(I),IS)
           IF (SSJ.GT.SS(I)) GO TO 542
           IF (SSI.GT.SSJ) GO TO 542
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 542
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 542
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 542
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           IIA=IA+WTX(I,IS+1)
           DO 543 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 543
                II=IIA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+2
                CJJ=C(JJ)*VAL5
                DO 544 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+CJJ*C(II)
  544           II=II+1
  543      JJ=JJ+1
  542 CONTINUE
C
C NO 120
C
      DO 9013 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSJ.GT.SSI) GO TO 9013
           IF (SSI.GT.SS(I)) GO TO 9013
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 9013
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 9013
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 9013
           JJA=JA+WTX(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+2
           DO 9014 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9014
                II=IA+WTX(J,IS+1)+WTY(M4)
                JJ=JJA+WTY(M4)
                Z=0.0D+00
                DO 9015 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9015           JJ=JJ+1
                INT(LAD+LADD(J+LOS))=INT(LAD+LADD(J+LOS))+VAL5*Z
 9014      CONTINUE
 9013 CONTINUE
      RETURN
C
C   XX   NO  11  ISEG=6
C
  305 IF (ISEG.EQ.4) GO TO 3051
      IF (ISEG.EQ.5) GO TO 3052
      IF (ISEG.EQ.8) GO TO 3054
      IF (ISEG.NE.6) GO TO 999
      XX=-VAL1/SQRT2+VAL3
      IF (DABS(XX).LT.1.0D-9) GO TO 550
      XX=XX+XX
C*********************************************************************
C IS THIS A BUG ???????????????????? TJL  11/08/84
C     VAL5=SQRTS*VAL1*2.0D+00
C*********************************************************************
      VAL5=SQRT2*VAL1*2.0D+00
      IF(IA.EQ.JA) GO TO 551
      IF (IS.NE.JS) GO TO 551
      DO 560 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 560
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 560
           IF (SSJ.EQ.SS(I)) M1=I-1
           II=IA+WTX(I,IS+1)+WTY(M2)
           JJ=JA+WTX(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD1=ARR+LADD(I+LOS)+KADD(I+AOS)
           Z1=0.0D+00
           DO 561 J=M2,M1
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
                Z=C(II)*C(JJ)
C******************  PROBLEM AREA ***** STAY AWAY ********************
C****************************************
      IF(LAD.EQ.136)WRITE(ITAPE3,*)' 561 LOOP,Z,XX=',Z,XX
                Z1=Z1+Z
C****************************************
      IF(LAD1.EQ.136)WRITE(ITAPE3,*)' 561 LOOP,Z1,XX=',Z1,XX
                INT(LAD+1)=INT(LAD+1)+Z*XX
                INT(LAD+2)=INT(LAD+2)+VAL5*Z
                JJ=JJ+1
  561      II=II+1
           INT(LAD1+1)=INT(LAD1+1)+Z1*XX
           INT(LAD1+2)=INT(LAD1+2)+VAL5*Z1
  560 CONTINUE
      GO TO 551
C
  550 VAL5=SQRT2*VAL1*2.0D+00
      IF(IA.EQ.JA) GO TO 555
      IF (IS.NE.JS ) GO TO 555
      DO 562 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 562
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 562
           IF (SSJ.EQ.SS(I)) M1=I-1
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
           Z1=0.0D+00
           JJ=JA+WTX(I,IS+1)+WTY(M2)
           II=IA+WTX(I,IS+1)+WTY(M2)
           DO 563 J=M2,M1
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+LADD(J+LOS)+2
                Z=VAL5*C(II)*C(JJ)
                Z1=Z1+Z
                INT(LAD)=INT(LAD)+Z
                JJ=JJ+1
  563      II=II+1
           INT(LAD1+2)=INT(LAD1+2)+Z1
  562 CONTINUE
      GO TO 555
C
C      NO 12
C
  551 DO 564 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SSI) GO TO 564
           IF (SSI.GT.SS(I)) GO TO 564
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 564
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 564
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 564
           II=IA+WTX(I,IS+1)+WTY(M2)
           JJA=JA+WTX(I,JS+1)
           DO 565 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GOTO 565
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                JJ=JJA+WTY(M4)
                LAD=ARR+KADD(J+AOS)
                CII=C(II)
                DO 566 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=C(JJ)*CII
                     INT(LAD1+1)=INT(LAD1+1)+Z*XX
                     INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
  566           JJ=JJ+1
  565      II=II+1
  564 CONTINUE
      GO TO 552
C
  555 DO 567 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SSI) GO TO 567
           IF (SSI.GT.SS(I)) GO TO 567
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 567
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 567
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 567
           II=IA+WTX(I,IS+1)+WTY(M2)
           JJA=JA+WTX(I,JS+1)
           DO 568 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 568
                JJ=JJA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+2
                CII=C(II)*VAL5
                DO 569 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+CII*C(JJ)
  569           JJ=JJ+1
  568      II=II+1
  567 CONTINUE
      GO TO 556
C
C    NO  13
C
  552 DO 570 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(SS(I),IS)
           IF (SSJ.GT.SS(I)) GO TO 570
           IF (SSI.GT.SSJ) GO TO 570
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 570
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 570
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 570
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           IIA=IA+WTX(I,IS+1)
           DO 571 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GOTO 571
                II=IIA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)
                CJJ=C(JJ)
                DO 572 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD1+3)=INT(LAD1+3)+Z*XX
                     INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
  572           II=II+1
  571      JJ=JJ+1
  570 CONTINUE
      GO TO 553
C
  556 DO 573 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(SS(I),IS)
           IF (SSJ.GT.SS(I)) GO TO 573
           IF (SSI.GT.SSJ) GO TO 573
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 573
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 573
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 573
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           IIA=IA+WTX(I,IS+1)
           DO 574 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 574
                II=IIA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+2
                CJJ=C(JJ)*VAL5
                DO 575 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))+CJJ*C(II)
  575           II=II+1
  574      JJ=JJ+1
  573 CONTINUE
      GO TO 557
C
C    NO  15
C
  553 DO 576 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(JS,SSI)
           IF (SSJ.GT.SSI) GO TO 576
           IF (SSI.GT.SS(I)) GO TO 576
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 576
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 576
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 576
           II=IA+WTX(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 577 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 577
                JJ=JA+WTX(J,JS+1)+WTY(M4)
                CII=C(II)
                DO 578 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=CII*C(JJ)
                     INT(LAD1+1)=INT(LAD1+1)-Z*XX
                     INT(LAD1+2)=INT(LAD1+2)-VAL5*Z
  578           JJ=JJ+1
  577      II=II+1
  576 CONTINUE
      GO TO 554
C
  557 DO 579 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(JS,SSI)
           IF (SSJ.GT.SSI) GO TO 579
           IF (SSI.GT.SS(I)) GO TO 579
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 579
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 579
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 579
           II=IA+WTX(I,IS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)+2
           DO 580 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 580
                JJ=JA+WTX(J,JS+1)+WTY(M4)
                CII=C(II)*VAL5
                DO 581 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))-CII*C(JJ)
  581           JJ=JJ+1
  580      II=II+1
  579 CONTINUE
      GO TO 558
C
C    NO   16
C
  554 DO 582 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SSJ) GO TO 582
           IF (SSJ.GT.SS(I)) GO TO 582
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 582
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 582
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 582
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 583 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 583
                II=IA+WTX(J,IS+1)+WTY(M4)
                CJJ=C(JJ)
                DO 584 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD1+3)=INT(LAD1+3)-Z*XX
                     INT(LAD1+2)=INT(LAD1+2)-VAL5*Z
  584           II=II+1
  583      JJ=JJ+1
  582 CONTINUE
C
C NO 120
C
      DO 9020 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSJ.GT.SSI) GO TO 9020
           IF (SSI.GT.SS(I)) GO TO 9020
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 9020
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 9020
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 9020
           JJA=JA+WTX(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 9021 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9021
                LAD1=LAD+LADD(J+LOS)
                Z=0.0D+00
                JJ=JJA+WTY(M4)
                II=IA+WTX(J,IS+1)+WTY(M4)
                DO 9022 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9022           JJ=JJ+1
                INT(LAD1+3)=INT(LAD1+3)+Z*XX
                INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
 9021      CONTINUE
 9020 CONTINUE
      IF (IA.EQ.JA) RETURN
C
C NO 121
C
      DO 9040 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SSI,JS)
           IF (SSI.GT.SSJ) GO TO 9040
           IF (SSJ.GT.SS(I)) GO TO 9040
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 9040
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 9040
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 9040
           IIA=IA+WTX(I,IS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 9041 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9041
                LAD1=LAD+LADD(J+LOS)
                Z=0.0D+00
                JJ=JA+WTX(J,JS+1)+WTY(M4)
                II=IIA+WTY(M4)
                DO 9042 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9042           JJ=JJ+1
                INT(LAD1+1)=INT(LAD1+1)+Z*XX
                INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
 9041      CONTINUE
 9040 CONTINUE
      RETURN
C
  558 DO 585 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SSJ) GO TO 585
           IF (SSJ.GT.SS(I)) GO TO 585
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 585
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 585
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 585
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=KADD(I+AOS)+ARR+2
           DO 586 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 586
                II=IA+WTX(J,IS+1)+WTY(M4)
                CJJ=C(JJ)*VAL5
                DO 587 K=M4,M3
                     INT(LAD+LADD(K+LOS))=INT(LAD+LADD(K+LOS))-CJJ*C(II)
  587           II=II+1
  586      JJ=JJ+1
  585 CONTINUE
C
C NO 120
C
      DO 9023 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSJ.GT.SSI) GO TO 9023
           IF (SSI.GT.SS(I)) GO TO 9023
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 9023
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 9023
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 9023
           JJA=JA+WTX(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 9024 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9024
                JJ=JJA+WTY(M4)
                II=IA+WTX(J,IS+1)+WTY(M4)
                LAD1=LAD+LADD(J+LOS)
                Z=0.0D+00
                DO 9025 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9025           JJ=JJ+1
                INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
 9024      CONTINUE
 9023 CONTINUE
C
C NO 121
C
      IF (IA.EQ.JA) RETURN
      DO 9043 I=3,N
           SSI=IEOR(SS(I),IS)
           SSJ=IEOR(SSI,JS)
           IF (SSI.GT.SSJ) GO TO 9043
           IF (SSJ.GT.SS(I)) GO TO 9043
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 9043
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 9043
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           IF (M1.LT.M2) GO TO 9043
           IIA=IA+WTX(I,IS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 9044 J=M2,M1
                LAD1=LAD+LADD(J+LOS)
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9044
                JJ=JA+WTX(J,JS+1)+WTY(M4)
                II=IIA+WTY(M4)
                Z=0.0D+00
                DO 9045 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9045           JJ=JJ+1
                INT(LAD1+2)=INT(LAD1+2)+VAL5*Z
 9044      CONTINUE
 9043 CONTINUE
      RETURN
C
C   XX   NO  11   ISEG=4
C
 3051 CONTINUE
      VAL4=-2.0D+00*VAL1/SQRT2
      IF(IA.EQ.JA) GO TO 8506
      IF (IS.NE.JS) GO TO 8506
      DO 501 I=2,N
           SSJ=IEOR(SS(I),JS)
           IF (SSJ.GT.SS(I)) GO TO 501
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 501
           IF (SSJ.EQ.SS(I)) M1=I-1
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD1=ARR+KADD(I+AOS)+LADD(I+LOS)
           Z1=0.0D+00
           II=IA+WTX(I,IS+1)+WTY(M2)
           JJ=JA+WTX(I,IS+1)+WTY(M2)
           DO 502 J=M2,M1
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)+LADD(J+LOS)
                Z=VAL4*C(II)*C(JJ)
                Z1=Z1+Z
                INT(LAD+1)=INT(LAD+1)+Z
                INT(LAD+2)=INT(LAD+2)-Z-Z
                JJ=JJ+1
  502      II=II+1
           INT(LAD1+1)=INT(LAD1+1)+Z1
           INT(LAD1+2)=INT(LAD1+2)-Z1-Z1
  501 CONTINUE
C
C    NO  13
C
 8506 DO 503 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(SS(I),IS)
           IF (SSJ.GT.SS(I)) GO TO 503
           IF (SSI.GT.SSJ) GO TO 503
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 503
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 503
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           IIA=IA+WTX(I,IS+1)
           DO 504 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 504
                II=IIA+WTY(M4)
                LSM=IEOR(ASM,SS(J))
                LOS=OS(LSM+1)
                LAD=ARR+KADD(J+AOS)
                CJJ=C(JJ)*VAL4
                DO 505 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=CJJ*C(II)
                     INT(LAD1+1)=INT(LAD1+1)+Z
                     INT(LAD1+2)=INT(LAD1+2)-Z-Z
  505           II=II+1
  504      JJ=JJ+1
  503 CONTINUE
C
C   NO  16
C
      DO 506 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSI.GT.SSJ) GO TO 506
           IF (SSJ.GT.SS(I)) GO TO 506
           M1=MX(SSJ+1)
           M2=MN(SSJ+1)
           IF (M2.GT.N) GO TO 506
           IF (SSJ.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           M3=MX(SSI+1)
           M4=MN(SSI+1)
           IF (M4.GT.N) GO TO 506
           JJ=JA+WTX(I,JS+1)+WTY(M2)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=KADD(I+AOS)+ARR
           DO 507 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 507
                II=IA+WTX(J,IS+1)+WTY(M4)
                CJJ=C(JJ)*VAL4
                DO 508 K=M4,M3
                     LAD1=LAD+LADD(K+LOS)
                     Z=C(II)*CJJ
                     INT(LAD1+1)=INT(LAD1+1)-Z
                     INT(LAD1+2)=INT(LAD1+2)+Z+Z
  508           II=II+1
  507      JJ=JJ+1
  506 CONTINUE
C
C NO 120
C
      DO 9034 I=3,N
           SSJ=IEOR(SS(I),JS)
           SSI=IEOR(IS,SSJ)
           IF (SSJ.GT.SSI) GO TO 9034
           IF (SSI.GT.SS(I)) GO TO 9034
           M1=MX(SSI+1)
           M2=MN(SSI+1)
           IF (M2.GT.N) GO TO 9034
           M3=MX(SSJ+1)
           M4=MN(SSJ+1)
           IF (M4.GT.N) GO TO 9034
           IF (SSI.EQ.SS(I)) M1=I-1
           IF (M2.LT.2) M2=2
           JJA=JA+WTX(I,JS+1)
           LSM=IEOR(ASM,SS(I))
           LOS=OS(LSM+1)
           LAD=ARR+KADD(I+AOS)
           DO 9035 J=M2,M1
                IF (SSI.EQ.SSJ) M3=J-1
                IF (M3.LT.M4) GO TO 9035
                LAD1=LAD+LADD(J+LOS)
                Z=0.0D+00
                JJ=JJA+WTY(M4)
                II=IA+WTX(J,IS+1)+WTY(M4)
                DO 9036 K=M4,M3
                     Z=Z+C(II)*C(JJ)
                     II=II+1
 9036           JJ=JJ+1
                Z=Z*VAL4
                INT(LAD1+1)=INT(LAD1+1)+Z
                INT(LAD1+2)=INT(LAD1+2)-Z-Z
 9035      CONTINUE
 9034 CONTINUE
  999 RETURN
      END
