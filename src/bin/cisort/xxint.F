      SUBROUTINE XXINT(INT,H,ORBSYM,IJGRP,IJADD,KADD,LADD,IINT,
     *ACRCY,ILBL,ILBM,IJXX,KLXX,IJWW,KLWW,iad,val,ival)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
CTJL  EXTENDED DUMMY INT,H,ORBSYM,IJGRP,IJADD,KADD,LADD,IINT
C
      integer nmax,nijvir,numij,norbs
      INTEGER ARR,SYMORB,IINT(1)
      INTEGER MAXB,LVFRM1,ORBSYM
      integer iad(1),ival(1)
      REAL*8 INT(NMAX),H(NIJVIR)
      real*8 val(1)
C
      DIMENSION KADD(SYMORB),LADD(SYMORB),IJADD(NUMIJ),ORBSYM(NORBS)
      DIMENSION IJGRP(NUMIJ)
C
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,NROW4O,NWKSOC,NLEVOC,NORBOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,MAXB,NROOTS,LVFRM1,NREFS
      COMMON /INTS/ NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP12,ITAP99,ITAP04
     *,             ITAPE3,ITAP05,ITAP06
      COMMON /SORT/   IBLOCK,LNBUF,MAX
      COMMON /XX4X/   NIJVIR
C
      SQRT2=SQRT(2.0D+00)
C
C     WRITE(ITAPE6,971) (INT(IQ),IQ=1,NMAX)
  971 FORMAT (//,(1X,5(1PE15.5)))
      DO 1 I=1,NORBS
        ORBSYM(I)=ORBSYM(I)-1
    1 CONTINUE
C
      IBLOCK=1
C
C     WRITE(ITAPE6,"(3E18.8)") INT
C
      DO 6 I=LVFRM1,2,-1
        IA=I*(I-1)/2
        II=IA+I
    2   IF (IJGRP(II).EQ.IBLOCK) GO TO 3
        IBLOCK=IBLOCK+1
        IF (IBLOCK.GT.NGROUP) THEN
          WRITE(ITAPE6,*) ' TOO MANY INT BLOCKS IN HXX,H'
          CALL MABORT
        END IF
        CALL NXTBLK(ACRCY,ILBL,ILBM,IINT)
        GO TO 2
    3   CONTINUE
        III=IJADD(IA+I)+KADD(I)
        IIIS=ORBSYM(I)
        IIIN=IIIS*NORBS
        DO 5 J=1,I-1
          IF (ORBSYM(J).NE.IIIS) GO TO 4
          H(IA+J)=INT(III+LADD(IIIN+J)+3)
    4     CONTINUE
    5   CONTINUE
    6 CONTINUE
C
C
C
      CALL SREW(ITAP52)
      CALL RSETSA(ITAP52,INTSRT)
      CALL SREAD(ITAP52,IINT,NMAX2)
C
C
      INTPT=INTSRT
      IBLOCK=1
C
      DO 110 I=NORBS,LEVFRM,-1
        IA=I*(I-1)/2
        IS=ORBSYM(I)
        DO 109 J=LVFRM1,2,-1
          IJ=IA+J
  101     IF (IJGRP(IJ).EQ.IBLOCK) GO TO 102
          IBLOCK=IBLOCK+1
          IF (IBLOCK.GT.NGROUP) THEN
            WRITE(ITAPE6,*) ' TOO MANY BLOCKS FOR 3X'
            CALL MABORT
          END IF
          CALL RWRIT(ITAP52,IINT,INTOWP(NMAX),INTPT)
          INTPT=INTPT+I2SEC(INTOWP(NMAX))
          CALL RREAD(ITAP52,IINT,INTOWP(NMAX),INTPT)
C     CALL NXTBLK(ACRCY,ILBL,ILBM,IINT)
          GO TO 101
  102     CONTINUE
          JA=J*(J-1)/2
          IJS=IEOR(IS,ORBSYM(J))
          IJN=IJS*NORBS
          IJP=IJADD(IJ)
C
          JJ=JA+J
          IJJ=IJP+KADD(IJN+J)
          IJJN=IS*NORBS
          DO 104 L=1,J-1
            IF (ORBSYM(L).NE.IS) GO TO 103
            IJJL=IJJ+LADD(IJJN+L)
            IL=IA+L
            JL=JA+L
            T1=INT(IJJL+1)
            T2=INT(IJJL+2)
            CALL OUTWY(SQRT2*T1,IL,JJ,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,
     *      IINT,iad,val,ival)
            CALL OUTWY(T1+T2,IJ,JL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
            CALL OUTXY(T1-T2,IJ,JL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
  103       CONTINUE
  104     CONTINUE
          DO 108 K=1,J-1
            IJK=IJP+KADD(IJN+K)
            IK=IA+K
            JK=JA+K
            KA=K*(K-1)/2
            IJKS=IEOR(IJS,ORBSYM(K))
            IJKN=IJKS*NORBS
            IF (ORBSYM(K).NE.IJKS) GO TO 105
            IJKK=IJK+LADD(IJKN+K)
            KK=KA+K
            T1=INT(IJKK+1)
            T2=INT(IJKK+2)
            CALL OUTWY(SQRT2*T1,IJ,KK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,
     *      IINT,
     #        iad,val,ival)
            CALL OUTWY(T2+T1,IK,JK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
            CALL OUTXY(T2-T1,IK,JK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
  105       CONTINUE
            DO 107 L=1,K-1
              IF (ORBSYM(L).NE.IJKS) GO TO 106
              IJKL=IJK+LADD(IJKN+L)
              IL=IA+L
              JL=JA+L
              KL=KA+L
              T1=INT(IJKL+1)
              T2=INT(IJKL+2)
              T3=INT(IJKL+3)
              CALL OUTWY(T1+T3,IJ,KL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
              CALL OUTXY(T1-T3,IJ,KL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
              CALL OUTWY(T2+T3,IK,JL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
              CALL OUTXY(T2-T3,IK,JL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
              CALL OUTWY(T2+T1,IL,JK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
              CALL OUTXY(T2-T1,IL,JK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
  106         CONTINUE
  107       CONTINUE
  108     CONTINUE
  109   CONTINUE
  110 CONTINUE
C
C     ----- 4-EXTERNAL INTEGRALS -----
C
      DO 23 I=LVFRM1,2,-1
        IA=I*(I-1)/2
        IS=ORBSYM(I)
        II=IA+I
    7   IF (IJGRP(II).EQ.IBLOCK) GO TO 8
        IBLOCK=IBLOCK+1
        IF (IBLOCK.GT.NGROUP) THEN
          WRITE(ITAPE6,*) ' TOO MANY INT BLOCKS IN HXX,H'
          CALL MABORT
        END IF
        CALL RWRIT(ITAP52,IINT,INTOWP(NMAX),INTPT)
        INTPT=INTPT+I2SEC(INTOWP(NMAX))
        CALL RREAD(ITAP52,IINT,INTOWP(NMAX),INTPT)
C     CALL NXTBLK(ACRCY,ILBL,ILBM,IINT)
        GO TO 7
    8   CONTINUE
        IIP=IJADD(II)
        III=IIP+KADD(I)
        IIIN=IS*NORBS
        DO 10 L=1,I-1
          IF (ORBSYM(L).NE.IS) GO TO 9
          IL=IA+L
          IIIL=III+LADD(IIIN+L)
          LL=L*(L+1)/2
          T=SQRT2*(INT(IIIL+1)+INT(IIIL+3))
          CALL OUTWW(T,II,IL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
          T=SQRT2*(INT(IIIL+2)+INT(IIIL+3))
          CALL OUTWW(T,IL,LL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
    9     CONTINUE
   10   CONTINUE
        DO 13 K=1,I-1
          IIK=IIP+KADD(K)
          IIKS=ORBSYM(K)
          IIKN=IIKS*NORBS
          KA=K*(K-1)/2
          IK=IA+K
          KK=KA+K
          IIKK=IIK+LADD(IIKN+K)
          T=INT(IIKK+1)
          CALL OUTWW(T,II,KK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
          DO 12 L=1,K-1
            IF (ORBSYM(L).NE.IIKS) GO TO 11
            IIKL=IIK+LADD(IIKN+L)
            IL=IA+L
            KL=KA+L
            T1=INT(IIKL+1)
            T2=INT(IIKL+2)
            HKL=H(KL)
            T=T2-T1+HKL
            CALL OUTXX(T,IK,IL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
            T=SQRT2*T1
            CALL OUTWW(T,II,KL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
            T=T1+T2+HKL
            CALL OUTWW(T,IK,IL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
   11       CONTINUE
   12     CONTINUE
   13   CONTINUE
C
        DO 22 J=I-1,2,-1
          IJ=IA+J
   14     IF (IJGRP(IJ).EQ.IBLOCK) GO TO 15
          IBLOCK=IBLOCK+1
          IF (IBLOCK.GT.NGROUP) THEN
            WRITE(ITAPE6,*) ' TOO MANY BLOCKS INT, HXX,IJ'
            CALL MABORT
          END IF
          CALL RWRIT(ITAP52,IINT,INTOWP(NMAX),INTPT)
          INTPT=INTPT+I2SEC(INTOWP(NMAX))
          CALL RREAD(ITAP52,IINT,INTOWP(NMAX),INTPT)
C     CALL NXTBLK(ACRCY,ILBL,ILBM,IINT)
          GO TO 14
   15     CONTINUE
          HIJ=H(IJ)
          IJP=IJADD(IJ)
          IJS=IEOR(IS,ORBSYM(J))
          IJN=IJS*NORBS
          JA=J*(J-1)/2
          JJ=JA+J
          IJJ=IJP+KADD(IJN+J)
          IJJS=IS
          IJJN=IJJS*NORBS
          DO 17 L=1,J-1
            IF (ORBSYM(L).NE.IJJS) GO TO 16
            IL=IA+L
            JL=JA+L
            IJJL=IJJ+LADD(IJJN+L)
            T1=INT(IJJL+1)
            T2=INT(IJJL+2)
            HIL=H(IL)
            T=T1-T2-HIL
            CALL OUTXX(T,IJ,JL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
            T=SQRT2*T1
            CALL OUTWW(T,IL,JJ,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
            T=T1+T2+HIL
            CALL OUTWW(T,IJ,JL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
   16       CONTINUE
   17     CONTINUE
C
          DO 21 K=1,J-1
            IJK=IJP+KADD(IJN+K)
            IJKS=IEOR(IJS,ORBSYM(K))
            IJKN=IJKS*NORBS
            IK=IA+K
            JK=JA+K
            KA=K*(K-1)/2
            KK=KA+K
            IF (ORBSYM(K).NE.IJKS) GO TO 18
            IJKK=IJK+LADD(IJKN+K)
            T1=INT(IJKK+1)
            T2=INT(IJKK+2)
            T=T2-T1+HIJ
            CALL OUTXX(T,IK,JK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
            T=SQRT2*T1
            CALL OUTWW(T,IJ,KK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
            T=T1+T2+HIJ
            CALL OUTWW(T,IK,JK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
   18       CONTINUE
C
            DO 20 L=1,K-1
              IF (ORBSYM(L).NE.IJKS) GO TO 19
              IJKL=IJK+LADD(IJKN+L)
              T1=INT(IJKL+1)
              T2=INT(IJKL+2)
              T3=INT(IJKL+3)
              IL=IA+L
              JL=JA+L
              KL=KA+L
              CALL OUTXX(T1-T3,IJ,KL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
              CALL OUTWW(T1+T3,IJ,KL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
              CALL OUTXX(T2-T3,IK,JL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
              CALL OUTWW(T2+T3,IK,JL,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
              CALL OUTXX(T2-T1,IL,JK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
              CALL OUTWW(T2+T1,IL,JK,INT,IJXX,KLXX,IJWW,KLWW,IJGRP,IINT,
     #        iad,val,ival)
   19         CONTINUE
   20       CONTINUE
   21     CONTINUE
   22   CONTINUE
   23 CONTINUE
C
      CALL RWRIT(ITAP52,IINT,INTOWP(NMAX),INTPT)
C
      RETURN
      END
