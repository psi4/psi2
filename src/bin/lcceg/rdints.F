C
C----------------------------------------------------------------------
C
      SUBROUTINE RDINTS(IDISK,MOINTS,LENINT,BUF,IBUF,LENBUF,NO,NV,
     .                  IOFF,IADD,IPP,IQQ,FOCK)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MOINTS(LENINT),BUF(LENBUF),FOCK(*),VAL
      COMMON/INTS/NSGOO,NSGVV,NSGOV,NSHOV,NSLOV,NSLVO
      INTEGER IBUF(2*LENBUF),
     .        IOFF(*),IADD(*),IPP(*),IQQ(*)
C
      IND(I,J)=I*(I-1)/2+J
C
      CALL SREW(IDISK)
      CALL ZERO(MOINTS,LENINT)
C
      ITEMP=255
      INTLEN=(INTOWP(LENBUF)-2)/INTOWP(1)
      MAXVAL=INTOWP(INTLEN)/(1+INTOWP(1))
      IBOFF=(MAXVAL+2)/INTOWP(1)
C
C***  A INTS  **** (OO|OO) *****
C
      IF (IDISK.NE.60) GOTO 110
  111 CONTINUE
      CALL SREAD(IDISK,IBUF,INTOWP(LENBUF))
      IFLG=IBUF(1)
      MBUF=IBUF(2)
      DO 12 XI=1,MBUF
         IJKL=IBUF(XI+2)
         I=ISHFT(IJKL,-24)
         J=IAND(ITEMP,ISHFT(IJKL,-16))
         IJ=IPP(I)+J
         K=IAND(ITEMP,ISHFT(IJKL,-8))
         L=IAND(ITEMP,IJKL)
         KL=IQQ(K)+L
         IJKL=IOFF(IJ)+IADD(KL)
         VAL=BUF(XI+IBOFF)
         IF(IJKL.LE.NSGOO) THEN
            MOINTS(IJKL)=VAL
C           WRITE(6,19)I,J,K,L,ITYPE(I,J,K,L)
19          FORMAT(4I6,5X,I6)
         ELSE
C           WRITE(6,29)I,J,K,L,ITYPE(I,J,K,L)
29          FORMAT(4I6,' F',3X,I6)
            TY=ITYPE(I,J,K,L)
            IF(TY.EQ.1.OR.TY.EQ.2.OR.TY.EQ.3)THEN
               IL=IND(I,L)
               FOCK(IL)=FOCK(IL)+VAL
               GO TO 912
            ELSE IF(TY.EQ.4)THEN
               II=IND(I,I)
               FOCK(II)=FOCK(II)+VAL+VAL
               LL=IND(L,L)
               FOCK(LL)=FOCK(LL)+VAL+VAL
               GO TO 912
            ELSE IF(TY.EQ.5)THEN
               II=IND(I,I)
               FOCK(II)=FOCK(II)-VAL
               LL=IND(L,L)
               FOCK(LL)=FOCK(LL)-VAL
               GO TO 912
            ELSE IF(TY.EQ.6)THEN
               KL=IND(K,L)
               FOCK(KL)=FOCK(KL)+VAL+VAL
               GO TO 912
            ELSE IF(TY.EQ.7)THEN
               JL=IND(J,L)
               FOCK(JL)=FOCK(JL)-VAL
               GO TO 912
            ELSE IF(TY.EQ.8)THEN
               IL=IND(I,L)
               FOCK(IL)=FOCK(IL)-VAL
               GO TO 912
            ELSE IF(TY.EQ.9)THEN
               IK=IND(I,K)
               FOCK(IK)=FOCK(IK)-VAL
               GO TO 912
            ELSE IF(TY.EQ.10.OR.TY.EQ.11)THEN
               IJ=IND(I,J)
               FOCK(IJ)=FOCK(IJ)+VAL+VAL
               GO TO 912
            ENDIF
  912       CONTINUE
         ENDIF
   12 CONTINUE
      IF (IFLG.EQ.0) GOTO 111
  110 CONTINUE
C
C***  B INTS  **** (VV|VV) *****
C
      IF (IDISK.NE.61) GOTO 120
  112 CONTINUE
      CALL SREAD(IDISK,IBUF,INTOWP(LENBUF))
      IFLG=IBUF(1)
      MBUF=IBUF(2)
      DO 22 II=1,MBUF
         IJKL=IBUF(II+2)
         I=ISHFT(IJKL,-24)-NO
         J=IAND(ITEMP,ISHFT(IJKL,-16))-NO
         IJ=IPP(I)+J
         K=IAND(ITEMP,ISHFT(IJKL,-8))-NO
         L=IAND(ITEMP,IJKL)-NO
         KL=IPP(K)+L
         IJKL=IOFF(IJ)+IADD(KL)
         IF(IJKL.GT.NSGVV)GO TO 22
         MOINTS(IJKL)=BUF(II+IBOFF)
   22 CONTINUE
      IF (IFLG.EQ.0) GOTO 112
  120 CONTINUE
C
C***  C INTS  **** (VV|OO) *****
C
      IF (IDISK.NE.62) GOTO 130
  113 CONTINUE
      CALL SREAD(IDISK,IBUF,INTOWP(LENBUF))
      IFLG=IBUF(1)
      MBUF=IBUF(2)
      DO 32 XI=1,MBUF
         IJKL=IBUF(XI+2)
         I=ISHFT(IJKL,-24)-NO
         J=IAND(ITEMP,ISHFT(IJKL,-16))-NO
         IJ=IPP(I)+J
         K=IAND(ITEMP,ISHFT(IJKL,-8))
         L=IAND(ITEMP,IJKL)
         KL=IQQ(K)+L
         IJKL=IOFF(IJ)+IADD(KL)
         VAL=BUF(XI+IBOFF)
         IF(IJKL.LE.NSGOV)THEN
            MOINTS(IJKL)=VAL
         ELSE IF(K.EQ.L) THEN
            IJ=IND(I+NO,J+NO)
            FOCK(IJ)=FOCK(IJ)+VAL+VAL
         ENDIF
   32 CONTINUE
      IF (IFLG.EQ.0) GOTO 113
  130 CONTINUE
C
C***  D INTS  **** (VO|VO) *****
C
      IF (IDISK.NE.63) GOTO 140
  114 CONTINUE
      CALL SREAD(IDISK,IBUF,INTOWP(LENBUF))
      IFLG=IBUF(1)
      MBUF=IBUF(2)
      DO 42 XI=1,MBUF
         IJKL=IBUF(XI+2)
         I=ISHFT(IJKL,-24)-NO
         J=IAND(ITEMP,ISHFT(IJKL,-16))
         IJ=IPP(I)+J
         K=IAND(ITEMP,ISHFT(IJKL,-8))-NO
         L=IAND(ITEMP,IJKL)
         KL=IQQ(K)+L
         IJKL=IOFF(IJ)+IADD(KL)
         VAL=BUF(XI+IBOFF)
         IF(IJKL.LE.NSHOV)THEN
            MOINTS(IJKL)=VAL
         ELSE IF(J.EQ.L) THEN
            IK=IND(I+NO,K+NO)
            FOCK(IK)=FOCK(IK)-VAL
         ENDIF
   42 CONTINUE
      IF (IFLG.EQ.0) GOTO 114
  140 CONTINUE
C
C***  E INTS  **** (VO|OO) *****
C
      IF (IDISK.NE.64) GOTO 150
  115 CONTINUE
      CALL SREAD(IDISK,IBUF,INTOWP(LENBUF))
      IFLG=IBUF(1)
      MBUF=IBUF(2)
      DO 52 XI=1,MBUF
         IJKL=IBUF(XI+2)
         I=ISHFT(IJKL,-24)-NO
         J=IAND(ITEMP,ISHFT(IJKL,-16))
         IJ=IPP(I)+J
         K=IAND(ITEMP,ISHFT(IJKL,-8))
         L=IAND(ITEMP,IJKL)
         KL=IQQ(K)+L
         IJKL=IOFF(IJ)+IADD(KL)
         VAL=BUF(XI+IBOFF)
         IF(IJKL.LE.NSLOV)THEN
            MOINTS(IJKL)=VAL
         ENDIF
         IF(IJKL.GT.NSLOV)THEN
            IF(K.EQ.L) THEN
               IJ=IND(I+NO,J)
               FOCK(IJ)=FOCK(IJ)+VAL+VAL
            ENDIF
            IF(K.EQ.J) THEN
               IL=IND(I+NO,L)
               FOCK(IL)=FOCK(IL)-VAL
            ELSE IF(J.EQ.L) THEN
               IK=IND(I+NO,K)
               FOCK(IK)=FOCK(IK)-VAL
            ENDIF
         ENDIF
   52 CONTINUE
      IF (IFLG.EQ.0) GOTO 115
  150 CONTINUE
C
      IF (IDISK.NE.65) GOTO 160
      CALL SREW(IDISK+1)
C
C***  F INTS  **** (VV|VO) *****
C
  116 CONTINUE
      CALL SREAD(IDISK,IBUF,INTOWP(LENBUF))
      IFLG=IBUF(1)
      MBUF=IBUF(2)
      DO 62 II=1,MBUF
         IJKL=IBUF(II+2)
         I=ISHFT(IJKL,-24)-NO
         J=IAND(ITEMP,ISHFT(IJKL,-16))-NO
         IJ=IPP(I)+J
         K=IAND(ITEMP,ISHFT(IJKL,-8))-NO
         L=IAND(ITEMP,IJKL)
         KL=IQQ(K)+L
         IJKL=IADD(KL)+IOFF(IJ)
         IF(IJKL.GT.NSLVO)GO TO 62
         MOINTS(IJKL)=BUF(II+IBOFF)
   62 CONTINUE
      IF (IFLG.EQ.0) GOTO 116
C
C***  F INTS  **** (VO|VV) *****
C
  117 CONTINUE
      CALL SREAD(IDISK+1,BUF,INTOWP(LENBUF))
      IFLG=IBUF(1)
      MBUF=IBUF(2)
      DO 72 II=1,MBUF
         IJKL=IBUF(II+2)
         I=ISHFT(IJKL,-24)-NO
         J=IAND(ITEMP,ISHFT(IJKL,-16))
         IJ=IQQ(I)+J
         K=IAND(ITEMP,ISHFT(IJKL,-8))-NO
         L=IAND(ITEMP,IJKL)-NO
         KL=IPP(K)+L
         IJKL=IADD(IJ)+IOFF(KL)
         IF(IJKL.GT.NSLVO)GO TO 72
         MOINTS(IJKL)=BUF(II+IBOFF)
   72 CONTINUE
      IF (IFLG.EQ.0) GOTO 117
  160 CONTINUE
C
      RETURN
      END
