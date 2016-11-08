C
C-----------------------------------------------------------------------
C
      SUBROUTINE WRTQQ(BUF,IBUF,GBUF,NBF,NO,NV,ITAP68,QOO,NTRO,QVV,NTRV,
     .                 QOV,IOFF,NTR,JOUT,ICNT,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      COMMON/BUFFER/TOL,IBFLEN,MAXBUF,IBOFF
      REAL*8 BUF(GBUF),TOL,QOO(NTRO),QVV(NTRV),QOV(NO,NV)
      INTEGER IBUF(GBUF*2),IOFF(NTR),FZO(NO),FZV(NV)
C
      IBUF(1)=1
      CALL SWRIT(ITAP68,BUF,INTOWP(GBUF))
      CALL ZERO(BUF,GBUF)
      ICNT=1
C
      DO 122 I=1,NO
         DO 120 J=1,I
            IJ=IOFF(I)+J
            IF(DABS(QOO(IJ)).GT.TOL) THEN
               LBL=IOR(J,ISHFT(I,8))
               IF(ICNT.GT.MAXBUF) THEN
                  CALL SWRIT(ITAP68,IBUF,INTOWP(GBUF))
                  CALL ZERO(BUF,GBUF)
                  ICNT=1
               ENDIF
               IBUF(1)=0
               IBUF(2)=ICNT
               IBUF(2+ICNT)=LBL
               BUF(ICNT+IBOFF)=QOO(IJ)
               ICNT=ICNT+1
            ENDIF
  120    CONTINUE
  122 CONTINUE
C
      DO 222 A=1,NV
         DO 220 B=1,A
            AB=IOFF(A)+B
            IF(DABS(QVV(AB)).GT.TOL) THEN
               LBL=IOR(B+NO,ISHFT(A+NO,8))
               IF(ICNT.GT.MAXBUF) THEN
                  CALL SWRIT(ITAP68,IBUF,INTOWP(GBUF))
                  CALL ZERO(BUF,GBUF)
                  ICNT=1
               ENDIF
               IBUF(1)=0
               IBUF(2)=ICNT
               IBUF(2+ICNT)=LBL
               BUF(ICNT+IBOFF)=QVV(AB)
               ICNT=ICNT+1
            ENDIF
  220    CONTINUE
  222 CONTINUE
C
      DO 322 I=1,NO
         DO 320 A=1,NV
            AI=(A-1)*NO+I
            IF(DABS(QOV(I,A)).GT.TOL) THEN
               LBL=IOR(I,ISHFT(A+NO,8))
               IF(ICNT.GT.MAXBUF) THEN
                  CALL SWRIT(ITAP68,IBUF,INTOWP(GBUF))
                  CALL ZERO(BUF,GBUF)
                  ICNT=1
               ENDIF
               IBUF(1)=0
               IBUF(2)=ICNT
               IBUF(2+ICNT)=LBL
               BUF(ICNT+IBOFF)=QOV(I,A)
               ICNT=ICNT+1
            ENDIF
  320    CONTINUE
  322 CONTINUE
C
      IBUF(1)=1
      CALL SWRIT(ITAP68,BUF,INTOWP(GBUF))
C
      RETURN
      END
