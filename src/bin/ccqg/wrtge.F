C
C----------------------------------------------------------------------
C
      SUBROUTINE WRTGE(BUF,IBUF,GBUF,NBF,NO,NV,ITAP68,GMAT,NDIMG,
     .                 IOFF,NTR,JOUT,ICNT,EOFF,EADD,INO,NOV,
     .                 FLOV,ORBSYM,NIRRED,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      COMMON/BUFFER/TOL,IBFLEN,MAXBUF,IBOFF
      REAL*8 BUF(GBUF),GMAT(NDIMG),TOL
      INTEGER IBUF(GBUF*2),IOFF(NTR),EOFF(NO*NV),EADD(NO*NO),INO(NOV),
     .        FLOV(NIRRED,4),ORBSYM(NBF),FZO(NO),FZV(NV)
C
      DO 227 A=1,NV
         IF(FZV(A).EQ.1)GO TO 227
         ASYM=ORBSYM(A+NO)
         DO 226 I=1,NO
            IF(FZO(I).EQ.1)GO TO 226
            ISYM=ORBSYM(I)
            AISYM=IEOR(ASYM,ISYM)
            AI=INO(A)+I
            DO 225 J=1,NO
               IF(FZO(J).EQ.1)GO TO 225
               JSYM=ORBSYM(J)
               KSYM=IEOR(AISYM,JSYM)+1
               FK=FLOV(KSYM,1)
               LK=FLOV(KSYM,2)
               IF(LK.GT.J)LK=J
               DO 224 K=FK,LK
                  JK=INO(J)+K
                  AIJK=EOFF(AI)+EADD(JK)
COUT              WRITE (JOUT,1) A+NO,I,J,K,GMAT(AIJK)/8.0D+00
                  IF(DABS(GMAT(AIJK)).GT.TOL) THEN
                     LBL=IOR(I,ISHFT(A+NO,8))
                     LBL=IOR(J,ISHFT(LBL,8))
                     LBL=IOR(K,ISHFT(LBL,8))
                     IF(ICNT.GT.MAXBUF) THEN
                        CALL SWRIT(ITAP68,IBUF,INTOWP(GBUF))
                        CALL ZERO(BUF,GBUF)
                        ICNT=1
                     ENDIF
                     IBUF(1)=0
                     IBUF(2)=ICNT
                     IBUF(2+ICNT)=LBL
                     BUF(ICNT+IBOFF)=GMAT(AIJK)
                     ICNT=ICNT+1
                  ENDIF
  224          CONTINUE
  225       CONTINUE
  226    CONTINUE
  227 CONTINUE
C
    1 FORMAT ('GEMAT(',I3,',',I3,',',I3,',',I3,')=',F20.12)
      RETURN
      END
