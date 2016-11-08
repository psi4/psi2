C
C----------------------------------------------------------------------
C
      SUBROUTINE FORMGB(BUF,IBUF,GBUF,TAU,Z2,NO,NV,NO2,NV2,NTRO,NTRV,
     .                  ICNT,IOFF,NTR,ZLX,JOUT,ITAP68,FLOV,ORBSYM,
     .                  NIRRED,NDIMT2,TOFF,TADD,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      COMMON/BUFFER/TOL,IBFLEN,MAXBUF,IBOFF
      REAL*8 VAL,TOL
      REAL*8 BUF(GBUF),Z2(NDIMT2),TAU(NDIMT2)
      INTEGER IOFF(NTR),ZLX(NV,NV),IBUF(GBUF*2),ORBSYM(NO+NV),
     .        FLOV(NIRRED,4),TOFF(NO,NO,2),TADD(NV,NV),
     .        FZO(NO),FZV(NV)
C
      DO 10 A=1,NV
         IF(FZV(A).EQ.1)GO TO 10
         ASYM=ORBSYM(A+NO)
         DO 20 B=1,A
            IF(FZV(B).EQ.1)GO TO 20
            BSYM=ORBSYM(B+NO)
            ABSYM=IEOR(ASYM,BSYM)
            DO 30 C=1,A
               IF(FZV(C).EQ.1)GO TO 30
               CSYM=ORBSYM(C+NO)
               ACSYM=IEOR(ASYM,CSYM)
               BCSYM=IEOR(BSYM,CSYM)
               DSYM=IEOR(ABSYM,CSYM)+1
               FD=FLOV(DSYM,3)-NO
               LD=FLOV(DSYM,4)-NO
               ZLBC=ZLX(B,C)
               MAXD=C
               IF (C.EQ.A) MAXD=B
               IF (LD.GT.MAXD) LD=MAXD
               DO 40 D=FD,LD
                  ZLBD=ZLX(B,D)
                  VAL=0.0D+00
                  DO 50 I=1,NO
                     IF(FZO(I).EQ.1)GO TO 50
                     ISYM=ORBSYM(I)
                     JSYM=IEOR(ISYM,ACSYM)+1
                     FJ=FLOV(JSYM,1)
                     LJ=FLOV(JSYM,2)
                     DO 60 J=FJ,LJ
                        IJBD=TOFF(I,J,ZLBD)+TADD(B,D)
                        IJAC=TOFF(I,J,1  )+TADD(A,C)
                        VAL=VAL+TAU(IJAC)*Z2(IJBD)+TAU(IJBD)*Z2(IJAC)
   60                CONTINUE
                     JSYM=IEOR(ISYM,BCSYM)+1
                     FJ=FLOV(JSYM,1)
                     LJ=FLOV(JSYM,2)
                     DO 61 J=FJ,LJ
                        IJAD=TOFF(I,J,1  )+TADD(A,D)
                        IJBC=TOFF(I,J,ZLBC)+TADD(B,C)
                        VAL=VAL+TAU(IJBC)*Z2(IJAD)+TAU(IJAD)*Z2(IJBC)
   61                CONTINUE
   50             CONTINUE
COUT              WRITE (JOUT,1) A+NO,B+NO,C+NO,D+NO,VAL/8.0D+00
                  IF(DABS(VAL).GT.TOL) THEN
                     LBL=IOR(B+NO,ISHFT(A+NO,8))
                     LBL=IOR(C+NO,ISHFT(LBL,8))
                     LBL=IOR(D+NO,ISHFT(LBL,8))
                     IF(ICNT.GT.MAXBUF) THEN
                        CALL SWRIT(ITAP68,IBUF,INTOWP(GBUF))
                        CALL ZERO(BUF,GBUF)
                        ICNT=1
                     ENDIF
                     IBUF(1)=0
                     IBUF(2)=ICNT
                     IBUF(2+ICNT)=LBL
                     BUF(ICNT+IBOFF)=VAL
                     ICNT=ICNT+1
                  ENDIF
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
    1 FORMAT ('GBMAT(',I3,',',I3,',',I3,',',I3,')=',F20.12)
      RETURN
      END
