C
C----------------------------------------------------------------------
C
      SUBROUTINE WRTGD(BUF,IBUF,GBUF,NBF,NO,NV,ITAP68,GMAT,NDIMG,
     .                 IOFF,NTR,JOUT,ICNT,DOFF,DADD,ITR,
     .                 FLOV,ORBSYM,NIRRED,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      COMMON/BUFFER/TOL,IBFLEN,MAXBUF,IBOFF
      REAL*8 BUF(GBUF),GMAT(NDIMG),TOL
      INTEGER IBUF(GBUF*2),IOFF(NTR),ITR(NV),DOFF(NO*NV),DADD(NO*NV),
     .        FLOV(NIRRED,4),ORBSYM(NBF),FZO(NO),FZV(NV)
C
C     INITIALIZE FILE68 AND DEFINE CONSTANTS
C
      TOL=1.0D-12
      IBFLEN=(INTOWP(GBUF)-2)/INTOWP(1)
      MAXBUF=INTOWP(IBFLEN)/(1+INTOWP(1))
      IBOFF=(MAXBUF+3)/INTOWP(1)
      WRITE(JOUT,*) ' BUFFER SIZE IS',GBUF,' REAL WORDS'
C     WRITE(JOUT,*) ' TO HOLD ',MAXBUF,' TERMS AT A TIME'
C
      CALL ZERO(BUF,GBUF)
      CALL SREW(ITAP68)
      ICNT=1
C
      DO 229 A=1,NV
         IF(FZV(A).EQ.1)GO TO 229
         ASYM=ORBSYM(A+NO)
         DO 228 I=1,NO
            IF(FZO(I).EQ.1)GO TO 228
            ISYM=ORBSYM(I)
            AISYM=IEOR(ASYM,ISYM)
            AI=ITR(A)+I
            DO 227 B=1,A
               IF(FZV(B).EQ.1)GO TO 227
               BSYM=ORBSYM(B+NO)
               JSYM=IEOR(AISYM,BSYM)+1
               FJ=FLOV(JSYM,1)
               LJ=FLOV(JSYM,2)
               LJ2=LJ
               IF(B.EQ.A)LJ2=I
               DO 226 J=FJ,LJ2
                  BJ=ITR(B)+J
                  AIBJ=DOFF(AI)+DADD(BJ)
c                 WRITE (JOUT,1) A+NO,I,B+NO,J,AIBJ,GMAT(AIBJ)/8.0D+00
                  IF(DABS(GMAT(AIBJ)).GT.TOL) THEN
                     LBL=IOR(I,ISHFT(A+NO,8))
                     LBL=IOR(B+NO,ISHFT(LBL,8))
                     LBL=IOR(J,ISHFT(LBL,8))
                     IF(ICNT.GT.MAXBUF) THEN
                        CALL SWRIT(ITAP68,IBUF,INTOWP(GBUF))
                        CALL ZERO(BUF,GBUF)
                        ICNT=1
                     ENDIF
                     IBUF(1)=0
                     IBUF(2)=ICNT
                     IBUF(2+ICNT)=LBL
                     BUF(ICNT+IBOFF)=GMAT(AIBJ)
                     ICNT=ICNT+1
                  ENDIF
  226          CONTINUE
  227       CONTINUE
  228    CONTINUE
  229 CONTINUE
C
    1 FORMAT ('GMAT(',I3,',',I3,',',I3,',',I3,')',I6,2X,F20.12)
      RETURN
      END
