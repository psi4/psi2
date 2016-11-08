C
C----------------------------------------------------------------------
C
      SUBROUTINE WRTGC(BUF,IBUF,GBUF,NBF,NO,NV,ITAP68,GMAT,NDIMG,
     .                 IOFF,NTR,JOUT,ICNT,COFF,CADD,INO,INV,NOV,
     .                 FLOV,ORBSYM,NIRRED,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      COMMON/BUFFER/TOL,IBFLEN,MAXBUF,IBOFF
      REAL*8 BUF(GBUF),GMAT(NDIMG),TOL
      INTEGER IBUF(GBUF*2),IOFF(NTR),COFF(NV*NV),CADD(NO*NO),
     .        INO(NOV),INV(NOV),FLOV(NIRRED,4),ORBSYM(NBF),FZO(*),FZV(*)
C
      DO 225 A=1,NV
         IF(FZV(A).EQ.1)GO TO 225
         ASYM=ORBSYM(A+NO)
         DO 224 B=1,A
            IF(FZV(B).EQ.1)GO TO 224
            BSYM=ORBSYM(B+NO)
            ABSYM=IEOR(ASYM,BSYM)
            AB=INV(A)+B
            DO 223 I=1,NO
               IF(FZO(I).EQ.1)GO TO 223
               ISYM=ORBSYM(I)
               JSYM=IEOR(ABSYM,ISYM)+1
               FJ=FLOV(JSYM,1)
               LJ=FLOV(JSYM,2)
               IF(LJ.GT.I)LJ=I
               DO 222 J=FJ,LJ
                  IJ=INO(I)+J
                  ABIJ=COFF(AB)+CADD(IJ)
c                 WRITE(JOUT,1) A+NO,B+NO,I,J,GMAT(ABIJ)/8.0D+00
                  IF(DABS(GMAT(ABIJ)).GT.TOL) THEN
                     LBL=IOR(B+NO,ISHFT(A+NO,8))
                     LBL=IOR(I,ISHFT(LBL,8))
                     LBL=IOR(J,ISHFT(LBL,8))
                     IF(ICNT.GT.MAXBUF) THEN
                        CALL SWRIT(ITAP68,IBUF,INTOWP(GBUF))
                        CALL ZERO(BUF,GBUF)
                        ICNT=1
                     ENDIF
                     IBUF(1)=0
                     IBUF(2)=ICNT
                     IBUF(2+ICNT)=LBL
                     BUF(ICNT+IBOFF)=GMAT(ABIJ)
                     ICNT=ICNT+1
                  ENDIF
  222          CONTINUE
  223       CONTINUE
  224    CONTINUE
  225 CONTINUE
C
    1 FORMAT ('GMAT(',I3,',',I3,',',I3,',',I3,')=',F20.12)
      RETURN
      END
