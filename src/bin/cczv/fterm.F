C
C-----------------------------------------------------------------------
C
      SUBROUTINE FTERM(Z1,Z2O,FMAT,NO,NV,NO2,NV2,NTRO,NTRV,NTR,IOFF,
     .                 JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,ZLX,UOFF,
     .                 FLOV,ORBSYM,NIRRED,TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 Z1(NO,NV),Z2O(NDIMT2  ),FMAT(NV*NV)
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER ZLX(NV,NV),UOFF(NO,NO,2)
      INTEGER FLOV(NIRRED,4),ORBSYM(NO+NV)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
C
      CALL RFILE(TPFMAT)
      POINT=1
      DO 9 TSYM =1,NIRRED
      DO 10 BE=1,NV
         IF(FZV(BE).EQ.1)GO TO 10
         BESYM=ORBSYM(BE+NO)
         DO 20 A=1,NV
            IF(FZV(A).EQ.1)GO TO 20
            ASYM=ORBSYM(A+NO)
            BEASYM=IEOR(BESYM,ASYM)
            IF(BEASYM.NE.TSYM-1)GO TO 20
            CALL WREADW(TPFMAT,IDIM,1,POINT,POINT)
            CALL WREADW(TPFMAT,FMAT,INTOWP(IDIM),POINT,POINT)
            ICNT=0
            DO 30 B=1,NV
               IF(FZV(B).EQ.1)GO TO 30
               BSYM=ORBSYM(B+NO)
               ISYM=IEOR(BEASYM,BSYM)+1
               FI=FLOV(ISYM,1)
               LI=FLOV(ISYM,2)
               ZL=ZLX(A,B)
               DO 40 I=FI,LI
                  ICNT=ICNT+1
                     DO 50 U=1,NO
                        IF(FZO(U).EQ.1)GO TO 50
                     USYM=ORBSYM(U)
                     IF(USYM.NE.BESYM)GO TO 50
                     UIAB=TOFF(U,I,ZL)+TADD(A,B)
                     Z1(U,BE)=Z1(U,BE)+FMAT(ICNT)*Z2O(UIAB)
CF                   Z1(U,BE)=Z1(U,BE)+FMAT(BE,A,B,I)*Z2O(UI,AB)
   50             CONTINUE
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
    9 CONTINUE
C
      CALL RCLOSE(TPFMAT,3)
      RETURN
      END
