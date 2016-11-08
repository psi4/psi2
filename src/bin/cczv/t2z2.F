C
C-----------------------------------------------------------------------
C
      SUBROUTINE T2Z2(TZV,TZO,Z2O,T2,NO,NV,NO2,NV2,NTRO,NTRV,NTR,IOFF,
     .                JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,ZLX,UOFF,
     .                FLOV,ORBSYM,NIRRED,
     .                TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 Z2O(NDIMT2),TZV(NV,NV),TZO(NO,NO),T2(NDIMT2),VALO,VALV
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV),
     .        ZLX(NV,NV),UOFF(NO,NO,2),
     .        FLOV(NIRRED,4),ORBSYM(NO+NV)
C
      CALL ZERO(TZO,NO2)
      CALL ZERO(TZV,NV2)
      DO 10 I=1,NO
         IF(FZO(I).EQ.1)GO TO 10
         ISYM=ORBSYM(I)
         DO 20 J=1,NO
            IF(FZO(J).EQ.1)GO TO 20
            JSYM=ORBSYM(J)
            IJSYM=IEOR(ISYM,JSYM)
            DO 30 C=1,NV
               IF(FZV(C).EQ.1)GO TO 30
               CSYM=ORBSYM(C+NO)
               ASYM=IEOR(IJSYM,CSYM)+1
               FA=FLOV(ASYM,3)-NO
               LA=FLOV(ASYM,4)-NO
               DO 40 A=FA,LA
                  ZLA=ZLX(A,C)
                  IJAC=TOFF(I,J,ZLA)+TADD(A,C)
                  DO 50 B=FA,LA
                     ZLB=ZLX(B,C)
                     IJBC=TOFF(I,J,ZLB)+TADD(B,C)
                     VALV=T2(IJAC)*Z2O(IJBC)
                     TZV(A,B)=TZV(A,B)+VALV
   50             CONTINUE
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      DO 12 A=1,NV
         IF(FZV(A).EQ.1)GO TO 12
         ASYM=ORBSYM(A+NO)
         DO 22 B=1,NV
            IF(FZV(B).EQ.1)GO TO 22
            BSYM=ORBSYM(B+NO)
            ZL=ZLX(A,B)
            ABSYM=IEOR(ASYM,BSYM)
            DO 32 K=1,NO
               IF(FZO(K).EQ.1)GO TO 32
               KSYM=ORBSYM(K)
               ISYM=IEOR(KSYM,ABSYM)+1
               FI=FLOV(ISYM,1)
               LI=FLOV(ISYM,2)
               DO 42 I=FI,LI
                  DO 52 J=FI,LI
                     IKAB=TOFF(I,K,ZL)+TADD(A,B)
                     JKAB=TOFF(J,K,ZL)+TADD(A,B)
                     VALO=T2(IKAB)*Z2O(JKAB)
                     TZO(I,J)=TZO(I,J)+VALO
   52             CONTINUE
   42          CONTINUE
   32       CONTINUE
   22    CONTINUE
   12 CONTINUE
C
      RETURN
      END
