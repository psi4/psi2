C
C-----------------------------------------------------------------------
C
      SUBROUTINE ETERM(Z1,Z2O,EMAT,NO,NV,NO2,NV2,NTRO,NTRV,NTR,IOFF,
     .                 JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,ZLX,UOFF,
     .                 FLOV,ORBSYM,NIRRED,TOFF,TADD,NDIMT2,
     .                 EMOFF,EMADD,LEMAT,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 Z1(NO,NV),Z2O(NDIMT2),EMAT(LEMAT)
      INTEGER IOFF(NTR),FZO(*),FZV(*),EMOFF(NO*NV),EMADD(NO*NO)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER ZLX(NV,NV),UOFF(NO,NO,2)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
      INTEGER FLOV(NIRRED,4),ORBSYM(NO+NV)
C
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         USYM=ORBSYM(U)
         DO 20 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 20
            BESYM=ORBSYM(BE+NO)
            IF(USYM.NE.BESYM)GO TO 20
            DO 30 I=1,NO
               IF(FZO(I).EQ.1)GO TO 30
               UIX=INO(U)+I
               ISYM=ORBSYM(I)
               IUSYM=IEOR(ISYM,USYM)
               DO 40 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 40
                  JSYM=ORBSYM(J)
                  ASYM=IEOR(IUSYM,JSYM)+1
                  FA=FLOV(ASYM,3)-NO
                  LA=FLOV(ASYM,4)-NO
                  DO 50 A=FA,LA
                     ZL=ZLX(BE,A)
                     IJBEA=TOFF(I,J,ZL)+TADD(BE,A)
                     JAX=INO(A)+J
                     UIJA=EMOFF(JAX)+EMADD(UIX)
                     Z1(U,BE)=Z1(U,BE)+EMAT(UIJA)*Z2O(IJBEA)
   50             CONTINUE
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END
