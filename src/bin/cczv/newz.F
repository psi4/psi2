C
C-----------------------------------------------------------------------
C
      SUBROUTINE NEWZ(Z2,Z1,Z2O,Z1O,NO,NV,NTR,NO2,NTRV,IOFF,JOUT,INO,
     .                INTRO,INTRV,NOV,INV,N2OV,
     .                ORBSYM,FLOV,NIRRED,ZLX,TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 Z2(NDIMT2  ),Z1(NO,NV),Z2O(NDIMT2  ),Z1O(NO,NV)
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER ORBSYM(NO+NV),FLOV(NIRRED,4),
     .        TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV)
C
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         USYM=ORBSYM(U)
         DO 20 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 20
            BESYM=ORBSYM(BE+NO)
            UBESYM=IEOR(USYM,BESYM)
            Z1O(U,BE)=Z1(U,BE)
            DO 30 V=1,NO
               IF(FZO(V).EQ.1)GO TO 30
               VSYM=ORBSYM(V)
               GASYM=IEOR(VSYM,UBESYM)+1
               FGA=FLOV(GASYM,3)-NO
               LGA=FLOV(GASYM,4)-NO
               IF(LGA.GT.BE)LGA=BE
               DO 40 GA=FGA,LGA
                  UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                  Z2O(UVBEGA)=Z2(UVBEGA)
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END
