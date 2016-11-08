C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZDIFF(Z2,Z1,Z2O,Z1O,DELZ,NO,NV,NTR,NO2,NTRV,IOFF,JOUT,
     .                 INO,INTRO,INTRV,NOV,INV,N2OV,DELZ2,
     .                 ZLX,FLOV,ORBSYM,NIRRED,
     .                 TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 Z2(NDIMT2),Z1(NO,NV),Z2O(NDIMT2),Z1O(NO,NV),DELZ,DSQRT
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER ORBSYM(NO+NV),FLOV(NIRRED,4)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV)
      REAL*8 DELZ2,VAL1,VAL2
C
      DELZ=0.0D+00
      DELZ2=0.0D+00
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         USYM=ORBSYM(U)
         DO 20 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 20
            BESYM=ORBSYM(BE+NO)
            UBESYM=IEOR(USYM,BESYM)
            VAL1=(Z1O(U,BE)-Z1(U,BE))*(Z1O(U,BE)-Z1(U,BE))
            DELZ=DELZ+DSQRT(VAL1)
            DELZ2=DELZ2+VAL1
            DO 30 V=1,NO
               IF(FZO(V).EQ.1)GO TO 30
               VSYM=ORBSYM(V)
               GASYM=IEOR(VSYM,UBESYM)+1
               FGA=FLOV(GASYM,3)-NO
               LGA=FLOV(GASYM,4)-NO
               IF(LGA.GT.BE)LGA=BE
               DO 40 GA=FGA,LGA
                  UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                  VAL2=
     .   +      (Z2O(UVBEGA)-Z2(UVBEGA))*(Z2O(UVBEGA)-Z2(UVBEGA))
                  DELZ=DELZ+DSQRT(VAL2)
                  DELZ2=DELZ2+VAL2
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
      DELZ2=DSQRT(DELZ2)
C
      RETURN
      END
