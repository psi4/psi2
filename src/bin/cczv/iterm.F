C
C-----------------------------------------------------------------------
C
      SUBROUTINE ITERM(Z2,Z1O,I1MAT,NO,NV,NO2,NV2,NTRO,NTRV,NTR,IOFF,
     .                 JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,
     .                 ZLX,FLOV,ORBSYM,NIRRED,
     .                 TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 Z2(NDIMT2),Z1O(NO,NV),I1MAT(NO,NV)
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER ORBSYM(NO+NV),FLOV(NIRRED,4)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV)
C
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         USYM=ORBSYM(U)
         DO 20 V=1,NO
            IF(FZO(V).EQ.1)GO TO 20
            VSYM=ORBSYM(V)
            UVSYM=IEOR(USYM,VSYM)
            DO 30 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 30
               BESYM=ORBSYM(BE+NO)
               GASYM=IEOR(BESYM,UVSYM)+1
               FGA=FLOV(GASYM,3)-NO
               LGA=FLOV(GASYM,4)-NO
               IF(LGA.GT.BE)LGA=BE
               DO 40 GA=FGA,LGA
                  UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                  Z2(UVBEGA)=Z2(UVBEGA)+I1MAT(U,GA)*Z1O(V,BE)+
     .                        I1MAT(V,BE)*Z1O(U,GA)-I1MAT(V,GA)*
     .                        Z1O(U,BE)-I1MAT(V,GA)*Z1O(U,BE)-
     .                        I1MAT(U,BE)*Z1O(V,GA)-I1MAT(U,BE)*
     .                        Z1O(V,GA)
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END
