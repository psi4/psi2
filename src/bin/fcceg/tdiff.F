      SUBROUTINE TDIFF(T2,T1,T2O,T1O,T1NORM,DELT,NO,NV,NDIMT2,
     .                 ORBSYM,FLOV,NIRRED,TOFF,TADD,ZLX,FZO,FZV,XE,XEO)
      IMPLICIT INTEGER(A-Z)
      REAL*8 XE,XEO,DELE,T1NORM,DELT
      REAL*8 T2(NDIMT2),T1(NO,NV),T2O(NDIMT2),T1O(NO,NV)
      INTEGER ORBSYM(NO+NV),FLOV(NIRRED,4),TOFF(NO,NO,2),TADD(NV,NV),
     .        ZLX(NV,NV),FZO(NO),FZV(NV)
C
      T1NORM=0.0D+00
      DELT=0.0D+00
C
      ncnt=0
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         ncnt=ncnt+1
         USYM=ORBSYM(U)
         DO 20 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 20
            BESYM=ORBSYM(BE+NO)
            IF(USYM.EQ.BESYM)THEN
               DELT=DELT+(T1(U,BE)-T1O(U,BE))*(T1(U,BE)-T1O(U,BE))
               T1NORM=T1NORM+T1(U,BE)*T1(U,BE)
            ENDIF
            UBESYM=IEOR(USYM,BESYM)
C           DO 30 V=1,NO
            DO 30 V=1,U
               IF(FZO(V).EQ.1)GO TO 30
               VSYM=ORBSYM(V)
               GASYM=IEOR(UBESYM,VSYM)+1
               FGA=FLOV(GASYM,3)-NO
               LGA=FLOV(GASYM,4)-NO
               LGA2=LGA
               IF(V.EQ.U)LGA2=BE
               DO 40 GA=FGA,LGA2
                  ZL=ZLX(BE,GA)
                  UVBEGA=TOFF(U,V,ZL)+TADD(BE,GA)
                  DELT=DELT
     .     +      (T2O(UVBEGA)-T2(UVBEGA))*(T2O(UVBEGA)-T2(UVBEGA))
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      DELT=DSQRT(DELT)
      T1NORM=DSQRT(T1NORM/(2*ncnt))
C
      RETURN
      END
