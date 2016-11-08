C
C-----------------------------------------------------------------------
C
      SUBROUTINE NORMZ(Z2,Z1,EIGVAL,NO,NV,NBF,NTR,NO2,NTRV,IOFF,JOUT,
     .                 INO,INTRO,INTRV,NOV,INV,N2OV,
     .                 ZLX,FLOV,ORBSYM,NIRRED,
     .                 TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 Z2(NDIMT2  ),Z1(NO,NV),EIGVAL(NBF)
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER ORBSYM(NO+NV),FLOV(NIRRED,4)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV)
C
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         USYM=ORBSYM(U)
         UU  =IOFF(U)+U
         DO 20 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 20
            BESYM=ORBSYM(BE+NO)
            UBESYM=IEOR(USYM,BESYM)
            BEBE=IOFF(BE+NO)+BE+NO
            Z1(U,BE)=Z1(U,BE)/(EIGVAL(BEBE)-EIGVAL(UU))
            DO 30 V=1,NO
               IF(FZO(V).EQ.1)GO TO 30
               VSYM=ORBSYM(V)
               VV  =IOFF(V)+V
               GASYM=IEOR(VSYM,UBESYM)+1
               FGA=FLOV(GASYM,3)-NO
               LGA=FLOV(GASYM,4)-NO
               IF(LGA.GT.BE)LGA=BE
               DO 40 GA=FGA,LGA
                  GAGA=IOFF(GA+NO)+GA+NO
                  UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                  Z2(UVBEGA)=Z2(UVBEGA)/(EIGVAL(BEBE)+EIGVAL(GAGA)-
     .                        EIGVAL(UU)-EIGVAL(VV))
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END
