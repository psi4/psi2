C
C-----------------------------------------------------------------------
C
      SUBROUTINE GTERM(Z2,Z1,Z1O,Z2O,GOMAT,GVMAT,NO,NV,NO2,NV2,NTRO,NTRV
     .                 ,NTR,IOFF,JOUT,FOCK,T1,INO,INTRO,INTRV,NOV,INV,
     .                 N2OV,ZLX,UOFF,FLOV,ORBSYM,NIRRED,TOFF,TADD,
     .                 NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 Z2(NDIMT2  ),Z1(NO,NV),Z2O(NDIMT2  ),Z1O(NO,NV),
     .       GOMAT(NO,NO),GVMAT(NV,NV)
      REAL * 8 FOCK(NTR),VAL,T1(NO,NV)
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER ZLX(NV,NV),UOFF(NO,NO,2)
      INTEGER FLOV(NIRRED,4),ORBSYM(NO+NV)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
C
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
      USYM=ORBSYM(U)
      FI=FLOV(USYM+1,1)
      LI=FLOV(USYM+1,2)
         DO 20 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 20
         BESYM=ORBSYM(BE+NO)
         FA2=FLOV(BESYM+1,3)-NO
         LA2=FLOV(BESYM+1,4)-NO
         UBESYM=IEOR(USYM,BESYM)
         UBE=IOFF(BE+NO)+U
         IF(USYM.NE.BESYM) GO TO 21
         Z1(U,BE)=Z1(U,BE)-2.D0*FOCK(UBE)
            DO 30 A=1,NV
               IF(FZV(A).EQ.1)GO TO 30
            Z1(U,BE)=Z1(U,BE)+GVMAT(A,BE)*Z1O(U,A)
   30       CONTINUE
            DO 40 I=1,NO
               IF(FZO(I).EQ.1)GO TO 40
            Z1(U,BE)=Z1(U,BE)+GOMAT(U,I)*Z1O(I,BE)
   40       CONTINUE
   21       CONTINUE
            DO 50 GA=1,BE
               IF(FZV(GA).EQ.1)GO TO 50
            GASYM=ORBSYM(GA+NO)
            FA1=FLOV(GASYM+1,3)-NO
            LA1=FLOV(GASYM+1,4)-NO
            VSYM=IEOR(GASYM,UBESYM)+1
            FV=FLOV(VSYM,1)
            LV=FLOV(VSYM,2)
               DO 60 V=FV,LV
               UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                  DO 70 A=FA1,LA1
                     ZLB=ZLX(BE,A)
                     UVBEA=TOFF(U,V,ZLB)+TADD(BE,A)
                     Z2(UVBEGA)=Z2(UVBEGA)+GVMAT(A,GA)*Z2O(UVBEA)
   70             CONTINUE
                  DO 71 A=FA2,LA2
                     ZLG=ZLX(GA,A)
                     VUGAA=TOFF(V,U,ZLG)+TADD(GA,A)
                     Z2(UVBEGA)=Z2(UVBEGA)+GVMAT(A,BE)*Z2O(VUGAA)
   71             CONTINUE
                  DO 100 I=FV,LV
                     UIBEGA=TOFF(U,I,1)+TADD(BE,GA)
                     Z2(UVBEGA)=Z2(UVBEGA)+GOMAT(V,I)*Z2O(UIBEGA)
  100             CONTINUE
                  DO 101 I=FI,LI
                     IVBEGA=TOFF(I,V,1)+TADD(BE,GA)
                     Z2(UVBEGA)=Z2(UVBEGA)+GOMAT(U,I)*Z2O(IVBEGA)
  101             CONTINUE
   60          CONTINUE
   50       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END
