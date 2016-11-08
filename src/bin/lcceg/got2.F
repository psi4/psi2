C
C ----------------------------------------------------------------------
C
      SUBROUTINE GOT2  (T2,T2O,GO,XVAL,AUX2,
     .                  ORBSYM,FLOV,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  TMP
      REAL*8  T2(*),T2O(*),GO(NO,NO),XVAL(*),AUX2(NO,NO)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV)
C
      DO 191 BE=1,NV
         IF(FZV(BE).EQ.1)GO TO 191
         BESYM=ORBSYM(BE+NO)
         DO 181 GA=1,NV
            IF(FZV(GA).EQ.1)GO TO 181
            GASYM=ORBSYM(GA+NO)
            BEGAS=IEOR(BESYM,GASYM)
            ZL=ZLX(BE,GA)
            DO 121 V=1,NO
               IF(FZO(V).EQ.1)GO TO 121
               VSYM=ORBSYM(V)
               ISYM=IEOR(VSYM,BEGAS)+1
               FI=FLOV(ISYM,1)
               LI=FLOV(ISYM,2)
               DO 111 I=FI,LI
                  IVBEGA=TOFF(I,V,ZL)+TADD(BE,GA)
                  AUX2(I,V)=T2O(IVBEGA)
  111          CONTINUE
  121       CONTINUE
            DO 151 U=1,NO
               IF(FZO(U).EQ.1)GO TO 151
               USYM=ORBSYM(U)
               FI=FLOV(USYM+1,1)
               LI=FLOV(USYM+1,2)
               DO 125 I=FI,LI
                  XVAL(I)=-GO(U,I)
  125          CONTINUE
               VSYM=IEOR(USYM,BEGAS)+1
               FV=FLOV(VSYM,1)
               LV=FLOV(VSYM,2)
               DO 141 V=FV,LV
                  TMP=0.0D0
                  DO 131 I=FI,LI
                     TMP=TMP+AUX2(I,V)*XVAL(I)
  131             CONTINUE
                  IF(BE.GE.GA)THEN
                     UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                     T2(UVBEGA)=T2(UVBEGA)+TMP
                  ENDIF
                  IF(BE.LE.GA)THEN
                     VUGABE=TOFF(V,U,1)+TADD(GA,BE)
                     T2(VUGABE)=T2(VUGABE)+TMP
                  ENDIF
  141          CONTINUE
  151       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
      RETURN
      END
