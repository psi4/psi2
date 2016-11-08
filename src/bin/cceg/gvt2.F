C
C ----------------------------------------------------------------------
C
      SUBROUTINE GVT2  (T2,T2O,GV,XVAL,AUX3,
     .                  ORBSYM,FLOV,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  TMP
      REAL*8  T2(*),T2O(*),GV(NV,NV),XVAL(*),AUX3(NV,NV)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV)
C
      DO 191 U=1,NO
         IF(FZO(U).EQ.1)GO TO 191
         USYM=ORBSYM(U)
         DO 181 V=1,NO
            IF(FZO(V).EQ.1)GO TO 181
            VSYM=ORBSYM(V)
            VUSYM=IEOR(VSYM,USYM)
            DO 121 GA=1,NV
               IF(FZV(GA).EQ.1)GO TO 121
               GASYM=ORBSYM(GA+NO)
               ASYM=IEOR(GASYM,VUSYM)+1
               FA=FLOV(ASYM,3)-NO
               LA=FLOV(ASYM,4)-NO
               DO 111 A=FA,LA
                  ZL=ZLX(A,GA)
                  UVAGA=TOFF(U,V,ZL)+TADD(A,GA)
                  AUX3(A,GA)=T2O(UVAGA)
  111          CONTINUE
  121       CONTINUE
            DO 151 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 151
               BESYM=ORBSYM(BE+NO)
               FA=FLOV(BESYM+1,3)-NO
               LA=FLOV(BESYM+1,4)-NO
               DO 125 A=FA,LA
                  XVAL(A)=GV(BE,A)
  125          CONTINUE
               GASYM=IEOR(BESYM,VUSYM)+1
               FGA=FLOV(GASYM,3)-NO
               LGA=FLOV(GASYM,4)-NO
               DO 141 GA=FGA,LGA
                  TMP=0.0D0
                  DO 131 A=FA,LA
                     TMP=TMP+AUX3(A,GA)*XVAL(A)
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
