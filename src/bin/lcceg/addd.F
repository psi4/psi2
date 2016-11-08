C
C ----------------------------------------------------------------------
C
      SUBROUTINE ADDD  (T2,DIN,XVAL,
     .                  ORBSYM,FLOV,DOFF,DADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  TMP
      REAL*8  T2(*),DIN(*),XVAL(*)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        DOFF(*),DADD(*)
C
      DO 191 U=1,NO
         IF(FZO(U).EQ.1)GO TO 191
         USYM=ORBSYM(U)
         DO 181 V=1,NO
            IF(FZO(V).EQ.1)GO TO 181
            VSYM=ORBSYM(V)
            VUSYM=IEOR(VSYM,USYM)
            DO 171 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 171
               BESYM=ORBSYM(BE+NO)
               GASYM=IEOR(BESYM,VUSYM)+1
               FGA=FLOV(GASYM,3)-NO
               LGA=FLOV(GASYM,4)-NO
               BEU=ITR(BE)+U
               IF(LGA.GT.BE)LGA=BE
               DO 161 GA=FGA,LGA
                  GAV=ITR(GA)+V
                  BEUGAV=DOFF(MAX0(BEU,GAV))+DADD(MIN0(BEU,GAV))
                  XVAL(GA)=DIN(BEUGAV)
  161          CONTINUE
               DO 162 GA=FGA,LGA
                  TMP = XVAL(GA)
                  UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                  T2(UVBEGA)=T2(UVBEGA)+TMP
  162          CONTINUE
  171       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
      RETURN
      END
