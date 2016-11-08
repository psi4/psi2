C
C ----------------------------------------------------------------------
C
      SUBROUTINE CDT1  (T1O,T2,CIN,DIN,XVAL,
     .                  ORBSYM,FLOV,COFF,CADD,DOFF,DADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL
      REAL*8  T1O(NO,NV),T2(*),CIN(*),DIN(*),XVAL(*)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        COFF(*),CADD(*),DOFF(*),DADD(*)
C
      DO 191 U=1,NO
         IF(FZO(U).EQ.1)GO TO 191
         USYM=ORBSYM(U)
         FA1=FLOV(USYM+1,3)-NO
         LA1=FLOV(USYM+1,4)-NO
         DO 181 V=1,NO
            IF(FZO(V).EQ.1)GO TO 181
            VSYM=ORBSYM(V)
            FA2=FLOV(VSYM+1,3)-NO
            LA2=FLOV(VSYM+1,4)-NO
            VUSYM=IEOR(VSYM,USYM)
            DO 171 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 171
               BESYM=ORBSYM(BE+NO)
               ISYM=IEOR(BESYM,VUSYM)+1
               FI=FLOV(ISYM,1)
               LI=FLOV(ISYM,2)
               FGA=FLOV(ISYM,3)-NO
               LGA=FLOV(ISYM,4)-NO
               BEU=ITR(BE)+U
               DO 134 I=FI,LI
                  XVAL(I)=0.0D0
                  DO 131 A=FA2,LA2
                     IA=ITR(A)+I
                     BEUIA=DOFF(MAX0(BEU,IA))+DADD(MIN0(BEU,IA))
                     XVAL(I)=XVAL(I)-DIN(BEUIA)*T1O(V,A)
  131             CONTINUE
                  IV=ITR(MAX0(I,V))+MIN0(I,V)
                  DO 133 A=FA1,LA1
                     BEA=ITV(MAX0(BE,A))+MIN0(BE,A)
                     IVBEA=COFF(BEA)+CADD(IV)
                     XVAL(I)=XVAL(I)-CIN(IVBEA)*T1O(U,A)
  133             CONTINUE
  134          CONTINUE
               DO 161 GA=FGA,LGA
                  VAL=0.0D0
                  DO 143 I=FI,LI
                     VAL=VAL+XVAL(I)*T1O(I,GA)
  143             CONTINUE
                  IF(BE.GE.GA)THEN
                     UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                     T2(UVBEGA)=T2(UVBEGA)+VAL
                  ENDIF
                  IF(BE.LE.GA)THEN
                     VUGABE=TOFF(V,U,1)+TADD(GA,BE)
                     T2(VUGABE)=T2(VUGABE)+VAL
                  ENDIF
  161          CONTINUE
  171       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
      RETURN
      END
