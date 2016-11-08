C
C-----------------------------------------------------------------------
C
      SUBROUTINE MTERM(Z2,TZO,TZV,NO,NV,NO2,NV2,NTRO,NTRV,NTR,IOFF,
     .                 MOINTS,LENINT,JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,
     .                 FLOV,ORBSYM,NIRRED,
     .                 DOFF,NOFF,DADD,NADD,TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 Z2(NDIMT2  ),TZO(NO,NO),TZV(NV,NV),MOINTS(LENINT),VALO,VALV
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER DOFF(NOFF),DADD(NADD)
      INTEGER FLOV(NIRRED,4),ORBSYM(NO+NV)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
C
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         USYM=ORBSYM(U)
         FIU=FLOV(USYM+1,1)
         LIU=FLOV(USYM+1,2)
         DO 20 V=1,NO
            IF(FZO(V).EQ.1)GO TO 20
            VSYM=ORBSYM(V)
            FIV=FLOV(VSYM+1,1)
            LIV=FLOV(VSYM+1,2)
            UVSYM=IEOR(USYM,VSYM)
            DO 30 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 30
               BESYM=ORBSYM(BE+NO)
               FA=FLOV(BESYM+1,3)-NO
               LA=FLOV(BESYM+1,4)-NO
               BEV=INO(BE)+V
               BEU=INO(BE)+U
               GASYM=IEOR(UVSYM,BESYM)+1
               FGA=FLOV(GASYM,3)-NO
               LGA=FLOV(GASYM,4)-NO
               LGA2=LGA
               IF(LGA.GT.BE)LGA2=BE
               DO 40 GA=FGA,LGA2
                  GAV=INO(GA)+V
                  GAU=INO(GA)+U
                  UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                  DO 50 A=FA,LA
                     AU=INO(A)+U
                     AV=INO(A)+V
                     GAVAU=DOFF(MAX0(GAV,AU))+DADD(MIN0(GAV,AU))
                     GAUAV=DOFF(MAX0(GAU,AV))+DADD(MIN0(GAU,AV))
                     VALV=(MOINTS(GAVAU)+MOINTS(GAVAU)-MOINTS(GAUAV))*
     .                    TZV(A,BE)
                     Z2(UVBEGA)=Z2(UVBEGA)+VALV
   50             CONTINUE
                  DO 51 A=FGA,LGA
                     AU=INO(A)+U
                     AV=INO(A)+V
                     BEVAU=DOFF(MAX0(BEV,AU))+DADD(MIN0(BEV,AU))
                     BEUAV=DOFF(MAX0(BEU,AV))+DADD(MIN0(BEU,AV))
                     VALV=(MOINTS(BEUAV)+MOINTS(BEUAV)-MOINTS(BEVAU))
     .                    *TZV(A,GA)
                     Z2(UVBEGA)=Z2(UVBEGA)+VALV
   51             CONTINUE
                  DO 60 I=FIU,LIU
                     GAI=INO(GA)+I
                     BEI=INO(BE)+I
                     GAVBEI=DOFF(MAX0(GAV,BEI))+DADD(MIN0(GAV,BEI))
                     GAIBEV=DOFF(MAX0(GAI,BEV))+DADD(MIN0(GAI,BEV))
                     VALO=(MOINTS(GAVBEI)+MOINTS(GAVBEI)-MOINTS(GAIBEV))
     .                    *TZO(I,U)
                     Z2(UVBEGA)=Z2(UVBEGA)+VALO
   60             CONTINUE
                  DO 61 I=FIV,LIV
                     GAI=INO(GA)+I
                     BEI=INO(BE)+I
                     BEUGAI=DOFF(MAX0(BEU,GAI))+DADD(MIN0(BEU,GAI))
                     BEIGAU=DOFF(MAX0(BEI,GAU))+DADD(MIN0(BEI,GAU))
                     VALO=          (MOINTS(BEUGAI)+MOINTS(BEUGAI)-
     .                    MOINTS(BEIGAU))*TZO(I,V)
                     Z2(UVBEGA)=Z2(UVBEGA)+VALO
   61             CONTINUE
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END
