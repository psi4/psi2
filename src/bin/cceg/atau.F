C
C ----------------------------------------------------------------------
C
      SUBROUTINE ATAU  (T1O,T2O,T2,AIN,DIN,EIN,AR3,AR4,AR5,AR6,
     .                  SUX1,SFOO,SUX3,SFVV,AUX1,
     .                  ORBSYM,FLOV,AOFF,AADD,DOFF,DADD,EOFF,EADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,MAXOO,MAXVV,ITAP91)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL1,VAL2,TMP,TMV,VALP,VALM
      REAL*8  T1O(NO,NV),T2O(*),T2(*),AIN(*),DIN(*),EIN(*),
     .        SUX1(*),SFOO(*),SUX3(*),SFVV(*),AUX1(NV,NO),
     .        AR3(MAXVV,MAXOO),AR4(MAXVV,MAXOO),
     .        AR5(MAXOO,MAXVV),AR6(MAXOO,MAXVV)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        AOFF(*),AADD(*),DOFF(*),DADD(*),EOFF(*),EADD(*)
C
      CALL TRNPS(AUX1,T1O,NO,NV)
C
      CALL RFILE(ITAP91)
      PON=1
      DO 393 TSYM=1,NIRRED
         XJI=0
         DO 190 I=1,NO
            IF(FZO(I).EQ.1)GO TO 190
            ISYM=ORBSYM(I)
            DO 180 J=1,I
               IF(FZO(J).EQ.1)GO TO 180
               JSYM=ORBSYM(J)
               JISYM=IEOR(JSYM,ISYM)
               IF(JISYM.NE.TSYM-1)GO TO 180
               XJI=XJI+1
               XBA=0
               DO 170 A=1,NV
                  IF(FZV(A).EQ.1)GO TO 170
                  IA=ITR(A)+I
                  JA=ITR(A)+J
                  ASYM=ORBSYM(A+NO)
                  BSYM=IEOR(ASYM,JISYM)+1
                  FB=FLOV(BSYM,3)-NO
                  LB=FLOV(BSYM,4)-NO
                  IF(LB.GT.A)LB=A
                  DO 160 B=FB,LB
                     XBA=XBA+1
                     IB=ITR(B)+I
                     JB=ITR(B)+J
                     IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                     IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                     VAL1=DIN(IAJB)
                     VAL2=DIN(IBJA)
                     AR3(XBA,XJI)=VAL1+VAL2
                     AR4(XBA,XJI)=VAL1-VAL2
  160             CONTINUE
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
         LIMO=XJI
         LIMV=XBA
C
         DO 390 U=1,NO
            IF(FZO(U).EQ.1)GO TO 390
            USYM=ORBSYM(U)
            FA1=FLOV(USYM+1,3)-NO
            LA1=FLOV(USYM+1,4)-NO
            DO 380 V=1,U
               IF(FZO(V).EQ.1)GO TO 380
               VSYM=ORBSYM(V)
               FA2=FLOV(VSYM+1,3)-NO
               LA2=FLOV(VSYM+1,4)-NO
               VUSYM=IEOR(VSYM,USYM)
               IF(VUSYM.NE.TSYM-1)GO TO 380
               BA=0
               DO 205 A=1,NV
                  IF(FZV(A).EQ.1)GO TO 205
                  ASYM=ORBSYM(A+NO)
                  BSYM=IEOR(ASYM,VUSYM)
                  FB=FLOV(BSYM+1,3)-NO
                  LB=FLOV(BSYM+1,4)-NO
                  IF(LB.GT.A)LB=A
                  DO 204 B=FB,LB
                     BA=BA+1
                     UVAB=TOFF(U,V,1)+TADD(A,B)
                     VUAB=TOFF(V,U,1)+TADD(A,B)
                     TMP=T2O(UVAB)+T1O(U,A)*T1O(V,B)
                     TMV=T2O(VUAB)+T1O(V,A)*T1O(U,B)
                     SUX3(BA)=TMP+TMV
                     SFVV(BA)=TMP-TMV
  204             CONTINUE
                  IF(ASYM.EQ.BSYM)SUX3(BA)=SUX3(BA)*0.5D0
  205          CONTINUE
C
               JI=0
               DO 340 I=1,NO
                  IF(FZO(I).EQ.1)GO TO 340
                  ISYM=ORBSYM(I)
                  IU=ITR(MAX0(I,U))+MIN0(I,U)
                  IV=ITR(MAX0(I,V))+MIN0(I,V)
                  JSYM=IEOR(ISYM,VUSYM)
                  FJ=FLOV(JSYM+1,1)
                  LJ=FLOV(JSYM+1,2)
                  IF(LJ.GT.I)LJ=I
                  DO 330 J=FJ,LJ
                     JI=JI+1
                     JV=ITR(MAX0(J,V))+MIN0(J,V)
                     JU=ITR(MAX0(J,U))+MIN0(J,U)
                     IUJV=AOFF(MAX0(IU,JV))+AADD(MIN0(IU,JV))
                     IVJU=AOFF(MAX0(IV,JU))+AADD(MIN0(IV,JU))
                     TMP=AIN(IUJV)+AIN(IVJU)
                     TMV=AIN(IUJV)-AIN(IVJU)
                     VAL1=0.0D0
                     VAL2=0.0D0
                     DO 315 A=FA2,LA2
                        JA=ITR(A)+J
                        IA=ITR(A)+I
                        JUIA=EOFF(IA)+EADD(JU)
                        IUJA=EOFF(JA)+EADD(IU)
                        VAL1=VAL1+EIN(IUJA)*AUX1(A,V)
                        VAL2=VAL2+EIN(JUIA)*AUX1(A,V)
  315                CONTINUE
                     DO 318 A=FA1,LA1
                        IA=ITR(A)+I
                        JA=ITR(A)+J
                        IAJV=EOFF(IA)+EADD(JV)
                        JAIV=EOFF(JA)+EADD(IV)
                        VAL1=VAL1+EIN(IAJV)*AUX1(A,U)
                        VAL2=VAL2+EIN(JAIV)*AUX1(A,U)
  318                CONTINUE
                     VALP=0.0D0
                     VALM=0.0D0
                     DO 320 BA=1,LIMV
                        VALP=VALP+AR3(BA,JI)*SUX3(BA)
                        VALM=VALM+AR4(BA,JI)* SFVV(BA)
  320                CONTINUE
                     SFOO(JI)=(TMP+VAL1+VAL2+VALP)*0.5D0
                     SUX1(JI)=(TMV+VAL1-VAL2+VALM)*0.5D0
  330             CONTINUE
                  IF(ISYM.EQ.JSYM) SFOO(JI)= SFOO(JI)*0.5D0
  340          CONTINUE
               CALL WWRITW(ITAP91,SFOO,INTOWP(MAXOO),PON,PON)
               CALL WWRITW(ITAP91,SUX1,INTOWP(MAXOO),PON,PON)
  380       CONTINUE
  390    CONTINUE
  393 CONTINUE
C
      PON=1
      DO 892 TSYM=1,NIRRED
         XAB=0
         DO 490 A=1,NV
            IF(FZV(A).EQ.1)GO TO 490
            ASYM=ORBSYM(A+NO)
            DO 480 B=1,A
               IF(FZV(B).EQ.1)GO TO 480
               BSYM=ORBSYM(B+NO)
               ABSYM=IEOR(ASYM,BSYM)
               IF(ABSYM.NE.TSYM-1)GO TO 480
               XAB=XAB+1
               XIJ=0
               DO 470 I=1,NO
                  IF(FZO(I).EQ.1)GO TO 470
                  ISYM=ORBSYM(I)
                  JSYM=IEOR(ISYM,ABSYM)+1
                  FJ=FLOV(JSYM,1)
                  LJ=FLOV(JSYM,2)
                  IF(LJ.GT.I)LJ=I
                  DO 460 J=FJ,LJ
                     XIJ=XIJ+1
                     IJAB=TOFF(I,J,1)+TADD(A,B)
                     JIAB=TOFF(J,I,1)+TADD(A,B)
                     VAL1=T2O(IJAB)+T1O(I,A)*T1O(J,B)
                     VAL2=T2O(JIAB)+T1O(I,B)*T1O(J,A)
                     AR5(XIJ,XAB)=VAL1+VAL2
                     AR6(XIJ,XAB)=VAL1-VAL2
  460             CONTINUE
  470          CONTINUE
  480       CONTINUE
  490    CONTINUE
         LIMO=XIJ
         LIMV=XAB
C
         DO 891 U=1,NO
            IF(FZO(U).EQ.1)GO TO 891
            USYM=ORBSYM(U)
            FA1=FLOV(USYM+1,3)-NO
            LA1=FLOV(USYM+1,4)-NO
            DO 881 V=1,U
               IF(FZO(V).EQ.1)GO TO 881
               VSYM=ORBSYM(V)
               FA2=FLOV(VSYM+1,3)-NO
               LA2=FLOV(VSYM+1,4)-NO
               VUSYM=IEOR(VSYM,USYM)
               IF(VUSYM.NE.TSYM-1)GO TO 881
               CALL WREADW(ITAP91,SFOO,INTOWP(MAXOO),PON,PON)
               CALL WREADW(ITAP91,SUX1,INTOWP(MAXOO),PON,PON)
C
               DO 570 BEGA=1,LIMV
                  TMP=0.0D0
                  TMV=0.0D0
                  DO 550 IJ=1,LIMO
                     TMP=TMP+SFOO(IJ)* AR5(IJ,BEGA)
                     TMV=TMV+SUX1(IJ)* AR6(IJ,BEGA)
  550             CONTINUE
                  SUX3(BEGA)=TMP
                  SFVV(BEGA)=TMV
  570          CONTINUE
C
               BEGA=0
               DO 674 BE=1,NV
                  IF(FZV(BE).EQ.1)GO TO 674
                  BESYM=ORBSYM(BE+NO)
                  GASYM=IEOR(BESYM,VUSYM)+1
                  FGA=FLOV(GASYM,3)-NO
                  LGA=FLOV(GASYM,4)-NO
                  IF(LGA.GT.BE)LGA=BE
                  DO 664 GA=FGA,LGA
                     BEGA=BEGA+1
                     UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                     T2(UVBEGA)=T2(UVBEGA)+SUX3(BEGA)+SFVV(BEGA)
                     IF(U.NE.V)THEN
                        VUBEGA=TOFF(V,U,1)+TADD(BE,GA)
                        T2(VUBEGA)=T2(VUBEGA)+SUX3(BEGA)-SFVV(BEGA)
                     ENDIF
  664             CONTINUE
  674          CONTINUE
C
  881       CONTINUE
  891    CONTINUE
  892 CONTINUE
      CALL RCLOSE(ITAP91,3)
C
      RETURN
      END
