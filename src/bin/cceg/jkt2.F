C
C ----------------------------------------------------------------------
C
      SUBROUTINE JKT2  (T1O,T2O,T2,CIN,DIN,EIN,FIN,AR1,AR2,
     .                  SUX1,SUX2,JAR,KAR,AUX1,
     .                  ORBSYM,FLOV,COFF,CADD,DOFF,DADD,
     .                  EOFF,EADD,FOFF,FADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,MAXOV,ITAP91)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,VALJ,VALK,VAL1,VAL2
      REAL*8  T1O(NO,NV),T2O(*),T2(*),CIN(*),DIN(*),EIN(*),FIN(*),
     .        SUX1(*),SUX2(*),JAR(*),KAR(*),AUX1(NV,NO),
     .        AR1(MAXOV,MAXOV),AR2(MAXOV,MAXOV)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        COFF(*),CADD(*),DOFF(*),DADD(*),
     .        EOFF(*),EADD(*),FOFF(*),FADD(*)
C
      CALL TRNPS(AUX1,T1O,NO,NV)
C
      ITAP91=91
      CALL RFILE(ITAP91)
      PON=1
      DO 695 TSYM=1,NIRRED
C
         XAI=0
         DO 190 I=1,NO
            IF(FZO(I).EQ.1)GO TO 190
            ISYM=ORBSYM(I)
            DO 170 A=1,NV
               IF(FZV(A).EQ.1)GO TO 170
               ASYM=ORBSYM(A+NO)
               IASYM=IEOR(ISYM,ASYM)
               IF(IASYM.NE.TSYM-1)GO TO 170
               IA=ITR(A)+I
               XAI=XAI+1
               XBJ=0
               DO 150 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 150
                  JA=ITR(A)+J
                  JSYM=ORBSYM(J)
                  BSYM=IEOR(JSYM,IASYM)+1
                  FB=FLOV(BSYM,3)-NO
                  LB=FLOV(BSYM,4)-NO
                  DO 130 B=FB,LB
                     XBJ=XBJ+1
                     IB=ITR(B)+I
                     JB=ITR(B)+J
                     IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                     IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                     AR1(XBJ,XAI)=-DIN(IAJB)-DIN(IAJB)+DIN(IBJA)
                     AR2(XBJ,XAI)=-DIN(IBJA)
  130             CONTINUE
  150          CONTINUE
  170       CONTINUE
  190    CONTINUE
C
         DO 690 U=1,NO
            IF(FZO(U).EQ.1)GO TO 690
            USYM=ORBSYM(U)
            DO 680 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 680
               BESYM=ORBSYM(BE+NO)
               BEUSYM=IEOR(BESYM,USYM)
               IF(BEUSYM.NE.TSYM-1)GO TO 680
C
               BEU=ITR(BE)+U
               XAI=0
               DO 490 I=1,NO
                  IF(FZO(I).EQ.1)GO TO 490
                  ISYM=ORBSYM(I)
                  ASYM=IEOR(ISYM,BEUSYM)+1
                  FA=FLOV(ASYM,3)-NO
                  LA=FLOV(ASYM,4)-NO
                  IU=ITR(MAX0(I,U))+MIN0(I,U)
                  EADDIU=EADD(IU)
                  DO 470 A=FA,LA
                     IA=ITR(A)+I
                     FADDIA=FADD(IA)
                     EOFFIA=EOFF(IA)
                     BEUIA=DOFF(MAX0(BEU,IA))+DADD(MIN0(BEU,IA))
                     VALJ=DIN(BEUIA)
                     BEA=ITV(MAX0(BE,A))+MIN0(BE,A)
                     FOFBEA=FOFF(BEA)
                     BEAIU=COFF(BEA)+CADD(IU)
                     VALK=CIN(BEAIU)
C
                     VAL=0.0D0
                     FB=FLOV(USYM+1,3)-NO
                     LB=FLOV(USYM+1,4)-NO
                     FB2=FB
                     LB2=LB
                     IF(LB.GT.BE)LB2=BE
                     DO 250 B=FB2,LB2
                        BEB=ITV(BE)+B
                        BEBIA=FOFF(BEB)+FADDIA
                        IB=ITR(B)+I
                        BEAIB=FOFBEA+FADD(IB)
                        VAL =VAL +FIN(BEBIA)*AUX1(B,U)
                        VALK=VALK+FIN(BEAIB)*AUX1(B,U)
  250                CONTINUE
                     FB2=FB
                     LB2=LB
                     IF(FB.LE.BE)FB2=BE+1
                     DO 251 B=FB2,LB2
                        BEB=ITV(B)+BE
                        BEBIA=FOFF(BEB)+FADDIA
                        IB=ITR(B)+I
                        BEAIB=FOFBEA+FADD(IB)
                        VAL =VAL +FIN(BEBIA)*AUX1(B,U)
                        VALK=VALK+FIN(BEAIB)*AUX1(B,U)
  251                CONTINUE
C
                     FJ=FLOV(BESYM+1,1)
                     LJ=FLOV(BESYM+1,2)
                     FJ2=FJ
                     LJ2=LJ
                     IF(LJ.GT.U)LJ2=U
                     DO 255 J=FJ2,LJ2
                        JU=ITR(U)+J
                        JUIA=EOFFIA+EADD(JU)
                        JA=ITR(A)+J
                        IUJA=EOFF(JA)+EADDIU
                        VAL =VAL -EIN(JUIA)*T1O(J,BE)
                        VALK=VALK-EIN(IUJA)*T1O(J,BE)
  255                CONTINUE
                     FJ2=FJ
                     LJ2=LJ
                     IF(FJ.LE.U)FJ2=U+1
                     DO 256 J=FJ2,LJ2
                        JU=ITR(J)+U
                        JUIA=EOFFIA+EADD(JU)
                        JA=ITR(A)+J
                        IUJA=EOFF(JA)+EADDIU
                        VAL =VAL -EIN(JUIA)*T1O(J,BE)
                        VALK=VALK-EIN(IUJA)*T1O(J,BE)
  256                CONTINUE
C
                     XAI=XAI+1
                     VALJ=VALJ+VAL
                     JAR(XAI)=VALJ+VALJ-VALK
                     KAR(XAI)=VALK
C
                     IF(BE.GE.A)THEN
                        UIBEA=TOFF(U,I,1)+TADD(BE,A)
                        T2(UIBEA)=T2(UIBEA)+VAL
                     ENDIF
                     IF(BE.LE.A)THEN
                        IUABE=TOFF(I,U,1)+TADD(A,BE)
                        T2(IUABE)=T2(IUABE)+VAL
                     ENDIF
  470             CONTINUE
  490          CONTINUE
 
               XBJ=0
               DO 580 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 580
                  JSYM=ORBSYM(J)
                  BSYM=IEOR(JSYM,BEUSYM)+1
                  FB=FLOV(BSYM,3)-NO
                  LB=FLOV(BSYM,4)-NO
                  DO 570 B=FB,LB
                     ZL=ZLX(B,BE)
                     UJBBE=TOFF(U,J,ZL)+TADD(B,BE)
                     JUBBE=TOFF(J,U,ZL)+TADD(B,BE)
                     XBJ=XBJ+1
                     SUX1(XBJ)=0.5D0*T2O(UJBBE)+T1O(U,B)*T1O(J,BE)
                     SUX2(XBJ)=SUX1(XBJ)-T2O(JUBBE)
  570             CONTINUE
  580          CONTINUE
               LIM=XBJ
C
               DO 642 AI=1,LIM
                  VALJ=0.0D0
                  VALK=0.0D0
                  DO 640 BJ=1,LIM
                     VALJ=VALJ+AR1(BJ,AI)*SUX2(BJ)
                     VALK=VALK+AR2(BJ,AI)*SUX1(BJ)
  640             CONTINUE
                  JAR(AI)=JAR(AI)+VALJ
                  KAR(AI)=KAR(AI)+VALK
  642          CONTINUE
C
               CALL WWRITW(ITAP91,JAR,INTOWP(LIM),PON,PON)
               CALL WWRITW(ITAP91,KAR,INTOWP(LIM),PON,PON)
C
  680       CONTINUE
  690    CONTINUE
  695 CONTINUE
C
      PON=1
      DO 995 TSYM=1,NIRRED
         XVG=0
         DO 790 V=1,NO
            IF(FZO(V).EQ.1)GO TO 790
            VSYM=ORBSYM(V)
            DO 770 GA=1,NV
               IF(FZV(GA).EQ.1)GO TO 770
               GASYM=ORBSYM(GA+NO)
               VGASYM=IEOR(VSYM,GASYM)
               IF(VGASYM.NE.TSYM-1)GO TO 770
               XVG=XVG+1
               XAI=0
               DO 750 I=1,NO
                  IF(FZO(I).EQ.1)GO TO 750
                  ISYM=ORBSYM(I)
                  ASYM=IEOR(ISYM,VGASYM)+1
                  FA=FLOV(ASYM,3)-NO
                  LA=FLOV(ASYM,4)-NO
                  DO 730 A=FA,LA
                     XAI=XAI+1
                     ZL=ZLX(A,GA)
                     IVAGA=TOFF(I,V,ZL)+TADD(A,GA)
                     VIAGA=TOFF(V,I,ZL)+TADD(A,GA)
                     AR1(XAI,XVG)=T2O(IVAGA)+T2O(IVAGA)-T2O(VIAGA)
                     AR2(XAI,XVG)=T2O(VIAGA)
  730             CONTINUE
  750          CONTINUE
  770       CONTINUE
  790    CONTINUE
         LIM=XVG
C
         DO 990 U=1,NO
            IF(FZO(U).EQ.1)GO TO 990
            USYM=ORBSYM(U)
            DO 980 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 980
               BESYM=ORBSYM(BE+NO)
               BEUSYM=IEOR(BESYM,USYM)
               IF(BEUSYM.NE.TSYM-1)GO TO 980
               CALL WREADW(ITAP91,JAR,INTOWP(LIM),PON,PON)
               CALL WREADW(ITAP91,KAR,INTOWP(LIM),PON,PON)
C
               DO 812 VGA=1,LIM
                  VAL1=0.0D0
                  VAL2=0.0D0
                  DO 811 AI=1,LIM
                     VAL1=VAL1+JAR(AI)*AR1(AI,VGA)
                     VAL2=VAL2+KAR(AI)*AR2(AI,VGA)
  811             CONTINUE
                  SUX1(VGA)= 0.5D0*VAL1
                  SUX2(VGA)=-0.5D0*VAL2
  812          CONTINUE
C
               VGA=0
               DO 860 V=1,NO
                  IF(FZO(V).EQ.1)GO TO 860
                  VSYM=ORBSYM(V)
                  GASYM=IEOR(VSYM,BEUSYM)+1
                  FGA=FLOV(GASYM,3)-NO
                  LGA=FLOV(GASYM,4)-NO
                  DO 850 GA=FGA,LGA
                     VGA=VGA+1
                     ZL=ZLX(BE,GA)
                     IF(BE.GE.GA)THEN
                        UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                        T2(UVBEGA)=T2(UVBEGA)+SUX1(VGA)+SUX2(VGA)
                        VUBEGA=TOFF(V,U,1)+TADD(BE,GA)
                        T2(VUBEGA)=T2(VUBEGA)+SUX2(VGA)+SUX2(VGA)
                     ENDIF
                     IF(BE.LE.GA)THEN
                        VUGABE=TOFF(V,U,1)+TADD(GA,BE)
                        T2(VUGABE)=T2(VUGABE)+SUX1(VGA)+SUX2(VGA)
                        UVGABE=TOFF(U,V,1)+TADD(GA,BE)
                        T2(UVGABE)=T2(UVGABE)+SUX2(VGA)+SUX2(VGA)
                     ENDIF
  850             CONTINUE
  860          CONTINUE
C
  980       CONTINUE
  990    CONTINUE
  995 CONTINUE
      CALL RCLOSE(ITAP91,3)
C
      RETURN
      END
