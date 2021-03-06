C
C-----------------------------------------------------------------------
C
      SUBROUTINE ABTERM(Z2,Z2O,AMAT,BMAT,NO,NV,NO2,NV2,NTRO,NTRV,NTR,
     .                  IOFF,JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,ZLX,UOFF,
     .                  ZLXO,UOFFV,FLOV,ORBSYM,NIRRED,Z1,T1,
     .                  AR3,AR4,AR5,AR6,ABIN,SFV,SUX1,SUX2,SUX3,
     .                  TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
      REAL*8 Z2(NDIMT2  ),Z2O(NDIMT2  ),AMAT(NO*NO),BMAT(NV2),VAL,
     .       Z1(NO,NV),T1(NO,NV)
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER ZLX(NV,NV),UOFF(NO,NO,2)
      INTEGER ZLXO(NO,NO),UOFFV(NV,NV,2)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
      DIMENSION FLOV(NIRRED,4),ORBSYM(NO+NV)
      REAL*8 AR3(MAXVV,MAXOO),AR4(MAXVV,MAXOO),ABIN(NOV,NOV),
     .       AR5(MAXOO,MAXVV),AR6(MAXOO,MAXVV),
     .       SFV(NV*NV),SUX3(NV*NV),SUX1(NO*NO),SUX2(NO*NO),
     .       TMP,TMV,VAL1,VAL2
C
C     A TERM
C
      CALL RFILE(TPAMAT)
      POINT=1
      DO 12392 TSYM=1,NIRRED
      XAB=0
      DO 33590 A=1,NV
      IF(FZV(A).EQ.1)GO TO 33590
      ASYM=ORBSYM(A+NO)
      DO 33580 B=1,A
      IF(FZV(B).EQ.1)GO TO 33580
      BSYM=ORBSYM(B+NO)
      ABSYM=IEOR(ASYM,BSYM)
      IF(ABSYM.NE.TSYM-1)GO TO 33580
      XAB=XAB+1
      XIJ=0
      DO 33570 I=1,NO
      IF(FZO(I).EQ.1)GO TO 33570
      ISYM=ORBSYM(I)
      JSYM=IEOR(ISYM,ABSYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      IF(LJ.GT.I)LJ=I
      DO 33560 J=FJ,LJ
      XIJ=XIJ+1
      IJAB=TOFF(I,J,1)+TADD(A,B)
      JIAB=TOFF(J,I,1)+TADD(A,B)
      VAL1=Z2O(IJAB)
      VAL2=Z2O(JIAB)
      AR5(XIJ,XAB)=VAL1+VAL2
      AR6(XIJ,XAB)=VAL1-VAL2
33560 CONTINUE
33570 CONTINUE
33580 CONTINUE
33590 CONTINUE
      LIMO=XIJ
      LIMV=XAB
C
      DO 12391 U=1,NO
      IF(FZO(U).EQ.1)GO TO 12391
      USYM=ORBSYM(U)
      DO 12381 V=1,U
      IF(FZO(V).EQ.1)GO TO 12381
      VSYM=ORBSYM(V)
      VUSYM=IEOR(VSYM,USYM)
      IF(VUSYM.NE.TSYM-1)GO TO 12381
C
        IDIM=NIJ(TSYM)
        CALL WREADW(TPAMAT,AMAT,INTOWP(IDIM),POINT,POINT)
        NADD=0
        DO 15925 I=1,NO
          IF(FZO(I).EQ.1)GO TO 15925
          ISYM=ORBSYM(I)
          JSYM=IEOR(ISYM,VUSYM)
          FJ=FLOV(JSYM+1,1)
          LJ=FLOV(JSYM+1,2)
          DO 15920 J=FJ,LJ
             NADD=NADD+1
             ABIN(I,J)=AMAT(NADD)
15920        CONTINUE
15925     CONTINUE
C
          XJI=0
          DO 15025 I=1,NO
          IF(FZO(I).EQ.1)GO TO 15025
          ISYM=ORBSYM(I)
          JSYM=IEOR(ISYM,VUSYM)
          FJ=FLOV(JSYM+1,1)
          LJ=FLOV(JSYM+1,2)
          IF(LJ.GT. I)LJ = I
             DO 15020 J=FJ,LJ
             XJI=XJI+1
             TMP=ABIN(I,J)+ABIN(J,I)
             TMV=ABIN(I,J)-ABIN(J,I)
            SUX1(XJI)= TMP*0.5D0
            SUX2(XJI)= TMV*0.5D0
15020        CONTINUE
             IF(ISYM.EQ.JSYM)SUX1(XJI)=SUX1(XJI)*0.5D0
15025 CONTINUE
      DO 14370 BEGA=1,LIMV
      TMP=0.0D0
      TMV=0.0D0
      DO 14350 IJ=1,LIMO
      TMP=TMP+SUX1(IJ)* AR5(IJ,BEGA)
      TMV=TMV+SUX2(IJ)* AR6(IJ,BEGA)
14350 CONTINUE
      SUX3(BEGA)=TMP
       SFV(BEGA)=TMV
14370 CONTINUE
C
      BEGA=0
      DO 52374 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 52374
      BESYM=ORBSYM(BE+NO)
      GASYM=IEOR(BESYM,VUSYM)+1
      FGA=FLOV(GASYM,3)-NO
      LGA=FLOV(GASYM,4)-NO
      IF(LGA.GT.BE)LGA=BE
      DO 52364 GA=FGA,LGA
      BEGA=BEGA+1
      UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
      VUBEGA=TOFF(V,U,1)+TADD(BE,GA)
 
      VAL1=SUX3(BEGA)+SFV(BEGA)
      VAL2=SUX3(BEGA)-SFV(BEGA)
           Z2(UVBEGA)=Z2(UVBEGA)+VAL1
           IF(BE.EQ.GA)VAL1=VAL1*0.5D0
           Z1(U,BE)=Z1(U,BE)+VAL1*T1(V,GA)
           Z1(V,GA)=Z1(V,GA)+VAL1*T1(U,BE)
      IF(U.NE.V)THEN
           Z2(VUBEGA)=Z2(VUBEGA)+VAL2
           IF(BE.EQ.GA)VAL2=VAL2*0.5D0
           Z1(V,BE)=Z1(V,BE)+VAL2*T1(U,GA)
           Z1(U,GA)=Z1(U,GA)+VAL2*T1(V,BE)
      ENDIF
52364 CONTINUE
52374 CONTINUE
C
12381 CONTINUE
12391 CONTINUE
12392 CONTINUE
      CALL RCLOSE(TPAMAT,3)
C
C     B TERM
C
      CALL RFILE(TPBMAT)
      POINT=1
      DO 14092 TSYM=1,NIRRED
C
      XAB=0
      DO 19490 A=1,NV
      IF(FZV(A).EQ.1)GO TO 19490
      ASYM=ORBSYM(A+NO)
      DO 19480 B=1,A
      IF(FZV(B).EQ.1)GO TO 19480
      BSYM=ORBSYM(B+NO)
      BASYM=IEOR(BSYM,ASYM)
      IF(BASYM.NE.TSYM-1)GO TO 19480
      XAB=XAB+1
      XIJ=0
      DO 19470 I=1,NO
      IF(FZO(I).EQ.1)GO TO 19470
      ISYM=ORBSYM(I)
      JSYM=IEOR(ISYM,BASYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      IF(LJ.GT.I)LJ=I
      DO 19460 J=FJ,LJ
      XIJ=XIJ+1
      IJAB=TOFF(I,J,1)+TADD(A,B)
      JIAB=TOFF(J,I,1)+TADD(A,B)
      VAL1=Z2O(IJAB)
      VAL2=Z2O(JIAB)
      AR3(XAB,XIJ)=VAL1+VAL2
      AR4(XAB,XIJ)=VAL1-VAL2
19460 CONTINUE
19470 CONTINUE
19480 CONTINUE
19490 CONTINUE
      LIMO=XIJ
      LIMV=XAB
C
      DO 14091 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 14091
      BESYM=ORBSYM(BE+NO)
      DO 14081 GA=1,BE
      IF(FZV(GA).EQ.1)GO TO 14081
      GASYM=ORBSYM(GA+NO)
      BEGAS=IEOR(GASYM,BESYM)
      IF(BEGAS.NE.TSYM-1)GO TO 14081
C
        IDIM=NAB(TSYM)
        CALL WREADW(TPBMAT,BMAT,INTOWP(IDIM),POINT,POINT)
        NADD=0
        DO 14925 A=1,NV
          IF(FZV(A).EQ.1)GO TO 14925
          ASYM=ORBSYM(A+NO)
          BSYM=IEOR(ASYM,BEGAS)
          FB=FLOV(BSYM+1,3)-NO
          LB=FLOV(BSYM+1,4)-NO
          DO 14920 B=FB,LB
             NADD=NADD+1
             ABIN(A,B)=BMAT(NADD)
14920        CONTINUE
14925     CONTINUE
C
          XBA=0
          DO 14025 A=1,NV
          IF(FZV(A).EQ.1)GO TO 14025
          ASYM=ORBSYM(A+NO)
          BSYM=IEOR(ASYM,BEGAS)
          FB=FLOV(BSYM+1,3)-NO
          LB=FLOV(BSYM+1,4)-NO
          IF(LB.GT. A)LB = A
             DO 14020 B=FB,LB
             XBA=XBA+1
             TMP=ABIN(A,B)+ABIN(B,A)
             TMV=ABIN(A,B)-ABIN(B,A)
             SFV(XBA)= TMP*0.5D0
            SUX3(XBA)= TMV*0.5D0
14020        CONTINUE
             IF(ASYM.EQ.BSYM)SFV(XBA)=SFV(XBA)*0.5D0
14025 CONTINUE
      DO 14077 UV=1,LIMO
      TMP = 0.0D0
      TMV = 0.0D0
      DO 14029 AB=1,LIMV
      TMP=TMP+ SFV(AB)* AR3(AB,UV)
      TMV=TMV+SUX3(AB)* AR4(AB,UV)
14029 CONTINUE
      SUX1(UV)=TMP
      SUX2(UV)=TMV
14077 CONTINUE
C
      BEGA=IOFF(BE)+GA
      UV=0
      DO 14068 U=1,NO
      IF(FZO(U).EQ.1)GO TO 14068
      USYM=ORBSYM(U)
      VSYM=IEOR(USYM,BEGAS)+1
      SV=FLOV(VSYM,1)
      LV=FLOV(VSYM,2)
      IF(LV.GT.U )LV=U
      DO 14063 V=SV,LV
      UV=UV+1
      VAL1=SUX1(UV)+SUX2(UV)
      VAL2=SUX1(UV)-SUX2(UV)
      UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
      Z2(UVBEGA)=Z2(UVBEGA)+VAL1
                  IF(BE.EQ.GA)VAL1=VAL1*0.5D0
                  Z1(U,BE)=Z1(U,BE)+VAL1*T1(V,GA)
                  Z1(V,GA)=Z1(V,GA)+VAL1*T1(U,BE)
      IF(U.NE.V)THEN
      VUBEGA=TOFF(V,U,1)+TADD(BE,GA)
      Z2(VUBEGA)=Z2(VUBEGA)+VAL2
                  IF(BE.EQ.GA)VAL2=VAL2*0.5D0
                  Z1(V,BE)=Z1(V,BE)+VAL2*T1(U,GA)
                  Z1(U,GA)=Z1(U,GA)+VAL2*T1(V,BE)
      ENDIF
C
14063 CONTINUE
14068 CONTINUE
14081 CONTINUE
14091 CONTINUE
14092 CONTINUE
      CALL RCLOSE(TPBMAT,3)
C
      RETURN
      END
