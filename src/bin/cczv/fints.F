C
C-----------------------------------------------------------------------
C
      SUBROUTINE FINTS(GVMAT,EMAT,BMAT,MOINTS,NO,NV,
     .           T1,T2,LENINT,NO2,NV2,NTRO,NTRV,NTR,IOFF,JOUT,FOCK,INO,
     .           INTRO,INTRV,NOV,INV,N2OV,ZLX,UOFF,FLOV,ORBSYM,NIRRED,
     .                 FOFF,NVNV,FADD,NM,
     .                 TOFF,TADD,NDIMT2,
     .                 AR3,AR4,SUX1,SUX2,SUX3,SFV,OMG1,
     .                 EMOFF,EMADD,LEMAT,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
      REAL*8 VAL,VAL1,VAL2,VAL3,VAL4,TMP,TMV
      REAL*8 T2(NDIMT2),T1(NO,NV),GVMAT(NV,NV),
     .       EMAT(LEMAT),FOCK(NTR),
     .       BMAT(NV2),MOINTS(LENINT),
     .       AR3(MAXVV,MAXOO),AR4(MAXVV,MAXOO),
     .       SUX1(NO*NO),SUX2(NO*NO),SUX3(NV*NV),SFV(NV*NV),OMG1(NO,NO)
      INTEGER IOFF(NTR),FZO(*),FZV(*),
     .        FOFF(NVNV),FADD(NM),EMOFF(NM),EMADD(NO*NO),
     .        TOFF(NO,NO,2),TADD(NV,NV),
     .        INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV),
     .        ZLX(NV,NV),UOFF(NO,NO,2),
     .        FLOV(NIRRED,4),ORBSYM(NO+NV)
C
      DO 9 IRI=1,NIRRED
      ISYM=IRI-1
      FI=FLOV(IRI,1)
      LI=FLOV(IRI,2)
      DO 10 I=FI,LI
         FA=FLOV(ISYM+1,3)-NO
         LA=FLOV(ISYM+1,4)-NO
         DO 11 B=1,NV
            IF(FZV(B).EQ.1)GO TO 11
            BSYM=ORBSYM(B+NO)
            FC=FLOV(BSYM+1,3)-NO
            LC=FLOV(BSYM+1,4)-NO
            DO 12 A=FA,LA
               AI=INO(A)+I
               DO 13 C=FC,LC
                  BC=INV(MAX0(B,C))+MIN0(B,C)
                  AIBC=FADD(AI)+FOFF(BC)
                  VAL=MOINTS(AIBC)
                  GVMAT(C,B)=GVMAT(C,B)-VAL*(T1(I,A)+T1(I,A))
   13          CONTINUE
   12       CONTINUE
   11    CONTINUE
   10 CONTINUE
      DO 80 I=FI,LI
         FB=FLOV(ISYM+1,3)-NO
         LB=FLOV(ISYM+1,4)-NO
         DO 81 B=FB,LB
            DO 82 A=1,NV
               IF(FZV(A).EQ.1)GO TO 82
               ASYM=ORBSYM(A+NO)
               FC=FLOV(ASYM+1,3)-NO
               LC=FLOV(ASYM+1,4)-NO
               AI=INO(A)+I
               DO 83 C=FC,LC
                  BC=INV(MAX0(B,C))+MIN0(B,C)
                  AIBC=FADD(AI)+FOFF(BC)
                  VAL=MOINTS(AIBC)
                  GVMAT(C,A)=GVMAT(C,A)+VAL*T1(I,B)
   83          CONTINUE
   82       CONTINUE
   81    CONTINUE
   80 CONTINUE
    9 CONTINUE
C
C     EMAT CONTRIBUTION
C
      DO 14092 TSYM=1,NIRRED
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
      VAL1=T2(IJAB)+T1(I,A)*T1(J,B)
      VAL2=T2(JIAB)+T1(J,A)*T1(I,B)
      AR3(XAB,XIJ)=VAL1+VAL2
      AR4(XAB,XIJ)=VAL1-VAL2
19460 CONTINUE
19470 CONTINUE
19480 CONTINUE
19490 CONTINUE
      LIMO=XIJ
      LIMV=XAB
C
      DO 14091 U=1,NO
      IF(FZO(U).EQ.1)GO TO 14091
      USYM=ORBSYM(U)
      DO 14081 D=1,NV
      IF(FZV(D).EQ.1)GO TO 14081
      DSYM=ORBSYM(D+NO)
      UDSYM=IEOR(USYM,DSYM)
      IF(UDSYM.NE.TSYM-1)GO TO 14081
          XBA=0
          DO 14025 A=1,NV
          IF(FZV(A).EQ.1)GO TO 14025
             AD=INV(MAX0(A,D))+MIN0(A,D)
             AU=INO(A)+U
          ASYM=ORBSYM(A+NO)
          BSYM=IEOR(ASYM,UDSYM)
          FB=FLOV(BSYM+1,3)-NO
          LB=FLOV(BSYM+1,4)-NO
          IF(LB.GT. A)LB = A
             DO 14020 B=FB,LB
             BD=INV(MAX0(B,D))+MIN0(B,D)
             BU=INO(B)+U
             AUBD=FOFF(BD)+FADD(AU)
             BUAD=FOFF(AD)+FADD(BU)
             XBA=XBA+1
             TMP=MOINTS(AUBD)+MOINTS(BUAD)
             TMV=MOINTS(AUBD)-MOINTS(BUAD)
             SFV(XBA)= TMP*0.5D0
            SUX3(XBA)= TMV*0.5D0
14020        CONTINUE
             IF(ASYM.EQ.BSYM)SFV(XBA)=SFV(XBA)*0.5D0
14025 CONTINUE
      DO 14077 IJ=1,LIMO
      TMP = 0.0D0
      TMV = 0.0D0
      DO 14029 AB=1,LIMV
      TMP=TMP+ SFV(AB)* AR3(AB,IJ)
      TMV=TMV+SUX3(AB)* AR4(AB,IJ)
14029 CONTINUE
      SUX1(IJ)=TMP
      SUX2(IJ)=TMV
14077 CONTINUE
C
      IJ=0
      DO 14068 I=1,NO
      IF(FZO(I).EQ.1)GO TO 14068
      ISYM=ORBSYM(I)
      JSYM=IEOR(ISYM,UDSYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      IF(LJ.GT.I )LJ=I
      DO 14063 J=FJ,LJ
      IJ=IJ+1
      OMG1(I,J)=SUX1(IJ)+SUX2(IJ)
      OMG1(J,I)=SUX1(IJ)-SUX2(IJ)
14063 CONTINUE
14068 CONTINUE
C
      DO 14768 I=1,NO
      IF(FZO(I).EQ.1)GO TO 14768
      ISYM=ORBSYM(I)
      JSYM=IEOR(ISYM,UDSYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      UIX=INO(U)+I
      DO 14763 J=FJ,LJ
      JDX=INO(D)+J
      UIJD=EMOFF(JDX)+EMADD(UIX)
      EMAT(UIJD)=EMAT(UIJD)+OMG1(I,J)
14763 CONTINUE
14768 CONTINUE
C
14081 CONTINUE
14091 CONTINUE
14092 CONTINUE
C
      DO 98 A=1,NV
      IF(FZV(A).EQ.1)GO TO 98
      DO 95 B=1,NV
      IF(FZV(B).EQ.1)GO TO 95
      AB=IOFF(MAX0(A,B)+NO)+MIN0(A,B)+NO
      IF(A.EQ.B)GO TO 93
      GVMAT(A,B)=GVMAT(A,B)-FOCK(AB)
   93 CONTINUE
      VAL=0.0D0
      DO 94 I=1,NO
      IF(FZO(I).EQ.1)GO TO 94
      IB=IOFF(B+NO)+I
      VAL=VAL+FOCK(IB)*T1(I,A)
   94 CONTINUE
      GVMAT(A,B)=GVMAT(A,B)+VAL
   95 CONTINUE
   98 CONTINUE
C
      CALL RFILE(TPBMAT)
      POINT=1
      DO 200 TSYM=1,NIRRED
      LAB=NAB(TSYM)
      DO 28 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 28
      BESYM=ORBSYM(BE+NO)
         DO 29 GA=1,BE
         IF(FZV(GA).EQ.1)GO TO 29
         GASYM=ORBSYM(GA+NO)
         BEGAS=IEOR(BESYM,GASYM)
         IF(BEGAS.NE.TSYM-1)GO TO 29
         BEGA=IOFF(BE)+GA
         CALL WREADW(TPBMAT,BMAT,INTOWP(LAB),POINT,IDUM )
         ICNT=0
            DO 30 A=1,NV
               IF(FZV(A).EQ.1)GO TO 30
               ASYM=ORBSYM(A+NO)
               BEA=INV (MAX0(BE,A))+MIN0(BE,A)
               BSYM=IEOR(BEGAS,ASYM)+1
               FB=FLOV(BSYM,3)-NO
               LB=FLOV(BSYM,4)-NO
               DO 31 B=FB,LB
                  GAB=INV (MAX0(GA,B))+MIN0(GA,B)
                  AB=INV(A)+B
                  ICNT=ICNT+1
                  FI=FLOV(ASYM+1,1)
                  LI=FLOV(ASYM+1,2)
                  DO 32 I=FI,LI
                     BEI=INO(BE)+I
                     BEIGAB=FADD(BEI)+FOFF(GAB)
CB                   BMAT(AB,BEGA)
                     BMAT(ICNT)=BMAT(ICNT)+MOINTS(BEIGAB)*T1(I,A)
   32             CONTINUE
                  FI=FLOV(BSYM,1)
                  LI=FLOV(BSYM,2)
                  DO 332 I=FI,LI
                     GAI=INO(GA)+I
                     GAIBEA=FADD (GAI)+FOFF(BEA)
CB                   BMAT(AB,BEGA)
                     BMAT(ICNT)=BMAT(ICNT)+MOINTS(GAIBEA)*T1(I,B)
  332             CONTINUE
   31          CONTINUE
   30       CONTINUE
      CALL WWRITW(TPBMAT,BMAT,INTOWP(LAB),POINT,POINT)
   29    CONTINUE
   28 CONTINUE
  200 CONTINUE
      CALL RCLOSE(TPBMAT,3)
C
      RETURN
      END
