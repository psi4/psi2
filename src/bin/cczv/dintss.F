C
C-----------------------------------------------------------------------
C
      SUBROUTINE DINTSS(I1MAT,GVMAT,GOMAT,JMAT,AMAT,OMG3,
     .                 FMAT,EMAT,T1,T2,OMG1,OMG2,MOINTS,LENINT,NO,NV,
     .                 NTRO,NTRV,NO2,NV2,NTR,IOFF,JOUT,INO,INTRO,INTRV,
     .                 NOV,INV,N2OV,ZLX,UOFF,FLOV,ORBSYM,NIRRED,
     .                 DOFF,NOFF,DADD,NADD,
     .                 TOFF,TADD,NDIMT2,AR3,AR4,SUX1,SUX2,SFV,SUX3,
     .                 EMOFF,EMADD,LEMAT,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 T2(NDIMT2),T1(NO,NV),GVMAT(NV,NV),VAL,VAL1,VAL2,
     .       EMAT(LEMAT),
     .       FMAT(NV*NV),MOINTS(LENINT),VALK,
     .       I1MAT(NO,NV),GOMAT(NO,NO),JMAT(NV*NO),AMAT(NO*NO),
     .       OMG1(NO,NOV),OMG2(NO,NV),OMG3(NV,NV)
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER DOFF(NOFF),DADD(NADD),EMOFF(NO*NV),EMADD(NO*NO)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
      INTEGER ZLX(NV,NV),UOFF(NO,NO,2)
      DIMENSION FLOV(NIRRED,4),ORBSYM(NO+NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
      REAL*8 AR3(MAXVV,MAXOO),AR4(MAXVV,MAXOO),
     .       SUX1(NO*NO),SUX2(NO*NO),SFV(NV*NV),SUX3(NV*NV),
     .       TMP,TMV
C
      CALL ZERO(I1MAT,NO*NV)
      DO 710 IRA=1,NIRRED
      ASYM=IRA-1
      FA=FLOV(IRA,3)-NO
      LA=FLOV(IRA,4)-NO
      FE=FLOV(IRA,3)-NO
      LE=FLOV(IRA,4)-NO
      DO 10 A=FA,LA
         DO 720 IRI=1,NIRRED
         ISYM=IRI-1
         FI=FLOV(IRI,1)
         LI=FLOV(IRI,2)
         FM=FLOV(IRI,1)
         LM=FLOV(IRI,2)
         DO 20 I=FI,LI
            AISYM=IEOR(ASYM,ISYM)
            AI=INO(A)+I
            DO 730 IRB=1,NIRRED
            BSYM=IRB-1
            FB=FLOV(IRB,3)-NO
            LB=FLOV(IRB,4)-NO
            DO 30 B=FB,LB
               JSYM=IEOR(AISYM,BSYM)+1
               FJ=FLOV(JSYM,1)
               LJ=FLOV(JSYM,2)
               DO 40 J=FJ,LJ
                  BJ=INO(B)+J
                  AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  VAL=MOINTS(AIBJ)
                  I1MAT(I,A)=I1MAT(I,A)+(VAL+VAL)*T1(J,B)
                  I1MAT(J,A)=I1MAT(J,A)-VAL*T1(I,B)
   40          CONTINUE
   30       CONTINUE
  730       CONTINUE
            DO 80 E=FE,LE
               DO 790 IRC=1,NIRRED
               CSYM=IRC-1
               FC=FLOV(IRC,3)-NO
               LC=FLOV(IRC,4)-NO
               DO 90 C=FC,LC
                  ZL=ZLX(E,C)
                  JSYM=IEOR(AISYM,CSYM)+1
                  FJ=FLOV(JSYM,1)
                  LJ=FLOV(JSYM,2)
                  DO 100 J=FJ,LJ
                     CJ=INO(C)+J
                     IJEC=TOFF(I,J,ZL)+TADD(E,C)
                     JIEC=TOFF(J,I,ZL)+TADD(E,C)
                     AICJ=DOFF(MAX0(AI,CJ))+DADD(MIN0(AI,CJ))
                     VAL=MOINTS(AICJ)
                     GVMAT(E,A)=GVMAT(E,A)+VAL*(T2(IJEC)+T2(IJEC)
     .                          -T2(JIEC)+T1(I,E)*T1(J,C)+T1(I,E)*T1(
     .                          J,C)-T1(J,E)*T1(I,C))
  100             CONTINUE
   90          CONTINUE
  790          CONTINUE
   80       CONTINUE
            DO 131 M=FM,LM
               DO 7141 IRK=1,NIRRED
               KSYM=IRK-1
               FK=FLOV(IRK,1)
               LK=FLOV(IRK,2)
               DO 141 K=FK,LK
                  BSYM=IEOR(AISYM,KSYM)+1
                  FB=FLOV(BSYM,3)-NO
                  LB=FLOV(BSYM,4)-NO
                  DO 151 B=FB,LB
                     BK=INO(B)+K
                     AIBK=DOFF(MAX0(AI,BK))+DADD(MIN0(AI,BK))
                     VAL=MOINTS(AIBK)
                     ZL=ZLX(A,B)
                     MKAB=TOFF(M,K,ZL)+TADD(A,B)
                     KMAB=TOFF(K,M,ZL)+TADD(A,B)
                     GOMAT(I,M)=GOMAT(I,M)+VAL*(T2(MKAB)+T2(MKAB)
     .                          -T2(KMAB)+T1(M,A)*T1(K,B)+T1(M,A)*T1(
     .                          K,B)-T1(M,B)*T1(K,A))
  151             CONTINUE
  141          CONTINUE
 7141          CONTINUE
  131       CONTINUE
            DO 300 J=1,NO
               IF(FZO(J).EQ.1)GO TO 300
               JSYM=ORBSYM(J)
               KSYM=IEOR(AISYM,JSYM)+1
               FK=FLOV(KSYM,1)
               LK=FLOV(KSYM,2)
               FB=FLOV(KSYM,3)-NO
               LB=FLOV(KSYM,4)-NO
               DO 310 K=FK,LK
                  JKX=INO(J)+K
                  JKIA=EMOFF(AI)+EMADD(JKX)
                  DO 320 B=FB,LB
                     BJ=INO(B)+J
                     AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                     EMAT(JKIA)=EMAT(JKIA)+MOINTS(AIBJ)*T1(K,B)
  320             CONTINUE
  310          CONTINUE
  300       CONTINUE
   20    CONTINUE
  720    CONTINUE
   10 CONTINUE
  710 CONTINUE
C
C     A TERM
C
      CALL RFILE(TPAMAT)
      PON=1
      DO 14092 TSYM=1,NIRRED
      IDIM=NIJ(TSYM)
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
      DO 14081 V=1,U
      IF(FZO(V).EQ.1)GO TO 14081
      VSYM=ORBSYM(V)
      UVSYM=IEOR(USYM,VSYM)
      IF(UVSYM.NE.TSYM-1)GO TO 14081
C
C
          XBA=0
          DO 14025 A=1,NV
             IF(FZV(A).EQ.1)GO TO 14025
             AV=INO(A)+V
             AU=INO(A)+U
             ASYM=ORBSYM(A+NO)
             BSYM=IEOR(ASYM,UVSYM)
             FB=FLOV(BSYM+1,3)-NO
             LB=FLOV(BSYM+1,4)-NO
             IF(LB.GT. A)LB = A
             DO 14020 B=FB,LB
                BV=INO(B)+V
                BU=INO(B)+U
                AUBV=DOFF(MAX0(AU,BV))+DADD(MIN0(AU,BV))
                BUAV=DOFF(MAX0(BU,AV))+DADD(MIN0(BU,AV))
                XBA=XBA+1
                TMP=MOINTS(AUBV)+MOINTS(BUAV)
                TMV=MOINTS(AUBV)-MOINTS(BUAV)
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
      JSYM=IEOR(ISYM,UVSYM)+1
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
      CALL WREADW(TPAMAT,AMAT,INTOWP(IDIM),PON,IDUM)
      IJ=0
      DO 14768 I=1,NO
      IF(FZO(I).EQ.1)GO TO 14768
      ISYM=ORBSYM(I)
      JSYM=IEOR(ISYM,UVSYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      DO 14763 J=FJ,LJ
      IJ=IJ+1
      AMAT(IJ)=AMAT(IJ)-OMG1(I,J)
14763 CONTINUE
14768 CONTINUE
      CALL WWRITW(TPAMAT,AMAT,INTOWP(IDIM),PON,PON)
C
14081 CONTINUE
14091 CONTINUE
14092 CONTINUE
      CALL RCLOSE(TPAMAT,3)
C
C
C
      CALL RFILE(TPJMAT)
      POINT=1
      DO 130 K=1,NO
         IF(FZO(K).EQ.1)GO TO 130
         KSYM=ORBSYM(K)
         DO 134 J=1,NO
         IF(FZO(J).EQ.1)GO TO 134
         JSYM=ORBSYM(J)
            KJSYM=IEOR(KSYM,JSYM)
            CALL WREADW(TPJMAT,IDIM,1,POINT,POINT)
            CALL WREADW(TPJMAT,JMAT,INTOWP(IDIM),POINT,IDUM )
            ICNT=0
            DO 132 I=1,NO
               IF(FZO(I).EQ.1)GO TO 132
               ISYM=ORBSYM(I)
               ASYM=IEOR(KJSYM,ISYM)+1
               FA=FLOV(ASYM,3)-NO
               LA=FLOV(ASYM,4)-NO
               FB=FLOV(JSYM+1,3)-NO
               LB=FLOV(JSYM+1,4)-NO
               DO 140 A=FA,LA
                  AK=INO(A)+K
                  AI=INO(A)+I
                  ICNT=ICNT+1
                  DO 150 B=FB,LB
                     BI=INO(B)+I
                     BK=INO(B)+K
                     AIBK=DOFF(MAX0(AI,BK))+DADD(MIN0(AI,BK))
                     AKBI=DOFF(MAX0(AK,BI))+DADD(MIN0(AK,BI))
                  VAL=MOINTS(AIBK)
                  VAL1=MOINTS(AKBI)
CJ                JMAT(A,I,K,J)=JMAT(A,I,K,J)+(VAL+VAL-VAL1)*T1(J,B)
                  JMAT(ICNT   )=JMAT(ICNT   )+(VAL+VAL-VAL1)*T1(J,B)
  150             CONTINUE
  140          CONTINUE
  132       CONTINUE
      CALL WWRITW(TPJMAT,JMAT,INTOWP(IDIM),POINT,POINT)
  134       CONTINUE
  130       CONTINUE
      CALL RCLOSE(TPJMAT,3)
C
      DO 711 IRI=1,NIRRED
      ISYM=IRI-1
      FI=FLOV(IRI,1)
      LI=FLOV(IRI,2)
      DO 11 I=FI,LI
         DO 721 IRJ=1,NIRRED
         JSYM=IRJ-1
         FJ=FLOV(IRJ,1)
         LJ=FLOV(IRJ,2)
         DO 21 J=FJ,LJ
            IJSYM=IEOR(ISYM,JSYM)
            CALL ZERO(OMG1,NO*NV)
            CALL ZERO(OMG2,NO*NV)
            DO 31 B=1,NV
               IF(FZV(B).EQ.1)GO TO 31
               BSYM=ORBSYM(B+NO)
               ASYM=IEOR(IJSYM,BSYM)+1
               FA=FLOV(ASYM,3)-NO
               LA=FLOV(ASYM,4)-NO
               FL=FLOV(ASYM,1)
               LL=FLOV(ASYM,2)
               BI=INO(B)+I
               BJ=INO(B)+J
               DO 41 L=FL,LL
                  VAL=0.0D+00
                  VALK=0.0D+00
                  DO 51 A=FA,LA
                     AI=INO(A)+I
                     AJ=INO(A)+J
                     AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                     BIAJ=DOFF(MAX0(AJ,BI))+DADD(MIN0(AJ,BI))
                     VAL=VAL-MOINTS(AIBJ)*T1(L,A)
                     VALK=VALK-MOINTS(BIAJ)*T1(L,A)
   51             CONTINUE
                  OMG1(L,B)=OMG1(L,B)+VAL
                  OMG2(L,B)=OMG2(L,B)+VALK
   41          CONTINUE
   31       CONTINUE
            DO 761 IRC=1,NIRRED
            CSYM=IRC-1
            FC=FLOV(IRC,3)-NO
            LC=FLOV(IRC,4)-NO
            DO 61 C=FC,LC
               ICSYM=IEOR(ISYM,CSYM)
               JCSYM=IEOR(JSYM,CSYM)
               DO 71 L=1,NO
                  IF(FZO(L).EQ.1)GO TO 71
                  ILX=INO(I)+L
                  LSYM=ORBSYM(L)
                  KSYM=IEOR(ICSYM,LSYM)+1
                  FK=FLOV(KSYM,1)
                  LK=FLOV(KSYM,2)
                  BSYM=IEOR(JCSYM,KSYM-1)+1
                  FB=FLOV(BSYM,3)-NO
                  LB=FLOV(BSYM,4)-NO
                  DO 81 K=FK,LK
                     KCX=INO(C)+K
                     ILKC=EMOFF(KCX)+EMADD(ILX)
                     DO 91 B=FB,LB
                       ZL=ZLX(C,B)
                       JKCB=TOFF(J,K,ZL)+TADD(C,B)
                       EMAT(ILKC)=EMAT(ILKC)+OMG1(L,B)*T2(JKCB)
   91                CONTINUE
   81             CONTINUE
                  BSYM=IEOR(JCSYM,LSYM)+1
                  FB=FLOV(BSYM,3)-NO
                  LB=FLOV(BSYM,4)-NO
                  DO 181 K=FK,LK
                     KCX=INO(C)+K
                     ILKC=EMOFF(KCX)+EMADD(ILX)
                     DO 191 B=FB,LB
                       ZL=ZLX(C,B)
                       JLCB=TOFF(J,L,ZL)+TADD(C,B)
                       EMAT(ILKC)=EMAT(ILKC)+OMG2(K,B)*T2(JLCB)
  191                CONTINUE
  181             CONTINUE
   71          CONTINUE
   61       CONTINUE
  761       CONTINUE
   21    CONTINUE
  721    CONTINUE
   11 CONTINUE
  711 CONTINUE
C
C
      CALL RFILE(TPFMAT)
      POINT=1
      DO 3328 TSYM=1,NIRRED
            DO 3330 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 3330
               BESYM=ORBSYM(BE+NO)
            DO 3331 A=1,NV
               IF(FZV(A).EQ.1)GO TO 3331
               ASYM=ORBSYM(A+NO)
               BEASYM=IEOR(BESYM,ASYM)
               IF(BEASYM.NE.TSYM-1)GO TO 3331
      CALL WREADW(TPFMAT,IDIM,1,POINT,POINT)
      CALL WREADW(TPFMAT,FMAT,INTOWP(IDIM),POINT,IDUM )
      ICNT1=0
      ICNT2=0
      ICNT3=0
      ICNT4=0
      ICNT5=0
      ICNT6=0
C
            DO 332 B=1,NV
               IF(FZV(B).EQ.1)GO TO 332
               BSYM=ORBSYM(B+NO)
               ISYM=IEOR(BEASYM,BSYM)+1
               FI=FLOV(ISYM,1)
               LI=FLOV(ISYM,2)
               FJ=FLOV(ASYM+1,1)
               LJ=FLOV(ASYM+1,2)
               DO 340 I=FI,LI
                     BI=INO(B)+I
      ICNT1=ICNT1+1
                  DO 350 J=FJ,LJ
                     BEJ=INO(BE)+J
                     BIBEJ=DOFF(MAX0(BI,BEJ))+DADD(MIN0(BI,BEJ))
CF                   FMAT(BE,A,B,I)=FMAT(BE,A,B,I)+MOINTS(BIBEJ)*T1(J,A)
      FMAT(ICNT1)=FMAT(ICNT1)+MOINTS(BIBEJ)*T1(J,A)
  350             CONTINUE
  340          CONTINUE
  332       CONTINUE
C
C
         FK=FLOV(BESYM+1,1)
         LK=FLOV(BESYM+1,2)
            DO 34 B=1,NV
               IF(FZV(B).EQ.1)GO TO 34
               BSYM=ORBSYM(B+NO)
               ISYM=IEOR(BEASYM,BSYM)+1
               FI=FLOV(ISYM,1)
               LI=FLOV(ISYM,2)
               DO 44 I=FI,LI
                  ICNT6=ICNT6+1
                  ZL=ZLX(B,A)
                  VAL=0.0D+00
                  DO 54 K=FK,LK
                     IKBA=TOFF(I,K,ZL)+TADD(B,A)
                     VAL=VAL+I1MAT(K,BE)*T2(IKBA)
   54             CONTINUE
CF                FMAT(BE,A,B,I)=FMAT(BE,A,B,I)+VAL
                  FMAT(ICNT6   )=FMAT(ICNT6   )+VAL
   44          CONTINUE
   34       CONTINUE
C
      CALL WWRITW(TPFMAT,FMAT,INTOWP(IDIM),POINT,POINT)
 3331       CONTINUE
 3330       CONTINUE
 3328       CONTINUE
      CALL RCLOSE(TPFMAT,3)
C
      RETURN
      END
