C
C-----------------------------------------------------------------------
C
      SUBROUTINE EINTS(JMAT,GOMAT,EMAT,FMAT,AMAT,MOINTS,NO,NV,T1,T2,
     .                 LENINT,NO2,NV2,NTRO,NTRV,NTR,IOFF,JOUT,FOCK,INO,
     .                 INTRO,INTRV,NOV,INV,N2OV,ZLX,UOFF,FLOV,ORBSYM,
     .                 NIRRED,EOFF,NOFF,EADD,NADD,
     .                 TOFF,TADD,NDIMT2,EMOFF,EMADD,LEMAT,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 VAL,VAL2,VAL3,VAL4
      REAL*8 T2(NDIMT2),T1(NO,NV),JMAT(NV*NO),GOMAT(NO,NO),
     .       EMAT(LEMAT),VAL1,FMAT(NV*NV),AMAT(NO*NO),
     .       MOINTS(LENINT),FOCK(NTR)
      INTEGER IOFF(NTR),FZO(*),FZV(*),
     .        INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV),
     .        EOFF(NOFF),EADD(NADD),EMOFF(NO*NV),EMADD(NO*NO),
     .        TOFF(NO,NO,2),TADD(NV,NV),
     .        ZLX(NV,NV),UOFF(NO,NO,2),
     .        FLOV(NIRRED,4),ORBSYM(NO+NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
C
      DO 9 IRA=1,NIRRED
         ASYM=IRA-1
         FA=FLOV(IRA,3)-NO
         LA=FLOV(IRA,4)-NO
         DO 10 A=FA,LA
            DO 11 J=1,NO
               IF(FZO(J).EQ.1)GO TO 11
               JSYM=ORBSYM(J)
               AJSYM=IEOR(ASYM,JSYM)
               DO 12 I=1,NO
                  IF(FZO(I).EQ.1)GO TO 12
                  ISYM=ORBSYM(I)
                  AI=INO(A)+I
                  KSYM=IEOR(AJSYM,ISYM)+1
                  FK=FLOV(KSYM,1)
                  LK=FLOV(KSYM,2)
                  DO 13 K=FK,LK
                     JK=INO (MAX0(J,K))+MIN0(J,K)
                     JKX=INO(J)+K
                     AIJK=EOFF(AI)+EADD(JK)
                     JKIA=EMOFF(AI)+EMADD(JKX)
                     VAL=MOINTS(AIJK)
                     GOMAT(J,K)=GOMAT(J,K)+VAL*(T1(I,A)+T1(I,A))
                     GOMAT(I,K)=GOMAT(I,K)-VAL*T1(J,A)
                     EMAT(JKIA)=EMAT(JKIA)+VAL
   13             CONTINUE
   12          CONTINUE
               DO 14 I=1,NO
                  IF(FZO(I).EQ.1)GO TO 14
                  ISYM=ORBSYM(I)
                  AI=INO(A)+I
                  JI=INO (MAX0(J,I))+MIN0(J,I)
                  JIX=INO(J)+I
                  LSYM=IEOR(AJSYM,ISYM)+1
                  ALSYM=IEOR(ASYM,LSYM-1)
                  FL=FLOV(LSYM,1)
                  LL=FLOV(LSYM,2)
                  DO 15 L=FL,LL
                     JL=INO (MAX0(J,L))+MIN0(J,L)
                     AL=INO(A)+L
                     AIJL=EOFF (AI)+EADD(JL)
                     ALJI=EOFF (AL)+EADD(JI)
                     VAL3=MOINTS(AIJL)
                     VAL4=MOINTS(ALJI)
                     DO 16 B=1,NV
                        IF(FZV(B).EQ.1)GO TO 16
                        BSYM=ORBSYM(B+NO)
                        ZL=ZLX(A,B)
                        JBX=INO(B)+J
                        KSYM=IEOR(ALSYM,BSYM)+1
                        FK=FLOV(KSYM,1)
                        LK=FLOV(KSYM,2)
                        DO 17 K=FK,LK
                           KLAB=TOFF(K,L,ZL)+TADD(A,B)
                           KBX=INO(B)+K
                           IKX=INO(I)+K
                           IKJB=EMOFF(JBX)+EMADD(IKX)
                           JIKB=EMOFF(KBX)+EMADD(JIX)
                           EMAT(IKJB)=EMAT(IKJB)-VAL3*T2(KLAB)
                           EMAT(JIKB)=EMAT(JIKB)-VAL4*T2(KLAB)
   17                   CONTINUE
   16                CONTINUE
   15             CONTINUE
   14          CONTINUE
   11       CONTINUE
   10    CONTINUE
    9 CONTINUE
C
      CALL RFILE(TPJMAT)
      POINT=1
      DO 44 J=1,NO
         IF(FZO(J).EQ.1)GO TO 44
         JSYM=ORBSYM(J)
         DO 43 K=1,NO
            IF(FZO(K).EQ.1)GO TO 43
            ICNT=0
            KSYM=ORBSYM(K)
            JKSYM=IEOR(KSYM,JSYM)
            JK=INO (MAX0(J,K))+MIN0(J,K)
            DO 42 I=1,NO
               IF(FZO(I).EQ.1)GO TO 42
               ISYM=ORBSYM(I)
               IK=INO (MAX0(I,K))+MIN0(I,K)
               ASYM=IEOR(JKSYM,ISYM)+1
               FA=FLOV(ASYM,3)-NO
               LA=FLOV(ASYM,4)-NO
               DO 41 A=FA,LA
                  ICNT=ICNT+1
                  AJ=INO(A)+J
                  AI=INO(A)+I
                  AIJK=EOFF(AI)+EADD(JK)
                  AJIK=EOFF(AJ)+EADD(IK)
                  VAL=MOINTS(AIJK)
                  VAL1=MOINTS(AJIK)
CJ                JMAT(A,I,J,K)
                  JMAT(ICNT)=VAL+VAL-VAL1
   41          CONTINUE
   42       CONTINUE
            CALL WWRITW(TPJMAT,ICNT,1,POINT,POINT)
            CALL WWRITW(TPJMAT,JMAT,INTOWP(ICNT),POINT,POINT)
   43    CONTINUE
   44 CONTINUE
      CALL RCLOSE(TPJMAT,3)
C
      CALL RFILE(TPFMAT)
      POINT=1
      DO 267 TSYM=1,NIRRED
         DO 26 A=1,NV
            IF(FZV(A).EQ.1)GO TO 26
            ASYM=ORBSYM(A+NO)
            DO 25 C=1,NV
               IF(FZV(C).EQ.1)GO TO 25
               CSYM=ORBSYM(C+NO)
               ACSYM=IEOR(ASYM,CSYM)
               IF(ACSYM.NE.TSYM-1)GO TO 25
               CALL WREADW(TPFMAT,IDIM,1,POINT,POINT)
               CALL WREADW(TPFMAT,FMAT,INTOWP(IDIM),POINT,IDUM)
               ICNT=0
               DO 21 B=1,NV
                  IF(FZV(B).EQ.1)GO TO 21
                  BSYM=ORBSYM(B+NO)
                  BCSYM=IEOR(BSYM,CSYM)
                  ZL=ZLX(B,C)
                  JSYM=IEOR(ACSYM,BSYM)+1
                  FJ=FLOV(JSYM,1)
                  LJ=FLOV(JSYM,2)
                  DO 22 J=FJ,LJ
                     VAL=0.0D+00
                     ICNT=ICNT+1
                     DO 823 IRI=1,NIRRED
                        ISYM=IRI-1
                        FI=FLOV(IRI,1)
                        LI=FLOV(IRI,2)
                        DO 23 I=FI,LI
                           AI=INO(A)+I
                           KSYM=IEOR(ISYM,BCSYM)+1
                           FK=FLOV(KSYM,1)
                           LK=FLOV(KSYM,2)
                           DO 24 K=FK,LK
                              KIBC=TOFF(K,I,ZL)+TADD(B,C)
                              JK=INO (MAX0(J,K))+MIN0(J,K)
                              AIJK=EOFF (AI)+EADD(JK)
                              VAL=VAL-MOINTS(AIJK)*T2(KIBC)
   24                      CONTINUE
   23                   CONTINUE
  823                CONTINUE
CF                   FMAT(A,C,B,J)
                     FMAT(ICNT)=FMAT(ICNT)+VAL
                     FK=FLOV(BSYM+1,1)
                     LK=FLOV(BSYM+1,2)
                     FI=FLOV(CSYM+1,1)
                     LI=FLOV(CSYM+1,2)
                     VAL=0.0D+00
                     DO 723 I=FI,LI
                        AI=INO(A)+I
                        DO 724 K=FK,LK
                           KIRC=UOFF(K,I,ZL)
                           JK=INO (MAX0(J,K))+MIN0(J,K)
                           AIJK=EOFF (AI)+EADD(JK)
                           VAL=VAL-MOINTS(AIJK)*T1(K,B)*T1(I,C)
  724                   CONTINUE
  723                CONTINUE
CF                   FMAT(A,C,B,J)
                     FMAT(ICNT)=FMAT(ICNT)+VAL
   22             CONTINUE
   21          CONTINUE
               CALL WWRITW(TPFMAT,FMAT,INTOWP(IDIM),POINT,POINT)
   25       CONTINUE
   26    CONTINUE
  267 CONTINUE
      CALL RCLOSE(TPFMAT,3)
C
      DO 98 I=1,NO
         IF(FZO(I).EQ.1)GO TO 98
         DO 95 J=1,NO
            IF(FZO(J).EQ.1)GO TO 95
            IJ=IOFF(MAX0(I,J))+MIN0(I,J)
            IF(I.NE.J)GOMAT(I,J)=GOMAT(I,J)+FOCK(IJ)
            VAL=0.0D0
            DO 94 A=1,NV
               IF(FZV(A).EQ.1)GO TO 94
               IA=IOFF(A+NO)+I
               VAL=VAL+FOCK(IA)*T1(J,A)
   94       CONTINUE
            GOMAT(I,J)=GOMAT(I,J)+VAL
   95    CONTINUE
   98 CONTINUE
C
 
      CALL RFILE(TPAMAT)
      PON=1
      DO 200 TSYM=1,NIRRED
         LIJ=NIJ(TSYM)
         DO 28 U=1,NO
            IF(FZO(U).EQ.1)GO TO 28
            USYM=ORBSYM(U)
            DO 29 V=1,U
               IF(FZO(V).EQ.1)GO TO 29
               VSYM=ORBSYM(V)
               UVSYM=IEOR(USYM,VSYM)
               IF(UVSYM.NE.TSYM-1)GO TO 29
               CALL WREADW(TPAMAT,AMAT,INTOWP(LIJ),PON,IDUM)
               ICNT=0
               DO 30 I=1,NO
                  IF(FZO(I).EQ.1)GO TO 30
                  ISYM=ORBSYM(I)
                  UI=INO (MAX0(U,I))+MIN0(U,I)
                  JSYM=IEOR(UVSYM,ISYM)+1
                  FJ=FLOV(JSYM,1)
                  LJ=FLOV(JSYM,2)
                  DO 31 J=FJ,LJ
                     VJ=INO (MAX0(V,J))+MIN0(V,J)
                     FA=FLOV(ISYM+1,3)-NO
                     LA=FLOV(ISYM+1,4)-NO
                     ICNT=ICNT+1
                     DO 32 A=FA,LA
                        AU=INO(A)+U
                        AUVJ=EOFF (AU)+EADD(VJ)
CA                      AMAT(UV,IJ)
                        AMAT(ICNT )=AMAT(ICNT )-MOINTS(AUVJ)*T1(I,A)
   32                CONTINUE
                     FA=FLOV(JSYM,3)-NO
                     LA=FLOV(JSYM,4)-NO
                     DO 33 A=FA,LA
                        AV=INO(A)+V
                        AVUI=EOFF (AV)+EADD(UI)
CA                      AMAT(UV,IJ)
                        AMAT(ICNT )=AMAT(ICNT )-MOINTS(AVUI)*T1(J,A)
   33                CONTINUE
   31             CONTINUE
   30          CONTINUE
               CALL WWRITW(TPAMAT,AMAT,INTOWP(LIJ),PON,PON)
   29       CONTINUE
   28    CONTINUE
  200 CONTINUE
      CALL RCLOSE(TPAMAT,3)
C
      RETURN
      END
