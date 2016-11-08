C
C-----------------------------------------------------------------------
C
      SUBROUTINE NFINTS(FMAT,MOINTS,T1,T2,NO2,NV2,NTRO,NTRV,NO,NV,
     .           LENINT,NTR,IOFF,JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,ZLX,
     .           UOFF,FLOV,ORBSYM,NIRRED,
     .           FOFF,NVNV,FADD,NM,
     .           TOFF,TADD,NDIMT2,DINTS,DOFF,DADD,
     .           AR1,AR2,JAR,KAR,SUX1,SUX2,AUX1,
     .           BUF,IBUF,LENBUF,NSHOV,NSLVO,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
      REAL*8 T2(NDIMT2),T1(NO,NV),VAL,VAL1,VAL2,VAL3,MOINTS(LENINT),
     .       VAL4,FMAT(NV*NV),VALL,VALJ,VALK,DINTS(LENINT)
      REAL*8 AR1(MAXOV,MAXOV),AR2(MAXOV,MAXOV),
     .       JAR(NO*NV),KAR(NO*NV),
     .       SUX1(NO*NV),SUX2(NO*NV),AUX1(NO,NO),BUF(*)
C     REAL*8 XMAT(9,9,5)
      INTEGER IOFF(NTR),FZO(*),FZV(*),NLIM(8),IBUF(*)
      INTEGER FOFF(NVNV),FADD(NM),DOFF(NM),DADD(NM)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER ZLX(NV,NV),UOFF(NO,NO,2)
      DIMENSION FLOV(NIRRED,4),ORBSYM(NO+NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
C
C     NEW F INTERMEDIATE
C
C     MODIFIED NOV 15, 88:
C     MOINTS CONTAINS F INTS. AFTER 1ST TERM IS STORE IN 89, D INTS
C     ARE READ IN MOINTS AND USED. F INTS ARE RESTORED IN MOINTS BEFORE
C     RETURNING TO MAIN.
C
C
      CALL RFILE(TPFMAT)
      itap89=89
      CALL RFILE(itap89)
C
      PON=1
      DO 14695 TSYM=1,NIRRED
      DO 14690 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 14690
      BESYM=ORBSYM(BE+NO)
      DO 14680 A=1,NV
      IF(FZV(A).EQ.1)GO TO 14680
      ASYM=ORBSYM(A+NO)
      BEASYM=IEOR(BESYM,ASYM)
      IF(BEASYM.NE.TSYM-1)GO TO 14680
      BEA=INV(MAX0(BE,A))+MIN0(BE,A)
C
      XCJ=0
      DO 14490 J=1,NO
      IF(FZO(J).EQ.1)GO TO 14490
      BEJ=INO(BE)+J
      JSYM=ORBSYM(J)
      CSYM=IEOR(JSYM,BEASYM)+1
      FC=FLOV(CSYM,3)-NO
      LC=FLOV(CSYM,4)-NO
      DO 14470 C=FC,LC
      CJ=INO(C)+J
      BEACJ=FOFF(BEA)+FADD(CJ)
      VALJ=MOINTS(BEACJ)
      CA=INV(MAX0(C,A))+MIN0(C,A)
      BEJCA=FOFF(CA)+FADD(BEJ)
      VALK=MOINTS(BEJCA)
C
C     FK=FLOV(ASYM+1,1)
C     LK=FLOV(ASYM+1,2)
C     DO 15450 K=FK,LK
C     BEK=INO(BE)+K
C     BEKCJ=DOFF(MAX0(BEK,CJ))+DADD(MIN0(BEK,CJ))
C     CK=INO(C)+K
C     BEJCK=DOFF(MAX0(BEJ,CK))+DADD(MIN0(BEJ,CK))
C     VALJ=VALJ-DINTS(BEKCJ)*T1(K,A)
C     VALK=VALK-DINTS(BEJCK)*T1(K,A)
14450 CONTINUE
C
      XCJ=XCJ+1
      JAR(XCJ)=VALJ+VALJ-VALK
      KAR(XCJ)=VALK
14470 CONTINUE
14490 CONTINUE
      LIM=XCJ
C
      CALL WWRITW(itap89,JAR,INTOWP(LIM),PON,PON)
      CALL WWRITW(itap89,KAR,INTOWP(LIM),PON,PON)
14680 CONTINUE
14690 CONTINUE
      NLIM(TSYM)=LIM
14695 CONTINUE
C
      CALL RDINTS(63,DINTS,NSHOV,
     .            BUF,IBUF,LENBUF,NO,NV,
     .            INO,INO,DOFF,DADD)
C
      PON=1
      DO 15695 TSYM=1,NIRRED
      LIM=NLIM(TSYM)
      DO 15690 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 15690
      BESYM=ORBSYM(BE+NO)
      DO 15680 A=1,NV
      IF(FZV(A).EQ.1)GO TO 15680
      ASYM=ORBSYM(A+NO)
      BEASYM=IEOR(BESYM,ASYM)
      IF(BEASYM.NE.TSYM-1)GO TO 15680
      CALL WREADW(itap89,JAR,INTOWP(LIM),PON,CON)
      CALL WREADW(itap89,KAR,INTOWP(LIM),CON,DUM)
C
      XCJ=0
      DO 15490 J=1,NO
      IF(FZO(J).EQ.1)GO TO 15490
      BEJ=INO(BE)+J
      JSYM=ORBSYM(J)
      CSYM=IEOR(JSYM,BEASYM)+1
      FC=FLOV(CSYM,3)-NO
      LC=FLOV(CSYM,4)-NO
      DO 15470 C=FC,LC
      CJ=INO(C)+J
      VALJ=0.0D0
      VALK=0.0D0
C
      FK=FLOV(ASYM+1,1)
      LK=FLOV(ASYM+1,2)
      DO 15450 K=FK,LK
      BEK=INO(BE)+K
      BEKCJ=DOFF(MAX0(BEK,CJ))+DADD(MIN0(BEK,CJ))
      CK=INO(C)+K
      BEJCK=DOFF(MAX0(BEJ,CK))+DADD(MIN0(BEJ,CK))
      VALJ=VALJ-DINTS(BEKCJ)*T1(K,A)
      VALK=VALK-DINTS(BEJCK)*T1(K,A)
15450 CONTINUE
C
      XCJ=XCJ+1
      JAR(XCJ)=JAR(XCJ)+VALJ+VALJ-VALK
      KAR(XCJ)=KAR(XCJ)+VALK
15470 CONTINUE
15490 CONTINUE
C
      CALL WWRITW(itap89,JAR,INTOWP(LIM),PON,PON)
      CALL WWRITW(itap89,KAR,INTOWP(LIM),PON,PON)
15680 CONTINUE
15690 CONTINUE
15695 CONTINUE
C
C   GO FOR THE EXCHANGE TERM
C
      POINT=1
      DO 985 TSYM=1,NIRRED
         DO 980 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 980
            BESYM=ORBSYM(BE+NO)
            DO 940 A=1,NV
               IF(FZV(A).EQ.1)GO TO 940
               ASYM=ORBSYM(A+NO)
               BEASYM=IEOR(BESYM,ASYM)
               IF(BEASYM.NE.TSYM-1)GO TO 940
               CALL WREADW(TPFMAT,IDIM,1,POINT,POINT)
               CALL WREADW(TPFMAT,FMAT,INTOWP(IDIM),POINT,IDUM )
               ICNT3=0
C              ICNT=0
C              DO 908 B=1,NV
C                 BSYM=ORBSYM(B+NO)
C                 ISYM=IEOR(BEASYM,BSYM)+1
C                 FI=FLOV(ISYM,1)
C                 LI=FLOV(ISYM,2)
C                 DO 906 I=FI,LI
C                    ICNT=ICNT+1
C                    VAL=0.0D0
C                    DO 904 C=1,NV
C                       CSYM=ORBSYM(C+NO)
C                       ZL=ZLX(C,A)
C                       BCSYM=IEOR(BSYM,CSYM)
C                       BC=INV(MAX0(B,C))+MIN0(B,C)
C                       JSYM=IEOR(BESYM,BCSYM)+1
C                       FJ=FLOV(JSYM,1)
C                       LJ=FLOV(JSYM,2)
C                       DO 902 J=FJ,LJ
C                          BEJ=INO(BE)+J
C                          BEJBC=FADD(BEJ)+FOFF(BC)
C                          IJCA=TOFF(I,J,ZL)+TADD(C,A)
C                          VAL=VAL+MOINTS(BEJBC)*T2(IJCA)
C 902                   CONTINUE
C 904                CONTINUE
C                    FMAT(ICNT)=FMAT(ICNT)+VAL
CF                   FMAT(BE,A,B,I)=FMAT(BE,A,B,I)+VAL
C 906             CONTINUE
C 908          CONTINUE
               CALL ZERO(AUX1,NO*NO)
               DO 918 I=1,NO
                  IF(FZO(I).EQ.1)GO TO 918
                  ISYM=ORBSYM(I)
                  JSYM=IEOR(BEASYM,ISYM)+1
                  FJ=FLOV(JSYM,1)
                  LJ=FLOV(JSYM,2)
                  DO 916 J=FJ,LJ
                     DO 914 C=1,NV
                        IF(FZV(C).EQ.1)GO TO 914
                        CJ=INO(C)+J
                        ZLA=ZLX(C,A)
                        CSYM=ORBSYM(C+NO)
                        CASYM=IEOR(CSYM,ASYM)
                        KSYM=IEOR(ISYM,CASYM)+1
                        FK=FLOV(KSYM,1)
                        LK=FLOV(KSYM,2)
                        DO 912 K=FK,LK
                           BEK=INO(BE)+K
                           IKCA=TOFF(I,K,ZLA)+TADD(C,A)
                           IK=UOFF(I,K,ZLA)
                           BEKCJ=DOFF(MAX0(BEK,CJ))+DADD(MIN0(BEK,CJ))
                           AUX1(J,I)=AUX1(J,I)-DINTS(BEKCJ)*T2(IKCA)
  912                   CONTINUE
  914                CONTINUE
  916             CONTINUE
  918          CONTINUE
               DO 926 B=1,NV
                  IF(FZV(B).EQ.1)GO TO 926
                  BSYM=ORBSYM(B+NO)
                  FJ=FLOV(BSYM+1,1)
                  LJ=FLOV(BSYM+1,2)
                  ISYM=IEOR(BEASYM,BSYM)+1
                  FI=FLOV(ISYM,1)
                  LI=FLOV(ISYM,2)
                  DO 924 I=FI,LI
                     ICNT3=ICNT3+1
                     DO 922 J=FJ,LJ
CF                      FMAT(BE,A,B,I)=FMAT(BE,A,B,I)+OMG1(J,I)*T1(J,B)
                        FMAT(ICNT3   )=FMAT(ICNT3   )+AUX1(J,I)*T1(J,B)
  922                CONTINUE
  924             CONTINUE
  926          CONTINUE
               CALL WWRITW(TPFMAT,FMAT,INTOWP(IDIM),POINT,POINT)
  940       CONTINUE
  980    CONTINUE
  985 CONTINUE
C
      CALL RDINTS(65,MOINTS,NSLVO,
     .            BUF,IBUF,LENBUF,NO,NV,
     .            INV,INO,FOFF,FADD)
C
      POINT=1
      DO 885 TSYM=1,NIRRED
         DO 880 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 880
            BESYM=ORBSYM(BE+NO)
            DO 840 A=1,NV
               IF(FZV(A).EQ.1)GO TO 840
               ASYM=ORBSYM(A+NO)
               BEASYM=IEOR(BESYM,ASYM)
               IF(BEASYM.NE.TSYM-1)GO TO 840
               CALL WREADW(TPFMAT,IDIM,1,POINT,POINT)
               CALL WREADW(TPFMAT,FMAT,INTOWP(IDIM),POINT,IDUM )
               ICNT=0
               DO 808 B=1,NV
                  IF(FZV(B).EQ.1)GO TO 808
                  BSYM=ORBSYM(B+NO)
                  ISYM=IEOR(BEASYM,BSYM)+1
                  FI=FLOV(ISYM,1)
                  LI=FLOV(ISYM,2)
                  DO 806 I=FI,LI
                     ICNT=ICNT+1
                     VAL=0.0D0
                     DO 804 C=1,NV
                        IF(FZV(C).EQ.1)GO TO 804
                        CSYM=ORBSYM(C+NO)
                        ZL=ZLX(C,A)
                        BCSYM=IEOR(BSYM,CSYM)
                        BC=INV(MAX0(B,C))+MIN0(B,C)
                        JSYM=IEOR(BESYM,BCSYM)+1
                        FJ=FLOV(JSYM,1)
                        LJ=FLOV(JSYM,2)
                        DO 802 J=FJ,LJ
                           BEJ=INO(BE)+J
                           BEJBC=FADD(BEJ)+FOFF(BC)
                           IJCA=TOFF(I,J,ZL)+TADD(C,A)
                           VAL=VAL+MOINTS(BEJBC)*T2(IJCA)
  802                   CONTINUE
  804                CONTINUE
                     FMAT(ICNT)=FMAT(ICNT)+VAL
  806             CONTINUE
  808          CONTINUE
               CALL WWRITW(TPFMAT,FMAT,INTOWP(IDIM),POINT,POINT)
  840       CONTINUE
  880    CONTINUE
  885 CONTINUE
C
C
C
      PON=1
      POINT=1
      DO 16695 TSYM=1,NIRRED
      LIM=NLIM(TSYM)
      XIB=0
      DO 16790 B=1,NV
      IF(FZV(B).EQ.1)GO TO 16790
      BSYM=ORBSYM(B+NO)
      DO 16770 I=1,NO
      IF(FZO(I).EQ.1)GO TO 16770
      ISYM=ORBSYM(I)
      IBSYM=IEOR(ISYM,BSYM)
      IF(IBSYM.NE.TSYM-1)GO TO 16770
      XIB=XIB+1
      XCJ=0
      DO 16750 J=1,NO
      IF(FZO(J).EQ.1)GO TO 16750
      JSYM=ORBSYM(J)
      CSYM=IEOR(JSYM,IBSYM)+1
      FC=FLOV(CSYM,3)-NO
      LC=FLOV(CSYM,4)-NO
      DO 16730 C=FC,LC
      XCJ=XCJ+1
      ZL=ZLX(C,B)
      JICB=TOFF(J,I,ZL)+TADD(C,B)
      IJCB=TOFF(I,J,ZL)+TADD(C,B)
      AR1(XCJ,XIB)=T2(JICB)+T2(JICB)-T2(IJCB)
      AR2(XCJ,XIB)=T2(IJCB)
16730 CONTINUE
16750 CONTINUE
16770 CONTINUE
16790 CONTINUE
C
      DO 16690 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 16690
      BESYM=ORBSYM(BE+NO)
      DO 16680 A=1,NV
      IF(FZV(A).EQ.1)GO TO 16680
      ASYM=ORBSYM(A+NO)
      BEASYM=IEOR(BESYM,ASYM)
      IF(BEASYM.NE.TSYM-1)GO TO 16680
      CALL WREADW(TPFMAT,IDIM,1,POINT,POINT)
      CALL WREADW(TPFMAT,FMAT,INTOWP(IDIM),POINT,DUM)
      ICNT=0
      BEA=INV(MAX0(BE,A))+MIN0(BE,A)
      CALL WREADW(itap89,JAR,INTOWP(LIM),PON,PON)
      CALL WREADW(itap89,KAR,INTOWP(LIM),PON,PON)
C
      DO 16661 BI=1,LIM
      VAL1=0.0D0
      VAL2=0.0D0
      DO 16391 CJ=1,LIM
      VAL1=VAL1+JAR(CJ)*AR1(CJ,BI)
      VAL2=VAL2+KAR(CJ)*AR2(CJ,BI)
16391 CONTINUE
      SUX1(BI)=-0.5D0*VAL1
      SUX2(BI)= 0.5D0*VAL2
16661 CONTINUE
C
      XBI=0
      DO 16660 B=1,NV
      IF(FZV(B).EQ.1)GO TO 16660
      BSYM=ORBSYM(B+NO)
      ISYM=IEOR(BSYM,BEASYM)+1
      FI=FLOV(ISYM,1)
      LI=FLOV(ISYM,2)
      DO 16410  I=FI,LI
      BI=INO(B)+I
      BIBEA=FOFF(BEA)+FADD(BI)
      XBI=XBI+1
      ICNT=ICNT+1
      FMAT(ICNT)=FMAT(ICNT)+SUX1(XBI)+SUX2(XBI)-MOINTS(BIBEA)
16410 CONTINUE
16660 CONTINUE
      CALL WWRITW(TPFMAT,FMAT,INTOWP(IDIM),POINT,POINT)
16680 CONTINUE
16690 CONTINUE
16695 CONTINUE
C
      CALL RCLOSE(itap89,4)
      CALL RCLOSE(TPFMAT,3)
C
      RETURN
      END
