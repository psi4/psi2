C
C-----------------------------------------------------------------------
C
      SUBROUTINE RQDENT(QMAT,RMAT,MOINTS,NO2,NV2,T1,T2,
     .                  NTRO,NTRV,NO,NV,
     .                  LENINT,NTR,IOFF,JOUT,INO,INTRO,INTRV,NOV,INV,
     .                  N2OV,ZLX,UOFF,FLOV,ORBSYM,NIRRED,
     .                  DOFF,NOFF,DADD,NADD,
     .                  TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 T2(NDIMT2),T1(NO,NV),QMAT(NV*NO),RMAT(NO*NV),MOINTS(LENINT)
     .       ,VAL,VAL2,VAL3,VALR,VALQ
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER DOFF(NOFF),DADD(NADD)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
      INTEGER ZLX(NV,NV),UOFF(NO,NO,2)
      DIMENSION FLOV(NIRRED,4),ORBSYM(NO+NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
C
      CALL RFILE(TPQMAT)
      CALL RFILE(TPRMAT)
      POINTR=1
      POINTQ=1
      DO 100 TSYM=1,NIRRED
      IDIM=NAI(TSYM)
      DO 10 BE=1,NV
         IF(FZV(BE).EQ.1)GO TO 10
         BESYM=ORBSYM(BE+NO)
         DO 20 U=1,NO
            IF(FZO(U).EQ.1)GO TO 20
            USYM=ORBSYM(U)
            BEUSYM=IEOR(BESYM,USYM)
            IF(BEUSYM.NE.TSYM-1)GO TO 20
            CALL WREADW(TPQMAT,QMAT,INTOWP(IDIM),POINTQ,IDUM )
            CALL WREADW(TPRMAT,RMAT,INTOWP(IDIM),POINTR,IDUM )
            ICNT=0
            BEU=INO(BE)+U
            DO 730 IRA=1,NIRRED
            ASYM=IRA-1
            FA=FLOV(IRA,3)-NO
            LA=FLOV(IRA,4)-NO
            ISYM=IEOR(BEUSYM,ASYM)+1
            FI=FLOV(ISYM,1)
            LI=FLOV(ISYM,2)
            DO 30 A=FA,LA
               DO 40 I=FI,LI
                  ICNT=ICNT+1
                  AI=INO(A)+I
                  BEUAI=DOFF(MAX0(BEU,AI))+DADD(MIN0(BEU,AI))
                  VAL=0.0D+00
                  VALR=0.0D+00
                  VAL2=0.0D+00
                  VAL3=0.0D+00
                  DO 750 IRJ=1,NIRRED
                  JSYM=IRJ-1
                  FJ=FLOV(IRJ,1)
                  LJ=FLOV(IRJ,2)
                  BSYM=IEOR(BEUSYM,JSYM)+1
                  FB=FLOV(BSYM,3)-NO
                  LB=FLOV(BSYM,4)-NO
                  DO 50 J=FJ,LJ
                     BEJ=INO(BE)+J
                     DO 60 B=FB,LB
                        ZL=ZLX(A,B)
                        JIAB=TOFF(J,I,ZL)+TADD(A,B)
                        IJAB=TOFF(I,J,ZL)+TADD(A,B)
                        BJ=INO(B)+J
                        BU=INO(B)+U
                        BEJBU=DOFF(MAX0(BEJ,BU))+DADD(MIN0(BEJ,BU))
                        BEUBJ=DOFF(MAX0(BEU,BJ))+DADD(MIN0(BEU,BJ))
                        VAL=VAL+MOINTS(BEUBJ)*T2(JIAB)
                        VAL2=VAL2+MOINTS(BEJBU)*T2(IJAB)
                        VAL3=VAL3+MOINTS(BEUBJ)*T2(IJAB)
                        VALR=VALR+MOINTS(BEJBU)*T2(JIAB)
   60                CONTINUE
   50             CONTINUE
  750             CONTINUE
                  VALQ=-VAL-VAL2+VAL3+VAL3+MOINTS(BEUAI)
                  QMAT(ICNT)=QMAT(ICNT)+VALQ+VALQ+VALR
                  RMAT(ICNT)=RMAT(ICNT)-VALR
                  VALQ=0.0D+00
                  VALR=0.0D+00
                  FJ=FLOV(IRA,1)
                  LJ=FLOV(IRA,2)
                  FB=FLOV(ISYM,3)-NO
                  LB=FLOV(ISYM,4)-NO
                  DO 51 J=FJ,LJ
                     BEJ=INO(BE)+J
                     DO 61 B=FB,LB
                        BU=INO(B)+U
                        BJ=INO(B)+J
                        BEUBJ=DOFF(MAX0(BEU,BJ))+DADD(MIN0(BEU,BJ))
                        BEJBU=DOFF(MAX0(BEJ,BU))+DADD(MIN0(BEJ,BU))
                        VALQ=VALQ+MOINTS(BEUBJ)*T1(J,A)*T1(I,B)
                        VALR=VALR+MOINTS(BEJBU)*T1(J,A)*T1(I,B)
   61                CONTINUE
   51             CONTINUE
                  QMAT(ICNT)=QMAT(ICNT)-VALQ-VALQ+VALR
                  RMAT(ICNT)=RMAT(ICNT)-VALR
   40          CONTINUE
   30       CONTINUE
  730       CONTINUE
            CALL WWRITW(TPQMAT,QMAT,INTOWP(IDIM),POINTQ,POINTQ)
            CALL WWRITW(TPRMAT,RMAT,INTOWP(IDIM),POINTR,POINTR)
   20    CONTINUE
   10 CONTINUE
  100 CONTINUE
      CALL RCLOSE(TPQMAT,3)
      CALL RCLOSE(TPRMAT,3)
C
      RETURN
      END
