C
C-----------------------------------------------------------------------
C
      SUBROUTINE BDENTS(BMAT,MOINTS,NO2,NV2,T1,T2,NTRO,NTRV,NO,NV,LENINT
     .                 ,NTR,IOFF,JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,FLOV,
     .                  ORBSYM,NIRRED,ZLX,UOFF,
     .                 DOFF,NOFF,DADD,NADD,
     .                 TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
      REAL*8 T1(NO,NV),BMAT(NV2),VAL,VAL2,T2(NDIMT2),MOINTS(LENINT)
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER DOFF(NOFF),DADD(NADD)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
      DIMENSION FLOV(NIRRED,4),ORBSYM(NO+NV)
      DIMENSION ZLX(NV,NV),UOFF(NO,NO,2)
    1 FORMAT (' B(',I2,',',I2,',',I2,',',I2,')=',F20.12)
C
      CALL RFILE(TPBMAT)
      POINT=1
      DO 100 TSYM=1,NIRRED
      LAB=NAB(TSYM)
      DO 10 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 10
      BESYM=ORBSYM(BE+NO)
          DO 20 GA=1,BE
          IF(FZV(GA).EQ.1)GO TO 20
          GASYM=ORBSYM(GA+NO)
          BEGAS=IEOR(BESYM,GASYM)
          BEGA=IOFF(BE)+GA
          IF(BEGAS.NE.TSYM-1)GO TO 20
      CALL WREADW(TPBMAT,BMAT,INTOWP(LAB),POINT,IDUM )
      ICNT=0
      DO 30 A=1,NV
      IF(FZV(A).EQ.1)GO TO 30
      ASYM=ORBSYM(A+NO)
      BSYM=IEOR(BEGAS,ASYM)+1
      FB=FLOV(BSYM,3)-NO
      LB=FLOV(BSYM,4)-NO
         DO 40 B=FB,LB
            AB=INV(A)+B
            ZLAB=ZLX(A,B)
            ICNT=ICNT+1
                  VAL=0.0D+00
                  DO 750 IRI=1,NIRRED
                  ISYM=IRI-1
                  JSYM=IEOR(BEGAS,ISYM)+1
                  FJ=FLOV(JSYM,1)
                  LJ=FLOV(JSYM,2)
                  FI=FLOV(IRI,1)
                  LI=FLOV(IRI,2)
                  DO 50 I=FI,LI
                     BEI=INO(BE)+I
                     DO 60 J=FJ,LJ
                        IJAB=TOFF(I,J,ZLAB)+TADD(A,B)
                        GAJ=INO(GA)+J
                        BEIGAJ=DOFF(MAX0(BEI,GAJ))+DADD(MIN0(BEI,GAJ))
                        VAL=VAL+MOINTS(BEIGAJ)*T2(IJAB)
   60                CONTINUE
   50             CONTINUE
  750             CONTINUE
CB                BMAT(AB,BEGA)=BMAT(AB,BEGA)-VAL
                  BMAT(ICNT)=BMAT(ICNT)-VAL
                  VAL=0.0D+00
                  FI=FLOV(ASYM+1,1)
                  LI=FLOV(ASYM+1,2)
                  FJ=FLOV(BSYM,1)
                  LJ=FLOV(BSYM,2)
                  DO 51 I=FI,LI
                     BEI=INO(BE)+I
                     DO 61 J=FJ,LJ
                        GAJ=INO(GA)+J
                        BEIGAJ=DOFF(MAX0(BEI,GAJ))+DADD(MIN0(BEI,GAJ))
                        VAL=VAL+MOINTS(BEIGAJ)*T1(I,A)*T1(J,B)
   61                CONTINUE
   51             CONTINUE
CB                BMAT(AB,BEGA)=BMAT(AB,BEGA)-VAL
                  BMAT(ICNT)=BMAT(ICNT)-VAL
   40          CONTINUE
   30       CONTINUE
      CALL WWRITW(TPBMAT,BMAT,INTOWP(LAB),POINT,POINT)
   20    CONTINUE
   10 CONTINUE
  100 CONTINUE
      CALL RCLOSE(TPBMAT,3)
      RETURN
      END
