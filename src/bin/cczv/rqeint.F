C
C-----------------------------------------------------------------------
C
      SUBROUTINE RQEINT(RMAT,QMAT,MOINTS,T1,
     .                  LENINT,NO,NV,NO2,NV2,NTRO,NTRV,
     .                  NTR,IOFF,JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,FLOV,
     .                  ORBSYM,NIRRED,EOFF,NOFF,EADD,NADD,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 VALR,VALQ
      REAL*8 T1(NO,NV),RMAT(NV*NO),QMAT(NV*NO),MOINTS(LENINT)
      INTEGER IOFF(NTR),FZO(*),FZV(*),
     .        INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV),
     .        EOFF(NOFF),EADD(NADD),
     .        FLOV(NIRRED,4),ORBSYM(NO+NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
C
      CALL RFILE(TPRMAT)
      CALL RFILE(TPQMAT)
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
            CALL WREADW(TPRMAT,RMAT,INTOWP(IDIM),POINTR,IDUM )
            ICNT=0
            BEU=INO(BE)+U
            DO 30 A=1,NV
               IF(FZV(A).EQ.1)GO TO 30
               ASYM=ORBSYM(A+NO)
               ISYM=IEOR(BEUSYM,ASYM)+1
               FI=FLOV(ISYM,1)
               LI=FLOV(ISYM,2)
               DO 40 I=FI,LI
                  ICNT=ICNT+1
                  IU=INO (MAX0(I,U))+MIN0(I,U)
                  VALR=0.0D+00
                  VALQ=0.0D+00
                  FJ=FLOV(ASYM+1,1)
                  LJ=FLOV(ASYM+1,2)
                  DO 50 J=FJ,LJ
                     BEJ=INO(BE)+J
                     BEJIU=EOFF (BEJ)+EADD(IU)
                     VALR=VALR+MOINTS(BEJIU)*T1(J,A)
                     IJ=INO (MAX0(I,J))+MIN0(I,J)
                     BEUIJ=EOFF (BEU)+EADD(IJ)
                     VALQ=VALQ+MOINTS(BEUIJ)*T1(J,A)
   50             CONTINUE
                  QMAT(ICNT)=-VALQ-VALQ -RMAT(ICNT)+VALR
                  RMAT(ICNT)=RMAT(ICNT)-VALR
   40          CONTINUE
   30       CONTINUE
            CALL WWRITW(TPRMAT,RMAT,INTOWP(IDIM),POINTR,POINTR)
            CALL WWRITW(TPQMAT,QMAT,INTOWP(ICNT),POINTQ,POINTQ)
   20    CONTINUE
   10 CONTINUE
  100 CONTINUE
      CALL RCLOSE(TPRMAT,3)
      CALL RCLOSE(TPQMAT,3)
      RETURN
      END
