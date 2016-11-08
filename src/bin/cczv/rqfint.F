C
C-----------------------------------------------------------------------
C
      SUBROUTINE RQFINT(QMAT,RMAT,MOINTS,T1,
     .                  LENINT,NO,NV,NO2,NV2,NTRO,NTRV,
     .                  NTR,IOFF,JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,FLOV,
     .                  ORBSYM,NIRRED,
     .                  FOFF,NVNV,FADD,NM,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 T1(NO,NV),QMAT(NV*NO),RMAT(NV*NO),MOINTS(LENINT),VALR,VALQ
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER FOFF(NVNV),FADD(NM)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
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
            DO 30 A=1,NV
               IF(FZV(A).EQ.1)GO TO 30
               BEA=INV (MAX0(BE,A))+MIN0(BE,A)
               ASYM=ORBSYM(A+NO)
               ISYM=IEOR(BEUSYM,ASYM)+1
               FI=FLOV(ISYM,1)
               LI=FLOV(ISYM,2)
               FB=FLOV(ISYM,3)-NO
               LB=FLOV(ISYM,4)-NO
               DO 40 I=FI,LI
                  VALQ=0.0D+00
                  VALR=0.0D+00
                  ICNT=ICNT+1
                  DO 50 B=FB,LB
                     BA=INV (MAX0(B,A))+MIN0(B,A)
                     BEUBA=FADD (BEU)+FOFF(BA)
                     VALQ=VALQ+MOINTS(BEUBA)*T1(I,B)
                     BU=INO(B)+U
                     BUBEA=FADD (BU)+FOFF(BEA)
                     VALR=VALR+MOINTS(BUBEA)*T1(I,B)
   50             CONTINUE
                  QMAT(ICNT)=QMAT(ICNT)+VALQ+VALQ-VALR
                  RMAT(ICNT)=RMAT(ICNT)+VALR
   40          CONTINUE
   30       CONTINUE
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
