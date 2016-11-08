C
C-----------------------------------------------------------------------
C
      SUBROUTINE RCINTS(RMAT,MOINTS,NO,NV,LENINT,NTRO,
     .                  NTR,IOFF,JOUT,INO,INTRO,INTRV,INV,NOV,N2OV,
     .                  FLOV,ORBSYM,NIRRED,COFF,NOFF,CADD,NADD,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 VAL
      REAL*8 MOINTS(LENINT),RMAT(NV*NO)
      INTEGER IOFF(NTR),FZO(*),FZV(*),
     .        INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV),
     .        COFF(NOFF),CADD(NADD),
     .        FLOV(NIRRED,4),ORBSYM(NO+NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
C
      CALL RFILE(TPRMAT)
      POINT=1
      DO 50 TSYM=1,NIRRED
         DO 10 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 10
            BESYM=ORBSYM(BE+NO)
            DO 20 U=1,NO
               IF(FZO(U).EQ.1)GO TO 20
               USYM=ORBSYM(U)
               BEUSYM=IEOR(BESYM,USYM)
               IF(BEUSYM.NE.TSYM-1)GO TO 20
               ICNT=0
               DO 30 A=1,NV
                  IF(FZV(A).EQ.1)GO TO 30
                  ASYM=ORBSYM(A+NO)
                  BEA=INV (MAX0(BE,A))+MIN0(BE,A)
                  ISYM=IEOR(BEUSYM,ASYM)+1
                  FI=FLOV(ISYM,1)
                  LI=FLOV(ISYM,2)
                  DO 40 I=FI,LI
                     ICNT=ICNT+1
                     IU=INO (MAX0(I,U))+MIN0(I,U)
                     BEAIU=COFF(BEA)+CADD(IU)
CR                   RMAT(BE,U,A,I)
                     RMAT(ICNT)=MOINTS(BEAIU)
   40             CONTINUE
   30          CONTINUE
               CALL WWRITW(TPRMAT,RMAT,INTOWP(ICNT),POINT,POINT)
   20       CONTINUE
   10    CONTINUE
   50 CONTINUE
      CALL RCLOSE(TPRMAT,3)
C
      RETURN
      END
