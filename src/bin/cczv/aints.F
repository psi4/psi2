C
C----------------------------------------------------------------------
C
      SUBROUTINE AINTS(AMAT,MOINTS,NTRO,NO2,LENINT,NO,NTR,IOFF,JOUT,
     .                 INO,NOV,AOFF,NOFF,AADD,NADD,
     .                 FLOV,NIRRED,ORBSYM,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 VAL
      REAL*8 AMAT(NO2),MOINTS(LENINT)
      INTEGER IOFF(NTR),FLOV(NIRRED,4),ORBSYM(NO),FZO(*),FZV(*),
     .        AOFF(NOFF),AADD(NADD),INO(NOV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
C
      CALL RFILE(TPAMAT)
      PON=1
      DO 100 TSYM=1,NIRRED
      LIJ=0
C
      DO 90 I=1,NO
         IF(FZO(I).EQ.1) GO TO 90
         ISYM=ORBSYM(I)
         JSYM=IEOR(ISYM,TSYM-1)+1
         FJ=FLOV(JSYM,1)
         LJ=FLOV(JSYM,2)
         DO 89 J=FJ,LJ
            LIJ=LIJ+1
   89    CONTINUE
   90 CONTINUE
      NIJ(TSYM)=LIJ
C     WRITE(6,*)'TSYM=',TSYM,' NIJ FOR AMAT',LIJ
C
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         USYM=ORBSYM(U)
         DO 11 V=1,U
            IF(FZO(V).EQ.1)GO TO 11
            VSYM=ORBSYM(V)
            UVSYM=IEOR(USYM,VSYM)
            IF(UVSYM.NE.TSYM-1)GO TO 11
            ICNT=0
            DO 12 I=1,NO
               IF(FZO(I).EQ.1)GO TO 12
               ISYM=ORBSYM(I)
               JSYM=IEOR(ISYM,UVSYM)+1
               FJ=FLOV(JSYM,1)
               LJ=FLOV(JSYM,2)
               UI=INO(MAX0(U,I))+MIN0(U,I)
               DO 13 J=FJ,LJ
                  VJ=INO(MAX0(V,J))+MIN0(V,J)
                  UIVJ=AOFF(MAX0(UI,VJ))+AADD(MIN0(UI,VJ))
                  VAL=MOINTS(UIVJ)
                  ICNT=ICNT+1
C                 AMAT(UV,IJ)
                  AMAT(ICNT)=-VAL
   13          CONTINUE
   12       CONTINUE
         CALL WWRITW(TPAMAT,AMAT,INTOWP(LIJ),PON,PON)
   11    CONTINUE
   10 CONTINUE
  100 CONTINUE
      CALL RCLOSE(TPAMAT,3)
C
      RETURN
      END
