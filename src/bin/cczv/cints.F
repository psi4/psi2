C
C-----------------------------------------------------------------------
C
      SUBROUTINE CINTS(EMAT,FMAT,MOINTS,NO,NV,T1,LENINT,NTRO,
     .                 NTR,IOFF,JOUT,INO,INTRO,INTRV,INV,NOV,N2OV,FLOV,
     .                 ORBSYM,NIRRED,COFF,NOFF,CADD,NADD,
     .                 EMOFF,EMADD,LEMAT,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 VAL
      REAL*8 MOINTS(LENINT),EMAT(LEMAT),FMAT(NV*NV),T1(NO,NV)
      INTEGER IOFF(NTR),FZO(*),FZV(*),
     .        INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV),
     .        COFF(NOFF),CADD(NADD),EMOFF(NO*NV),EMADD(NO*NO),
     .        FLOV(NIRRED,4),ORBSYM(NO+NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
C
      CALL ZERO(EMAT,LEMAT)
C
      DO 10 A=1,NV
         IF(FZV(A).EQ.1)GO TO 10
         ASYM=ORBSYM(A+NO)
         DO 11 I=1,NO
            IF(FZO(I).EQ.1)GO TO 11
            ISYM=ORBSYM(I)
            AISYM=IEOR(ASYM,ISYM)
            DO 14 J=1,NO
               IF(FZO(J).EQ.1)GO TO 14
               JSYM=ORBSYM(J)
               IJ=INO (MAX0(I,J))+MIN0(I,J)
               IJX=INO(I)+J
               KSYM=IEOR(AISYM,JSYM)+1
               FK=FLOV(KSYM,1)
               LK=FLOV(KSYM,2)
               FB=FLOV(KSYM,3)-NO
               LB=FLOV(KSYM,4)-NO
               DO 15 K=FK,LK
                  KAX=INO(A)+K
                  IJKA=EMOFF(KAX)+EMADD(IJX)
                  DO 16 B=FB,LB
                     AB=INV (MAX0(A,B))+MIN0(A,B)
                     ABIJ=COFF(AB)+CADD(IJ)
                     VAL=MOINTS(ABIJ)
                     EMAT(IJKA)=EMAT(IJKA)+VAL*T1(K,B)
   16             CONTINUE
   15          CONTINUE
   14       CONTINUE
   11    CONTINUE
   10 CONTINUE
C
      CALL RFILE(TPFMAT)
      POINT=1
      DO 23 TSYM=1,NIRRED
         DO 20 A=1,NV
            IF(FZV(A).EQ.1)GO TO 20
            ASYM=ORBSYM(A+NO)
            DO 21 B=1,NV
               IF(FZV(B).EQ.1)GO TO 21
               BSYM=ORBSYM(B+NO)
               ABSYM=IEOR(ASYM,BSYM)
               IF(ABSYM.NE.TSYM-1)GO TO 21
               AB=INV (MAX0(A,B))+MIN0(A,B)
               ICNT=0
               DO 17 C=1,NV
                  IF(FZV(C).EQ.1)GO TO 17
                  CSYM=ORBSYM(C+NO)
                  ISYM=IEOR(ABSYM,CSYM)+1
                  FJ=FLOV(CSYM+1,1)
                  LJ=FLOV(CSYM+1,2)
                  FI=FLOV(ISYM,1)
                  LI=FLOV(ISYM,2)
                  DO 18 I=FI,LI
                     ICNT=ICNT+1
                     VAL=0.0D0
                     DO 19 J=FJ,LJ
                        IJ=INO (MAX0(I,J))+MIN0(I,J)
                        ABIJ=COFF(AB)+CADD(IJ)
                        VAL=VAL+MOINTS(ABIJ)*T1(J,C)
   19                CONTINUE
CF                   FMAT(A,B,C,I)=FMAT(A,B,C,I)+VAL
                     FMAT(ICNT)=+VAL
   18             CONTINUE
   17          CONTINUE
               CALL WWRITW(TPFMAT,ICNT,1,POINT,POINT)
               CALL WWRITW(TPFMAT,FMAT,INTOWP(ICNT),POINT,POINT)
   21       CONTINUE
   20    CONTINUE
   23 CONTINUE
      CALL RCLOSE(TPFMAT,3)
C
      RETURN
      END
