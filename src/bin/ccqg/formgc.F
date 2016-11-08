C THIS SHOULD BE COMPILED AT LEVEL 2 TO AVOID DUMMY ERROR MESSAGES! on ibm
      SUBROUTINE FORMGC(GCMAT,TAU,T2,T1,Z2,Z1,NO,NV,NO2,NV2,NTRO,NTRV,
     .                  IOFF,NTR,JOUT,ZLX,FLOV,ORBSYM,NIRRED,
     .                  OPTION,NVV,NDIMT2,TOFF,TADD,
     .                  COFF,CADD,NSGOV,INO,INV,NOV,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 OPTION
      REAL*8 VAL,VAL2,FACIJ,FACAB
      REAL*8 GCMAT(NSGOV),Z2(NDIMT2),T2(NDIMT2),T1(NO,NV),
     .       Z1(NO,NV),TAU(NDIMT2),NVV(NV,NV)
      INTEGER IOFF(NTR),ZLX(NV,NV),FLOV(NIRRED,4),ORBSYM(NO+NV),
     .        TOFF(NO,NO,2),TADD(NV,NV),COFF(NV*NV),CADD(NO*NO),
     .        INO(NOV),INV(NOV),FZO(NO),FZV(NV)
C
      CALL ZERO(GCMAT,NSGOV)
C
      DO 10 A=1,NV
         IF(FZV(A).EQ.1)GO TO 10
         ASYM=ORBSYM(A+NO)
         DO 20 B=1,NV
            IF(FZV(B).EQ.1)GO TO 20
            BSYM=ORBSYM(B+NO)
            FACAB=1.0D+00
            IF (A.EQ.B) FACAB=2.0D+00
            AB=INV(MAX0(A,B))+MIN0(A,B)
            ABSYM=IEOR(ASYM,BSYM)
            DO 30 I=1,NO
               IF(FZO(I).EQ.1)GO TO 30
               ISYM=ORBSYM(I)
               IASYM=IEOR(ISYM,ASYM)
               JSYM=IEOR(ABSYM,ISYM)+1
               FJ=FLOV(JSYM,1)
               TJ=FLOV(JSYM,2)
               DO 40 J=FJ,TJ
                  FACIJ=1.0D+00
                  IF (I.EQ.J) FACIJ=2.0D+00
                  IJ=INO(MAX0(I,J))+MIN0(I,J)
                  ABIJ=COFF(AB)+CADD(IJ)
                  GCMAT(ABIJ)=GCMAT(ABIJ)-Z1(I,A)*T1(J,B)*FACIJ*FACAB
            IF(I.EQ.J.AND.(OPTION.EQ.'SDT1'.or.option.eq.'CCT '))
     .               GCMAT(ABIJ)=GCMAT(ABIJ)+NVV(A,B)*4.D0*FACAB
                  VAL=0.0D+00
                  DO 50 L=1,NO
                     IF(FZO(L).EQ.1)GO TO 50
                     LSYM=ORBSYM(L)
                     DSYM=IEOR(IASYM,LSYM)+1
                     FD=FLOV(DSYM,3)-NO
                     TD=FLOV(DSYM,4)-NO
                     DO 60 D=FD,TD
                        ZLAD=ZLX(A,D)
                        ZLBD=ZLX(B,D)
                        LJBD=TOFF(L,J,ZLBD)+TADD(B,D)
                        LIAD=TOFF(L,I,ZLAD)+TADD(A,D)
                        JLBD=TOFF(J,L,ZLBD)+TADD(B,D)
                        ILAD=TOFF(I,L,ZLAD)+TADD(A,D)
                        VAL=VAL-TAU(LIAD)*Z2(LJBD)-T2(ILAD)*Z2(JLBD)
   60                CONTINUE
   50             CONTINUE
                  GCMAT(ABIJ)=GCMAT(ABIJ)+VAL*FACIJ*FACAB
   40          CONTINUE
               IF(ASYM.NE.BSYM)GO TO 48
               DO 47 J=I,I
                  IF(FZO(J).EQ.1)GO TO 47
                  IJ=INO(I)+J
                  ABIJ=COFF(AB)+CADD(IJ)
                  VAL=0.0D+00
                  VAL2=0.0D+00
                  DO 57 K=1,NO
                     IF(FZO(K).EQ.1)GO TO 57
                     KSYM=ORBSYM(K)
                     VAL=VAL+T1(K,A)*Z1(K,B)
                     DO 67 L=1,NO
                        IF(FZO(L).EQ.1)GO TO 67
                        LSYM=ORBSYM(L)
                        KLSYM=IEOR(KSYM,LSYM)
                        CSYM=IEOR(KLSYM,ASYM)+1
                        FC=FLOV(CSYM,3)-NO
                        TC=FLOV(CSYM,4)-NO
                        DO 77 C=FC,TC
                           ZLAC=ZLX(A,C)
                           ZLBC=ZLX(B,C)
                           KLAC=TOFF(K,L,ZLAC)+TADD(A,C)
                           KLBC=TOFF(K,L,ZLBC)+TADD(B,C)
                           VAL2=VAL2+T2(KLAC)*Z2(KLBC)
   77                   CONTINUE
   67                CONTINUE
   57             CONTINUE
                  GCMAT(ABIJ)=GCMAT(ABIJ)+(VAL+VAL+VAL+VAL+
     .                                     VAL2+VAL2+VAL2+VAL2)*FACAB
   47          CONTINUE
   48          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END
