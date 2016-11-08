      SUBROUTINE FORMQ (QOO,QVV,QOV,TAU,T2,T1,Z2,Z1,NO,NV,NO2,NV2,NTRO,
     .                  NTRV,IOFF,NTR,JOUT,OPTION,
     .                  ORBSYM,FLOV,ZLX,ITR,ITV,NT3,NIRRED,NTAU,NOV,
     .                  T3OFF,T3ADD,W3,MOO,NVV,NDIMT2,TOFF,TADD,
     .                  FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 OPTION
      COMMON/TAPES/ITAP81,ITAP82,ITAP83,ITAP99,itap76,i76,itap74
      REAL*8 VAL,VAL1,VAL2,VAL3,VAL4,VAL5,VAL6,VAL7,VAL8
      REAL*8 QOO(NTRO),QVV(NTRV),QOV(NO,NV),TAU(NDIMT2),T2(NDIMT2),
     .       T1(NO,NV),Z2(NDIMT2),Z1(NO,NV),W3(NTAU),
     .       MOO(NO,NO),NVV(NV,NV)
      INTEGER IOFF(NTR),ORBSYM(NO+NV),FLOV(NIRRED,4),ZLX(NV,NV),
     .        ITR(NO),ITV(NOV),NT3(NIRRED),FZO(NO),FZV(NV),
     .        T3OFF(NO,NO,2,NIRRED),T3ADD(NV,NV,NIRRED),
     .        TOFF(NO,NO,2),TADD(NV,NV)
C
      CALL ZERO(QOO,NTRO)
      CALL ZERO(QVV,NTRV)
      CALL ZERO(QOV,NO*NV)
C
      DO 10 I=1,NO
         IF(FZO(I).EQ.1)GO TO 10
         ISYM=ORBSYM(I)
         FA=FLOV(ISYM+1,3)-NO
         TA=FLOV(ISYM+1,4)-NO
         DO 11 A=FA,TA
            VAL1=0.0D+00
            DO 12 J=1,NO
               IF(FZO(J).EQ.1)GO TO 12
               JSYM=ORBSYM(J)
               FB=FLOV(JSYM+1,3)-NO
               TB=FLOV(JSYM+1,4)-NO
               DO 13 B=FB,TB
                  ZLAB=ZLX(A,B)
                  IJAB=TOFF(I,J,ZLAB) +TADD(A,B)
                  JIAB=TOFF(J,I,ZLAB) +TADD(A,B)
                  VAL1=VAL1+Z1(J,B)*(T2(IJAB)+T2(IJAB)-TAU(JIAB))
   13          CONTINUE
   12       CONTINUE
C
            VAL2=0.0D+00
            DO 28 J=1,NO
               IF(FZO(J).EQ.1)GO TO 28
               JSYM=ORBSYM(J)
               DO 29 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 29
                  KSYM=ORBSYM(K)
                  JKSYM=IEOR(JSYM,KSYM)
                  CSYM=IEOR(JKSYM,ISYM)+1
                  FC=FLOV(CSYM,3)-NO
                  TC=FLOV(CSYM,4)-NO
                  DO 30 B=FA,TA
                     DO 31 C=FC,TC
                        ZLAC=ZLX(A,C)
                        ZLBC=ZLX(B,C)
                        JKBC=TOFF(J,K,ZLBC)+TADD(B,C)
                        JKAC=TOFF(J,K,ZLAC) +TADD(A,C)
                        VAL2=VAL2-T1(I,B)*T2(JKAC)*Z2(JKBC)
   31                CONTINUE
   30             CONTINUE
   29          CONTINUE
   28       CONTINUE
C
            FJ=FLOV(ISYM+1,1)
            TJ=FLOV(ISYM+1,2)
            DO 128 J=FJ,TJ
               DO 129 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 129
                  KSYM=ORBSYM(K)
                  IKSYM=IEOR(ISYM,KSYM)
                  DO 130 B=1,NV
                     IF(FZV(B).EQ.1)GO TO 130
                     BSYM=ORBSYM(B+NO)
                     CSYM=IEOR(IKSYM,BSYM)+1
                     FC=FLOV(CSYM,3)-NO
                     TC=FLOV(CSYM,4)-NO
                     DO 131 C=FC,TC
                        ZLBC=ZLX(B,C)
                        JKBC=TOFF(J,K,ZLBC) +TADD(B,C)
                        IKBC=TOFF(I,K,ZLBC) +TADD(B,C)
                        VAL2=VAL2-T1(J,A)*T2(IKBC)*Z2(JKBC)
  131                CONTINUE
  130             CONTINUE
  129          CONTINUE
  128       CONTINUE
            QOV(I,A)=QOV(I,A)+VAL1+T1(I,A)+T1(I,A)+Z1(I,A)+VAL2
   11    CONTINUE
   10 CONTINUE
C
C     IF(OPTION.EQ.'SDT1')THEN
C     I81=1
C     DO 190 C=1,NV
C     IF(FZV(C).EQ.1)GO TO 190
C     CSYM=ORBSYM(C+NO)
C     DO 185 K=1,NO
C     IF(FZO(K).EQ.1)GO TO 185
C     KSYM=ORBSYM(K)
C     KCSYM=IEOR(KSYM,CSYM)+1
C     DIM=NT3(KCSYM)
C     CALL WREADW(ITAP81,W3,INTOWP(DIM),I81,I81)
C     DO 180 A=1,NV
C     IF(FZV(A).EQ.1)GO TO 180
C     ASYM=ORBSYM(A+NO)
C     FI=FLOV(ASYM+1,1)
C     LI=FLOV(ASYM+1,2)
C     DO 175 I=FI,LI
C     DO 170 J=1,NO
C     IF(FZO(J).EQ.1)GO TO 170
C     JSYM=ORBSYM(J)
C     JKSYM=IEOR(JSYM,KSYM)
C     BSYM=IEOR(JKSYM,CSYM)+1
C     FB=FLOV(BSYM,3)-NO
C     LB=FLOV(BSYM,4)-NO
C     DO 160 B=FB,LB
C     BC=IOFF(MAX0(B,C))+MIN0(B,C)
C     ZLBC=ZLX(B,C)
C     ZLAB=ZLX(A,B)
C     JK=ZOFF(J,K,ZLBC)
C     IJBA=T3OFF(J,I,ZLAB,KCSYM)+T3ADD(B,A,KCSYM)
C     IJAB=T3OFF(I,J,ZLAB,KCSYM)+T3ADD(A,B,KCSYM)
CC!   QOV(I,A)=QOV(I,A)+(W3(IJAB)-W3(IJBA))*Z2(JK,BC)
C 160 CONTINUE
C 170 CONTINUE
C 175 CONTINUE
C 180 CONTINUE
C 185 CONTINUE
C 190 CONTINUE
C     ENDIF
C
      DO 15 A=1,NV
         IF(FZV(A).EQ.1)GO TO 15
         ASYM=ORBSYM(A+NO)
         FB=FLOV(ASYM+1,3)-NO
         TB=FLOV(ASYM+1,4)-NO
         IF(TB.GT.A)TB=A
         DO 16 B=FB,TB
            AB=IOFF(A)+B
            VAL1=0.0D+00
            VAL2=0.0D+00
            DO 17 I=1,NO
               IF(FZO(I).EQ.1)GO TO 17
               ISYM=ORBSYM(I)
               VAL1=VAL1+T1(I,A)*Z1(I,B)+T1(I,B)*Z1(I,A)
               DO 18 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 18
                  JSYM=ORBSYM(J)
                  IJSYM=IEOR(ISYM,JSYM)
                  CSYM=IEOR(IJSYM,ASYM)+1
                  FC=FLOV(CSYM,3)-NO
                  TC=FLOV(CSYM,4)-NO
                  DO 19 C=FC,TC
                     ZLBC=ZLX(B,C)
                     ZLAC=ZLX(A,C)
                     IJAC=TOFF(I,J,ZLAC)+TADD(A,C)
                     IJBC=TOFF(I,J,ZLBC)+TADD(B,C)
                     VAL2=VAL2+T2(IJAC)*Z2(IJBC)+T2(IJBC)*Z2(IJAC)
   19             CONTINUE
   18          CONTINUE
   17       CONTINUE
            QVV(AB)=QVV(AB)+(VAL1+VAL2)
            IF(OPTION.EQ.'SDT1'.or.option.eq.'CCT ')
     .         QVV(AB)=QVV(AB)+NVV(A,B)+NVV(B,A)
   16    CONTINUE
   15 CONTINUE
C
      DO 22 I=1,NO
         IF(FZO(I).EQ.1)GO TO 22
         ISYM=ORBSYM(I)
         FJ=FLOV(ISYM+1,1)
         TJ=FLOV(ISYM+1,2)
         IF(TJ.GT.I)TJ=I
         DO 23 J=FJ,TJ
            IJ=IOFF(I)+J
            VAL1=0.0D+00
            VAL2=0.0D+00
            DO 24 A=1,NV
               IF(FZV(A).EQ.1)GO TO 24
               ASYM=ORBSYM(A+NO)
               IASYM=IEOR(ISYM,ASYM)
               VAL1=VAL1-T1(I,A)*Z1(J,A)-T1(J,A)*Z1(I,A)
               DO 25 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 25
                  KSYM=ORBSYM(K)
                  BSYM=IEOR(IASYM,KSYM)+1
                  FB=FLOV(BSYM,3)-NO
                  TB=FLOV(BSYM,4)-NO
                  DO 26 B=FB,TB
                     ZLAB=ZLX(A,B)
                     IKAB=TOFF(I,K,ZLAB)+TADD(A,B)
                     JKAB=TOFF(J,K,ZLAB)+TADD(A,B)
                     VAL2=VAL2-T2(IKAB)*Z2(JKAB)-T2(JKAB)*Z2(IKAB)
   26             CONTINUE
   25          CONTINUE
   24       CONTINUE
            QOO(IJ)=QOO(IJ)+(VAL1+VAL2)
            IF(OPTION.EQ.'SDT1'.or.option.eq.'CCT ')
     .         QOO(IJ)=QOO(IJ)-MOO(I,J)-MOO(J,I)
   23    CONTINUE
   22 CONTINUE
C
 6044 FORMAT (////,2X,' ******   QOO MATRIX   ******')
 6045 FORMAT (////,2X,' ******   QVV MATRIX   ******')
 6046 FORMAT (////,2X,' ******   QOV MATRIX   ******')
C
c     WRITE (JOUT,6046)
c     CALL MATOUT(QOV,NO,NV,NO,NV,JOUT)
c     WRITE (JOUT,6044)
c     CALL PRINT(QOO,NTRO,NO,JOUT)
c     WRITE (JOUT,6045)
c     CALL PRINT(QVV,NTRV,NV,JOUT)
C
      RETURN
      END
