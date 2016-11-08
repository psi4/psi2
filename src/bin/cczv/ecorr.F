C
C-----------------------------------------------------------------------
C
      SUBROUTINE ECORR(T2,T1,MOINTS,LENINT,NO,NV,NO2,NV2,NTRO,NTRV,IOFF,
     .                 NTR,JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,ZLX,UOFF,
     .                 FLOV,ORBSYM,NIRRED,DOFF,NOFF,DADD,NADD,
     .                 TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      REAL*8 T2(NDIMT2),T1(NO,NV),MOINTS(LENINT),EC
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NV)
      INTEGER DOFF(NOFF),DADD(NADD)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV)
      INTEGER ZLX(NV,NV),UOFF(NO,NO,2)
      DIMENSION FLOV(NIRRED,4),ORBSYM(NO+NV)
C
      EC=0.0D+00
      DO 10 I=1,NO
         IF(FZO(I).EQ.1)GO TO 10
         ISYM=ORBSYM(I)
         DO 20 A=1,NV
            IF(FZV(A).EQ.1)GO TO 20
            ASYM=ORBSYM(A+NO)
            IASYM=IEOR(ISYM,ASYM)
            AI=INO(A)+I
            DO 30 J=1,NO
               IF(FZO(J).EQ.1)GO TO 30
               JSYM=ORBSYM(J)
               BSYM=IEOR(IASYM,JSYM)+1
               FB=FLOV(BSYM,3)-NO
               LB=FLOV(BSYM,4)-NO
               DO 40 B=FB,LB
                  BJ=INO(B)+J
                  ZL=ZLX(A,B)
                  IJAB=TOFF(I,J,ZL)+TADD(A,B)
                  JIAB=TOFF(J,I,ZL)+TADD(A,B)
                  AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  EC=EC+MOINTS(AIBJ)*(T2(IJAB)+T2(IJAB)-T2(JIAB)+
     .               T1(I,A)*T1(J,B)+T1(I,A)*T1(J,B)-T1(J,A)*T1(I,B))
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      WRITE(6,*) 'CORRELATION ENERGY FROM ZCCSD  = ',EC
      WRITE(4,*) 'CORRELATION ENERGY FROM ZCCSD  = ',EC
      RETURN
      END
