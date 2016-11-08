C
C-----------------------------------------------------------------------
C
      SUBROUTINE PCORR(T2,T1,II,XE,UOFF,VADD,IPQ,ORBSYM,FLOV,
     .                 NIT,NO,NV,NT,NT2,NTIN,NDIMT2,NIRRED,FZO,FZV,FOCK,
     .                 DOFF,DADD,NM,ITR)
      IMPLICIT INTEGER (A-Z)
      REAL*8 XE,VAL1,VAL2
      REAL*8 T2(*),T1(NO,NV),FOCK(*),II(*)
      INTEGER ORBSYM(NT),FLOV(NIRRED,4),IPQ(NT2),UOFF(NO,NO,2),
     .        VADD(NV,NV),FZO(NO),FZV(NV),DOFF(NM),DADD(NM),ITR(NV)
C
C     CALCULATE THE CORRELATION ENERGY
C
      XE=0.0D0
      DO 81 I=1,NO
         IF(FZO(I).EQ.1)GO TO 81
         ISYM=ORBSYM(I)
         DO 80 J=1,NO
            IF(FZO(J).EQ.1)GO TO 80
            JSYM=ORBSYM(J)
            IJSYM=IEOR(ISYM,JSYM)
            DO 78 A=1,NV
               IF(FZV(A).EQ.1)GO TO 78
               ASYM=ORBSYM(A+NO)
               BSYM=IEOR(IJSYM,ASYM)+1
               FB=FLOV(BSYM,3)-NO
               LB=FLOV(BSYM,4)-NO
               IF(LB.GT.A-1)LB=A-1
               IA=ITR(A)+I
               DO 76 B=FB,LB
                  JB=ITR(B)+J
                  IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                  IJAB=UOFF(I,J,1)+VADD(A,B)
                  JIAB=UOFF(J,I,1)+VADD(A,B)
                  VAL1=T2(IJAB)
                  VAL2=T2(JIAB)
                  XE=XE+II(IAJB)*(VAL1+VAL1-VAL2)
   76          CONTINUE
               FB=FLOV(BSYM,3)-NO
               LB=FLOV(BSYM,4)-NO
               IF(FB.LT.A)FB=A
               DO 77 B=FB,LB
                  JB=ITR(B)+J
                  IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                  IJAB=UOFF(J,I,1)+VADD(B,A)
                  JIAB=UOFF(I,J,1)+VADD(B,A)
                  VAL1=T2(IJAB)
                  VAL2=T2(JIAB)
                  XE=XE+II(IAJB)*(VAL1+VAL1-VAL2)
   77          CONTINUE
   78       CONTINUE
   80    CONTINUE
   81 CONTINUE
C
C     ADD FOCK MATRIX DEPENDENT TERMS
C
      DO 95 I=1,NO
         IF(FZO(I).EQ.1)GO TO 95
         DO 93 A=1,NV
            IF(FZV(A).EQ.1)GO TO 93
            IA=IPQ(A+NO)+I
            XE=XE+2.D0*FOCK(IA)*T1(I,A)
   93    CONTINUE
   95 CONTINUE
C
C     WRITE(6,9000)
      WRITE(6,9010) NIT,XE,XE
      WRITE(4,9010) NIT,XE,XE
 9000 FORMAT(/)
 9010 FORMAT('NIT =',I3,5X,'CORREL.ENERGY =',F16.12,3X,D23.16)
      END
