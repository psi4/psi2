C
C ----------------------------------------------------------------------
C
      SUBROUTINE INITT (T1,T2,DIN,FOCK,
     .                  ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .                  NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL
      REAL*8  T1(NO,NV),T2(*),DIN(*),FOCK(*)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),DOFF(*),DADD(*),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),TOFF(NO,NO,2),TADD(NV,NV)
C
      CALL ZERO(T1,NO*NV)
C
      DO 146 I=1,NO
         IF(FZO(I).EQ.1)GO TO 146
         ISYM=ORBSYM(I)
         IY=IPQ(I)+I
         DO 136 J=1,NO
            IF(FZO(J).EQ.1)GO TO 136
            JSYM=ORBSYM(J)
            JISYM=IEOR(JSYM,ISYM)
            JJ=IPQ(J)+J
            DO 126 A=1,NV
               IF(FZV(A).EQ.1)GO TO 126
               ASYM=ORBSYM(A+NO)
               BSYM=IEOR(ASYM,JISYM)+1
               FB=FLOV(BSYM,3)-NO
               LB=FLOV(BSYM,4)-NO
               AA=IPQ(A+NO)+A+NO
               IA=ITR(A)+I
               IF(LB.GT.A)LB=A
               DO 116 B=FB,LB
                  BY=IPQ(B+NO)+B+NO
                  JB=ITR(B)+J
                  IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                  IJAB=TOFF(I,J,1)+TADD(A,B)
                  VAL=FOCK(IY)+FOCK(JJ)-FOCK(AA)-FOCK(BY)
                  T2(IJAB)=DIN(IAJB)/VAL
  116          CONTINUE
  126       CONTINUE
  136    CONTINUE
  146 CONTINUE
C
      RETURN
      END
