C
C ----------------------------------------------------------------------
C
      SUBROUTINE HIAINT(T1O,T2O,DIN,FOCK,HIA,XVAL,
     .                  ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .                  ZLX,NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,TMP
      REAL*8  T1O(NO,NV),T2O(*),DIN(*),FOCK(*),HIA(NO,NV),XVAL(*)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),DOFF(*),DADD(*),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV)
C
      DO 137 A=1,NV
         IF(FZV(A).EQ.1)GO TO 137
         ASYM=ORBSYM(A+NO)
         FI=FLOV(ASYM+1,1)
         LI=FLOV(ASYM+1,2)
         DO 136 I=FI,LI
            VAL=0.0D0
            IA=ITR(A)+I
            FIA=IPQ(A+NO)+I
            DO 134 J=1,NO
               IF(FZO(J).EQ.1)GO TO 134
               JSYM=ORBSYM(J)
               FB=FLOV(JSYM+1,3)-NO
               LB=FLOV(JSYM+1,4)-NO
               JA=ITR(A)+J
               DO 130 B=FB,LB
                  JB=ITR(B)+J
                  IB=ITR(B)+I
                  IAJB=DOFF(MAX0(JB,IA))+DADD(MIN0(JB,IA))
                  IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                  VAL=VAL+(DIN(IAJB)+DIN(IAJB)-DIN(IBJA))*T1O(J,B)
  130          CONTINUE
  134       CONTINUE
            HIA(I,A)=VAL
     .                   +FOCK(FIA)
  136    CONTINUE
  137 CONTINUE
C
      RETURN
      END
