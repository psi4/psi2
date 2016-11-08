C
C ----------------------------------------------------------------------
C
      SUBROUTINE HOINT (T1O,T2O,DIN,FOCK,HO,XVAL,
     .                  ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .                  ZLX,NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,TMP
      REAL*8  T1O(NO,NV),T2O(*),DIN(*),FOCK(*),HO(NO,NO),XVAL(*)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),DOFF(*),DADD(*),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV)
C
      DO 132 U=1,NO
         IF(FZO(U).EQ.1)GO TO 132
         USYM=ORBSYM(U)
         FI=FLOV(USYM+1,1)
         LI=FLOV(USYM+1,2)
         DO 131 I=FI,LI
            VAL=0.0D0
            FUI=IPQ(MAX0(U,I))+MIN0(U,I)
            IF(I.NE.U)VAL=VAL+FOCK(FUI)
            HO(U,I)=VAL
  131    CONTINUE
  132 CONTINUE
C
      DO 193 A=1,NV
         IF(FZV(A).EQ.1)GO TO 193
         ASYM=ORBSYM(A+NO)
         DO 183 B=1,NV
            IF(FZV(B).EQ.1)GO TO 183
            BSYM=ORBSYM(B+NO)
            BASYM=IEOR(BSYM,ASYM)
            ZL=ZLX(A,B)
            DO 173 J=1,NO
               IF(FZO(J).EQ.1)GO TO 173
               JSYM=ORBSYM(J)
               USYM=IEOR(JSYM,BASYM)+1
               FU=FLOV(USYM,1)
               LU=FLOV(USYM,2)
               JA=ITR(A)+J
               JB=ITR(B)+J
               DO 143 U=FU,LU
                  UJAB=TOFF(U,J,ZL)+TADD(A,B)
                  XVAL(U)=T2O(UJAB)+T1O(U,A)*T1O(J,B)
  143          CONTINUE
               DO 163 I=FU,LU
                  IA=ITR(A)+I
                  IB=ITR(B)+I
                  IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                  IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                  TMP=DIN(IAJB)+DIN(IAJB)-DIN(IBJA)
                  DO 153 U=FU,LU
                     HO(U,I)=HO(U,I)+TMP*XVAL(U)
  153             CONTINUE
  163          CONTINUE
  173       CONTINUE
  183    CONTINUE
  193 CONTINUE
C
      RETURN
      END
