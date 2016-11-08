C
C ----------------------------------------------------------------------
C
      SUBROUTINE HVINT (T1O,T2O,DIN,FOCK,HV,XVAL,
     .                  ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .                  ZLX,NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,TMP
      REAL*8  T1O(NO,NV),T2O(*),DIN(*),FOCK(*),HV(NV,NV),XVAL(*)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),DOFF(*),DADD(*),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV)
C
      DO 132 BE=1,NV
         IF(FZV(BE).EQ.1)GO TO 132
         BESYM=ORBSYM(BE+NO)
         FA=FLOV(BESYM+1,3)-NO
         LA=FLOV(BESYM+1,4)-NO
         DO 131 A=FA,LA
            VAL=0.0D0
            FABE=IPQ(MAX0(A,BE)+NO)+MIN0(A,BE)+NO
            IF(A.NE.BE)VAL=VAL+FOCK(FABE)
            HV(BE,A)=VAL
  131    CONTINUE
  132 CONTINUE
C
      DO 139 B=1,NV
         IF(FZV(B).EQ.1)GO TO 139
         BSYM=ORBSYM(B+NO)
         DO 138 J=1,NO
            IF(FZO(J).EQ.1)GO TO 138
            JSYM=ORBSYM(J)
            JBSYM=IEOR(JSYM,BSYM)
            JB=ITR(B)+J
            DO 137 I=1,NO
               IF(FZO(I).EQ.1)GO TO 137
               ISYM=ORBSYM(I)
               BESYM=IEOR(ISYM,JBSYM)+1
               FBE=FLOV(BESYM,3)-NO
               LBE=FLOV(BESYM,4)-NO
               IB=ITR(B)+I
               DO 133 BE=FBE,LBE
                  ZL=ZLX(BE,B)
                  IJBEB=TOFF(I,J,ZL)+TADD(BE,B)
                  XVAL(BE)=T2O(IJBEB)+T1O(I,BE)*T1O(J,B)
  133          CONTINUE
               DO 135 A=FBE,LBE
                  IA=ITR(A)+I
                  JA=ITR(A)+J
                  IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                  IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                  TMP=DIN(IAJB)+DIN(IAJB)-DIN(IBJA)
                  DO 134 BE=FBE,LBE
                     HV(BE,A)=HV(BE,A)-TMP*XVAL(BE)
  134             CONTINUE
  135          CONTINUE
  137       CONTINUE
  138    CONTINUE
  139 CONTINUE
C
      RETURN
      END
