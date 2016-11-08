C
C ----------------------------------------------------------------------
C
      SUBROUTINE FT2   (T1,T1O,T2O,FIN,AUX4,XVAL,
     .                  ORBSYM,FLOV,FOFF,FADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL1,VAL2,TMP,TMV
      REAL*8  T1(NO,NV),T1O(NO,NV),T2O(*),FIN(*),XVAL(*),AUX4(NO,NV)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        FOFF(*),FADD(*)
C
      CALL ZERO(XVAL,NO)
C
      DO 191 I=1,NO
         IF(FZO(I).EQ.1)GO TO 191
         ISYM=ORBSYM(I)
         DO 181 A=1,NV
            IF(FZV(A).EQ.1)GO TO 181
            ASYM=ORBSYM(A+NO)
            IA=IPQ(A+NO)+I
            AISYM=IEOR(ASYM,ISYM)
            IA=ITR(A)+I
            DO 171 B=1,NV
               IF(FZV(B).EQ.1)GO TO 171
               BSYM=ORBSYM(B+NO)
               USYM=IEOR(BSYM,AISYM)+1
               FU=FLOV(USYM,1)
               LU=FLOV(USYM,2)
               FBE=FLOV(USYM,3)-NO
               LBE=FLOV(USYM,4)-NO
               ZL=ZLX(A,B)
               DO 111 U=FU,LU
                  IUAB=TOFF(I,U,ZL)+TADD(A,B)
                  UIAB=TOFF(U,I,ZL)+TADD(A,B)
                  VAL1=T2O(IUAB)+T1O(I,A)*T1O(U,B)
                  VAL2=T2O(UIAB)+T1O(U,A)*T1O(I,B)
                  XVAL(U)=VAL1+VAL1-VAL2
  111          CONTINUE
               DO 161 BE=FBE,LBE
                  BEB=ITV(MAX0(BE,B))+MIN0(BE,B)
                  IABEB=FOFF(BEB)+FADD(IA)
                  DO 151 U=FU,LU
                     TMV=AUX4(U,BE)
                     TMV=TMV+FIN(IABEB)*XVAL(U)
                     AUX4(U,BE)=TMV
  151             CONTINUE
  161          CONTINUE
  171       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
C     ADD INTO T1
C
      DO 291 U=1,NO
         IF(FZO(U).EQ.1)GO TO 291
         USYM=ORBSYM(U)
         FBE=FLOV(USYM+1,3)-NO
         LBE=FLOV(USYM+1,4)-NO
         DO 271 BE=FBE,LBE
            T1(U,BE)=T1(U,BE)+AUX4(U,BE)
  271    CONTINUE
  291 CONTINUE
C
      RETURN
      END
