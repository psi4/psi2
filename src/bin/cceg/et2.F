C
C ----------------------------------------------------------------------
C
      SUBROUTINE ET2   (T1O,T2O,EIN,AUX1,XVAL,
     .                  ORBSYM,FLOV,EOFF,EADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL1,VAL2,TMP,TMV
      REAL*8  T1O(NO,NV),T2O(*),EIN(*),XVAL(*),AUX1(NV,NO)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        EOFF(*),EADD(*)
C
      CALL ZERO(AUX1,NO*NV)
      CALL ZERO(XVAL,NV)
C
      DO 191 I=1,NO
         IF(FZO(I).EQ.1)GO TO 191
         ISYM=ORBSYM(I)
         DO 181 J=1,NO
            IF(FZO(J).EQ.1)GO TO 181
            JSYM=ORBSYM(J)
            JISYM=IEOR(JSYM,ISYM)
            DO 171 A=1,NV
               IF(FZV(A).EQ.1)GO TO 171
               ASYM=ORBSYM(A+NO)
               IA=ITR(A)+I
               BESYM=IEOR(ASYM,JISYM)+1
               FBE=FLOV(BESYM,3)-NO
               LBE=FLOV(BESYM,4)-NO
               FU=FLOV(BESYM,1)
               LU=FLOV(BESYM,2)
               DO 101 BE=FBE,LBE
                  ZL=ZLX(A,BE)
                  IJABE=TOFF(I,J,ZL)+TADD(A,BE)
                  JIABE=TOFF(J,I,ZL)+TADD(A,BE)
                  VAL1=T2O(IJABE)+T1O(I,A)*T1O(J,BE)
                  VAL2=T2O(JIABE)+T1O(J,A)*T1O(I,BE)
                  XVAL(BE)=VAL1+VAL1-VAL2
  101          CONTINUE
               DO 161 U=FU,LU
                  JU=ITR(MAX0(J,U))+MIN0(J,U)
                  IAJU=EOFF(IA)+EADD(JU)
                  TMP = -EIN(IAJU)
                  DO 151 BE=FBE,LBE
                     TMV = AUX1(BE,U)
                     TMV = TMV + TMP *XVAL(BE)
                     AUX1(BE,U) = TMV
  151             CONTINUE
  161          CONTINUE
  171       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
      RETURN
      END
