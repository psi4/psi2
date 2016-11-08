C
C----------------------------------------------------------------------
C
      SUBROUTINE WRTEM(GEMAT,ITR,ORBSYM,FLOV,NIRRED,NO,NV,EOFF,EADD)
      IMPLICIT INTEGER (A-Z)
      REAL*8 GEMAT(*),VAL
      INTEGER ITR(*),ORBSYM(*),FLOV(NIRRED,4),EOFF(*),EADD(*)
      DO 190 A=1,NV
      ASYM=ORBSYM(A+NO)
      DO 180 I=1,NO
      ISYM=ORBSYM(I)
      AISYM=IEOR(ASYM,ISYM)
      DO 170 J=1,NO
      JSYM=ORBSYM(J)
      DO 160 K=1,J
      KSYM=ORBSYM(K)
      JKSYM=IEOR(JSYM,KSYM)
      IF(AISYM.NE.JKSYM)GO TO 160
      AI=ITR(A)+I
      JK=ITR(J)+K
      KJ=ITR(K)+J
      AIJK=EOFF(AI)+EADD(JK)
C     AIKJ=EOFF(AI)+EADD(KJ)
C     VAL=(GEMAT(AIJK)+GEMAT(AIKJ))*0.5D0
      VAL= GEMAT(AIJK)             *0.125D0
                  WRITE (6,1) A+NO,I,J,K,VAL
    1 FORMAT (' GEMAT(',I3,',',I3,',',I3,',',I3,')=',F20.12,F20.12)
  160 CONTINUE
  170 CONTINUE
  180 CONTINUE
  190 CONTINUE
C
      RETURN
      END
