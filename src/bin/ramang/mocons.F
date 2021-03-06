      SUBROUTINE MOCONS(SA,NNA,SM,NNM,EAO,U,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SA(NNA),SM(NNM)
      DIMENSION EAO(NBFAO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)
      COMMON/BASIS/NBASIS,NBFAO,NTRI,NTRI2,NBATRI,NBASQ
      COMMON/SIGNS/IOFF(466),IPRNT
      DATA ZERO / 0.0D+00 /
C
C   TRANSFORM INTEGRALS FROM AO TO MO BASIS
      DO 101 II=1,NBASIS
      DO 101 JJ=1,NBASIS
      IIJJ=IOFF(MAX0(II,JJ))+MIN0(II,JJ)
      T(II,JJ)=SA(IIJJ)
  101 CONTINUE
      DO 103 II=1,NBASIS
      DO 103 J=1,NBASIS
      SUM=ZERO
      DO 102 JJ=1,NBASIS
      SUM=SUM+EAO(JJ,J)*T(II,JJ)
  102 CONTINUE
      U(II,J)=SUM
  103 CONTINUE
      DO 105 I=1,NBASIS
      DO 105 J=1,I
      IJ=IOFF(I)+J
      SUM=ZERO
      DO 104 II=1,NBASIS
      SUM=SUM+EAO(II,I)*U(II,J)
  104 CONTINUE
      SM(IJ)=SUM
  105 CONTINUE
C
      RETURN
      END
