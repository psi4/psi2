      SUBROUTINE MOCONV(SA,NNA,SM,NNM,ESO,U,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer nna,nnm,nbfso,nbasis,nbfao
      DIMENSION SA(NNA),SM(NNM)
      DIMENSION ESO(NBFSO,NBASIS),U(NBFAO,NBFAO),T(NBFAO,NBFAO)
      COMMON/BASIS/NBASIS,NTRI,NSYM,NBFAO,NBFSO,NBATRI
      COMMON/SIGNS/IOFF(256),IPRNT
      DATA ZERO / 0.0D+00 /
      DATA XLIMIT / 1.0D-35 /
C
C    TRANSFORMATION FROM SO TO MO BASIS
C
      DO 101 II=1,NBFSO
      DO 101 JJ=1,II
      IIJJ=IOFF(II)+JJ
      SAX=SA(IIJJ)
      IF(DABS(SAX).LT.XLIMIT) SAX=ZERO
      T(II,JJ)=SAX
      T(JJ,II)=SAX
  101 CONTINUE
      DO 103 II=1,NBFSO
      DO 103 J=1,NBASIS
      SUM=ZERO
      DO 102 JJ=1,NBFSO
      SUM=SUM+ESO(JJ,J)*T(II,JJ)
  102 CONTINUE
      U(II,J)=SUM
  103 CONTINUE
      DO 105 I=1,NBASIS
      DO 105 J=1,I
      IJ=IOFF(I)+J
      SUM=ZERO
      DO 104 II=1,NBFSO
      SUM=SUM+ESO(II,I)*U(II,J)
  104 CONTINUE
      SM(IJ)=SUM
  105 CONTINUE
C
      RETURN
      END
