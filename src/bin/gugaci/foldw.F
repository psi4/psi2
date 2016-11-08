      SUBROUTINE FOLDW(SA,S,N)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 sqrt2,sqrt3,sqt1p5,asqrt2,hsqrt3
      common /sqrts/ sqrt2,sqrt3,sqt1p5,asqrt2,hsqrt3
C
C     DIMENSION SA(1),S(N,N)
C
C     SQRT2=DSQRT(2.0D+00)
C     IJ=0
C     DO 2 I=1,N
C     DO 1 J=1,I-1
C     IJ=IJ+1
C     SA(IJ)=SA(IJ)+S(I,J)+S(J,I)
C   1 CONTINUE
C     IJ=IJ+1
C     SA(IJ)=SA(IJ)+SQRT2*S(I,I)
C   2 CONTINUE
C     RETURN
C
      DIMENSION SA(1),S(1)
C
cets111090      SQRT2=DSQRT(2.0D+00)
      K=1
      JISV=1
      DO 2 I=1,N
      JI=JISV
      IJ=I
      DO 1 J=1,I-1
      SA(K)=SA(K)+S(IJ)+S(JI)
      K=K+1
      IJ=IJ+N
      JI=JI+1
    1 CONTINUE
      SA(K)=SA(K)+SQRT2*S(IJ)
      K=K+1
      JISV=JISV+N
    2 CONTINUE
      RETURN
      END
