      SUBROUTINE FOLD (SA,S,N)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     DIMENSION SA(1),S(N,N)
C
C     IJ=0
C     DO 2 I=2,N
C     DO 1 J=1,I-1
C     IJ=IJ+1
C     SA(IJ)=SA(IJ)+S(I,J)-S(J,I)
C   1 CONTINUE
C   2 CONTINUE
C     RETURN
C
      DIMENSION SA(1),S(1)
C
      K=1
      JISV=1
      DO 2 I=1,N
      JI=JISV
      IJ=I
      DO 1 J=1,I-1
      SA(K)=SA(K)+S(IJ)-S(JI)
      K=K+1
      IJ=IJ+N
      JI=JI+1
    1 CONTINUE
      JISV=JISV+N
    2 CONTINUE
      RETURN
      END
