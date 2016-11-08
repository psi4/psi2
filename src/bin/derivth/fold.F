      SUBROUTINE FOLD(A,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N)
C
      DO 101 J=1,N
      DO 101 I=1,J-1
      A(I,J)=A(I,J)+A(J,I)
      A(J,I)=A(I,J)
  101 CONTINUE
C
      RETURN
      END
