      SUBROUTINE VADD(A,B,C,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N),B(N),C(N)
C
      DO 101 I=1,N
      A(I)=B(I)+C(I)
  101 CONTINUE
      RETURN
      END
