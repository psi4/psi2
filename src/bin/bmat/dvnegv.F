      SUBROUTINE DVNEGV(A,B,C,N)
      REAL*8 A(1),B(1),C(1)
C
C  A = B - C
C
      DO 10 I=1,N
      A(I)=B(I)-C(I)
   10 CONTINUE
C
      RETURN
      END
