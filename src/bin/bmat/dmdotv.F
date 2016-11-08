      SUBROUTINE DMDOTV(A,B,C,N)
      REAL*8 A(1),B(50,1),C(1)
      DATA ZERO/0.0D0/
C
C  A = B * C
C
      DO 10 I=1,N
      A(I)=ZERO
      DO 10 J=1,N
      A(I)=A(I)+B(I,J)*C(J)
   10 CONTINUE
C
      RETURN
      END
