      SUBROUTINE VEWPSX(A,B,C,Z,N)
      REAL*8 A(N),B(N),C(N),Z
      DO 1 I=1,N
         A(I)=B(I)+Z*C(I)
    1 CONTINUE
      RETURN
      END
