      SUBROUTINE VMOVE(A,B,N)
      REAL*8 A(N),B(N)
      DO 10 I=1,N
      A(I)=B(I)
   10 CONTINUE
      RETURN
      END
