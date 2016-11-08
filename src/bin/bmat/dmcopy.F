      SUBROUTINE DMCOPY(A,B,N)
      REAL*8 A(50,1),B(50,1)
C
C  COPIES MATRIX B INTO A
C
      DO 10 I=1,N
      DO 10 J=1,N
      A(I,J)=B(I,J)
   10 CONTINUE
C
      RETURN
      END
