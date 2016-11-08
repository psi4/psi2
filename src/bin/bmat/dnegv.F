      SUBROUTINE DNEGV(A,N)
      REAL*8 A(1)
C
C  A = -A
C
      DO 10 I=1,N
      A(I)=-A(I)
   10 CONTINUE
C
      RETURN
      END
