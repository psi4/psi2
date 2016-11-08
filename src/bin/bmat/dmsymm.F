      SUBROUTINE DMSYMM(A,N)
      REAL*8 A(50,1)
      DATA HALF/0.5D0/
C
C  SYMMETRIZES A MATRIX A
C
      DO 10 I=2,N
      DO 10 J=1,I-1
      A(I,J)=HALF*(A(I,J)+A(J,I))
      A(J,I)=A(I,J)
   10 CONTINUE
C
      RETURN
      END
