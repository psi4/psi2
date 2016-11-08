C     ////////////////////////////////////////////////////
      SUBROUTINE NOM (U)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3)
      Y=SCALAR(U,U)
      X = 1.0D0/DSQRT(Y)
      DO 10 I=1,3
         U(I) = U(I)*X
   10 CONTINUE
      RETURN
      END
