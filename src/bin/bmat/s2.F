C     /////////////////////////////////////////////////
      REAL*8 FUNCTION S2(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      Y = 1.0D0-X*X
      S2 = DSQRT(Y)
      RETURN
      END
