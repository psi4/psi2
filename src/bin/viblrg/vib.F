       SUBROUTINE VIB(N,X)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION X(N)
       FACTOR = 1302.79426D+00
       DO 10 I=1,N
       V = X(I)
       IF(V.LT.0) V = -V
       V = DSQRT(V) * FACTOR
       IF(X(I).LT.0) V = -V
 10    X(I) = V
       RETURN
       END
