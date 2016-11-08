       SUBROUTINE PRINMA(N,F)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION F(N,N)
       K = N
       MIN = 1
  5    MAX = MIN + K -1
       IF((MAX-MIN).GT.10) MAX = MIN + 9
       DO 10 I=1,N
       WRITE(6,1)(F(I,J),J=MIN,MAX)
10     CONTINUE
 1     FORMAT(1X,10F12.7)
       MIN = MAX + 1
       K = K- 10
       WRITE(6,2)
       IF(N.GE.MIN) GO TO 5
 2     FORMAT(/,13(10(1H*)))
       RETURN
       END
