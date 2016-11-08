      SUBROUTINE SHMIDT(N)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /INOUT/  IOUT,ICHECK,ITAP30
      COMMON/SCHMDT/ A(18,18),DD(18)
C
      PARAMETER (ZERO=0.0D+00, TOL=1.0D-10)
    1 FORMAT(2X,' D (FIRST TIME)= ',D20.10)
    2 FORMAT(2X,' D (SECOND TIME)= ',D20.10)
      D=ZERO
CCCCC WRITE(IOUT,*) ' INSIDE SHMIDT '
C     CALL MATOUT(A,18,18,N,N,IOUT)
      DO 10 I=1,N
C     WRITE(IOUT,1)D
   10 D=D+A(I,1)**2
C     IF(D .LT. TOL) GO TO 20
      D=DSQRT(D)
      IF(D.LT.TOL) GO TO 20
      DO 15 I=1,N
   15 A(I,1)=A(I,1)/D
   20 DD(1)=D
      IF(N.EQ.1) RETURN
      DO 100 J=2,N
         KMAX=J-1
         DO 30 K=1,KMAX
            DUM=ZERO
            DO 25 L=1,N
   25       DUM=DUM+A(L,J)*A(L,K)
            DO 30 I=1,N
   30    A(I,J)=A(I,J)-DUM*A(I,K)
         D=ZERO
         DO 40 I=1,N
         IF (DABS(A(I,J)).LT. TOL ) A(I,J)=ZERO
   40    D=D+A(I,J)**2
C     IF(D .LT. TOL) GO TO 90
         D=DSQRT(D)
         IF(D.LT.TOL) GO TO 90
         DO 50 I=1,N
   50    A(I,J)=A(I,J)/D
   90    DD(J)=D
  100 CONTINUE
      RETURN
      END
