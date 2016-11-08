      SUBROUTINE BRINT(A,NAD,M,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*6 LINE(21)
      DIMENSION A(NAD)
    1 FORMAT(2X,10(7X,I5))
    2 FORMAT(2X,21A6)
    3 FORMAT(2X,I2,2X,10F12.7)
    4 FORMAT(/)
      DATA LINE / 21*'------' /
      II=0
      JJ=0
  200 II=II+1
      JJ=JJ+1
      KK=10*JJ
      NN=KK+KK*(KK-1)/2
      MM=M
      IF(M.GT.KK) MM=KK
      LL=2*(MM-II+1)+1
      WRITE (IOUT,1) (I,I=II,MM)
      WRITE (IOUT,2) (LINE(I),I=1,LL)
      DO 101 I=II,M
      I1=I*(I-1)/2+II
      I2=I+I*(I-1)/2
      IF(I2.GT.NN) I2=I1+9
      WRITE (IOUT,3) I,(A(J),J=I1,I2)
  101 CONTINUE
      IF(M.LE.KK) GO TO 201
      WRITE (IOUT,4)
      II=KK
      GO TO 200
  201 RETURN
      END
