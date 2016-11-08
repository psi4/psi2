      SUBROUTINE PRINTS(A,NAD,M,IOUT2)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      CHARACTER*1 LINE(126)
      DIMENSION A(NAD)
      DATA LINE /126*'-'/
    1 FORMAT(2X,10(7X,I5))
    2 FORMAT(2X,21A6)
    3 FORMAT(2X,I2,2X,10F12.6)
    4 FORMAT(/)
      IOUT=IOUT2
      IF (IOUT .EQ. 0) IOUT=6
      II=0
      JJ=0
  200 II=II+1
      JJ=JJ+1
      KK=10*JJ
      NN=KK*(KK+1)/2
      MM=M
      IF(M.GT.KK) MM=KK
      LL=2*(MM-II+1)+1
      WRITE(IOUT,1) (I,I=II,MM)
      WRITE(IOUT,2) (LINE(I),I=1,LL)
      DO 101 I=II,M
         I1=II+I*(I-1)/2
         I2=I*(I+1)/2
         IF(I2.GT.NN) I2=I1+9
         WRITE(IOUT,3) I,(A(J),J=I1,I2)
  101 CONTINUE
      IF(M.LE.KK) GO TO 201
      WRITE(IOUT,4)
      II=KK
      GO TO 200
  201 RETURN
      END
