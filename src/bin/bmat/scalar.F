C     ///////////////////////////////////////
      REAL*8 FUNCTION SCALAR(U,V)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION U(3) , V(3)
      SCALAR = 0.0D0
      DO 10 I=1,3
          SCALAR = SCALAR + U(I)*V(I)
   10 CONTINUE
      RETURN
      END
