      SUBROUTINE FNDNEG(EIGVAL,NEG,NVAR)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 EIGVAL(1)
      DATA ZERO/0.0D0/
C
C  FIND THE NUMBER OF NEGATIVE EIGENVALUES
C  THIS WILL BE RETURNED IN NEG
C
      DO 9 I=1,NVAR
      IF(EIGVAL(I).GT.ZERO) GO TO 10
    9 CONTINUE
   10 NEG=I-1
C
      RETURN
      END
