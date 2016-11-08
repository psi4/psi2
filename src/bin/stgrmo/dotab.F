      SUBROUTINE DOTAB(XAB,A,B,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N),B(N)
      DATA A00 / 0.0D+00 /
C
      XAB=A00
      DO 101 I=1,N
      XAB=XAB+A(I)*B(I)
  101 CONTINUE
C
      RETURN
      END
