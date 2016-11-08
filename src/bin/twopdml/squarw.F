      SUBROUTINE SQUARW(CA,C,N)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C      EXTENDED DUMMY CA,C
C
C     DIMENSION CA(1),C(N,N)
C
C     SQRT2=DSQRT(2.0D+00)
C     IJ=0
C     DO 2 I=1,N
C     DO 1 J=1,I-1
C     IJ=IJ+1
C     T=CA(IJ)
C     C(I,J)=T
C     C(J,I)=-T
C   1 CONTINUE
C     IJ=IJ+1
C     C(I,I)=SQRT2*CA(IJ)
C   2 CONTINUE
C     RETURN
C
      DIMENSION CA(1),C(1)
C
      SQRT2=DSQRT(2.0D+00)
      K=1
      JISV=1
      DO 2 I=1,N
           IJ=I
           JI=JISV
           DO 1 J=1,I-1
                T=CA(K)
                C(IJ)=T
                C(JI)=T
                K=K+1
                IJ=IJ+N
                JI=JI+1
    1      CONTINUE
           C(IJ)=SQRT2*CA(K)
           K=K+1
           JISV=JISV+N
    2 CONTINUE
      RETURN
      END
