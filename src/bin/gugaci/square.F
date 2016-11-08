      SUBROUTINE SQUARE(CA,C,N)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     DIMENSION CA(1),C(N,N)
C
C     IJ=0
C     C(1,1)=0.0
C     DO 2 I=2,N
C     DO 1 J=1,I-1
C     IJ=IJ+1
C     T=CA(IJ)
C     C(I,J)=T
C     C(J,I)=-T
C   1 CONTINUE
C     C(I,I)=0.0
C   2 CONTINUE
C     RETURN
C
      DIMENSION CA(1),C(1)
C
      K=1
      JISV=1
      DO 2 I=1,N
      IJ=I
      JI=JISV
      DO 1 J=1,I-1
      T=CA(K)
      C(IJ)=T
      C(JI)=-T
      K=K+1
      IJ=IJ+N
      JI=JI+1
    1 CONTINUE
      C(IJ)=0.0D+00
      JISV=JISV+N
    2 CONTINUE
      RETURN
      END
