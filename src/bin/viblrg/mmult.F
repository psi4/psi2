      SUBROUTINE MMULT(MAT1,MAT2,PROD,I,K,J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 MAT1(I,K),MAT2(K,J),PROD(I,J)
      DO 10 I1=1,I
      DO 10 J1=1,J
        PROD(I1,J1)=0.0D+00
        DO 10 K1=1,K
   10     PROD(I1,J1)=PROD(I1,J1)+MAT1(I1,K1)*MAT2(K1,J1)
      RETURN
      END
