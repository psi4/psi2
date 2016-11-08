      FUNCTION ENUC(ZAN,C)
C***********************************************************************
C     THIS FUNCTION CALCULATES THE NUCLEAR REPULSION ENERGY.           *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /INFOA/  NAT,NUM,NQW,NX
      PARAMETER (ZERO=0.0D+00)
      DIMENSION ZAN(50),C(3,50)
      ENUC=ZERO
      IF(NAT.EQ.1) RETURN
      DO 20 I=2,NAT
        NI=I-1
        DO 20 J=1,NI
          RR=ZERO
          DO 10 K=1,3
   10     RR=RR+(C(K,I)-C(K,J))**2
   20 ENUC=ENUC+ZAN(I)*ZAN(J)/DSQRT(RR)
C
      RETURN
      END
