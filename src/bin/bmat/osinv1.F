C     //////////////////////////////////////////////////////////////
      SUBROUTINE OSINV1 (A,N,D,TOL,L,M)
C
C     PARAMETERS   A - INPUT MATRIX , DESTROYED IN COMPUTATION AND
C                      REPLACED BY RESULTANT INVERSE
C                  N - ORDER OF MATRIX A
C                  D - RESULTANT DETERMINANT
C            L AND M - WORK VECTORS OF LENGHT N
C                TOL - IF PIVOT ELEMENT IS LESS THAN THIS PARAMETER THE
C                      MATRIX IS TAKEN FOR SINGULAR (USUALLY 1.0E-8)
C     A DETERMINANT OF ZERO INDICATES THAT THE MATRIX IS SINGULAR
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1), M(1), L(1)
      PARAMETER(ONE=1.0D0,ZERO=0.0D0)
      D=ONE
      NK=-N
      DO 180 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
      IF (DABS(BIGA)-DABS(A(IJ))) 10,20,20
   10 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
      J=L(K)
      IF(J-K) 50,50,30
   30 KI=K-N
      DO 40 I=1,N
      KI=KI+N
      HOLO=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   40 A(JI)=HOLO
   50 I=M(K)
      IF(I-K) 80,80,60
   60 JP=N*(I-1)
      DO 70 J=1,N
      JK=NK+J
      JI=JP+J
      HOLO=-A(JK)
      A(JK)=A(JI)
   70 A(JI)=HOLO
   80 IF(DABS(BIGA)-TOL) 90,100,100
   90 D=ZERO
      RETURN
  100 DO 120 I=1,N
      IF(I-K) 110,120,110
  110 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
  120 CONTINUE
      DO 150 I=1,N
      IK=NK+I
      IJ=I-N
      DO 150 J=1,N
      IJ=IJ+N
      IF(I-K) 130,150,130
  130 IF(J-K) 140,150,140
  140 KJ=IJ-I+K
      A(IJ)=A(IK)*A(KJ)+A(IJ)
  150 CONTINUE
      KJ=K-N
      DO 170 J=1,N
      KJ=KJ+N
      IF(J-K) 160,170,160
  160 A(KJ)=A(KJ)/BIGA
  170 CONTINUE
      D=D*BIGA
      A(KK)=ONE/BIGA
  180 CONTINUE
      K=N
  190 K=K-1
      IF(K) 260,260,200
  200 I=L(K)
      IF(I-K) 230,230,210
  210 JQ=N*(K-1)
      JR=N*(I-1)
      DO 220 J=1,N
      JK=JQ+J
      HOLO=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  220 A(JI)=HOLO
  230 J=M(K)
      IF (J-K) 190,190,240
  240 KI=K-N
      DO 250 I=1,N
      KI=KI+N
      HOLO=A(KI)
      JI=KI+J-K
      A(KI)=-A(JI)
  250 A(JI)=HOLO
      GO TO 190
  260 RETURN
      END
