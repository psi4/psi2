C     /////////////////////////////////////////////////////////////
      SUBROUTINE AHX3(NC,NSX,L1,L2,L3,H11,H21,H31,H22,H32,H33,A,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NC,NC),X(NC,NC)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)
      DO 10  N=1,NSX
      DO 10  M=1,N
          DO 15  I=1,3
          DO 15  J=1,3
          W1=A(L1+I,M)*A(L1+J,N)
          W2=A(L2+I,M)*A(L2+J,N)
          W3=A(L3+I,M)*A(L3+J,N)
          X(M,N)=X(M,N)+W1*H11(I,J)+W2*H22(I,J)+W3*H33(I,J)
          W1=A(L2+I,M)*A(L1+J,N)+A(L1+J,M)*A(L2+I,N)
          W2=A(L3+I,M)*A(L1+J,N)+A(L1+J,M)*A(L3+I,N)
          W3=A(L3+I,M)*A(L2+J,N)+A(L2+J,M)*A(L3+I,N)
  15      X(M,N)=X(M,N)+W1*H21(I,J)+W2*H31(I,J)+W3*H32(I,J)
  10  CONTINUE
      DO 20  N=1,NSX
      DO 20  M=1,N
  20  X(N,M)=X(M,N)
      RETURN
      END
