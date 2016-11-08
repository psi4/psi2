      SUBROUTINE SORT(A,MORD,NORD,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N),MORD(N),NORD(N)
C
      DO 101 I=1,N
      NORD(I)=I
  101 CONTINUE
      IF(N.LE.1) GO TO 201
      N1=N-1
      DO 103 I=1,N1
      IR=I+1
      DO 102 J=IR,N
      IF( DABS(A(I)) .GT. DABS(A(J)) ) GO TO 102
      AA=A(I)
      A(I)=A(J)
      A(J)=AA
      MM=MORD(I)
      MORD(I)=MORD(J)
      MORD(J)=MM
      NN=NORD(I)
      NORD(I)=NORD(J)
      NORD(J)=NN
  102 CONTINUE
  103 CONTINUE
C
  201 CONTINUE
      RETURN
      END
