      SUBROUTINE ORDER(A,NORD,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N),NORD(N)
C
      DO 101 I=1,N
      NORD(I)=I
  101 CONTINUE
      IF(N.LE.1) GO TO 201
      N1=N-1
      DO 103 I=1,N1
      IR=I+1
      DO 102 J=IR,N
      IF(A(I).LT.A(J)) GO TO 102
      AA=A(I)
      A(I)=A(J)
      A(J)=AA
      NN=NORD(I)
      NORD(I)=NORD(J)
      NORD(J)=NN
  102 CONTINUE
  103 CONTINUE
  201 RETURN
C
      ENTRY REORDR(A,NORD,N)
      DO 105 I=1,N
  105 NORD(I)=I
      IF(N.LE.1) GO TO 205
      N1=N-1
      DO 107 I=1,N1
      IR=I+1
      DO 106 J=IR,N
      IF(A(I).GE.A(J)) GO TO 106
      AA=A(I)
      A(I)=A(J)
      A(J)=AA
      NN=NORD(I)
      NORD(I)=NORD(J)
      NORD(J)=NN
  106 CONTINUE
  107 CONTINUE
  205 CONTINUE
C
      RETURN
      END
