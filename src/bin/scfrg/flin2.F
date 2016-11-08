      SUBROUTINE FLIN2(A,IDIM,IN,IM,DET)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
C     LINEAR SIMULTANEOUS EQUATION
C
C     A(IN*IN) * X(IN*IM) = B(IN*IM)
C
C     A & B SHOULD BE STORED ON A(IN*(IN+IM))
C     SOLUTION X WILL BE STORED ON B PART IN DIMENSION A.
C
      DIMENSION A(IDIM,1)
C
      N=IN
      NR=IM
      JMAX=N+NR
CP    CALL MATOUT(A,N,JMAX,N,JMAX,ICHECK)
      SIGN=1.0
C M IS THE STAGE OF ELIMINATION
      DO 49 M=1,N
         TEMP=0.0
         DO 41 I=M,N
            IF(M.GT.1)A(I,M)=A(I,M)-DOTG(A(I,1),IDIM,A(1,M),1,M-1)
            IF(ABS(A(I,M)).LE.TEMP)GOTO 41
            TEMP=ABS(A(I,M))
            IMAX=I
 41      CONTINUE
         IF(TEMP.LE.0.0)GOTO 999
         IF(IMAX.EQ.M)GOTO 45
         SIGN=-SIGN
         DO 44 J=1,JMAX
            STOR=A(M,J)
            A(M,J)=A(IMAX,J)
            A(IMAX,J)=STOR
 44      CONTINUE
 45      CONTINUE
         JJ=M+1
         IF(JJ.GT.JMAX)GOTO 49
         IF(M.GT.1)GOTO 47
         DO 46 J=JJ,JMAX
            A(1,J)=A(1,J)/A(1,1)
 46      CONTINUE
         D=A(1,1)
         GOTO 49
 47      CONTINUE
         DO 48 J=JJ,JMAX
            A(M,J)=(A(M,J)-DOTG(A(M,1),IDIM,A(1,J),1,M-1))/A(M,M)
 48      CONTINUE
         D=D*A(M,M)
 49   CONTINUE
      IF(NR.EQ.0)RETURN
      DO 59 I=1,NR
         NPI=N+I
         DO 58 K=2,N
            J=N+1-K
            A(J,NPI)=A(J,NPI)-DOTG(A(J,J+1),IDIM,A(J+1,NPI),1,K-1)
 58      CONTINUE
 59   CONTINUE
      DET=D*SIGN
      RETURN
C ON ZERO PIVOT, SET DET=0.AND RETURN TO CALLING PROGRAM NOV 1972
 999  DET=0.0
      RETURN
      END
