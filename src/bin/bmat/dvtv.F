      REAL*8 FUNCTION DVTV(A,B,N)
      REAL*8 A(1),B(1)
C
C  DVTV = A' * B
C
      DVTV=0.0D0
      DO 10 I=1,N
      DVTV=DVTV+A(I)*B(I)
   10 CONTINUE
C
      RETURN
      END
