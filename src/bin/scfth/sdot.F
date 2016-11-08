      FUNCTION SDOT(N,A,IX,B,IY)
      REAL*8 A(IX,N),B(IY,N),SDOT
      SDOT=0.0D+00
      DO 10 I=1,N
      SDOT=SDOT+A(1,I)*B(1,I)
   10 CONTINUE
      RETURN
      END
