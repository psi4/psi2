      FUNCTION DOTG(A,NA,B,NB,N)
      REAL*8 A,B,D,DOTG
C
      DIMENSION A(1),B(1)
C
      IAPT=1
      IBPT=1
      D=0.0
      DO 10 I=1,N
         D=D+A(IAPT)*B(IBPT)
         IAPT=IAPT+NA
         IBPT=IBPT+NB
 10   CONTINUE
      DOTG=D
      RETURN
      END
