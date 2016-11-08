      double precision FUNCTION DOT(A,B,N)
C
      REAL*8 A(N),B(N),T,t1

C
      T=0.0D+00
      t1=0.0d0
      DO 1 I=1,N-1,2
           T=T+A(I)*B(I)
           t1=t1+a(i+1)*b(i+1)
    1 CONTINUE
      if(iand(n,1).ne.0) t=t+a(n)*b(n)
      DOT=T+t1
      RETURN
      END
