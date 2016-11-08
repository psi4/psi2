C
C-----------------------------------------------------------------------
C
      SUBROUTINE DSINV(A,N,EPS,IER,DET)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N)
      EPS = 1.D-08
      IF(N-1) 12,1,1
   1  IER=0
      KPIV=0
      DO 11 K=1,N
      KPIV=KPIV+K
      IND=KPIV
      LEND=K-1
      TOL=DABS(EPS*A(KPIV))
      DO 11 I=K,N
      DSUM=0.D0
      IF(LEND) 2,4,2
   2  DO 3 L=1,LEND
      LANF=KPIV-L
      LIND=IND-L
   3  DSUM=DSUM+A(LANF)*A(LIND)
   4  DSUM=A(IND)-DSUM
      IF(I-K) 10,5,10
   5  IF(DSUM-TOL) 6,6,9
   6  IF(DSUM) 12,12,7
   7  IF(IER) 8,8,9
   8  IER=K-1
   9  DPIV=DSQRT(DSUM)
      A(KPIV)=DPIV
      DPIV=1.D0/DPIV
      GO TO 11
  10  A(IND)=DSUM*DPIV
  11  IND=IND+I
      GO TO 13
  12  IER=-1
  13  IF(IER) 23,14,14
  14  DET=1.D0
      DO 15 LL=1,N
      IX= (LL*(LL+1))/2
      AUX=A(IX)
  15  DET = DET*(AUX)**2
      IPIV=(N*(N+1))/2
      IND=IPIV
      DO 20 I=1,N
      DIN=1.D0/A(IPIV)
      A(IPIV)=DIN
      MIN=N
      KEND=I-1
      LANF=N-KEND
      IF(KEND) 19,19,16
  16  J=IND
      DO 18 K=1,KEND
      WORK=0.D0
      MIN=MIN-1
      LHOR=IPIV
      LVER=J
      DO 17 L=LANF,MIN
      LVER=LVER+1
      LHOR=LHOR+L
  17  WORK=WORK+A(LVER)*A(LHOR)
      A(J)=-WORK*DIN
  18  J=J-MIN
  19  IPIV=IPIV-MIN
  20  IND=IND-1
      DO 22 I=1,N
      IPIV=IPIV+I
      J=IPIV
      DO 22 K=I,N
      WORK=0.D0
      LHOR=J
      DO 21 L=K,N
      LVER=LHOR+K-I
      WORK=WORK+A(LHOR)*A(LVER)
  21  LHOR=LHOR+L
      A(J)=WORK
  22  J=J+K
  23  CONTINUE
      IF(IER) 54,56,54
   54 WRITE( 6,603) IER
   56 CONTINUE
  603 FORMAT(//,' ****** UNSUCCESFUL INVERSION  ********* ',I4,
     .        /,' ****** DETERMINANT LESS THAN ZERO ***** ')
      RETURN
      END
