      SUBROUTINE TRITRN(N,KI1,KL1,X,B,C,DI)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      COMMON /TIMX/ITIME,ICONT0,JCOUNT,ITIME2,ICONT2,ITIME3,ICONT3,
     1ITIME4,ICONT4,ITIME1,ICONT1
      INTEGER START,END
      DIMENSION X(1), B(1), C(1), DI(1)
      DATA A0 /0.0D0/
C
      ICONT0=ICONT0+1
CSEL  CALL ETIME2(START)
      IF(N.NE.1) GO TO 16
      JCOUNT=JCOUNT+1
      X(1)=B(KL1+1)*C(KI1+1)*C(KI1+1)
      RETURN
  16  KI=KI1
      KL1P1=KL1+1
      IJ=0
      DO 64 I=1,N
         LIK=KI
         KL=KL1P1
         KI=KI+1
         DI(1)=C(KI)*B(KL)
C START OF CA BLOCK
C     DO 48 K=2,N
C     KM1=K-1
C     DI(K)=A0
C     KI=KI+1
C     DO 32 L=1,KM1
C     KL=KL+1
C     DI(L)=DI(L)+C(KI)*B(KL)
C32   DI(K)=DI(K)+C(LIK+L)*B(KL)
C     KL=KL+1
C  48 DI(K)=DI(K)+C(KI)*B(KL)
         DO 48 K=2,N
            KI=KI+1
            DIK=0.0D+00
            CKI=C(KI)
            DO 32 L=1,K-1
               BKL=B(KL+L)
               DI(L)=DI(L)+CKI*BKL
   32       DIK=DIK+C(LIK+L)*BKL
            KL=KL+L
   48    DI(K)=DIK+C(KI)*B(KL)
C
C
C
         LJ=KI1
C START OF CB BLOCK
C     DO 64 J=1,I
C     IJ=IJ+1
C     X(IJ)=A0
C     DO 64 L=1,N
C     LJ=LJ+1
C  64 X(IJ)=X(IJ)+DI(L)*C(LJ)
         DO 2 J=1,I
            XIJ=0.0D+00
            DO 1 L=1,N
               XIJ=XIJ+DI(L)*C(LJ+L)
    1       CONTINUE
            LJ=LJ+N
            X(IJ+J)=XIJ
    2    CONTINUE
         IJ=IJ+I
   64 CONTINUE
C
C
CSEL  CALL ETIME2(END)
      ITIME=ITIME+END-START
      RETURN
      END
