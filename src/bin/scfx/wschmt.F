C     //////////////////////////////////////////////////////
      SUBROUTINE WSCHMT(NQ,NDOC,QALP,RALP,IXALP,SUC,SUQ,SUR,SCRTCH,
     $                 SEX,SUND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SUC(NQ,NQ),SUQ(NQ),SUR(NQ),SCRTCH(NQ,NQ),SEX(NQ,NQ)
      DIMENSION SUND(NQ,1)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
         IF(IXALP.EQ.QALP) THEN
            IYALP=RALP
         ELSE
            IYALP=QALP
         END IF
      DO 4 J=NDOC+1,NQ
      DO 4 I=1,NQ
 4    SUC(I,J)=SUND(I,J-NDOC)
      DO 5 J=1,NQ
      DO 5 I=1,NQ
 5    SCRTCH(I,J)=SUC(I,J)
      DO 10  N=1,NDOC
         DO 20 M=1,N-1
         VN=ZERO
            DO 25  K=1,NQ
            DO 25  L=1,NQ
 25         VN=VN+SEX(K,L)*SUC(L,M)*SUC(K,N)
            DO 30  I=1,NQ
 30         SUC(I,N)=SUC(I,N)-VN*SUC(I,M)
 20      CONTINUE
         VN=ZERO
         DO 35  K=1,NQ
         DO 35  L=1,NQ
 35      VN=VN+SEX(K,L)*SUC(K,N)*SUC(L,N)
         VN=ONE/DSQRT(VN)
         DO 40  I=1,NQ
 40      SUC(I,N)=SUC(I,N)*VN
 10   CONTINUE
      DO 60  M=1,NDOC
         VN1=ZERO
         VN2=ZERO
         VN3=ZERO
         DO 65  K=1,NQ
         DO 65  L=1,NQ
         XX=SEX(K,L)*SUC(L,M)
         VN1=VN1+SUQ(K)*XX
         VN2=VN2+SUR(K)*XX
 65      VN3=VN3+SUC(K,IYALP)*XX
         DO 70  I=1,NQ
         SUQ(I)=SUQ(I)-VN1*SUC(I,M)
         SUR(I)=SUR(I)-VN2*SUC(I,M)
 70      SUC(I,IYALP)=SUC(I,IYALP)-VN3*SUC(I,M)
 60   CONTINUE
      VN1=ZERO
      VN2=ZERO
      DO 75  K=1,NQ
      DO 75  L=1,NQ
      VN1=VN1+SEX(K,L)*SUQ(K)*SUQ(L)
 75   VN2=VN2+SEX(K,L)*SUR(K)*SUR(L)
      VN1=ONE/DSQRT(VN1)
      VN2=ONE/DSQRT(VN2)
      DO 80  I=1,NQ
      SUQ(I)=SUQ(I)*VN1
      SUR(I)=SUR(I)*VN2
 80   SUC(I,IXALP)=SUR(I)
      VN1=ZERO
      VN2=ZERO
      DO 90  K=1,NQ
      DO 90  L=1,NQ
 90   VN1=VN1+SUC(L,IXALP)*SUC(K,IYALP)*SEX(K,L)
      DO 95  I=1,NQ
 95   SUC(I,IYALP)=SUC(I,IYALP)-VN1*SUC(I,IXALP)
      DO 100  K=1,NQ
      DO 100  L=1,NQ
 100  VN2=VN2+SUC(K,IYALP)*SUC(L,IYALP)*SEX(K,L)
      VN2=ONE/DSQRT(VN2)
      DO 105  I=1,NQ
 105  SUC(I,IYALP)=VN2*SUC(I,IYALP)
      DO 120  N=NDOC+1,NQ
      IF(N.EQ.IXALP.OR.N.EQ.IYALP) GO TO 120
        DO 130  M=1,NQ
        IF(M.GE.N.AND.M.NE.IXALP.AND.M.NE.IYALP)GO TO 130
        VN=ZERO
        DO 135  K=1,NQ
        DO 135  L=1,NQ
  135   VN=VN+SEX(K,L)*SUC(K,N)*SUC(L,M)
        DO 140  I=1,NQ
  140   SUC(I,N)=SUC(I,N)-VN*SUC(I,M)
  130   CONTINUE
        VN=ZERO
        DO 145  K=1,NQ
        DO 145  L=1,NQ
  145   VN=VN+SEX(K,L)*SUC(K,N)*SUC(L,N)
        VN=ONE/DSQRT(VN)
        DO 150  I=1,NQ
  150   SUC(I,N)=SUC(I,N)*VN
  120 CONTINUE
      RETURN
      END
