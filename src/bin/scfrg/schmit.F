      SUBROUTINE SCHMIT(NSYM,NV,ISWTCH,THRSH,NLAMDA,SU,S,V)
C
C
C
C  ISWTCH.EQ.1   ORTHONORMALIZE ALL VECTORS
C  ISWTCH.EQ.0   ORTHONORMALIZE OCCUPIED VECTORS ONLY
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      COMMON /TAPES/ NIN,IR,IOUT,ICHECK,IFLAG3,ITAPE,I30
      COMMON /TIMX/ITIME,ICONT0,JCOUNT,ITIME2,ICONT2,ITIME3,ICONT3,
     1ITIME4,ICONT4,ITIME1,ICONT1
      INTEGER START,END
      DIMENSION NV(2,1),V(1),NLAMDA(1),S(1),SU(1)
      DATA A0,A1/0.0D0,1.0D0/
      ICONT1=ICONT1+1
      KI=0
      KL=0
C
C     FOR I=1,NSYM
C     PRINT '(' IRREP #',I2,' HAS',I3,' MOS')', I,NLAMDA(I)
C     END FOR
C
      DO 64 IIRP=1,NSYM
         KLIIRP=KL
         N=NLAMDA(IIRP)
         IF(ISWTCH.EQ.0) GO TO 12
         NNV=N
         GO TO 16
   12    NNV=NV(1,IIRP)+NV(2,IIRP)
         IF(NNV.EQ.0) NNV=1
   16    DO 60 I=1,NNV
            LIK=KI
            KI=KI+1
            KL=KLIIRP+1
            V(1)=SU(KI)*S(KL)
            IF(N.EQ.1) GO TO 28
C     DO 24 K=2,N
C     V(K)=A0
C     KI=KI+1
C     KM1=K-1
C     DO 20 L=1,KM1
C     KL=KL+1
C     V(K)=V(K)+SU(LIK+L)*S(KL)
C  20 V(L)=V(L)+SU(KI)*S(KL)
C     KL=KL+1
C  24 V(K)=V(K)+SU(KI)*S(KL)
            DO 2 K=2,N
               KI=KI+1
               SUKI=SU(KI)
               VK=0.0
               DO 1 L=1,K-1
                  SKL=S(KL+L)
                  VK=VK+SU(LIK+L)*SKL
                  V(L)=V(L)+SUKI*SKL
    1          CONTINUE
               KL=KL+K
               V(K)=VK+SU(KI)*S(KL)
    2       CONTINUE
C
C
C
C SECTION TO NORMALIZE  VECTORS
   28       TEMP=A0
            DO 32 L=1,N
   32       TEMP=TEMP+V(L)*SU(LIK+L)
            IF(TEMP.GT.THRSH) GO TO 40
            IF(TEMP.GT.A0) GO TO 36
            WRITE (IOUT,1024) I, IIRP, TEMP
            IF(IFLAG3.GT.0) WRITE(ICHECK,1024) I,IIRP,TEMP
            STOP
   36       WRITE (IOUT,1032) I, IIRP, TEMP
            IF(IFLAG3.GT.0) WRITE(ICHECK,1032) I, IIRP, TEMP
            TEMP=THRSH
   40       TEMP=A1/ DSQRT(TEMP)
C
C
            DO 44 L=1,N
               V(L)=TEMP*V(L)
               SU(LIK+L)=TEMP*SU(LIK+L)
   44       CONTINUE
            IF(I.EQ.NNV) GO TO 60
C SECTION TO ORTHOGONALIZE VECTORS
            IP1=I+1
            LJ=KI
            DO 56 J=IP1,NNV
               TEMP=A0
               DO 52 L=1,N
   52          TEMP=TEMP+V(L)*SU(LJ+L)
               DO 56 L=1,N
                  LJ=LJ+1
                  SU(LJ)=SU(LJ)-TEMP*SU(LIK+L)
   56       CONTINUE
   60    CONTINUE
   64 KI=KI+(N-NNV)*N
CSEL  CALL ETIME2(END)
      ITIME1=ITIME1+END-START
      RETURN
 1024 FORMAT(3H MO,I4,3H OF,I3,28H IRREP CANNOT BE NORMALIZED,,E18.8/)
 1032 FORMAT(20H NORMALIZATION OF MO,I4,3H OF,I3,23H IRREP BELOW THRESHO
     XLD,,E18.8/)
      END
