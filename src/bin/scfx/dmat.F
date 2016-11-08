C     ////////////////////////////////////////////////////////////////
      SUBROUTINE DMAT(NLAMDA,CM,SU,ICALL,D)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      COMMON /INPT/  NSYM,IPMAX,NMAT,ITMAX,NSFMAX,NTOTAL,MXCOEF,MSFRU
      COMMON/FLAGS/ IDEBUG,ISTO,IPUNCH,ISAVE,IRORDR
      COMMON /timex/ITIME,ICONT0,JCOUNT,ITIME2,ICONT2,ITIME3,ICONT3,
     1ITIME4,ICONT4,ITIME1,ICONT1
      DIMENSION D(NMAT)
      DIMENSION NLAMDA(NSYM),CM(2,MSFRU),SU(MXCOEF)
      KL=0
      MAX =0
      KK=0
      DO 100 LAM=1,NSYM
         NN=NLAMDA(LAM)
         IHBD = (NN*NN+NN)/2
         DO 40 I=1,IHBD
  40     D(KK+I)=0.0D0
         DO 90 I=1,NN
            KL=KL+1
            MIN = MAX + 1
            MAX = MAX + NN
            IF(CM(ICALL,KL).EQ.0.0D0)GO TO 90
            CM1=CM(ICALL,KL)
            K=0
            DO 80 IQ = MIN,MAX
               TMP1 = CM1 * SU(IQ)
               TMP2 = TMP1 + TMP1
               DO 80 IP = MIN, IQ
                  K = K+1
                  IF(IP .EQ. IQ)GOTO60
  70              D(KK+K)=D(KK+K)+TMP2*SU(IP)
                  GO TO 80
  60              D(KK+K)=D(KK+K)+TMP1*SU(IP)
  80        CONTINUE
  90     CONTINUE
         KK=KK + IHBD
  100 CONTINUE
      IF(IDEBUG.EQ.0) GO TO 110
      IF(ICALL.EQ.2) GO TO 4
      WRITE(ICHECK,1)(D(I),I=1,NMAT)
      RETURN
   4  WRITE(ICHECK,2)(D(I),I=1,NMAT)
  110 RETURN
   1  FORMAT(21H CLOSED-SHELL DENSITY//(7E18.8))
   2  FORMAT(21H   OPEN-SHELL DENSITY//(7E18.8))
      END
