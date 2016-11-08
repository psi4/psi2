      SUBROUTINE DMAT1(NLAMDA,CM,SU,ICALL,D1,D2,IOPEN,SQDIF,NEXT,A1,
     1ICOUNT)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
C
C     FIRST PART OF THIS FORMS THE DENSITY MATRICES . THE SECOND PART
C     ATTEMPS TO EXTRAPOLATE THEM
C
      integer nmat,nsym,msfru,mxcoef
      DIMENSION D1(NMAT),D2(NMAT),A1(1)
      DIMENSION NLAMDA(NSYM),CM(2,MSFRU),SU(MXCOEF)
      LOGICAL IPOPLE
      COMMON /POP/ IPOPLE
      COMMON /INPT/  NSYM,IPMAX,NMAT,ITMAX,NSFMAX,NTOTAL,MXCOEF,MSFRU
      COMMON /TAPES/ NIN,IR,IOUT,ICHECK,IFLAG3,ITAPE,I30
      COMMON/FLAGS/ IDEBUG,ISTO,IPUNCH,ISAVE,IRORDR
      COMMON /EXTR/ NMAT2,NMAT3,NMAT4,NMAT5,NMAT6,NMAT7
      COMMON/TIMX/ ITIME,NCOUNT,JCOUNT,ITIME2,ICONT2,ITIME3,ICONT3,
     1ITIME4,ICONT4,ITIME1,ICONT1
      INTEGER START,END
      ICONT2=ICONT2+1
      IF (ICALL.EQ.2) GO TO 111
C     FORM CLOSED SHELL DENSITY MATRIX
      KL=0
      MAX =0
      KK=0
      DO 100 LAM=1,NSYM
         NN=NLAMDA(LAM)
         IHBD = (NN*NN+NN)/2
         DO 40 I=1,IHBD
  40     D1(KK+I)=0.0D0
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
  70              D1(KK+K)=D1(KK+K)+TMP2*SU(IP)
                  GO TO 80
  60              D1(KK+K)=D1(KK+K)+TMP1*SU(IP)
  80        CONTINUE
  90     CONTINUE
         KK=KK + IHBD
  100 CONTINUE
      GO TO 110
C
C     FORM THE OPEN-SHELL DENSITY MATRIX
C
111   KL=0
      MAX =0
      KK=0
      DO 1100 LAM=1,NSYM
         NN=NLAMDA(LAM)
         IHBD = (NN*NN+NN)/2
         DO 140 I=1,IHBD
 140     D2(KK+I)=0.0D0
         DO 190 I=1,NN
            KL=KL+1
            MIN = MAX + 1
            MAX = MAX + NN
            IF(CM(ICALL,KL).EQ.0.0D0)GO TO 190
            CM1=CM(ICALL,KL)
            K=0
            DO 180 IQ = MIN,MAX
               TMP1 = CM1 * SU(IQ)
               TMP2 = TMP1 + TMP1
               DO 180 IP = MIN, IQ
                  K = K+1
                  IF(IP .EQ. IQ)GOTO160
 170              D2(KK+K)=D2(KK+K)+TMP2*SU(IP)
                  GO TO 180
 160              D2(KK+K)=D2(KK+K)+TMP1*SU(IP)
 180        CONTINUE
 190     CONTINUE
         KK=KK + IHBD
 1100 CONTINUE
  110 IF(IDEBUG.EQ.0) GO TO 112
      IF(ICALL.EQ.2) GO TO 4
      WRITE(ICHECK,1)(D1(I),I=1,NMAT)
      GO TO 112
   4  WRITE(ICHECK,2)(D2(I),I=1,NMAT)
C
C     TEST TO SEE IF WE WANT TO DO THE EXTRAPOLATION
C
112   IF (.NOT.IPOPLE) THEN
         RETURN
      ELSE IF (ICALL.EQ.1.AND.IOPEN.GT.0) THEN
         RETURN
      ELSE IF (SQDIF.LT.1.0D-09) THEN
         RETURN
      END IF
      ICOUNT=ICOUNT+1
      IF (ICOUNT.EQ.1) THEN
         IF (IOPEN.EQ.0) THEN
            NMAT2=2*NMAT
            NMAT3=3*NMAT
            NMAT4=4*NMAT
            NMAT5=5*NMAT
            NMAT6=6*NMAT
         ELSE
            NMAT2=2*NMAT
            NMAT3=3*NMAT
            NMAT4=4*NMAT
            NMAT5=5*NMAT
            NMAT6=6*NMAT
            NMAT7=7*NMAT
            NMAT8=8*NMAT
C     NMAT9=9*NMAT
C     NMAT10=10*NMAT
C     NMAT11=11*NMAT
C     NMAT12=12*NMAT
         END IF
      END IF
C
C     NOW STORE THE DENSITY MATRICES PROPERLY
C
      IF (ICOUNT.EQ.1) THEN
         IF (IOPEN.EQ.0) THEN
            DO 200 I=1,NMAT
               A1(I)=D1(I)
200         CONTINUE
         ELSE
            DO 210 I=1,NMAT
               A1(I)=D1(I)
               A1(NMAT+I)=D2(I)
210         CONTINUE
         END IF
      ELSE IF (ICOUNT.EQ.2) THEN
         IF (IOPEN.EQ.0) THEN
            DO 220 I=1,NMAT
               A1(NMAT+I)=D1(I)
220         CONTINUE
         ELSE
            DO 230 I=1,NMAT
               A1(NMAT2+I)=D1(I)
               A1(NMAT3+I)=D2(I)
230         CONTINUE
         END IF
      ELSE IF (ICOUNT.EQ.3) THEN
         IF (IOPEN.EQ.0) THEN
            DO 240 I=1,NMAT
               A1(NMAT2+I)=D1(I)
240         CONTINUE
         ELSE
            DO 250 I=1,NMAT
               A1(NMAT4+I)=D1(I)
               A1(NMAT5+I)=D2(I)
250         CONTINUE
         END IF
      ELSE IF (ICOUNT.EQ.4) THEN
         IF (IOPEN.EQ.0) THEN
            DO 260 I=1,NMAT
               A1(NMAT3+I)=D1(I)
260         CONTINUE
         ELSE
            DO 270 I=1,NMAT
               A1(NMAT6+I)=D1(I)
               A1(NMAT7+I)=D2(I)
270         CONTINUE
         END IF
      ELSE
         IF (IOPEN.EQ.0) THEN
            DO 280 I=1,NMAT
               A1(I)      =A1(NMAT+I)
               A1(NMAT+I) =A1(NMAT2+I)
               A1(NMAT2+I)=A1(NMAT3+I)
               A1(NMAT3+I)=D1(I)
280         CONTINUE
         ELSE IF (ICOUNT .GT. 4) THEN
            DO 290 I=1,NMAT2
               A1(I)=A1(NMAT2+I)
               A1(NMAT2+I)=A1(NMAT4+I)
               A1(NMAT4+I)=A1(NMAT6+I)
290         CONTINUE
            DO 300 I=1,NMAT
               A1(NMAT6+I)=D1(I)
               A1(NMAT7+I)=D2(I)
300         CONTINUE
         END IF
      END IF
C
C     LETS FORM THE DENSITY MATRICES DIFFERENCE BUT WE ONLY EXTRAPOLTE A
C     4 ITERATIONS OR MORE DEPENDING ON THE LINEARITY OF THE DENSITY MAT
C     DIFFERENCES
C
      IF (ICOUNT.LT.4) RETURN
C     IF (IOPEN.EQ.0) THEN
C     DO 305 I=1,NMAT
C     A1(NMAT4+I)=A1(NMAT3+I)-A1(NMAT2+I)
C     A1(NMAT5+I)=A1(NMAT2+I)-A1(NMAT+I)
C     A1(NMAT6+I)=A1(NMAT+I )-A1(I)
C305  CONTINUE
C     ELSE
C     DO 310 I=1,NMAT2
C     A1(NMAT8+I)=A1(NMAT6+I)-A1(NMAT4+I)
C     A1(NMAT10+I)=A1(NMAT4+I)-A1(NMAT2+I)
C     A1(NMAT12+I)=A1(NMAT2+I)-A1(I)
C310  CONTINUE
C     END IF
C
C     NOW LETS MULTIPLY THESE DIFFERENCES TOGETHER
C
      IF (IOPEN.EQ.0) THEN
         SP11=TRACLO(A1(NMAT3+1),A1(NMAT2+1),A1(NMAT3+1),A1(NMAT2+1),
     *   ICALL,NLAMDA)
         SP12=TRACLO(A1(NMAT2+1),A1(NMAT+1),A1(NMAT3+1),A1(NMAT2+1),
     *   ICALL,NLAMDA)
         SP13=TRACLO(A1(NMAT+1),A1(1),A1(NMAT3+1),A1(NMAT2+1),ICALL,
     *   NLAMDA)
         SP22=TRACLO(A1(NMAT2+1),A1(NMAT+1),A1(NMAT2+1),A1(NMAT+1),
     *   ICALL,NLAMDA)
         SP23=TRACLO(A1(NMAT+1),A1(1),A1(NMAT2+1),A1(NMAT+1),ICALL,
     *   NLAMDA)
         SP33=TRACLO(A1(NMAT+1),A1(1),A1(NMAT+1),A1(1),ICALL,NLAMDA)
      ELSE
         SP11=TRACLO(A1(NMAT6+1),A1(NMAT4+1),A1(NMAT6+1),A1(NMAT4+1),
     *   ICALL,NLAMDA)
         SP12=TRACLO(A1(NMAT4+1),A1(NMAT2+1),A1(NMAT6+1),A1(NMAT4+1),
     *   ICALL,NLAMDA)
         SP13=TRACLO(A1(NMAT2+1),A1(1),A1(NMAT6+1),A1(NMAT4+1),ICALL,
     *   NLAMDA)
         SP22=TRACLO(A1(NMAT4+1),A1(NMAT2+1),A1(NMAT4+1),A1(NMAT2+1),
     *   ICALL,NLAMDA)
         SP23=TRACLO(A1(NMAT2+1),A1(1),A1(NMAT4+1),A1(NMAT2+1),ICALL,
     *   NLAMDA)
         SP33=TRACLO(A1(NMAT2+1),A1(1),A1(NMAT2+1),A1(1),ICALL,NLAMDA)
      END IF
      IF (IOPEN.EQ.0) THEN
         DP1=DSQRT(SP11)
         DP2=DSQRT(SP22)
         DP3=DSQRT(SP33)
      ELSE
         DP1=DSQRT(SP11/2.0D0)
         DP2=DSQRT(SP22/2.0D0)
         DP3=DSQRT(SP33/2.0D0)
      END IF
C
C     FIND COSINE OF ANGLE BETWEEN DISPLACEMENTS
C
      IF (IOPEN.EQ.0) THEN
         COSPHI=SP12/(DP1*DP2)
      ELSE
         COSPHI=SP12/(2.0D0*DP1*DP2)
      END IF
C
C     FIND COSINE OF ANGLE DP(3) AND PLANE OF DP(1) AND DP(2)
C
      P=SP11*SP22-SP12*SP12
      X=(SP13*SP22-SP12*SP23)/(SP11*SP22-SP12*SP12)
      Y=(SP23*SP11-SP12*SP13)/(SP11*SP22-SP12*SP12)
      IF (IOPEN.EQ.0) THEN
         COSPSI=DSQRT(X*X*SP11+Y*Y*SP22+2.0D0*X*Y*SP12)/DP3
      ELSE
         COSPSI=DSQRT((X*X*SP11+Y*Y*SP22+2.0D0*X*Y*SP12)/2.0D0)/DP3
      END IF
C
C     TEST TO SEE IF EXTRAPOLATION IS POSSIBLE
C
      IF ((COSPSI-0.99).LT.0.0) RETURN
      Y=-Y/X
      X=1.0/X
      XY=Y*Y+4.0*X
      IF (XY.LT.0) RETURN
      XY=ABS(Y)+DSQRT(XY)
C
C     TEST IF WE DO A THREE POINT OR FOUR POINT EXTRAPOLATION
C
      IF (XY.LT.1.9) THEN
         XXX=X/(1.0-X-Y)
         YYY=(X+Y)/(1.0-X-Y)
         IF (IOPEN.EQ.0) THEN
            DO 350 I=1,NMAT
               D1(I)=D1(I)+XXX*(A1(NMAT2+I)-A1(NMAT+I))+YYY*(A1(
     *         NMAT3+I)-A1(NMAT2+I))
350         CONTINUE
         ELSE
            DO 360 I=1,NMAT
               D1(I)=D1(I)+XXX*(A1(NMAT4+I)-A1(NMAT2+I))+YYY*(A1(
     *         NMAT6+I)-A1(NMAT4+I))
               D2(I)=D2(I)+XXX*(A1(NMAT5+I)-A1(NMAT3+I))+YYY*(A1(
     *         NMAT7+I)-A1(NMAT5+I))
360         CONTINUE
         END IF
         ICOUNT=0
      ELSE
         IF ((ABS(COSPHI)-0.995).LE.0.0) RETURN
         X=DP1/(DP2*COSPHI-DP1)
         IF (IOPEN.EQ.0) THEN
            DO 370 I=1,NMAT
               D1(I)=D1(I)+X*(A1(NMAT3+I)-A1(NMAT2+I))
370         CONTINUE
            ICOUNT=0
         ELSE
            DO 380 I=1,NMAT
               D1(I)=D1(I)+X*(A1(NMAT6+I)-A1(NMAT4+I))
               D2(I)=D2(I)+X*(A1(NMAT7+I)-A1(NMAT5+I))
380         CONTINUE
            ICOUNT=0
         END IF
      END IF
C
      WRITE (IOUT,3)
      IF (IFLAG3.GT.0) WRITE (IFLAG3,3)
C
   1  FORMAT(21H CLOSED-SHELL DENSITY//(7E18.8))
   2  FORMAT(21H   OPEN-SHELL DENSITY//(7E18.8))
   3  FORMAT (' EXTRAPOLATION OF THE DENSITY MATRICES')
      RETURN
      END
