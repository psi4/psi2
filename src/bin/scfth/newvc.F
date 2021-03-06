      SUBROUTINE NEWVC(C,H,SMHALF,T1,T2,U,NUMSYM,NSYM,NNPBLK,
     #                  NUMBLK,EIGVAL,NUM,IOUT,TRIANG)
C
      IMPLICIT INTEGER (A-Z)
C
      REAL*8 C(NUMBLK),H(NNPBLK),SMHALF(NNPBLK),T1(*),T2(*)
      REAL*8 EIGVAL(NUM),TRIANG(*),U(NUMBLK)
      INTEGER NUMSYM(NSYM)
C
C     ----- DIAGONALIZE H -----
C
      PTTRI=1
      PTSQR=1
      PTEIG=1
      DO 5 SYM=1,NSYM
         NUMS=NUMSYM(SYM)
         IF (NUMS.LE.0) GO TO 5
         NNPS=(NUMS+1)*NUMS/2
         CALL VMOVE(TRIANG,H(PTTRI),NNPS)
         CALL RSP(NUMS,NUMS,NNPS,TRIANG,EIGVAL(PTEIG),1,U(PTSQR),T1,T2)
         CALL TRTOSQ(T1,SMHALF(PTTRI),NUMS,NNPS)
         CALL EBC(C(PTSQR),T1,U(PTSQR),NUMS,NUMS,NUMS)
         PTTRI=PTTRI+NNPS
         PTSQR=PTSQR+NUMS**2
         PTEIG=PTEIG+NUMS
    5 CONTINUE
C
C
      RETURN
      END
