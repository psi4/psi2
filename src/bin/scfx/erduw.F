C     /////////////////////////////////////////////////////////////////
      SUBROUTINE ERDUW(A,B,NA,EPSLON,CM,IRORDR)
C
C
C
C TAKEN MAINLY FROM VERSION 4, AUG., 1971, OF JACSCF, U. OF WA.
C   MATRIX DIAGONALIZATION BY THE JACOBI METHOD.
C     A = REAL SYMMETRIC MATRIX TO BE DIAGONALIZED. IT IS STORED
C       BY COLUMNS WITH ALL SUB-DIAGONAL ELEMENTS OMITTED, SO A(I,J) IS
C       STORED AS A((J*(J-1))/2+I).
C     B = MATRIX TO BE MULTIPLIED BY THE MATRIX OF EIGENVECTORS.
C     NA = DIMENSION OF THE MATRICES.
C      EPSLON IS THE CONVERGENCE CRITERION FOR OFF-DIAGONAL ELEMENTS.
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      DIMENSION CM(2,1)
      DIMENSION A(1), B(1)
      COMMON /DIM21/ IPQ(256)
      COMMON /timex/ ITIME,ICONT0,JCOUNT,ITIME2,ICONT2,ITIME3,ICONT3,
     1ITIME4,ICONT4,ITIME1,ICONT1
      COMMON /TAPES/NIN,IR,IOUT,ICHECK,IFLAG3,ITAPE,ITERML
      INTEGER START,END
      DATA A0, A1S2 /0.0D0, 0.5D0/
      DATA FUZZ /1.0D-05/
      ICONT3=ICONT3+1
      IF(NA.EQ.1) RETURN
      LOOPC=128
      NN=(NA*(NA-1))/16
      AMAX=A0
      SUM=A0
      IJ=1
      DO 24 I=1,NA
         DO 16 J=1,I
            IF(I.EQ.J) GO TO 8
            AMAX=MAX(AMAX,ABS(A(IJ)))
    8       TERM=A(IJ)*A(IJ)
            SUM=SUM+TERM+TERM
   16    IJ=IJ+1
   24 SUM=SUM-TERM
      THRSHG= DSQRT(SUM/NA)*EPSLON
      IF(THRSHG.GE.AMAX) RETURN
      THRESH=MAX(A1S2*AMAX,THRSHG)
   32 N=0
      IJ=2
      JJ=1
      DO 112 J=2,NA
         JJ=JJ+J
         JM1=J-1
         II=0
         DO 104 I=1,JM1
            II=II+I
            IF(ABS(A(IJ)).LT.THRESH) GO TO 104
            N=N+1
            SUM=A1S2*(A(JJ)+A(II))
            TERM=A1S2*(A(JJ)-A(II))
            AMAX= SIGN( DSQRT(TERM*TERM+A(IJ)*A(IJ)),TERM)
            C= DSQRT((AMAX+TERM)/(AMAX+AMAX))
            S=A(IJ)/(C*(AMAX+AMAX))
            A(II)=SUM-AMAX
            A(JJ)=SUM+AMAX
            A(IJ)=A0
            IM1=I-1
            IF(IM1) 40,56,40
   40       KI=II-I
            KJ=JJ-J
            DO 48 K=1,IM1
               KI=KI+1
               KJ=KJ+1
               TERM=C*A(KI)-S*A(KJ)
               A(KJ)=S*A(KI)+C*A(KJ)
   48       A(KI)=TERM
   56       IF(JM1.EQ.I) GO TO 72
            IP1=I+1
            IK=II+I
            KJ=IJ
            DO 64 K=IP1,JM1
               KJ=KJ+1
               TERM=C*A(IK)-S*A(KJ)
               A(KJ)=S*A(IK)+C*A(KJ)
               A(IK)=TERM
   64       IK=IK+K
   72       IF(J.EQ.NA) GO TO 88
            JP1=J+1
            IK=JJ+I
            JK=JJ+J
            DO 80 K=JP1,NA
               TERM=C*A(IK)-S*A(JK)
               A(JK)=S*A(IK)+C*A(JK)
               A(IK)=TERM
               IK=IK+K
   80       JK=JK+K
   88       KI=IM1*NA
            KJ=JM1*NA
            DO 96 K=1,NA
               KI=KI+1
               KJ=KJ+1
               TERM=C*B(KI)-S*B(KJ)
               B(KJ)=S*B(KI)+C*B(KJ)
   96       B(KI)=TERM
  104    IJ=IJ+1
  112 IJ=IJ+1
      LOOPC=LOOPC-1
      IF(LOOPC) 120,1024,120
  120 IF(N.GT.NN) GO TO 32
      IF(THRESH.EQ.THRSHG) GO TO 128
      THRESH=MAX(A1S2*THRESH,THRSHG)
      GO TO 32
  128 IF(N) 32,136,32
  136 IF(IRORDR.EQ.1) GO TO 160
      DO 150 I=1,NA
         II=IPQ(I+1)
         IF(ABS(CM(2,I)).GT.FUZZ) GO TO 150
         DO 145 J=I,NA
            JJ=IPQ(J+1)
            IF(ABS(CM(2,J)).GT.FUZZ) GO TO 145
            IF(A(II).GT.A(JJ)) THEN
               TEMP=A(II)
               A(II)=A(JJ)
               A(JJ)=TEMP
               DO 140 K=1,NA
                  KI=NA*(I-1)+K
                  KJ=NA*(J-1)+K
                  TEMP=B(KI)
                  B(KI)=B(KJ)
  140          B(KJ)=TEMP
            END IF
  145    CONTINUE
  150 CONTINUE
  160 CONTINUE
      ITIME3=ITIME3+END-START
      RETURN
 1024 WRITE (ICHECK,2048)
      STOP
 2048 FORMAT(12H ERROR ERDUW)
      END
