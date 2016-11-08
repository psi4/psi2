C======================================================================-
      SUBROUTINE LIGEN3(N,NTWO,NSQUAR,A,VEC,EIG)
C     SUBROUTINE LIGEN3(N,NTWO,NSQUAR,A,VEC,EIG,
C    .                  GAMMA,BETA,BETASQ,W,IPOSV,IVPOS,IORD,P,Q)
      IMPLICIT REAL*8 (A-H,O-Z)
C        *****  A GIVENS HOUSHOLDER MATRIX DIAGONALIZATION  *****
C        *****  ROUTINE SAME AS EIGEN BUT WORKS WITH A      *****
C        *****  LINEAR ARRAY.                               *****
      COMMON/SIGNS/IA(500),IPRNT
C     DIMENSION       IA(NT2)
C     DIMENSION IN(N),IA(NT2)
      DIMENSION W(200),GAMMA(200),BETA(200),BETASQ(200)
      DIMENSION P(200),Q(200),IPOSV(400),IVPOS(400),IORD(400)
      DIMENSION A(NTWO),VEC(NSQUAR),EIG(N)
      EQUIVALENCE (IPOSV(1),GAMMA(1)),(IVPOS(1),BETA(1)),
     1 (IORD(1),BETASQ(1)),(P(1),BETA(1)),(Q(1),BETA(1))
      RHOSQ=1.D-18
      IF(N.EQ.0) GO TO 560
      N1=N-1
      N2=N-2
      GAMMA(1)=A(1)
      IF(N2) 280,270,120
  120 DO 260 NR=1,N2
      IK=IA(NR+1)+NR
      B=A(IK)
      S=0.0D+00
      DO 130 I=NR,N2
      IJ=IA(I+2)+NR
  130 S=S+A(IJ)**2
C        *****  PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION ****
      A(IK)=0.0D+00
      IF(S.LE.0.0D+00) GO TO 250
      S=S+B*B
      SGN=+1.0D+00
      IF(B.GE.0.0D+00) GO TO 160
      SGN=-1.0D+00
  160 SQRTS= DSQRT(S)
      D=SGN/(SQRTS+SQRTS)
      TEMP= DSQRT(.5D+00+B*D)
      W(NR)=TEMP
      A(IK)=TEMP
      D=D/TEMP
      B=-SGN*SQRTS
C        *****  -D- IS FACTOR OF PROPORTIONALITY. NOW       *****
C        *****  COMPUTE AND SAVE -W- VECTOR. EXTRA SINGLY   *****
C        *****  SUBSCRIPTED -W- VECTOR FOR SPEED.           *****
      DO 170 I=NR,N2
      IJ=IA(I+2)+NR
      TEMP=D*A(IJ)
      W(I+1)=TEMP
  170 A(IJ)=TEMP
C        *****  PREMULTIPLY VECTOR -W- BY MATRIX -A- TO     *****
C        *****  OBTAIN -P- VECTOR. SIMULTANEOUSLY ACCUMULATE ****
C        *****  DOT PRODUCT -WP- -- SCALR -K-.              *****
      WTAW=0.0D+00
      DO 220 I=NR,N1
      SUM=0.0D+00
      II=IA(I+1)
      DO 180 J=NR,I
      IJ=II+J+1
  180 SUM=SUM+A(IJ)*W(J)
      I1=I+1
      IF(N1.LT.I1) GO TO 210
      DO 200 J=I1,N1
      IJ=IA(J+1)+I+1
  200 SUM=SUM+A(IJ)*W(J)
  210 P(I)=SUM
  220 WTAW=WTAW+SUM*W(I)
      DO 230 I=NR,N1
  230 Q(I)=P(I)-WTAW*W(I)
C        *****  NOW FORM -PAP- MATRIX, REQUIRED PART        *****
      DO 240 J=NR,N1
      QJ=Q(J)
      WJ=W(J)
      JJ=J+1
      DO 240 I=J,N1
      IJ=IA(I+1)+JJ
  240 A(IJ)=A(IJ)-2.0D+00*(W(I)*QJ+WJ*Q(I))
  250 BETA(NR)=B
      BETASQ(NR)=B*B
      IL=IK+1
  260 GAMMA(NR+1)=A(IL)
  270 IJ=IA(N)+N-1
      B=A(IJ)
      BETA(N-1)=B
      BETASQ(N-1)=B*B
      IJ=IJ+1
      GAMMA(N)=A(IJ)
  280 BETASQ(N)=0.0D+00
C        *****  ADJOIN AN IDENTYTY MATRIX TO BE POST-       *****
C        *****  MULTIPLIED BY ROTATIONS                     *****
      NN=N*N
      DO 299 I=1,NN
  299 VEC(I)=0.0D+00
      DO 300 I=1,N
      IJ=I+(I-1)*N
  300 VEC(IJ)=1.0D+00
      M=N
      SUM=0.0D+00
      NPAS=1
      GO TO 400
  310 SUM=SUM+SHIFT
      COSA=1.0D+00
      G=GAMMA(1)-SHIFT
      PP=G
      PPBS=PP*PP+BETASQ(1)
      PPBR= DSQRT(PPBS)
      DO 370 J=1,M
      COSAP=COSA
      IF(PPBS.NE.0.0D+00) GO TO 320
      SINA=0.0D+00
      SINA2=0.0D+00
      COSA=1.0D+00
      GO TO 350
  320 SINA=BETA(J)/PPBR
      SINA2=BETASQ(J)/PPBS
      COSA=PP/PPBR
C        *****  POSTMULTIPLY IDENTITY BY -P- TRANSPOSE .    *****
      NT=J+NPAS
      IF(NT.LT.N) GO TO 330
      NT=N
  330 IJ=(J-1)*N
      IK=IJ+N
      DO 340 I=1,NT
      II=I+IJ
      IL=I+IK
      TEMP=COSA*VEC(II)+SINA*VEC(IL)
      VEC(IL)=-SINA*VEC(II)+COSA*VEC(IL)
  340 VEC(II)=TEMP
  350 DIA=GAMMA(J+1)-SHIFT
      U=SINA2*(G+DIA)
      GAMMA(J)=G+U
      G=DIA-U
      PP=DIA*COSA-SINA*COSAP*BETA(J)
      IF(J.NE.M) GO TO 360
      BETA(J)=SINA*PP
      BETASQ(J)=SINA2*PP*PP
      GO TO 380
  360 PPBS=PP*PP+BETASQ(J+1)
      PPBR= DSQRT(PPBS)
      BETA(J)=SINA*PPBR
  370 BETASQ(J)=SINA2*PPBS
  380 GAMMA(M+1)=G
C        *****  TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT ****
      NPAS=NPAS+1
      IF(BETASQ(M).GT.RHOSQ) GO TO 410
  390 EIG(M+1)=GAMMA(M+1)+SUM
  400 BETA(M)=0.0D+00
      BETASQ(M)=0.0D+00
      M=M-1
      IF(M.EQ.0) GO TO 430
      IF(BETASQ(M).LE.RHOSQ) GO TO 390
C        *****  TAKE ROOT OF CORMER 2 BY 2 NEAREST TO       *****
C        *****  LOWER DIAGONAL IN VALUE AS ESTIMATE OF      *****
C        *****  EIGENVALUE TO USE FOR SHIFT                 *****
  410 A2=GAMMA(M+1)
      R2=.5D+00*A2
      R1=.5D+00*GAMMA(M)
      R12=R1+R2
      DIF=R1-R2
      TEMP= DSQRT(DIF*DIF+BETASQ(M))
      R1=R12+TEMP
      R2=R12-TEMP
      DIF= DABS(A2-R1)- DABS(A2-R2)
      IF(DIF.LT.0.0D+00) GO TO 420
      SHIFT=R2
      GO TO 310
  420 SHIFT=R1
      GO TO 310
  430 EIG(1)=GAMMA(1)+SUM
      DO 440 J=1,N
      IPOSV(J)=J
      IVPOS(J)=J
  440  IORD(J)=J
      M=N
      GO TO 470
  450 DO 460 J=1,M
      IF(EIG(J).LE.EIG(J+1)) GO TO 460
      TEMP=EIG(J)
      EIG(J)=EIG(J+1)
      EIG(J+1)=TEMP
      ITEMP=IORD(J)
      IORD(J)=IORD(J+1)
      IORD(J+1)=ITEMP
  460 CONTINUE
  470 M=M-1
      IF(M.NE.0) GO TO 450
      IF(N1.EQ.0) GO TO 500
      DO 490 L=1,N1
      NV=IORD(L)
      NP=IPOSV(NV)
      IF(NP.EQ.L) GO TO 490
      LV=IVPOS(L)
      IVPOS(NP)=LV
      IPOSV(LV)=NP
      IL=(L-1)*N
      IK=(NP-1)*N
      DO 480 I=1,N
      II=I+IL
      JJ=I+IK
      TEMP=VEC(II)
      VEC(II)=VEC(JJ)
  480 VEC(JJ)=TEMP
  490 CONTINUE
C        *****  BACK TRANSFORM THE VECTORS OF THE TRIPLE    *****
C        *****  DIAGONAL MATRIX.                            *****
  500 DO 550 NRR=1,N
      K=N1
  510 K=K-1
      IF(K.LE.0) GO TO 550
      SUM=0.0D+00
      DO 520 I=K,N1
      IJ=IA(I+1)+K
      IK=I+1+(NRR-1)*N
  520 SUM=SUM+VEC(IK)*A(IJ)
      SUM=SUM+SUM
      DO 530 I=K,N1
      IJ=IA(I+1)+K
      IK=I+1+(NRR-1)*N
  530 VEC(IK)=VEC(IK)-SUM*A(IJ)
      GO TO 510
  550 CONTINUE
  560 CONTINUE
      RETURN
      END

