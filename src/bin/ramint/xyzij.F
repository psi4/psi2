      SUBROUTINE XYZIJ
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/EXPNT/EIIA,EJIA,EKIB,EIJIA,EKLIB,RHOTSQ,EI,EJ,EK,A,B
      COMMON/GEOIN/XIMXJ,YIMYJ,ZIMZJ,XKMXL,YKMYL,ZKMZL,
     1             XAMXB,YAMYB,ZAMZB
      COMMON/SETD /BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,
     1             NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,IOFFST
      COMMON/TWOD0/XIN(2304),YIN(2304),ZIN(2304)
      COMMON/TWODI/DXI(2304),DYI(2304),DZI(2304)
      COMMON/TWODJ/DXJ(2304),DYJ(2304),DZJ(2304)
      DATA ZERO,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /
C
      RTIASQ = RHOTSQ / (A*A)
      RTIAB  = RHOTSQ / (A*B)
      DIC00  = EI/A - RTIASQ*EI - 1
      DJC00  = EJ/A - RTIASQ*EJ
      DICP00 = EI * RTIAB
      DJCP00 = EJ * RTIAB
      RHO2T  = RHOTSQ + RHOTSQ
      RHO2TI = RHO2T*EIIA
      RHO2TJ = RHO2T*EJIA
      XICNST =  -XIMXJ*EIJIA*TWO - XAMXB*RHO2TI
      YICNST =  -YIMYJ*EIJIA*TWO - YAMYB*RHO2TI
      ZICNST = (-ZIMZJ*EIJIA*TWO - ZAMZB*RHO2TI)
      XJCNST =   XIMXJ*EIJIA*TWO - XAMXB*RHO2TJ
      YJCNST =   YIMYJ*EIJIA*TWO - YAMYB*RHO2TJ
      ZJCNST = ( ZIMZJ*EIJIA*TWO - ZAMZB*RHO2TJ)
C
C
C     ----- I(0,0) -----
C
      I = IOFFST + 1
      XIN(I) = ONE
      YIN(I) = ONE
      ZIN(I) = F00
      DXI(I) = ZERO
      DYI(I) = ZERO
      DZI(I) = ZERO
      DXJ(I) = ZERO
      DYJ(I) = ZERO
      DZJ(I) = ZERO
C
C     ----- I(1,0)=C00 -----
C
      IF (NMAX.NE.0) THEN
        I2 = IOFFST + 2
        XIN(I2) = XC00
        YIN(I2) = YC00
        ZIN(I2) = ZC00 * F00
        DXI(I2) = DIC00
        DYI(I2) = DIC00
        DZI(I2) = DIC00 * F00
        DXJ(I2) = DJC00
        DYJ(I2) = DJC00
        DZJ(I2) = DJC00 * F00
      END IF
C
C     ----- I(0,1)=CP00 -----
C
      IF (MMAX.NE.0) THEN
        I = IOFFST + 33
        XIN(I) = XCP00
        YIN(I) = YCP00
        ZIN(I) = ZCP00 * F00
        DXI(I) = DICP00
        DYI(I) = DICP00
        DZI(I) = DICP00 * F00
        DXJ(I) = DJCP00
        DYJ(I) = DJCP00
        DZJ(I) = DJCP00 * F00
C
C     ----- I(1,1)=B00 * I(0,0) + CP00 * I(1,0) -----
C
        IF (NMAX.NE.0) THEN
          I = IOFFST + 34
          XIN(I) = XCP00 * XIN(I2) + B00
          YIN(I) = YCP00 * YIN(I2) + B00
          ZIN(I) = ZCP00 * ZIN(I2) + B00 * F00
          DXI(I) = DICP00 * XIN(I2) + XCP00 * DXI(I2)
          DYI(I) = DICP00 * YIN(I2) + YCP00 * DYI(I2)
          DZI(I) = DICP00 * ZIN(I2) + ZCP00 * DZI(I2)
          DXJ(I) = DJCP00 * XIN(I2) + XCP00 * DXJ(I2)
          DYJ(I) = DJCP00 * YIN(I2) + YCP00 * DYJ(I2)
          DZJ(I) = DJCP00 * ZIN(I2) + ZCP00 * DZJ(I2)
        END IF
      END IF
C
C     ----- I(N,0)=(N-1)B10 * I(N-2,0) + C00 * I(N-1,0) -----
C
      IF (NMAX.GT.1) THEN
        DO 1 N=2,NMAX
          NL10=N + IOFFST
          N0=NL10+1
          NL20=NL10 - 1
          B10NL1=(N-1) * B10
          XIN(N0) = B10NL1 * XIN(NL20) + XC00 * XIN(NL10)
          YIN(N0) = B10NL1 * YIN(NL20) + YC00 * YIN(NL10)
          ZIN(N0) = B10NL1 * ZIN(NL20) + ZC00 * ZIN(NL10)
          DXI(N0) = B10NL1 * DXI(NL20) + DIC00 * XIN(NL10) + XC00* DXI
     *    (NL10)
          DYI(N0) = B10NL1 * DYI(NL20) + DIC00 * YIN(NL10) + YC00* DYI
     *    (NL10)
          DZI(N0) = B10NL1 * DZI(NL20) + DIC00 * ZIN(NL10) + ZC00* DZI
     *    (NL10)
          DXJ(N0) = B10NL1 * DXJ(NL20) + DJC00 * XIN(NL10) + XC00* DXJ
     *    (NL10)
          DYJ(N0) = B10NL1 * DYJ(NL20) + DJC00 * YIN(NL10) + YC00* DYJ
     *    (NL10)
          DZJ(N0) = B10NL1 * DZJ(NL20) + DJC00 * ZIN(NL10) + ZC00* DZJ
     *    (NL10)
C
C     ----- I(N,1)=N * B00 * I(N-1,0) + CP00 * I(N,0) -----
C
          IF (MMAX.NE.0) THEN
            N1=N0+32
            B00N = B00 * N
            XIN(N1) = B00N * XIN(NL10) + XCP00 * XIN(N0)
            YIN(N1) = B00N * YIN(NL10) + YCP00 * YIN(N0)
            ZIN(N1) = B00N * ZIN(NL10) + ZCP00 * ZIN(N0)
            DXI(N1) = B00N * DXI(NL10) + DICP00 * XIN(N0) + XCP00 *
     *      DXI(N0)
            DYI(N1) = B00N * DYI(NL10) + DICP00 * YIN(N0) + YCP00 *
     *      DYI(N0)
            DZI(N1) = B00N * DZI(NL10) + DICP00 * ZIN(N0) + ZCP00 *
     *      DZI(N0)
            DXJ(N1) = B00N * DXJ(NL10) + DJCP00 * XIN(N0) + XCP00 *
     *      DXJ(N0)
            DYJ(N1) = B00N * DYJ(NL10) + DJCP00 * YIN(N0) + YCP00 *
     *      DYJ(N0)
            DZJ(N1) = B00N * DZJ(NL10) + DJCP00 * ZIN(N0) + ZCP00 *
     *      DZJ(N0)
          END IF
    1   CONTINUE
      END IF
C
C     ----- I(0,M)=(M-1)BP01 * I(0,M-2) + CP00 * I(0,M-1) -----
C
      IF (MMAX.GT.1) THEN
        DO 2 M=2,MMAX
          I0M = M*32 +  1 + IOFFST
          I0ML1 = I0M - 32
          I0ML2 = I0ML1 - 32
          I1M   = I0M + 1
          BPML1 = (M-1) * BP01
          XIN(I0M) = BPML1 * XIN(I0ML2) + XCP00 * XIN(I0ML1)
          YIN(I0M) = BPML1 * YIN(I0ML2) + YCP00 * YIN(I0ML1)
          ZIN(I0M) = BPML1 * ZIN(I0ML2) + ZCP00 * ZIN(I0ML1)
          DXI(I0M) = BPML1 * DXI(I0ML2) + DICP00 * XIN(I0ML1)
     1    +                                XCP00 * DXI(I0ML1)
          DYI(I0M) = BPML1 * DYI(I0ML2) + DICP00 * YIN(I0ML1)
     1    +                                YCP00 * DYI(I0ML1)
          DZI(I0M) = BPML1 * DZI(I0ML2) + DICP00 * ZIN(I0ML1)
     1    +                                ZCP00 * DZI(I0ML1)
          DXJ(I0M) = BPML1 * DXJ(I0ML2) + DJCP00 * XIN(I0ML1)
     1    +                                XCP00 * DXJ(I0ML1)
          DYJ(I0M) = BPML1 * DYJ(I0ML2) + DJCP00 * YIN(I0ML1)
     1    +                                YCP00 * DYJ(I0ML1)
          DZJ(I0M) = BPML1 * DZJ(I0ML2) + DJCP00 * ZIN(I0ML1)
     1    +                                ZCP00 * DZJ(I0ML1)
C
C     ----- I(1,M)=M * B00 * I(0,M-1) + C00 * I(0,M) -----
C
          IF (NMAX.NE.0) THEN
            B00M = M * B00
            XIN(I1M) = B00M * XIN(I0ML1) + XC00 * XIN(I0M)
            YIN(I1M) = B00M * YIN(I0ML1) + YC00 * YIN(I0M)
            ZIN(I1M) = B00M * ZIN(I0ML1) + ZC00 * ZIN(I0M)
            DXI(I1M) = B00M * DXI(I0ML1) + DIC00 * XIN(I0M) + XC00 *
     *      DXI(I0M)
            DYI(I1M) = B00M * DYI(I0ML1) + DIC00 * YIN(I0M) + YC00 *
     *      DYI(I0M)
            DZI(I1M) = B00M * DZI(I0ML1) + DIC00 * ZIN(I0M) + ZC00 *
     *      DZI(I0M)
            DXJ(I1M) = B00M * DXJ(I0ML1) + DJC00 * XIN(I0M) + XC00 *
     *      DXJ(I0M)
            DYJ(I1M) = B00M * DYJ(I0ML1) + DJC00 * YIN(I0M) + YC00 *
     *      DYJ(I0M)
            DZJ(I1M) = B00M * DZJ(I0ML1) + DJC00 * ZIN(I0M) + ZC00 *
     *      DZJ(I0M)
          END IF
    2   CONTINUE
      END IF
C
C     ----- I(N,M)=(N-1)B10 * I(N-2,M) + M * B00 * I(N-1,M-1) +
C                                            C00 * I(N-1,M)   -----
C
      IF (NMAX.GT.1.AND.MMAX.GT.1) THEN
        DO 4 M=2,MMAX
          M8=M*32 + IOFFST
          DO 3 N=2,NMAX
            NM=N + 1 + M8
            NL1M=NM - 1
            NL1ML1= NL1M - 32
            NL2M = NL1M - 1
            B10NL1 = (N-1) * B10
            B00M = M * B00
            XIN(NM) = B10NL1 * XIN(NL2M) + B00M * XIN(NL1ML1) +XC00 *
     *      XIN(NL1M)
            YIN(NM) = B10NL1 * YIN(NL2M) + B00M * YIN(NL1ML1) +YC00 *
     *      YIN(NL1M)
            ZIN(NM) = B10NL1 * ZIN(NL2M) + B00M * ZIN(NL1ML1) +ZC00 *
     *      ZIN(NL1M)
            DXI(NM) = B10NL1 * DXI(NL2M) + B00M * DXI(NL1ML1)
     1      +          DIC00 * XIN(NL1M) + XC00  * DXI(NL1M)
            DYI(NM) = B10NL1 * DYI(NL2M) + B00M * DYI(NL1ML1)
     1      +          DIC00 * YIN(NL1M) + YC00  * DYI(NL1M)
            DZI(NM) = B10NL1 * DZI(NL2M) + B00M * DZI(NL1ML1)
     1      +          DIC00 * ZIN(NL1M) + ZC00  * DZI(NL1M)
            DXJ(NM) = B10NL1 * DXJ(NL2M) + B00M * DXJ(NL1ML1)
     1      +          DJC00 * XIN(NL1M) + XC00  * DXJ(NL1M)
            DYJ(NM) = B10NL1 * DYJ(NL2M) + B00M * DYJ(NL1ML1)
     1      +          DJC00 * YIN(NL1M) + YC00  * DYJ(NL1M)
            DZJ(NM) = B10NL1 * DZJ(NL2M) + B00M * DZJ(NL1ML1)
     1      +          DJC00 * ZIN(NL1M) + ZC00  * DZJ(NL1M)
    3     CONTINUE
    4   CONTINUE
      END IF
C
C     ----- I(I,J,M)=I(I+1,J-1,M) + (QI-QJ) * I(I,J-1,M) -----
C
      DO 7 M=0,MMAX
        MPT = M * 32 + IOFFST
        DO 6 J=1,NJMAX
          JM = J * 8 + MPT
          DO 5 I=1,NIMAX + NJMAX - J + 1
            IJM = I + JM
            IJL1M = IJM - 8
            I1JL1M = IJL1M + 1
            XIN(IJM) = XIN(I1JL1M) + XIMXJ * XIN(IJL1M)
            YIN(IJM) = YIN(I1JL1M) + YIMYJ * YIN(IJL1M)
            ZIN(IJM) = ZIN(I1JL1M) + ZIMZJ * ZIN(IJL1M)
            DXI(IJM) = DXI(I1JL1M) + XIMXJ * DXI(IJL1M) + XIN(IJL1M)
            DYI(IJM) = DYI(I1JL1M) + YIMYJ * DYI(IJL1M) + YIN(IJL1M)
            DZI(IJM) = DZI(I1JL1M) + ZIMZJ * DZI(IJL1M) + ZIN(IJL1M)
            DXJ(IJM) = DXJ(I1JL1M) + XIMXJ * DXJ(IJL1M) - XIN(IJL1M)
            DYJ(IJM) = DYJ(I1JL1M) + YIMYJ * DYJ(IJL1M) - YIN(IJL1M)
            DZJ(IJM) = DZJ(I1JL1M) + ZIMZJ * DZJ(IJL1M) - ZIN(IJL1M)
    5     CONTINUE
    6   CONTINUE
    7 CONTINUE
C
C     ----- I(I,J,M,0) = I(I,J,M) -----
C
      DO 10 M=MMAX,0,-1
        MNEW = M * 64 + IOFFST
        MOLD = M * 32 + IOFFST
        DO 9 J=NJMAX,0,-1
          JMNEW = MNEW + J * 16
          JMOLD = MOLD + J * 8
          DO 8 I=NIMAX+1,1,-1
            IJMNEW = JMNEW + I
            IJMOLD = JMOLD + I
            XIN(IJMNEW) = XIN(IJMOLD)
            YIN(IJMNEW) = YIN(IJMOLD)
            ZIN(IJMNEW) = ZIN(IJMOLD)
            DXI(IJMNEW) = DXI(IJMOLD)
            DYI(IJMNEW) = DYI(IJMOLD)
            DZI(IJMNEW) = DZI(IJMOLD)
            DXJ(IJMNEW) = DXJ(IJMOLD)
            DYJ(IJMNEW) = DYJ(IJMOLD)
            DZJ(IJMNEW) = DZJ(IJMOLD)
    8     CONTINUE
    9   CONTINUE
   10 CONTINUE
C
C     ----- I(I,J,K,L)=I(I,J,K+1,L-1) + (QK-QL) * I(I,J,K,L-1) -----
C
      DO 14 L=1,NLMAX
        LPT = 4 * L + IOFFST
        DO 13 K=0,NKMAX + NLMAX - L
          KL = K * 64 + LPT
          DO 12 J=0,NJMAX
            JKL = J * 16 + KL
            DO 11 I=1,NIMAX+1
              IJKL = I + JKL
              IJKLM1 = IJKL - 4
              KP1LM1 = IJKLM1 + 64
              XIN(IJKL) = XIN(KP1LM1) + XKMXL * XIN(IJKLM1)
              YIN(IJKL) = YIN(KP1LM1) + YKMYL * YIN(IJKLM1)
              ZIN(IJKL) = ZIN(KP1LM1) + ZKMZL * ZIN(IJKLM1)
              DXI(IJKL) = DXI(KP1LM1) + XKMXL * DXI(IJKLM1)
              DYI(IJKL) = DYI(KP1LM1) + YKMYL * DYI(IJKLM1)
              DZI(IJKL) = DZI(KP1LM1) + ZKMZL * DZI(IJKLM1)
              DXJ(IJKL) = DXJ(KP1LM1) + XKMXL * DXJ(IJKLM1)
              DYJ(IJKL) = DYJ(KP1LM1) + YKMYL * DYJ(IJKLM1)
              DZJ(IJKL) = DZJ(KP1LM1) + ZKMZL * DZJ(IJKLM1)
   11       CONTINUE
   12     CONTINUE
   13   CONTINUE
   14 CONTINUE
C
C     ----- ADD DERIVATIVE OF EXP(-D) AND EXP(-G) -----
C
      DO 18 L=0,NLMAX
        LPT = L*4 + IOFFST
        DO 17 K=0,NKMAX
          KL = K*64 + LPT
          DO 16 J=0,NJMAX
            JKL = J*16 + KL
            DO 15 I=1,NIMAX+1
              IJKL = I + JKL
              DXI(IJKL) = DXI(IJKL) + XICNST * XIN(IJKL)
              DYI(IJKL) = DYI(IJKL) + YICNST * YIN(IJKL)
              DZI(IJKL) = DZI(IJKL) + ZICNST * ZIN(IJKL)
              DXJ(IJKL) = DXJ(IJKL) + XJCNST * XIN(IJKL)
              DYJ(IJKL) = DYJ(IJKL) + YJCNST * YIN(IJKL)
              DZJ(IJKL) = DZJ(IJKL) + ZJCNST * ZIN(IJKL)
   15       CONTINUE
   16     CONTINUE
   17   CONTINUE
   18 CONTINUE
      RETURN
      END
