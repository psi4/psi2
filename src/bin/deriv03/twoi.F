      SUBROUTINE TWOI
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/EXPNT/EIIA,EJIA,EKIB,EIJIA,EKLIB,RHOTSQ,EI,EJ,EK,A,B
      COMMON/GEOIN/XIMXJ,YIMYJ,ZIMZJ,XKMXL,YKMYL,ZKMZL,
     1             XAMXB,YAMYB,ZAMZB
      COMMON/SETD /BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,
     1             NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,IOFFST
      COMMON/TOLER/DLIM
ctph  COMMON/TWOD0/XIN(2304),YIN(2304),ZIN(2304)
ctph  COMMON/TWODI/DXI(2304),DYI(2304),DZI(2304)
ctph  COMMON/TWDII/DXIXI(2304),DYIYI(2304),DZIZI(2304)
      COMMON/TWOD0/XIN(6250),YIN(6250),ZIN(6250)
      COMMON/TWODI/DXI(6250),DYI(6250),DZI(6250)
      COMMON/TWDII/DXIXI(6250),DYIYI(6250),DZIZI(6250)
      DATA ZERO,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /
C
      RTIASQ = RHOTSQ / (A*A)
      RTIBSQ = RHOTSQ / (B*B)
      RTIAB  = RHOTSQ / (A*B)
      DIC00  = EI/A - RTIASQ*EI - 1
      DICP00 = EI * RTIAB
      RHO2T  = RHOTSQ + RHOTSQ
      RHO2TI = RHO2T*EIIA
      XICNST =  -XIMXJ*EIJIA*TWO - XAMXB*RHO2TI
      YICNST =  -YIMYJ*EIJIA*TWO - YAMYB*RHO2TI
      ZICNST = (-ZIMZJ*EIJIA*TWO - ZAMZB*RHO2TI)
      DIDICN = -EIJIA*TWO - EIIA * RHO2TI
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
      DXIXI(I) = ZERO
      DYIYI(I) = ZERO
      DZIZI(I) = ZERO
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
        DXIXI(I2) = ZERO
        DYIYI(I2) = ZERO
        DZIZI(I2) = ZERO
      END IF
C
C     ----- I(0,1)=CP00 -----
C
      IF (MMAX.NE.0) THEN
ctph    I = IOFFST + 33
        I = IOFFST + 51
        XIN(I) = XCP00
        YIN(I) = YCP00
        ZIN(I) = ZCP00 * F00
        DXI(I) = DICP00
        DYI(I) = DICP00
        DZI(I) = DICP00 * F00
        DXIXI(I) = ZERO
        DYIYI(I) = ZERO
        DZIZI(I) = ZERO
C
C     ----- I(1,1)=B00 * I(0,0) + CP00 * I(1,0) -----
C
        IF (NMAX.NE.0) THEN
ctph      I = IOFFST + 34
          I = IOFFST + 52
          XIN(I) = XCP00 * XIN(I2) + B00
          YIN(I) = YCP00 * YIN(I2) + B00
          ZIN(I) = ZCP00 * ZIN(I2) + B00 * F00
          DXI(I) = DICP00 * XIN(I2) + XCP00 * DXI(I2)
          DYI(I) = DICP00 * YIN(I2) + YCP00 * DYI(I2)
          DZI(I) = DICP00 * ZIN(I2) + ZCP00 * DZI(I2)
          DXIXI(I) = DICP00 * DXI(I2) + DICP00 * DXI(I2)
          DYIYI(I) = DICP00 * DYI(I2) + DICP00 * DYI(I2)
          DZIZI(I) = DICP00 * DZI(I2) + DICP00 * DZI(I2)
        END IF
      END IF
C
C     ----- I(N,0)=(N-1)B10 * I(N-2,0) + C00 * I(N-1,0) -----
C
      IF (NMAX.GT.1) THEN
        DO 101 N=2,NMAX
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
          DXIXI(N0) = B10NL1 * DXIXI(NL20) + DIC00 * DXI(NL10)
     #    +             XC00 * DXIXI(NL10) + DIC00 * DXI(NL10)
          DYIYI(N0) = B10NL1 * DYIYI(NL20) + DIC00 * DYI(NL10)
     #    +             YC00 * DYIYI(NL10) + DIC00 * DYI(NL10)
          DZIZI(N0) = B10NL1 * DZIZI(NL20) + DIC00 * DZI(NL10)
     #    +             ZC00 * DZIZI(NL10) + DIC00 * DZI(NL10)
C
C     ----- I(N,1)=N * B00 * I(N-1,0) + CP00 * I(N,0) -----
C
          IF (MMAX.NE.0) THEN
ctph        N1=N0+32
            N1=N0+50
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
            DXIXI(N1) = B00N  * DXIXI(NL10) + DICP00 * DXI(N0)
     #      +           XCP00 * DXIXI(N0)   + DICP00 * DXI(N0)
            DYIYI(N1) = B00N  * DYIYI(NL10) + DICP00 * DYI(N0)
     #      +           YCP00 * DYIYI(N0)   + DICP00 * DYI(N0)
            DZIZI(N1) = B00N  * DZIZI(NL10) + DICP00 * DZI(N0)
     #      +           ZCP00 * DZIZI(N0)   + DICP00 * DZI(N0)
          END IF
  101   CONTINUE
      END IF
C
C     ----- I(0,M)=(M-1)BP01 * I(0,M-2) + CP00 * I(0,M-1) -----
C
      IF (MMAX.GT.1) THEN
        DO 102 M=2,MMAX
ctph      I0M = M*32 +  1 + IOFFST
ctph      I0ML1 = I0M - 32
ctph      I0ML2 = I0ML1 - 32
          I0M = M*50 +  1 + IOFFST
          I0ML1 = I0M - 50
          I0ML2 = I0ML1 - 50
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
          DXIXI(I0M) = BPML1 * DXIXI(I0ML2) + DICP00 * DXI(I0ML1)
     #    +            XCP00 * DXIXI(I0ML1) + DICP00 * DXI(I0ML1)
          DYIYI(I0M) = BPML1 * DYIYI(I0ML2) + DICP00 * DYI(I0ML1)
     #    +            YCP00 * DYIYI(I0ML1) + DICP00 * DYI(I0ML1)
          DZIZI(I0M) = BPML1 * DZIZI(I0ML2) + DICP00 * DZI(I0ML1)
     #    +            ZCP00 * DZIZI(I0ML1) + DICP00 * DZI(I0ML1)
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
            DXIXI(I1M) = B00M * DXIXI(I0ML1) + DIC00 * DXI(I0M)
     #      +            XC00 * DXIXI(I0M)   + DIC00 * DXI(I0M)
            DYIYI(I1M) = B00M * DYIYI(I0ML1) + DIC00 * DYI(I0M)
     #      +            YC00 * DYIYI(I0M)   + DIC00 * DYI(I0M)
            DZIZI(I1M) = B00M * DZIZI(I0ML1) + DIC00 * DZI(I0M)
     #      +            ZC00 * DZIZI(I0M)   + DIC00 * DZI(I0M)
          END IF
  102   CONTINUE
      END IF
C
C     ----- I(N,M)=(N-1)B10 * I(N-2,M) + M * B00 * I(N-1,M-1) +
C                                            C00 * I(N-1,M)   -----
C
      IF (NMAX.GT.1.AND.MMAX.GT.1) THEN
        DO 104 M=2,MMAX
ctph      M8=M*32 + IOFFST
          M8=M*50 + IOFFST
          DO 103 N=2,NMAX
            NM=N + 1 + M8
            NL1M=NM - 1
ctph        NL1ML1= NL1M - 32
            NL1ML1= NL1M - 50
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
            DXIXI(NM) = B10NL1 * DXIXI(NL2M) + DIC00 * DXI(NL1M)
     #      +             XC00 * DXIXI(NL1M) + DIC00 * DXI(NL1M)
     #      +             B00M * DXIXI(NL1ML1)
            DYIYI(NM) = B10NL1 * DYIYI(NL2M) + DIC00 * DYI(NL1M)
     #      +             YC00 * DYIYI(NL1M) + DIC00 * DYI(NL1M)
     #      +             B00M * DYIYI(NL1ML1)
            DZIZI(NM) = B10NL1 * DZIZI(NL2M) + DIC00 * DZI(NL1M)
     #      +             ZC00 * DZIZI(NL1M) + DIC00 * DZI(NL1M)
     #      +             B00M * DZIZI(NL1ML1)
  103     CONTINUE
  104   CONTINUE
      END IF
C
C     ----- I(I,J,M)=I(I+1,J-1,M) + (QI-QJ) * I(I,J-1,M) -----
C
      DO 107 M=0,MMAX
ctph    MPT = M * 32 + IOFFST
        MPT = M * 50 + IOFFST
        DO 106 J=1,NJMAX
ctph      JM = J * 8 + MPT
          JM = J *10 + MPT
          DO 105 I=1,NIMAX + NJMAX - J + 1
            IJM = I + JM
ctph        IJL1M = IJM - 8
            IJL1M = IJM -10
            I1JL1M = IJL1M + 1
            XIN(IJM) = XIN(I1JL1M) + XIMXJ * XIN(IJL1M)
            YIN(IJM) = YIN(I1JL1M) + YIMYJ * YIN(IJL1M)
            ZIN(IJM) = ZIN(I1JL1M) + ZIMZJ * ZIN(IJL1M)
            DXI(IJM) = DXI(I1JL1M) + XIMXJ * DXI(IJL1M) + XIN(IJL1M)
            DYI(IJM) = DYI(I1JL1M) + YIMYJ * DYI(IJL1M) + YIN(IJL1M)
            DZI(IJM) = DZI(I1JL1M) + ZIMZJ * DZI(IJL1M) + ZIN(IJL1M)
            DXIXI(IJM) =        DXIXI(I1JL1M) + DXI(IJL1M)
     #      +            XIMXJ * DXIXI(IJL1M) + DXI(IJL1M)
            DYIYI(IJM) =        DYIYI(I1JL1M) + DYI(IJL1M)
     #      +            YIMYJ * DYIYI(IJL1M) + DYI(IJL1M)
            DZIZI(IJM) =        DZIZI(I1JL1M) + DZI(IJL1M)
     #      +            ZIMZJ * DZIZI(IJL1M) + DZI(IJL1M)
  105     CONTINUE
  106   CONTINUE
  107 CONTINUE
C
C     ----- I(I,J,M,0) = I(I,J,M) -----
C
      DO 110 M=MMAX,0,-1
ctph    MNEW = M * 64 + IOFFST
ctph    MOLD = M * 32 + IOFFST
        MNEW = M *125 + IOFFST
        MOLD = M * 50 + IOFFST
        DO 109 J=NJMAX,0,-1
ctph      JMNEW = MNEW + J * 16
ctph      JMOLD = MOLD + J * 8
          JMNEW = MNEW + J * 25
          JMOLD = MOLD + J * 10
          DO 108 I=NIMAX+1,1,-1
            IJMNEW = JMNEW + I
            IJMOLD = JMOLD + I
            XIN(IJMNEW) = XIN(IJMOLD)
            YIN(IJMNEW) = YIN(IJMOLD)
            ZIN(IJMNEW) = ZIN(IJMOLD)
            DXI(IJMNEW) = DXI(IJMOLD)
            DYI(IJMNEW) = DYI(IJMOLD)
            DZI(IJMNEW) = DZI(IJMOLD)
            DXIXI(IJMNEW) = DXIXI(IJMOLD)
            DYIYI(IJMNEW) = DYIYI(IJMOLD)
            DZIZI(IJMNEW) = DZIZI(IJMOLD)
  108     CONTINUE
  109   CONTINUE
  110 CONTINUE
C
C     ----- I(I,J,K,L)=I(I,J,K+1,L-1) + (QK-QL) * I(I,J,K,L-1) -----
C
      DO 114 L=1,NLMAX
ctph    LPT = 4 * L + IOFFST
        LPT = 5 * L + IOFFST
        DO 113 K=0,NKMAX + NLMAX - L
ctph      KL = K * 64 + LPT
          KL = K *125 + LPT
          DO 112 J=0,NJMAX
ctph        JKL = J * 16 + KL
            JKL = J * 25 + KL
            DO 111 I=1,NIMAX+1
              IJKL = I + JKL
ctph          IJKLM1 = IJKL - 4
ctph          KP1LM1 = IJKLM1 + 64
              IJKLM1 = IJKL - 5
              KP1LM1 = IJKLM1 + 125
              XIN(IJKL) = XIN(KP1LM1) + XKMXL * XIN(IJKLM1)
              YIN(IJKL) = YIN(KP1LM1) + YKMYL * YIN(IJKLM1)
              ZIN(IJKL) = ZIN(KP1LM1) + ZKMZL * ZIN(IJKLM1)
              DXI(IJKL) = DXI(KP1LM1) + XKMXL * DXI(IJKLM1)
              DYI(IJKL) = DYI(KP1LM1) + YKMYL * DYI(IJKLM1)
              DZI(IJKL) = DZI(KP1LM1) + ZKMZL * DZI(IJKLM1)
              DXIXI(IJKL) =          DXIXI(KP1LM1)
     #        +              XKMXL * DXIXI(IJKLM1)
              DYIYI(IJKL) =          DYIYI(KP1LM1)
     #        +              YKMYL * DYIYI(IJKLM1)
              DZIZI(IJKL) =          DZIZI(KP1LM1)
     #        +              ZKMZL * DZIZI(IJKLM1)
  111       CONTINUE
  112     CONTINUE
  113   CONTINUE
  114 CONTINUE
C
C     ----- ADD DERIVATIVE OF EXP(-D) AND EXP(-G) -----
C
      DO 118 L=0,NLMAX
ctph    LPT = L*4 + IOFFST
        LPT = L*5 + IOFFST
        DO 117 K=0,NKMAX
ctph      KL = K*64 + LPT
          KL = K*125+ LPT
          DO 116 J=0,NJMAX
ctph        JKL = J*16 + KL
            JKL = J*25 + KL
            DO 115 I=1,NIMAX+1
              IJKL = I + JKL
              DXIXI(IJKL) = DXIXI(IJKL) + (XICNST*XICNST + DIDICN) *
     *        XIN(IJKL)+             XICNST * DXI(IJKL) + XICNST *
     *        DXI(IJKL)
              DYIYI(IJKL) = DYIYI(IJKL) + (YICNST*YICNST + DIDICN) *
     *        YIN(IJKL)+             YICNST * DYI(IJKL) + YICNST *
     *        DYI(IJKL)
              DZIZI(IJKL) = DZIZI(IJKL) + (ZICNST*ZICNST + DIDICN) *
     *        ZIN(IJKL)+             ZICNST * DZI(IJKL) + ZICNST *
     *        DZI(IJKL)
              DXI(IJKL) = DXI(IJKL) + XICNST * XIN(IJKL)
              DYI(IJKL) = DYI(IJKL) + YICNST * YIN(IJKL)
              DZI(IJKL) = DZI(IJKL) + ZICNST * ZIN(IJKL)
  115       CONTINUE
  116     CONTINUE
  117   CONTINUE
  118 CONTINUE
      RETURN
      END
