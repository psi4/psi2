      SUBROUTINE TWOIK
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/EXPNT/EIIA,EJIA,EKIB,EIJIA,EKLIB,RHOTSQ,EI,EJ,EK,A,B
      COMMON/GEOIN/XIMXJ,YIMYJ,ZIMZJ,XKMXL,YKMYL,ZKMZL,
     1             XAMXB,YAMYB,ZAMZB
      COMMON/IKTWO/DIC00,DKC00,DICP00,DKCP00,
     1             XICNST,YICNST,ZICNST,XKCNST,YKCNST,ZKCNST,
     2             DIDICN,DIDKCN,DKDKCN
      COMMON/SETD /BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,
     1             NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,IOFFST
      COMMON/TOLER/DLIM
ctph  COMMON/TWOD0/XIN(2304),YIN(2304),ZIN(2304)
ctph  COMMON/TWODI/DXI(2304),DYI(2304),DZI(2304)
ctph  COMMON/TWODK/DXK(2304),DYK(2304),DZK(2304)
ctph  COMMON/TWDII/DXIXI(2304),DYIYI(2304),DZIZI(2304)
ctph  COMMON/TWDKK/DXKXK(2304),DYKYK(2304),DZKZK(2304)
ctph  COMMON/TWDIK/DXIXK(2304),DYIYK(2304),DZIZK(2304)
      COMMON/TWOD0/XIN(6250),YIN(6250),ZIN(6250)
      COMMON/TWODI/DXI(6250),DYI(6250),DZI(6250)
      COMMON/TWODK/DXK(6250),DYK(6250),DZK(6250)
      COMMON/TWDII/DXIXI(6250),DYIYI(6250),DZIZI(6250)
      COMMON/TWDKK/DXKXK(6250),DYKYK(6250),DZKZK(6250)
      COMMON/TWDIK/DXIXK(6250),DYIYK(6250),DZIZK(6250)
      DATA ZERO,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /
C
      RTIASQ = RHOTSQ / (A*A)
      RTIBSQ = RHOTSQ / (B*B)
      RTIAB  = RHOTSQ / (A*B)
      DIC00  = EI/A - RTIASQ*EI - 1
      DKC00  = EK * RTIAB
      DICP00 = EI * RTIAB
      DKCP00 = EK/B - RTIBSQ*EK - 1
      RHO2T  = RHOTSQ + RHOTSQ
      RHO2TI = RHO2T*EIIA
      RHO2TK = RHO2T*EKIB
      XICNST =  -XIMXJ*EIJIA*TWO - XAMXB*RHO2TI
      YICNST =  -YIMYJ*EIJIA*TWO - YAMYB*RHO2TI
      ZICNST = (-ZIMZJ*EIJIA*TWO - ZAMZB*RHO2TI)
      XKCNST =  -XKMXL*EKLIB*TWO + XAMXB*RHO2TK
      YKCNST =  -YKMYL*EKLIB*TWO + YAMYB*RHO2TK
      ZKCNST = (-ZKMZL*EKLIB*TWO + ZAMZB*RHO2TK)
      DIDICN = -EIJIA*TWO - EIIA * RHO2TI
      DIDKCN =              EKIB * RHO2TI
      DKDKCN = -EKLIB*TWO - EKIB * RHO2TK
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
      DXK(I) = ZERO
      DYK(I) = ZERO
      DZK(I) = ZERO
      DXIXI(I) = ZERO
      DXIXK(I) = ZERO
      DXKXK(I) = ZERO
      DYIYI(I) = ZERO
      DYIYK(I) = ZERO
      DYKYK(I) = ZERO
      DZIZI(I) = ZERO
      DZIZK(I) = ZERO
      DZKZK(I) = ZERO
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
        DXK(I2) = DKC00
        DYK(I2) = DKC00
        DZK(I2) = DKC00 * F00
        DXIXI(I2) = ZERO
        DXIXK(I2) = ZERO
        DXKXK(I2) = ZERO
        DYIYI(I2) = ZERO
        DYIYK(I2) = ZERO
        DYKYK(I2) = ZERO
        DZIZI(I2) = ZERO
        DZIZK(I2) = ZERO
        DZKZK(I2) = ZERO
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
        DXK(I) = DKCP00
        DYK(I) = DKCP00
        DZK(I) = DKCP00 * F00
        DXIXI(I) = ZERO
        DXIXK(I) = ZERO
        DXKXK(I) = ZERO
        DYIYI(I) = ZERO
        DYIYK(I) = ZERO
        DYKYK(I) = ZERO
        DZIZI(I) = ZERO
        DZIZK(I) = ZERO
        DZKZK(I) = ZERO
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
          DXK(I) = DKCP00 * XIN(I2) + XCP00 * DXK(I2)
          DYK(I) = DKCP00 * YIN(I2) + YCP00 * DYK(I2)
          DZK(I) = DKCP00 * ZIN(I2) + ZCP00 * DZK(I2)
          DXIXI(I) = DICP00 * DXI(I2) + DICP00 * DXI(I2)
          DXIXK(I) = DICP00 * DXK(I2) + DKCP00 * DXI(I2)
          DXKXK(I) = DKCP00 * DXK(I2) + DKCP00 * DXK(I2)
          DYIYI(I) = DICP00 * DYI(I2) + DICP00 * DYI(I2)
          DYIYK(I) = DICP00 * DYK(I2) + DKCP00 * DYI(I2)
          DYKYK(I) = DKCP00 * DYK(I2) + DKCP00 * DYK(I2)
          DZIZI(I) = DICP00 * DZI(I2) + DICP00 * DZI(I2)
          DZIZK(I) = DICP00 * DZK(I2) + DKCP00 * DZI(I2)
          DZKZK(I) = DKCP00 * DZK(I2) + DKCP00 * DZK(I2)
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
          DXK(N0) = B10NL1 * DXK(NL20) + DKC00 * XIN(NL10) + XC00* DXK
     *    (NL10)
          DYK(N0) = B10NL1 * DYK(NL20) + DKC00 * YIN(NL10) + YC00* DYK
     *    (NL10)
          DZK(N0) = B10NL1 * DZK(NL20) + DKC00 * ZIN(NL10) + ZC00* DZK
     *    (NL10)
          DXIXI(N0) = B10NL1 * DXIXI(NL20) + DIC00 * DXI(NL10)
     #    +             XC00 * DXIXI(NL10) + DIC00 * DXI(NL10)
          DXIXK(N0) = B10NL1 * DXIXK(NL20) + DIC00 * DXK(NL10)
     #    +             XC00 * DXIXK(NL10) + DKC00 * DXI(NL10)
          DXKXK(N0) = B10NL1 * DXKXK(NL20) + DKC00 * DXK(NL10)
     #    +             XC00 * DXKXK(NL10) + DKC00 * DXK(NL10)
          DYIYI(N0) = B10NL1 * DYIYI(NL20) + DIC00 * DYI(NL10)
     #    +             YC00 * DYIYI(NL10) + DIC00 * DYI(NL10)
          DYIYK(N0) = B10NL1 * DYIYK(NL20) + DIC00 * DYK(NL10)
     #    +             YC00 * DYIYK(NL10) + DKC00 * DYI(NL10)
          DYKYK(N0) = B10NL1 * DYKYK(NL20) + DKC00 * DYK(NL10)
     #    +             YC00 * DYKYK(NL10) + DKC00 * DYK(NL10)
          DZIZI(N0) = B10NL1 * DZIZI(NL20) + DIC00 * DZI(NL10)
     #    +             ZC00 * DZIZI(NL10) + DIC00 * DZI(NL10)
          DZIZK(N0) = B10NL1 * DZIZK(NL20) + DIC00 * DZK(NL10)
     #    +             ZC00 * DZIZK(NL10) + DKC00 * DZI(NL10)
          DZKZK(N0) = B10NL1 * DZKZK(NL20) + DKC00 * DZK(NL10)
     #    +             ZC00 * DZKZK(NL10) + DKC00 * DZK(NL10)
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
            DXK(N1) = B00N * DXK(NL10) + DKCP00 * XIN(N0) + XCP00 *
     *      DXK(N0)
            DYK(N1) = B00N * DYK(NL10) + DKCP00 * YIN(N0) + YCP00 *
     *      DYK(N0)
            DZK(N1) = B00N * DZK(NL10) + DKCP00 * ZIN(N0) + ZCP00 *
     *      DZK(N0)
            DXIXI(N1) = B00N  * DXIXI(NL10) + DICP00 * DXI(N0)
     #      +           XCP00 * DXIXI(N0)   + DICP00 * DXI(N0)
            DXIXK(N1) = B00N  * DXIXK(NL10) + DICP00 * DXK(N0)
     #      +           XCP00 * DXIXK(N0)   + DKCP00 * DXI(N0)
            DXKXK(N1) = B00N  * DXKXK(NL10) + DKCP00 * DXK(N0)
     #      +           XCP00 * DXKXK(N0)   + DKCP00 * DXK(N0)
            DYIYI(N1) = B00N  * DYIYI(NL10) + DICP00 * DYI(N0)
     #      +           YCP00 * DYIYI(N0)   + DICP00 * DYI(N0)
            DYIYK(N1) = B00N  * DYIYK(NL10) + DICP00 * DYK(N0)
     #      +           YCP00 * DYIYK(N0)   + DKCP00 * DYI(N0)
            DYKYK(N1) = B00N  * DYKYK(NL10) + DKCP00 * DYK(N0)
     #      +           YCP00 * DYKYK(N0)   + DKCP00 * DYK(N0)
            DZIZI(N1) = B00N  * DZIZI(NL10) + DICP00 * DZI(N0)
     #      +           ZCP00 * DZIZI(N0)   + DICP00 * DZI(N0)
            DZIZK(N1) = B00N  * DZIZK(NL10) + DICP00 * DZK(N0)
     #      +           ZCP00 * DZIZK(N0)   + DKCP00 * DZI(N0)
            DZKZK(N1) = B00N  * DZKZK(NL10) + DKCP00 * DZK(N0)
     #      +           ZCP00 * DZKZK(N0)   + DKCP00 * DZK(N0)
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
          DXK(I0M) = BPML1 * DXK(I0ML2) + DKCP00 * XIN(I0ML1)
     1    +                                XCP00 * DXK(I0ML1)
          DYK(I0M) = BPML1 * DYK(I0ML2) + DKCP00 * YIN(I0ML1)
     1    +                                YCP00 * DYK(I0ML1)
          DZK(I0M) = BPML1 * DZK(I0ML2) + DKCP00 * ZIN(I0ML1)
     1    +                                ZCP00 * DZK(I0ML1)
          DXIXI(I0M) = BPML1 * DXIXI(I0ML2) + DICP00 * DXI(I0ML1)
     #    +            XCP00 * DXIXI(I0ML1) + DICP00 * DXI(I0ML1)
          DXIXK(I0M) = BPML1 * DXIXK(I0ML2) + DICP00 * DXK(I0ML1)
     #    +            XCP00 * DXIXK(I0ML1) + DKCP00 * DXI(I0ML1)
          DXKXK(I0M) = BPML1 * DXKXK(I0ML2) + DKCP00 * DXK(I0ML1)
     #    +            XCP00 * DXKXK(I0ML1) + DKCP00 * DXK(I0ML1)
          DYIYI(I0M) = BPML1 * DYIYI(I0ML2) + DICP00 * DYI(I0ML1)
     #    +            YCP00 * DYIYI(I0ML1) + DICP00 * DYI(I0ML1)
          DYIYK(I0M) = BPML1 * DYIYK(I0ML2) + DICP00 * DYK(I0ML1)
     #    +            YCP00 * DYIYK(I0ML1) + DKCP00 * DYI(I0ML1)
          DYKYK(I0M) = BPML1 * DYKYK(I0ML2) + DKCP00 * DYK(I0ML1)
     #    +            YCP00 * DYKYK(I0ML1) + DKCP00 * DYK(I0ML1)
          DZIZI(I0M) = BPML1 * DZIZI(I0ML2) + DICP00 * DZI(I0ML1)
     #    +            ZCP00 * DZIZI(I0ML1) + DICP00 * DZI(I0ML1)
          DZIZK(I0M) = BPML1 * DZIZK(I0ML2) + DICP00 * DZK(I0ML1)
     #    +            ZCP00 * DZIZK(I0ML1) + DKCP00 * DZI(I0ML1)
          DZKZK(I0M) = BPML1 * DZKZK(I0ML2) + DKCP00 * DZK(I0ML1)
     #    +            ZCP00 * DZKZK(I0ML1) + DKCP00 * DZK(I0ML1)
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
            DXK(I1M) = B00M * DXK(I0ML1) + DKC00 * XIN(I0M) + XC00 *
     *      DXK(I0M)
            DYK(I1M) = B00M * DYK(I0ML1) + DKC00 * YIN(I0M) + YC00 *
     *      DYK(I0M)
            DZK(I1M) = B00M * DZK(I0ML1) + DKC00 * ZIN(I0M) + ZC00 *
     *      DZK(I0M)
            DXIXI(I1M) = B00M * DXIXI(I0ML1) + DIC00 * DXI(I0M)
     #      +            XC00 * DXIXI(I0M)   + DIC00 * DXI(I0M)
            DXIXK(I1M) = B00M * DXIXK(I0ML1) + DIC00 * DXK(I0M)
     #      +            XC00 * DXIXK(I0M)   + DKC00 * DXI(I0M)
            DXKXK(I1M) = B00M * DXKXK(I0ML1) + DKC00 * DXK(I0M)
     #      +            XC00 * DXKXK(I0M)   + DKC00 * DXK(I0M)
            DYIYI(I1M) = B00M * DYIYI(I0ML1) + DIC00 * DYI(I0M)
     #      +            YC00 * DYIYI(I0M)   + DIC00 * DYI(I0M)
            DYIYK(I1M) = B00M * DYIYK(I0ML1) + DIC00 * DYK(I0M)
     #      +            YC00 * DYIYK(I0M)   + DKC00 * DYI(I0M)
            DYKYK(I1M) = B00M * DYKYK(I0ML1) + DKC00 * DYK(I0M)
     #      +            YC00 * DYKYK(I0M)   + DKC00 * DYK(I0M)
            DZIZI(I1M) = B00M * DZIZI(I0ML1) + DIC00 * DZI(I0M)
     #      +            ZC00 * DZIZI(I0M)   + DIC00 * DZI(I0M)
            DZIZK(I1M) = B00M * DZIZK(I0ML1) + DIC00 * DZK(I0M)
     #      +            ZC00 * DZIZK(I0M)   + DKC00 * DZI(I0M)
            DZKZK(I1M) = B00M * DZKZK(I0ML1) + DKC00 * DZK(I0M)
     #      +            ZC00 * DZKZK(I0M)   + DKC00 * DZK(I0M)
          END IF
  102   CONTINUE
      END IF
C
      CALL TWOIK1
C
      RETURN
      END
