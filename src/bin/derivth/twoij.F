      SUBROUTINE TWOIJ
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/EXPNT/EIIA,EJIA,EKIB,EIJIA,EKLIB,RHOTSQ,EI,EJ,EK,A,B
      COMMON/GEOIN/XIMXJ,YIMYJ,ZIMZJ,XKMXL,YKMYL,ZKMZL,
     1             XAMXB,YAMYB,ZAMZB
      COMMON/IJTWO/DIC00,DJC00,DICP00,DJCP00,
     1             XICNST,YICNST,ZICNST,XJCNST,YJCNST,ZJCNST,
     2             DIDICN,DIDJCN,DJDJCN
      COMMON/SETD /BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,
     1             NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,IOFFST
      COMMON/TOLER/DLIM
ctph  COMMON/TWOD0/XIN(2304),YIN(2304),ZIN(2304)
ctph  COMMON/TWODI/DXI(2304),DYI(2304),DZI(2304)
ctph  COMMON/TWODJ/DXJ(2304),DYJ(2304),DZJ(2304)
ctph  COMMON/TWDII/DXIXI(2304),DYIYI(2304),DZIZI(2304)
ctph  COMMON/TWDJJ/DXJXJ(2304),DYJYJ(2304),DZJZJ(2304)
ctph  COMMON/TWDIJ/DXIXJ(2304),DYIYJ(2304),DZIZJ(2304)
      COMMON/TWOD0/XIN(6250),YIN(6250),ZIN(6250)
      COMMON/TWODI/DXI(6250),DYI(6250),DZI(6250)
      COMMON/TWODJ/DXJ(6250),DYJ(6250),DZJ(6250)
      COMMON/TWDII/DXIXI(6250),DYIYI(6250),DZIZI(6250)
      COMMON/TWDJJ/DXJXJ(6250),DYJYJ(6250),DZJZJ(6250)
      COMMON/TWDIJ/DXIXJ(6250),DYIYJ(6250),DZIZJ(6250)
      DATA ZERO,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /
C
      RTIASQ = RHOTSQ / (A*A)
      RTIBSQ = RHOTSQ / (B*B)
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
      DIDICN = -EIJIA*TWO - EIIA * RHO2TI
      DIDJCN =  EIJIA*TWO - EJIA * RHO2TI
      DJDJCN = -EIJIA*TWO - EJIA * RHO2TJ
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
      DXIXI(I) = ZERO
      DXIXJ(I) = ZERO
      DXJXJ(I) = ZERO
      DYIYI(I) = ZERO
      DYIYJ(I) = ZERO
      DYJYJ(I) = ZERO
      DZIZI(I) = ZERO
      DZIZJ(I) = ZERO
      DZJZJ(I) = ZERO
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
        DXIXI(I2) = ZERO
        DXIXJ(I2) = ZERO
        DXJXJ(I2) = ZERO
        DYIYI(I2) = ZERO
        DYIYJ(I2) = ZERO
        DYJYJ(I2) = ZERO
        DZIZI(I2) = ZERO
        DZIZJ(I2) = ZERO
        DZJZJ(I2) = ZERO
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
        DXJ(I) = DJCP00
        DYJ(I) = DJCP00
        DZJ(I) = DJCP00 * F00
        DXIXI(I) = ZERO
        DXIXJ(I) = ZERO
        DXJXJ(I) = ZERO
        DYIYI(I) = ZERO
        DYIYJ(I) = ZERO
        DYJYJ(I) = ZERO
        DZIZI(I) = ZERO
        DZIZJ(I) = ZERO
        DZJZJ(I) = ZERO
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
          DXJ(I) = DJCP00 * XIN(I2) + XCP00 * DXJ(I2)
          DYJ(I) = DJCP00 * YIN(I2) + YCP00 * DYJ(I2)
          DZJ(I) = DJCP00 * ZIN(I2) + ZCP00 * DZJ(I2)
          DXIXI(I) = DICP00 * DXI(I2) + DICP00 * DXI(I2)
          DXIXJ(I) = DICP00 * DXJ(I2) + DJCP00 * DXI(I2)
          DXJXJ(I) = DJCP00 * DXJ(I2) + DJCP00 * DXJ(I2)
          DYIYI(I) = DICP00 * DYI(I2) + DICP00 * DYI(I2)
          DYIYJ(I) = DICP00 * DYJ(I2) + DJCP00 * DYI(I2)
          DYJYJ(I) = DJCP00 * DYJ(I2) + DJCP00 * DYJ(I2)
          DZIZI(I) = DICP00 * DZI(I2) + DICP00 * DZI(I2)
          DZIZJ(I) = DICP00 * DZJ(I2) + DJCP00 * DZI(I2)
          DZJZJ(I) = DJCP00 * DZJ(I2) + DJCP00 * DZJ(I2)
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
          DXJ(N0) = B10NL1 * DXJ(NL20) + DJC00 * XIN(NL10) + XC00* DXJ
     *    (NL10)
          DYJ(N0) = B10NL1 * DYJ(NL20) + DJC00 * YIN(NL10) + YC00* DYJ
     *    (NL10)
          DZJ(N0) = B10NL1 * DZJ(NL20) + DJC00 * ZIN(NL10) + ZC00* DZJ
     *    (NL10)
          DXIXI(N0) = B10NL1 * DXIXI(NL20) + DIC00 * DXI(NL10)
     #    +             XC00 * DXIXI(NL10) + DIC00 * DXI(NL10)
          DXIXJ(N0) = B10NL1 * DXIXJ(NL20) + DIC00 * DXJ(NL10)
     #    +             XC00 * DXIXJ(NL10) + DJC00 * DXI(NL10)
          DXJXJ(N0) = B10NL1 * DXJXJ(NL20) + DJC00 * DXJ(NL10)
     #    +             XC00 * DXJXJ(NL10) + DJC00 * DXJ(NL10)
          DYIYI(N0) = B10NL1 * DYIYI(NL20) + DIC00 * DYI(NL10)
     #    +             YC00 * DYIYI(NL10) + DIC00 * DYI(NL10)
          DYIYJ(N0) = B10NL1 * DYIYJ(NL20) + DIC00 * DYJ(NL10)
     #    +             YC00 * DYIYJ(NL10) + DJC00 * DYI(NL10)
          DYJYJ(N0) = B10NL1 * DYJYJ(NL20) + DJC00 * DYJ(NL10)
     #    +             YC00 * DYJYJ(NL10) + DJC00 * DYJ(NL10)
          DZIZI(N0) = B10NL1 * DZIZI(NL20) + DIC00 * DZI(NL10)
     #    +             ZC00 * DZIZI(NL10) + DIC00 * DZI(NL10)
          DZIZJ(N0) = B10NL1 * DZIZJ(NL20) + DIC00 * DZJ(NL10)
     #    +             ZC00 * DZIZJ(NL10) + DJC00 * DZI(NL10)
          DZJZJ(N0) = B10NL1 * DZJZJ(NL20) + DJC00 * DZJ(NL10)
     #    +             ZC00 * DZJZJ(NL10) + DJC00 * DZJ(NL10)
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
            DXJ(N1) = B00N * DXJ(NL10) + DJCP00 * XIN(N0) + XCP00 *
     *      DXJ(N0)
            DYJ(N1) = B00N * DYJ(NL10) + DJCP00 * YIN(N0) + YCP00 *
     *      DYJ(N0)
            DZJ(N1) = B00N * DZJ(NL10) + DJCP00 * ZIN(N0) + ZCP00 *
     *      DZJ(N0)
            DXIXI(N1) = B00N  * DXIXI(NL10) + DICP00 * DXI(N0)
     #      +           XCP00 * DXIXI(N0)   + DICP00 * DXI(N0)
            DXIXJ(N1) = B00N  * DXIXJ(NL10) + DICP00 * DXJ(N0)
     #      +           XCP00 * DXIXJ(N0)   + DJCP00 * DXI(N0)
            DXJXJ(N1) = B00N  * DXJXJ(NL10) + DJCP00 * DXJ(N0)
     #      +           XCP00 * DXJXJ(N0)   + DJCP00 * DXJ(N0)
            DYIYI(N1) = B00N  * DYIYI(NL10) + DICP00 * DYI(N0)
     #      +           YCP00 * DYIYI(N0)   + DICP00 * DYI(N0)
            DYIYJ(N1) = B00N  * DYIYJ(NL10) + DICP00 * DYJ(N0)
     #      +           YCP00 * DYIYJ(N0)   + DJCP00 * DYI(N0)
            DYJYJ(N1) = B00N  * DYJYJ(NL10) + DJCP00 * DYJ(N0)
     #      +           YCP00 * DYJYJ(N0)   + DJCP00 * DYJ(N0)
            DZIZI(N1) = B00N  * DZIZI(NL10) + DICP00 * DZI(N0)
     #      +           ZCP00 * DZIZI(N0)   + DICP00 * DZI(N0)
            DZIZJ(N1) = B00N  * DZIZJ(NL10) + DICP00 * DZJ(N0)
     #      +           ZCP00 * DZIZJ(N0)   + DJCP00 * DZI(N0)
            DZJZJ(N1) = B00N  * DZJZJ(NL10) + DJCP00 * DZJ(N0)
     #      +           ZCP00 * DZJZJ(N0)   + DJCP00 * DZJ(N0)
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
          DXJ(I0M) = BPML1 * DXJ(I0ML2) + DJCP00 * XIN(I0ML1)
     1    +                                XCP00 * DXJ(I0ML1)
          DYJ(I0M) = BPML1 * DYJ(I0ML2) + DJCP00 * YIN(I0ML1)
     1    +                                YCP00 * DYJ(I0ML1)
          DZJ(I0M) = BPML1 * DZJ(I0ML2) + DJCP00 * ZIN(I0ML1)
     1    +                                ZCP00 * DZJ(I0ML1)
          DXIXI(I0M) = BPML1 * DXIXI(I0ML2) + DICP00 * DXI(I0ML1)
     #    +            XCP00 * DXIXI(I0ML1) + DICP00 * DXI(I0ML1)
          DXIXJ(I0M) = BPML1 * DXIXJ(I0ML2) + DICP00 * DXJ(I0ML1)
     #    +            XCP00 * DXIXJ(I0ML1) + DJCP00 * DXI(I0ML1)
          DXJXJ(I0M) = BPML1 * DXJXJ(I0ML2) + DJCP00 * DXJ(I0ML1)
     #    +            XCP00 * DXJXJ(I0ML1) + DJCP00 * DXJ(I0ML1)
          DYIYI(I0M) = BPML1 * DYIYI(I0ML2) + DICP00 * DYI(I0ML1)
     #    +            YCP00 * DYIYI(I0ML1) + DICP00 * DYI(I0ML1)
          DYIYJ(I0M) = BPML1 * DYIYJ(I0ML2) + DICP00 * DYJ(I0ML1)
     #    +            YCP00 * DYIYJ(I0ML1) + DJCP00 * DYI(I0ML1)
          DYJYJ(I0M) = BPML1 * DYJYJ(I0ML2) + DJCP00 * DYJ(I0ML1)
     #    +            YCP00 * DYJYJ(I0ML1) + DJCP00 * DYJ(I0ML1)
          DZIZI(I0M) = BPML1 * DZIZI(I0ML2) + DICP00 * DZI(I0ML1)
     #    +            ZCP00 * DZIZI(I0ML1) + DICP00 * DZI(I0ML1)
          DZIZJ(I0M) = BPML1 * DZIZJ(I0ML2) + DICP00 * DZJ(I0ML1)
     #    +            ZCP00 * DZIZJ(I0ML1) + DJCP00 * DZI(I0ML1)
          DZJZJ(I0M) = BPML1 * DZJZJ(I0ML2) + DJCP00 * DZJ(I0ML1)
     #    +            ZCP00 * DZJZJ(I0ML1) + DJCP00 * DZJ(I0ML1)
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
            DXIXI(I1M) = B00M * DXIXI(I0ML1) + DIC00 * DXI(I0M)
     #      +            XC00 * DXIXI(I0M)   + DIC00 * DXI(I0M)
            DXIXJ(I1M) = B00M * DXIXJ(I0ML1) + DIC00 * DXJ(I0M)
     #      +            XC00 * DXIXJ(I0M)   + DJC00 * DXI(I0M)
            DXJXJ(I1M) = B00M * DXJXJ(I0ML1) + DJC00 * DXJ(I0M)
     #      +            XC00 * DXJXJ(I0M)   + DJC00 * DXJ(I0M)
            DYIYI(I1M) = B00M * DYIYI(I0ML1) + DIC00 * DYI(I0M)
     #      +            YC00 * DYIYI(I0M)   + DIC00 * DYI(I0M)
            DYIYJ(I1M) = B00M * DYIYJ(I0ML1) + DIC00 * DYJ(I0M)
     #      +            YC00 * DYIYJ(I0M)   + DJC00 * DYI(I0M)
            DYJYJ(I1M) = B00M * DYJYJ(I0ML1) + DJC00 * DYJ(I0M)
     #      +            YC00 * DYJYJ(I0M)   + DJC00 * DYJ(I0M)
            DZIZI(I1M) = B00M * DZIZI(I0ML1) + DIC00 * DZI(I0M)
     #      +            ZC00 * DZIZI(I0M)   + DIC00 * DZI(I0M)
            DZIZJ(I1M) = B00M * DZIZJ(I0ML1) + DIC00 * DZJ(I0M)
     #      +            ZC00 * DZIZJ(I0M)   + DJC00 * DZI(I0M)
            DZJZJ(I1M) = B00M * DZJZJ(I0ML1) + DJC00 * DZJ(I0M)
     #      +            ZC00 * DZJZJ(I0M)   + DJC00 * DZJ(I0M)
          END IF
  102   CONTINUE
      END IF
C
      CALL TWOIJ1
C
      RETURN
      END
