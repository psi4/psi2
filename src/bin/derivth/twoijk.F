C***********************************************************************
C* For array size changes from/to  (50625) or (10000) search C11-14-88 *
C* (50625) is for running in 15Mb while (10000) is for 7Mb.    RBR     *
C***********************************************************************
C* this small core version will only do up to f functions.  change all *
C* abdens(50625) and derint(50625,9) to (10000) and (10000,9)          *
C***********************************************************************
C* LAST UPDATED NOVEMBER 11, 1988 BY TRACY HAMILTON TO DO D5,F7,G9     *
C* DERIVATIVES (OR D6,F10,G15).  SEARCH FOR CTPH.                      *
C***********************************************************************
C* LAST UPDATED NOVEMBER 11, 1987 BY RICHARD REMINGTON                 *
C* VARIABLE DLIM CHANGED FROM 1.0E-10 TO 1.0E-16 FOR BETTER PRECISION; *
C* ALSO PUT DLIM IN COMMON TOLER                                       *
C***********************************************************************
C* LAST UPDATED ON APRIL 25, 1985 BY YUKIO YAMAGUCHI                   *
C***********************************************************************
      SUBROUTINE TWOIJK
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/EXPNT/EIIA,EJIA,EKIB,EIJIA,EKLIB,RHOTSQ,EI,EJ,EK,A,B
      COMMON/GEOIN/XIMXJ,YIMYJ,ZIMZJ,XKMXL,YKMYL,ZKMZL,
     1             XAMXB,YAMYB,ZAMZB
      COMMON/SETD /BP01,B00,B10,XCP00,XC00,YCP00,YC00,ZCP00,ZC00,F00,
     1             NIMAX,NJMAX,NKMAX,NLMAX,NMAX,MMAX,IOFFST
      COMMON/TOLER/DLIM
ctph  COMMON/TWOD0/XIN(2304),YIN(2304),ZIN(2304)
ctph  COMMON/TWODI/DXI(2304),DYI(2304),DZI(2304)
ctph  COMMON/TWODJ/DXJ(2304),DYJ(2304),DZJ(2304)
ctph  COMMON/TWODK/DXK(2304),DYK(2304),DZK(2304)
ctph  COMMON/TWDII/DXIXI(2304),DYIYI(2304),DZIZI(2304)
ctph  COMMON/TWDJJ/DXJXJ(2304),DYJYJ(2304),DZJZJ(2304)
ctph  COMMON/TWDKK/DXKXK(2304),DYKYK(2304),DZKZK(2304)
ctph  COMMON/TWDIJ/DXIXJ(2304),DYIYJ(2304),DZIZJ(2304)
ctph  COMMON/TWDIK/DXIXK(2304),DYIYK(2304),DZIZK(2304)
ctph  COMMON/TWDJK/DXJXK(2304),DYJYK(2304),DZJZK(2304)
      COMMON/TWOD0/XIN(6250),YIN(6250),ZIN(6250)
      COMMON/TWODI/DXI(6250),DYI(6250),DZI(6250)
      COMMON/TWODJ/DXJ(6250),DYJ(6250),DZJ(6250)
      COMMON/TWODK/DXK(6250),DYK(6250),DZK(6250)
      COMMON/TWDII/DXIXI(6250),DYIYI(6250),DZIZI(6250)
      COMMON/TWDJJ/DXJXJ(6250),DYJYJ(6250),DZJZJ(6250)
      COMMON/TWDKK/DXKXK(6250),DYKYK(6250),DZKZK(6250)
      COMMON/TWDIJ/DXIXJ(6250),DYIYJ(6250),DZIZJ(6250)
      COMMON/TWDIK/DXIXK(6250),DYIYK(6250),DZIZK(6250)
      COMMON/TWDJK/DXJXK(6250),DYJYK(6250),DZJZK(6250)
      COMMON/VARIL/RTIASQ,RTIBSQ,RTIAB,DIC00,DJC00,DKC00
      DATA ZERO,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /
C
      RTIASQ = RHOTSQ / (A*A)
      RTIBSQ = RHOTSQ / (B*B)
      RTIAB  = RHOTSQ / (A*B)
      DIC00  = EI/A - RTIASQ*EI - 1
      DJC00  = EJ/A - RTIASQ*EJ
      DKC00  = EK * RTIAB
      DICP00 = EI * RTIAB
      DJCP00 = EJ * RTIAB
      DKCP00 = EK/B - RTIBSQ*EK - 1
      RHO2T  = RHOTSQ + RHOTSQ
      RHO2TI = RHO2T*EIIA
      RHO2TJ = RHO2T*EJIA
      RHO2TK = RHO2T*EKIB
      XICNST =  -XIMXJ*EIJIA*TWO - XAMXB*RHO2TI
      YICNST =  -YIMYJ*EIJIA*TWO - YAMYB*RHO2TI
      ZICNST = (-ZIMZJ*EIJIA*TWO - ZAMZB*RHO2TI)
      XJCNST =   XIMXJ*EIJIA*TWO - XAMXB*RHO2TJ
      YJCNST =   YIMYJ*EIJIA*TWO - YAMYB*RHO2TJ
      ZJCNST = ( ZIMZJ*EIJIA*TWO - ZAMZB*RHO2TJ)
      XKCNST =  -XKMXL*EKLIB*TWO + XAMXB*RHO2TK
      YKCNST =  -YKMYL*EKLIB*TWO + YAMYB*RHO2TK
      ZKCNST = (-ZKMZL*EKLIB*TWO + ZAMZB*RHO2TK)
      DIDICN = -EIJIA*TWO - EIIA * RHO2TI
      DIDJCN =  EIJIA*TWO - EJIA * RHO2TI
      DIDKCN =              EKIB * RHO2TI
      DJDJCN = -EIJIA*TWO - EJIA * RHO2TJ
      DJDKCN =              EKIB * RHO2TJ
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
      DXJ(I) = ZERO
      DYJ(I) = ZERO
      DZJ(I) = ZERO
      DXK(I) = ZERO
      DYK(I) = ZERO
      DZK(I) = ZERO
      DXIXI(I) = ZERO
      DXIXJ(I) = ZERO
      DXIXK(I) = ZERO
      DXJXJ(I) = ZERO
      DXJXK(I) = ZERO
      DXKXK(I) = ZERO
      DYIYI(I) = ZERO
      DYIYJ(I) = ZERO
      DYIYK(I) = ZERO
      DYJYJ(I) = ZERO
      DYJYK(I) = ZERO
      DYKYK(I) = ZERO
      DZIZI(I) = ZERO
      DZIZJ(I) = ZERO
      DZIZK(I) = ZERO
      DZJZJ(I) = ZERO
      DZJZK(I) = ZERO
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
        DXJ(I2) = DJC00
        DYJ(I2) = DJC00
        DZJ(I2) = DJC00 * F00
        DXK(I2) = DKC00
        DYK(I2) = DKC00
        DZK(I2) = DKC00 * F00
        DXIXI(I2) = ZERO
        DXIXJ(I2) = ZERO
        DXIXK(I2) = ZERO
        DXJXJ(I2) = ZERO
        DXJXK(I2) = ZERO
        DXKXK(I2) = ZERO
        DYIYI(I2) = ZERO
        DYIYJ(I2) = ZERO
        DYIYK(I2) = ZERO
        DYJYJ(I2) = ZERO
        DYJYK(I2) = ZERO
        DYKYK(I2) = ZERO
        DZIZI(I2) = ZERO
        DZIZJ(I2) = ZERO
        DZIZK(I2) = ZERO
        DZJZJ(I2) = ZERO
        DZJZK(I2) = ZERO
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
        DXJ(I) = DJCP00
        DYJ(I) = DJCP00
        DZJ(I) = DJCP00 * F00
        DXK(I) = DKCP00
        DYK(I) = DKCP00
        DZK(I) = DKCP00 * F00
        DXIXI(I) = ZERO
        DXIXJ(I) = ZERO
        DXIXK(I) = ZERO
        DXJXJ(I) = ZERO
        DXJXK(I) = ZERO
        DXKXK(I) = ZERO
        DYIYI(I) = ZERO
        DYIYJ(I) = ZERO
        DYIYK(I) = ZERO
        DYJYJ(I) = ZERO
        DYJYK(I) = ZERO
        DYKYK(I) = ZERO
        DZIZI(I) = ZERO
        DZIZJ(I) = ZERO
        DZIZK(I) = ZERO
        DZJZJ(I) = ZERO
        DZJZK(I) = ZERO
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
          DXJ(I) = DJCP00 * XIN(I2) + XCP00 * DXJ(I2)
          DYJ(I) = DJCP00 * YIN(I2) + YCP00 * DYJ(I2)
          DZJ(I) = DJCP00 * ZIN(I2) + ZCP00 * DZJ(I2)
          DXK(I) = DKCP00 * XIN(I2) + XCP00 * DXK(I2)
          DYK(I) = DKCP00 * YIN(I2) + YCP00 * DYK(I2)
          DZK(I) = DKCP00 * ZIN(I2) + ZCP00 * DZK(I2)
          DXIXI(I) = DICP00 * DXI(I2) + DICP00 * DXI(I2)
          DXIXJ(I) = DICP00 * DXJ(I2) + DJCP00 * DXI(I2)
          DXIXK(I) = DICP00 * DXK(I2) + DKCP00 * DXI(I2)
          DXJXJ(I) = DJCP00 * DXJ(I2) + DJCP00 * DXJ(I2)
          DXJXK(I) = DJCP00 * DXK(I2) + DKCP00 * DXJ(I2)
          DXKXK(I) = DKCP00 * DXK(I2) + DKCP00 * DXK(I2)
          DYIYI(I) = DICP00 * DYI(I2) + DICP00 * DYI(I2)
          DYIYJ(I) = DICP00 * DYJ(I2) + DJCP00 * DYI(I2)
          DYIYK(I) = DICP00 * DYK(I2) + DKCP00 * DYI(I2)
          DYJYJ(I) = DJCP00 * DYJ(I2) + DJCP00 * DYJ(I2)
          DYJYK(I) = DJCP00 * DYK(I2) + DKCP00 * DYJ(I2)
          DYKYK(I) = DKCP00 * DYK(I2) + DKCP00 * DYK(I2)
          DZIZI(I) = DICP00 * DZI(I2) + DICP00 * DZI(I2)
          DZIZJ(I) = DICP00 * DZJ(I2) + DJCP00 * DZI(I2)
          DZIZK(I) = DICP00 * DZK(I2) + DKCP00 * DZI(I2)
          DZJZJ(I) = DJCP00 * DZJ(I2) + DJCP00 * DZJ(I2)
          DZJZK(I) = DJCP00 * DZK(I2) + DKCP00 * DZJ(I2)
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
          DXJ(N0) = B10NL1 * DXJ(NL20) + DJC00 * XIN(NL10) + XC00* DXJ
     *    (NL10)
          DYJ(N0) = B10NL1 * DYJ(NL20) + DJC00 * YIN(NL10) + YC00* DYJ
     *    (NL10)
          DZJ(N0) = B10NL1 * DZJ(NL20) + DJC00 * ZIN(NL10) + ZC00* DZJ
     *    (NL10)
          DXK(N0) = B10NL1 * DXK(NL20) + DKC00 * XIN(NL10) + XC00* DXK
     *    (NL10)
          DYK(N0) = B10NL1 * DYK(NL20) + DKC00 * YIN(NL10) + YC00* DYK
     *    (NL10)
          DZK(N0) = B10NL1 * DZK(NL20) + DKC00 * ZIN(NL10) + ZC00* DZK
     *    (NL10)
          DXIXI(N0) = B10NL1 * DXIXI(NL20) + DIC00 * DXI(NL10)
     #    +             XC00 * DXIXI(NL10) + DIC00 * DXI(NL10)
          DXIXJ(N0) = B10NL1 * DXIXJ(NL20) + DIC00 * DXJ(NL10)
     #    +             XC00 * DXIXJ(NL10) + DJC00 * DXI(NL10)
          DXIXK(N0) = B10NL1 * DXIXK(NL20) + DIC00 * DXK(NL10)
     #    +             XC00 * DXIXK(NL10) + DKC00 * DXI(NL10)
          DXJXJ(N0) = B10NL1 * DXJXJ(NL20) + DJC00 * DXJ(NL10)
     #    +             XC00 * DXJXJ(NL10) + DJC00 * DXJ(NL10)
          DXJXK(N0) = B10NL1 * DXJXK(NL20) + DJC00 * DXK(NL10)
     #    +             XC00 * DXJXK(NL10) + DKC00 * DXJ(NL10)
          DXKXK(N0) = B10NL1 * DXKXK(NL20) + DKC00 * DXK(NL10)
     #    +             XC00 * DXKXK(NL10) + DKC00 * DXK(NL10)
          DYIYI(N0) = B10NL1 * DYIYI(NL20) + DIC00 * DYI(NL10)
     #    +             YC00 * DYIYI(NL10) + DIC00 * DYI(NL10)
          DYIYJ(N0) = B10NL1 * DYIYJ(NL20) + DIC00 * DYJ(NL10)
     #    +             YC00 * DYIYJ(NL10) + DJC00 * DYI(NL10)
          DYIYK(N0) = B10NL1 * DYIYK(NL20) + DIC00 * DYK(NL10)
     #    +             YC00 * DYIYK(NL10) + DKC00 * DYI(NL10)
          DYJYJ(N0) = B10NL1 * DYJYJ(NL20) + DJC00 * DYJ(NL10)
     #    +             YC00 * DYJYJ(NL10) + DJC00 * DYJ(NL10)
          DYJYK(N0) = B10NL1 * DYJYK(NL20) + DJC00 * DYK(NL10)
     #    +             YC00 * DYJYK(NL10) + DKC00 * DYJ(NL10)
          DYKYK(N0) = B10NL1 * DYKYK(NL20) + DKC00 * DYK(NL10)
     #    +             YC00 * DYKYK(NL10) + DKC00 * DYK(NL10)
          DZIZI(N0) = B10NL1 * DZIZI(NL20) + DIC00 * DZI(NL10)
     #    +             ZC00 * DZIZI(NL10) + DIC00 * DZI(NL10)
          DZIZJ(N0) = B10NL1 * DZIZJ(NL20) + DIC00 * DZJ(NL10)
     #    +             ZC00 * DZIZJ(NL10) + DJC00 * DZI(NL10)
          DZIZK(N0) = B10NL1 * DZIZK(NL20) + DIC00 * DZK(NL10)
     #    +             ZC00 * DZIZK(NL10) + DKC00 * DZI(NL10)
          DZJZJ(N0) = B10NL1 * DZJZJ(NL20) + DJC00 * DZJ(NL10)
     #    +             ZC00 * DZJZJ(NL10) + DJC00 * DZJ(NL10)
          DZJZK(N0) = B10NL1 * DZJZK(NL20) + DJC00 * DZK(NL10)
     #    +             ZC00 * DZJZK(NL10) + DKC00 * DZJ(NL10)
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
            DXJ(N1) = B00N * DXJ(NL10) + DJCP00 * XIN(N0) + XCP00 *
     *      DXJ(N0)
            DYJ(N1) = B00N * DYJ(NL10) + DJCP00 * YIN(N0) + YCP00 *
     *      DYJ(N0)
            DZJ(N1) = B00N * DZJ(NL10) + DJCP00 * ZIN(N0) + ZCP00 *
     *      DZJ(N0)
            DXK(N1) = B00N * DXK(NL10) + DKCP00 * XIN(N0) + XCP00 *
     *      DXK(N0)
            DYK(N1) = B00N * DYK(NL10) + DKCP00 * YIN(N0) + YCP00 *
     *      DYK(N0)
            DZK(N1) = B00N * DZK(NL10) + DKCP00 * ZIN(N0) + ZCP00 *
     *      DZK(N0)
            DXIXI(N1) = B00N  * DXIXI(NL10) + DICP00 * DXI(N0)
     #      +           XCP00 * DXIXI(N0)   + DICP00 * DXI(N0)
            DXIXJ(N1) = B00N  * DXIXJ(NL10) + DICP00 * DXJ(N0)
     #      +           XCP00 * DXIXJ(N0)   + DJCP00 * DXI(N0)
            DXIXK(N1) = B00N  * DXIXK(NL10) + DICP00 * DXK(N0)
     #      +           XCP00 * DXIXK(N0)   + DKCP00 * DXI(N0)
            DXJXJ(N1) = B00N  * DXJXJ(NL10) + DJCP00 * DXJ(N0)
     #      +           XCP00 * DXJXJ(N0)   + DJCP00 * DXJ(N0)
            DXJXK(N1) = B00N  * DXJXK(NL10) + DJCP00 * DXK(N0)
     #      +           XCP00 * DXJXK(N0)   + DKCP00 * DXJ(N0)
            DXKXK(N1) = B00N  * DXKXK(NL10) + DKCP00 * DXK(N0)
     #      +           XCP00 * DXKXK(N0)   + DKCP00 * DXK(N0)
            DYIYI(N1) = B00N  * DYIYI(NL10) + DICP00 * DYI(N0)
     #      +           YCP00 * DYIYI(N0)   + DICP00 * DYI(N0)
            DYIYJ(N1) = B00N  * DYIYJ(NL10) + DICP00 * DYJ(N0)
     #      +           YCP00 * DYIYJ(N0)   + DJCP00 * DYI(N0)
            DYIYK(N1) = B00N  * DYIYK(NL10) + DICP00 * DYK(N0)
     #      +           YCP00 * DYIYK(N0)   + DKCP00 * DYI(N0)
            DYJYJ(N1) = B00N  * DYJYJ(NL10) + DJCP00 * DYJ(N0)
     #      +           YCP00 * DYJYJ(N0)   + DJCP00 * DYJ(N0)
            DYJYK(N1) = B00N  * DYJYK(NL10) + DJCP00 * DYK(N0)
     #      +           YCP00 * DYJYK(N0)   + DKCP00 * DYJ(N0)
            DYKYK(N1) = B00N  * DYKYK(NL10) + DKCP00 * DYK(N0)
     #      +           YCP00 * DYKYK(N0)   + DKCP00 * DYK(N0)
            DZIZI(N1) = B00N  * DZIZI(NL10) + DICP00 * DZI(N0)
     #      +           ZCP00 * DZIZI(N0)   + DICP00 * DZI(N0)
            DZIZJ(N1) = B00N  * DZIZJ(NL10) + DICP00 * DZJ(N0)
     #      +           ZCP00 * DZIZJ(N0)   + DJCP00 * DZI(N0)
            DZIZK(N1) = B00N  * DZIZK(NL10) + DICP00 * DZK(N0)
     #      +           ZCP00 * DZIZK(N0)   + DKCP00 * DZI(N0)
            DZJZJ(N1) = B00N  * DZJZJ(NL10) + DJCP00 * DZJ(N0)
     #      +           ZCP00 * DZJZJ(N0)   + DJCP00 * DZJ(N0)
            DZJZK(N1) = B00N  * DZJZK(NL10) + DJCP00 * DZK(N0)
     #      +           ZCP00 * DZJZK(N0)   + DKCP00 * DZJ(N0)
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
          DXJ(I0M) = BPML1 * DXJ(I0ML2) + DJCP00 * XIN(I0ML1)
     1    +                                XCP00 * DXJ(I0ML1)
          DYJ(I0M) = BPML1 * DYJ(I0ML2) + DJCP00 * YIN(I0ML1)
     1    +                                YCP00 * DYJ(I0ML1)
          DZJ(I0M) = BPML1 * DZJ(I0ML2) + DJCP00 * ZIN(I0ML1)
     1    +                                ZCP00 * DZJ(I0ML1)
          DXK(I0M) = BPML1 * DXK(I0ML2) + DKCP00 * XIN(I0ML1)
     1    +                                XCP00 * DXK(I0ML1)
          DYK(I0M) = BPML1 * DYK(I0ML2) + DKCP00 * YIN(I0ML1)
     1    +                                YCP00 * DYK(I0ML1)
          DZK(I0M) = BPML1 * DZK(I0ML2) + DKCP00 * ZIN(I0ML1)
     1    +                                ZCP00 * DZK(I0ML1)
          DXIXI(I0M) = BPML1 * DXIXI(I0ML2) + DICP00 * DXI(I0ML1)
     #    +            XCP00 * DXIXI(I0ML1) + DICP00 * DXI(I0ML1)
          DXIXJ(I0M) = BPML1 * DXIXJ(I0ML2) + DICP00 * DXJ(I0ML1)
     #    +            XCP00 * DXIXJ(I0ML1) + DJCP00 * DXI(I0ML1)
          DXIXK(I0M) = BPML1 * DXIXK(I0ML2) + DICP00 * DXK(I0ML1)
     #    +            XCP00 * DXIXK(I0ML1) + DKCP00 * DXI(I0ML1)
          DXJXJ(I0M) = BPML1 * DXJXJ(I0ML2) + DJCP00 * DXJ(I0ML1)
     #    +            XCP00 * DXJXJ(I0ML1) + DJCP00 * DXJ(I0ML1)
          DXJXK(I0M) = BPML1 * DXJXK(I0ML2) + DJCP00 * DXK(I0ML1)
     #    +            XCP00 * DXJXK(I0ML1) + DKCP00 * DXJ(I0ML1)
          DXKXK(I0M) = BPML1 * DXKXK(I0ML2) + DKCP00 * DXK(I0ML1)
     #    +            XCP00 * DXKXK(I0ML1) + DKCP00 * DXK(I0ML1)
          DYIYI(I0M) = BPML1 * DYIYI(I0ML2) + DICP00 * DYI(I0ML1)
     #    +            YCP00 * DYIYI(I0ML1) + DICP00 * DYI(I0ML1)
          DYIYJ(I0M) = BPML1 * DYIYJ(I0ML2) + DICP00 * DYJ(I0ML1)
     #    +            YCP00 * DYIYJ(I0ML1) + DJCP00 * DYI(I0ML1)
          DYIYK(I0M) = BPML1 * DYIYK(I0ML2) + DICP00 * DYK(I0ML1)
     #    +            YCP00 * DYIYK(I0ML1) + DKCP00 * DYI(I0ML1)
          DYJYJ(I0M) = BPML1 * DYJYJ(I0ML2) + DJCP00 * DYJ(I0ML1)
     #    +            YCP00 * DYJYJ(I0ML1) + DJCP00 * DYJ(I0ML1)
          DYJYK(I0M) = BPML1 * DYJYK(I0ML2) + DJCP00 * DYK(I0ML1)
     #    +            YCP00 * DYJYK(I0ML1) + DKCP00 * DYJ(I0ML1)
          DYKYK(I0M) = BPML1 * DYKYK(I0ML2) + DKCP00 * DYK(I0ML1)
     #    +            YCP00 * DYKYK(I0ML1) + DKCP00 * DYK(I0ML1)
          DZIZI(I0M) = BPML1 * DZIZI(I0ML2) + DICP00 * DZI(I0ML1)
     #    +            ZCP00 * DZIZI(I0ML1) + DICP00 * DZI(I0ML1)
          DZIZJ(I0M) = BPML1 * DZIZJ(I0ML2) + DICP00 * DZJ(I0ML1)
     #    +            ZCP00 * DZIZJ(I0ML1) + DJCP00 * DZI(I0ML1)
          DZIZK(I0M) = BPML1 * DZIZK(I0ML2) + DICP00 * DZK(I0ML1)
     #    +            ZCP00 * DZIZK(I0ML1) + DKCP00 * DZI(I0ML1)
          DZJZJ(I0M) = BPML1 * DZJZJ(I0ML2) + DJCP00 * DZJ(I0ML1)
     #    +            ZCP00 * DZJZJ(I0ML1) + DJCP00 * DZJ(I0ML1)
          DZJZK(I0M) = BPML1 * DZJZK(I0ML2) + DJCP00 * DZK(I0ML1)
     #    +            ZCP00 * DZJZK(I0ML1) + DKCP00 * DZJ(I0ML1)
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
            DXJ(I1M) = B00M * DXJ(I0ML1) + DJC00 * XIN(I0M) + XC00 *
     *      DXJ(I0M)
            DYJ(I1M) = B00M * DYJ(I0ML1) + DJC00 * YIN(I0M) + YC00 *
     *      DYJ(I0M)
            DZJ(I1M) = B00M * DZJ(I0ML1) + DJC00 * ZIN(I0M) + ZC00 *
     *      DZJ(I0M)
            DXK(I1M) = B00M * DXK(I0ML1) + DKC00 * XIN(I0M) + XC00 *
     *      DXK(I0M)
            DYK(I1M) = B00M * DYK(I0ML1) + DKC00 * YIN(I0M) + YC00 *
     *      DYK(I0M)
            DZK(I1M) = B00M * DZK(I0ML1) + DKC00 * ZIN(I0M) + ZC00 *
     *      DZK(I0M)
            DXIXI(I1M) = B00M * DXIXI(I0ML1) + DIC00 * DXI(I0M)
     #      +            XC00 * DXIXI(I0M)   + DIC00 * DXI(I0M)
            DXIXJ(I1M) = B00M * DXIXJ(I0ML1) + DIC00 * DXJ(I0M)
     #      +            XC00 * DXIXJ(I0M)   + DJC00 * DXI(I0M)
            DXIXK(I1M) = B00M * DXIXK(I0ML1) + DIC00 * DXK(I0M)
     #      +            XC00 * DXIXK(I0M)   + DKC00 * DXI(I0M)
            DXJXJ(I1M) = B00M * DXJXJ(I0ML1) + DJC00 * DXJ(I0M)
     #      +            XC00 * DXJXJ(I0M)   + DJC00 * DXJ(I0M)
            DXJXK(I1M) = B00M * DXJXK(I0ML1) + DJC00 * DXK(I0M)
     #      +            XC00 * DXJXK(I0M)   + DKC00 * DXJ(I0M)
            DXKXK(I1M) = B00M * DXKXK(I0ML1) + DKC00 * DXK(I0M)
     #      +            XC00 * DXKXK(I0M)   + DKC00 * DXK(I0M)
            DYIYI(I1M) = B00M * DYIYI(I0ML1) + DIC00 * DYI(I0M)
     #      +            YC00 * DYIYI(I0M)   + DIC00 * DYI(I0M)
            DYIYJ(I1M) = B00M * DYIYJ(I0ML1) + DIC00 * DYJ(I0M)
     #      +            YC00 * DYIYJ(I0M)   + DJC00 * DYI(I0M)
            DYIYK(I1M) = B00M * DYIYK(I0ML1) + DIC00 * DYK(I0M)
     #      +            YC00 * DYIYK(I0M)   + DKC00 * DYI(I0M)
            DYJYJ(I1M) = B00M * DYJYJ(I0ML1) + DJC00 * DYJ(I0M)
     #      +            YC00 * DYJYJ(I0M)   + DJC00 * DYJ(I0M)
            DYJYK(I1M) = B00M * DYJYK(I0ML1) + DJC00 * DYK(I0M)
     #      +            YC00 * DYJYK(I0M)   + DKC00 * DYJ(I0M)
            DYKYK(I1M) = B00M * DYKYK(I0ML1) + DKC00 * DYK(I0M)
     #      +            YC00 * DYKYK(I0M)   + DKC00 * DYK(I0M)
            DZIZI(I1M) = B00M * DZIZI(I0ML1) + DIC00 * DZI(I0M)
     #      +            ZC00 * DZIZI(I0M)   + DIC00 * DZI(I0M)
            DZIZJ(I1M) = B00M * DZIZJ(I0ML1) + DIC00 * DZJ(I0M)
     #      +            ZC00 * DZIZJ(I0M)   + DJC00 * DZI(I0M)
            DZIZK(I1M) = B00M * DZIZK(I0ML1) + DIC00 * DZK(I0M)
     #      +            ZC00 * DZIZK(I0M)   + DKC00 * DZI(I0M)
            DZJZJ(I1M) = B00M * DZJZJ(I0ML1) + DJC00 * DZJ(I0M)
     #      +            ZC00 * DZJZJ(I0M)   + DJC00 * DZJ(I0M)
            DZJZK(I1M) = B00M * DZJZK(I0ML1) + DJC00 * DZK(I0M)
     #      +            ZC00 * DZJZK(I0M)   + DKC00 * DZJ(I0M)
            DZKZK(I1M) = B00M * DZKZK(I0ML1) + DKC00 * DZK(I0M)
     #      +            ZC00 * DZKZK(I0M)   + DKC00 * DZK(I0M)
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
            DXJ(NM) = B10NL1 * DXJ(NL2M) + B00M * DXJ(NL1ML1)
     1      +          DJC00 * XIN(NL1M) + XC00  * DXJ(NL1M)
            DYJ(NM) = B10NL1 * DYJ(NL2M) + B00M * DYJ(NL1ML1)
     1      +          DJC00 * YIN(NL1M) + YC00  * DYJ(NL1M)
            DZJ(NM) = B10NL1 * DZJ(NL2M) + B00M * DZJ(NL1ML1)
     1      +          DJC00 * ZIN(NL1M) + ZC00  * DZJ(NL1M)
            DXK(NM) = B10NL1 * DXK(NL2M) + B00M * DXK(NL1ML1)
     1      +          DKC00 * XIN(NL1M) + XC00  * DXK(NL1M)
            DYK(NM) = B10NL1 * DYK(NL2M) + B00M * DYK(NL1ML1)
     1      +          DKC00 * YIN(NL1M) + YC00  * DYK(NL1M)
            DZK(NM) = B10NL1 * DZK(NL2M) + B00M * DZK(NL1ML1)
     1      +          DKC00 * ZIN(NL1M) + ZC00  * DZK(NL1M)
            DXIXI(NM) = B10NL1 * DXIXI(NL2M) + DIC00 * DXI(NL1M)
     #      +             XC00 * DXIXI(NL1M) + DIC00 * DXI(NL1M)
     #      +             B00M * DXIXI(NL1ML1)
            DXIXJ(NM) = B10NL1 * DXIXJ(NL2M) + DIC00 * DXJ(NL1M)
     #      +             XC00 * DXIXJ(NL1M) + DJC00 * DXI(NL1M)
     #      +             B00M * DXIXJ(NL1ML1)
            DXIXK(NM) = B10NL1 * DXIXK(NL2M) + DIC00 * DXK(NL1M)
     #      +             XC00 * DXIXK(NL1M) + DKC00 * DXI(NL1M)
     #      +             B00M * DXIXK(NL1ML1)
            DXJXJ(NM) = B10NL1 * DXJXJ(NL2M) + DJC00 * DXJ(NL1M)
     #      +             XC00 * DXJXJ(NL1M) + DJC00 * DXJ(NL1M)
     #      +             B00M * DXJXJ(NL1ML1)
            DXJXK(NM) = B10NL1 * DXJXK(NL2M) + DJC00 * DXK(NL1M)
     #      +             XC00 * DXJXK(NL1M) + DKC00 * DXJ(NL1M)
     #      +             B00M * DXJXK(NL1ML1)
            DXKXK(NM) = B10NL1 * DXKXK(NL2M) + DKC00 * DXK(NL1M)
     #      +             XC00 * DXKXK(NL1M) + DKC00 * DXK(NL1M)
     #      +             B00M * DXKXK(NL1ML1)
            DYIYI(NM) = B10NL1 * DYIYI(NL2M) + DIC00 * DYI(NL1M)
     #      +             YC00 * DYIYI(NL1M) + DIC00 * DYI(NL1M)
     #      +             B00M * DYIYI(NL1ML1)
            DYIYJ(NM) = B10NL1 * DYIYJ(NL2M) + DIC00 * DYJ(NL1M)
     #      +             YC00 * DYIYJ(NL1M) + DJC00 * DYI(NL1M)
     #      +             B00M * DYIYJ(NL1ML1)
            DYIYK(NM) = B10NL1 * DYIYK(NL2M) + DIC00 * DYK(NL1M)
     #      +             YC00 * DYIYK(NL1M) + DKC00 * DYI(NL1M)
     #      +             B00M * DYIYK(NL1ML1)
            DYJYJ(NM) = B10NL1 * DYJYJ(NL2M) + DJC00 * DYJ(NL1M)
     #      +             YC00 * DYJYJ(NL1M) + DJC00 * DYJ(NL1M)
     #      +             B00M * DYJYJ(NL1ML1)
            DYJYK(NM) = B10NL1 * DYJYK(NL2M) + DJC00 * DYK(NL1M)
     #      +             YC00 * DYJYK(NL1M) + DKC00 * DYJ(NL1M)
     #      +             B00M * DYJYK(NL1ML1)
            DYKYK(NM) = B10NL1 * DYKYK(NL2M) + DKC00 * DYK(NL1M)
     #      +             YC00 * DYKYK(NL1M) + DKC00 * DYK(NL1M)
     #      +             B00M * DYKYK(NL1ML1)
            DZIZI(NM) = B10NL1 * DZIZI(NL2M) + DIC00 * DZI(NL1M)
     #      +             ZC00 * DZIZI(NL1M) + DIC00 * DZI(NL1M)
     #      +             B00M * DZIZI(NL1ML1)
            DZIZJ(NM) = B10NL1 * DZIZJ(NL2M) + DIC00 * DZJ(NL1M)
     #      +             ZC00 * DZIZJ(NL1M) + DJC00 * DZI(NL1M)
     #      +             B00M * DZIZJ(NL1ML1)
            DZIZK(NM) = B10NL1 * DZIZK(NL2M) + DIC00 * DZK(NL1M)
     #      +             ZC00 * DZIZK(NL1M) + DKC00 * DZI(NL1M)
     #      +             B00M * DZIZK(NL1ML1)
            DZJZJ(NM) = B10NL1 * DZJZJ(NL2M) + DJC00 * DZJ(NL1M)
     #      +             ZC00 * DZJZJ(NL1M) + DJC00 * DZJ(NL1M)
     #      +             B00M * DZJZJ(NL1ML1)
            DZJZK(NM) = B10NL1 * DZJZK(NL2M) + DJC00 * DZK(NL1M)
     #      +             ZC00 * DZJZK(NL1M) + DKC00 * DZJ(NL1M)
     #      +             B00M * DZJZK(NL1ML1)
            DZKZK(NM) = B10NL1 * DZKZK(NL2M) + DKC00 * DZK(NL1M)
     #      +             ZC00 * DZKZK(NL1M) + DKC00 * DZK(NL1M)
     #      +             B00M * DZKZK(NL1ML1)
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
            DXJ(IJM) = DXJ(I1JL1M) + XIMXJ * DXJ(IJL1M) - XIN(IJL1M)
            DYJ(IJM) = DYJ(I1JL1M) + YIMYJ * DYJ(IJL1M) - YIN(IJL1M)
            DZJ(IJM) = DZJ(I1JL1M) + ZIMZJ * DZJ(IJL1M) - ZIN(IJL1M)
            DXK(IJM) = DXK(I1JL1M) + XIMXJ * DXK(IJL1M)
            DYK(IJM) = DYK(I1JL1M) + YIMYJ * DYK(IJL1M)
            DZK(IJM) = DZK(I1JL1M) + ZIMZJ * DZK(IJL1M)
            DXIXI(IJM) =        DXIXI(I1JL1M) + DXI(IJL1M)
     #      +            XIMXJ * DXIXI(IJL1M) + DXI(IJL1M)
            DXIXJ(IJM) =        DXIXJ(I1JL1M) + DXJ(IJL1M)
     #      +            XIMXJ * DXIXJ(IJL1M) - DXI(IJL1M)
            DXIXK(IJM) =        DXIXK(I1JL1M) + DXK(IJL1M)
     #      +            XIMXJ * DXIXK(IJL1M)
            DXJXJ(IJM) =        DXJXJ(I1JL1M) - DXJ(IJL1M)
     #      +            XIMXJ * DXJXJ(IJL1M) - DXJ(IJL1M)
            DXJXK(IJM) =        DXJXK(I1JL1M) - DXK(IJL1M)
     #      +            XIMXJ * DXJXK(IJL1M)
            DXKXK(IJM) =        DXKXK(I1JL1M)
     #      +            XIMXJ * DXKXK(IJL1M)
            DYIYI(IJM) =        DYIYI(I1JL1M) + DYI(IJL1M)
     #      +            YIMYJ * DYIYI(IJL1M) + DYI(IJL1M)
            DYIYJ(IJM) =        DYIYJ(I1JL1M) + DYJ(IJL1M)
     #      +            YIMYJ * DYIYJ(IJL1M) - DYI(IJL1M)
            DYIYK(IJM) =        DYIYK(I1JL1M) + DYK(IJL1M)
     #      +            YIMYJ * DYIYK(IJL1M)
            DYJYJ(IJM) =        DYJYJ(I1JL1M) - DYJ(IJL1M)
     #      +            YIMYJ * DYJYJ(IJL1M) - DYJ(IJL1M)
            DYJYK(IJM) =        DYJYK(I1JL1M) - DYK(IJL1M)
     #      +            YIMYJ * DYJYK(IJL1M)
            DYKYK(IJM) =        DYKYK(I1JL1M)
     #      +            YIMYJ * DYKYK(IJL1M)
            DZIZI(IJM) =        DZIZI(I1JL1M) + DZI(IJL1M)
     #      +            ZIMZJ * DZIZI(IJL1M) + DZI(IJL1M)
            DZIZJ(IJM) =        DZIZJ(I1JL1M) + DZJ(IJL1M)
     #      +            ZIMZJ * DZIZJ(IJL1M) - DZI(IJL1M)
            DZIZK(IJM) =        DZIZK(I1JL1M) + DZK(IJL1M)
     #      +            ZIMZJ * DZIZK(IJL1M)
            DZJZJ(IJM) =        DZJZJ(I1JL1M) - DZJ(IJL1M)
     #      +            ZIMZJ * DZJZJ(IJL1M) - DZJ(IJL1M)
            DZJZK(IJM) =        DZJZK(I1JL1M) - DZK(IJL1M)
     #      +            ZIMZJ * DZJZK(IJL1M)
            DZKZK(IJM) =        DZKZK(I1JL1M)
     #      +            ZIMZJ * DZKZK(IJL1M)
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
            DXJ(IJMNEW) = DXJ(IJMOLD)
            DYJ(IJMNEW) = DYJ(IJMOLD)
            DZJ(IJMNEW) = DZJ(IJMOLD)
            DXK(IJMNEW) = DXK(IJMOLD)
            DYK(IJMNEW) = DYK(IJMOLD)
            DZK(IJMNEW) = DZK(IJMOLD)
            DXIXI(IJMNEW) = DXIXI(IJMOLD)
            DXIXJ(IJMNEW) = DXIXJ(IJMOLD)
            DXIXK(IJMNEW) = DXIXK(IJMOLD)
            DXJXJ(IJMNEW) = DXJXJ(IJMOLD)
            DXJXK(IJMNEW) = DXJXK(IJMOLD)
            DXKXK(IJMNEW) = DXKXK(IJMOLD)
            DYIYI(IJMNEW) = DYIYI(IJMOLD)
            DYIYJ(IJMNEW) = DYIYJ(IJMOLD)
            DYIYK(IJMNEW) = DYIYK(IJMOLD)
            DYJYJ(IJMNEW) = DYJYJ(IJMOLD)
            DYJYK(IJMNEW) = DYJYK(IJMOLD)
            DYKYK(IJMNEW) = DYKYK(IJMOLD)
            DZIZI(IJMNEW) = DZIZI(IJMOLD)
            DZIZJ(IJMNEW) = DZIZJ(IJMOLD)
            DZIZK(IJMNEW) = DZIZK(IJMOLD)
            DZJZJ(IJMNEW) = DZJZJ(IJMOLD)
            DZJZK(IJMNEW) = DZJZK(IJMOLD)
            DZKZK(IJMNEW) = DZKZK(IJMOLD)
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
              DXJ(IJKL) = DXJ(KP1LM1) + XKMXL * DXJ(IJKLM1)
              DYJ(IJKL) = DYJ(KP1LM1) + YKMYL * DYJ(IJKLM1)
              DZJ(IJKL) = DZJ(KP1LM1) + ZKMZL * DZJ(IJKLM1)
              DXK(IJKL) = DXK(KP1LM1) + XKMXL * DXK(IJKLM1) + XIN(
     *        IJKLM1)
              DYK(IJKL) = DYK(KP1LM1) + YKMYL * DYK(IJKLM1) + YIN(
     *        IJKLM1)
              DZK(IJKL) = DZK(KP1LM1) + ZKMZL * DZK(IJKLM1) + ZIN(
     *        IJKLM1)
              DXIXI(IJKL) =          DXIXI(KP1LM1)
     #        +              XKMXL * DXIXI(IJKLM1)
              DXIXJ(IJKL) =          DXIXJ(KP1LM1)
     #        +              XKMXL * DXIXJ(IJKLM1)
              DXIXK(IJKL) =          DXIXK(KP1LM1)
     #        +              XKMXL * DXIXK(IJKLM1) + DXI(IJKLM1)
              DXJXJ(IJKL) =          DXJXJ(KP1LM1)
     #        +              XKMXL * DXJXJ(IJKLM1)
              DXJXK(IJKL) =          DXJXK(KP1LM1)
     #        +              XKMXL * DXJXK(IJKLM1) + DXJ(IJKLM1)
              DXKXK(IJKL) =          DXKXK(KP1LM1) + DXK(IJKLM1)
     #        +              XKMXL * DXKXK(IJKLM1) + DXK(IJKLM1)
              DYIYI(IJKL) =          DYIYI(KP1LM1)
     #        +              YKMYL * DYIYI(IJKLM1)
              DYIYJ(IJKL) =          DYIYJ(KP1LM1)
     #        +              YKMYL * DYIYJ(IJKLM1)
              DYIYK(IJKL) =          DYIYK(KP1LM1)
     #        +              YKMYL * DYIYK(IJKLM1) + DYI(IJKLM1)
              DYJYJ(IJKL) =          DYJYJ(KP1LM1)
     #        +              YKMYL * DYJYJ(IJKLM1)
              DYJYK(IJKL) =          DYJYK(KP1LM1)
     #        +              YKMYL * DYJYK(IJKLM1) + DYJ(IJKLM1)
              DYKYK(IJKL) =          DYKYK(KP1LM1) + DYK(IJKLM1)
     #        +              YKMYL * DYKYK(IJKLM1) + DYK(IJKLM1)
              DZIZI(IJKL) =          DZIZI(KP1LM1)
     #        +              ZKMZL * DZIZI(IJKLM1)
              DZIZJ(IJKL) =          DZIZJ(KP1LM1)
     #        +              ZKMZL * DZIZJ(IJKLM1)
              DZIZK(IJKL) =          DZIZK(KP1LM1)
     #        +              ZKMZL * DZIZK(IJKLM1) + DZI(IJKLM1)
              DZJZJ(IJKL) =          DZJZJ(KP1LM1)
     #        +              ZKMZL * DZJZJ(IJKLM1)
              DZJZK(IJKL) =          DZJZK(KP1LM1)
     #        +              ZKMZL * DZJZK(IJKLM1) + DZJ(IJKLM1)
              DZKZK(IJKL) =          DZKZK(KP1LM1) + DZK(IJKLM1)
     #        +              ZKMZL * DZKZK(IJKLM1) + DZK(IJKLM1)
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
              DXIXJ(IJKL) = DXIXJ(IJKL) + (XICNST*XJCNST + DIDJCN) *
     *        XIN(IJKL)+             XICNST * DXJ(IJKL) + XJCNST *
     *        DXI(IJKL)
              DXIXK(IJKL) = DXIXK(IJKL) + (XICNST*XKCNST + DIDKCN) *
     *        XIN(IJKL)+             XICNST * DXK(IJKL) + XKCNST *
     *        DXI(IJKL)
              DXJXJ(IJKL) = DXJXJ(IJKL) + (XJCNST*XJCNST + DJDJCN) *
     *        XIN(IJKL)+             XJCNST * DXJ(IJKL) + XJCNST *
     *        DXJ(IJKL)
              DXJXK(IJKL) = DXJXK(IJKL) + (XJCNST*XKCNST + DJDKCN) *
     *        XIN(IJKL)+             XJCNST * DXK(IJKL) + XKCNST *
     *        DXJ(IJKL)
              DXKXK(IJKL) = DXKXK(IJKL) + (XKCNST*XKCNST + DKDKCN) *
     *        XIN(IJKL)+             XKCNST * DXK(IJKL) + XKCNST *
     *        DXK(IJKL)
              DYIYI(IJKL) = DYIYI(IJKL) + (YICNST*YICNST + DIDICN) *
     *        YIN(IJKL)+             YICNST * DYI(IJKL) + YICNST *
     *        DYI(IJKL)
              DYIYJ(IJKL) = DYIYJ(IJKL) + (YICNST*YJCNST + DIDJCN) *
     *        YIN(IJKL)+             YICNST * DYJ(IJKL) + YJCNST *
     *        DYI(IJKL)
              DYIYK(IJKL) = DYIYK(IJKL) + (YICNST*YKCNST + DIDKCN) *
     *        YIN(IJKL)+             YICNST * DYK(IJKL) + YKCNST *
     *        DYI(IJKL)
              DYJYJ(IJKL) = DYJYJ(IJKL) + (YJCNST*YJCNST + DJDJCN) *
     *        YIN(IJKL)+             YJCNST * DYJ(IJKL) + YJCNST *
     *        DYJ(IJKL)
              DYJYK(IJKL) = DYJYK(IJKL) + (YJCNST*YKCNST + DJDKCN) *
     *        YIN(IJKL)+             YJCNST * DYK(IJKL) + YKCNST *
     *        DYJ(IJKL)
              DYKYK(IJKL) = DYKYK(IJKL) + (YKCNST*YKCNST + DKDKCN) *
     *        YIN(IJKL)+             YKCNST * DYK(IJKL) + YKCNST *
     *        DYK(IJKL)
              DZIZI(IJKL) = DZIZI(IJKL) + (ZICNST*ZICNST + DIDICN) *
     *        ZIN(IJKL)+             ZICNST * DZI(IJKL) + ZICNST *
     *        DZI(IJKL)
              DZIZJ(IJKL) = DZIZJ(IJKL) + (ZICNST*ZJCNST + DIDJCN) *
     *        ZIN(IJKL)+             ZICNST * DZJ(IJKL) + ZJCNST *
     *        DZI(IJKL)
              DZIZK(IJKL) = DZIZK(IJKL) + (ZICNST*ZKCNST + DIDKCN) *
     *        ZIN(IJKL)+             ZICNST * DZK(IJKL) + ZKCNST *
     *        DZI(IJKL)
              DZJZJ(IJKL) = DZJZJ(IJKL) + (ZJCNST*ZJCNST + DJDJCN) *
     *        ZIN(IJKL)+             ZJCNST * DZJ(IJKL) + ZJCNST *
     *        DZJ(IJKL)
              DZJZK(IJKL) = DZJZK(IJKL) + (ZJCNST*ZKCNST + DJDKCN) *
     *        ZIN(IJKL)+             ZJCNST * DZK(IJKL) + ZKCNST *
     *        DZJ(IJKL)
              DZKZK(IJKL) = DZKZK(IJKL) + (ZKCNST*ZKCNST + DKDKCN) *
     *        ZIN(IJKL)+             ZKCNST * DZK(IJKL) + ZKCNST *
     *        DZK(IJKL)
              DXI(IJKL) = DXI(IJKL) + XICNST * XIN(IJKL)
              DYI(IJKL) = DYI(IJKL) + YICNST * YIN(IJKL)
              DZI(IJKL) = DZI(IJKL) + ZICNST * ZIN(IJKL)
              DXJ(IJKL) = DXJ(IJKL) + XJCNST * XIN(IJKL)
              DYJ(IJKL) = DYJ(IJKL) + YJCNST * YIN(IJKL)
              DZJ(IJKL) = DZJ(IJKL) + ZJCNST * ZIN(IJKL)
              DXK(IJKL) = DXK(IJKL) + XKCNST * XIN(IJKL)
              DYK(IJKL) = DYK(IJKL) + YKCNST * YIN(IJKL)
              DZK(IJKL) = DZK(IJKL) + ZKCNST * ZIN(IJKL)
  115       CONTINUE
  116     CONTINUE
  117   CONTINUE
  118 CONTINUE
      RETURN
      END
