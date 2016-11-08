      SUBROUTINE TWOIJ1
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
            DXIXI(NM) = B10NL1 * DXIXI(NL2M) + DIC00 * DXI(NL1M)
     #      +             XC00 * DXIXI(NL1M) + DIC00 * DXI(NL1M)
     #      +             B00M * DXIXI(NL1ML1)
            DXIXJ(NM) = B10NL1 * DXIXJ(NL2M) + DIC00 * DXJ(NL1M)
     #      +             XC00 * DXIXJ(NL1M) + DJC00 * DXI(NL1M)
     #      +             B00M * DXIXJ(NL1ML1)
            DXJXJ(NM) = B10NL1 * DXJXJ(NL2M) + DJC00 * DXJ(NL1M)
     #      +             XC00 * DXJXJ(NL1M) + DJC00 * DXJ(NL1M)
     #      +             B00M * DXJXJ(NL1ML1)
            DYIYI(NM) = B10NL1 * DYIYI(NL2M) + DIC00 * DYI(NL1M)
     #      +             YC00 * DYIYI(NL1M) + DIC00 * DYI(NL1M)
     #      +             B00M * DYIYI(NL1ML1)
            DYIYJ(NM) = B10NL1 * DYIYJ(NL2M) + DIC00 * DYJ(NL1M)
     #      +             YC00 * DYIYJ(NL1M) + DJC00 * DYI(NL1M)
     #      +             B00M * DYIYJ(NL1ML1)
            DYJYJ(NM) = B10NL1 * DYJYJ(NL2M) + DJC00 * DYJ(NL1M)
     #      +             YC00 * DYJYJ(NL1M) + DJC00 * DYJ(NL1M)
     #      +             B00M * DYJYJ(NL1ML1)
            DZIZI(NM) = B10NL1 * DZIZI(NL2M) + DIC00 * DZI(NL1M)
     #      +             ZC00 * DZIZI(NL1M) + DIC00 * DZI(NL1M)
     #      +             B00M * DZIZI(NL1ML1)
            DZIZJ(NM) = B10NL1 * DZIZJ(NL2M) + DIC00 * DZJ(NL1M)
     #      +             ZC00 * DZIZJ(NL1M) + DJC00 * DZI(NL1M)
     #      +             B00M * DZIZJ(NL1ML1)
            DZJZJ(NM) = B10NL1 * DZJZJ(NL2M) + DJC00 * DZJ(NL1M)
     #      +             ZC00 * DZJZJ(NL1M) + DJC00 * DZJ(NL1M)
     #      +             B00M * DZJZJ(NL1ML1)
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
            DXIXI(IJM) =        DXIXI(I1JL1M) + DXI(IJL1M)
     #      +            XIMXJ * DXIXI(IJL1M) + DXI(IJL1M)
            DXIXJ(IJM) =        DXIXJ(I1JL1M) + DXJ(IJL1M)
     #      +            XIMXJ * DXIXJ(IJL1M) - DXI(IJL1M)
            DXJXJ(IJM) =        DXJXJ(I1JL1M) - DXJ(IJL1M)
     #      +            XIMXJ * DXJXJ(IJL1M) - DXJ(IJL1M)
            DYIYI(IJM) =        DYIYI(I1JL1M) + DYI(IJL1M)
     #      +            YIMYJ * DYIYI(IJL1M) + DYI(IJL1M)
            DYIYJ(IJM) =        DYIYJ(I1JL1M) + DYJ(IJL1M)
     #      +            YIMYJ * DYIYJ(IJL1M) - DYI(IJL1M)
            DYJYJ(IJM) =        DYJYJ(I1JL1M) - DYJ(IJL1M)
     #      +            YIMYJ * DYJYJ(IJL1M) - DYJ(IJL1M)
            DZIZI(IJM) =        DZIZI(I1JL1M) + DZI(IJL1M)
     #      +            ZIMZJ * DZIZI(IJL1M) + DZI(IJL1M)
            DZIZJ(IJM) =        DZIZJ(I1JL1M) + DZJ(IJL1M)
     #      +            ZIMZJ * DZIZJ(IJL1M) - DZI(IJL1M)
            DZJZJ(IJM) =        DZJZJ(I1JL1M) - DZJ(IJL1M)
     #      +            ZIMZJ * DZJZJ(IJL1M) - DZJ(IJL1M)
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
            DXIXI(IJMNEW) = DXIXI(IJMOLD)
            DXIXJ(IJMNEW) = DXIXJ(IJMOLD)
            DXJXJ(IJMNEW) = DXJXJ(IJMOLD)
            DYIYI(IJMNEW) = DYIYI(IJMOLD)
            DYIYJ(IJMNEW) = DYIYJ(IJMOLD)
            DYJYJ(IJMNEW) = DYJYJ(IJMOLD)
            DZIZI(IJMNEW) = DZIZI(IJMOLD)
            DZIZJ(IJMNEW) = DZIZJ(IJMOLD)
            DZJZJ(IJMNEW) = DZJZJ(IJMOLD)
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
              DXIXI(IJKL) =          DXIXI(KP1LM1)
     #        +              XKMXL * DXIXI(IJKLM1)
              DXIXJ(IJKL) =          DXIXJ(KP1LM1)
     #        +              XKMXL * DXIXJ(IJKLM1)
              DXJXJ(IJKL) =          DXJXJ(KP1LM1)
     #        +              XKMXL * DXJXJ(IJKLM1)
              DYIYI(IJKL) =          DYIYI(KP1LM1)
     #        +              YKMYL * DYIYI(IJKLM1)
              DYIYJ(IJKL) =          DYIYJ(KP1LM1)
     #        +              YKMYL * DYIYJ(IJKLM1)
              DYJYJ(IJKL) =          DYJYJ(KP1LM1)
     #        +              YKMYL * DYJYJ(IJKLM1)
              DZIZI(IJKL) =          DZIZI(KP1LM1)
     #        +              ZKMZL * DZIZI(IJKLM1)
              DZIZJ(IJKL) =          DZIZJ(KP1LM1)
     #        +              ZKMZL * DZIZJ(IJKLM1)
              DZJZJ(IJKL) =          DZJZJ(KP1LM1)
     #        +              ZKMZL * DZJZJ(IJKLM1)
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
              DXJXJ(IJKL) = DXJXJ(IJKL) + (XJCNST*XJCNST + DJDJCN) *
     *        XIN(IJKL)+             XJCNST * DXJ(IJKL) + XJCNST *
     *        DXJ(IJKL)
              DYIYI(IJKL) = DYIYI(IJKL) + (YICNST*YICNST + DIDICN) *
     *        YIN(IJKL)+             YICNST * DYI(IJKL) + YICNST *
     *        DYI(IJKL)
              DYIYJ(IJKL) = DYIYJ(IJKL) + (YICNST*YJCNST + DIDJCN) *
     *        YIN(IJKL)+             YICNST * DYJ(IJKL) + YJCNST *
     *        DYI(IJKL)
              DYJYJ(IJKL) = DYJYJ(IJKL) + (YJCNST*YJCNST + DJDJCN) *
     *        YIN(IJKL)+             YJCNST * DYJ(IJKL) + YJCNST *
     *        DYJ(IJKL)
              DZIZI(IJKL) = DZIZI(IJKL) + (ZICNST*ZICNST + DIDICN) *
     *        ZIN(IJKL)+             ZICNST * DZI(IJKL) + ZICNST *
     *        DZI(IJKL)
              DZIZJ(IJKL) = DZIZJ(IJKL) + (ZICNST*ZJCNST + DIDJCN) *
     *        ZIN(IJKL)+             ZICNST * DZJ(IJKL) + ZJCNST *
     *        DZI(IJKL)
              DZJZJ(IJKL) = DZJZJ(IJKL) + (ZJCNST*ZJCNST + DJDJCN) *
     *        ZIN(IJKL)+             ZJCNST * DZJ(IJKL) + ZJCNST *
     *        DZJ(IJKL)
              DXI(IJKL) = DXI(IJKL) + XICNST * XIN(IJKL)
              DYI(IJKL) = DYI(IJKL) + YICNST * YIN(IJKL)
              DZI(IJKL) = DZI(IJKL) + ZICNST * ZIN(IJKL)
              DXJ(IJKL) = DXJ(IJKL) + XJCNST * XIN(IJKL)
              DYJ(IJKL) = DYJ(IJKL) + YJCNST * YIN(IJKL)
              DZJ(IJKL) = DZJ(IJKL) + ZJCNST * ZIN(IJKL)
  115       CONTINUE
  116     CONTINUE
  117   CONTINUE
  118 CONTINUE
      RETURN
      END
