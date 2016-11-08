      SUBROUTINE TWOIK1
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
            DXK(NM) = B10NL1 * DXK(NL2M) + B00M * DXK(NL1ML1)
     1      +          DKC00 * XIN(NL1M) + XC00  * DXK(NL1M)
            DYK(NM) = B10NL1 * DYK(NL2M) + B00M * DYK(NL1ML1)
     1      +          DKC00 * YIN(NL1M) + YC00  * DYK(NL1M)
            DZK(NM) = B10NL1 * DZK(NL2M) + B00M * DZK(NL1ML1)
     1      +          DKC00 * ZIN(NL1M) + ZC00  * DZK(NL1M)
            DXIXI(NM) = B10NL1 * DXIXI(NL2M) + DIC00 * DXI(NL1M)
     #      +             XC00 * DXIXI(NL1M) + DIC00 * DXI(NL1M)
     #      +             B00M * DXIXI(NL1ML1)
            DXIXK(NM) = B10NL1 * DXIXK(NL2M) + DIC00 * DXK(NL1M)
     #      +             XC00 * DXIXK(NL1M) + DKC00 * DXI(NL1M)
     #      +             B00M * DXIXK(NL1ML1)
            DXKXK(NM) = B10NL1 * DXKXK(NL2M) + DKC00 * DXK(NL1M)
     #      +             XC00 * DXKXK(NL1M) + DKC00 * DXK(NL1M)
     #      +             B00M * DXKXK(NL1ML1)
            DYIYI(NM) = B10NL1 * DYIYI(NL2M) + DIC00 * DYI(NL1M)
     #      +             YC00 * DYIYI(NL1M) + DIC00 * DYI(NL1M)
     #      +             B00M * DYIYI(NL1ML1)
            DYIYK(NM) = B10NL1 * DYIYK(NL2M) + DIC00 * DYK(NL1M)
     #      +             YC00 * DYIYK(NL1M) + DKC00 * DYI(NL1M)
     #      +             B00M * DYIYK(NL1ML1)
            DYKYK(NM) = B10NL1 * DYKYK(NL2M) + DKC00 * DYK(NL1M)
     #      +             YC00 * DYKYK(NL1M) + DKC00 * DYK(NL1M)
     #      +             B00M * DYKYK(NL1ML1)
            DZIZI(NM) = B10NL1 * DZIZI(NL2M) + DIC00 * DZI(NL1M)
     #      +             ZC00 * DZIZI(NL1M) + DIC00 * DZI(NL1M)
     #      +             B00M * DZIZI(NL1ML1)
            DZIZK(NM) = B10NL1 * DZIZK(NL2M) + DIC00 * DZK(NL1M)
     #      +             ZC00 * DZIZK(NL1M) + DKC00 * DZI(NL1M)
     #      +             B00M * DZIZK(NL1ML1)
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
            DXK(IJM) = DXK(I1JL1M) + XIMXJ * DXK(IJL1M)
            DYK(IJM) = DYK(I1JL1M) + YIMYJ * DYK(IJL1M)
            DZK(IJM) = DZK(I1JL1M) + ZIMZJ * DZK(IJL1M)
            DXIXI(IJM) =        DXIXI(I1JL1M) + DXI(IJL1M)
     #      +            XIMXJ * DXIXI(IJL1M) + DXI(IJL1M)
            DXIXK(IJM) =        DXIXK(I1JL1M) + DXK(IJL1M)
     #      +            XIMXJ * DXIXK(IJL1M)
            DXKXK(IJM) =        DXKXK(I1JL1M)
     #      +            XIMXJ * DXKXK(IJL1M)
            DYIYI(IJM) =        DYIYI(I1JL1M) + DYI(IJL1M)
     #      +            YIMYJ * DYIYI(IJL1M) + DYI(IJL1M)
            DYIYK(IJM) =        DYIYK(I1JL1M) + DYK(IJL1M)
     #      +            YIMYJ * DYIYK(IJL1M)
            DYKYK(IJM) =        DYKYK(I1JL1M)
     #      +            YIMYJ * DYKYK(IJL1M)
            DZIZI(IJM) =        DZIZI(I1JL1M) + DZI(IJL1M)
     #      +            ZIMZJ * DZIZI(IJL1M) + DZI(IJL1M)
            DZIZK(IJM) =        DZIZK(I1JL1M) + DZK(IJL1M)
     #      +            ZIMZJ * DZIZK(IJL1M)
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
            DXK(IJMNEW) = DXK(IJMOLD)
            DYK(IJMNEW) = DYK(IJMOLD)
            DZK(IJMNEW) = DZK(IJMOLD)
            DXIXI(IJMNEW) = DXIXI(IJMOLD)
            DXIXK(IJMNEW) = DXIXK(IJMOLD)
            DXKXK(IJMNEW) = DXKXK(IJMOLD)
            DYIYI(IJMNEW) = DYIYI(IJMOLD)
            DYIYK(IJMNEW) = DYIYK(IJMOLD)
            DYKYK(IJMNEW) = DYKYK(IJMOLD)
            DZIZI(IJMNEW) = DZIZI(IJMOLD)
            DZIZK(IJMNEW) = DZIZK(IJMOLD)
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
              DXK(IJKL) = DXK(KP1LM1) + XKMXL * DXK(IJKLM1) + XIN(
     *        IJKLM1)
              DYK(IJKL) = DYK(KP1LM1) + YKMYL * DYK(IJKLM1) + YIN(
     *        IJKLM1)
              DZK(IJKL) = DZK(KP1LM1) + ZKMZL * DZK(IJKLM1) + ZIN(
     *        IJKLM1)
              DXIXI(IJKL) =          DXIXI(KP1LM1)
     #        +              XKMXL * DXIXI(IJKLM1)
              DXIXK(IJKL) =          DXIXK(KP1LM1)
     #        +              XKMXL * DXIXK(IJKLM1) + DXI(IJKLM1)
              DXKXK(IJKL) =          DXKXK(KP1LM1) + DXK(IJKLM1)
     #        +              XKMXL * DXKXK(IJKLM1) + DXK(IJKLM1)
              DYIYI(IJKL) =          DYIYI(KP1LM1)
     #        +              YKMYL * DYIYI(IJKLM1)
              DYIYK(IJKL) =          DYIYK(KP1LM1)
     #        +              YKMYL * DYIYK(IJKLM1) + DYI(IJKLM1)
              DYKYK(IJKL) =          DYKYK(KP1LM1) + DYK(IJKLM1)
     #        +              YKMYL * DYKYK(IJKLM1) + DYK(IJKLM1)
              DZIZI(IJKL) =          DZIZI(KP1LM1)
     #        +              ZKMZL * DZIZI(IJKLM1)
              DZIZK(IJKL) =          DZIZK(KP1LM1)
     #        +              ZKMZL * DZIZK(IJKLM1) + DZI(IJKLM1)
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
              DXIXK(IJKL) = DXIXK(IJKL) + (XICNST*XKCNST + DIDKCN) *
     *        XIN(IJKL)+             XICNST * DXK(IJKL) + XKCNST *
     *        DXI(IJKL)
              DXKXK(IJKL) = DXKXK(IJKL) + (XKCNST*XKCNST + DKDKCN) *
     *        XIN(IJKL)+             XKCNST * DXK(IJKL) + XKCNST *
     *        DXK(IJKL)
              DYIYI(IJKL) = DYIYI(IJKL) + (YICNST*YICNST + DIDICN) *
     *        YIN(IJKL)+             YICNST * DYI(IJKL) + YICNST *
     *        DYI(IJKL)
              DYIYK(IJKL) = DYIYK(IJKL) + (YICNST*YKCNST + DIDKCN) *
     *        YIN(IJKL)+             YICNST * DYK(IJKL) + YKCNST *
     *        DYI(IJKL)
              DYKYK(IJKL) = DYKYK(IJKL) + (YKCNST*YKCNST + DKDKCN) *
     *        YIN(IJKL)+             YKCNST * DYK(IJKL) + YKCNST *
     *        DYK(IJKL)
              DZIZI(IJKL) = DZIZI(IJKL) + (ZICNST*ZICNST + DIDICN) *
     *        ZIN(IJKL)+             ZICNST * DZI(IJKL) + ZICNST *
     *        DZI(IJKL)
              DZIZK(IJKL) = DZIZK(IJKL) + (ZICNST*ZKCNST + DIDKCN) *
     *        ZIN(IJKL)+             ZICNST * DZK(IJKL) + ZKCNST *
     *        DZI(IJKL)
              DZKZK(IJKL) = DZKZK(IJKL) + (ZKCNST*ZKCNST + DKDKCN) *
     *        ZIN(IJKL)+             ZKCNST * DZK(IJKL) + ZKCNST *
     *        DZK(IJKL)
              DXI(IJKL) = DXI(IJKL) + XICNST * XIN(IJKL)
              DYI(IJKL) = DYI(IJKL) + YICNST * YIN(IJKL)
              DZI(IJKL) = DZI(IJKL) + ZICNST * ZIN(IJKL)
              DXK(IJKL) = DXK(IJKL) + XKCNST * XIN(IJKL)
              DYK(IJKL) = DYK(IJKL) + YKCNST * YIN(IJKL)
              DZK(IJKL) = DZK(IJKL) + ZKCNST * ZIN(IJKL)
  115       CONTINUE
  116     CONTINUE
  117   CONTINUE
  118 CONTINUE
      RETURN
      END
