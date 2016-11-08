C
C
C----------------------------------------------------------------------
C
      SUBROUTINE FRMTAU(ORBSYM,FLOV,NIRRED,NO,NV,T2,PKT2,NDIMT2,T1,
     .                  FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      INTEGER FLOV(NIRRED,4),ORBSYM(NO+NV),FZO(NO),FZV(NV)
      REAL*8 T2(NDIMT2),PKT2(NDIMT2),T1(NO,NV)
C
C      FORM TAU
C
       NOFF=0
       DO 150 TSYM=1,NIRRED
          DO 140 USYM=1,NIRRED
             DO 130 VSYM=1,NIRRED
                UVSYM=IEOR(USYM-1,VSYM-1)+1
                IF(UVSYM.NE.TSYM) GO TO 130
                DO 125 U=FLOV(USYM,1),FLOV(USYM,2)
                   DO 120 V=FLOV(VSYM,1),FLOV(VSYM,2)
                      DO 119 BE=1,NV
                         BESYM=ORBSYM(BE+NO)
                         IF(FZV(BE).EQ.1) GO TO 119
                         DO 118 GA=1,BE
                            IF(FZV(GA).EQ.1) GO TO 118
                            GASYM=ORBSYM(GA+NO)
                            BEGSYM=IEOR(BESYM,GASYM)+1
                            IF(BEGSYM.NE.TSYM)GO TO 118
                            NOFF=NOFF+1
                            T2(NOFF)=PKT2(NOFF)+T1(U,BE)*T1(V,GA)
  118                    CONTINUE
  119                 CONTINUE
  120              CONTINUE
  125           CONTINUE
  130        CONTINUE
  140     CONTINUE
  150  CONTINUE
C
      RETURN
      END
