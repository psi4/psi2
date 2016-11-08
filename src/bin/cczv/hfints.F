C
C-----------------------------------------------------------------------
C
      SUBROUTINE HFINTS(HMAT,MOINTS,NO2,NV2,NTRO,NTRV,NO,NV,
     .                  LENINT,NTR,IOFF,JOUT,INO,INTRO,INTRV,
     .                  NOV,INV,N2OV,FLOV,ORBSYM,NIRRED,
     .                  FOFF,NVNV,FADD,NM,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 VAL,VAL1,VAL2,VAL3,VAL4
      REAL*8 MOINTS(LENINT),HMAT(NV*NV)
      INTEGER IOFF(NTR),FZO(*),FZV(*),
     .        FOFF(NVNV),FADD(NM),
     .        INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV),
     .        FLOV(NIRRED,4),ORBSYM(NO+NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
C
      CALL RFILE(TPHMAT)
      POINT=1
      DO 10 V=1,NO
         IF(FZO(V).EQ.1)GO TO 10
         VSYM=ORBSYM(V)
         DO 20 GA=1,NV
            IF(FZV(GA).EQ.1)GO TO 20
            GASYM=ORBSYM(GA+NO)
            VGASYM=IEOR(VSYM,GASYM)
            VGA=INO(GA)+V
            ICNT=0
            DO 30 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 30
               BESYM=ORBSYM(BE+NO)
               ASYM=IEOR(VGASYM,BESYM)+1
               VBE=INO(BE)+V
               FA=FLOV(ASYM,3)-NO
               LA=FLOV(ASYM,4)-NO
               DO 40 A=FA,LA
                  ICNT=ICNT+1
                  BEA=INV (MAX0(BE,A))+MIN0(BE,A)
                  GAA=INV (MAX0(GA,A))+MIN0(GA,A)
                  VGABEA=FADD (VGA)+FOFF(BEA)
                  VBEGAA=FADD (VBE)+FOFF(GAA)
                  VAL1= MOINTS(VGABEA)
                  VAL2= MOINTS(VBEGAA)
CH                HMAT(V,GA,BE,A)
C      HMAT(ICNT)=-MOINTS(VGABEA)-MOINTS(VGABEA)+MOINTS(VBEGAA)
                  HMAT(ICNT)=-VAL1-VAL1+VAL2
   40          CONTINUE
   30       CONTINUE
            CALL WWRITW(TPHMAT,ICNT,1,POINT,POINT)
            CALL WWRITW(TPHMAT,HMAT,INTOWP(ICNT),POINT,POINT)
   20    CONTINUE
   10 CONTINUE
      CALL RCLOSE(TPHMAT,3)
C
      RETURN
      END
