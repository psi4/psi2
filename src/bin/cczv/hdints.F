C
C-----------------------------------------------------------------------
C
      SUBROUTINE HDINTS(HMAT,MOINTS,NO2,NV2,T1,NTRO,NTRV,NO,NV,LENINT,
     .                 NTR,IOFF,JOUT,INO,INTRO,INTRV,NOV,INV,N2OV,FLOV,
     .                 ORBSYM,NIRRED,
     .                 DOFF,NOFF,DADD,NADD,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 T1(NO,NV),HMAT(NV*NV),VAL,VAL2,MOINTS(LENINT),VALK
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER DOFF(NOFF),DADD(NADD)
      DIMENSION FLOV(NIRRED,4),ORBSYM(NO+NV)
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
            GAV=INO(GA)+V
      CALL WREADW(TPHMAT,IDIM,1,POINT,POINT)
      CALL WREADW(TPHMAT,HMAT,INTOWP(IDIM),POINT,IDUM )
      ICNT=0
            DO 30 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 30
               BESYM=ORBSYM(BE+NO)
               ISYM=IEOR(VGASYM,BESYM)+1
               FI=FLOV(ISYM,1)
               LI=FLOV(ISYM,2)
               FA=FLOV(ISYM,3)-NO
               LA=FLOV(ISYM,4)-NO
               BEV=INO(BE)+V
               DO 40 A=FA,LA
                  VAL=0.0D+00
                  VAL2=0.0D+00
                  ICNT=ICNT+1
                  DO 50 I=FI,LI
                     BEI=INO(BE)+I
                     GAI=INO(GA)+I
                     BEIGAV=DOFF(MAX0(BEI,GAV))+DADD(MIN0(BEI,GAV))
                     BEVGAI=DOFF(MAX0(BEV,GAI))+DADD(MIN0(BEV,GAI))
                     VAL=VAL+MOINTS(BEIGAV)*T1(I,A)
                     VAL2=VAL2+MOINTS(BEVGAI)*T1(I,A)
   50             CONTINUE
CH                HMAT(V,GA,BE,A)=HMAT(V,GA,BE,A)+VAL+VAL-VAL2
                  HMAT(ICNT     )=HMAT(ICNT     )+VAL+VAL-VAL2
   40         CONTINUE
   30      CONTINUE
           CALL WWRITW(TPHMAT,HMAT,INTOWP(IDIM),POINT,POINT)
   20    CONTINUE
   10 CONTINUE
      CALL RCLOSE(TPHMAT,3)
      RETURN
      END
