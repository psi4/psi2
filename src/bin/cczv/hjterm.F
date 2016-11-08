C
C-----------------------------------------------------------------------
C
      SUBROUTINE HJTERM(Z1,Z2,Z1O,TZO,TZV,HMAT,JMAT,NO,NV,NO2,NV2,NTRO,
     .                  NTRV,NTR,IOFF,JOUT,FOCK,INO,INTRO,INTRV,NOV,INV,
     .                  N2OV,FLOV,ORBSYM,NIRRED,
     .                  ZLX,TOFF,TADD,NDIMT2,FZO,FZV)
      IMPLICIT INTEGER(A-Z)
      REAL*8 Z2(NDIMT2  ),Z1(NO,NV),Z1O(NO,NV),TZO(NO,NO),TZV(NV,NV),
     .       HMAT(NV*NV),JMAT(NV*NO),FOCK(NTR),VAL
C    .       HMAT(NV*NV),JMAT(NV,NO,NO,NO),FOCK(NTR),VAL
      INTEGER IOFF(NTR),FZO(*),FZV(*)
      INTEGER INO(NOV),INTRO(N2OV),INTRV(NO*NV),INV(NOV)
      INTEGER ZLX(NV,NV),TOFF(NO,NO,2),TADD(NV,NV)
      DIMENSION FLOV(NIRRED,4),ORBSYM(NO+NV)
      COMMON/TAPES/TPBMAT,TPFMAT,TPHMAT,TPRMAT,TPQMAT,TPJMAT,ITAP71,
     .             TPAMAT
C
      CALL RFILE(TPHMAT)
      POINT=1
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         USYM=ORBSYM(U)
         DO 20 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 20
            BESYM=ORBSYM(BE+NO)
            UBESYM=IEOR(USYM,BESYM)
            CALL WREADW(TPHMAT,IDIM,1,POINT,POINT)
            CALL WREADW(TPHMAT,HMAT,INTOWP(IDIM),POINT,POINT)
            ICNT1=0
            ICNT2=0
            IF (USYM.NE.BESYM) GO TO 31
            DO 30 I=1,NO
               IF(FZO(I).EQ.1)GO TO 30
               IBE=IOFF(BE+NO)+I
               Z1(U,BE)=Z1(U,BE)+FOCK(IBE)*TZO(I,U)
   30       CONTINUE
   31       CONTINUE
            DO 50 A=1,NV
               IF(FZV(A).EQ.1)GO TO 50
               ASYM=ORBSYM(A+NO)
               UA=IOFF(A+NO)+U
               Z1(U,BE)=Z1(U,BE)+FOCK(UA)*TZV(A,BE)
               BSYM=IEOR(UBESYM,ASYM)+1
               FB=FLOV(BSYM,3)-NO
               LB=FLOV(BSYM,4)-NO
               DO 60 B=FB,LB
       ICNT1=ICNT1+1
                  Z1(U,BE)=Z1(U,BE)+HMAT(ICNT1   )*TZV(A,B)
CH                Z1(U,BE)=Z1(U,BE)+HMAT(U,BE,A,B)*TZV(A,B)
   60          CONTINUE
   50       CONTINUE
C
            UBE=IOFF(BE+NO)+U
            FI=FLOV(BESYM+1,1)
            LI=FLOV(BESYM+1,2)
               DO 80 GA=1,BE-1
               IF(FZV(GA).EQ.1)GO TO 80
               GASYM=ORBSYM(GA+NO)
            VSYM=IEOR(UBESYM,GASYM)+1
            FV=FLOV(VSYM,1)
            LV=FLOV(VSYM,2)
            FA=FLOV(VSYM,3)-NO
            LA=FLOV(VSYM,4)-NO
            DO 70 V=FV,LV
            VBE=IOFF(BE+NO)+V
            UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
            Z2(UVBEGA)=Z2(UVBEGA)
     .      -2.D0*FOCK(UBE)*Z1O(V,GA)+FOCK(VBE)*Z1O(U,GA)
   70       CONTINUE
                  DO 90 A=FA,LA
       ICNT2=ICNT2+1
                  DO 89 V=FV,LV
            UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                  Z2(UVBEGA)=Z2(UVBEGA)
     .                          +HMAT(ICNT2    )*Z1O(V,A)
CH   .                          +HMAT(U,BE,GA,A)*Z1O(V,A)
   89             CONTINUE
   90             CONTINUE
   80             CONTINUE
C
C      GO FOR GA=BE
C
               GA=BE
               GASYM=ORBSYM(GA+NO)
            VSYM=IEOR(UBESYM,GASYM)+1
            FV=FLOV(VSYM,1)
            LV=FLOV(VSYM,2)
            FA=FLOV(VSYM,3)-NO
            LA=FLOV(VSYM,4)-NO
            DO 75 V=FV,LV
            VBE=IOFF(BE+NO)+V
            VAL=  -2.D0*FOCK(UBE)*Z1O(V,GA)+FOCK(VBE)*Z1O(U,GA)
            UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
            VUBEGA=TOFF(V,U,1)+TADD(BE,GA)
            Z2(UVBEGA)=Z2(UVBEGA)+VAL
            Z2(VUBEGA)=Z2(VUBEGA)+VAL
   75       CONTINUE
                  DO 95 A=FA,LA
       ICNT2=ICNT2+1
                  DO 85 V=FV,LV
                  VAL=HMAT(ICNT2    )*Z1O(V,A)
CH                VAL=HMAT(U,BE,GA,A)*Z1O(V,A)
            UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
            VUBEGA=TOFF(V,U,1)+TADD(BE,GA)
                  Z2(UVBEGA)=Z2(UVBEGA)+VAL
                  Z2(VUBEGA)=Z2(VUBEGA)+VAL
   85             CONTINUE
   95             CONTINUE
C
            DO 81 GA=BE+1,NV
               IF(FZV(GA).EQ.1)GO TO 81
            GASYM=ORBSYM(GA+NO)
            VSYM=IEOR(UBESYM,GASYM)+1
            FV=FLOV(VSYM,1)
            LV=FLOV(VSYM,2)
            FA=FLOV(VSYM,3)-NO
            LA=FLOV(VSYM,4)-NO
            DO 71 V=FV,LV
            VBE=IOFF(BE+NO)+V
            VUGABE=TOFF(V,U,1)+TADD(GA,BE)
            Z2(VUGABE)=Z2(VUGABE)
     .                -2.D0*FOCK(UBE)*Z1O(V,GA)+FOCK(VBE)*Z1O(U,GA)
   71       CONTINUE
                  DO 91 A=FA,LA
       ICNT2=ICNT2+1
                  DO 72 V=FV,LV
            VUGABE=TOFF(V,U,1)+TADD(GA,BE)
                  Z2(VUGABE)=Z2(VUGABE)+HMAT(ICNT2    )*Z1O(V,A)
CH                Z2(VU,GABE)=Z2(VU,GABE)+HMAT(U,BE,GA,A)*Z1O(V,A)
   72             CONTINUE
   91             CONTINUE
   81          CONTINUE
   20    CONTINUE
   10 CONTINUE
      CALL RCLOSE(TPHMAT,3)
C
C
      CALL RFILE(TPJMAT)
      POINT=1
      DO 300 V=1,NO
      IF(FZO(V).EQ.1)GO TO 300
      VSYM=ORBSYM(V)
      DO 310 I=1,NO
      IF(FZO(I).EQ.1)GO TO 310
      ISYM=ORBSYM(I)
      VISYM=IEOR(VSYM,ISYM)
      CALL WREADW(TPJMAT,IDIM,1,POINT,POINT)
      CALL WREADW(TPJMAT,JMAT,INTOWP(IDIM),POINT,POINT)
      ICNT=0
      DO 270 U=1,NO
      IF(FZO(U).EQ.1)GO TO 270
      USYM=ORBSYM(U)
         BESYM=IEOR(VISYM,USYM)+1
         FBE=FLOV(BESYM,3)-NO
         LBE=FLOV(BESYM,4)-NO
         DO 320 BE=FBE,LBE
      ICNT=ICNT+1
C
C              Z1(U,BE)=Z1(U,BE)+JMAT(BE,U,I,J)*TZO(I,J)
         Z1(U,BE)=Z1(U,BE)+JMAT(ICNT    )*TZO(V,I)
CJ       Z1(U,BE)=Z1(U,BE)+JMAT(BE,U,V,I)*TZO(V,I)
           FGA=FLOV(ISYM+1,3)-NO
           LGA=FLOV(ISYM+1,4)-NO
           FGA2=FGA
           LGA2=LGA
           IF(LGA.GT.BE)LGA2=BE
           DO 280 GA=FGA,LGA2
           UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
           Z2(UVBEGA)=Z2(UVBEGA)+JMAT(ICNT    )*Z1O(I,GA)
  280      CONTINUE
           IF(FGA.LT.BE)FGA2=BE
           DO 281 GA=FGA2,LGA
           VUGABE=TOFF(V,U,1)+TADD(GA,BE)
           Z2(VUGABE)=Z2(VUGABE)+JMAT(ICNT)*Z1O(I,GA)
  281      CONTINUE
C
  320    CONTINUE
  270 CONTINUE
  310 CONTINUE
  300 CONTINUE
      CALL RCLOSE(TPJMAT,3)
      RETURN
      END
