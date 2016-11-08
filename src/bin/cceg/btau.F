C
C ----------------------------------------------------------------------
C
      SUBROUTINE BTAU  (T1O,T2O,T2,FIN,AR3,AR4,BIN,BUFINT,
     .                  SUX1,SUX2,SUX3,SUX4,
     .                  ORBSYM,FLOV,FOFF,FADD,NAB,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,MAXOO,MAXVV,ITAP71)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL1,VAL2,TMP,TMV
      REAL*8  T1O(NO,NV),T2O(*),T2(*),FIN(*),BIN(NV,NV),BUFINT(*),
     .        SUX1(*),SUX2(*),SUX3(*),SUX4(*),
     .        AR3(MAXVV,MAXOO),AR4(MAXVV,MAXOO)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO,2),TADD(NV,NV),ZLX(NV,NV),
     .        FOFF(*),FADD(*),NAB(*)
C
C     NAB MUST INCLUDE FROZEN PAIRS
C     BECAUSE FILE71 CONTAINS (VV,VV) INTEGRALS FOR ALL PAIRS
C     POINTER MUST BE INCREMENTED ALSO FOR FROZEN PAIRS!
C
      CALL RFILE(ITAP71)
      POINT=1
C
      DO 992 TSYM=1,NIRRED
         LAB=0
         DO 123 A=1,NV
C           IF(FZV(A).EQ.1)GO TO 123
            ASYM=ORBSYM(A+NO)
            DO 121 B=1,NV
C              IF(FZV(B).EQ.1)GO TO 121
               BSYM=ORBSYM(B+NO)
               ABSYM=IEOR(ASYM,BSYM)
               IF(ABSYM.NE.TSYM-1)GO TO 121
               LAB=LAB+1
  121       CONTINUE
  123    CONTINUE
         NAB(TSYM)=LAB
C
C        WRITE(6,*)'TSYM',TSYM, ' NAB(TSYM)=',NAB(TSYM)
C
         XAB=0
         DO 290 A=1,NV
            IF(FZV(A).EQ.1)GO TO 290
            ASYM=ORBSYM(A+NO)
            DO 280 B=1,A
               IF(FZV(B).EQ.1)GO TO 280
               BSYM=ORBSYM(B+NO)
               BASYM=IEOR(BSYM,ASYM)
               IF(BASYM.NE.TSYM-1)GO TO 280
               XAB=XAB+1
               XIJ=0
               DO 270 I=1,NO
                  IF(FZO(I).EQ.1)GO TO 270
                  ISYM=ORBSYM(I)
                  JSYM=IEOR(ISYM,BASYM)+1
                  FJ=FLOV(JSYM,1)
                  LJ=FLOV(JSYM,2)
                  IF(LJ.GT.I)LJ=I
                  DO 260 J=FJ,LJ
                     XIJ=XIJ+1
                     IJAB=TOFF(I,J,1)+TADD(A,B)
                     JIAB=TOFF(J,I,1)+TADD(A,B)
                     VAL1=T2O(IJAB)+T1O(I,A)*T1O(J,B)
                     VAL2=T2O(JIAB)+T1O(I,B)*T1O(J,A)
                     AR3(XAB,XIJ)=VAL1+VAL2
                     AR4(XAB,XIJ)=VAL1-VAL2
  260             CONTINUE
  270          CONTINUE
  280       CONTINUE
  290    CONTINUE
         LIMO=XIJ
         LIMV=XAB
C
C        WRITE(6,*)'TSYM=',TSYM,' XIJ=',XIJ,' XAB=',XAB,'   ',XIJ*XAB
C
         DO 991 BE=1,NV
C           IF(FZV(BE).EQ.1)GO TO 991
            BESYM=ORBSYM(BE+NO)
            FI1=FLOV(BESYM+1,1)
            LI1=FLOV(BESYM+1,2)
            DO 981 GA=1,BE
C              IF(FZV(GA).EQ.1)GO TO 981
               GASYM=ORBSYM(GA+NO)
               FI2=FLOV(GASYM+1,1)
               LI2=FLOV(GASYM+1,2)
               BEGAS=IEOR(GASYM,BESYM)
               IF(BEGAS.NE.TSYM-1)GO TO 981
               LAB=NAB(TSYM)
               IF(FZV(BE).EQ.1.OR.FZV(GA).EQ.1)THEN
                  POINT=POINT+INTOWP(LAB)
                  GO TO 981
               ENDIF
               CALL WREADW(ITAP71,BUFINT,INTOWP(LAB),POINT,POINT)
               CALL ZERO(BIN,NV*NV)
               NADD=0
               DO 325 A=1,NV
C                 IF(FZV(A).EQ.1)GO TO 325
                  ASYM=ORBSYM(A+NO)
                  DO 320 B=1,NV
C                    IF(FZV(B).EQ.1)GO TO 320
                     BSYM=ORBSYM(B+NO)
                     ABSYM=IEOR(ASYM,BSYM)
                     IF(ABSYM.NE.TSYM-1)GO TO 320
                     NADD=NADD+1
                     BIN(A,B)=BUFINT(NADD)
  320             CONTINUE
  325          CONTINUE
C
               XBA=0
               DO 425 A=1,NV
                  IF(FZV(A).EQ.1)GO TO 425
                  ASYM=ORBSYM(A+NO)
                  BEA=ITV(MAX0(BE,A))+MIN0(BE,A)
                  GAA=ITV(MAX0(GA,A))+MIN0(GA,A)
                  BSYM=IEOR(ASYM,BEGAS)
                  FB=FLOV(BSYM+1,3)-NO
                  LB=FLOV(BSYM+1,4)-NO
                  IF(LB.GT.A)LB=A
                  DO 420 B=FB,LB
                     XBA=XBA+1
                     GAB=ITV(MAX0(GA,B))+MIN0(GA,B)
                     BEB=ITV(MAX0(BE,B))+MIN0(BE,B)
                     TMP=BIN(A,B)+BIN(B,A)
                     TMV=BIN(A,B)-BIN(B,A)
                     VAL1=0.0D0
                     VAL2=0.0D0
                     DO 408 I=FI2,LI2
                        IA=ITR(A)+I
                        IABEB=FOFF(BEB)+FADD(IA)
                        VAL2=VAL2-FIN(IABEB)*T1O(I,GA)
                        IB=ITR(B)+I
                        IBBEA=FOFF(BEA)+FADD(IB)
                        VAL1=VAL1-FIN(IBBEA)*T1O(I,GA)
  408                CONTINUE
                     DO 410 I=FI1,LI1
                        IA=ITR(A)+I
                        IAGAB=FOFF(GAB)+FADD(IA)
                        VAL1=VAL1-FIN(IAGAB)*T1O(I,BE)
                        IB=ITR(B)+I
                        IBGAA=FOFF(GAA)+FADD(IB)
                        VAL2=VAL2-FIN(IBGAA)*T1O(I,BE)
  410                CONTINUE
                     SUX3(XBA)=(TMP+VAL1+VAL2)*0.5D0
                     SUX4(XBA)=(TMV+VAL1-VAL2)*0.5D0
  420             CONTINUE
                  IF(ASYM.EQ.BSYM)SUX3(XBA)=SUX3(XBA)*0.5D0
  425          CONTINUE
C
C              GO FOR THE N**6 LOOP
C
               DO 577 UV=1,LIMO
                  TMP = 0.0D0
                  TMV = 0.0D0
                  DO 529 AB=1,LIMV
                     TMP=TMP+SUX3(AB)* AR3(AB,UV)
                     TMV=TMV+SUX4(AB)* AR4(AB,UV)
  529             CONTINUE
                  SUX1(UV)=TMP
                  SUX2(UV)=TMV
  577          CONTINUE
C
               UV=0
               DO 668 U=1,NO
                  IF(FZO(U).EQ.1)GO TO 668
                  USYM=ORBSYM(U)
                  VSYM=IEOR(USYM,BEGAS)+1
                  SV=FLOV(VSYM,1)
                  LV=FLOV(VSYM,2)
                  IF(LV.GT.U )LV=U
                  DO 663 V=SV,LV
                     UV=UV+1
                     UVBEGA=TOFF(U,V,1)+TADD(BE,GA)
                     T2(UVBEGA)=T2(UVBEGA)+SUX1(UV)+SUX2(UV)
                     IF(U.NE.V)THEN
                        VUBEGA=TOFF(V,U,1)+TADD(BE,GA)
                        T2(VUBEGA)=T2(VUBEGA)+SUX1(UV)-SUX2(UV)
                     ENDIF
  663             CONTINUE
  668          CONTINUE
C
  981       CONTINUE
  991    CONTINUE
  992 CONTINUE
      CALL RCLOSE(ITAP71,3)
C
      RETURN
      END
