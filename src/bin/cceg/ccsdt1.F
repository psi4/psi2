C
C CCSDT-1A SUBROUTINE. WRITTEN BY GUSTAVO E. SCUSERIA. OCTOBER 7, 1988.
C
C***********************************************************************
C  NOTICE OF PROGRAM MODIFICATION
c 
C      BY:  gustavo e scuseria 
C    DATE:  august 31, 1989 
C  REASON:  get rid of file81                        
c
C      BY:  GUSTAVO E. SCUSERIA
C    DATE:  NOVEMBER 2, 1988
C  REASON:  RCLOSE (ITAP81,4) AFTER UPDATING T1 & T2
c
C***********************************************************************
      SUBROUTINE CCSDT1( S,T1N,T2 ,T1,TAU,NO,NV,NDIMT2,NIT,
     .                 ORBSYM,FLOV,NIRRED,TOFF,TADD,ZLX,FZO,FZV,
     .                 UOFF,VADD,OPTION,RSTR,NT3,DOFF,DADD,EOFF,EADD,
     .                 FOFF,FADD,ITV,ITR,IPQ,DIN,EIN,FIN,fl81,FOCK,
     .                 AUX1,AUX2,AUX3,SUX1,SUX2,SUX3,
     .                 XAR1,XAR2,NOACT,NVACT,MAXO,MAXV,MAXOV,
     .                 AR1,AR2,SUX11,SUX22,XAR3,XAR4,AUX4,AUX5)
      IMPLICIT INTEGER(A-Z)
      REAL*8  S(NDIMT2),T1N(NO,NV),T2(NDIMT2),T1(NO,NV),TAU(*),FOCK(*),
     .       VAL,VAL1,VAL2,DIN(*),EIN(*),FIN(*),
     .       AUX1(NV,NO),AUX2(NV,NO),AUX3(NV,NV),SUX11(NO*NV),
     .       AUX4(NO,NO),AUX5(NO,NV),
     .       SUX1(NV*NO),SUX2(NV*NO),SUX3(NV*NV),SUX22(NO*NV),
     .       XAR1(MAXV,MAXOV),XAR2(MAXV,MAXOV),
     .       XAR3(MAXO,MAXOV),XAR4(MAXO,MAXOV),
     .       AR1(MAXOV,MAXOV),AR2(MAXOV,MAXOV)
      CHARACTER*4 OPTION,RSTR,SORT,fl81
      INTEGER ORBSYM(*),FLOV(NIRRED,4),TOFF(NO,NO,2,NIRRED),ZLX(NV,NV),
     .        TADD(NV,NV,NIRRED),FZO(NO),FZV(NV),UOFF(NO,NO,2),NT3(8),
     .        VADD(NV,NV),DOFF(*),DADD(*),FOFF(*),FADD(*),EOFF(*),
     .        EADD(*),ITR(*),ITV(*),IPQ(*)
      integer prcntr
C
      SORT='NO  '
      IF(NIT.EQ.3)THEN
         if (prcntr('IS_ON BRIEF').eq.0)
     &     WRITE(6,*)'ENTERING CCSDT-1A SUBROUTINE'
C         if (prcntr('IS_ON BRIEF').eq.0) CALL TIMIT(1,6)
      END IF
c
      itap81=81
      call rfile(itap81)
      i81=1
C
      DO 9990 A=1,NV
         IF(FZV(A).EQ.1)GO TO 9990
         XA=IPQ(A+NO)+A+NO
         ASYM=ORBSYM(A+NO)
         DO 9985 I=1,NO
            IF(FZO(I).EQ.1)GO TO 9985
            XI=IPQ(I)+I
            ISYM=ORBSYM(I)
            AI=ITR(A)+I
            AISYM=IEOR(ISYM,ASYM)
            AIS=AISYM+1
            DIM=NT3(AIS)
            CALL ZERO (TAU,DIM)
C
C  F INTS contribution to t3.
C
            DO 3910 B=1,NV
               IF(FZV(B).EQ.1)GO TO 3910
               BSYM=ORBSYM(B+NO)
               ESYM=IEOR(BSYM,AISYM)+1
               FE=FLOV(ESYM,3)-NO
               LE=FLOV(ESYM,4)-NO
               DO 3905 E=FE,LE
                  BE=ITV(MAX0(B,E))+MIN0(B,E)
                  AIBE=FOFF(BE)+FADD(AI)
                  AUX3(E,B)=FIN(AIBE)
 3905          CONTINUE
 3910       CONTINUE
            DO 3980 C=1,NV
               IF(FZV(C).EQ.1)GO TO 3980
               CSYM=ORBSYM(C+NO)
               DO 3975 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 3975
                  KSYM=ORBSYM(K)
                  CKSYM=IEOR(CSYM,KSYM)
                  AICKS=IEOR(AISYM,CKSYM)
                  DO 3920 E=1,NV
                     IF(FZV(E).EQ.1)GO TO 3920
                     ZLCE=ZLX(C,E)
                     ESYM=ORBSYM(E+NO)
                     JSYM=IEOR(ESYM,CKSYM)+1
                     FJ=FLOV(JSYM,1)
                     LJ=FLOV(JSYM,2)
                     DO 3915 J=FJ,LJ
                        KJCE=UOFF(K,J,ZLCE)+VADD(C,E)
                        AUX1(E,J)=T2(KJCE)
 3915                CONTINUE
 3920             CONTINUE
                  DO 3935 B=1,C
                     IF(FZV(B).EQ.1)GO TO 3935
                     BSYM=ORBSYM(B+NO)
                     JSYM=IEOR(BSYM,AICKS)+1
                     FJ=FLOV(JSYM,1)
                     LJ=FLOV(JSYM,2)
                     ESYM=IEOR(BSYM,AISYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 3930 J=FJ,LJ
                     VAL=0.0D0
                        DO 3925 E=FE,LE
                        VAL=VAL+AUX3(E,B)*AUX1(E,J)
 3925                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)+VAL
 3930                CONTINUE
 3935             CONTINUE
                  DO 3950 B=C,NV
                     IF(FZV(B).EQ.1)GO TO 3950
                     BSYM=ORBSYM(B+NO)
                     JSYM=IEOR(BSYM,AICKS)+1
                     FJ=FLOV(JSYM,1)
                     LJ=FLOV(JSYM,2)
                     ESYM=IEOR(BSYM,AISYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 3945 J=FJ,LJ
                     VAL=0.0D0
                        DO 3940 E=FE,LE
                        VAL=VAL+AUX3(E,B)*AUX1(E,J)
 3940                   CONTINUE
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)+VAL
 3945                CONTINUE
 3950             CONTINUE
 3975          CONTINUE
 3980       CONTINUE
C
C           IF(SORT.NE.'YES ')THEN
c
            DO 4110 K=1,NO
               IF(FZO(K).EQ.1)GO TO 4110
               KSYM=ORBSYM(K)
               ESYM=IEOR(KSYM,AISYM)+1
               FE=FLOV(ESYM,3)-NO
               LE=FLOV(ESYM,4)-NO
               DO 4105 E=FE,LE
                  ZLAE=ZLX(A,E)
                  IKAE=UOFF(I,K,ZLAE)+VADD(A,E)
                  AUX1(E,K)=T2(IKAE)
 4105          CONTINUE
 4110       CONTINUE
            DO 4180 B=1,NV
               IF(FZV(B).EQ.1)GO TO 4180
               BSYM=ORBSYM(B+NO)
               DO 4175 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 4175
                  BJ=ITR(B)+J
                  JSYM=ORBSYM(J)
                  BJSYM=IEOR(BSYM,JSYM)
                  AIBJS=IEOR(AISYM,BJSYM)
                  DO 4120 C=1,NV
                     IF(FZV(C).EQ.1)GO TO 4120
                     CSYM=ORBSYM(C+NO)
                     ESYM=IEOR(CSYM,BJSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4115 E=FE,LE
                        CE=ITV(MAX0(C,E))+MIN0(C,E)
                        BJCE=FOFF(CE)+FADD(BJ)
                        AUX3(E,C)=FIN(BJCE)
 4115                CONTINUE
 4120             CONTINUE
                  DO 4135 C=1,B
                     IF(FZV(C).EQ.1)GO TO 4135
                     CSYM=ORBSYM(C+NO)
                     KSYM=IEOR(CSYM,AIBJS)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     ESYM=IEOR(CSYM,BJSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4130 K=FK,LK
                     VAL=0.0D0
                        DO 4125 E=FE,LE
                        VAL=VAL+AUX3(E,C)*AUX1(E,K)
 4125                   CONTINUE
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)+VAL
 4130                CONTINUE
 4135             CONTINUE
                  DO 4150 C=B,NV
                     IF(FZV(C).EQ.1)GO TO 4150
                     CSYM=ORBSYM(C+NO)
                     KSYM=IEOR(CSYM,AIBJS)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     ESYM=IEOR(CSYM,BJSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4145 K=FK,LK
                     VAL=0.0D0
                        DO 4140 E=FE,LE
                        VAL=VAL+AUX3(E,C)*AUX1(E,K)
 4140                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)+VAL
 4145                CONTINUE
 4150             CONTINUE
 4175          CONTINUE
 4180       CONTINUE
C
            DO 4292 XSYM=1,NIRRED
            BJSEM=IEOR(AISYM,XSYM-1)
            XCK=0
            DO 4204 C=1,NV
               IF(FZV(C).EQ.1)GO TO 4204
               CSYM=ORBSYM(C+NO)
               DO 4203 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 4203
                  KSYM=ORBSYM(K)
                  CKSYM=IEOR(CSYM,KSYM)
                  IF(CKSYM.NE.XSYM-1)GO TO 4203
                  CK=ITR(C)+K
                  XCK=XCK+1
                  ESYM=IEOR(CKSYM,ASYM)+1
                  FE=FLOV(ESYM,3)-NO
                  LE=FLOV(ESYM,4)-NO
                  DO 4202 E=FE,LE
                     AE=ITV(MAX0(A,E))+MIN0(A,E)
                     CKAE=FOFF(AE)+FADD(CK)
                     XAR1(E,XCK)=FIN(CKAE)
 4202             CONTINUE
 4203          CONTINUE
 4204       CONTINUE
            XBJ=0
            DO 4208 B=1,NV
               IF(FZV(B).EQ.1)GO TO 4208
               BSYM=ORBSYM(B+NO)
               DO 4207 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 4207
                  JSYM=ORBSYM(J)
                  BJSYM=IEOR(BSYM,JSYM)
                  IF(BJSYM.NE.BJSEM)GO TO 4207
                  XBJ=XBJ+1
                  ESYM=IEOR(BJSYM,ISYM)+1
                  FE=FLOV(ESYM,3)-NO
                  LE=FLOV(ESYM,4)-NO
                  DO 4206 E=FE,LE
                     ZLBE=ZLX(B,E)
                     JIBE=UOFF(J,I,ZLBE)+VADD(B,E)
                     XAR2(E,XBJ)=T2(JIBE)
 4206             CONTINUE
 4207          CONTINUE
 4208       CONTINUE
            XBJ=0
            DO 4280 B=1,NV
               IF(FZV(B).EQ.1)GO TO 4280
               BSYM=ORBSYM(B+NO)
               DO 4275 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 4275
                  JSYM=ORBSYM(J)
                  BJSYM=IEOR(BSYM,JSYM)
                  IF(BJSYM.NE.BJSEM)GO TO 4275
                  XBJ=XBJ+1
                  XCK=0
                  DO 4235 C=1,B-1
                     IF(FZV(C).EQ.1)GO TO 4235
                     CSYM=ORBSYM(C+NO)
                     KSYM=IEOR(CSYM,XSYM-1)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     ESYM=IEOR(ISYM,BJSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4230 K=FK,LK
                     VAL=0.0D0
                     XCK=XCK+1
                        DO 4225 E=FE,LE
                        VAL=VAL+XAR1(E,XCK)*XAR2(E,XBJ)
 4225                   CONTINUE
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)+VAL
 4230                CONTINUE
 4235             CONTINUE
                  DO 4239 C=B,B
                     IF(FZV(C).EQ.1)GO TO 4239
                     CSYM=ORBSYM(C+NO)
                     KSYM=IEOR(CSYM,XSYM-1)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     ESYM=IEOR(ISYM,BJSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4238 K=FK,LK
                     XCK=XCK+1
                     VAL=0.0D0
                        DO 4237 E=FE,LE
                        VAL=VAL+XAR1(E,XCK)*XAR2(E,XBJ)
 4237                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)+VAL
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)+VAL
 4238                CONTINUE
 4239             CONTINUE
                  DO 4250 C=B+1,NV
                     IF(FZV(C).EQ.1)GO TO 4250
                     CSYM=ORBSYM(C+NO)
                     KSYM=IEOR(CSYM,XSYM-1)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     ESYM=IEOR(ISYM,BJSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4245 K=FK,LK
                     XCK=XCK+1
                     VAL=0.0D0
                        DO 4240 E=FE,LE
                        VAL=VAL+XAR1(E,XCK)*XAR2(E,XBJ)
 4240                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)+VAL
 4245                CONTINUE
 4250             CONTINUE
 4275          CONTINUE
 4280       CONTINUE
 4292       CONTINUE
C
C  E INTS contribution to t3.
C
            DO 4910 J=1,NO
               IF(FZO(J).EQ.1)GO TO 4910
               JSYM=ORBSYM(J)
               MSYM=IEOR(JSYM,AISYM)+1
               FM=FLOV(MSYM,1)
               LM=FLOV(MSYM,2)
               DO 4905 M=FM,LM
                  JM=ITR(MAX0(J,M))+MIN0(J,M)
                  AIJM=EOFF(AI)+EADD(JM)
                  AUX4(M,J)=EIN(AIJM)
 4905          CONTINUE
 4910       CONTINUE
            DO 4980 C=1,NV
               IF(FZV(C).EQ.1)GO TO 4980
               CSYM=ORBSYM(C+NO)
               DO 4975 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 4975
                  KSYM=ORBSYM(K)
                  CKSYM=IEOR(CSYM,KSYM)
                  AICKS=IEOR(AISYM,CKSYM)
                  DO 4920 M=1,NO
                     IF(FZO(M).EQ.1)GO TO 4920
                     MSYM=ORBSYM(M)
                     BSYM=IEOR(MSYM,CKSYM)+1
                     FB=FLOV(BSYM,3)-NO
                     LB=FLOV(BSYM,4)-NO
                     DO 4915 B=FB,LB
                        ZLCB=ZLX(C,B)
                        KMCB=UOFF(K,M,ZLCB)+VADD(C,B)
                        AUX5(M,B)=T2(KMCB)
 4915                CONTINUE
 4920             CONTINUE
                  DO 4935 J=1,NO
                     IF(FZO(J).EQ.1)GO TO 4935
                     JSYM=ORBSYM(J)
                     BSYM=IEOR(JSYM,AICKS)+1
                     FB=FLOV(BSYM,3)-NO
                     LB=FLOV(BSYM,4)-NO
                     FB2=FB
                     LB2=LB
                     MSYM=IEOR(JSYM,AISYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     IF(LB.GT.C)LB2=C
                     DO 4930 B=FB2,LB2
                     VAL=0.0D0
                        DO 4925 M=FM,LM
                        VAL=VAL+AUX4(M,J)*AUX5(M,B)
 4925                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)-VAL
 4930                CONTINUE
                     FB2=FB
                     LB2=LB
                     IF(FB.LT.C)FB2=C
                     DO 4931 B=FB2,LB2
                     VAL=0.0D0
                        DO 4926 M=FM,LM
                        VAL=VAL+AUX4(M,J)*AUX5(M,B)
 4926                   CONTINUE
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)-VAL
 4931                CONTINUE
 4935             CONTINUE
 4950             CONTINUE
 4975          CONTINUE
 4980       CONTINUE
C
C           IF(SORT.NE.'YES ')THEN
c
            DO 5110 C=1,NV
               IF(FZV(C).EQ.1)GO TO 5110
               ZLAC=ZLX(A,C)
               CSYM=ORBSYM(C+NO)
               MSYM=IEOR(CSYM,AISYM)+1
               FM=FLOV(MSYM,1)
               LM=FLOV(MSYM,2)
               DO 5105 M=FM,LM
                  IMAC=UOFF(I,M,ZLAC)+VADD(A,C)
                  AUX5(M,C)=T2(IMAC)
 5105          CONTINUE
 5110       CONTINUE
            DO 5180 B=1,NV
               IF(FZV(B).EQ.1)GO TO 5180
               BSYM=ORBSYM(B+NO)
               DO 5175 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 5175
                  BJ=ITR(B)+J
                  JSYM=ORBSYM(J)
                  BJSYM=IEOR(BSYM,JSYM)
                  AIBJS=IEOR(AISYM,BJSYM)
                  DO 5120 K=1,NO
                     IF(FZO(K).EQ.1)GO TO 5120
                     KSYM=ORBSYM(K)
                     MSYM=IEOR(KSYM,BJSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5115 M=FM,LM
                        KM=ITR(MAX0(K,M))+MIN0(K,M)
                        BJKM=EOFF(BJ)+EADD(KM)
                        AUX4(M,K)=EIN(BJKM)
 5115                CONTINUE
 5120             CONTINUE
                  DO 5135 K=1,NO
                     IF(FZO(K).EQ.1)GO TO 5135
                     KSYM=ORBSYM(K)
                     CSYM=IEOR(KSYM,AIBJS)+1
                     FC=FLOV(CSYM,3)-NO
                     LC=FLOV(CSYM,4)-NO
                     FC2=FC
                     LC2=LC
                     MSYM=IEOR(KSYM,BJSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     IF(LC.GT.B)LC2=B
                     DO 5130 C=FC2,LC2
                     VAL=0.0D0
                        DO 5125 M=FM,LM
                        VAL=VAL+AUX4(M,K)*AUX5(M,C)
 5125                   CONTINUE
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)-VAL
 5130                CONTINUE
                     FC2=FC
                     LC2=LC
                     IF(FC.LT.B)FC2=B
                     DO 5131 C=FC2,LC2
                     VAL=0.0D0
                        DO 5126 M=FM,LM
                        VAL=VAL+AUX4(M,K)*AUX5(M,C)
 5126                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)-VAL
 5131                CONTINUE
 5135             CONTINUE
 5175          CONTINUE
 5180       CONTINUE
C
            DO 5292 XSYM=1,NIRRED
            BJSEM=IEOR(AISYM,XSYM-1)
            XCK=0
            DO 5204 C=1,NV
               IF(FZV(C).EQ.1)GO TO 5204
               CSYM=ORBSYM(C+NO)
               DO 5203 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 5203
                  KSYM=ORBSYM(K)
                  CKSYM=IEOR(CSYM,KSYM)
                  IF(CKSYM.NE.XSYM-1)GO TO 5203
                  CK=ITR(C)+K
                  XCK=XCK+1
                  MSYM=IEOR(CKSYM,ISYM)+1
                  FM=FLOV(MSYM,1)
                  LM=FLOV(MSYM,2)
                  DO 5202 M=FM,LM
                     IM=ITR(MAX0(I,M))+MIN0(I,M)
                     CKIM=EOFF(CK)+EADD(IM)
                     XAR3(M,XCK)=EIN(CKIM)
 5202             CONTINUE
 5203          CONTINUE
 5204       CONTINUE
            XBJ=0
            DO 5208 B=1,NV
               IF(FZV(B).EQ.1)GO TO 5208
               ZLBA=ZLX(B,A)
               BSYM=ORBSYM(B+NO)
               DO 5207 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 5207
                  JSYM=ORBSYM(J)
                  BJSYM=IEOR(BSYM,JSYM)
                  IF(BJSYM.NE.BJSEM)GO TO 5207
                  XBJ=XBJ+1
                  MSYM=IEOR(BJSYM,ASYM)+1
                  FM=FLOV(MSYM,1)
                  LM=FLOV(MSYM,2)
                  DO 5206 M=FM,LM
                     JMBA=UOFF(J,M,ZLBA)+VADD(B,A)
                     XAR4(M,XBJ)=T2(JMBA)
 5206             CONTINUE
 5207          CONTINUE
 5208       CONTINUE
            XBJ=0
            DO 5280 B=1,NV
               IF(FZV(B).EQ.1)GO TO 5280
               BSYM=ORBSYM(B+NO)
               DO 5275 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 5275
                  JSYM=ORBSYM(J)
                  BJSYM=IEOR(BSYM,JSYM)
                  IF(BJSYM.NE.BJSEM)GO TO 5275
                  XBJ=XBJ+1
                  XCK=0
                  DO 5235 C=1,B-1
                     IF(FZV(C).EQ.1)GO TO 5235
                     CSYM=ORBSYM(C+NO)
                     KSYM=IEOR(CSYM,XSYM-1)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     MSYM=IEOR(ASYM,BJSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5230 K=FK,LK
                     VAL=0.0D0
                     XCK=XCK+1
                        DO 5225 M=FM,LM
                        VAL=VAL+XAR3(M,XCK)*XAR4(M,XBJ)
 5225                   CONTINUE
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)-VAL
 5230                CONTINUE
 5235             CONTINUE
                  DO 5249 C=B,B
                     IF(FZV(C).EQ.1)GO TO 5249
                     CSYM=ORBSYM(C+NO)
                     KSYM=IEOR(CSYM,XSYM-1)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     MSYM=IEOR(ASYM,BJSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5238 K=FK,LK
                     XCK=XCK+1
                     VAL=0.0D0
                        DO 5237 M=FM,LM
                        VAL=VAL+XAR3(M,XCK)*XAR4(M,XBJ)
 5237                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)-VAL
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)-VAL
 5238                CONTINUE
 5249             CONTINUE
                  DO 5250 C=B+1,NV
                     IF(FZV(C).EQ.1)GO TO 5250
                     CSYM=ORBSYM(C+NO)
                     KSYM=IEOR(CSYM,XSYM-1)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     MSYM=IEOR(ASYM,BJSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5245 K=FK,LK
                     XCK=XCK+1
                     VAL=0.0D0
                        DO 5240 M=FM,LM
                        VAL=VAL+XAR3(M,XCK)*XAR4(M,XBJ)
 5240                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)-VAL
 5245                CONTINUE
 5250             CONTINUE
 5275          CONTINUE
 5280       CONTINUE
 5292       CONTINUE
C
C  DIVIDE BY MO ENERGIES
C
            DO 3570 J=1,NO
            IF(FZO(J).EQ.1)GO TO 3570
            XJ=IPQ(J)+J
            JSYM=ORBSYM(J)
            IJSYM=IEOR(ISYM,JSYM)
            DO 3560 B=1,NV
               IF(FZV(B).EQ.1)GO TO 3560
               XB=IPQ(B+NO)+B+NO
               BSYM=ORBSYM(B+NO)
               ABSYM=IEOR(ASYM,BSYM)
               DO 3550 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 3550
                  XK=IPQ(K)+K
                  KSYM=ORBSYM(K)
                  IJKSYM=IEOR(IJSYM,KSYM)
                  CSYM=IEOR(IJKSYM,ABSYM)+1
                  FC=FLOV(CSYM,3)-NO
                  LC=FLOV(CSYM,4)-NO
                  LC2=LC
                  IF(LC.GT.B)LC2=B
                  DO 3552 C=FC,LC2
                     IF(FZV(C).EQ.1)GO TO 3552
                     XC=IPQ(C+NO)+C+NO
                     ZLBC=ZLX(B,C)
                     JKBC=TOFF(J,K,ZLBC,AIS)+TADD(B,C,AIS)
                     VAL=FOCK(XI)+FOCK(XJ)+FOCK(XK)
     .                  -FOCK(XA)-FOCK(XB)-FOCK(XC)
                     TAU(JKBC)=TAU(JKBC)/VAL
 3552                CONTINUE
 3550             CONTINUE
 3560          CONTINUE
 3570       CONTINUE
c
      if(fl81.eq.'GRAD')then
         call wwritw(itap81,tau,intowp(dim),i81,i81)
         go to 9985
      endif
C
C >>> ADD CCSDT TERMS TO T1 EQUATION. D INTS.
C
      IA=ITR(A)+I
      DO 3180 U=1,NO
      IF(FZO(U).EQ.1)GO TO 3180
      USYM=ORBSYM(U)
      BESYM=USYM
      FBE=FLOV(BESYM+1,3)-NO
      LBE=FLOV(BESYM+1,4)-NO
      DO 3175 BE=FBE,LBE
      DO 3170 J=1,NO
      IF(FZO(J).EQ.1)GO TO 3170
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      JA=ITR(A)+J
      BSYM=IEOR(IJSYM,ASYM)+1
      FB=FLOV(BSYM,3)-NO
      LB=FLOV(BSYM,4)-NO
      DO 3165 B=FB,LB
      IB=ITR(B)+I
      JB=ITR(B)+J
      ZLBEB=ZLX(BE,B)
      ZLBBE=ZLX(B,BE)
      IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
      IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
      JUBEB=TOFF(J,U,ZLBEB,ais)+TADD(BE,B,ais)
      JUBBE=TOFF(J,U,ZLBBE,ais)+TADD(B,BE,ais)
      T1N(U,BE)=T1N(U,BE)+(DIN(IAJB)+DIN(IAJB)-DIN(IBJA))
     .                   *(TAU(JUBBE)-TAU(JUBEB))
 3165 CONTINUE
 3170 CONTINUE
 3175 CONTINUE
 3180 CONTINUE
c  
c ... END CCSDT D TERMS
c 
c >>> ADD CCSDT TERMS TO T2 EQUATION. F INTS.
C
         DO 6282 XSYM=1,NIRRED
            XBj=0
            DO 6220 B=1,NV
               IF(FZV(B).EQ.1)GO TO 6220
               BSYM=ORBSYM(B+NO)
               DO 6215 j=1,NO
                  IF(FZO(j).EQ.1)GO TO 6215
                  jSYM=ORBSYM(j)
                  BjSYM=IEOR(BSYM,jSYM)
                  IF(BjSYM.NE.XSYM-1)GO TO 6215
                  AiBjS=IEOR(AiSYM,XSYM-1)
                  XBj=XBj+1
                  XGAV=0
                  DO 6210 GA=1,NV
                     IF(FZV(GA).EQ.1)GO TO 6210
                     GASYM=ORBSYM(GA+NO)
                     VSYM=IEOR(GASYM,AiBjS)+1
                     FV=FLOV(VSYM,1)
                     LV=FLOV(VSYM,2)
                     ZLGAB=ZLX(GA,B)
                     DO 6205 V=FV,LV
                        XGAV=XGAV+1
                        VjGAB=TOFF(V,j,ZLGAB,AiS)+TADD(GA,B,AiS)
                        VjBGA=TOFF(j,V,ZLGAB,AiS)+TADD(GA,B,AiS)
                        AR1(XBj,XGAV)=TAU(VjGAB)
                        AR2(XBj,XGAV)=TAU(VjBGA)
 6205                CONTINUE
 6210             CONTINUE
 6215          CONTINUE
 6220       CONTINUE
            LIMBj=XBj
            LIMGAV=XGAV
            BESYM=IEOR(ASYM,XSYM-1)
            BEiSYM=IEOR(BESYM,iSYM)
            FBE=FLOV(BESYM+1,3)-NO
            LBE=FLOV(BESYM+1,4)-NO
            DO 6270 BE=FBE,LBE
               BEA=ITV(MAX0(BE,A))+MIN0(BE,A)
               XBj=0
               DO 6230 B=1,NV
                  IF(FZV(B).EQ.1)GO TO 6230
                  BEB=ITV(MAX0(BE,B))+MIN0(BE,B)
                  BSYM=ORBSYM(B+NO)
                  DO 6225 j=1,NO
                     IF(FZO(j).EQ.1)GO TO 6225
                     jSYM=ORBSYM(j)
                     BjSYM=IEOR(BSYM,jSYM)
                     IF(BjSYM.NE.XSYM-1)GO TO 6225
                     XBj=XBj+1
                     jA=ITR(A)+j
                     jB=ITR(B)+j
                     BEBjA=FOFF(BEB)+FADD(jA)
                     BEAjB=FOFF(BEA)+FADD(jB)
                     SUX1(XBj)= FIN(BEAjB)+FIN(BEAjB)-FIN(BEBjA)
                     SUX2(XBj)=-FIN(BEAjB)
 6225             CONTINUE
 6230          CONTINUE
               DO 6240 XGAV=1,LIMGAV
                  VAL1=0.0D0
                  VAL2=0.0D0
                  DO 6235 XBj=1,LIMBj
                     VAL1=VAL1+SUX1(XBj)*AR1(XBj,XGAV)
                     VAL2=VAL2+SUX2(XBj)*AR2(XBj,XGAV)
 6235             CONTINUE
                  SUX11(XGAV)=VAL1
                  SUX22(XGAV)=VAL2
 6240          CONTINUE
               XGAV=0
               DO 6250 GA=1,BE-1
                  IF(FZV(GA).EQ.1)GO TO 6250
                  GASYM=ORBSYM(GA+NO)
                  VSYM=IEOR(GASYM,BEiSYM)+1
                  FV=FLOV(VSYM,1)
                  LV=FLOV(VSYM,2)
                  DO 6245 V=FV,LV
                     XGAV=XGAV+1
                     iVBEGA=UOFF(i,V,1)+VADD(BE,GA)
                     S(iVBEGA)=S(iVBEGA)+SUX11(XGAV)+SUX22(XGAV)
 6245             CONTINUE
 6250          CONTINUE
               DO 6251 GA=BE,BE
                  IF(FZV(GA).EQ.1)GO TO 6251
                  GASYM=ORBSYM(GA+NO)
                  VSYM=IEOR(GASYM,BEiSYM)+1
                  FV=FLOV(VSYM,1)
                  LV=FLOV(VSYM,2)
                  DO 6246 V=FV,LV
                     XGAV=XGAV+1
                     iVBEGA=UOFF(i,V,1)+VADD(BE,GA)
                     S(iVBEGA)=S(iVBEGA)+SUX11(XGAV)+SUX22(XGAV)
                     ViGABE=UOFF(V,i,1)+VADD(GA,BE)
                     S(ViGABE)=S(ViGABE)+SUX11(XGAV)+SUX22(XGAV)
 6246             CONTINUE
 6251          CONTINUE
               DO 6252 GA=BE+1,NV
                  IF(FZV(GA).EQ.1)GO TO 6252
                  GASYM=ORBSYM(GA+NO)
                  VSYM=IEOR(GASYM,BEiSYM)+1
                  FV=FLOV(VSYM,1)
                  LV=FLOV(VSYM,2)
                  DO 6247 V=FV,LV
                     XGAV=XGAV+1
                     ViGABE=UOFF(V,i,1)+VADD(GA,BE)
                     S(ViGABE)=S(ViGABE)+SUX11(XGAV)+SUX22(XGAV)
 6247             CONTINUE
 6252          CONTINUE
 6270       CONTINUE
 6282    CONTINUE
c   
c ... END CCSDT F TERMS
c
c
C >>> ADD CCSDT TERMS TO T2 EQUATION. E INTS.
C
         DO 6482 XSYM=1,NIRRED
            XbJ=0
            DO 6420 b=1,NV
               IF(FZV(b).EQ.1)GO TO 6420
               bSYM=ORBSYM(b+NO)
               DO 6415 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 6415
                  JSYM=ORBSYM(J)
                  bJSYM=IEOR(bSYM,JSYM)
                  IF(bJSYM.NE.XSYM-1)GO TO 6415
                  aIbJS=IEOR(aISYM,XSYM-1)
                  XbJ=XbJ+1
                  XGAV=0
                  DO 6410 GA=1,NV
                     IF(FZV(GA).EQ.1)GO TO 6410
                     GASYM=ORBSYM(GA+NO)
                     VSYM=IEOR(GASYM,aIbJS)+1
                     FV=FLOV(VSYM,1)
                     LV=FLOV(VSYM,2)
                     ZLGAb=ZLX(GA,b)
                     DO 6405 V=FV,LV
                        XGAV=XGAV+1
                        VJGAb=TOFF(V,J,ZLGAb,aIS)+TADD(GA,b,aIS)
                        VJbGA=TOFF(J,V,ZLGAb,aIS)+TADD(GA,b,aIS)
                        AR1(XbJ,XGAV)=TAU(VJGAb)
                        AR2(XbJ,XGAV)=TAU(VJbGA)
 6405                CONTINUE
 6410             CONTINUE
 6415          CONTINUE
 6420       CONTINUE
            LIMbJ=XbJ
            LIMGAV=XGAV
            USYM=IEOR(ISYM,XSYM-1)
            aUSYM=IEOR(aSYM,USYM)
            FU=FLOV(USYM+1,1)
            LU=FLOV(USYM+1,2)
            DO 6470 U=FU,LU
               UI=ITR(MAX0(U,I))+MIN0(U,I)
               XbJ=0
               DO 6430 b=1,NV
                  IF(FZV(b).EQ.1)GO TO 6430
                  bI=ITR(b)+I
                  bSYM=ORBSYM(b+NO)
                  DO 6425 J=1,NO
                     IF(FZO(J).EQ.1)GO TO 6425
                     JSYM=ORBSYM(J)
                     bJSYM=IEOR(bSYM,JSYM)
                     IF(bJSYM.NE.XSYM-1)GO TO 6425
                     XbJ=XbJ+1
                     UJ=ITR(MAX0(U,J))+MIN0(U,J)
                     bJ=ITR(b)+J
                     UJbI=EOFF(bI)+EADD(UJ)
                     UIbJ=EOFF(bJ)+EADD(UI)
                     SUX1(XbJ)=-EIN(UIbJ)-EIN(UIbJ)+EIN(UJbI)
                     SUX2(XbJ)= EIN(UIbJ)
 6425             CONTINUE
 6430          CONTINUE
               DO 6440 XGAV=1,LIMGAV
                  VAL1=0.0D0
                  VAL2=0.0D0
                  DO 6435 XbJ=1,LIMbJ
                     VAL1=VAL1+SUX1(XbJ)*AR1(XbJ,XGAV)
                     VAL2=VAL2+SUX2(XbJ)*AR2(XbJ,XGAV)
 6435             CONTINUE
                  SUX11(XGAV)=VAL1
                  SUX22(XGAV)=VAL2
 6440          CONTINUE
               XGAV=0
               DO 6450 GA=1,a-1
                  IF(FZV(GA).EQ.1)GO TO 6450
                  GASYM=ORBSYM(GA+NO)
                  VSYM=IEOR(GASYM,aUSYM)+1
                  FV=FLOV(VSYM,1)
                  LV=FLOV(VSYM,2)
                  DO 6445 V=FV,LV
                     XGAV=XGAV+1
                     UVaGA=UOFF(U,V,1)+VADD(a,GA)
                     S(UVaGA)=S(UVaGA)+SUX11(XGAV)+SUX22(XGAV)
 6445             CONTINUE
 6450          CONTINUE
               DO 6451 GA=a,a
                  IF(FZV(GA).EQ.1)GO TO 6451
                  GASYM=ORBSYM(GA+NO)
                  VSYM=IEOR(GASYM,aUSYM)+1
                  FV=FLOV(VSYM,1)
                  LV=FLOV(VSYM,2)
                  DO 6446 V=FV,LV
                     XGAV=XGAV+1
                     UVaGA=UOFF(U,V,1)+VADD(a,GA)
                     S(UVaGA)=S(UVaGA)+SUX11(XGAV)+SUX22(XGAV)
                     VUGAa=UOFF(V,U,1)+VADD(GA,a)
                     S(VUGAa)=S(VUGAa)+SUX11(XGAV)+SUX22(XGAV)
 6446             CONTINUE
 6451          CONTINUE
               DO 6452 GA=a+1,NV
                  IF(FZV(GA).EQ.1)GO TO 6452
                  GASYM=ORBSYM(GA+NO)
                  VSYM=IEOR(GASYM,aUSYM)+1
                  FV=FLOV(VSYM,1)
                  LV=FLOV(VSYM,2)
                  DO 6447 V=FV,LV
                     XGAV=XGAV+1
                     VUGAa=UOFF(V,U,1)+VADD(GA,a)
                     S(VUGAa)=S(VUGAa)+SUX11(XGAV)+SUX22(XGAV)
 6447             CONTINUE
 6452          CONTINUE
 6470       CONTINUE
 6482    CONTINUE
c 
c ... END CCSDT E TERMS
c
 9985    CONTINUE
 9990 CONTINUE
c
      call rclose(itap81,3)
c
C     IF(SORT.EQ.'YES ')
C    .CALL SORTT(TAU,T2,NTAU,ORBSYM,NO,NV,FLOV,NIRRED,NT3,TOFF,TADD,
C    .           ITAP81,ITAP97,ZLX)
C     CALL FRTAU(TAU,T2,T1 ,UOFF,VADD,ORBSYM,FLOV,NO,NV,
C    .           NDIMT2,NIRRED,FZO,FZV,OPTION)
C
          RETURN
          END
