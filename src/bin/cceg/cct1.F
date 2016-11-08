ccsd(t) routine based on sdt1. will probably not work ...
      SUBROUTINE CCT1 (S,T1N,T2 ,T1,TAU,NO,NV,NDIMT2,NIT,
     .                 ORBSYM,FLOV,NIRRED,TOFF,TADD,ZLX,FZO,FZV,
     .                 UOFF,VADD,OPTION,RSTR,NT3,DOFF,DADD,EOFF,EADD,
     .                 FOFF,FADD,ITV,ITR,IPQ,DIN,EIN,FIN,fl81,FOCK,
     .                 AUX1,AUX2,AUX3,SUX1,SUX2,SUX3,
     .                 XAR1,XAR2,NOACT,NVACT,MAXO,MAXV,MAXOV,
     .                 AR1,AR2,SUX11,SUX22,XAR3,XAR4,AUX4,AUX5)
      IMPLICIT INTEGER(A-Z)
      REAL*8  S(NDIMT2),T1N(NO,NV),T2(NDIMT2),T1(NO,NV),TAU(*),FOCK(*),
     .       VAL,VAL1,VAL2,DIN(*),EIN(*),FIN(*),val4,
     .       AUX1(NV,NO),AUX2(NV,NO),AUX3(NV,NV),SUX11(NO*NV),
     .       AUX4(NO,NO),AUX5(NO,NV),
     .       SUX1(NV*NO),SUX2(NV*NO),SUX3(NV*NV),SUX22(NO*NV),
     .       XAR1(MAXV,MAXOV),XAR2(MAXV,MAXOV),
     .       XAR3(MAXO,MAXOV),XAR4(MAXO,MAXOV),
     .       AR1(MAXOV,MAXOV),AR2(MAXOV,MAXOV)
      real*8 w(5,5,5),r(5,5,5)
      CHARACTER*4 OPTION,RSTR,SORT,fl81
      INTEGER ORBSYM(*),FLOV(NIRRED,4),TOFF(NO,NO,2,NIRRED),ZLX(NV,NV),
     .        TADD(NV,NV,NIRRED),FZO(NO),FZV(NV),UOFF(NO,NO,2),NT3(8),
     .        VADD(NV,NV),DOFF(*),DADD(*),FOFF(*),FADD(*),EOFF(*),
     .        EADD(*),ITR(*),ITV(*),IPQ(*)
C
      SORT='NO  '
c
      call zero(w,no*no*no)
      call zero(r,no*no*no)
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
c           DO 3980 C=1,NV
            DO 3980 C=1,a 
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
                     val2=val+val
                     val4=val2+val2
                     w(i,k,j)=w(i,k,j)+val
                     r(i,k,j)=r(i,k,j)+val4
                     r(i,j,k)=r(i,j,k)-val2
                     r(k,i,j)=r(k,i,j)-val2
                     r(j,i,k)=r(j,i,k)+val
                     r(j,k,i)=r(j,k,i)-val2
                     r(k,j,i)=r(k,j,i)+val
c                    KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
c                    TAU(KJCB)=TAU(KJCB)+VAL
 3930                CONTINUE
 3935             CONTINUE
c                 DO 3950 B=C,NV
c                    IF(FZV(B).EQ.1)GO TO 3950
c                    BSYM=ORBSYM(B+NO)
c                    JSYM=IEOR(BSYM,AICKS)+1
c                    FJ=FLOV(JSYM,1)
c                    LJ=FLOV(JSYM,2)
c                    ESYM=IEOR(BSYM,AISYM)+1
c                    FE=FLOV(ESYM,3)-NO
c                    LE=FLOV(ESYM,4)-NO
c                    DO 3945 J=FJ,LJ
c                    VAL=0.0D0
c                       DO 3940 E=FE,LE
c                       VAL=VAL+AUX3(E,B)*AUX1(E,J)
c3940                   CONTINUE
c                    JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
c                    TAU(JKBC)=TAU(JKBC)+VAL
c3945                CONTINUE
c3950             CONTINUE
 3975          CONTINUE
 3980       CONTINUE
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
c           DO 4180 B=1,NV
            DO 4180 B=1,a 
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
                     val2=val+val
                     val4=val2+val2
                     w(i,j,k)=w(i,j,k)+val
                     r(i,j,k)=r(i,j,k)+val4
                     r(i,k,j)=r(i,k,j)-val2
                     r(j,i,k)=r(j,i,k)-val2
                     r(k,i,j)=r(k,i,j)+val
                     r(k,j,i)=r(k,j,i)-val2
                     r(j,k,i)=r(j,k,i)+val
c                    JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
c                    TAU(JKBC)=TAU(JKBC)+VAL
 4130                CONTINUE
 4135             CONTINUE
c                 DO 4150 C=B,NV
c                    IF(FZV(C).EQ.1)GO TO 4150
c                    CSYM=ORBSYM(C+NO)
c                    KSYM=IEOR(CSYM,AIBJS)+1
c                    FK=FLOV(KSYM,1)
c                    LK=FLOV(KSYM,2)
c                    ESYM=IEOR(CSYM,BJSYM)+1
c                    FE=FLOV(ESYM,3)-NO
c                    LE=FLOV(ESYM,4)-NO
c                    DO 4145 K=FK,LK
c                    VAL=0.0D0
c                       DO 4140 E=FE,LE
c                       VAL=VAL+AUX3(E,C)*AUX1(E,K)
c4140                   CONTINUE
c                    KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
c                    TAU(KJCB)=TAU(KJCB)+VAL
c4145                CONTINUE
c4150             CONTINUE
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
c                 DO 4235 C=1,B-1
                  DO 4235 C=1,B
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
c                    JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
c                    TAU(JKBC)=TAU(JKBC)+VAL
                     val2=val+val
                     val4=val2+val2
                     w(i,j,k)=w(i,j,k)+val
                     r(i,j,k)=r(i,j,k)+val4
                     r(i,k,j)=r(i,k,j)-val2
                     r(j,i,k)=r(j,i,k)-val2
                     r(k,i,j)=r(k,i,j)+val
                     r(k,j,i)=r(k,j,i)-val2
                     r(j,k,i)=r(j,k,i)+val
 4230                CONTINUE
 4235             CONTINUE
c                 DO 4239 C=B,B
c                    IF(FZV(C).EQ.1)GO TO 4239
c                    CSYM=ORBSYM(C+NO)
c                    KSYM=IEOR(CSYM,XSYM-1)+1
c                    FK=FLOV(KSYM,1)
c                    LK=FLOV(KSYM,2)
c                    ESYM=IEOR(ISYM,BJSYM)+1
c                    FE=FLOV(ESYM,3)-NO
c                    LE=FLOV(ESYM,4)-NO
c                    DO 4238 K=FK,LK
c                    XCK=XCK+1
c                    VAL=0.0D0
c                       DO 4237 E=FE,LE
c                       VAL=VAL+XAR1(E,XCK)*XAR2(E,XBJ)
c4237                   CONTINUE
c                    KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
c                    TAU(KJCB)=TAU(KJCB)+VAL
c                    JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
c                    TAU(JKBC)=TAU(JKBC)+VAL
c4238                CONTINUE
c4239             CONTINUE
c                 DO 4250 C=B+1,NV
c                    IF(FZV(C).EQ.1)GO TO 4250
c                    CSYM=ORBSYM(C+NO)
c                    KSYM=IEOR(CSYM,XSYM-1)+1
c                    FK=FLOV(KSYM,1)
c                    LK=FLOV(KSYM,2)
c                    ESYM=IEOR(ISYM,BJSYM)+1
c                    FE=FLOV(ESYM,3)-NO
c                    LE=FLOV(ESYM,4)-NO
c                    DO 4245 K=FK,LK
c                    XCK=XCK+1
c                    VAL=0.0D0
c                       DO 4240 E=FE,LE
c                       VAL=VAL+XAR1(E,XCK)*XAR2(E,XBJ)
c4240                   CONTINUE
c                    KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
c                    TAU(KJCB)=TAU(KJCB)+VAL
c4245                CONTINUE
c4250             CONTINUE
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
c           DO 4980 C=1,NV
            DO 4980 C=1,a 
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
c                    FB2=FB
c                    LB2=LB
c                    MSYM=IEOR(JSYM,AISYM)+1
c                    FM=FLOV(MSYM,1)
c                    LM=FLOV(MSYM,2)
c                    IF(LB.GT.C)LB2=C
c                    DO 4930 B=FB2,LB2
c                    VAL=0.0D0
c                       DO 4925 M=FM,LM
c                       VAL=VAL+AUX4(M,J)*AUX5(M,B)
c4925                   CONTINUE
c                    KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
c                    TAU(KJCB)=TAU(KJCB)-VAL
c4930                CONTINUE
                     FB2=FB
                     LB2=LB
                     IF(FB.LT.C)FB2=C
                     DO 4931 B=FB2,LB2
                     VAL=0.0D0
                        DO 4926 M=FM,LM
                        VAL=VAL+AUX4(M,J)*AUX5(M,B)
 4926                   CONTINUE
c                    JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
c                    TAU(JKBC)=TAU(JKBC)-VAL
                     val2=val+val
                     val4=val2+val2
                     w(i,k,j)=w(i,k,j)+val
                     r(i,k,j)=r(i,k,j)+val4
                     r(i,j,k)=r(i,j,k)-val2
                     r(k,i,j)=r(k,i,j)-val2
                     r(j,i,k)=r(j,i,k)+val
                     r(j,k,i)=r(j,k,i)-val2
                     r(k,j,i)=r(k,j,i)+val
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
c           DO 5180 B=1,NV
            DO 5180 B=1,a   
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
c                    JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
c                    TAU(JKBC)=TAU(JKBC)-VAL
                     val2=val+val
                     val4=val2+val2
                     w(i,j,k)=w(i,j,k)+val
                     r(i,j,k)=r(i,j,k)+val4
                     r(i,k,j)=r(i,k,j)-val2
                     r(j,i,k)=r(j,i,k)-val2
                     r(k,i,j)=r(k,i,j)+val
                     r(k,j,i)=r(k,j,i)-val2
                     r(j,k,i)=r(j,k,i)+val
 5130                CONTINUE
c                    FC2=FC
c                    LC2=LC
c                    IF(FC.LT.B)FC2=B
c                    DO 5131 C=FC2,LC2
c                    VAL=0.0D0
c                       DO 5126 M=FM,LM
c                       VAL=VAL+AUX4(M,K)*AUX5(M,C)
c5126                   CONTINUE
c                    KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
c                    TAU(KJCB)=TAU(KJCB)-VAL
c5131                CONTINUE
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
c           DO 5280 B=1,NV
            DO 5280 B=1,a 
               IF(FZV(B).EQ.1)GO TO 5280
               BSYM=ORBSYM(B+NO)
               DO 5275 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 5275
                  JSYM=ORBSYM(J)
                  BJSYM=IEOR(BSYM,JSYM)
                  IF(BJSYM.NE.BJSEM)GO TO 5275
                  XBJ=XBJ+1
                  XCK=0
c                 DO 5235 C=1,B-1
                  DO 5235 C=1,B
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
c                    JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
c                    TAU(JKBC)=TAU(JKBC)-VAL
                     val2=val+val
                     val4=val2+val2
                     w(i,j,k)=w(i,j,k)+val
                     r(i,j,k)=r(i,j,k)+val4
                     r(i,k,j)=r(i,k,j)-val2
                     r(j,i,k)=r(j,i,k)-val2
                     r(k,i,j)=r(k,i,j)+val
                     r(k,j,i)=r(k,j,i)-val2
                     r(j,k,i)=r(j,k,i)+val
 5230                CONTINUE
 5235             CONTINUE
c                 DO 5249 C=B,B
c                    IF(FZV(C).EQ.1)GO TO 5249
c                    CSYM=ORBSYM(C+NO)
c                    KSYM=IEOR(CSYM,XSYM-1)+1
c                    FK=FLOV(KSYM,1)
c                    LK=FLOV(KSYM,2)
c                    MSYM=IEOR(ASYM,BJSYM)+1
c                    FM=FLOV(MSYM,1)
c                    LM=FLOV(MSYM,2)
c                    DO 5238 K=FK,LK
c                    XCK=XCK+1
c                    VAL=0.0D0
c                       DO 5237 M=FM,LM
c                       VAL=VAL+XAR3(M,XCK)*XAR4(M,XBJ)
c5237                   CONTINUE
c                    KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
c                    TAU(KJCB)=TAU(KJCB)-VAL
c                    JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
c                    TAU(JKBC)=TAU(JKBC)-VAL
c5238                CONTINUE
c5249             CONTINUE
c                 DO 5250 C=B+1,NV
c                    IF(FZV(C).EQ.1)GO TO 5250
c                    CSYM=ORBSYM(C+NO)
c                    KSYM=IEOR(CSYM,XSYM-1)+1
c                    FK=FLOV(KSYM,1)
c                    LK=FLOV(KSYM,2)
c                    MSYM=IEOR(ASYM,BJSYM)+1
c                    FM=FLOV(MSYM,1)
c                    LM=FLOV(MSYM,2)
c                    DO 5245 K=FK,LK
c                    XCK=XCK+1
c                    VAL=0.0D0
c                       DO 5240 M=FM,LM
c                       VAL=VAL+XAR3(M,XCK)*XAR4(M,XBJ)
c5240                   CONTINUE
c                    KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
c                    TAU(KJCB)=TAU(KJCB)-VAL
c5245                CONTINUE
c5250             CONTINUE
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
c           DO 3560 B=1,NV
            DO 3560 B=1,a 
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
                     XC=IPQ(C+NO)+C+NO
c                    ZLBC=ZLX(B,C)
c                    JKBC=TOFF(J,K,ZLBC,AIS)+TADD(B,C,AIS)
                     VAL=FOCK(XI)+FOCK(XJ)+FOCK(XK)
     .                  -FOCK(XA)-FOCK(XB)-FOCK(XC)
c                    TAU(JKBC)=TAU(JKBC)/VAL
                     r(i,j,k)=r(i,j,k)/val
 3552                CONTINUE
 3550             CONTINUE
 3560          CONTINUE
 3570       CONTINUE
c
c
 9985    CONTINUE
 9990 CONTINUE
c
          RETURN
          END
