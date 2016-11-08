      subroutine yfterm(a,asym,i,isym,ai,aisym,ais,dim,tau,fin,
     .                  fzo,fzv,orbsym,flov,nirred,no,nv,
     .                  itr,itv,ipq,foff,fadd,toff,tadd,zlx,uoff,
     .                  vadd,nt3,z2o,aux1,aux3,xar1,xar2,xaro,xarv,
     .                  maxv)
      IMPLICIT INTEGER(A-Z)
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
      REAL*8 z2o(*),TAU(*),VAL,VAL1,VAL2,FIN(*),
     .       AUX1(NV,NO),AUX3(NV,NV),
     .       XAR1(MAXV,MAXOV),XAR2(MAXV,MAXOV),xaro(maxv,maxoo),
     .       xarv(maxv,maxvv) 
      INTEGER ORBSYM(*),FLOV(NIRRED,4),TOFF(NO,NO,2,NIRRED),ZLX(NV,NV),
     .        TADD(NV,NV,NIRRED),FZO(NO),FZV(NV),UOFF(NO,NO,2),NT3(8),
     .        VADD(NV,NV),FOFF(*),FADD(*),ITR(*),ITV(*),IPQ(*)
c
            DO 3910 B=1,NV
               IF(FZV(B).EQ.1)GO TO 3910
               BSYM=ORBSYM(B+NO)
               ESYM=IEOR(BSYM,AISYM)+1
               FE=FLOV(ESYM,3)-NO
               LE=FLOV(ESYM,4)-NO
               bI=ITR(b)+I
               DO 3905 E=FE,LE
                  BE=ITV(MAX0(B,E))+MIN0(B,E)
                  aE=ITV(MAX0(a,E))+MIN0(a,E)
                  AIBE=FOFF(BE)+FADD(AI)
                  bIaE=FOFF(aE)+FADD(bI)
                  AUX3(E,B)=-FIN(AIBE)-fin(aibe)+fin(biae)
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
                        AUX1(E,J)=z2o(KJCE)
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
                  AUX1(E,K)=z2o(IKAE)
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
                     cJ=ITR(c)+J
                     DO 4115 E=FE,LE
                        CE=ITV(MAX0(C,E))+MIN0(C,E)
                        bE=ITV(MAX0(b,E))+MIN0(b,E)
                        BJCE=FOFF(CE)+FADD(BJ)
                        cJbE=FOFF(bE)+FADD(cJ)
                        AUX3(E,C)=-FIN(BJCE)-fin(bjce)+fin(cjbe)
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
                  aK=ITR(a)+K
                  XCK=XCK+1
                  ESYM=IEOR(CKSYM,ASYM)+1
                  FE=FLOV(ESYM,3)-NO
                  LE=FLOV(ESYM,4)-NO
                  DO 4202 E=FE,LE
                     AE=ITV(MAX0(A,E))+MIN0(A,E)
                     cE=ITV(MAX0(c,E))+MIN0(c,E)
                     CKAE=FOFF(AE)+FADD(CK)
                     aKcE=FOFF(cE)+FADD(aK)
                     XAR1(E,XCK)=-FIN(CKAE)-fin(ckae)+fin(akce)
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
                     XAR2(E,XBJ)=z2o(JIBE)
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
C
            DO 4392 XSYM=1,NIRRED
            cJSEM=IEOR(AISYM,XSYM-1)
            XbK=0
            DO 4304 b=1,NV
               IF(FZV(b).EQ.1)GO TO 4304
               bSYM=ORBSYM(b+NO)
               DO 4303 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 4303
                  KSYM=ORBSYM(K)
                  bKSYM=IEOR(bSYM,KSYM)
                  IF(bKSYM.NE.XSYM-1)GO TO 4303
                  bK=ITR(b)+K
                  XbK=XbK+1
                  ESYM=IEOR(bKSYM,ASYM)+1
                  FE=FLOV(ESYM,3)-NO
                  LE=FLOV(ESYM,4)-NO
                  DO 4302 E=FE,LE
                     AE=ITV(MAX0(A,E))+MIN0(A,E)
                     bKAE=FOFF(AE)+FADD(bK)
                     XAR1(E,XbK)=FIN(bKAE)
 4302             CONTINUE
 4303          CONTINUE
 4304       CONTINUE
            XcJ=0
            DO 4308 c=1,NV
               IF(FZV(c).EQ.1)GO TO 4308
               cSYM=ORBSYM(c+NO)
               DO 4307 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 4307
                  JSYM=ORBSYM(J)
                  cJSYM=IEOR(cSYM,JSYM)
                  IF(cJSYM.NE.cJSEM)GO TO 4307
                  XcJ=XcJ+1
                  ESYM=IEOR(cJSYM,ISYM)+1
                  FE=FLOV(ESYM,3)-NO
                  LE=FLOV(ESYM,4)-NO
                  DO 4306 E=FE,LE
                     ZLcE=ZLX(c,E)
                     JIcE=UOFF(J,I,ZLcE)+VADD(c,E)
                     XAR2(E,XcJ)=z2o(JIcE)
 4306             CONTINUE
 4307          CONTINUE
 4308       CONTINUE
            XcJ=0
            DO 4380 c=1,NV
               IF(FZV(c).EQ.1)GO TO 4380
               cSYM=ORBSYM(c+NO)
               DO 4375 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 4375
                  JSYM=ORBSYM(J)
                  cJSYM=IEOR(cSYM,JSYM)
                  IF(cJSYM.NE.cJSEM)GO TO 4375
                  XcJ=XcJ+1
                  XbK=0
                  DO 4335 b=1,c-1
                     IF(FZV(b).EQ.1)GO TO 4335
                     bSYM=ORBSYM(b+NO)
                     KSYM=IEOR(bSYM,XSYM-1)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     ESYM=IEOR(ISYM,cJSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4330 K=FK,LK
                     VAL=0.0D0
                     XbK=XbK+1
                        DO 4325 E=FE,LE
                        VAL=VAL+XAR1(E,XbK)*XAR2(E,XcJ)
 4325                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)+VAL
 4330                CONTINUE
 4335             CONTINUE
                  DO 4339 b=c,c
                     IF(FZV(b).EQ.1)GO TO 4339
                     bSYM=ORBSYM(b+NO)
                     KSYM=IEOR(bSYM,XSYM-1)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     ESYM=IEOR(ISYM,cJSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4338 K=FK,LK
                     XbK=XbK+1
                     VAL=0.0D0
                        DO 4337 E=FE,LE
                        VAL=VAL+XAR1(E,XbK)*XAR2(E,XcJ)
 4337                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)+VAL
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)+VAL
 4338                CONTINUE
 4339             CONTINUE
                  DO 4350 b=c+1,NV
                     IF(FZV(b).EQ.1)GO TO 4350
                     bSYM=ORBSYM(b+NO)
                     KSYM=IEOR(bSYM,XSYM-1)+1
                     FK=FLOV(KSYM,1)
                     LK=FLOV(KSYM,2)
                     ESYM=IEOR(ISYM,cJSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4345 K=FK,LK
                     XbK=XbK+1
                     VAL=0.0D0
                        DO 4340 E=FE,LE
                        VAL=VAL+XAR1(E,XbK)*XAR2(E,XcJ)
 4340                   CONTINUE
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)+VAL
 4345                CONTINUE
 4350             CONTINUE
 4375          CONTINUE
 4380       CONTINUE
 4392       CONTINUE
C
            DO 4492 XSYM=1,NIRRED
            BkSEM=IEOR(AISYM,XSYM-1)
            XCj=0
            DO 4404 C=1,NV
               IF(FZV(C).EQ.1)GO TO 4404
               CSYM=ORBSYM(C+NO)
               DO 4403 j=1,NO
                  IF(FZO(j).EQ.1)GO TO 4403
                  jSYM=ORBSYM(j)
                  CjSYM=IEOR(CSYM,jSYM)
                  IF(CjSYM.NE.XSYM-1)GO TO 4403
                  aj=ITR(a)+j
                  XCj=XCj+1
                  ESYM=IEOR(CjSYM,ASYM)+1
                  FE=FLOV(ESYM,3)-NO
                  LE=FLOV(ESYM,4)-NO
                  DO 4402 E=FE,LE
                     cE=ITV(MAX0(c,E))+MIN0(c,E)
                     ajcE=FOFF(cE)+FADD(aj)
                     XAR1(E,XCj)=fin(ajce)
 4402             CONTINUE
 4403          CONTINUE
 4404       CONTINUE
            XBk=0
            DO 4408 B=1,NV
               IF(FZV(B).EQ.1)GO TO 4408
               BSYM=ORBSYM(B+NO)
               DO 4407 k=1,NO
                  IF(FZO(k).EQ.1)GO TO 4407
                  kSYM=ORBSYM(k)
                  BkSYM=IEOR(BSYM,kSYM)
                  IF(BkSYM.NE.BkSEM)GO TO 4407
                  XBk=XBk+1
                  ESYM=IEOR(BkSYM,ISYM)+1
                  FE=FLOV(ESYM,3)-NO
                  LE=FLOV(ESYM,4)-NO
                  DO 4406 E=FE,LE
                     ZLBE=ZLX(B,E)
                     ikBE=UOFF(i,k,ZLBE)+VADD(B,E)
                     XAR2(E,XBk)=z2o(ikBE)
 4406             CONTINUE
 4407          CONTINUE
 4408       CONTINUE
            XBk=0
            DO 4480 B=1,NV
               IF(FZV(B).EQ.1)GO TO 4480
               BSYM=ORBSYM(B+NO)
               DO 4475 k=1,NO
                  IF(FZO(k).EQ.1)GO TO 4475
                  kSYM=ORBSYM(k)
                  BkSYM=IEOR(BSYM,kSYM)
                  IF(BkSYM.NE.BkSEM)GO TO 4475
                  XBk=XBk+1
                  XCj=0
                  DO 4435 C=1,B-1
                     IF(FZV(C).EQ.1)GO TO 4435
                     CSYM=ORBSYM(C+NO)
                     jSYM=IEOR(CSYM,XSYM-1)+1
                     Fj=FLOV(jSYM,1)
                     Lj=FLOV(jSYM,2)
                     ESYM=IEOR(ISYM,BkSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4430 j=Fj,Lj
                     VAL=0.0D0
                     XCj=XCj+1
                        DO 4425 E=FE,LE
                        VAL=VAL+XAR1(E,XCj)*XAR2(E,XBk)
 4425                   CONTINUE
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)+VAL
 4430                CONTINUE
 4435             CONTINUE
                  DO 4439 C=B,B
                     IF(FZV(C).EQ.1)GO TO 4439
                     CSYM=ORBSYM(C+NO)
                     jSYM=IEOR(CSYM,XSYM-1)+1
                     Fj=FLOV(jSYM,1)
                     Lj=FLOV(jSYM,2)
                     ESYM=IEOR(ISYM,bkSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4438 j=Fj,Lj
                     XCj=XCj+1
                     VAL=0.0D0
                        DO 4437 E=FE,LE
                        VAL=VAL+XAR1(E,XCj)*XAR2(E,XBk)
 4437                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)+VAL
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)+VAL
 4438                CONTINUE
 4439             CONTINUE
                  DO 4450 C=B+1,NV
                     IF(FZV(C).EQ.1)GO TO 4450
                     CSYM=ORBSYM(C+NO)
                     jSYM=IEOR(CSYM,XSYM-1)+1
                     Fj=FLOV(jSYM,1)
                     Lj=FLOV(jSYM,2)
                     ESYM=IEOR(ISYM,bkSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4445 j=Fj,Lj
                     XCj=XCj+1
                     VAL=0.0D0
                        DO 4440 E=FE,LE
                        VAL=VAL+XAR1(E,XCj)*XAR2(E,XBk)
 4440                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)+VAL
 4445                CONTINUE
 4450             CONTINUE
 4475          CONTINUE
 4480       CONTINUE
 4492       CONTINUE
C
            DO 4592 XSYM=1,NIRRED
            jkSEM=IEOR(AISYM,XSYM-1)
            Xcb=0
            DO 4504 C=1,NV
               IF(FZV(C).EQ.1)GO TO 4504
               CSYM=ORBSYM(C+NO)
               ci=itr(c)+i
               DO 4503 b=1,nv
                  IF(FZv(b).EQ.1)GO TO 4503
                  bSYM=ORBSYM(b+no)
                  CbSYM=IEOR(CSYM,bSYM)
                  IF(CbSYM.NE.XSYM-1)GO TO 4503
                  XCb=XCb+1
                  ESYM=IEOR(CbSYM,iSYM)+1
                  FE=FLOV(ESYM,3)-NO
                  LE=FLOV(ESYM,4)-NO
                  DO 4502 E=FE,LE
                     bE=ITV(MAX0(b,E))+MIN0(b,E)
                     cibe=FOFF(bE)+FADD(ci)
                     XARv(E,XCb)=fin(cibe)
 4502             CONTINUE
 4503          CONTINUE
 4504       CONTINUE
            Xjk=0
            DO 4508 j=1,No
               IF(FZo(j).EQ.1)GO TO 4508
               jSYM=ORBSYM(j)
               DO 4507 k=1,NO
                  IF(FZO(k).EQ.1)GO TO 4507
                  kSYM=ORBSYM(k)
                  jkSYM=IEOR(jSYM,kSYM)
                  IF(jkSYM.NE.jkSEM)GO TO 4507
                  Xjk=Xjk+1
                  ESYM=IEOR(jkSYM,aSYM)+1
                  FE=FLOV(ESYM,3)-NO
                  LE=FLOV(ESYM,4)-NO
                  DO 4506 E=FE,LE
                     ZLaE=ZLX(a,E)
                     kjaE=UOFF(k,j,ZLaE)+VADD(a,E)
                     XARo(E,Xjk)=z2o(kjaE)
 4506             CONTINUE
 4507          CONTINUE
 4508       CONTINUE
            Xjk=0
            DO 4580 j=1,No
               IF(FZo(j).EQ.1)GO TO 4580
               jSYM=ORBSYM(j)
               DO 4575 k=1,NO
                  IF(FZO(k).EQ.1)GO TO 4575
                  kSYM=ORBSYM(k)
                  jkSYM=IEOR(jSYM,kSYM)
                  IF(jkSYM.NE.jkSEM)GO TO 4575
                  Xjk=Xjk+1
                  XCb=0
                  DO 4535 C=1,nv
                     IF(FZV(C).EQ.1)GO TO 4535
                     CSYM=ORBSYM(C+NO)
                     bSYM=IEOR(CSYM,XSYM-1)+1
                     Fb=FLOV(bSYM,3)-no
                     Lb=FLOV(bSYM,4)-no
                     ESYM=IEOR(aSYM,jkSYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 4530 b=fb,lb
                        VAL=0.0D0
                        XCb=XCb+1
                        zlbc=zlx(b,c)
                        DO 4525 E=FE,LE
                           VAL=VAL+XARv(E,XCb)*XARo(E,Xjk)
 4525                   CONTINUE
                        if(b.ne.c)then
                           JKBC=TOFF(J,K,zlbc,AIS)+TADD(B,C,AIS)
                           TAU(JKBC)=TAU(JKBC)+VAL
                        else
                           KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                           TAU(KJCB)=TAU(KJCB)+VAL
                           JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                           TAU(JKBC)=TAU(JKBC)+VAL
                        endif
 4530                CONTINUE
 4535             CONTINUE
 4575       CONTINUE
 4580       CONTINUE
 4592       CONTINUE
C
      return
      end
