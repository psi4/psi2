      subroutine yeterm(a,asym,i,isym,ai,aisym,ais,dim,tau,ein,
     .                  fzo,fzv,orbsym,flov,nirred,no,nv,
     .                  itr,itv,ipq,eoff,eadd,toff,tadd,zlx,uoff,
     .                  vadd,nt3,z2o,aux4,aux5,xar3,xar4,yaro,yarv,
     .                  maxo)
      IMPLICIT INTEGER(A-Z)
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
      REAL*8 z2o(*),TAU(*),VAL,VAL1,VAL2,EIN(*),
     .       AUX4(NO,NO),AUX5(NO,NV),
     .       XAR3(MAXO,MAXOV),XAR4(MAXO,MAXOV),
     .       yaro(maxo,maxoo),yarv(maxo,maxvv)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),TOFF(NO,NO,2,NIRRED),ZLX(NV,NV),
     .        TADD(NV,NV,NIRRED),FZO(NO),FZV(NV),UOFF(NO,NO,2),NT3(8),
     .        VADD(NV,NV),EOFF(*),EADD(*),ITR(*),ITV(*),IPQ(*)
c
            DO 4910 J=1,NO
               IF(FZO(J).EQ.1)GO TO 4910
               aj=itr(a)+j
               JSYM=ORBSYM(J)
               MSYM=IEOR(JSYM,AISYM)+1
               FM=FLOV(MSYM,1)
               LM=FLOV(MSYM,2)
               DO 4905 M=FM,LM
                  JM=ITR(MAX0(J,M))+MIN0(J,M)
                  iM=ITR(MAX0(i,M))+MIN0(i,M)
                  AIJM=EOFF(AI)+EADD(JM)
                  AjiM=EOFF(Aj)+EADD(iM)
                  AUX4(M,J)=EIN(AIJM)+ein(aijm)-ein(ajim)
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
                        AUX5(M,B)=z2o(KMCB)
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
                     TAU(KJCB)=TAU(KJCB)+VAL
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
                     TAU(JKBC)=TAU(JKBC)+VAL
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
                  AUX5(M,C)=z2o(IMAC)
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
                     bk=itr(b)+k
                     KSYM=ORBSYM(K)
                     MSYM=IEOR(KSYM,BJSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5115 M=FM,LM
                        KM=ITR(MAX0(K,M))+MIN0(K,M)
                        jM=ITR(MAX0(j,M))+MIN0(j,M)
                        BJKM=EOFF(BJ)+EADD(KM)
                        BkjM=EOFF(Bk)+EADD(jM)
                        AUX4(M,K)=EIN(BJKM)+ein(bjkm)-ein(bkjm)
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
                     TAU(JKBC)=TAU(JKBC)+VAL
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
                     TAU(KJCB)=TAU(KJCB)+VAL
 5131                CONTINUE
 5135             CONTINUE
 5175          CONTINUE
 5180       CONTINUE
c
c
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
                  Ci=ITR(C)+i
                  XCK=XCK+1
                  MSYM=IEOR(CKSYM,ISYM)+1
                  FM=FLOV(MSYM,1)
                  LM=FLOV(MSYM,2)
                  DO 5202 M=FM,LM
                     IM=ITR(MAX0(I,M))+MIN0(I,M)
                     kM=ITR(MAX0(k,M))+MIN0(k,M)
                     CKIM=EOFF(CK)+EADD(IM)
                     CikM=EOFF(Ci)+EADD(kM)
                     XAR3(M,XCK)=EIN(CKIM)+ein(ckim)-ein(cikm)
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
                     XAR4(M,XBJ)=z2o(JMBA)
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
                     TAU(JKBC)=TAU(JKBC)+VAL
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
                     TAU(KJCB)=TAU(KJCB)+VAL
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)+VAL
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
                     TAU(KJCB)=TAU(KJCB)+VAL
 5245                CONTINUE
 5250             CONTINUE
 5275          CONTINUE
 5280       CONTINUE
 5292       CONTINUE
C
c
            DO 5392 XSYM=1,NIRRED
            cJSEM=IEOR(AISYM,XSYM-1)
            XbK=0
            DO 5304 b=1,NV
               IF(FZV(b).EQ.1)GO TO 5304
               bSYM=ORBSYM(b+NO)
               DO 5303 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 5303
                  KSYM=ORBSYM(K)
                  bKSYM=IEOR(bSYM,KSYM)
                  IF(bKSYM.NE.XSYM-1)GO TO 5303
                  bK=ITR(b)+K
                  XbK=XbK+1
                  MSYM=IEOR(bKSYM,ISYM)+1
                  FM=FLOV(MSYM,1)
                  LM=FLOV(MSYM,2)
                  DO 5302 M=FM,LM
                     IM=ITR(MAX0(I,M))+MIN0(I,M)
                     bKIM=EOFF(bK)+EADD(IM)
                     XAR3(M,XbK)=ein(bkim)
 5302             CONTINUE
 5303          CONTINUE
 5304       CONTINUE
            XcJ=0
            DO 5308 c=1,NV
               IF(FZV(c).EQ.1)GO TO 5308
               ZLcA=ZLX(c,A)
               cSYM=ORBSYM(c+NO)
               DO 5307 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 5307
                  JSYM=ORBSYM(J)
                  cJSYM=IEOR(cSYM,JSYM)
                  IF(cJSYM.NE.cJSEM)GO TO 5307
                  XcJ=XcJ+1
                  MSYM=IEOR(cJSYM,ASYM)+1
                  FM=FLOV(MSYM,1)
                  LM=FLOV(MSYM,2)
                  DO 5306 M=FM,LM
                     JMcA=UOFF(J,M,ZLcA)+VADD(c,A)
                     XAR4(M,XcJ)=z2o(JMcA)
 5306             CONTINUE
 5307          CONTINUE
 5308       CONTINUE
            XBk=0
            DO 5380 B=1,NV
               IF(FZV(B).EQ.1)GO TO 5380
               BSYM=ORBSYM(B+NO)
               DO 5375 k=1,NO
                  IF(FZO(k).EQ.1)GO TO 5375
                  kSYM=ORBSYM(k)
                  BkSYM=IEOR(BSYM,kSYM)
                  IF(BkSYM.NE.xsym-1)GO TO 5375
                  XBk=XBk+1
                  XCj=0
                  DO 5335 C=1,B-1
                     IF(FZV(C).EQ.1)GO TO 5335
                     CSYM=ORBSYM(C+NO)
                     jsym=ieor(csym,cjsem)+1
                     fj=FLOV(jSYM,1)
                     Lj=FLOV(jSYM,2)
                     MSYM=IEOR(iSYM,BkSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5330 j=Fj,Lj
                     VAL=0.0D0
                     XCj=XCj+1
                        DO 5325 M=FM,LM
                        VAL=VAL+XAR4(M,XCj)*XAR3(M,XBk)
 5325                   CONTINUE
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)-VAL
 5330                CONTINUE
 5335             CONTINUE
                  DO 5349 C=B,B
                     IF(FZV(C).EQ.1)GO TO 5349
                     CSYM=ORBSYM(C+NO)
                     jSYM=IEOR(CSYM,cjsem)+1
                     Fj=FLOV(jSYM,1)
                     Lj=FLOV(jSYM,2)
                     MSYM=IEOR(iSYM,BkSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5338 j=Fj,Lj
                     XCj=XCj+1
                     VAL=0.0D0
                        DO 5337 M=FM,LM
                        VAL=VAL+XAR4(M,XCj)*XAR3(M,XBk)
 5337                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)-VAL
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)-VAL
 5338                CONTINUE
 5349             CONTINUE
                  DO 5350 C=B+1,NV
                     IF(FZV(C).EQ.1)GO TO 5350
                     CSYM=ORBSYM(C+NO)
                     jSYM=IEOR(CSYM,cjsem)+1
                     Fj=FLOV(jSYM,1)
                     Lj=FLOV(jSYM,2)
                     MSYM=IEOR(iSYM,BkSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5345 j=Fj,Lj
                     XCj=XCj+1
                     VAL=0.0D0
                        DO 5340 M=FM,LM
                        VAL=VAL+XAR4(M,XCj)*XAR3(M,XBk)
 5340                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)-VAL
 5345                CONTINUE
 5350             CONTINUE
 5375          CONTINUE
 5380       CONTINUE
 5392       CONTINUE
C
c
            DO 5492 XSYM=1,NIRRED
            bkSEM=IEOR(AISYM,XSYM-1)
            Xcj=0
            DO 5404 c=1,NV
               IF(FZV(c).EQ.1)GO TO 5404
               cSYM=ORBSYM(c+NO)
               ci=ITR(c)+i
               DO 5403 j=1,NO
                  IF(FZO(j).EQ.1)GO TO 5403
                  jSYM=ORBSYM(j)
                  cjSYM=IEOR(cSYM,jSYM)
                  IF(cjSYM.NE.XSYM-1)GO TO 5403
                  Xcj=Xcj+1
                  MSYM=IEOR(cjSYM,ISYM)+1
                  FM=FLOV(MSYM,1)
                  LM=FLOV(MSYM,2)
                  DO 5402 M=FM,LM
                     jM=ITR(MAX0(j,M))+MIN0(j,M)
                     cijM=EOFF(ci)+EADD(jM)
                     XAR3(M,Xcj)=ein(cijm)
 5402             CONTINUE
 5403          CONTINUE
 5404       CONTINUE
            Xbk=0
            DO 5408 b=1,NV
               IF(FZV(b).EQ.1)GO TO 5408
               ZLbA=ZLX(b,A)
               bSYM=ORBSYM(b+NO)
               DO 5407 k=1,NO
                  IF(FZO(k).EQ.1)GO TO 5407
                  kSYM=ORBSYM(k)
                  bkSYM=IEOR(bSYM,kSYM)
                  IF(bkSYM.NE.bkSEM)GO TO 5407
                  Xbk=Xbk+1
                  MSYM=IEOR(bkSYM,ASYM)+1
                  FM=FLOV(MSYM,1)
                  LM=FLOV(MSYM,2)
                  DO 5406 M=FM,LM
                     mkbA=UOFF(m,k,ZLbA)+VADD(b,A)
                     XAR4(M,Xbk)=z2o(mkbA)
 5406             CONTINUE
 5407          CONTINUE
 5408       CONTINUE
            XBk=0
            DO 5480 B=1,NV
               IF(FZV(B).EQ.1)GO TO 5480
               BSYM=ORBSYM(B+NO)
               DO 5475 k=1,NO
                  IF(FZO(k).EQ.1)GO TO 5475
                  kSYM=ORBSYM(k)
                  BkSYM=IEOR(BSYM,kSYM)
                  IF(BkSYM.NE.BkSEM)GO TO 5475
                  XBk=XBk+1
                  XCj=0
                  DO 5435 C=1,B-1
                     IF(FZV(C).EQ.1)GO TO 5435
                     CSYM=ORBSYM(C+NO)
                     jSYM=IEOR(CSYM,XSYM-1)+1
                     Fj=FLOV(jSYM,1)
                     Lj=FLOV(jSYM,2)
                     MSYM=IEOR(ASYM,BkSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5430 j=Fj,Lj
                     VAL=0.0D0
                     XCj=XCj+1
                        DO 5425 M=FM,LM
                        VAL=VAL+XAR3(M,XCj)*XAR4(M,XBk)
 5425                   CONTINUE
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)-VAL
 5430                CONTINUE
 5435             CONTINUE
                  DO 5449 C=B,B
                     IF(FZV(C).EQ.1)GO TO 5449
                     CSYM=ORBSYM(C+NO)
                     jSYM=IEOR(CSYM,XSYM-1)+1
                     Fj=FLOV(jSYM,1)
                     Lj=FLOV(jSYM,2)
                     MSYM=IEOR(ASYM,BkSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5438 j=Fj,Lj
                     XCj=XCj+1
                     VAL=0.0D0
                        DO 5437 M=FM,LM
                        VAL=VAL+XAR3(M,XCj)*XAR4(M,XBk)
 5437                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)-VAL
                     JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                     TAU(JKBC)=TAU(JKBC)-VAL
 5438                CONTINUE
 5449             CONTINUE
                  DO 5450 C=B+1,NV
                     IF(FZV(C).EQ.1)GO TO 5450
                     CSYM=ORBSYM(C+NO)
                     jSYM=IEOR(CSYM,XSYM-1)+1
                     Fj=FLOV(jSYM,1)
                     Lj=FLOV(jSYM,2)
                     MSYM=IEOR(ASYM,BkSYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 5445 j=Fj,Lj
                     XCj=XCj+1
                     VAL=0.0D0
                        DO 5440 M=FM,LM
                        VAL=VAL+XAR3(M,XCj)*XAR4(M,XBk)
 5440                   CONTINUE
                     KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                     TAU(KJCB)=TAU(KJCB)-VAL
 5445                CONTINUE
 5450             CONTINUE
 5475          CONTINUE
 5480       CONTINUE
 5492       CONTINUE
C
            DO 5592 XSYM=1,NIRRED
            bcSEM=IEOR(AISYM,XSYM-1)
            Xkj=0
            DO 5504 k=1,No
               IF(FZo(k).EQ.1)GO TO 5504
               kSYM=ORBSYM(k)
               DO 5503 j=1,no
                  IF(FZo(j).EQ.1)GO TO 5503
                  ja=itr(a)+j
                  jSYM=ORBSYM(j)
                  kjSYM=IEOR(kSYM,jSYM)
                  IF(kjSYM.NE.XSYM-1)GO TO 5503
                  Xkj=Xkj+1
                  mSYM=IEOR(kjSYM,aSYM)+1
                  Fm=FLOV(mSYM,1)
                  Lm=FLOV(mSYM,2)
                  DO 5502 m=Fm,Lm
                     km=ITr(MAX0(k,m))+MIN0(k,m)
                     kmja=eOFF(ja)+eADD(km)
                     yARo(m,Xkj)=ein(kmja)
 5502             CONTINUE
 5503          CONTINUE
 5504       CONTINUE
            Xbc=0
            DO 5508 b=1,Nv
               IF(FZv(b).EQ.1)GO TO 5508
               bSYM=ORBSYM(b+no)
               DO 5507 c=1,Nv
                  IF(FZv(c).EQ.1)GO TO 5507
                  cSYM=ORBSYM(c+no)
                  BCSYM=IEOR(bSYM,cSYM)
                  IF(bcSYM.NE.bcSEM)GO TO 5507
                  ZLbc=ZLX(b,c)
                  Xbc=Xbc+1
                  mSYM=IEOR(bcSYM,iSYM)+1
                  Fm=FLOV(mSYM,1)
                  Lm=FLOV(mSYM,2)
                  DO 5506 m=Fm,Lm
                     imbc=UOFF(i,m,zlbc)+VADD(b,c)
                     yARv(m,Xbc)=z2o(imbc)
 5506             CONTINUE
 5507          CONTINUE
 5508       CONTINUE
            Xjk=0
            DO 5580 j=1,No
               IF(FZo(j).EQ.1)GO TO 5580
               jSYM=ORBSYM(j)
               DO 5575 k=1,NO
                  IF(FZO(k).EQ.1)GO TO 5575
                  kSYM=ORBSYM(k)
                  jkSYM=IEOR(jSYM,kSYM)
                  IF(jkSYM.NE.xsym-1)GO TO 5575
                  Xjk=Xjk+1
                  XCb=0
                  DO 5535 C=1,nv
                     IF(FZV(C).EQ.1)GO TO 5535
                     CSYM=ORBSYM(C+NO)
                     bSYM=IEOR(CSYM,bcsem)+1
                     Fb=FLOV(bSYM,3)-no
                     Lb=FLOV(bSYM,4)-no
                     mSYM=IEOR(aSYM,jkSYM)+1
                     Fm=FLOV(mSYM,1)
                     Lm=FLOV(mSYM,2)
                     DO 5530 b=fb,lb
                        VAL=0.0D0
                        XCb=XCb+1
                        zlbc=zlx(b,c)
                        DO 5525 m=Fm,Lm
                           VAL=VAL+yARv(m,XCb)*yARo(m,Xjk)
 5525                   CONTINUE
                        if(b.ne.c)then
                           JKBC=TOFF(J,K,zlbc,AIS)+TADD(B,C,AIS)
                           TAU(JKBC)=TAU(JKBC)-VAL
                        else
                           KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
                           TAU(KJCB)=TAU(KJCB)-VAL
                           JKBC=TOFF(J,K,1,AIS)+TADD(B,C,AIS)
                           TAU(JKBC)=TAU(JKBC)-VAL
                        endif
 5530                CONTINUE
 5535             CONTINUE
 5575       CONTINUE
 5580       CONTINUE
 5592       CONTINUE
C
           return
           end
