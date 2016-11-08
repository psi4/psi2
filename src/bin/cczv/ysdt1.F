      SUBROUTINE YSDT1(z2,z1,z2o,z1o,TAU,NO,NV,NDIMT2,NIT,
     .                 ORBSYM,FLOV,NIRRED,TOFF,TADD,ZLX,FZO,FZV,
     .                 UOFF,VADD,NT3,DOFF,DADD,EOFF,EADD,
     .                 FOFF,FADD,ITV,ITR,IPQ,DIN,EIN,FIN,fl82,FOCK,
     .                 AUX1,AUX3,SUX1,
     .                 XAR1,XAR2,MAXO,MAXV,
     .                 AR1,SUX11,XAR3,XAR4,AUX4,AUX5,
     .                 xaro,xarv,yaro,yarv,
     .                 buf,ibuf,intbuf,nshov,nslov,nslvo)
      IMPLICIT INTEGER(A-Z)
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
      REAL*8  z2(NDIMT2),z1(NO,NV),z2o(NDIMT2),z1o(NO,NV),
     .       TAU(*),FOCK(*),
     .       VAL,VAL1,VAL2,DIN(*),EIN(*),FIN(*),
     .       AUX1(NV,NO),AUX3(NV,NV),SUX11(NO*NV),
     .       AUX4(NO,NO),AUX5(NO,NV),SUX1(NV*NO),
     .       XAR1(MAXV,MAXOV),XAR2(MAXV,MAXOV),xaro(maxv,maxoo),
     .       XAR3(MAXO,MAXOV),XAR4(MAXO,MAXOV),xarv(maxv,maxvv),
     .       AR1(MAXOV,MAXOV),yaro(maxo,maxoo),yarv(maxo,maxvv)
      CHARACTER*4 SORT,fl82
      INTEGER ORBSYM(*),FLOV(NIRRED,4),TOFF(NO,NO,2,NIRRED),ZLX(NV,NV),
     .        TADD(NV,NV,NIRRED),FZO(NO),FZV(NV),UOFF(NO,NO,2),NT3(8),
     .        VADD(NV,NV),DOFF(*),DADD(*),FOFF(*),FADD(*),EOFF(*),
     .        EADD(*),ITR(*),ITV(*),IPQ(*)
C
c     SORT='NO  '
c
      call rdints(63,din,nshov,buf,ibuf,intbuf,no,nv,itr,itr,doff,dadd)
      call rdints(64,ein,nslov,buf,ibuf,intbuf,no,nv,itr,itr,eoff,eadd)
      call rdints(65,fin,nslvo,buf,ibuf,intbuf,no,nv,itv,itr,foff,fadd)
c
      itap82=82
      call rfile(itap82)
      i82=1
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
C  d INTS contribution to y3.
C
      call ydterm(a,asym,i,isym,ai,aisym,ais,dim,tau,din,
     .            fzo,fzv,orbsym,flov,nirred,no,nv,
     .            itr,itv,ipq,doff,dadd,toff,tadd,zlx,uoff,
     .            vadd,nt3,z1o)
c
C  F INTS contribution to y3.
C
      call yfterm(a,asym,i,isym,ai,aisym,ais,dim,tau,fin,
     .            fzo,fzv,orbsym,flov,nirred,no,nv,
     .            itr,itv,ipq,foff,fadd,toff,tadd,zlx,uoff,
     .            vadd,nt3,z2o,aux1,aux3,xar1,xar2,xaro,xarv,
     .            maxv)
C
C  E INTS contribution to t3.
C
      call yeterm(a,asym,i,isym,ai,aisym,ais,dim,tau,ein,
     .            fzo,fzv,orbsym,flov,nirred,no,nv,
     .            itr,itv,ipq,eoff,eadd,toff,tadd,zlx,uoff,
     .            vadd,nt3,z2o,aux4,aux5,xar3,xar4,yaro,yarv,
     .            maxo)
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
      if(fl82.eq.'YES ')then
         call wwritw(itap82,tau,intowp(dim),i82,i82)
         go to 9985
      endif
c 
c >>> ADD CCSDT TERMS TO z2o EQUATION. F INTS.
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
                        AR1(XBj,XGAV)=TAU(VjGAB)
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
                  BSYM=ORBSYM(B+NO)
                  DO 6225 j=1,NO
                     IF(FZO(j).EQ.1)GO TO 6225
                     jSYM=ORBSYM(j)
                     BjSYM=IEOR(BSYM,jSYM)
                     IF(BjSYM.NE.XSYM-1)GO TO 6225
                     XBj=XBj+1
                     jB=ITR(B)+j
                     BEAjB=FOFF(BEA)+FADD(jB)
                     SUX1(XBj)= FIN(BEAjB)
 6225             CONTINUE
 6230          CONTINUE
               DO 6240 XGAV=1,LIMGAV
                  VAL1=0.0D0
                  VAL2=0.0D0
                  DO 6235 XBj=1,LIMBj
                     VAL1=VAL1+SUX1(XBj)*AR1(XBj,XGAV)
 6235             CONTINUE
                  SUX11(XGAV)=VAL1
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
                     z2(iVBEGA)=z2(iVBEGA)+SUX11(XGAV)
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
                     z2(iVBEGA)=z2(iVBEGA)+SUX11(XGAV)
                     ViGABE=UOFF(V,i,1)+VADD(GA,BE)
                     z2(ViGABE)=z2(ViGABE)+SUX11(XGAV)
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
                     z2(ViGABE)=z2(ViGABE)+SUX11(XGAV)
 6247             CONTINUE
 6252          CONTINUE
 6270       CONTINUE
 6282    CONTINUE
c   
C >>> ADD CCSDT TERMS TO z2o EQUATION. E INTS.
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
                        AR1(XbJ,XGAV)=TAU(VJGAb)
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
                     bJ=ITR(b)+J
                     UIbJ=EOFF(bJ)+EADD(UI)
                     SUX1(XbJ)=-EIN(UIbJ)
 6425             CONTINUE
 6430          CONTINUE
               DO 6440 XGAV=1,LIMGAV
                  VAL1=0.0D0
                  VAL2=0.0D0
                  DO 6435 XbJ=1,LIMbJ
                     VAL1=VAL1+SUX1(XbJ)*AR1(XbJ,XGAV)
 6435             CONTINUE
                  SUX11(XGAV)=VAL1
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
                     z2(UVaGA)=z2(UVaGA)+SUX11(XGAV)
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
                     z2(UVaGA)=z2(UVaGA)+SUX11(XGAV)
                     VUGAa=UOFF(V,U,1)+VADD(GA,a)
                     z2(VUGAa)=z2(VUGAa)+SUX11(XGAV)
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
                     z2(VUGAa)=z2(VUGAa)+SUX11(XGAV)
 6447             CONTINUE
 6452          CONTINUE
 6470       CONTINUE
 6482    CONTINUE
c 
 9985    CONTINUE
 9990 CONTINUE
c
      call rclose(itap82,3)
c
C
          RETURN
          END
