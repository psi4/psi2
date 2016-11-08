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
C***********************************************************************
      SUBROUTINE coret3( S,T1N,T2 ,T1,TAU,NO,NV,NDIMT2,NIT,
     .                 ORBSYM,FLOV,NIRRED,TOFF,TADD,ZLX,FZO,FZV,
     .                 UOFF,VADD,OPTION,RSTR,NT3,DOFF,DADD,EOFF,EADD,
     .                 FOFF,FADD,ITV,ITR,IPQ,DIN,EIN,FIN,fl81,FOCK,
     .                 AUX1,AUX2,AUX3,SUX1,SUX2,SUX3,
     .                 XAR1,XAR2,NOACT,NVACT,MAXO,MAXV,MAXOV,
     .                 AR1,AR2,SUX11,SUX22,XAR3,XAR4,AUX4,AUX5,
     .                 t3arr,i1off,i2off,i3off,ndimt3,trsym)
      IMPLICIT INTEGER(A-Z)
      common/sheet/xnm,xirred
      REAL*8  S(NDIMT2),T1N(NO,NV),T2(NDIMT2),T1(NO,NV),TAU(*),FOCK(*),
     .       VAL,VAL1,VAL2,DIN(*),EIN(*),FIN(*),tnorm,
     .       AUX1(NV,NO),AUX2(NV,NO),AUX3(NV,NV),SUX11(NO*NV),
     .       AUX4(NO,NO),AUX5(NO,NV),
     .       SUX1(NV*NO),SUX2(NV*NO),SUX3(NV*NV),SUX22(NO*NV),
     .       XAR1(MAXV,MAXOV),XAR2(MAXV,MAXOV),
     .       XAR3(MAXO,MAXOV),XAR4(MAXO,MAXOV),
     .       AR1(MAXOV,MAXOV),AR2(MAXOV,MAXOV),
     .       t3arr(*)
      CHARACTER*4 OPTION,RSTR,fl81
      INTEGER ORBSYM(*),FLOV(NIRRED,4),TOFF(NO,NO,2,NIRRED),ZLX(NV,NV),
     .        TADD(NV,NV,NIRRED),FZO(NO),FZV(NV),UOFF(NO,NO,2),NT3(8),
     .        VADD(NV,NV),DOFF(*),DADD(*),FOFF(*),FADD(*),EOFF(*),
     .        EADD(*),ITR(*),ITV(*),IPQ(*),
     .        i1off(*),i2off(no*nv,nirred),i3off(*),trsym(*)
C
      xirred=nirred
      xnm=no*nv
c
      do 100 a=1,nv
         asym=orbsym(a+no)
         do 101 i=1,no
            isym=orbsym(i)
            ais=ieor(asym,isym)+1
            ai=itr(a)+i
            trsym(ai)=ais
  101    continue
  100 continue
c

      CALL ZERO (t3arr,ndimt3)
C
C  F INTS contribution to t3.
C
      DO 3990 A=1,NV
         IF(FZV(A).EQ.1)GO TO 3990
         XA=IPQ(A+NO)+A+NO
         ASYM=ORBSYM(A+NO)
         DO 3985 I=1,NO
            IF(FZO(I).EQ.1)GO TO 3985
            XI=IPQ(I)+I
            ISYM=ORBSYM(I)
            AI=ITR(A)+I
            AISYM=IEOR(ISYM,ASYM)
            AIS=AISYM+1
c
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
                  ck=itr(c)+k
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
C                 DO 3935 B=1,C
                  DO 3935 B=1,nv
                     IF(FZV(B).EQ.1)GO TO 3935
                     BSYM=ORBSYM(B+NO)
                     JSYM=IEOR(BSYM,AICKS)+1
                     FJ=FLOV(JSYM,1)
                     LJ=FLOV(JSYM,2)
                     ESYM=IEOR(BSYM,AISYM)+1
                     FE=FLOV(ESYM,3)-NO
                     LE=FLOV(ESYM,4)-NO
                     DO 3930 J=FJ,LJ
                     bj=itr(b)+j
                     VAL=0.0D0
                        DO 3925 E=FE,LE
                        VAL=VAL+AUX3(E,B)*AUX1(E,J)
 3925                   CONTINUE
c                    call indt3(ai,bj,ck,aibjck,i1off,i2off,i3off,
c    .                          trsym)
                     k1=max0(ai,bj,ck)
                     k3=min0(ai,bj,ck)
                     k3s=trsym(k3)
                     k2=ai+bj+ck-k1-k3
                     aibjck=i1off(k1)+i2off(k2,k3s)+i3off(k3)
                     t3arr(aibjck)=t3arr(aibjck)+val
 3930                CONTINUE
 3935             CONTINUE
 3975          CONTINUE
 3980       CONTINUE
 3985       CONTINUE
 3990       CONTINUE
C
C  E INTS contribution to t3.
C
      DO 4990 A=1,NV
         IF(FZV(A).EQ.1)GO TO 4990
         XA=IPQ(A+NO)+A+NO
         ASYM=ORBSYM(A+NO)
         DO 4985 I=1,NO
            IF(FZO(I).EQ.1)GO TO 4985
            XI=IPQ(I)+I
            ISYM=ORBSYM(I)
            AI=ITR(A)+I
            AISYM=IEOR(ISYM,ASYM)
            AIS=AISYM+1
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
                  ck=itr(c)+k
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
                     MSYM=IEOR(JSYM,AISYM)+1
                     FM=FLOV(MSYM,1)
                     LM=FLOV(MSYM,2)
                     DO 4930 B=FB,LB
                     bj=itr(b)+j
                     VAL=0.0D0
                        DO 4925 M=FM,LM
                        VAL=VAL+AUX4(M,J)*AUX5(M,B)
 4925                   CONTINUE
c                    KJCB=TOFF(K,J,1,AIS)+TADD(C,B,AIS)
c                    TAU(KJCB)=TAU(KJCB)-VAL
c                    call indt3(ai,bj,ck,aibjck,i1off,i2off,i3off,
c    .                          trsym)
                     k1=max0(ai,bj,ck)
                     k3=min0(ai,bj,ck)
                     k3s=trsym(k3)
                     k2=ai+bj+ck-k1-k3
                     aibjck=i1off(k1)+i2off(k2,k3s)+i3off(k3)
                     t3arr(aibjck)=t3arr(aibjck)-val
 4930                CONTINUE
 4935             CONTINUE
 4950             CONTINUE
 4975          CONTINUE
 4980       CONTINUE
c
 4985    CONTINUE
 4990 CONTINUE
c
C
C  DIVIDE BY MO ENERGIES
C
      ind=0
      DO 3590 A=1,NV
         IF(FZV(A).EQ.1)GO TO 3590
         XA=IPQ(A+NO)+A+NO
         ASYM=ORBSYM(A+NO)
         DO 3585 I=1,NO
            IF(FZO(I).EQ.1)GO TO 3585
            XI=IPQ(I)+I
            ISYM=ORBSYM(I)
            AI=ITR(A)+I
            AISYM=IEOR(ISYM,ASYM)
            DO 3570 B=1,NV
               IF(FZV(B).EQ.1)GO TO 3570
               XB=IPQ(B+NO)+B+NO
               BSYM=ORBSYM(B+NO)
               DO 3560 J=1,NO
                  IF(FZO(J).EQ.1)GO TO 3560
                  XJ=IPQ(J)+J
                  JSYM=ORBSYM(J)
                  bJSYM=IEOR(bSYM,JSYM)
                  aibjs=IEOR(AiSYM,BjSYM)
                  bj=itr(b)+j
                  if(bj.gt.ai)go to 3560
                  DO 3550 c=1,Nv
                     IF(FZv(c).EQ.1)GO TO 3550
                     Xc=IPQ(c+no)+c+no
                     cSYM=ORBSYM(c+no)
                     kSYM=IEOR(aibjs,csYM)+1
                     Fk=FLOV(kSYM,1)
                     Lk=FLOV(kSYM,2)
                     DO 3540 k=Fk,Lk
                        Xk=IPQ(k)+k
                        ck=itr(c)+k
                        if(ck.gt.bj)go to 3540
                        VAL=FOCK(XI)+FOCK(XJ)+FOCK(XK)
     .                     -FOCK(XA)-FOCK(XB)-FOCK(XC)
                        ind=ind+1
                        t3arr(ind)=t3arr(ind)/val
c           call indt3(ai,bj,ck,aibjck,i1off,i2off,i3off,trsym)
c           t3arr(aibjck)=t3arr(aibjck)/val
c           write(6,675)a+no,i,b+no,j,c+no,k,t3arr(aibjck)
  675       format('div',6i3,f12.8)
 3540                CONTINUE
 3550             CONTINUE
 3560          CONTINUE
 3570       CONTINUE
 3585    CONTINUE
 3590 CONTINUE
c
      if(fl81.eq.'GRAD')then
      itap81=81
      call rfile(itap81)
      i81=1
C
      DO 9990 A=1,NV
         IF(FZV(A).EQ.1)GO TO 9990
         ASYM=ORBSYM(A+NO)
         DO 9985 I=1,NO
            IF(FZO(I).EQ.1)GO TO 9985
            ISYM=ORBSYM(I)
            ai=itr(a)+i
            AISYM=IEOR(ISYM,ASYM)
            AIS=AISYM+1
            DIM=NT3(AIS)
            DO 9970 J=1,NO
            IF(FZO(J).EQ.1)GO TO 9970
            JSYM=ORBSYM(J)
            IJSYM=IEOR(ISYM,JSYM)
            DO 9960 B=1,NV
               IF(FZV(B).EQ.1)GO TO 9960
               BSYM=ORBSYM(B+NO)
               ABSYM=IEOR(ASYM,BSYM)
               bj=itr(b)+j
               DO 9950 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 9950
                  KSYM=ORBSYM(K)
                  IJKSYM=IEOR(IJSYM,KSYM)
                  CSYM=IEOR(IJKSYM,ABSYM)+1
                  FC=FLOV(CSYM,3)-NO
                  LC=FLOV(CSYM,4)-NO
                  LC2=LC
                  IF(LC.GT.B)LC2=B
                  DO 9952 C=FC,LC2
                     IF(FZV(C).EQ.1)GO TO 9952
                     ck=itr(c)+k
                     ZLBC=ZLX(B,C)
                     JKBC=TOFF(J,K,ZLBC,AIS)+TADD(B,C,AIS)
                     k1=max0(ai,bj,ck)
                     k3=min0(ai,bj,ck)
                     k3s=trsym(k3)
                     k2=ai+bj+ck-k1-k3
                     aibjck=i1off(k1)+i2off(k2,k3s)+i3off(k3)
                     TAU(JKBC)=t3arr(aibjck)   
 9952                CONTINUE
 9950             CONTINUE
 9960          CONTINUE
 9970       CONTINUE
            call wwritw(itap81,tau,intowp(dim),i81,i81)
 9985    CONTINUE
 9990 CONTINUE
      call rclose(itap81,3)
      endif
C
      DO 3690 A=1,NV
         IF(FZV(A).EQ.1)GO TO 3690
         XA=IPQ(A+NO)+A+NO
         ASYM=ORBSYM(A+NO)
         DO 3685 I=1,NO
            IF(FZO(I).EQ.1)GO TO 3685
            XI=IPQ(I)+I
            ISYM=ORBSYM(I)
            AI=ITR(A)+I
            if(asym.ne.isym)go to 3685
            Xa=IPQ(a+NO)+a+NO
            aiaiai=i1off(ai)+i2off(ai,1)+i3off(ai)
            t3arr(aiaiai)=t3arr(aiaiai)*6.0d0
c           write(6,675)a+no,i,a+no,i,a+no,i,t3arr(aiaiai)
 3685    CONTINUE
 3690 CONTINUE
C
C
      DO 3790 A=1,NV
         IF(FZV(A).EQ.1)GO TO 3790
         XA=IPQ(A+NO)+A+NO
         ASYM=ORBSYM(A+NO)
         DO 3785 I=1,NO
            IF(FZO(I).EQ.1)GO TO 3785
            XI=IPQ(I)+I
            ISYM=ORBSYM(I)
            AI=ITR(A)+I
            ais=ieor(asym,isym)+1
            do 3780 b=1,nv
               IF(FZV(b).EQ.1)GO TO 3780
               Xb=IPQ(b+NO)+b+NO
               bSYM=ORBSYM(b+NO)
               DO 3775 j=1,NO
                  IF(FZO(j).EQ.1)GO TO 3775
                  jSYM=ORBSYM(j)
                  if(bsym.ne.jsym)go to 3775
                  bj=ITR(b)+j
                  if(ai.eq.bj)go to 3775
                  Xj=IPQ(j)+j
                  bj=ITR(b)+j
                  VAL=(FOCK(XI)+fock(xi)+fock(xj)
     .               - FOCK(XA)-fock(xa)-fock(xb))*0.5d0
                  if(ai.gt.bj)then
                     aiaibj=i1off(ai)+i2off(ai,1)+i3off(bj)
                     t3arr(aiaibj)=t3arr(aiaibj)*2.0d0
c                  write(6,675)a+no,i,a+no,i,b+no,j,t3arr(aiaibj)
                  else
                     bjaiai=i1off(bj)+i2off(ai,ais)+i3off(ai)
                     t3arr(bjaiai)=t3arr(bjaiai)*2.0d0
c                  write(6,675)b+no,j,a+no,i,a+no,i,t3arr(bjaiai)
                  endif
 3775          CONTINUE
 3780       CONTINUE
 3785    CONTINUE
 3790 CONTINUE
C
C
C >>> ADD CCSDT TERMS TO T1 EQUATION. D INTS.
C
      DO 3190 A=1,NV
         IF(FZV(A).EQ.1)GO TO 3190
         XA=IPQ(A+NO)+A+NO
         ASYM=ORBSYM(A+NO)
         DO 3185 I=1,NO
            IF(FZO(I).EQ.1)GO TO 3185
            XI=IPQ(I)+I
            ISYM=ORBSYM(I)
            AI=ITR(A)+I
            AISYM=IEOR(ISYM,ASYM)
            AIS=AISYM+1
c
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
c     ZLBEB=ZLX(BE,B)
c     ZLBBE=ZLX(B,BE)
      IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
      IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
c     JUBEB=TOFF(J,U,ZLBEB,ais)+TADD(BE,B,ais)
c     JUBBE=TOFF(J,U,ZLBBE,ais)+TADD(B,BE,ais)
      bej=itr(be)+j
      beu=itr(be)+u
      bj=itr(b)+j
      bu=itr(b)+u
c
c     call indt3(ai,bej,bu,jubeb,i1off,i2off,i3off,
c    .                          trsym)
         k1=max0(ai,bej,bu)
         k3=min0(ai,bej,bu)
         k3s=trsym(k3)
         k2=ai+bej+bu-k1-k3
         jubeb=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c
c     call indt3(ai,beu,bj,jubbe,i1off,i2off,i3off,
c    .                          trsym)
         k1=max0(ai,beu,bj)
         k3=min0(ai,beu,bj)
         k3s=trsym(k3)
         k2=ai+beu+bj-k1-k3
         jubbe=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      T1N(U,BE)=T1N(U,BE)+(DIN(IAJB)+DIN(IAJB)-DIN(IBJA))
     .                   *(t3arr(JUBBE)-t3arr(JUBEB))
c    .                   *(TAU(JUBBE)-TAU(JUBEB))
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
                  bj=itr(b)+j
                  XGAV=0
                  DO 6210 GA=1,NV
                     IF(FZV(GA).EQ.1)GO TO 6210
                     GASYM=ORBSYM(GA+NO)
                     VSYM=IEOR(GASYM,AiBjS)+1
                     FV=FLOV(VSYM,1)
                     LV=FLOV(VSYM,2)
                     ZLGAB=ZLX(GA,B)
                     gaj=itr(ga)+j
                     DO 6205 V=FV,LV
                        gav=itr(ga)+v
                        bv=itr(b)+v
                        XGAV=XGAV+1
         k1=max0(ai,gav,bj)
         k3=min0(ai,gav,bj)
         k3s=trsym(k3)
         k2=ai+gav+bj-k1-k3
         vjgab=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c                       call indt3(ai,gav,bj,vjgab,i1off,i2off,i3off,
c    .                          trsym)
         k1=max0(ai,gaj,bv)
         k3=min0(ai,gaj,bv)
         k3s=trsym(k3)
         k2=ai+gaj+bv-k1-k3
         vjbga=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c                       call indt3(ai,gaj,bv,vjbga,i1off,i2off,i3off,
c    .                          trsym)
                        AR1(XBj,XGAV)=t3arr(VjGAB)
                        AR2(XBj,XGAV)=t3arr(VjBGA)
c                       VjGAB=TOFF(V,j,ZLGAB,AiS)+TADD(GA,B,AiS)
c                       VjBGA=TOFF(j,V,ZLGAB,AiS)+TADD(GA,B,AiS)
c                       AR1(XBj,XGAV)=TAU(VjGAB)
c                       AR2(XBj,XGAV)=TAU(VjBGA)
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
                  bj=itr(b)+j
                  XGAV=0
                  DO 6410 GA=1,NV
                     IF(FZV(GA).EQ.1)GO TO 6410
                     GASYM=ORBSYM(GA+NO)
                     VSYM=IEOR(GASYM,aIbJS)+1
                     FV=FLOV(VSYM,1)
                     LV=FLOV(VSYM,2)
                     gaj=itr(ga)+j
                     ZLGAb=ZLX(GA,b)
                     DO 6405 V=FV,LV
                        gav=itr(ga)+v
                        bv=itr(b)+v
                        XGAV=XGAV+1
         k1=max0(ai,gaj,bv)
         k3=min0(ai,gaj,bv)
         k3s=trsym(k3)
         k2=ai+gaj+bv-k1-k3
         vjbga=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c                       call indt3(ai,gaj,bv,vjbga,i1off,i2off,i3off,
c    .                          trsym)
         k1=max0(ai,gav,bj)
         k3=min0(ai,gav,bj)
         k3s=trsym(k3)
         k2=ai+gav+bj-k1-k3
         vjgab=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c                       call indt3(ai,gav,bj,vjgab,i1off,i2off,i3off,
c    .                          trsym)
                        AR1(XbJ,XGAV)=t3arr(VJGAb)
                        AR2(XbJ,XGAV)=t3arr(VJbGA)
c                       VJGAb=TOFF(V,J,ZLGAb,aIS)+TADD(GA,b,aIS)
c                       VJbGA=TOFF(J,V,ZLGAb,aIS)+TADD(GA,b,aIS)
c                       AR1(XbJ,XGAV)=TAU(VJGAb)
c                       AR2(XbJ,XGAV)=TAU(VJbGA)
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
 3185    CONTINUE
 3190 CONTINUE
c
          RETURN
          END
