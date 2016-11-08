      SUBROUTINE t3tots(ITR,itv,din,ein,fin,
     .                  doff,dadd,eoff,eadd,foff,fadd,
     .                  NO,NV,NT,NT2,NDIMT2,ndimt3,nirred,level,
     .                  t1o,t1,t2o,t2,tau,t3o,t3,
     .                  ORBSYM,FLOV,UOFF,VADD,ZLX,FZO,FZV,FOCK,
     .                  i1off,i2off,i3off,trsym)
      IMPLICIT INTEGER (A-Z)
      REAL*8 VAL,VAL1,VAL2,VAL3,xfac
      real*8 t3o(*),t3(*),din(*),ein(*),fin(*), 
     .       t2(NDIMT2),t2o(NDIMT2),TAU(NDIMT2),T1O(NO,NV),t1(NO,NV),
     .       fock(*)
      integer ORBSYM(NT),FLOV(NIRRED,4),
     .        UOFF(NO,NO,2),VADD(NV,NV),ZLX(NV,NV),
     .        ITR(NV),FZO(NO),FZV(NV),itv(*),
     .        i1off(*),i2off(no*nv,nirred),i3off(*),trsym(*),
     .        doff(*),dadd(*),eoff(*),eadd(*),foff(*),fadd(*)
C
      DO 3190 U=1,NO
      USYM=ORBSYM(U)
      IF(FZO(U).EQ.1)GO TO 3190
      FBE=FLOV(USYM+1,3)-NO
      LBE=FLOV(USYM+1,4)-NO
      DO 3185 BE=FBE,LBE
      beu=itr(be)+u
      DO 3180 I=1,NO
      IF(FZO(I).EQ.1)GO TO 3180
      ISYM=ORBSYM(I)
      DO 3175 J=1,NO
      IF(FZO(J).EQ.1)GO TO 3175
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      bej=itr(be)+j
      DO 3170 A=1,NV
      IF(FZV(A).EQ.1)GO TO 3170
      ASYM=ORBSYM(A+NO)
      ia=itr(a)+i
      ja=itr(a)+j
      ai=itr(a)+i
      BSYM=IEOR(IJSYM,ASYM)+1
      FB=FLOV(BSYM,3)-NO
      LB=FLOV(BSYM,4)-NO
      DO 3165 B=FB,LB
      ib=itr(b)+i
      jb=itr(b)+j
      IAJB=doff(MAX0(IA,JB))+dadd(MIN0(IA,JB))
      IBJA=doff(MAX0(IB,JA))+dadd(MIN0(IB,JA))
      bj=itr(b)+j
      bu=itr(b)+u
c     IDX1=X3OFF(I,J,U)+X3ADD(A,BE,B)
           k1=max0(ai,bej,bu)
           k3=min0(ai,bej,bu)
           k3s=trsym(k3)
           k2=ai+bej+bu-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(I,J,U)+X3ADD(A,B,BE)
           k1=max0(ai,bj,beu)
           k3=min0(ai,bj,beu)
           k3s=trsym(k3)
           k2=ai+bj+beu-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t1(U,BE)=t1(U,BE)
     .         -(din(IAJB)+din(IAJB)-din(IBJA))*(t3o(IDX1)-t3o(IDX2))         
 3165 CONTINUE
 3170 CONTINUE
 3175 CONTINUE
 3180 CONTINUE
 3185 CONTINUE
 3190 CONTINUE
C
      DO 3290 U=1,NO
      IF(FZO(U).EQ.1)GO TO 3290
      USYM=ORBSYM(U)
      DO 3285 V=1,NO
      IF(FZO(V).EQ.1)GO TO 3285
      VSYM=ORBSYM(V)
      UVSYM=IEOR(USYM,VSYM)
      DO 3280 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 3280
      BESYM=ORBSYM(BE+NO)
      GASYM=IEOR(UVSYM,BESYM)+1
      FGA=FLOV(GASYM,3)-NO
      LGA=FLOV(GASYM,4)-NO
C
      LGA2=LGA
      IF(LGA.GT.BE)LGA2=BE
      DO 3276 GA=FGA,LGA2
      UVBEGA=UOFF(U,V,1)+VADD(BE,GA)
      gav=itr(ga)+v
      DO 3271 I=1,NO
      IF(FZO(I).EQ.1)GO TO 3271
      ISYM=ORBSYM(I)
      bei=itr(be)+i
      gai=itr(ga)+i
      IU=itr(MAX0(I,U))+MIN0(I,U)
      DO 3266 A=1,NV
      IF(FZV(A).EQ.1)GO TO 3266
      ASYM=ORBSYM(A+NO)
      IASYM=IEOR(ISYM,ASYM)
      IA=itr(A)+I
      av=itr(a)+v
      au=itr(a)+u
      BEA=itv(MAX0(BE,A))+MIN0(BE,A)
      JSYM=IEOR(IASYM,USYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      DO 3241 J=FJ,LJ
      JU=itr(MAX0(J,U))+MIN0(J,U)
      JA=itr(A)+J
      IUJA=eoff(JA)+eadd(IU)
      IAJU=eoff(IA)+eadd(JU)
      aj=itr(a)+j
      gaj=itr(ga)+j
c     IDX1=X3OFF(I,J,V)+X3ADD(BE,A,GA)
           k1=max0(aj,bei,gav)
           k3=min0(aj,bei,gav)
           k3s=trsym(k3)
           k2=aj+bei+gav-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(I,J,V)+X3ADD(BE,GA,A)
           k1=max0(gaj,bei,av)
           k3=min0(gaj,bei,av)
           k3s=trsym(k3)
           k2=gaj+bei+av-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t2(UVBEGA)=t2(UVBEGA)-(ein(IUJA)+ein(IUJA)-ein(IAJU))*t3o(IDX1)
     .                   + ein(IUJA) * t3o(IDX2)
 3241 CONTINUE
c
      BSYM=IEOR(IASYM,BESYM)+1
      FB=FLOV(BSYM,3)-NO
      LB=FLOV(BSYM,4)-NO
      DO 3251 B=FB,LB
      IB=itr(B)+I
      bi=itr(b)+i
      bv=itr(b)+v
      BEB=Itv(MAX0(BE,B))+MIN0(BE,B)
      BEBIA=foff(BEB)+fadd(IA)
      BAIB=foff(BEA)+fadd(IB)
c     IDX1=X3OFF(U,I,V)+X3ADD(A,B,GA)
           k1=max0(au,bi,gav)
           k3=min0(au,bi,gav)
           k3s=trsym(k3)
           k2=au+bi+gav-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(U,I,V)+X3ADD(A,GA,B)
           k1=max0(gai,bv,au)
           k3=min0(gai,bv,au)
           k3s=trsym(k3)
           k2=gai+bv+au-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t2(UVBEGA)=t2(UVBEGA)+(fin(BAIB)+fin(BAIB)-fin(BEBIA))*t3o(IDX1)
     .                   - fin(BAIB) * t3o(IDX2)
 3251 CONTINUE
 3266 CONTINUE
 3271 CONTINUE
 3276 CONTINUE
C
      FGA2=FGA
      IF(FGA.LT.BE)FGA2=BE
      DO 3376 GA=FGA2,LGA
      gav=itr(ga)+v
      VUGABE=UOFF(V,U,1)+VADD(GA,BE)
      DO 3371 I=1,NO
      IF(FZO(I).EQ.1)GO TO 3371
      ISYM=ORBSYM(I)
      bei=itr(be)+i
      gai=itr(ga)+i
      IU=itr(MAX0(I,U))+MIN0(I,U)
      DO 3366 A=1,NV
      IF(FZV(A).EQ.1)GO TO 3366
      ASYM=ORBSYM(A+NO)
      IASYM=IEOR(ISYM,ASYM)
      IA=itr(A)+I
      av=itr(a)+v
      au=itr(a)+u
      BEA=itv(MAX0(BE,A))+MIN0(BE,A)
      JSYM=IEOR(IASYM,USYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      DO 3341 J=FJ,LJ
      JU=itr(MAX0(J,U))+MIN0(J,U)
      JA=itr(A)+J
      IUJA=eoff(JA)+eadd(IU)
      IAJU=eoff(IA)+eadd(JU)
      aj=itr(a)+j
      gaj=itr(ga)+j
c     IDX1=X3OFF(I,J,V)+X3ADD(BE,A,GA)
           k1=max0(aj,bei,gav)
           k3=min0(aj,bei,gav)
           k3s=trsym(k3)
           k2=aj+bei+gav-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(I,J,V)+X3ADD(BE,GA,A)
           k1=max0(gaj,bei,av)
           k3=min0(gaj,bei,av)
           k3s=trsym(k3)
           k2=gaj+bei+av-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t2(VUGABE)=t2(VUGABE)-(ein(IUJA)+ein(IUJA)-ein(IAJU))*t3o(IDX1)
     .                   + ein(IUJA) * t3o(IDX2)
 3341 CONTINUE
      BSYM=IEOR(IASYM,BESYM)+1
      FB=FLOV(BSYM,3)-NO
      LB=FLOV(BSYM,4)-NO
      DO 3351 B=FB,LB
      IB=itr(B)+I
      BEB=Itv(MAX0(BE,B))+MIN0(BE,B)
      BEBIA=foff(BEB)+fadd(IA)
      BAIB=foff(BEA)+fadd(IB)
      bi=itr(b)+i
      bv=itr(b)+v
c     IDX1=X3OFF(U,I,V)+X3ADD(A,B,GA)
           k1=max0(au,bi,gav)
           k3=min0(au,bi,gav)
           k3s=trsym(k3)
           k2=au+bi+gav-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(U,I,V)+X3ADD(A,GA,B)
           k1=max0(gai,bv,au)
           k3=min0(gai,bv,au)
           k3s=trsym(k3)
           k2=gai+bv+au-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t2(VUGABE)=t2(VUGABE)+(fin(BAIB)+fin(BAIB)-fin(BEBIA))*t3o(IDX1)
     .                   - fin(BAIB) * t3o(IDX2)
 3351 CONTINUE
 3366 CONTINUE
 3371 CONTINUE
 3376 CONTINUE
C
 3280 CONTINUE
 3285 CONTINUE
 3290 CONTINUE
C
C     ADD NON-LINEAR T3*T1 CONTRIBUTION TO T2 EQUATION (CCSDT-1B)
C
      DO 6398 I=1,NO
      IF(FZO(I).EQ.1)GO TO 6398
      ISYM=ORBSYM(I)
      FA=FLOV(ISYM+1,3)-NO
      LA=FLOV(ISYM+1,4)-NO
      DO 6395 A=FA,LA
      IA=itr(A)+I
      ai=itr(a)+i
      VAL=0.0D0
      DO 6350 J=1,NO
      IF(FZO(J).EQ.1)GO TO 6350
      JSYM=ORBSYM(J)
      FB=FLOV(JSYM+1,3)-NO
      LB=FLOV(JSYM+1,4)-NO
      JA=itr(A)+J
      DO 6345 B=FB,LB
      JB=itr(B)+J
      IAJB=doff(MAX0(IA,JB))+dadd(MIN0(IA,JB))
      IB=itr(B)+I
      IBJA=doff(MAX0(IB,JA))+dadd(MIN0(IB,JA))
      VAL=VAL+(din(IAJB)+din(IAJB)-din(IBJA))*T1O(J,B)
 6345 CONTINUE
 6350 CONTINUE
      DO 6390 U=1,NO
      IF(FZO(U).EQ.1)GO TO 6390
      au=itr(a)+u
      USYM=ORBSYM(U)
      DO 6385 V=1,NO
      IF(FZO(V).EQ.1)GO TO 6385
      VSYM=ORBSYM(V)
      UVSYM=IEOR(USYM,VSYM)
      DO 6380 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 6380
      beu=itr(be)+u
      bei=itr(be)+i
      BESYM=ORBSYM(BE+NO)
      GASYM=IEOR(UVSYM,BESYM)+1
      FGA=FLOV(GASYM,3)-NO
      LGA=FLOV(GASYM,4)-NO
C
      LGA2=LGA
      IF(LGA.GT.BE)LGA2=BE
      DO 6375 GA=FGA,LGA2
      UVBEGA=UOFF(U,V,1)+VADD(BE,GA)
      gav=itr(ga)+v
c     IDX1=X3OFF(U,V,I)+X3ADD(BE,GA,A)
           k1=max0(beu,ai,gav)
           k3=min0(beu,ai,gav)
           k3s=trsym(k3)
           k2=beu+ai+gav-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(U,V,I)+X3ADD(A,GA,BE)
           k1=max0(bei,au,gav)
           k3=min0(bei,au,gav)
           k3s=trsym(k3)
           k2=bei+au+gav-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t2(UVBEGA)=t2(UVBEGA)+VAL* (t3o(IDX1)-t3o(IDX2))
 6375 CONTINUE
C
      FGA2=FGA
      IF(FGA.LT.BE)FGA2=BE
      DO 6377 GA=FGA2,LGA
      VUGABE=UOFF(V,U,1)+VADD(GA,BE)
      gav=itr(ga)+v
c     IDX1=X3OFF(U,V,I)+X3ADD(BE,GA,A)
           k1=max0(beu,ai,gav)
           k3=min0(beu,ai,gav)
           k3s=trsym(k3)
           k2=beu+ai+gav-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(U,V,I)+X3ADD(A,GA,BE)
           k1=max0(bei,au,gav)
           k3=min0(bei,au,gav)
           k3s=trsym(k3)
           k2=bei+au+gav-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t2(VUGABE)=t2(VUGABE)+VAL*(t3o(IDX1)-t3o(IDX2))
 6377 CONTINUE
C
 6380 CONTINUE
 6385 CONTINUE
 6390 CONTINUE
 6395 CONTINUE
 6398 CONTINUE
C
      DO 6498 I=1,NO
      IF(FZO(I).EQ.1)GO TO 6498
      ISYM=ORBSYM(I)
      DO 6495 A=1,NV
      IF(FZV(A).EQ.1)GO TO 6495
      ASYM=ORBSYM(A+NO)
      IA=itr(A)+I
      ai=itr(a)+i
      IASYM=IEOR(ISYM,ASYM)
      DO 6490 J=1,NO
      IF(FZO(J).EQ.1)GO TO 6490
      JA=itr(A)+J
      JSYM=ORBSYM(J)
      VSYM=IEOR(IASYM,JSYM)
      IV=FLOV(VSYM+1,1)
      LV=FLOV(VSYM+1,2)
      FB=FLOV(VSYM+1,3)-NO
      LB=FLOV(VSYM+1,4)-NO
      DO 6480 V=IV,LV
      VAL1=0.0D0
      VAL2=0.0D0
      DO 6445 B=FB,LB
      JB=itr(B)+J
      IAJB=doff(MAX0(IA,JB))+dadd(min0(IA,JB))
      IB=itr(B)+I
      IBJA=doff(MAX0(IB,JA))+dadd(min0(IB,JA))
      VAL1=VAL1+(din(IAJB)+din(IAJB)-din(IBJA))*T1O(V,B)
      VAL2=VAL2+ din(IAJB)                   *T1O(V,B)
 6445 CONTINUE
      DO 6470 U=1,NO
      IF(FZO(U).EQ.1)GO TO 6470
      au=itr(a)+u
      USYM=ORBSYM(U)
      UVSYM=IEOR(USYM,VSYM)
      DO 6460 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 6460
      beu=itr(be)+u
      bei=itr(be)+i
      BESYM=ORBSYM(BE+NO)
      GASYM=IEOR(UVSYM,BESYM)+1
      FGA=FLOV(GASYM,3)-NO
      LGA=FLOV(GASYM,4)-NO
C
      LGA2=LGA
      IF(LGA.GT.BE)LGA2=BE
      DO 6455 GA=FGA,LGA2
      UVBEGA=UOFF(U,V,1)+VADD(BE,GA)
      gaj=itr(ga)+j
c     IDX1=X3OFF(U,I,J)+X3ADD(BE,A,GA)
           k1=max0(beu,ai,gaj)
           k3=min0(beu,ai,gaj)
           k3s=trsym(k3)
           k2=beu+ai+gaj-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(U,I,J)+X3ADD(A,BE,GA)
           k1=max0(bei,au,gaj)
           k3=min0(bei,au,gaj)
           k3s=trsym(k3)
           k2=bei+au+gaj-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t2(UVBEGA)=t2(UVBEGA)-VAL1* t3o(IDX1)+VAL2* t3o(IDX2)
 6455 CONTINUE
C
      FGA2=FGA
      IF(FGA.LT.BE)FGA2=BE
      DO 6457 GA=FGA2,LGA
      VUGABE=UOFF(V,U,1)+VADD(GA,BE)
      gaj=itr(ga)+j
c     IDX1=X3OFF(U,I,J)+X3ADD(BE,A,GA)
           k1=max0(beu,ai,gaj)
           k3=min0(beu,ai,gaj)
           k3s=trsym(k3)
           k2=beu+ai+gaj-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(U,I,J)+X3ADD(A,BE,GA)
           k1=max0(bei,au,gaj)
           k3=min0(bei,au,gaj)
           k3s=trsym(k3)
           k2=bei+au+gaj-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t2(VUGABE)=t2(VUGABE)-VAL1* t3o(IDX1)+VAL2* t3o(IDX2)
 6457 CONTINUE
C
 6460 CONTINUE
 6470 CONTINUE
 6480 CONTINUE
 6490 CONTINUE
 6495 CONTINUE
 6498 CONTINUE
C
      DO 6598 I=1,NO
      IF(FZO(I).EQ.1)GO TO 6598
      ISYM=ORBSYM(I)
      DO 6595 A=1,NV
      IF(FZV(A).EQ.1)GO TO 6595
      ASYM=ORBSYM(A+NO)
      IA=itr(A)+I
      ai=itr(a)+i
      IASYM=IEOR(ISYM,ASYM)
      DO 6590 B=1,NV
      IF(FZV(B).EQ.1)GO TO 6590
      IB=itr(B)+I
      BSYM=ORBSYM(B+NO)
      GASYM=IEOR(IASYM,BSYM)
      FJ=FLOV(GASYM+1,1)
      LJ=FLOV(GASYM+1,2)
      FGA=FLOV(GASYM+1,3)-NO
      LGA=FLOV(GASYM+1,4)-NO
      DO 6580 GA=FGA,LGA
      VAL1=0.0D0
      VAL2=0.0D0
      DO 6545 J=FJ,LJ
      JB=itr(B)+J
      JA=itr(A)+J
      IAJB=doff(MAX0(IA,JB))+dadd(MIN0(IA,JB))
      IBJA=doff(MAX0(IB,JA))+dadd(MIN0(IB,JA))
      VAL1=VAL1+(din(IAJB)+din(IAJB)-din(IBJA))*T1O(J,GA)
      VAL2=VAL2+ din(IAJB)                   *T1O(J,GA)
 6545 CONTINUE
      DO 6570 U=1,NO
      IF(FZO(U).EQ.1)GO TO 6570
      USYM=ORBSYM(U)
      au=itr(a)+u
      DO 6560 V=1,NO
      IF(FZO(V).EQ.1)GO TO 6560
      bv=itr(b)+v
      VSYM=ORBSYM(V)
      UVSYM=IEOR(USYM,VSYM)
      BESYM=IEOR(UVSYM,GASYM)+1
      FBE=FLOV(BESYM,3)-NO
      LBE=FLOV(BESYM,4)-NO
C
      LBE2=LBE
      IF(LBE.GT.GA)LBE2=GA
      DO 6555 BE=FBE,LBE2
      VUGABE=UOFF(V,U,1)+VADD(GA,BE)
      beu=itr(be)+u
      bei=itr(be)+i
c     IDX1=X3OFF(U,V,I)+X3ADD(BE,B,A)
           k1=max0(beu,ai,bv)
           k3=min0(beu,ai,bv)
           k3s=trsym(k3)
           k2=beu+ai+bv-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(U,V,I)+X3ADD(A,B,BE)
           k1=max0(bei,au,bv)
           k3=min0(bei,au,bv)
           k3s=trsym(k3)
           k2=bei+au+bv-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t2(VUGABE)=t2(VUGABE)-VAL1* t3o(IDX1) +VAL2* t3o(IDX2)
 6555 CONTINUE
C
      FBE2=FBE
      IF(FBE.LT.GA)FBE2=GA
      DO 6557 BE=FBE2,LBE
      UVBEGA=UOFF(U,V,1)+VADD(BE,GA)
      beu=itr(be)+u
      bei=itr(be)+i
c     IDX1=X3OFF(U,V,I)+X3ADD(BE,B,A)
           k1=max0(beu,ai,bv)
           k3=min0(beu,ai,bv)
           k3s=trsym(k3)
           k2=beu+ai+bv-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(U,V,I)+X3ADD(A,B,BE)
           k1=max0(bei,au,bv)
           k3=min0(bei,au,bv)
           k3s=trsym(k3)
           k2=bei+au+bv-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t2(UVBEGA)=t2(UVBEGA)-VAL1* t3o(IDX1) +VAL2* t3o(IDX2)
 6557 CONTINUE
C
 6560 CONTINUE
 6570 CONTINUE
 6580 CONTINUE
 6590 CONTINUE
 6595 CONTINUE
 6598 CONTINUE
c
      return  
      END
