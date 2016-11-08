      SUBROUTINE fccsdt(IPQ,ITR,itv,ain,bin,cin,din,ein,fin,
     .                  aoff,aadd,boff,badd,coff,cadd,
     .                  doff,dadd,eoff,eadd,foff,fadd,
     .                  NO,NV,NT,NT2,NDIMT2,ndimt3,nirred,level,
     .                  t1o,t1,t2o,t2,tau,t3o,t3,
     .                  ORBSYM,FLOV,UOFF,VADD,ZLX,FZO,FZV,FOCK,
     .                  F1,F2,F3,F4,F5,F6,
     .                  CHI,CHI1,CHI2,CHIO,CHIV,FFF,
     .                  i1off,i2off,i3off,trsym)
      IMPLICIT INTEGER (A-Z)
      REAL*8 VAL,VAL1,VAL2,VAL3,xfac
      real*8 t3o(*),t3(*),
     .       ain(*),bin(*),cin(*),din(*),ein(*),fin(*)  
c     INTEGER X3OFF(5,5,5),X3ADD(9,9,9) 
      REAL*8 CHI1(NV,NV,NO,NV),CHI2(NV,NO,NO,NO),CHI(NT,NT),
     .       F1(NV,NO,NV,NV),F2(NV,NO,NO,NO),F3(NV,NO,NO,NV),
     .       F4(NV,NO,NV,NO),F5(NO,NV,NV,NO),F6(NO,NV,NO,NV),
     .       CHIO(NO,NO),CHIV(NV,NV),FFF(NO,NV) 
      real*8 t2(NDIMT2),t2o(NDIMT2),TAU(NDIMT2),T1O(NO,NV),t1(NO,NV),
     .       fock(*)
      integer ORBSYM(NT),FLOV(NIRRED,4),
     .        UOFF(NO,NO,2),VADD(NV,NV),ZLX(NV,NV),
     .        IPQ(NT2),ITR(NV),FZO(NO),FZV(NV),itv(*),
     .        i1off(*),i2off(no*nv,nirred),i3off(*),trsym(*),
     .        aoff(*),aadd(*),boff(*),badd(*),coff(*),cadd(*),
     .        doff(*),dadd(*),eoff(*),eadd(*),foff(*),fadd(*)
C
c     call x3ARR(ORBSYM,FLOV,NIRRED,X3OFF,X3ADD,no)
c
c >>> read integrals from 78
c
c     call zero(ii,6000)
c     call rd78(ii)
c     WRITE(6,*)'integrals read '                 
c     WRITE(6,'(5d14.7)')(ii (IDX),IDX=1,100   )
c
c >>> initialize t3
c
c     call zero(t3,ndimt3)
C
C >>> form trsym
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
c >>> form tau !!
c
      DO 1290 U=1,NO
      IF(FZO(U).EQ.1)GO TO 1290
      USYM=ORBSYM(U)
      DO 1285 V=1,NO
      IF(FZO(V).EQ.1)GO TO 1285
      VSYM=ORBSYM(V)
      UVSYM=IEOR(USYM,VSYM)
      DO 1280 BE=1,NV
      IF(FZV(BE).EQ.1)GO TO 1280
      BESYM=ORBSYM(BE+NO)
      GASYM=IEOR(UVSYM,BESYM)+1
      FGA=FLOV(GASYM,3)-NO
      LGA=FLOV(GASYM,4)-NO
C
      LGA2=LGA
      IF(LGA.GT.BE)LGA2=BE
      DO 1276 GA=FGA,LGA2
      UVBEGA=UOFF(U,V,1)+VADD(BE,GA)
      tau(UVBEGA)=t2o(UVBEGA)+t1o(u,be)*t1o(v,ga)
 1276 CONTINUE
C
      FGA2=FGA
      IF(FGA.LT.BE)FGA2=BE
      DO 1277 GA=FGA2,LGA
      VUGABE=UOFF(V,U,1)+VADD(GA,BE)
      tau(VUGABE)=t2o(VUGABE)+t1o(u,be)*t1o(v,ga)
 1277 CONTINUE
C
 1280 CONTINUE
 1285 CONTINUE
 1290 CONTINUE
C
C >>> ADD T3 TERMS TO T1 AND T2 (do not divide by eig)
C
       call t3tots(ITR,itv,din,ein,fin,
     .             doff,dadd,eoff,eadd,foff,fadd,
     .             NO,NV,NT,NT2,NDIMT2,ndimt3,nirred,level,
     .             t1o,t1,t2o,t2,tau,t3o,t3,
     .             ORBSYM,FLOV,UOFF,VADD,ZLX,FZO,FZV,FOCK,
     .             i1off,i2off,i3off,trsym)
C
C >>> GO INTO CCSDT 
C
c  this routine splitted only because of compiler optimization
c
            call fsdt1 (IPQ,ITR,itv,ain,bin,cin,din,ein,fin,
     .                  aoff,aadd,boff,badd,coff,cadd,
     .                  doff,dadd,eoff,eadd,foff,fadd,
     .                  NO,NV,NT,NT2,NDIMT2,ndimt3,nirred,level,
     .                  t1o,t1,t2o,t2,tau,t3o,t3,
     .                  ORBSYM,FLOV,UOFF,VADD,ZLX,FZO,FZV,FOCK,
     .                  F1,F2,F3,F4,F5,F6,
     .                  CHI,CHI1,CHI2,CHIO,CHIV,FFF,
     .                  i1off,i2off,i3off,trsym)
c
c     call wwritw(90,t3,intowp(ndimt3),1,iend)
c     call zero(t3,ndimt3)
c     
C
C     GO FOR THE REST OF CCSDT
C
      IF(LEVEL.GT.3)THEN
C
C     >>CHIO<<
C
      CALL ZERO(CHIO,NO*NO)
      DO 4090 I=1,NO
      IF(FZO(I).EQ.1)GO TO 4090
      ISYM=ORBSYM(I)
      DO 4080 M=1,NO
      IF(FZO(M).EQ.1)GO TO 4080
      MSYM=ORBSYM(M)
      IF(MSYM.NE.ISYM)GO TO 4080
      IM=itr(MAX0(I,M))+MIN0(I,M)
CJN   IF(I.NE.M)  CHIO(I,M)=CHIO(I,M)+EIG(I)
      fim=ipq(max0(i,m))+min0(i,m)
      if(i.ne.m)chio(i,m)=chio(i,m)+fock(fim)
C
      IF(LEVEL.GT.4)THEN
      fe=flov(isym+1,3)-no
      le=flov(isym+1,4)-no
      do 4071 e=fe,le
      fme=ipq(e+no)+m
      chio(i,m)=chio(i,m)+fock(fme)*t1(i,e)
 4071 continue
      VAL=0.0D0
      DO 4072 N=1,NO
      IF(FZO(N).EQ.1)GO TO 4072
      IN=itr(MAX0(I,N))+MIN0(I,N)
      NSYM=ORBSYM(N)
      ESYM=NSYM
      FE=FLOV(ESYM+1,3)-NO
      LE=FLOV(ESYM+1,4)-NO
      DO 4070 E=FE,LE
      EM=itr(E)+M
      EN=itr(E)+N
      IMEN=eoff(EN)+eadd(IM)
      INEM=eoff(EM)+eadd(IN)
      VAL=VAL+(ein(IMEN)+ein(IMEN)-ein(INEM))*T1O(N,E)
 4070 CONTINUE
 4072 CONTINUE
      CHIO(I,M)=CHIO(I,M)+VAL
C
      VAL=0.0D0
      DO 4079 N=1,NO
      IF(FZO(N).EQ.1)GO TO 4079
      NSYM=ORBSYM(N)
      MNSYM=IEOR(MSYM,NSYM)
      DO 4077 E=1,NV
      IF(FZV(E).EQ.1)GO TO 4077
      ME=itr(E)+M
      NE=itr(E)+N
      ESYM=ORBSYM(E+NO)
      FSYM=IEOR(MNSYM,ESYM)+1
      FF=FLOV(FSYM,3)-NO
      LF=FLOV(FSYM,4)-NO
      DO 4075 F=FF,LF
      MF=itr(F)+M
      NF=itr(F)+N
      MENF=doff(MAX0(ME,NF))+dadd(MIN0(ME,NF))
      MFNE=doff(MAX0(MF,NE))+dadd(MIN0(MF,NE))
      ZL=ZLX(E,F)
      INEF=UOFF(I,N,ZL)+VADD(E,F)
      VAL=VAL+(din(MENF)+din(MENF)-din(MFNE))*TAU(INEF)
 4075 CONTINUE
 4077 CONTINUE
 4079 CONTINUE
      CHIO(I,M)=CHIO(I,M)+VAL
C
      ENDIF
c
 4080 CONTINUE
 4090 CONTINUE
C
C     >>CHIV<<
C
      CALL ZERO(CHIV,NV*NV)
      DO 4190 A=1,NV
      IF(FZV(A).EQ.1)GO TO 4190
      ASYM=ORBSYM(A+NO)
      ESYM=ASYM
      FE=FLOV(ESYM+1,3)-NO
      LE=FLOV(ESYM+1,4)-NO
      DO 4180 E=FE,LE
      AE=itv(MAX0(A,E))+MIN0(A,E)
CJN   IF(A.EQ.E)  CHIV(A,E)=CHIV(A,E)+EIG(A+NO)
      fae=ipq(max0(a,e)+no)+min0(a,e)+no
      if(a.ne.e)chiv(a,e)=chiv(a,e)+fock(fae)
C
      IF(LEVEL.GT.4)THEN
      fm=flov(asym+1,1)
      lm=flov(asym+1,2)
      do 4171 m=fm,lm
      fme=ipq(e+no)+m
      chiv(a,e)=chiv(a,e)-fock(fme)*t1(m,a)
 4171 continue
      VAL=0.0D0
      DO 4172 M=1,NO
      IF(FZO(M).EQ.1)GO TO 4172
      ME=itr(E)+M
      MSYM=ORBSYM(M)
      FSYM=MSYM
      FF=FLOV(FSYM+1,3)-NO
      LF=FLOV(FSYM+1,4)-NO
      DO 4170 F=FF,LF
      MF=itr(F)+M
      AF=itv(MAX0(A,F))+MIN0(A,F)
      AEMF=foff(AE)+fadd(MF)
      AFME=foff(AF)+fadd(ME)
      VAL=VAL+(fin(AEMF)+fin(AEMF)-fin(AFME))*T1O(M,F)
 4170 CONTINUE
 4172 CONTINUE
      CHIV(A,E)=CHIV(A,E)+VAL
C
      VAL=0.0D0
      DO 4179 N=1,NO
      IF(FZO(N).EQ.1)GO TO 4179
      NE=itr(E)+N
      NSYM=ORBSYM(N)
      DO 4177 M=1,NO
      IF(FZO(M).EQ.1)GO TO 4177
      ME=itr(E)+M
      MSYM=ORBSYM(M)
      MNSYM=IEOR(MSYM,NSYM)
      FSYM=IEOR(MNSYM,ESYM)+1
      FF=FLOV(FSYM,3)-NO
      LF=FLOV(FSYM,4)-NO
      DO 4175 F=FF,LF
      MF=itr(F)+M
      NF=itr(F)+N
      MENF=doff(MAX0(ME,NF))+dadd(MIN0(ME,NF))
      MFNE=doff(MAX0(MF,NE))+dadd(MIN0(MF,NE))
      ZL=ZLX(A,F)
      MNAF=UOFF(M,N,ZL)+VADD(A,F)
      VAL=VAL+(din(MENF)+din(MENF)-din(MFNE))*TAU(MNAF)
 4175 CONTINUE
 4177 CONTINUE
 4179 CONTINUE
      CHIV(A,E)=CHIV(A,E)-VAL
C
      ENDIF
c
 4180 CONTINUE
 4190 CONTINUE
C
C     >>CHI3<<    O O O O
C
      CALL ZERO(CHI,NT*NT)
      DO 4291 J=1,NO
      IF(FZO(J).EQ.1)GO TO 4291
      JSYM=ORBSYM(J)
      DO 4286 K=1,NO
      IF(FZO(K).EQ.1)GO TO 4286
      KSYM=ORBSYM(K)
      JKSYM=IEOR(JSYM,KSYM)
      DO 4231 M=1,NO
      IF(FZO(M).EQ.1)GO TO 4231
      JM=itr(MAX0(J,M))+MIN0(J,M)
      MSYM=ORBSYM(M)
      NSYM=IEOR(JKSYM,MSYM)+1
      FN=FLOV(NSYM,1)
      LN=FLOV(NSYM,2)
      DO 4226 N=FN,LN
      KN=itr(MAX0(K,N))+MIN0(K,N)
      JMKN=aoff(MAX0(JM,KN))+aadd(MIN0(JM,KN))
      CHI(M,N)=ain(JMKN)
C
      IF(LEVEL.GT.4)THEN
      VAL=0.0D0
      DO 4217 E=1,NV
      IF(FZV(E).EQ.1)GO TO 4217
      ME=itr(E)+M
      ESYM=ORBSYM(E+NO)
      FSYM=IEOR(JKSYM,ESYM)+1
      FF=FLOV(FSYM,3)-NO
      LF=FLOV(FSYM,4)-NO
      DO 4216 F=FF,LF
      NF=itr(F)+N
      MENF=doff(MAX0(ME,NF))+dadd(MIN0(ME,NF))
      ZL=ZLX(E,F)
      JKEF=UOFF(J,K,ZL)+VADD(E,F)
      VAL=VAL+din(MENF)*TAU(JKEF)
 4216 CONTINUE
 4217 CONTINUE
      CHI(M,N)=CHI(M,N)+VAL
C
      VAL=0.0D0
      DO 4221 E=1,NV
      IF(FZV(E).EQ.1)GO TO 4221
      EN=itr(E)+N
      JMEN=eoff(EN)+eadd(JM)
      EM=itr(E)+M
      EMKN=eoff(EM)+eadd(KN)
      VAL=VAL+ein(JMEN)*T1O(K,E)+ein(EMKN)*T1O(J,E)
 4221 CONTINUE
      CHI(M,N)=CHI(M,N)+VAL
C
      ENDIF
c
 4226 CONTINUE
 4231 CONTINUE
C
      DO 4280 I=1,NO
      IF(FZO(I).EQ.1)GO TO 4280
      ISYM=ORBSYM(I)
      IJKSYM=IEOR(ISYM,JKSYM)
      DO 4265 A=1,NV
      IF(FZV(A).EQ.1)GO TO 4265
      ASYM=ORBSYM(A+NO)
      ai=itr(a)+i
      DO 4260 B=1,NV
      IF(FZV(B).EQ.1)GO TO 4260
      bj=itr(b)+j
      BSYM=ORBSYM(B+NO)
      ABSYM=IEOR(ASYM,BSYM)
      CSYM=IEOR(ABSYM,IJKSYM)+1
      FC=FLOV(CSYM,3)-NO
      LC=FLOV(CSYM,4)-NO
      DO 4255 C=FC,LC
      VAL=0.0D0
      ck=itr(c)+k
      iflag=0
      if(ai.ge.bj.and.bj.ge.ck)iflag=1
      if(bj.ge.ai.and.ai.ge.ck)iflag=2
      if(ck.ge.bj.and.bj.ge.ai)iflag=3
      if(iflag.eq.0)go to 4255
      xfac=1.0d0
      if(ai.eq.bj.and.iflag.ne.3)xfac=2.0d0
      if(ai.eq.ck.and.iflag.ne.2)xfac=2.0d0
      if(ai.eq.bj.and.ai.eq.ck)xfac=3.0d0
      DO 4240 M=1,NO
      IF(FZO(M).EQ.1)GO TO 4240
      bm=itr(b)+m
      MSYM=ORBSYM(M)
      NSYM=IEOR(MSYM,JKSYM)
      FN=FLOV(NSYM+1,1)
      LN=FLOV(NSYM+1,2)
      DO 4235 N=FN,LN
      cn=itr(c)+n
c     IDX=X3OFF(I,M,N)+X3ADD(A,B,C)
           k1=max0(bm,ai,cn)
           k3=min0(bm,ai,cn)
           k3s=trsym(k3)
           k2=bm+ai+cn-k1-k3
           idx=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      VAL=VAL+CHI(M,N)*t3o(IDX)
 4235 CONTINUE
 4240 CONTINUE
c     IDX1=X3OFF(I,J,K)+X3ADD(A,B,C)
           k1=max0(ai,bj,ck)
           k3=min0(ai,bj,ck)
           k3s=trsym(k3)
           k2=ai+bj+ck-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     t3(IDX1)=t3(IDX1)+VAL*0.5d0
      t3(IDX1)=t3(IDX1)+VAL*xfac
 4255 CONTINUE
 4260 CONTINUE
 4265 CONTINUE
 4280 CONTINUE
 4286 CONTINUE
 4291 CONTINUE
C
C     >>CHI4<<    V V V V
C
      CALL ZERO(CHI,NT*NT)
      DO 4390 B=1,NV
      IF(FZV(B).EQ.1)GO TO 4390
      BSYM=ORBSYM(B+NO)
      DO 4385 C=1,NV
      IF(FZV(C).EQ.1)GO TO 4385
      CSYM=ORBSYM(C+NO)
      BCSYM=IEOR(BSYM,CSYM)
      DO 4340 E=1,NV
      IF(FZV(E).EQ.1)GO TO 4340
      BE=itv(MAX0(B,E))+MIN0(B,E)
      ESYM=ORBSYM(E+NO)
      FSYM=IEOR(BCSYM,ESYM)+1
      FF=FLOV(FSYM,3)-NO
      LF=FLOV(FSYM,4)-NO
      DO 4335 F=FF,LF
      CF=itv(MAX0(C,F))+MIN0(C,F)
      BECF=boff(MAX0(BE,CF))+badd(MIN0(BE,CF))
      CHI(E,F)=bin(BECF)
C
      IF(LEVEL.GT.4)THEN
      VAL=0.0D0
      DO 4317 M=1,NO
      IF(FZO(M).EQ.1)GO TO 4317
      ME=itr(E)+M
      MSYM=ORBSYM(M)
      NSYM=IEOR(BCSYM,MSYM)+1
      FN=FLOV(NSYM,1)
      LN=FLOV(NSYM,2)
      DO 4315 N=FN,LN
      NF=itr(F)+N
      MENF=doff(MAX0(ME,NF))+dadd(MIN0(ME,NF))
      ZL=ZLX(B,C)
      MNBC=UOFF(M,N,ZL)+VADD(B,C)
      VAL=VAL+din(MENF)*TAU(MNBC)
 4315 CONTINUE
 4317 CONTINUE
      CHI(E,F)=CHI(E,F)+VAL
C
      VAL=0.0D0
      fm=flov(csym+1,1)
      lm=flov(csym+1,2)
      DO 4320 M=fm,lm
      MF=itr(F)+M
      BEMF=foff(BE)+fadd(MF)
      VAL=VAL+fin(BEMF)*T1O(M,C)
 4320 CONTINUE
      CHI(E,F)=CHI(E,F)-VAL
C
      VAL=0.0D0
      fm=flov(bsym+1,1)
      lm=flov(bsym+1,2)
      DO 4325 M=fm,lm
      ME=itr(E)+M
      MECF=foff(CF)+fadd(ME)
      VAL=VAL+fin(MECF)*T1O(M,B)
 4325 CONTINUE
      CHI(E,F)=CHI(E,F)-VAL
C
      ENDIF
 4335 CONTINUE
 4340 CONTINUE
C
      DO 4380 A=1,NV
      IF(FZV(A).EQ.1)GO TO 4380
      ASYM=ORBSYM(A+NO)
      ABCSYM=IEOR(ASYM,BCSYM)
      DO 4365 I=1,NO
      IF(FZO(I).EQ.1)GO TO 4365
      ai=itr(a)+i
      ISYM=ORBSYM(I)
      DO 4360 J=1,NO
      IF(FZO(J).EQ.1)GO TO 4360
      bj=itr(b)+j
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      KSYM=IEOR(IJSYM,ABCSYM)+1
      FK=FLOV(KSYM,1)
      LK=FLOV(KSYM,2)
      DO 4355 K=FK,LK
      VAL=0.0D0
      ck=itr(c)+k
      iflag=0
      if(ai.ge.bj.and.bj.ge.ck)iflag=1
      if(bj.ge.ai.and.ai.ge.ck)iflag=2
      if(ck.ge.bj.and.bj.ge.ai)iflag=3
      if(iflag.eq.0)go to 4355
      xfac=1.0d0
      if(ai.eq.bj.and.iflag.ne.3)xfac=2.0d0
      if(ai.eq.ck.and.iflag.ne.2)xfac=2.0d0
      if(ai.eq.bj.and.ai.eq.ck)xfac=3.0d0
      DO 4350 E=1,NV
      IF(FZV(E).EQ.1)GO TO 4350
      ej=itr(e)+j
      ESYM=ORBSYM(E+NO)
      FSYM=IEOR(ESYM,BCSYM)
      IF=FLOV(FSYM+1,3)-NO
      LF=FLOV(FSYM+1,4)-NO
      DO 4348 F=IF,LF
      fk=itr(f)+k
c     IDX=X3OFF(I,J,K)+X3ADD(A,E,F)
           k1=max0(ej,ai,fk)
           k3=min0(ej,ai,fk)
           k3s=trsym(k3)
           k2=ej+ai+fk-k1-k3
           idx=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      VAL=VAL+CHI(E,F)*t3o(IDX)
 4348 CONTINUE
 4350 CONTINUE
c     IDX1=X3OFF(I,J,K)+X3ADD(A,B,C)
           k1=max0(ai,bj,ck)
           k3=min0(ai,bj,ck)
           k3s=trsym(k3)
           k2=ai+bj+ck-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      t3(IDX1)=t3(IDX1)+VAL *xfac
c     t3(IDX1)=t3(IDX1)+VAL *0.5d0
 4355 CONTINUE
 4360 CONTINUE
 4365 CONTINUE
 4380 CONTINUE
 4385 CONTINUE
 4390 CONTINUE
C
C     >>CHI5<<    V O V O
C
      CALL ZERO(CHI,NT*NT)
      DO 4490 A=1,NV
      IF(FZV(A).EQ.1)GO TO 4490
      ASYM=ORBSYM(A+NO)
      DO 4485 I=1,NO
      IF(FZO(I).EQ.1)GO TO 4485
      ISYM=ORBSYM(I)
      ai=itr(a)+i
      AISYM=IEOR(ASYM,ISYM)
      DO 4435 M=1,NO
      IF(FZO(M).EQ.1)GO TO 4435
      MSYM=ORBSYM(M)
      MI=itr(MAX0(M,I))+MIN0(M,I)
      ESYM=IEOR(AISYM,MSYM)
      FE=FLOV(ESYM+1,3)-NO
      LE=FLOV(ESYM+1,4)-NO
      DO 4430 E=FE,LE
      AE=itv(MAX0(A,E))+MIN0(A,E)
      AEMI=coff(AE)+cadd(MI)
      CHI(M,E)=cin(AEMI)
C
      IF(LEVEL.GT.4)THEN
      VAL=0.0D0
      MESYM=IEOR(MSYM,ESYM)
      DO 4417 N=1,NO
      IF(FZO(N).EQ.1)GO TO 4417
      NE=itr(E)+N
      NSYM=ORBSYM(N)
      FSYM=IEOR(MESYM,NSYM)+1
      FF=FLOV(FSYM,3)-NO
      LF=FLOV(FSYM,4)-NO
      DO 4415 F=FF,LF
      MF=itr(F)+M
      MFNE=doff(MAX0(MF,NE))+dadd(MIN0(MF,NE))
      ZL=ZLX(F,A)
      INFA=UOFF(I,N,ZL)+VADD(F,A)
      VAL=VAL+din(MFNE)*TAU(INFA)
 4415 CONTINUE
 4417 CONTINUE
      CHI(M,E)=CHI(M,E)-VAL
C
      VAL=0.0D0
      fn=flov(asym+1,1)
      ln=flov(asym+1,2)
      DO 4420 N=fn,ln
      NE=itr(E)+N
      NEMI=eoff(NE)+eadd(MI)
      VAL=VAL+ein(NEMI)*T1O(N,A)
 4420 CONTINUE
      CHI(M,E)=CHI(M,E)-VAL
C
      VAL=0.0D0
      ff=flov(isym+1,3)-no
      lf=flov(isym+1,4)-no
      DO 4422 F=ff,lf
      MF=itr(F)+M
      AEMF=foff(AE)+fadd(MF)
      VAL=VAL+fin(AEMF)*T1O(I,F)
 4422 CONTINUE
      CHI(M,E)=CHI(M,E)+VAL
C
      ENDIF
 4430 CONTINUE
 4435 CONTINUE
C
      DO 4480 B=1,NV
      IF(FZV(B).EQ.1)GO TO 4480
      BSYM=ORBSYM(B+NO)
      AIBSYM=IEOR(BSYM,AISYM)
      DO 4465 J=1,NO
      IF(FZO(J).EQ.1)GO TO 4465
      JSYM=ORBSYM(J)
      bj=itr(b)+j
      DO 4460 C=1,NV
      IF(FZV(C).EQ.1)GO TO 4460
      CSYM=ORBSYM(C+NO)
      JCSYM=IEOR(CSYM,JSYM)
      KSYM=IEOR(JCSYM,AIBSYM)+1
      FK=FLOV(KSYM,1)
      LK=FLOV(KSYM,2)
      DO 4455 K=FK,LK
      VAL1=0.0D0
c     VAL2=0.0D0
c     VAL3=0.0D0
      ck=itr(c)+k
      iflag=0
      if(ai.ge.bj.and.bj.ge.ck)iflag=1
      if(bj.ge.ai.and.ai.ge.ck)iflag=2
      if(ck.ge.bj.and.bj.ge.ai)iflag=3
      if(iflag.eq.0)go to 4455
      xfac=1.0d0
      if(ai.eq.bj.and.iflag.ne.3)xfac=2.0d0
      if(ai.eq.ck.and.iflag.ne.2)xfac=2.0d0
      if(ai.eq.bj.and.ai.eq.ck)xfac=3.0d0
      DO 4440 E=1,NV
      IF(FZV(E).EQ.1)GO TO 4440
      ESYM=ORBSYM(E+NO)
      ej=itr(e)+j
      ek=itr(e)+k
      MSYM=IEOR(ESYM,AISYM)
      FM=FLOV(MSYM+1,1)
      LM=FLOV(MSYM+1,2)
      DO 4438 M=FM,LM
      em=itr(e)+m
      bm=itr(b)+m
      cm=itr(c)+m
c     IDX1=X3OFF(M,J,K)+X3ADD(E,B,C)
           k1=max0(em,bj,ck)
           k3=min0(em,bj,ck)
           k3s=trsym(k3)
           k2=em+bj+ck-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(M,J,K)+X3ADD(B,E,C)
c          k1=max0(bm,ej,ck)
c          k3=min0(bm,ej,ck)
c          k3s=trsym(k3)
c          k2=bm+ej+ck-k1-k3
c          idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX3=X3OFF(M,J,K)+X3ADD(C,B,E)
c          k1=max0(cm,ek,bj)
c          k3=min0(cm,ek,bj)
c          k3s=trsym(k3)
c          k2=cm+ek+bj-k1-k3
c          idx3=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      VAL1=VAL1-CHI(M,E)* t3o(IDX1)
c     VAL2=VAL2-CHI(M,E)* t3o(IDX2)
c     VAL3=VAL3-CHI(M,E)* t3o(IDX3)
 4438 CONTINUE
 4440 CONTINUE
c
c     IDX1=X3OFF(I,J,K)+X3ADD(A,B,C)
           k1=max0(ai,bj,ck)
           k3=min0(ai,bj,ck)
           k3s=trsym(k3)
           k2=ai+bj+ck-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     t3(IDX1)=t3(IDX1)+VAL1 *0.5d0
      t3(IDX1)=t3(IDX1)+VAL1 *xfac
c
c     IDX1=X3OFF(I,J,K)+X3ADD(B,A,C)
c          bi=itr(b)+i
c          aj=itr(a)+j
c          k1=max0(bi,aj,ck)
c          k3=min0(bi,aj,ck)
c          k3s=trsym(k3)
c          k2=bi+aj+ck-k1-k3
c          idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     t3(IDX1)=t3(IDX1)+VAL2 
c
c     IDX1=X3OFF(I,J,K)+X3ADD(C,B,A)
c          ci =itr(c)+i
c          ak =itr(a)+k
c          k1=max0(ci,bj,ak)
c          k3=min0(ci,bj,ak)
c          k3s=trsym(k3)
c          k2=ci+bj+ak-k1-k3
c          idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     t3(IDX1)=t3(IDX1)+VAL3
c
 4455 CONTINUE
 4460 CONTINUE
 4465 CONTINUE
 4480 CONTINUE
c
      DO 5480 B=1,NV
      IF(FZV(B).EQ.1)GO TO 5480
      BSYM=ORBSYM(B+NO)
      AIBSYM=IEOR(BSYM,AISYM)
      DO 5465 J=1,NO
      IF(FZO(J).EQ.1)GO TO 5465
      JSYM=ORBSYM(J)
      bj=itr(b)+j
      aj=itr(a)+j
      bi=itr(b)+i
      DO 5460 C=1,NV
      IF(FZV(C).EQ.1)GO TO 5460
      CSYM=ORBSYM(C+NO)
      JCSYM=IEOR(CSYM,JSYM)
      KSYM=IEOR(JCSYM,AIBSYM)+1
      FK=FLOV(KSYM,1)
      LK=FLOV(KSYM,2)
      DO 5455 K=FK,LK
      VAL2=0.0D0
      ck=itr(c)+k
      iflag=0
      if(bi.ge.aj.and.aj.ge.ck)iflag=1
      if(aj.ge.bi.and.bi.ge.ck)iflag=2
      if(ck.ge.aj.and.aj.ge.bi)iflag=3
      if(iflag.eq.0)go to 5455
      xfac=1.0d0
      if(bi.eq.aj.and.iflag.ne.3)xfac=2.0d0
      if(bi.eq.ck.and.iflag.ne.2)xfac=2.0d0
      if(bi.eq.aj.and.bi.eq.ck)xfac=3.0d0
      DO 5440 E=1,NV
      IF(FZV(E).EQ.1)GO TO 5440
      ESYM=ORBSYM(E+NO)
      ej=itr(e)+j
      ek=itr(e)+k
      MSYM=IEOR(ESYM,AISYM)
      FM=FLOV(MSYM+1,1)
      LM=FLOV(MSYM+1,2)
      DO 5438 M=FM,LM
      em=itr(e)+m
      bm=itr(b)+m
      cm=itr(c)+m
c     IDX2=X3OFF(M,J,K)+X3ADD(B,E,C)
           k1=max0(bm,ej,ck)
           k3=min0(bm,ej,ck)
           k3s=trsym(k3)
           k2=bm+ej+ck-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      VAL2=VAL2-CHI(M,E)* t3o(IDX2)
 5438 CONTINUE
 5440 CONTINUE
c
c     IDX1=X3OFF(I,J,K)+X3ADD(B,A,C)
           k1=max0(bi,aj,ck)
           k3=min0(bi,aj,ck)
           k3s=trsym(k3)
           k2=bi+aj+ck-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     t3(IDX1)=t3(IDX1)+VAL2 *0.5d0
      t3(IDX1)=t3(IDX1)+VAL2*xfac
 5455 CONTINUE
 5460 CONTINUE
 5465 CONTINUE
 5480 CONTINUE
c

      DO 5680 B=1,NV
      IF(FZV(B).EQ.1)GO TO 5680
      BSYM=ORBSYM(B+NO)
      AIBSYM=IEOR(BSYM,AISYM)
      DO 5665 J=1,NO
      IF(FZO(J).EQ.1)GO TO 5665
      JSYM=ORBSYM(J)
      bj=itr(b)+j
      DO 5660 C=1,NV
      IF(FZV(C).EQ.1)GO TO 5660
      CSYM=ORBSYM(C+NO)
      JCSYM=IEOR(CSYM,JSYM)
      KSYM=IEOR(JCSYM,AIBSYM)+1
      FK=FLOV(KSYM,1)
      LK=FLOV(KSYM,2)
      ci =itr(c)+i
      DO 5655 K=FK,LK
      VAL3=0.0D0
      ck=itr(c)+k
      ak =itr(a)+k
      iflag=0
      if(ci.ge.bj.and.bj.ge.ak)iflag=1
      if(bj.ge.ci.and.ci.ge.ak)iflag=2
      if(ak.ge.bj.and.bj.ge.ci)iflag=3
      if(iflag.eq.0)go to 5655
      xfac=1.0d0
      if(ci.eq.bj.and.iflag.ne.3)xfac=2.0d0
      if(ci.eq.ak.and.iflag.ne.2)xfac=2.0d0
      if(ci.eq.bj.and.ci.eq.ak)xfac=3.0d0
      DO 5640 E=1,NV
      IF(FZV(E).EQ.1)GO TO 5640
      ESYM=ORBSYM(E+NO)
      ej=itr(e)+j
      ek=itr(e)+k
      MSYM=IEOR(ESYM,AISYM)
      FM=FLOV(MSYM+1,1)
      LM=FLOV(MSYM+1,2)
      DO 5638 M=FM,LM
      em=itr(e)+m
      bm=itr(b)+m
      cm=itr(c)+m
c     IDX3=X3OFF(M,J,K)+X3ADD(C,B,E)
           k1=max0(cm,ek,bj)
           k3=min0(cm,ek,bj)
           k3s=trsym(k3)
           k2=cm+ek+bj-k1-k3
           idx3=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      VAL3=VAL3-CHI(M,E)* t3o(IDX3)
 5638 CONTINUE
 5640 CONTINUE
c
c     IDX1=X3OFF(I,J,K)+X3ADD(C,B,A)
           k1=max0(ci,bj,ak)
           k3=min0(ci,bj,ak)
           k3s=trsym(k3)
           k2=ci+bj+ak-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     t3(IDX1)=t3(IDX1)+VAL3*0.5d0
      t3(IDX1)=t3(IDX1)+VAL3*xfac
 5655 CONTINUE
 5660 CONTINUE
 5665 CONTINUE
 5680 CONTINUE
c
 4485 CONTINUE
 4490 CONTINUE
C
C     >>CHI6<<    V O O V
C
      CALL ZERO(CHI,NT*NT)
      DO 4590 A=1,NV
      IF(FZV(A).EQ.1)GO TO 4590
      ASYM=ORBSYM(A+NO)
      DO 4585 I=1,NO
      IF(FZO(I).EQ.1)GO TO 4585
      ISYM=ORBSYM(I)
      AISYM=IEOR(ASYM,ISYM)
      AI=itr(A)+I
      DO 4530 E=1,NV
      IF(FZV(E).EQ.1)GO TO 4530
      ESYM=ORBSYM(E+NO)
      MSYM=IEOR(AISYM,ESYM)
      FM=FLOV(MSYM+1,1)
      LM=FLOV(MSYM+1,2)
      DO 4525 M=FM,LM
      ME=itr(E)+M
      AIME=doff(MAX0(AI,ME))+dadd(MIN0(AI,ME))
      CHI(M,E)=din(AIME)
C
      IF(LEVEL.GT.4)THEN
      ai=itr(a)+i
      VAL=0.0D0
      MESYM=IEOR(MSYM,ESYM)
      DO 4517 N=1,NO
      IF(FZO(N).EQ.1)GO TO 4517
      NE=itr(E)+N
      ME=itr(E)+M
      NSYM=ORBSYM(N)
      FSYM=IEOR(MESYM,NSYM)+1
      FF=FLOV(FSYM,3)-NO
      LF=FLOV(FSYM,4)-NO
      DO 4515 F=FF,LF
      MF=itr(F)+M
      NF=itr(F)+N
      MFNE=doff(MAX0(MF,NE))+dadd(MIN0(MF,NE))
      MENF=doff(MAX0(ME,NF))+dadd(MIN0(ME,NF))
      ZL=ZLX(F,A)
      INFA=UOFF(I,N,ZL)+VADD(F,A)
      INAF=UOFF(N,I,ZL)+VADD(F,A)
      VAL=VAL+(din(MENF)+din(MENF)-din(MFNE))*t2o(INAF)
     .                          -din(MENF) *TAU(INFA)
 4515 CONTINUE
 4517 CONTINUE
      CHI(M,E)=CHI(M,E)+VAL
C
      VAL=0.0D0
      fn=flov(asym+1,1)
      ln=flov(asym+1,2)
      DO 4520 N=fn,ln
      NI=itr(MAX0(N,I))+MIN0(N,I)
      NIME=eoff(ME)+eadd(NI)
      VAL=VAL+ein(NIME)*T1O(N,A)
 4520 CONTINUE
      CHI(M,E)=CHI(M,E)-VAL
C
      VAL=0.0D0
      EM=itr(E)+M
      ff=flov(isym+1,3)-no
      lf=flov(isym+1,4)-no
      DO 4522 F=ff,lf
      AF=itv(MAX0(A,F))+MIN0(A,F)
      AFEM=foff(AF)+fadd(EM)
      VAL=VAL+fin(AFEM)*T1O(I,F)
 4522 CONTINUE
      CHI(M,E)=CHI(M,E)+VAL
C
      ENDIF
c
 4525 CONTINUE
 4530 CONTINUE
C
      DO 4580 B=1,NV
      IF(FZV(B).EQ.1)GO TO 4580
      BSYM=ORBSYM(B+NO)
      AIBSYM=IEOR(BSYM,AISYM)
      DO 4565 J=1,NO
      IF(FZO(J).EQ.1)GO TO 4565
      JSYM=ORBSYM(J)
      bj=itr(b)+j
      DO 4560 C=1,NV
      IF(FZV(C).EQ.1)GO TO 4560
      CSYM=ORBSYM(C+NO)
      JCSYM=IEOR(CSYM,JSYM)
      KSYM=IEOR(JCSYM,AIBSYM)+1
      FK=FLOV(KSYM,1)
      LK=FLOV(KSYM,2)
      DO 4555 K=FK,LK
      VAL=0.0D0
      ck=itr(c)+k
      iflag=0
      if(ai.ge.bj.and.bj.ge.ck)iflag=1
      if(bj.ge.ai.and.ai.ge.ck)iflag=2
      if(ck.ge.bj.and.bj.ge.ai)iflag=3
      if(iflag.eq.0)go to 4555
      xfac=1.0d0
      if(ai.eq.bj.and.iflag.ne.3)xfac=2.0d0
      if(ai.eq.ck.and.iflag.ne.2)xfac=2.0d0
      if(ai.eq.bj.and.ai.eq.ck)xfac=3.0d0
      DO 4540 E=1,NV
      IF(FZV(E).EQ.1)GO TO 4540
      ESYM=ORBSYM(E+NO)
      ej=itr(e)+j
      ek=itr(e)+k
      MSYM=IEOR(ESYM,AISYM)
      FM=FLOV(MSYM+1,1)
      LM=FLOV(MSYM+1,2)
      DO 4538 M=FM,LM
      em=itr(e)+m
      bm=itr(b)+m
      cm=itr(c)+m
c     IDX1=X3OFF(M,J,K)+X3ADD(E,B,C)
           k1=max0(em,bj,ck)
           k3=min0(em,bj,ck)
           k3s=trsym(k3)
           k2=em+bj+ck-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(M,J,K)+X3ADD(B,E,C)
           k1=max0(bm,ej,ck)
           k3=min0(bm,ej,ck)
           k3s=trsym(k3)
           k2=bm+ej+ck-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX3=X3OFF(M,J,K)+X3ADD(C,B,E)
           k1=max0(cm,ek,bj)
           k3=min0(cm,ek,bj)
           k3s=trsym(k3)
           k2=cm+ek+bj-k1-k3
           idx3=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      VAL=VAL+CHI(M,E)*(t3o(IDX1)+t3o(IDX1)-t3o(IDX2)-t3o(IDX3))
 4538 CONTINUE
 4540 CONTINUE
c     IDX1=X3OFF(I,J,K)+X3ADD(A,B,C)
           k1=max0(ai,bj,ck)
           k3=min0(ai,bj,ck)
           k3s=trsym(k3)
           k2=ai+bj+ck-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     t3(IDX1)=t3(IDX1)+VAL *0.5d0
      t3(IDX1)=t3(IDX1)+VAL *xfac
C
 4555 CONTINUE
 4560 CONTINUE
 4565 CONTINUE
 4580 CONTINUE
 4585 CONTINUE
 4590 CONTINUE
C
C     ADD THE 2ND SET OF CONTRIBUTIONS
C
      DO 5590 I=1,NO
      IF(FZO(I).EQ.1)GO TO 5590
      ISYM=ORBSYM(I)
      DO 5580 J=1,NO
      IF(FZO(J).EQ.1)GO TO 5580
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      DO 5570 K=1,NO
      IF(FZO(K).EQ.1)GO TO 5570
      KSYM=ORBSYM(K)
      JKSYM=IEOR(JSYM,KSYM)
      IJKSYM=IEOR(IJSYM,KSYM)
      DO 5560 A=1,NV
      IF(FZV(A).EQ.1)GO TO 5560
      ai=itr(a)+i
      ASYM=ORBSYM(A+NO)
      AISYM=IEOR(ASYM,ISYM)
      DO 5550 B=1,NV
      IF(FZV(B).EQ.1)GO TO 5550
      BSYM=ORBSYM(B+NO)
      bj=itr(b)+j
      BISYM=IEOR(BSYM,ISYM)
      ABSYM=IEOR(ASYM,BSYM)
      CSYM=IEOR(IJKSYM,ABSYM)
      BCSYM=IEOR(BSYM,CSYM)
      CISYM=IEOR(CSYM,ISYM)
      FC=FLOV(CSYM+1,3)-NO
      LC=FLOV(CSYM+1,4)-NO
      DO 5540 C=FC,LC
      VAL=0.0D0
      ck=itr(c)+k
      iflag=0
      if(ai.ge.bj.and.bj.ge.ck)iflag=1
      if(bj.ge.ai.and.ai.ge.ck)iflag=2
      if(ck.ge.bj.and.bj.ge.ai)iflag=3
      if(iflag.eq.0)go to 5540
      xfac=1.0d0
      if(ai.eq.bj.and.iflag.ne.3)xfac=2.0d0
      if(ai.eq.ck.and.iflag.ne.2)xfac=2.0d0
      if(ai.eq.bj.and.ai.eq.ck)xfac=3.0d0
C
      FE=FLOV(ASYM+1,3)-NO
      LE=FLOV(ASYM+1,4)-NO
      DO 5000 E=FE,LE
      ei=itr(e)+i
c     IDX=X3OFF(I,J,K)+X3ADD(E,B,C)
           k1=max0(ck,ei,bj)
           k3=min0(ck,ei,bj)
           k3s=trsym(k3)
           k2=ck+ei+bj-k1-k3
           idx=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      VAL=VAL+CHIV(A,E)*t3o(IDX)
 5000 CONTINUE
C
      FM=FLOV(ISYM+1,1)
      LM=FLOV(ISYM+1,2)
      DO 5100 M=FM,LM
      am=itr(a)+m
c     IDX=X3OFF(M,J,K)+X3ADD(A,B,C)
           k1=max0(ck,bj,am)
           k3=min0(ck,bj,am)
           k3s=trsym(k3)
           k2=ck+bj+am-k1-k3
           idx=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      VAL=VAL-CHIO(I,M)*t3o(IDX)
 5100 CONTINUE
C
c     IDX1=X3OFF(I,J,K)+X3ADD(A,B,C)
           k1=max0(ai,bj,ck)
           k3=min0(ai,bj,ck)
           k3s=trsym(k3)
           k2=ai+bj+ck-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c          t3(IDX1)=t3(IDX1)+VAL *0.5d0
           t3(IDX1)=t3(IDX1)+VAL *xfac
 5540 CONTINUE
 5550 CONTINUE
 5560 CONTINUE
 5570 CONTINUE
 5580 CONTINUE
 5590 CONTINUE
C
C
C     ENDING LEVEL.GT.3 (2ND SET OF CONTRIBUTIONS)
      ENDIF
C
c   renormalize triple and double diagonal elements
c
c     DO 3690 A=1,NV
c        IF(FZV(A).EQ.1)GO TO 3690
c        ASYM=ORBSYM(A+NO)
c        DO 3685 I=1,NO
c           IF(FZO(I).EQ.1)GO TO 3685
c           ISYM=ORBSYM(I)
c           AI=ITR(A)+I
c           if(asym.ne.isym)go to 3685
c           aiaiai=i1off(ai)+i2off(ai,1)+i3off(ai)
c           t3(aiaiai)=t3(aiaiai)*3.0d0
c3685    CONTINUE
c3690 CONTINUE
c
c     DO 3790 A=1,NV
c        IF(FZV(A).EQ.1)GO TO 3790
c        ASYM=ORBSYM(A+NO)
c        DO 3785 I=1,NO
c           IF(FZO(I).EQ.1)GO TO 3785
c           ISYM=ORBSYM(I)
c           AI=ITR(A)+I
c           ais=ieor(asym,isym)+1
c           do 3780 b=1,nv
c              IF(FZV(b).EQ.1)GO TO 3780
c              bSYM=ORBSYM(b+NO)
c              DO 3775 j=1,NO
c                 IF(FZO(j).EQ.1)GO TO 3775
c                 jSYM=ORBSYM(j)
c                 if(bsym.ne.jsym)go to 3775
c                 bj=ITR(b)+j
c                 if(ai.eq.bj)go to 3775
c                 if(ai.gt.bj)then
c                    aiaibj=i1off(ai)+i2off(ai,1)+i3off(bj)
c                    t3(aiaibj)=t3(aiaibj)*2.0d0
c                 else
c                    bjaiai=i1off(bj)+i2off(ai,ais)+i3off(ai)
c                    t3(bjaiai)=t3(bjaiai)*2.0d0
c                 endif
c3775          CONTINUE
c3780       CONTINUE
c3785    CONTINUE
c3790 CONTINUE
c
c     call wreadw(90,t3o,intowp(ndimt3),1,iend)
c
C     DIVIDE BY MO ENERGIES
C
      DO 3590 I=1,NO
      IF(FZO(I).EQ.1)GO TO 3590
      ISYM=ORBSYM(I)
      IY=IPQ(I)+I
      DO 3580 J=1,NO
      IF(FZO(J).EQ.1)GO TO 3580
      JSYM=ORBSYM(J)
      IJSYM=IEOR(ISYM,JSYM)
      JJ=IPQ(J)+J
      DO 3570 K=1,NO
      IF(FZO(K).EQ.1)GO TO 3570
      KSYM=ORBSYM(K)
      IJKSYM=IEOR(IJSYM,KSYM)
      KK=IPQ(K)+K
      DO 3560 A=1,NV
      IF(FZV(A).EQ.1)GO TO 3560
      ASYM=ORBSYM(A+NO)
      AA=IPQ(A+NO)+A+NO
      ai=itr(a)+i
      DO 3550 B=1,NV
      IF(FZV(B).EQ.1)GO TO 3550
      bj=itr(b)+j
      if(bj.gt.ai)go to 3550
      BSYM=ORBSYM(B+NO)
      ABSYM=IEOR(ASYM,BSYM)
      BB=IPQ(B+NO)+B+NO
      CSYM=IEOR(IJKSYM,ABSYM)+1
      FC=FLOV(CSYM,3)-NO
      LC=FLOV(CSYM,4)-NO
      DO 3540 C=FC,LC
      CC=IPQ(C+NO)+C+NO
      ck=itr(c)+k
      if(ck.gt.bj)go to 3540
      VAL=FOCK(IY)+FOCK(JJ)+FOCK(KK)-FOCK(AA)-FOCK(BB)-FOCK(CC)
c     IDX=X3OFF(I,J,K)+X3ADD(A,B,C)
           k1=max0(ai,bj,ck)
           k3=min0(ai,bj,ck)
           k3s=trsym(k3)
           k2=ai+bj+ck-k1-k3
           idx=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     t3(IDX)=(t3(IDX)+t3o(idx))/VAL
      t3(IDX)=t3(IDX)/VAL
 3540 CONTINUE
 3550 CONTINUE
 3560 CONTINUE
 3570 CONTINUE
 3580 CONTINUE
 3590 CONTINUE
C
c    transfer t3 back into t3o
c
c     do 3890 i=1,ndimt3
c        t3o(i)=t3(i)
c3890 continue
c 
c     CALL SWRIT(ITAP69,NDIMT3,1)
c     CALL SWRIT(ITAP69,t3,INTOWP(NDIMT3))
C
c     CALL TDIFF(S,t1,T2,T1O,T1NORM,DELT,NO,NV,NDIMT2,
c    .           ORBSYM,FLOV,NIRRED,UOFF,VADD,ZLX,FZO,FZV,
c    .           t3,t3o,NDIMT3,X3OFF,X3ADD)
C
c     CALL TWRIT(t1,S,ORBSYM,FLOV,NIRRED,UOFF,VADD,IVAL,TVAL,
c    .           ZLX,FZO,FZV,NO,NV,NIRRED,ITYP,NXT1,NXT2,NCSF,
c    .           t3,X3OFF,X3ADD,NXT3)
c
c     WRITE(6,*)' T3 AMPLITUDES new code'
c     WRITE(6,'(5d14.7)')(t3(IDX),IDX=1,100   )
c     WRITE(6,*)'FINAL T2 AMPLITUDES'
c     WRITE(6,'(10F7.4)')( t2(IDX),IDX=1,NDIMT2)
c
c     CALL DIISD3(T1O,t1,T2,S,NO,NV,NDIMT2,IT,MAXDIM,
c    .            CCC,BBB,ITAP98,ITAP99,ITC,BB2,
c    .            t3,t3o,NDIMT3)
c
      return  
      END
