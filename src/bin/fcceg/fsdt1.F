      SUBROUTINE fsdt1(IPQ,ITR,itv,ain,bin,cin,din,ein,fin,
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
C >>> GO INTO CCSDT 
C
      CALL ZERO (FFF,NO*NV)
      CALL ZERO (F1,NV*NO*NV*NV)
      CALL ZERO (F2,NV*NO*NO*NO)
      CALL ZERO (F3,NV*NO*NO*NV)
      CALL ZERO (F4,NV*NO*NV*NO)
      CALL ZERO (F5,NO*NV*NV*NO)
      CALL ZERO (F6,NO*NV*NO*NV)
      IF(LEVEL.GT.1)THEN
c
C    WATCH OUT! THERE IS A NON-SCF TERM FOR CCSDT2, NOT CODED.
c  term has now been coded. may23, 1990.
c
c     IF(LEVEL.GT.2)THEN
      DO 2085 M=1,NO
      IF(FZO(M).EQ.1)GO TO 2085
      msym=orbsym(m)
      fe=flov(msym+1,3)-no
      le=flov(msym+1,4)-no
      DO 2080 E=fe,le
      ME=itr(E)+M
      fme=ipq(e+no)+m
      val=fock(fme)
c     VAL=0.0D0
      if(level.gt.2)then
      DO 2075 F=1,NV
      IF(FZV(F).EQ.1)GO TO 2075
      fsym=orbsym(f+no)
      fn=flov(fsym+1,1)
      ln=flov(fsym+1,2)
      MF=itr(F)+M
      DO 2070 N=fn,ln
      NF=itr(F)+N
      NE=itr(E)+N
      MENF=doff(MAX0(ME,NF))+dadd(MIN0(ME,NF))
      MFNE=doff(MAX0(MF,NE))+dadd(MIN0(MF,NE))
      VAL=VAL+(din(MENF)+din(MENF)-din(MFNE))*T1O(N,F)
 2070 CONTINUE
 2075 CONTINUE
      endif
      FFF(M,E)=VAL
 2080 CONTINUE
 2085 CONTINUE
c     ENDIF
C
C     >F1<
C
      DO 2190 A=1,NV
      IF(FZV(A).EQ.1)GO TO 2190
      asym=orbsym(a+no)
      fn=flov(asym+1,1)
      ln=flov(asym+1,2)
      DO 2185 M=1,NO
      IF(FZO(M).EQ.1)GO TO 2185
      msym=orbsym(m)
      amsym=ieor(asym,msym)
      DO 2180 E=1,NV
      IF(FZV(E).EQ.1)GO TO 2180
      esym=orbsym(e+no)
      AE=itv(MAX0(A,E))+MIN0(A,E)
      fsym=ieor(esym,amsym)+1
      ff=flov(fsym,3)-no
      lf=flov(fsym,4)-no
      DO 2175 F=ff,lf
      MF=itr(F)+M
      AEMF=foff(AE)+fadd(mf)
      F1(A,M,E,F)=fin(AEMF)
C
      IF(LEVEL.GT.2)THEN
      VAL=0.0D0
      DO 2170 N=fn,ln
      NE=itr(E)+N
      NEMF=doff(MAX0(NE,MF))+dadd(MIN0(NE,MF))
      VAL=VAL+din(NEMF)*T1O(N,A)
 2170 CONTINUE
      F1(A,M,E,F)=F1(A,M,E,F)-VAL
      ENDIF
C
 2175 CONTINUE
 2180 CONTINUE
 2185 CONTINUE
 2190 CONTINUE
C
C     >F2<
C
      DO 2290 E=1,NV
      IF(FZV(E).EQ.1)GO TO 2290
      esym=orbsym(e+no)
      DO 2285 I=1,NO
      IF(FZO(I).EQ.1)GO TO 2285
      isym=orbsym(i)
      ff=flov(isym+1,3)-no
      lf=flov(isym+1,4)-no
      eisym=ieor(esym,isym)
      DO 2280 M=1,NO
      IF(FZO(M).EQ.1)GO TO 2280
      msym=orbsym(m)
      nsym=ieor(msym,eisym)+1
      fn=flov(nsym,1)
      ln=flov(nsym,2)
      ME=itr(E)+M
      DO 2275 N=fn,ln
      IN=itr(MAX0(I,N))+MIN0(I,N)
      EMIN=eoff(ME)+eadd(IN)
      F2(E,I,M,N)=ein(EMIN)
C
      IF(LEVEL.GT.2)THEN
      VAL=0.0D0
      DO 2270 F=ff,lf
      NF=itr(F)+N
      MENF=doff(MAX0(ME,NF))+dadd(MIN0(ME,NF))
      VAL=VAL+din(MENF)*T1O(I,F)
 2270 CONTINUE
      F2(E,I,M,N)=F2(E,I,M,N)+VAL
      ENDIF
C
 2275 CONTINUE
 2280 CONTINUE
 2285 CONTINUE
 2290 CONTINUE
C
C     >F3<
C
      DO 2390 A=1,NV
      IF(FZV(A).EQ.1)GO TO 2390
      asym=orbsym(a+no)
      DO 2385 M=1,NO
      IF(FZO(M).EQ.1)GO TO 2385
      msym=orbsym(m)
      amsym=ieor(asym,msym)
      DO 2380 I=1,NO
      IF(FZO(I).EQ.1)GO TO 2380
      isym=orbsym(i)
      ff=flov(isym+1,3)-no
      lf=flov(isym+1,4)-no
      AI=itr(A)+I
      esym=ieor(isym,amsym)+1
      fe=flov(esym,3)-no
      le=flov(esym,4)-no
      DO 2375 E=fe,le
      ME=itr(E)+M
      AIME=doff(MAX0(ME,AI))+dadd(MIN0(ME,AI))
      F3(A,M,I,E)=din(AIME)
C
      IF(LEVEL.GT.2)THEN
      VAL=0.0D0
      DO 2370 F=ff,lf
      AF=itv(MAX0(A,F))+MIN0(A,F)
      MEAF=foff(AF)+fadd(ME)
      VAL=VAL+fin(MEAF)*T1O(I,F)
 2370 CONTINUE
      F3(A,M,I,E)=F3(A,M,I,E)+VAL
      ENDIF
C
 2375 CONTINUE
 2380 CONTINUE
 2385 CONTINUE
 2390 CONTINUE
C
C     >F4<
C
      DO 2490 A=1,NV
      IF(FZV(A).EQ.1)GO TO 2490
      asym=orbsym(a+no)
      DO 2485 M=1,NO
      IF(FZO(M).EQ.1)GO TO 2485
      msym=orbsym(m)
      amsym=ieor(asym,msym)
      DO 2480 E=1,NV
      IF(FZV(E).EQ.1)GO TO 2480
      esym=orbsym(e+no)
      isym=ieor(esym,amsym)+1
      fi=flov(isym,1)
      li=flov(isym,2)
      ff=flov(isym,3)-no
      lf=flov(isym,4)-no
      AE=itv(MAX0(A,E))+MIN0(A,E)
      DO 2475 I=fi,li
      IM=itr(MAX0(I,M))+MIN0(I,M)
      AEIM=coff(AE)+cadd(IM)
      F4(A,M,E,I)=cin(AEIM)
C
      IF(LEVEL.GT.2)THEN
      VAL=0.0D0
      DO 2470 F=ff,lf
      MF=itr(F)+M
      AEMF=foff(AE)+fadd(MF)
      VAL=VAL+fin(AEMF)*T1O(I,F)
 2470 CONTINUE
      F4(A,M,E,I)=F4(A,M,E,I)+VAL
      ENDIF
C
 2475 CONTINUE
 2480 CONTINUE
 2485 CONTINUE
 2490 CONTINUE
C
C     >F5<
C
      DO 2590 I=1,NO
      IF(FZO(I).EQ.1)GO TO 2590
      isym=orbsym(i)
      DO 2585 E=1,NV
      IF(FZV(E).EQ.1)GO TO 2585
      esym=orbsym(e+no)
      iesym=ieor(esym,isym)
      DO 2580 A=1,NV
      IF(FZV(A).EQ.1)GO TO 2580
      asym=orbsym(a+no)
      fn=flov(asym+1,1)
      ln=flov(asym+1,2)
      msym=ieor(asym,iesym)+1
      fm=flov(msym,1)
      lm=flov(msym,2)
      IA=itr(A)+I
      DO 2575 M=fm,lm
      EM=itr(E)+M
      IAEM=doff(MAX0(IA,EM))+dadd(MIN0(IA,EM))
      F5(I,E,A,M)=din(IAEM)
C
      IF(LEVEL.GT.2)THEN
      VAL=0.0D0
      DO 2570 N=fn,ln
      IN=itr(MAX0(I,N))+MIN0(I,N)
      INEM=eoff(EM)+eadd(IN)
      VAL=VAL+ein(INEM)*T1O(N,A)
 2570 CONTINUE
      F5(I,E,A,M)=F5(I,E,A,M)-VAL
      ENDIF
C
 2575 CONTINUE
 2580 CONTINUE
 2585 CONTINUE
 2590 CONTINUE
C
C     >F6<
C
      DO 2690 I=1,NO
      IF(FZO(I).EQ.1)GO TO 2690
      isym=orbsym(i)
      DO 2685 E=1,NV
      IF(FZV(E).EQ.1)GO TO 2685
      esym=orbsym(e+no)
      iesym=ieor(isym,esym)
      DO 2680 M=1,NO
      IF(FZO(M).EQ.1)GO TO 2680
      msym=orbsym(m)
      asym=ieor(msym,iesym)+1
      fa=flov(asym,3)-no
      la=flov(asym,4)-no
      fn=flov(asym,1)
      ln=flov(asym,2)
      IM=itr(MAX0(I,M))+MIN0(I,M)
      DO 2675 A=fa,la
      EA=itv(MAX0(E,A))+MIN0(E,A)
      IMEA=coff(EA)+cadd(IM)
      F6(I,E,M,A)=cin(IMEA)
C
      IF(LEVEL.GT.2)THEN
      VAL=0.0D0
      DO 2670 N=fn,ln
      EN=itr(E)+N
      IMEN=eoff(EN)+eadd(IM)
      VAL=VAL+ein(IMEN)*T1O(N,A)
 2670 CONTINUE
      F6(I,E,M,A)=F6(I,E,M,A)-VAL
      ENDIF
C
 2675 CONTINUE
 2680 CONTINUE
 2685 CONTINUE
 2690 CONTINUE
      ENDIF
C
C     >>CHI1<<
C
      CALL ZERO(CHI1,NV*NV*NO*NV)
      DO 2790 b=1,NV
      IF(FZV(b).EQ.1)GO TO 2790
      bSYM=ORBSYM(b+NO)
      DO 2785 e=1,NV
      IF(FZV(e).EQ.1)GO TO 2785
      eSYM=ORBSYM(e+NO)
      beSYM=IEOR(bSYM,eSYM)
      BE=itv(MAX0(B,E))+MIN0(B,E)
      DO 2780 I=1,NO
      IF(FZO(I).EQ.1)GO TO 2780
      ISYM=ORBSYM(I)
      aSYM=IEOR(ISYM,beSYM)+1
      IASYM=IEOR(ISYM,ASYM-1)
      IBSYM=IEOR(ISYM,BSYM)
      absym=ieor(asym-1,bsym)
      Fa=FLOV(aSYM,3)-NO
      La=FLOV(aSYM,4)-NO
      DO 2775 a=Fa,La
      ZLAB=ZLX(A,B)
      BE=itv(MAX0(B,E))+MIN0(B,E)
      AI=itr(a)+i
      AIBE=foff(be)+fadd(ai)
      CHI1(A,B,I,E)=fin(AIBE)
C
      IF(LEVEL.GT.1)THEN
      VAL=0.0D0
      DO 2760 N=1,NO
      IF(FZO(N).EQ.1)GO TO 2760
      NSYM=ORBSYM(N)
      MSYM=IEOR(ABSYM,NSYM)+1
      FM=FLOV(MSYM,1)
      LM=FLOV(MSYM,2)
      DO 2759 M=FM,LM
      NMAB=UOFF(N,M,ZLAB)+VADD(A,B)
      VAL=VAL+F2(E,I,M,N)*t2o(NMAB)
 2759 CONTINUE
 2760 CONTINUE
      CHI1(A,B,I,E)=CHI1(A,B,I,E)+VAL
C
      VAL=0.0D0
      DO 2750 F=1,NV
      IF(FZV(F).EQ.1)GO TO 2750
      ZLAF=ZLX(A,F)
      ZLBF=ZLX(B,F)
      FSYM=ORBSYM(F+NO)
      MSYM=IEOR(IASYM,FSYM)+1
      FM=FLOV(MSYM,1)
      LM=FLOV(MSYM,2)
      DO 2748 M=FM,LM
      IMAF=UOFF(I,M,ZLAF)+VADD(A,F)
      MIAF=UOFF(M,I,ZLAF)+VADD(A,F)
      VAL=VAL+(F1(B,M,E,F)+F1(B,M,E,F)-F1(B,M,F,E))*t2o(IMAF)
      VAL=VAL- F1(B,M,E,F)*t2o(MIAF)
 2748 CONTINUE
      MSYM=IEOR(IBSYM,FSYM)+1
      FM=FLOV(MSYM,1)
      LM=FLOV(MSYM,2)
      DO 2749 M=FM,LM
      MIBF=UOFF(M,I,ZLBF)+VADD(B,F)
      VAL=VAL- F1(A,M,F,E) *t2o(MIBF)
 2749 CONTINUE
 2750 CONTINUE
      CHI1(A,B,I,E)=CHI1(A,B,I,E)+VAL
C
      IF(LEVEL.GT.2)THEN
      VAL=0.0D0
      DO 2762 N=1,NO
      IF(FZO(N).EQ.1)GO TO 2762
      NSYM=ORBSYM(N)
      MSYM=IEOR(ABSYM,NSYM)+1
      FM=FLOV(MSYM,1)
      LM=FLOV(MSYM,2)
      DO 2769 M=FM,LM
      VAL=VAL+F2(E,I,M,N)*T1O(N,A)*T1O(M,B)
 2769 CONTINUE
 2762 CONTINUE
      CHI1(A,B,I,E)=CHI1(A,B,I,E)+VAL
C
      VAL=0.0D0
      DO 2745 M=1,NO
      IF(FZO(M).EQ.1)GO TO 2745
      VAL=VAL+F3(A,M,I,E)*T1O(M,B)+F4(B,M,E,I)*T1O(M,A)
 2745 CONTINUE
      CHI1(A,B,I,E)=CHI1(A,B,I,E)-VAL
C
      VAL=0.0D0
      BE=itv(MAX0(B,E))+MIN0(B,E)
      ff=flov(isym+1,3)-no
      lf=flov(isym+1,4)-no
      DO 2740 F=ff,lf
      AF=itv(MAX0(A,F))+MIN0(A,F)
      AFBE=boff(MAX0(AF,BE))+badd(MIN0(AF,BE))
      VAL=VAL+bin(AFBE)*T1O(I,F)
 2740 CONTINUE
      CHI1(A,B,I,E)=CHI1(A,B,I,E)+VAL
      ENDIF
C
      IF(LEVEL.GT.4)THEN
      ai=itr(a)+i
      VAL=0.0D0
      DO 2763 M=1,NO
      IF(FZO(M).EQ.1)GO TO 2763
      bm=itr(b)+m
      MSYM=ORBSYM(M)
      ME=itr(E)+M
      DO 2768 N=1,NO
      IF(FZO(N).EQ.1)GO TO 2768
      an=itr(a)+n
      bn=itr(b)+n
      NSYM=ORBSYM(N)
      MNSYM=IEOR(MSYM,NSYM)
      FSYM=IEOR(MNSYM,ESYM)+1
      FF=FLOV(FSYM,3)-NO
      LF=FLOV(FSYM,4)-NO
      DO 2767 F=FF,LF
      NF=itr(F)+N
      MENF=doff(MAX0(ME,NF))+dadd(MIN0(ME,NF))
      fn=itr(f)+n
      fi=itr(f)+i
      fm=itr(f)+m
c     IDX1=X3OFF(I,M,N)+X3ADD(A,B,F)
           k1=max0(fn,ai,bm)
           k3=min0(fn,ai,bm)
           k3s=trsym(k3)
           k2=fn+ai+bm-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(I,M,N)+X3ADD(F,B,A)
           k1=max0(fi,an,bm)
           k3=min0(fi,an,bm)
           k3s=trsym(k3)
           k2=fi+an+bm-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX3=X3OFF(I,M,N)+X3ADD(A,F,B)
           k1=max0(fm,ai,bn)
           k3=min0(fm,ai,bn)
           k3s=trsym(k3)
           k2=fm+ai+bn-k1-k3
           idx3=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      VAL=VAL+din(MENF)*(t3o(IDX1)+t3o(IDX1)-t3o(IDX2)-t3o(IDX3))
 2767 CONTINUE
 2768 CONTINUE
 2763 CONTINUE
      CHI1(A,B,I,E)=CHI1(A,B,I,E)-VAL
      ENDIF
C
      ENDIF
 2775 CONTINUE
 2780 CONTINUE
 2785 CONTINUE
 2790 CONTINUE
C
C     >>CHI2<<
C
      CALL ZERO(CHI2,NO*NO*NO*NV)
      DO 2890 A=1,NV
      IF(FZV(A).EQ.1)GO TO 2890
      ASYM=ORBSYM(A+NO)
      DO 2885 M=1,NO
      IF(FZO(M).EQ.1)GO TO 2885
      MSYM=ORBSYM(M)
      AMSYM=IEOR(ASYM,MSYM)
      DO 2880 I=1,NO
      IF(FZO(I).EQ.1)GO TO 2880
      ISYM=ORBSYM(I)
      IASYM=IEOR(ISYM,ASYM)
      IMSYM=IEOR(ISYM,MSYM)
      JSYM=IEOR(AMSYM,ISYM)+1
      IJSYM=IEOR(ISYM,JSYM-1)
      AJSYM=IEOR(ASYM,JSYM-1)
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      DO 2875 J=FJ,LJ
      AI=itr(A)+I
      JM=itr(MAX0(J,M))+MIN0(J,M)
      AIJM=eoff(AI)+eadd(JM)
      CHI2(A,M,I,J)=ein(AIJM)
C
      IF(LEVEL.GT.1)THEN
      VAL=0.0D0
      DO 2860 E=1,NV
      IF(FZV(E).EQ.1)GO TO 2860
      ESYM=ORBSYM(E+NO)
      FSYM=IEOR(IJSYM,ESYM)+1
      FF=FLOV(FSYM,3)-NO
      LF=FLOV(FSYM,4)-NO
      DO 2859 F=FF,LF
      ZL=ZLX(E,F)
      IJEF=UOFF(I,J,ZL)+VADD(E,F)
      VAL=VAL+F1(A,M,E,F)*t2o (IJEF)
 2859 CONTINUE
 2860 CONTINUE
      CHI2(A,M,I,J)=CHI2(A,M,I,J)+VAL
C
      VAL=0.0D0
      DO 2850 E=1,NV
      IF(FZV(E).EQ.1)GO TO 2850
      ZLAE=ZLX(A,E)
      ESYM=ORBSYM(E+NO)
      NSYM=IEOR(IASYM,ESYM)+1
      FN=FLOV(NSYM,1)
      LN=FLOV(NSYM,2)
      DO 2848 N=FN,LN
      INAE=UOFF(I,N,ZLAE)+VADD(A,E)
      NIAE=UOFF(N,I,ZLAE)+VADD(A,E)
      VAL=VAL+(F2(E,J,N,M)+F2(E,J,N,M)-F2(E,J,M,N))*t2o(INAE)
      VAL=VAL- F2(E,J,N,M)*t2o(NIAE)
 2848 CONTINUE
      NSYM=IEOR(AJSYM,ESYM)+1
      FN=FLOV(NSYM,1)
      LN=FLOV(NSYM,2)
      DO 2849 N=FN,LN
      NJAE=UOFF(N,J,ZLAE)+VADD(A,E)
      VAL=VAL- F2(E,I,M,N) *t2o(NJAE)
 2849 CONTINUE
 2850 CONTINUE
      CHI2(A,M,I,J)=CHI2(A,M,I,J)+VAL
C
      VAL=0.0D0
      EESYM=IEOR(JSYM-1,IASYM)+1
      FE=FLOV(EESYM,3)-NO
      LE=FLOV(EESYM,4)-NO
      DO 2865 E=FE,LE
      ZL=ZLX(A,E)
      IJAE=UOFF(I,J,ZL)+VADD(A,E)
      VAL=VAL+FFF(M,E)*t2o(IJAE)
 2865 CONTINUE
      CHI2(A,M,I,J)=CHI2(A,M,I,J)+VAL
C
      IF(LEVEL.GT.2)THEN
      VAL=0.0D0
      DO 2861 E=1,NV
      IF(FZV(E).EQ.1)GO TO 2861
      ESYM=ORBSYM(E+NO)
      FSYM=IEOR(IJSYM,ESYM)+1
      FF=FLOV(FSYM,3)-NO
      LF=FLOV(FSYM,4)-NO
      DO 2862 F=FF,LF
      VAL=VAL+F1(A,M,E,F)*T1O(I,E)*T1O(J,F)
 2862 CONTINUE
 2861 CONTINUE
      CHI2(A,M,I,J)=CHI2(A,M,I,J)+VAL
C
      VAL=0.0D0
      DO 2845 E=1,NV
      IF(FZV(E).EQ.1)GO TO 2845
      VAL=VAL+F5(I,E,A,M)*T1O(J,E)+F6(J,E,M,A)*T1O(I,E)
 2845 CONTINUE
      CHI2(A,M,I,J)=CHI2(A,M,I,J)+VAL
C
      VAL=0.0D0
      fn=flov(asym+1,1)
      ln=flov(asym+1,2)
      DO 2840 N=fn,ln
      IN=itr(MAX0(I,N))+MIN0(I,N)
      INJM=aoff(MAX0(IN,JM))+aadd(MIN0(IN,JM))
      VAL=VAL+ain(INJM)*T1O(N,A)
 2840 CONTINUE
      CHI2(A,M,I,J)=CHI2(A,M,I,J)-VAL
      ENDIF
C
      IF(LEVEL.GT.4)THEN
      ai=itr(a)+i
      VAL=0.0D0
      DO 2870 N=1,NO
      IF(FZO(N).EQ.1)GO TO 2870
      NSYM=ORBSYM(N)
      an=itr(a)+n
      DO 2869 E=1,NV
      IF(FZV(E).EQ.1)GO TO 2869
      ej=itr(e)+j
      en=itr(e)+n
      ESYM=ORBSYM(E+NO)
      ME=itr(E)+M
      NESYM=IEOR(NSYM,ESYM)
      FSYM=IEOR(NESYM,MSYM)+1
      FF=FLOV(FSYM,3)-NO
      LF=FLOV(FSYM,4)-NO
      DO 2868 F=FF,LF
      NF=itr(F)+N
      MENF=doff(MAX0(ME,NF))+dadd(MIN0(ME,NF))
      fi=itr(f)+i
      fj=itr(f)+j
      fn=itr(f)+n
c     IDX1=X3OFF(I,J,N)+X3ADD(A,E,F)
           k1=max0(fn,ai,ej)
           k3=min0(fn,ai,ej)
           k3s=trsym(k3)
           k2=fn+ai+ej-k1-k3
           idx1=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX2=X3OFF(I,J,N)+X3ADD(F,E,A)
           k1=max0(fi,an,ej)
           k3=min0(fi,an,ej)
           k3s=trsym(k3)
           k2=fi+an+ej-k1-k3
           idx2=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c     IDX3=X3OFF(I,J,N)+X3ADD(A,F,E)
           k1=max0(fj,ai,en)
           k3=min0(fj,ai,en)
           k3s=trsym(k3)
           k2=fj+ai+en-k1-k3
           idx3=i1off(k1)+i2off(k2,k3s)+i3off(k3)
      VAL=VAL+din(MENF)*(t3o(IDX1)+t3o(IDX1)-t3o(IDX2)-t3o(IDX3))
 2868 CONTINUE
 2869 CONTINUE
 2870 CONTINUE
      CHI2(A,M,I,J)=CHI2(A,M,I,J)+VAL
      ENDIF
C
      ENDIF
 2875 CONTINUE
 2880 CONTINUE
 2885 CONTINUE
 2890 CONTINUE
C
C     CALCULATE FIRST SET OF CONTRIBUTIONS
C
      DO 2990 I=1,NO
      IF(FZO(I).EQ.1)GO TO 2990
      ISYM=ORBSYM(I)
      DO 2980 J=1,NO
      IF(FZO(J).EQ.1)GO TO 2980
      JSYM=ORBSYM(J)
      DO 2970 K=1,NO
      IF(FZO(K).EQ.1)GO TO 2970
      KSYM=ORBSYM(K)
      JKSYM=IEOR(JSYM,KSYM)
      IJKSYM=IEOR(ISYM,JKSYM)
      DO 2960 A=1,NV
      IF(FZV(A).EQ.1)GO TO 2960
      ai=itr(a)+i
      ASYM=ORBSYM(A+NO)
      DO 2950 B=1,NV
      IF(FZV(B).EQ.1)GO TO 2950
      bj=itr(b)+j
c     if(bj.gt.ai)go to 2960
      BSYM=ORBSYM(B+NO)
      ABSYM=IEOR(ASYM,BSYM)
      CSYM=IEOR(IJKSYM,ABSYM)
      FC=FLOV(CSYM+1,3)-NO
      LC=FLOV(CSYM+1,4)-NO
      BCSYM=IEOR(BSYM,CSYM)
      ESYM=IEOR(JKSYM,CSYM)+1
      FE=FLOV(ESYM,3)-NO
      LE=FLOV(ESYM,4)-NO
      MSYM=IEOR(BCSYM,KSYM)+1
      FM=FLOV(MSYM,1)
      LM=FLOV(MSYM,2)
      DO 2940 C=FC,LC
      ck=itr(c)+k
c     if(ck.gt.bj)go to 2940
      VAL=0.0D0
      ZLBC=ZLX(B,C)
C
      DO 2930 E=FE,LE
      ZLCE=ZLX(C,E)
      KJCE=UOFF(K,J,ZLCE)+VADD(C,E)
      VAL=VAL+CHI1(A,B,I,E)*t2o(KJCE)
 2930 CONTINUE
C
      DO 2920 M=FM,LM
      MKBC=UOFF(M,K,ZLBC)+VADD(B,C)
      VAL=VAL-CHI2(A,M,I,J)*t2o(MKBC)
 2920 CONTINUE
C
c     IDX1=X3OFF(I,J,K)+X3ADD(A,B,C)
c     IDX2=X3OFF(J,K,I)+X3ADD(B,C,A)
c     IDX3=X3OFF(K,I,J)+X3ADD(C,A,B)
c     IDX4=X3OFF(I,K,J)+X3ADD(A,C,B)
c     IDX5=X3OFF(K,J,I)+X3ADD(C,B,A)
c     IDX6=X3OFF(J,I,K)+X3ADD(B,A,C)
c     t3(IDX1)=t3(IDX1)+VAL
c     t3(IDX2)=t3(IDX2)+VAL
c     t3(IDX3)=t3(IDX3)+VAL
c     t3(IDX4)=t3(IDX4)+VAL
c     t3(IDX5)=t3(IDX5)+VAL
c     t3(IDX6)=t3(IDX6)+VAL
c       xfac=1.0d0
c       if(ai.eq.bj.or.ai.eq.ck.or.bj.eq.ck)xfac=2.0d0
c       if(ai.eq.bj.and.ai.eq.ck)xfac=6.0d0
           k1=max0(ai,bj,ck)
           k3=min0(ai,bj,ck)
           k3s=trsym(k3)
           k2=ai+bj+ck-k1-k3
           idxn=i1off(k1)+i2off(k2,k3s)+i3off(k3)
c          t3(IDXn)=t3(idxn)+VAL*xfac
           t3(IDXn)=t3(idxn)+VAL
 2940 CONTINUE
 2950 CONTINUE
 2960 CONTINUE
 2970 CONTINUE
 2980 CONTINUE
 2990 CONTINUE
C
c   renormalize triple and double diagonal elements
c
      DO 4690 A=1,NV
         IF(FZV(A).EQ.1)GO TO 4690
         ASYM=ORBSYM(A+NO)
         DO 4685 I=1,NO
            IF(FZO(I).EQ.1)GO TO 4685
            ISYM=ORBSYM(I)
            AI=ITR(A)+I
            if(asym.ne.isym)go to 4685
            aiaiai=i1off(ai)+i2off(ai,1)+i3off(ai)
            t3(aiaiai)=t3(aiaiai)*6.0d0
 4685    CONTINUE
 4690 CONTINUE
c
      DO 4790 A=1,NV
         IF(FZV(A).EQ.1)GO TO 4790
         ASYM=ORBSYM(A+NO)
         DO 4785 I=1,NO
            IF(FZO(I).EQ.1)GO TO 4785
            ISYM=ORBSYM(I)
            AI=ITR(A)+I
            ais=ieor(asym,isym)+1
            do 4780 b=1,nv
               IF(FZV(b).EQ.1)GO TO 4780
               bSYM=ORBSYM(b+NO)
               DO 4775 j=1,NO
                  IF(FZO(j).EQ.1)GO TO 4775
                  jSYM=ORBSYM(j)
                  if(bsym.ne.jsym)go to 4775
                  bj=ITR(b)+j
                  if(ai.eq.bj)go to 4775
                  if(ai.gt.bj)then
                     aiaibj=i1off(ai)+i2off(ai,1)+i3off(bj)
                     t3(aiaibj)=t3(aiaibj)*2.0d0
                  else
                     bjaiai=i1off(bj)+i2off(ai,ais)+i3off(ai)
                     t3(bjaiai)=t3(bjaiai)*2.0d0
                  endif
 4775          CONTINUE
 4780       CONTINUE
 4785    CONTINUE
 4790 CONTINUE
c
      return
      end
