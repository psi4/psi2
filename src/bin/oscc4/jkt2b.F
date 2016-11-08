      SUBROUTINE JKT2b (T1O,T2O,T2,CIN,DIN,EIN,FIN,AR1,AR2,
     .                  SUX1,SUX2,JAR,KAR,AUX1,
     .                  ORBSYM,FLOV,COFF,CADD,DOFF,DADD,
     .                  EOFF,EADD,FOFF,FADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,MAXOV,ITAP91,
     .                  map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,VALJ,VALK,VAL1,VAL2,xval,sgn,tval
      REAL*8  T1O(NO,NV),T2O(*),T2(*),CIN(*),DIN(*),EIN(*),FIN(*),
     .        SUX1(*),SUX2(*),JAR(*),KAR(*),AUX1(NV,NO),
     .        AR1(MAXOV,MAXOV),AR2(MAXOV,MAXOV)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        COFF(*),CADD(*),DOFF(*),DADD(*),
     .        EOFF(*),EADD(*),FOFF(*),FADD(*),
     .        map(*),spin(*),ssov(2,4),xspin(6),yspin(6),zspin(6)
C
      CALL TRNPS(AUX1,T1O,NO,NV)
C
      xspin(1)=1
      yspin(1)=2

      xspin(2)=2
      yspin(2)=1
c
c     write(6,*)' in jkt2b,  maxov=',maxov
c     write(6,*)' in jkt2b,  cints='
c     write(6,897)(cin(ikj),ikj=1,21)
  897 format(f20.12)
      ITAP91=91
      CALL RFILE(ITAP91)
c
      do 998 case=1,2
      spin1=xspin(case)
      spin2=yspin(case)
c     write(6,*)'B CASE =',case,spin1,spin2
c
      PON=1
      DO 695 TSYM=1,NIRRED
         XAI=0
         DO 190 I=1,NO
            IF(FZO(I).EQ.1)GO TO 190
            ispn=spin(i)
            if(ispn.ne.spin1)go to 190
            xi=map(i)
            ISYM=ORBSYM(I)
            fa=ssov(spin2,3)-no
            la=ssov(spin2,4)-no
            DO 170 A=fa,la
               if(fzv(a).eq.1)go to 170
               xa=map(a+no)
               ASYM=ORBSYM(A+NO)
               IASYM=IEOR(ISYM,ASYM)
               IF(IASYM.NE.TSYM-1)GO TO 170
               XAI=XAI+1
               XBJ=0
               fj=ssov(spin2,1)
               lj=ssov(spin2,2)
               DO 150 J=fj,lj
                  IF(FZO(J).EQ.1)GO TO 150
                  xj=map(j)
                  ja=ipq(max0(xj,xa))+min0(xj,xa)
                  JSYM=ORBSYM(J)
                  BSYM=IEOR(JSYM,IASYM)+1
                  FB=FLOV(BSYM,spin1,3)-NO
                  LB=FLOV(BSYM,spin1,4)-NO
                  DO 130 B=FB,LB
                     XBJ=XBJ+1
                     xb=map(b+no)
                     ib=ipq(max0(xi,xb))+min0(xi,xb)
                     IbJa=DOFF(MAX0(Ib,Ja))+DADD(MIN0(Ib,Ja))
                     xval=din(ibja)
                     AR1(XBJ,XAI)=-xval
  130             CONTINUE
  150          CONTINUE
c          write(6,*)'xbj=',xbj
  170       CONTINUE
  190    CONTINUE
c          write(6,*)'xai=',xai
c
         DO 690 U=1,NO
            IF(FZO(U).EQ.1)GO TO 690
            uspn=spin(u)
            if(uspn.ne.spin1)go to 690
            xu=map(u)
            USYM=ORBSYM(U)
            DO 680 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 680
               bespn=spin(be+no)
               if(bespn.ne.spin2)go to 680
               xbe=map(be+no)
               BESYM=ORBSYM(BE+NO)
               BEUSYM=IEOR(BESYM,USYM)
               IF(BEUSYM.NE.TSYM-1)GO TO 680
C
               BEU=Ipq(max0(xbe,xu))+min0(xbe,xu)
               XAI=0
               fi=ssov(spin1,1)
               li=ssov(spin1,2)
               DO 490 I=fi,li
                  IF(FZO(i).EQ.1)GO TO 490
                  xi=map(i)
                  ISYM=ORBSYM(I)
                  ASYM=IEOR(ISYM,BEUSYM)+1
                  FA=FLOV(ASYM,spin2,3)-NO
                  LA=FLOV(ASYM,spin2,4)-NO
                  IU=ipq(MAX0(xI,xU))+MIN0(xI,xU)
                  BEi=Ipq(max0(xbe,xi))+min0(xbe,xi)
                  DO 470 A=FA,LA
                     xa=map(a+no)
                     uA=ipq(MAX0(xu,xA))+MIN0(xu,xA)
                     BEA=ipq(MAX0(xBE,xA))+MIN0(xBE,xA)
                     BEAIU=COFF(max0(BEA,iu))+CADD(min0(bea,IU))
                     valj=CIN(BEAIU)
C
                     FB=FLOV(USYM+1,uspn,3)-NO
                     LB=FLOV(USYM+1,uspn,4)-NO
                     DO 250 B=FB,LB
                        xb=map(b+no)
                        ib=ipq(max0(xi,xb))+min0(xi,xb)
                        BEaIb=FOFF(max0(BEa,ib))+FADD(min0(bea,Ib))
                        xval=fin(beaib)
                        VALj=VALj+ xval*AUX1(B,U)
  250                CONTINUE
c
                     FJ=FLOV(BESYM+1,bespn,1)
                     LJ=FLOV(BESYM+1,bespn,2)
                     DO 255 J=FJ,LJ
                        xj=map(j)
                        ja=ipq(max0(xj,xa))+min0(xj,xa)
                        iUjA=EOFF(max0(iu,jA))+EADD(min0(iU,ja))
                        xval=ein(iuja)
                        VALj=VALj- xval    *T1O(J,BE)
  255                CONTINUE
C
                     XAI=XAI+1
                     lim=xai
                     JAR(XAI)= +VALJ
C
  470             CONTINUE
  490          CONTINUE
C
               fi=ssov(spin2,1)
               li=ssov(spin2,2)
               DO 491 I=fi,li
                  IF(FZO(i).EQ.1)GO TO 491
                  xi=map(i)
                  ISYM=ORBSYM(I)
                  ASYM=IEOR(ISYM,BEUSYM)+1
                  FA=FLOV(ASYM,spin1,3)-NO
                  LA=FLOV(ASYM,spin1,4)-NO
                  BEi=Ipq(max0(xbe,xi))+min0(xbe,xi)
                  DO 471 A=FA,LA
                     xa=map(a+no)
                     uA=ipq(MAX0(xu,xA))+MIN0(xu,xA)
c
                     VAL=0.0D0
                     FB=FLOV(USYM+1,uspn,3)-NO
                     LB=FLOV(USYM+1,uspn,4)-NO
                     DO 251 B=FB,LB
                        xb=map(b+no)
                        aB=ipq(max0(xa,xb))+min0(xa,xb)
                        BEiab=FOFF(max0(BEi,ab))+FADD(min0(bei,ab))
                        VAL =VAL + fin(beiab)*AUX1(B,U)
  251                CONTINUE
c
                     FJ=FLOV(BESYM+1,bespn,1)
                     LJ=FLOV(BESYM+1,bespn,2)
                     DO 256 J=FJ,LJ
                        xj=map(j)
                        Ji=ipq(max0(xj,xi))+min0(xj,xi)
                        JiuA=EOFF(max0(ji,uA))+EADD(min0(Ji,ua))
                        VAL =VAL -ein(jiua)*T1O(J,BE)
  256                CONTINUE
C
               if(u.ne.i.and.be.ne.a)then
                     UIBEA=TOFF(U,I)+TADD(BE,A)
                     T2(UIBEA)=T2(UIBEA)+ VAL*0.5d0
               endif
  471             CONTINUE
  491          CONTINUE
c
               XBJ=0
               fj=ssov(spin2,1)
               lj=ssov(spin2,2)
               DO 580 J=fj,lj
                  IF(FZO(j).EQ.1)GO TO 580
                  JSYM=ORBSYM(J)
                  BSYM=IEOR(JSYM,BEUSYM)+1
                  FB=FLOV(BSYM,spin1,3)-NO
                  LB=FLOV(BSYM,spin1,4)-NO
                  DO 570 B=FB,LB
                     sgn=1.0d0
                     if(u.eq.j.or.be.eq.b)sgn=0.0d0
                     if(be.gt.b)sgn=-sgn
                     if(j.gt.u)sgn=-sgn
                     UJBBE=TOFF(U,J)+TADD(B,BE)
                     XBJ=XBJ+1
                     SUX1(XBJ)=0.5D0*T2O(UJBBE)*sgn+T1O(U,B)*T1O(J,BE)
  570             CONTINUE
  580          CONTINUE
cbug !!!       lim=xbj
               lim2=xbj
c              write(6,*)'lim.lim2',lim,lim2
c
               DO 642 AI=1,LIM
                  VALJ=0.0D0
                  DO 640 BJ=1,LIM2
                     VALJ=VALJ+AR1(BJ,AI)*SUX1(BJ)
  640             CONTINUE
                  JAR(AI)=JAR(AI)+VALJ
  642          CONTINUE
c
               CALL WWRITW(ITAP91,JAR,INTOWP(LIM),PON,PON)
c              write(6,*)'b wr tsym lim pon',tsym,lim,pon
C
  680       CONTINUE
  690    CONTINUE
  695 CONTINUE
C
      PON=1
      DO 995 TSYM=1,NIRRED
         XVG=0
         DO 790 V=1,no
            IF(FZO(V).EQ.1)GO TO 790
            VSYM=ORBSYM(V)
            vspn=spin(v)
            if(vspn.ne.spin2)go to 790
            DO 770 GA=1,nv
               IF(FZV(GA).EQ.1)GO TO 770
               gaspn=spin(ga+no)
               if(gaspn.ne.spin1)go to 770
               GASYM=ORBSYM(GA+NO)
               VGASYM=IEOR(VSYM,GASYM)
               IF(VGASYM.NE.TSYM-1)GO TO 770
               XVG=XVG+1
               XAI=0
               fi=ssov(spin1,1)
               li=ssov(spin1,2)
               DO 750 I=fi,li
                  IF(FZO(i).EQ.1)GO TO 750
                  ISYM=ORBSYM(I)
                  ASYM=IEOR(ISYM,VGASYM)+1
                  FA=FLOV(ASYM,spin2,3)-NO
                  LA=FLOV(ASYM,spin2,4)-NO
                  DO 730 A=FA,LA
                     sgn=1.0d0
                     if(i.eq.v.or.a.eq.ga)sgn=0.0d0
                     if(a.lt.ga)sgn=-sgn 
                     if(v.gt.i)sgn=-sgn 
                     XAI=XAI+1
                     IVAGA=TOFF(I,V)+TADD(A,GA)
                     AR1(XAI,XVG)=T2O(IVAGA)*sgn
  730             CONTINUE
  750          CONTINUE
  770       CONTINUE
  790    CONTINUE
cbug!!!!!!!        LIM=XVG
         LIM=Xai
c        write(6,*)'xai,xvg',xai,xvg
C
         DO 990 U=1,NO
            IF(FZO(U).EQ.1)GO TO 990
            USYM=ORBSYM(U)
            uspn=spin(u)
            if(uspn.ne.spin1)go to 990
            DO 980 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 980
               bespn=spin(be+no)
               if(bespn.ne.spin2)go to 980
               BESYM=ORBSYM(BE+NO)
               BEUSYM=IEOR(BESYM,USYM)
               IF(BEUSYM.NE.TSYM-1)GO TO 980
               CALL WREADW(ITAP91,JAR,INTOWP(LIM),PON,PON)
c              write(6,*)'b RD tsym lim pon',tsym,lim,pon
C
               VGA=0
               DO 860 V=1,no
                  IF(FZO(V).EQ.1)GO TO 860
                  VSYM=ORBSYM(V)
                  GASYM=IEOR(VSYM,BEUSYM)+1
                  vspn=spin(v)
                  if(vspn.ne.spin2)go to 860
                  FGA=FLOV(GASYM,spin1,3)-NO
                  LGA=FLOV(GASYM,spin1,4)-NO
                  DO 850 GA=FGA,LGA
                     VGA=VGA+1
               if(u.ne.v.and.be.ne.ga)then
                     VAL1=0.0D0
                     DO 810 AI=1,LIM
                        VAL1=VAL1+JAR(AI)*AR1(AI,VGA)
c             write(6,678)map(u),map(v),map(be+no),map(ga+no),spin(u),
c    .       spin(v), spin(be+no),spin(ga+no),jar(ai),ar1(ai,vga),val1
  678      format(4i3,3x,4i2,3x,3f12.8)
  810                CONTINUE
                     UVBEGA=TOFF(U,V)+TADD(BE,GA)
                     T2(UVBEGA)=T2(UVBEGA)+val1
               endif
  850             CONTINUE
  860          CONTINUE
C
  980       CONTINUE
  990    CONTINUE
  995 CONTINUE
c
  998 continue
C
      CALL RCLOSE(ITAP91,3)
c
      RETURN
      END
