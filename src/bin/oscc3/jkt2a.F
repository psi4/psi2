      SUBROUTINE JKT2a (T1O,T2O,T2,CIN,DIN,EIN,FIN,AR1,AR2,
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
     .        map(*),spin(*),ssov(2,4),xspin(8),yspin(8),zspin(8)
C
      CALL TRNPS(AUX1,T1O,NO,NV)
C
      xspin(1)=1
      yspin(1)=1
      zspin(1)=1

      xspin(2)=1
      yspin(2)=2
      zspin(2)=1

      xspin(3)=1
      yspin(3)=1
      zspin(3)=2

      xspin(4)=1
      yspin(4)=2
      zspin(4)=2

      xspin(5)=2
      yspin(5)=1
      zspin(5)=1

      xspin(6)=2
      yspin(6)=2
      zspin(6)=1

      xspin(7)=2
      yspin(7)=1
      zspin(7)=2

      xspin(8)=2
      yspin(8)=2
      zspin(8)=2
c
      ITAP91=91
      CALL RFILE(ITAP91)
c
      do 998 case=1,8
      spin1=xspin(case)
      spin2=yspin(case)
      spin3=zspin(case)
c     write(6,*)' this is A CASE =',case,'    SPINS=',spin1,spin2,spin3
      if(case.eq.2.or.case.eq.4.or.case.eq.6.or.case.eq.8)go to 700
c
      PON=1
      DO 1695 TSYM=1,NIRRED
         XAI=0
         DO 1190 I=1,NO
            IF(FZO(I).EQ.1)GO TO 1190
            ispn=spin(i)
            if(ispn.ne.spin3)go to 1190
            xi=map(i)
            ISYM=ORBSYM(I)
            fa=ssov(spin3,3)-no
            la=ssov(spin3,4)-no
            DO 1170 A=fa,la
               if(fzv(a).eq.1)go to 1170
               xa=map(a+no)
               ASYM=ORBSYM(A+NO)
               IASYM=IEOR(ISYM,ASYM)
               IF(IASYM.NE.TSYM-1)GO TO 1170
               ia=ipq(max0(xi,xa))+min0(xi,xa)
               XAI=XAI+1
               XBJ=0
               fj=ssov(1,1)
               lj=ssov(1,2)
               DO 1150 J=fj,lj
                  IF(FZO(J).EQ.1)GO TO 1150
                  xj=map(j)
                  ja=ipq(max0(xj,xa))+min0(xj,xa)
                  JSYM=ORBSYM(J)
                  BSYM=IEOR(JSYM,IASYM)+1
                  FB=FLOV(BSYM,1,3)-NO
                  LB=FLOV(BSYM,1,4)-NO
                  DO 1130 B=FB,LB
                     XBJ=XBJ+1
                     xb=map(b+no)
                     jb=ipq(max0(xj,xb))+min0(xj,xb)
                     IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                     xval=din(iajb)
                     if(spin3.eq.1)then
                        ib=ipq(max0(xi,xb))+min0(xi,xb)
                        IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                        xval=xval-din(ibja)
                     endif
                     AR1(XBJ,XAI)=-xval
 1130             CONTINUE
 1150          CONTINUE
 1170       CONTINUE
 1190    CONTINUE
         limai=xai
c        write(6,*)'limai=',limai
C
         DO 1690 U=1,NO
            IF(FZO(U).EQ.1)GO TO 1690
            uspn=spin(u)
            if(uspn.ne.spin1)go to 1690
            xu=map(u)
            USYM=ORBSYM(U)
            DO 1680 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 1680
               bespn=spin(be+no)
               if(bespn.ne.spin1)go to 1680
               xbe=map(be+no)
               BESYM=ORBSYM(BE+NO)
               BEUSYM=IEOR(BESYM,USYM)
               IF(BEUSYM.NE.TSYM-1)GO TO 1680
c
               XBJ=0
               fj=ssov(1,1)
               lj=ssov(1,2)
               DO 1580 J=fj,lj
                  IF(FZO(J).EQ.1)GO TO 1580
                  JSYM=ORBSYM(J)
                  BSYM=IEOR(JSYM,BEUSYM)+1
                  FB=FLOV(BSYM,1,3)-NO
                  LB=FLOV(BSYM,1,4)-NO
                  DO 1570 B=FB,LB
                     sgn=1.0d0
                     if(j.eq.u.or.be.eq.b)sgn=0.0d0
                     if(j.gt.u)sgn=-sgn
                     if(be.gt.b)sgn=-sgn
                     UJBBE=TOFF(U,J)+TADD(B,BE)
                     XBJ=XBJ+1
                     SUX1(XBJ)=0.5D0*T2O(UJBBE)*sgn+T1O(U,B)*T1O(J,BE)
 1570             CONTINUE
 1580          CONTINUE
               LIMbj=XBJ
c        write(6,*)'limbj=',limbj
c
               DO 1642 AI=1,LIMai
                  VALJ=0.0D0
                  DO 1640 BJ=1,LIMbj
                     VALJ=VALJ+AR1(BJ,AI)*SUX1(BJ)
 1640             CONTINUE
                  JAR(AI)=+VALJ
 1642          CONTINUE
c
c             write(6,*)'WR1,tsym,limai,pon',tsym,limai,pon
               CALL WWRITW(ITAP91,JAR,INTOWP(LIMai),PON,PON)
c
 1680       CONTINUE
 1690    CONTINUE
 1695 CONTINUE
c
      PON=1
      DO 695 TSYM=1,NIRRED
         XAI=0
         DO 190 I=1,NO
            IF(FZO(I).EQ.1)GO TO 190
            ispn=spin(i)
            if(ispn.ne.spin3)go to 190
            xi=map(i)
            ISYM=ORBSYM(I)
            fa=ssov(spin3,3)-no
            la=ssov(spin3,4)-no
            DO 170 A=fa,la
               if(fzv(a).eq.1)go to 170
               xa=map(a+no)
               ASYM=ORBSYM(A+NO)
               IASYM=IEOR(ISYM,ASYM)
               IF(IASYM.NE.TSYM-1)GO TO 170
               ia=ipq(max0(xi,xa))+min0(xi,xa)
               XAI=XAI+1
               XBJ=0
               fj=ssov(2,1)
               lj=ssov(2,2)
               DO 150 J=fj,lj
                  IF(FZO(J).EQ.1)GO TO 150
                  xj=map(j)
                  ja=ipq(max0(xj,xa))+min0(xj,xa)
                  JSYM=ORBSYM(J)
                  BSYM=IEOR(JSYM,IASYM)+1
                  FB=FLOV(BSYM,2,3)-NO
                  LB=FLOV(BSYM,2,4)-NO
                  DO 130 B=FB,LB
                     XBJ=XBJ+1
                     xb=map(b+no)
                     jb=ipq(max0(xj,xb))+min0(xj,xb)
                     IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                     xval=din(iajb)
                     if(spin3.eq.2)then
                        ib=ipq(max0(xi,xb))+min0(xi,xb)
                        IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                        xval=xval-din(ibja)
                     endif
                     AR1(XBJ,XAI)=-xval
  130             CONTINUE
  150          CONTINUE
  170       CONTINUE
  190    CONTINUE
         LIMai=XAI
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
               if(bespn.ne.spin1)go to 680
               xbe=map(be+no)
               BESYM=ORBSYM(BE+NO)
               BEUSYM=IEOR(BESYM,USYM)
               IF(BEUSYM.NE.TSYM-1)GO TO 680
C
c             write(6,*)'RD1,tsym,limai,pon',tsym,limai,pon
               CALL Wreadw(ITAP91,JAR,INTOWP(limai),PON,PON1)
C
               BEU=Ipq(max0(xbe,xu))+min0(xbe,xu)
               XAI=0
               fi=ssov(spin3,1)
               li=ssov(spin3,2)
               DO 490 I=fi,li
                  IF(FZO(i).EQ.1)GO TO 490
                  xi=map(i)
                  ISYM=ORBSYM(I)
                  ASYM=IEOR(ISYM,BEUSYM)+1
                  FA=FLOV(ASYM,spin3,3)-NO
                  LA=FLOV(ASYM,spin3,4)-NO
                  IU=ipq(MAX0(xI,xU))+MIN0(xI,xU)
                  BEi=Ipq(max0(xbe,xi))+min0(xbe,xi)
                  DO 470 A=FA,LA
                     xa=map(a+no)
                     ia=ipq(max0(xi,xa))+min0(xi,xa)
                     ua=ipq(max0(xu,xa))+min0(xu,xa)
                     BEUIA=DOFF(MAX0(BEU,IA))+DADD(MIN0(BEU,IA))
                     VALJ=DIN(BEUIA)
                     if(spin1.eq.spin3)then
                        BEA=ipq(MAX0(xBE,xA))+MIN0(xBE,xA)
                        BEAIU=COFF(max0(BEA,iu))+CADD(min0(bea,IU))
                        valj=valj-CIN(BEAIU)
                     endif
C
                     VAL=0.0D0
                     FB=FLOV(USYM+1,uspn,3)-NO
                     LB=FLOV(USYM+1,uspn,4)-NO
                     DO 250 B=FB,LB
                        xb=map(b+no)
                        BEB=ipq(max0(xbe,xb))+min0(xbe,xb)
                        BEBIA=FOFF(max0(BEB,ia))+FADD(min0(beb,IA))
                        tval=fin(bebia)
                        xval=fin(bebia)
                        if(spin1.eq.spin3)then
                           ib=ipq(max0(xi,xb))+min0(xi,xb)
                           BEaIb=FOFF(max0(BEa,ib))+FADD(min0(bea,Ib))
                           xval=xval-fin(beaib)
                           aB=ipq(max0(xa,xb))+min0(xa,xb)
                           BEiab=FOFF(max0(BEi,ab))+FADD(min0(bei,ab))
                           tval=tval-fin(beiab)
                        endif
                        VAL =VAL + tval     *AUX1(B,U)
                        VALj=VALj+ xval     *AUX1(B,U)
  250                CONTINUE
c
                     FJ=FLOV(BESYM+1,bespn,1)
                     LJ=FLOV(BESYM+1,bespn,2)
                     DO 255 J=FJ,LJ
                        xj=map(j)
                        JU=ipq(max0(xj,xu))+min0(xj,xu)
                        JUIA=EOFF(max0(ju,IA))+EADD(min0(JU,ia))
                        tval=ein(juia)
                        xval=ein(juia)
                        if(spin1.eq.spin3)then
                           Ji=ipq(max0(xj,xi))+min0(xj,xi)
                           JiuA=EOFF(max0(ji,uA))+EADD(min0(Ji,ua))
                           tval=tval-ein(jiua)
                           ja=ipq(max0(xj,xa))+min0(xj,xa)
                           iUjA=EOFF(max0(iu,jA))+EADD(min0(iU,ja))
                           xval=xval-ein(iuja)
                        endif
                        VAL =VAL - tval    *T1O(J,BE)
                        VALj=VALj- xval    *T1O(J,BE)
  255                CONTINUE
C
                     XAI=XAI+1
                     limai=xai
                     JAR(XAI)=jar(xai)+VALJ
C
c     if(case.eq.2.or.case.eq.4.or.case.eq.6.or.case.eq.8)then      
             if(u.eq.i.or.be.eq.a)val=0.0d0
         if((u.lt.i.and.be.gt.a).or.(u.gt.i.and.be.lt.a))val=-val
                        UIBEA=TOFF(U,I)+TADD(BE,A)
                        T2(UIBEA)=T2(UIBEA)+ VAL*0.5d0
c             endif
  470             CONTINUE
  490          CONTINUE
c
               XBJ=0
               fj=ssov(2,1)
               lj=ssov(2,2)
               DO 580 J=fj,lj
                  IF(FZO(j).EQ.1)GO TO 580
                  JSYM=ORBSYM(J)
                  BSYM=IEOR(JSYM,BEUSYM)+1
                  FB=FLOV(BSYM,2,3)-NO
                  LB=FLOV(BSYM,2,4)-NO
                  DO 570 B=FB,LB
                     sgn=1.0d0
                     if(j.gt.u)sgn=-sgn
                     if(be.gt.b)sgn=-sgn
                     if(j.eq.u.or.be.eq.b)sgn=0.0d0
                     UJBBE=TOFF(U,J)+TADD(B,BE)
                     XBJ=XBJ+1
                     SUX1(XBJ)=0.5D0*T2O(UJBBE)*sgn+T1O(U,B)*T1O(J,BE)
  570             CONTINUE
  580          CONTINUE
               limbj=XBJ
c
               DO 642 AI=1,limai
                  VALJ=0.0D0
                  DO 640 BJ=1,limbj
                     VALJ=VALJ+AR1(BJ,AI)*SUX1(BJ)
  640             CONTINUE
                  JAR(AI)=JAR(AI)+VALJ
  642          CONTINUE
C
c             write(6,*)'WR2,tsym,limai,pon',tsym,limai,pon
               CALL WWRITW(ITAP91,JAR,INTOWP(limai),PON,PON)
C
  680       CONTINUE
  690    CONTINUE
  695 CONTINUE
c
  700 continue
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
               if(gaspn.ne.spin2)go to 770
               GASYM=ORBSYM(GA+NO)
               VGASYM=IEOR(VSYM,GASYM)
               IF(VGASYM.NE.TSYM-1)GO TO 770
               XVG=XVG+1
               XAI=0
               fi=ssov(spin3,1)
               li=ssov(spin3,2)
               DO 750 I=fi,li
                  IF(FZO(i).EQ.1)GO TO 750
                  ISYM=ORBSYM(I)
                  ASYM=IEOR(ISYM,VGASYM)+1
                  FA=FLOV(ASYM,spin3,3)-NO
                  LA=FLOV(ASYM,spin3,4)-NO
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
c        LIM=XVG
         limai=xai
C
         DO 990 U=1,NO
            IF(FZO(U).EQ.1)GO TO 990
            USYM=ORBSYM(U)
            uspn=spin(u)
            if(uspn.ne.spin1)go to 990
            DO 980 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 980
               bespn=spin(be+no)
               if(bespn.ne.spin1)go to 980
               BESYM=ORBSYM(BE+NO)
               BEUSYM=IEOR(BESYM,USYM)
               IF(BEUSYM.NE.TSYM-1)GO TO 980
c             write(6,*)'RD2,tsym,limai,pon',tsym,limai,pon
               CALL WREADW(ITAP91,JAR,INTOWP(limai),PON,PON)
C
               VGA=0
               DO 860 V=1,no
                  IF(FZO(V).EQ.1)GO TO 860
                  VSYM=ORBSYM(V)
                  GASYM=IEOR(VSYM,BEUSYM)+1
                  vspn=spin(v)
                  if(vspn.ne.spin2)go to 860
                  FGA=FLOV(GASYM,spin2,3)-NO
                  LGA=FLOV(GASYM,spin2,4)-NO
                  DO 850 GA=FGA,LGA
                     VGA=VGA+1
                  if(u.ne.v.and.be.ne.ga)then     
                     VAL1=0.0D0
                  DO 810 AI=1,limai
                     VAL1=VAL1+JAR(AI)*AR1(AI,VGA)
c             write(6,678)map(u),map(v),map(be+no),map(ga+no),spin(u),
c    .       spin(v), spin(be+no),spin(ga+no),jar(ai),ar1(ai,vga),val1
  810             CONTINUE
                  UVBEGA=TOFF(U,V)+TADD(BE,GA)
         if((u.lt.v.and.be.gt.ga).or.(u.gt.v.and.be.lt.ga))val1=-val1
c                 write(6,*)'val1=',val1
                  T2(UVBEGA)=T2(UVBEGA)+val1
  678      format(4i3,3x,4i2,3x,3(f12.8))
                  endif
  850          CONTINUE
  860          CONTINUE
C
  980       CONTINUE
  990    CONTINUE
  995 CONTINUE
c
  998 CONTINUE
      CALL RCLOSE(ITAP91,3)
c
      RETURN
      END
