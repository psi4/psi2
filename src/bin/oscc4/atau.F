      SUBROUTINE ATAU  (T1O,T2O,T2,AIN,DIN,EIN,AR3,AR4,AR5,AR6,
     .                  SUX1,SFOO,SUX3,SFVV,AUX1,
     .                  ORBSYM,FLOV,AOFF,AADD,DOFF,DADD,EOFF,EADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,MAXOO,MAXVV,ITAP91,
     .                  map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL1,VAL2,TMP,TMV,VALP,VALM
      REAL*8  T1O(NO,NV),T2O(*),T2(*),AIN(*),DIN(*),EIN(*),
     .        SUX1(*),SFOO(*),SUX3(*),SFVV(*),AUX1(NV,NO),
     .        AR3(MAXVV,MAXOO),AR4(MAXVV,MAXOO),
     .        AR5(MAXOO,MAXVV),AR6(MAXOO,MAXVV)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        AOFF(*),AADD(*),DOFF(*),DADD(*),EOFF(*),EADD(*),
     .        map(*),spin(*),ssov(2,4),xspin(3),yspin(3)
C
      CALL TRNPS(AUX1,T1O,NO,NV)
C
      xspin(1)=1
      yspin(1)=1
      xspin(2)=2
      yspin(2)=1
      xspin(3)=2
      yspin(3)=2
c        
c     write(6,*)'in atau'
c     write(6,*)'intowp(0)=',intowp(0)
      do 999 case=1,3
      spin1=xspin(case)
      spin2=yspin(case)
c     write(6,*)'CASE, spins',case,spin1,spin2
c
      CALL RFILE(ITAP91)
      PON=1
      DO 393 TSYM=1,NIRRED
         XJI=0
         fi=ssov(spin1,1)
         li=ssov(spin1,2)
         DO 190 I=fi,li
            if(fzo(i).eq.1)go to 190
            ISYM=ORBSYM(I)
            xi=map(i)
            fj=ssov(spin2,1)
            lj=ssov(spin2,2)
            if(lj.ge.i)lj=i-1
            DO 180 J=fj,lj
               if(fzo(j).eq.1)go to 180
               xj=map(j)
               JSYM=ORBSYM(J)
               JISYM=IEOR(JSYM,ISYM)
               IF(JISYM.NE.TSYM-1)GO TO 180
               XJI=XJI+1
               XBA=0
               fa=ssov(spin1,3)-no
               la=ssov(spin1,4)-no
               DO 170 A=fa,la
                  if(fzv(a).eq.1)go to 170
                  xa=map(a+no)
                  IA=ipq(max0(xi,xa))+min0(xi,xa)
                  jA=ipq(max0(xj,xa))+min0(xj,xa)
                  ASYM=ORBSYM(A+NO)
                  BSYM=IEOR(ASYM,JISYM)+1
                  FB=FLOV(BSYM,spin2,3)-NO
                  LB=FLOV(BSYM,spin2,4)-NO
                  IF(LB.Ge.A)LB=A-1
                  DO 160 B=FB,LB
                     xb=map(b+no)
                     XBA=XBA+1
                     jb=ipq(max0(xj,xb))+min0(xj,xb)
                     IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                     VAL1=DIN(IAJB)
                     if(case.ne.2)then
                        ib=ipq(max0(xi,xb))+min0(xi,xb)
                        IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                        val1=val1-din(ibja)
                     endif
                     AR3(XBA,XJI)=VAL1
  160             CONTINUE
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
         LIMO=XJI
         LIMV=XBA
c        write(6,*)'limo,limv',limo,limv
C
         DO 390 U=1,NO
            IF(FZO(U).EQ.1)GO TO 390
            uspn=spin(u)
            xu=map(u)
            if(uspn.ne.spin1)go to 390
            USYM=ORBSYM(U)
            FA1=FLOV(USYM+1,spin1,3)-NO
            LA1=FLOV(USYM+1,spin1,4)-NO
            DO 380 V=1,U-1
               IF(FZO(V).EQ.1)GO TO 380
               vspn=spin(v)
               xv=map(v)
               if(vspn.ne.spin2)go to 380
               VSYM=ORBSYM(V)
               FA2=FLOV(VSYM+1,spin2,3)-NO
               LA2=FLOV(VSYM+1,spin2,4)-NO
               VUSYM=IEOR(VSYM,USYM)
               IF(VUSYM.NE.TSYM-1)GO TO 380
               BA=0
               fa=ssov(spin1,3)-no
               la=ssov(spin1,4)-no
               DO 205 A=fa,la
                  if(fzv(a).eq.1)go to 205
                  ASYM=ORBSYM(A+NO)
                  BSYM=IEOR(ASYM,VUSYM)
                  FB=FLOV(BSYM+1,spin2,3)-NO
                  LB=FLOV(BSYM+1,spin2,4)-NO
                  IF(LB.Ge.A)LB=A-1
                  DO 204 B=FB,LB
                     BA=BA+1
                     UVAB=TOFF(U,V)+TADD(A,B)
                     TMP=T2O(UVAB)+T1O(U,A)*T1O(V,B)
                     if(case.ne.2)tmp=tmp-t1o(u,b)*t1o(v,a)
                     SUX3(BA)=TMP
  204             CONTINUE
  205          CONTINUE
C
               JI=0
               fi=ssov(spin1,1)
               li=ssov(spin1,2)
               DO 340 I=fi,li
                  if(fzo(i).eq.1)go to 340
                  xi=map(i)
                  ISYM=ORBSYM(I)
                  IU=ipq(MAX0(xI,xU))+MIN0(xI,xU)
                  IV=ipq(MAX0(xI,xV))+MIN0(xI,xV)
                  JSYM=IEOR(ISYM,VUSYM)
                  FJ=FLOV(JSYM+1,spin2,1)
                  LJ=FLOV(JSYM+1,spin2,2)
                  IF(LJ.Ge.I)LJ=I-1
                  DO 330 J=FJ,LJ
                     xj=map(j)
                     JI=JI+1
                     JV=ipq(MAX0(xJ,xV))+MIN0(xJ,xV)
                     IUJV=AOFF(MAX0(IU,JV))+AADD(MIN0(IU,JV))
                     TMP=AIN(IUJV)
                     if(case.ne.2)then
                        JU=ipq(MAX0(xJ,xU))+MIN0(xJ,xU)
                        IVJU=AOFF(MAX0(IV,JU))+AADD(MIN0(IV,JU))
                        tmp=tmp-ain(ivju)
                     endif
                     VAL1=0.0D0
                     DO 315 A=FA2,LA2
                        xa=map(a+no)
                        ja=ipq(max0(xj,xa))+min0(xj,xa)
                        iUjA=EOFF(max0(iu,jA))+EADD(min0(iU,ja))
                        tmv=ein(iuja)
                        if(case.ne.2)then
                           ia=ipq(max0(xi,xa))+min0(xi,xa)
                           jUiA=EOFF(max0(ju,iA))+EADD(min0(jU,ia))
                           tmv=tmv-ein(juia)
                        endif
                        VAL1=VAL1+tmv*AUX1(A,V)
  315                CONTINUE
                     DO 318 A=FA1,LA1
                        xa=map(a+no)
                        ia=ipq(max0(xi,xa))+min0(xi,xa)
                        IAJV=EOFF(max0(IA,jv))+EADD(min0(ia,JV))
                        tmv=ein(iajv)
                        if(case.ne.2)then
                           ja=ipq(max0(xj,xa))+min0(xj,xa)
                           jAiV=EOFF(max0(jA,iv))+EADD(min0(ja,iV))
                           tmv=tmv-ein(jaiv)
                        endif
                        VAL1=VAL1+tmv*AUX1(A,U)
  318                CONTINUE
                     VALP=0.0D0
                     DO 320 BA=1,LIMV
                        VALP=VALP+AR3(BA,JI)*SUX3(BA)
  320                CONTINUE
                     SFOO(JI)= TMP+VAL1+VALP
  330             CONTINUE
  340          CONTINUE
               CALL WWRITW(ITAP91,SFOO,INTOWP(limo),PON,PON)
c              write(6,*)'wr,limo,pon',limo,pon
  380       CONTINUE
  390    CONTINUE
  393 CONTINUE
C
      PON=1
      DO 892 TSYM=1,NIRRED
         XAB=0
cges     if xab remains zero (ie does not go into the inner loops)
cges     xij is undefined and causes problems
         xij=0
         fa=ssov(spin1,3)-no
         la=ssov(spin1,4)-no
         DO 490 A=fa,la
            if(fzv(a).eq.1)go to 490
            ASYM=ORBSYM(A+NO)
            fb=ssov(spin2,3)-no
            lb=ssov(spin2,4)-no
            if(lb.ge.a)lb=a-1
            DO 480 B=fb,lb
               if(fzv(b).eq.1)go to 480
               BSYM=ORBSYM(B+NO)
               ABSYM=IEOR(ASYM,BSYM)
               IF(ABSYM.NE.TSYM-1)GO TO 480
               XAB=XAB+1
               XIJ=0
               fi=ssov(spin1,1)
               li=ssov(spin1,2)
               DO 470 I=fi,li
                  if(fzo(i).eq.1)go to 470
                  ISYM=ORBSYM(I)
                  JSYM=IEOR(ISYM,ABSYM)+1
                  FJ=FLOV(JSYM,spin2,1)
                  LJ=FLOV(JSYM,spin2,2)
                  IF(LJ.Ge.I)LJ=I-1
                  DO 460 J=FJ,LJ
                     XIJ=XIJ+1
                     IJAB=TOFF(I,J)+TADD(A,B)
                     VAL1=T2O(IJAB)+T1O(I,A)*T1O(J,B)
                     if(case.ne.2)val1=val1-t1o(i,b)*t1o(j,a)
                     AR5(XIJ,XAB)=VAL1
  460             CONTINUE
  470          CONTINUE
  480       CONTINUE
  490    CONTINUE
         LIMO=XIJ
         LIMV=XAB
C
         DO 891 U=1,no
            if(fzo(u).eq.1)go to 891
            USYM=ORBSYM(U)
            uspn=spin(u)
            if(uspn.ne.spin1)go to 891
            DO 881 V=1,u-1
               if(fzo(v).eq.1)go to 881
               vspn=spin(v)
               if(vspn.ne.spin2)go to 881
               VSYM=ORBSYM(V)
               VUSYM=IEOR(VSYM,USYM)
               IF(VUSYM.NE.TSYM-1)GO TO 881
               CALL WREADW(ITAP91,SFOO,INTOWP(limo),PON,PON)
c              write(6,*)'rd,limo,pon',limo,pon
C
               DO 570 BEGA=1,LIMV
                  TMP=0.0D0
                  DO 550 IJ=1,LIMO
                     TMP=TMP+SFOO(IJ)* AR5(IJ,BEGA)
  550             CONTINUE
                  SUX3(BEGA)=TMP
  570          CONTINUE
C
               BEGA=0
               fbe=ssov(spin1,3)-no
               lbe=ssov(spin1,4)-no
               DO 674 BE=fbe,lbe
                  if(fzv(be).eq.1)go to 674
                  BESYM=ORBSYM(BE+NO)
                  GASYM=IEOR(BESYM,VUSYM)+1
                  FGA=FLOV(GASYM,spin2,3)-NO
                  LGA=FLOV(GASYM,spin2,4)-NO
                  IF(LGA.Ge.BE)LGA=BE-1
                  DO 664 GA=FGA,LGA
                     BEGA=BEGA+1
                     UVBEGA=TOFF(U,V)+TADD(BE,GA)
                     T2(UVBEGA)=T2(UVBEGA)+SUX3(BEGA)
  664             CONTINUE
  674          CONTINUE
C
  881       CONTINUE
  891    CONTINUE
  892 CONTINUE
      CALL RCLOSE(ITAP91,3)
C
  999 continue
c
      RETURN
      END
