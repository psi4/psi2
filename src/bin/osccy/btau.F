      SUBROUTINE BTAU  (T1O,T2O,T2,xIN,AR3,BIN,BUFINT,
     .                  SUX1,SUX3,
     .                  ORBSYM,FLOV,iOFF,iADD,NAB,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,
     .                  NO,NV,NIRRED,MAXOO,MAXVV,ITAP71,
     .                  map,spin,ssov,ndimt2)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL1,VAL2,TMP
      REAL*8  T1O(NO,NV),T2O(*),T2(ndimt2),xIN(*),BIN(NV,NV),BUFINT(*),
     .        SUX1(*),SUX3(*),AR3(MAXVV,MAXOO)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),
     .        iOFF(*),iADD(*),NAB(*),
     .        map(*),spin(*),ssov(2,4),xspin(3),yspin(3)
c 
      xspin(1)=1
      yspin(1)=1
      xspin(2)=2
      yspin(2)=1
      xspin(3)=2
      yspin(3)=2
c        
      do 999 tsym=1,nirred
      do 998 case=1,3
      spin1=xspin(case)
      spin2=yspin(case)
c
         XAB=0
         fa=ssov(spin1,3)-no
         la=ssov(spin1,4)-no
         DO 290 A=fa,la
            if(fzv(a).eq.1)go to 290
            ASYM=ORBSYM(A+NO)
            bsym=ieor(asym,tsym-1)
            fb=flov(bsym+1,spin2,3)-no
            lb=flov(bsym+1,spin2,4)-no
            if(lb.ge.a)lb=a-1
            DO 280 B=fb,lb
               XAB=XAB+1
               XIJ=0
               fi=ssov(spin1,1)
               li=ssov(spin1,2)
               DO 270 I=fi,li
                  if(fzo(i).eq.1)go to 270
                  ISYM=ORBSYM(I)
                  JSYM=IEOR(ISYM,tsym-1)
                  FJ=FLOV(JSYM+1,spin2,1)
                  LJ=FLOV(JSYM+1,spin2,2)
                  IF(LJ.Ge.I)LJ=I-1
                  DO 260 J=FJ,LJ
                     XIJ=XIJ+1
                     IJAB=TOFF(I,J)+TADD(A,B)
                     VAL1=T2O(IJAB)+T1O(I,A)*T1O(J,B)
                     if(case.ne.2)val1=val1-t1o(i,b)*t1o(j,a)
                     AR3(XAB,XIJ)=VAL1
  260             CONTINUE
  270          CONTINUE
  280       CONTINUE
  290    CONTINUE
         LIMO=XIJ
         LIMV=XAB
C
c        WRITE(6,*)'TSYM=',TSYM,' XIJ=',XIJ,' XAB=',XAB,'   ',XIJ*XAB
C
         DO 991 BE=1,NV
            IF(FZV(BE).EQ.1)GO TO 991
            BESYM=ORBSYM(BE+NO)
            bespn=spin(be+no)
            if(bespn.ne.spin1)go to 991
            xbe=map(be+no)
            FI1=FLOV(BESYM+1,bespn,1)
            LI1=FLOV(BESYM+1,bespn,2)
            DO 981 GA=1,BE-1
               IF(FZV(GA).EQ.1)GO TO 981
               GASYM=ORBSYM(GA+NO)
               gaspn=spin(ga+no)
               if(gaspn.ne.spin2)go to 981
               xga=map(ga+no)
               FI2=FLOV(GASYM+1,gaspn,1)
               LI2=FLOV(GASYM+1,gaspn,2)
               BEGAS=IEOR(GASYM,BESYM)
               IF(BEGAS.NE.TSYM-1)GO TO 981
C
               fa=ssov(bespn,3)-no
               la=ssov(bespn,4)-no
               xba=0
               DO 425 A=fa,la
                  if(fzv(a).eq.1)go to 425
                  ASYM=ORBSYM(A+NO)
                  xa=map(a+no)
                  BEA=ipq(MAX0(xBE,xA))+MIN0(xBE,xA)
                  GAA=ipq(MAX0(xGA,xA))+MIN0(xGA,xA)
                  BSYM=IEOR(ASYM,BEGAS)
                  FB=FLOV(BSYM+1,gaspn,3)-NO
                  LB=FLOV(BSYM+1,gaspn,4)-NO
                  IF(LB.Ge.A)LB=A-1
                  DO 420 B=FB,LB
                     xb=map(b+no)
                     XBA=XBA+1
                     GAB=ipq(MAX0(xGA,xB))+MIN0(xGA,xB)
                     BEB=ipq(MAX0(xBE,xB))+MIN0(xBE,xB)
                     beagab=ioff(max0(bea,gab))+iadd(min0(bea,gab))
                     TMP=xin(beagab)
                     if(case.ne.2)then
                        bebgaa=ioff(max0(beb,gaa))+iadd(min0(beb,gaa))
                        TMP=TMP-xin(bebgaa)
                     endif
                     VAL1=0.0D0
                     DO 410 I=FI2,LI2
                        xi=map(i)
                        IB=ipq(max0(xi,xb))+min0(xi,xb)
                        ibbea=ioff(max0(ib,bea))+iadd(min0(ib,bea))
                        val2=xin(ibbea)
                        if(case.ne.2)then
                        Ia=ipq(max0(xi,xa))+min0(xi,xa)
                        iabeb=ioff(max0(ia,beb))+iadd(min0(ia,beb))
                        val2=val2-xin(iabeb)
                        endif
                        VAL1=VAL1-val2*T1O(I,GA)
  410                CONTINUE
                     DO 408 I=FI1,LI1
                        xi=map(i)
                        Ia=ipq(max0(xi,xa))+min0(xi,xa)
                        iagab=ioff(max0(ia,gab))+iadd(min0(ia,gab))
                        val2=xin(iagab)
                        if(case.ne.2)then
                        IB=ipq(max0(xi,xb))+min0(xi,xb)
                        ibgaa=ioff(max0(ib,gaa))+iadd(min0(ib,gaa))
                        val2=val2-xin(ibgaa)
                        endif
                        VAL1=VAL1-val2*T1O(I,BE)
  408                CONTINUE
                     sux3(xba)=TMP+VAL1
  420             CONTINUE
  425          CONTINUE
C
C              GO FOR THE N**6 LOOP
C
               DO 577 UV=1,LIMO
                  TMP = 0.0D0
                  DO 529 AB=1,LIMV
                     TMP=TMP+SUX3(AB)* AR3(AB,UV)
  529             CONTINUE
                  SUX1(UV)=TMP
  577          CONTINUE
C
               UV=0
               DO 668 U=1,NO
                  IF(FZO(U).EQ.1)GO TO 668
                  USYM=ORBSYM(U)
                  VSYM=IEOR(USYM,BEGAS)+1
                  uspn=spin(u)
                  if(uspn.ne.spin1)go to 668
                  SV=FLOV(VSYM,spin2,1)
                  LV=FLOV(VSYM,spin2,2)
                  IF(LV.Ge.U )LV=U-1
                  DO 663 V=SV,LV
                     UV=UV+1
                     UVBEGA=TOFF(U,V)+TADD(BE,GA)
                     T2(UVBEGA)=T2(UVBEGA)+SUX1(UV)
  663             CONTINUE
  668          CONTINUE
C
  981       CONTINUE
  991    CONTINUE
c
  998 continue
  999 CONTINUE
C
      RETURN
      END
