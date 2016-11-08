      SUBROUTINE BTAU  (T1O,T2O,T2,xIN,AR3,BIN,BUFINT,
     .                  SUX1,SUX3,
     .                  ORBSYM,FLOV,iOFF,iADD,NAB,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,
     .                  NO,NV,NIRRED,MAXOO,MAXVV,ITAP71,
     .                  map,spin,ssov,ndimt2,nocs,ipont,omap)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL1,VAL2,TMP,tmr
      REAL*8  T1O(NO,NV),T2O(*),T2(ndimt2),xIN(*),BIN(NV,NV),
     .        BUFINT(nv*nv),
     .        SUX1(*),SUX3(*),AR3(MAXVV,MAXOO)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),
     .        iOFF(*),iADD(*),NAB(*),
     .        map(*),spin(*),ssov(2,4),xspin(3),yspin(3),
     .        omap(*)
      integer ipont(nv,nv)
C
C >>> NAB MUST INCLUDE FROZEN PAIRS
C     BECAUSE FILE71 CONTAINS (VV,VV) INTEGRALS FOR ALL PAIRS
C     POINTER MUST BE INCREMENTED ALSO FOR FROZEN PAIRS!
c     calculate how many ints in a beta-gama pair by using beta  spin
C
      itap71=71
      CALL RFILE(ITAP71)
C
      do 992 tsym=1,nirred
         LAB=0
         DO 123 A=1,nv 
C           IF(FZV(A).EQ.1)GO TO 123
            aspn=spin(a+no)
            xa=map(a+no)
            if(aspn.ne.2)go to 123
            ASYM=ORBSYM(A+NO)
            DO 121 B=1,nv 
C              IF(FZV(B).EQ.1)GO TO 121
               bspn=spin(b+no)
               xb=map(b+no)
               if(bspn.ne.2)go to 121
               bSYM=ORBSYM(b+NO)
               absym=ieor(asym,bsym)
               if(absym.ne.tsym-1)go to 121
               LAB=LAB+1
  121       CONTINUE
  123    CONTINUE
         NAB(TSYM)=LAB
c        WRITE(6,*)'TSYM',TSYM, ' NAB(TSYM)=',NAB(TSYM)
  992 continue
c     call flush(6)
c
c >>> construct a pointer to find the required batch of integrals in 71
c     note that reading is still done forwardly only for cases=1,3
c     however not all batches are needed
c     the C lines are important because ALL integrals are kept in 71
c
      icnt=1
      do 993 tsym=1,nirred
         lab=nab(tsym)
         do 126 a=1,nv
Cimportant  IF(FZV(A).EQ.1)GO TO 126
            aspn=spin(a+no)
            xa=map(a+no)
            if(aspn.ne.2)go to 126
            ASYM=ORBSYM(A+NO)
            DO 124 B=1,a 
Cimportant     IF(FZV(B).EQ.1)GO TO 124
               bspn=spin(b+no)
               xb=map(b+no)
               if(bspn.ne.2)go to 124
               bSYM=ORBSYM(b+NO)
               absym=ieor(asym,bsym)
               if(absym.ne.tsym-1)go to 124
ccccccccc      write(6,*)'xa,xb,tsym',xa,xb,tsym
               ipont(xa,xb)=icnt
               ipont(xb,xa)=icnt
               icnt=icnt+intowp(lab)
  124       continue
  126    continue
c           
  993 continue
C
      xspin(1)=1
      yspin(1)=1
      xspin(2)=2
      yspin(2)=1
      xspin(3)=2
      yspin(3)=2
c        
      do 998 case=1,3
      spin1=xspin(case)
      spin2=yspin(case)
      POINT=1
c
      do 999 tsym=1,nirred
c
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
            oxbe=omap(be+no)
            FI1=FLOV(BESYM+1,bespn,1)
            LI1=FLOV(BESYM+1,bespn,2)
            DO 981 GA=1,BE-1
               IF(FZV(GA).EQ.1)GO TO 981
               GASYM=ORBSYM(GA+NO)
               gaspn=spin(ga+no)
               if(gaspn.ne.spin2)go to 981
               xga=map(ga+no)
               oxga=omap(ga+no)
c
               FI2=FLOV(GASYM+1,gaspn,1)
               LI2=FLOV(GASYM+1,gaspn,2)
               BEGAS=IEOR(GASYM,BESYM)
               IF(BEGAS.NE.TSYM-1)GO TO 981
c         write(6,618)case,tsym,xbe,xga
  618    format('case=',i4,'  tsym=',i4,'  xbe=',i4,'  xga=',i4)
               LAB=NAB(TSYM)
               point=ipont(xbe,xga)
c
c              IF(FZV(BE).EQ.1.OR.FZV(GA).EQ.1.or.be.eq.ga)THEN
c                 POINT=POINT+INTOWP(LAB)
c                 GO TO 981
c              ENDIF
c
c              write(6,*)'reading 71 at',point
               CALL WREADW(ITAP71,BUFINT,INTOWP(LAB),POINT,POINT)
c              write(6,*)'bufint'
c              write(6,'(4f20.12)')(bufint(ijk),ijk=1,lab)
               CALL ZERO(BIN,nv*nv)
               NADD=0
               DO 325 A=1,nv 
Cimportant        IF(FZV(A).EQ.1)GO TO 325
                  aspn=spin(a+no)
                  if(aspn.ne.2)go to 325
                  oxa=omap(a+no)-nocs
                  ASYM=ORBSYM(A+NO)
                  DO 320 B=1,nv 
Cimportant           IF(FZV(B).EQ.1)GO TO 320
                     bspn=spin(b+no)
                     if(bspn.ne.2)go to 320
                     oxb=omap(b+no)-nocs
                     BSYM=ORBSYM(B+NO)
                     ABSYM=IEOR(ASYM,BSYM)
                     IF(ABSYM.NE.TSYM-1)GO TO 320
                     NADD=NADD+1
                     if(oxbe.ge.oxga)then
                     BIN(oxA,oxB)=BUFINT(NADD)
                     else
                     BIN(oxb,oxa)=BUFINT(NADD)
                     endif
  320             CONTINUE
  325          CONTINUE
C
               fa=ssov(bespn,3)-no
               la=ssov(bespn,4)-no
               xba=0
               DO 425 A=fa,la
                  if(fzv(a).eq.1)go to 425
                  ASYM=ORBSYM(A+NO)
                  xa=map(a+no)
                  oxa=omap(a+no)
                  BEA=ipq(MAX0(xBE,xA))+MIN0(xBE,xA)
                  GAA=ipq(MAX0(xGA,xA))+MIN0(xGA,xA)
                  BSYM=IEOR(ASYM,BEGAS)
                  FB=FLOV(BSYM+1,gaspn,3)-NO
                  LB=FLOV(BSYM+1,gaspn,4)-NO
                  IF(LB.Ge.A)LB=A-1
                  DO 420 B=FB,LB
                     xb=map(b+no)
                     oxb=omap(b+no)
                     XBA=XBA+1
                     GAB=ipq(MAX0(xGA,xB))+MIN0(xGA,xB)
                     BEB=ipq(MAX0(xBE,xB))+MIN0(xBE,xB)
c                    beagab=ioff(max0(bea,gab))+iadd(min0(bea,gab))
cc                    TMr=xin(beagab)
                     tmp=bin(oxa-nocs,oxb-nocs)
c                    write(6,'(2f20.12)')tmr,tmp
                     if(case.ne.2)then
c                       bebgaa=ioff(max0(beb,gaa))+iadd(min0(beb,gaa))
cc                       TMr=TMr-xin(bebgaa)
                        TMP=TMP-bin(oxb-nocs,oxa-nocs)
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
  999 CONTINUE
  998 continue
      CALL RCLOSE(ITAP71,3)
C
      RETURN
      END
C
