      SUBROUTINE CDT1  (T1O,T2,CIN,DIN,XVAL,
     .                  ORBSYM,FLOV,COFF,CADD,DOFF,DADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,tmp,sgn
      REAL*8  T1O(NO,NV),T2(*),CIN(*),DIN(*),XVAL(*)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        COFF(*),CADD(*),DOFF(*),DADD(*),
     .        map(*),spin(*),ssov(2,4)
C
      DO 191 U=1,NO
         IF(FZO(U).EQ.1)GO TO 191
         USYM=ORBSYM(U)
         uspn=spin(u)
         xu=map(u)
         DO 181 V=1,u-1
            IF(FZO(V).EQ.1)GO TO 181
            vspn=spin(v)
            xv=map(v)
            VSYM=ORBSYM(V)
            FA=FLOV(VSYM+1,vspn,3)-NO
            LA=FLOV(VSYM+1,vspn,4)-NO
            VUSYM=IEOR(VSYM,USYM)
            fbe=ssov(uspn,3)-no
            lbe=ssov(uspn,4)-no
            DO 171 BE=fbe,lbe
               if(fzv(be).eq.1)go to 171
               xbe=map(be+no)
               BESYM=ORBSYM(BE+NO)
               ISYM=IEOR(BESYM,VUSYM)+1
               FI=FLOV(ISYM,vspn,1)
               LI=FLOV(ISYM,vspn,2)
               FGA=FLOV(ISYM,vspn,3)-NO
               LGA=FLOV(ISYM,vspn,4)-NO
               if(lga.ge.be)lga=be-1
               BEU=ipq(max0(xbe,xu))+min0(xbe,xu)
               DO 134 I=FI,LI
                  xi=map(i)
                  iU=ipq(max0(xi,xu))+min0(xi,xu)
                  XVAL(I)=0.0D0
                  DO 131 A=FA,LA
                     xa=map(a+no)
                     IA=ipq(max0(xi,xa))+min0(xi,xa)
                     beA=ipq(max0(xbe,xa))+min0(xbe,xa)
                     BEUIA=DOFF(MAX0(BEU,IA))+DADD(MIN0(BEU,IA))
                     BEaIu=cOFF(MAX0(BEa,Iu))+cADD(MIN0(BEa,Iu))
                     tmp=din(beuia)
                     if(uspn.eq.vspn)tmp=tmp-cin(beaiu)
                     XVAL(I)=XVAL(I)-tmp*T1O(V,A)
  131             CONTINUE
  134          CONTINUE
               DO 161 GA=FGA,LGA
                  VAL=0.0D0
                  DO 143 I=FI,LI
                     VAL=VAL+XVAL(I)*T1O(I,GA)
  143             CONTINUE
                  UVBEGA=TOFF(U,v)+TADD(be,ga)
                  T2(UVBEGA)=T2(UVBEGA)+VAL
  161          CONTINUE
  171       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
      DO 291 U=1,NO
         IF(FZO(U).EQ.1)GO TO 291
         USYM=ORBSYM(U)
         uspn=spin(u)
         xu=map(u)
         DO 281 V=u+1,no
            IF(FZO(V).EQ.1)GO TO 281
            vspn=spin(v)
            xv=map(v)
            VSYM=ORBSYM(V)
            FA=FLOV(VSYM+1,vspn,3)-NO
            LA=FLOV(VSYM+1,vspn,4)-NO
            VUSYM=IEOR(VSYM,USYM)
            fbe=ssov(uspn,3)-no
            lbe=ssov(uspn,4)-no
            DO 271 BE=fbe,lbe
               if(fzv(be).eq.1)go to 271
               xbe=map(be+no)
               BESYM=ORBSYM(BE+NO)
               ISYM=IEOR(BESYM,VUSYM)+1
               FI=FLOV(ISYM,vspn,1)
               LI=FLOV(ISYM,vspn,2)
               FGA=FLOV(ISYM,vspn,3)-NO
               LGA=FLOV(ISYM,vspn,4)-NO
               if(fga.le.be)fga=be+1
               BEU=ipq(max0(xbe,xu))+min0(xbe,xu)
               DO 234 I=FI,LI
                  xi=map(i)
                  iU=ipq(max0(xi,xu))+min0(xi,xu)
                  XVAL(I)=0.0D0
                  DO 231 A=FA,LA
                     xa=map(a+no)
                     IA=ipq(max0(xi,xa))+min0(xi,xa)
                     beA=ipq(max0(xbe,xa))+min0(xbe,xa)
                     BEUIA=DOFF(MAX0(BEU,IA))+DADD(MIN0(BEU,IA))
                     BEaIu=cOFF(MAX0(BEa,Iu))+cADD(MIN0(BEa,Iu))
                     tmp=din(beuia)
                     if(uspn.eq.vspn)tmp=tmp-cin(beaiu)
                     XVAL(I)=XVAL(I)-tmp*T1O(V,A)
  231             CONTINUE
  234          CONTINUE
               DO 261 GA=FGA,LGA
                  VAL=0.0D0
                  DO 243 I=FI,LI
                     VAL=VAL+XVAL(I)*T1O(I,GA)
  243             CONTINUE
c                 UVBEGA=TOFF(v,u)+TADD(ga,be)
                  UVBEGA=TOFF(U,v)+TADD(be,ga)
                  T2(UVBEGA)=T2(UVBEGA)+VAL
  261          CONTINUE
  271       CONTINUE
  281    CONTINUE
  291 CONTINUE
c
      DO 391 U=1,NO
         IF(FZO(U).EQ.1)GO TO 391
         USYM=ORBSYM(U)
         uspn=spin(u)
         xu=map(u)
         DO 381 V=1,u-1
            IF(FZO(V).EQ.1)GO TO 381
            vspn=spin(v)
            xv=map(v)
            VSYM=ORBSYM(V)
            FA=FLOV(VSYM+1,vspn,3)-NO
            LA=FLOV(VSYM+1,vspn,4)-NO
            VUSYM=IEOR(VSYM,USYM)
            fbe=ssov(vspn,3)-no
            lbe=ssov(vspn,4)-no
            DO 371 BE=fbe,lbe
               if(fzv(be).eq.1)go to 371
               xbe=map(be+no)
               BESYM=ORBSYM(BE+NO)
               ISYM=IEOR(BESYM,VUSYM)+1
               FI=FLOV(ISYM,uspn,1)
               LI=FLOV(ISYM,uspn,2)
               FGA=FLOV(ISYM,uspn,3)-NO
               LGA=FLOV(ISYM,uspn,4)-NO
               if(fga.le.be)fga=be+1
               BEU=ipq(max0(xbe,xu))+min0(xbe,xu)
               DO 334 I=FI,LI
                  xi=map(i)
                  iU=ipq(max0(xi,xu))+min0(xi,xu)
                  XVAL(I)=0.0D0
                  DO 331 A=FA,LA
                     xa=map(a+no)
                     IA=ipq(max0(xi,xa))+min0(xi,xa)
                     beA=ipq(max0(xbe,xa))+min0(xbe,xa)
                     BEUIA=DOFF(MAX0(BEU,IA))+DADD(MIN0(BEU,IA))
                     BEaIu=cOFF(MAX0(BEa,Iu))+cADD(MIN0(BEa,Iu))
                     tmp=cin(beaiu)
                     if(uspn.eq.vspn)tmp=tmp-din(beuia)
                     XVAL(I)=XVAL(I)-tmp*T1O(V,A)
  331             CONTINUE
  334          CONTINUE
               DO 361 GA=FGA,LGA
                  VAL=0.0D0
                  DO 343 I=FI,LI
                     VAL=VAL+XVAL(I)*T1O(I,GA)
  343             CONTINUE
c                 UVBEGA=TOFF(u,v)+TADD(ga,be)
                  UVBEGA=TOFF(U,v)+TADD(be,ga)
                  T2(UVBEGA)=T2(UVBEGA)+VAL
  361          CONTINUE
  371       CONTINUE
  381    CONTINUE
  391 CONTINUE
c
      DO 491 U=1,NO
         IF(FZO(U).EQ.1)GO TO 491
         USYM=ORBSYM(U)
         uspn=spin(u)
         xu=map(u)
         DO 481 V=u+1,no
            IF(FZO(V).EQ.1)GO TO 481
            vspn=spin(v)
            xv=map(v)
            VSYM=ORBSYM(V)
            FA=FLOV(VSYM+1,vspn,3)-NO
            LA=FLOV(VSYM+1,vspn,4)-NO
            VUSYM=IEOR(VSYM,USYM)
            fbe=ssov(vspn,3)-no
            lbe=ssov(vspn,4)-no
            DO 471 BE=fbe,lbe
               if(fzv(be).eq.1)go to 471
               xbe=map(be+no)
               BESYM=ORBSYM(BE+NO)
               ISYM=IEOR(BESYM,VUSYM)+1
               FI=FLOV(ISYM,uspn,1)
               LI=FLOV(ISYM,uspn,2)
               FGA=FLOV(ISYM,uspn,3)-NO
               LGA=FLOV(ISYM,uspn,4)-NO
               if(lga.ge.be)lga=be-1
               BEU=ipq(max0(xbe,xu))+min0(xbe,xu)
               DO 434 I=FI,LI
                  xi=map(i)
                  iU=ipq(max0(xi,xu))+min0(xi,xu)
                  XVAL(I)=0.0D0
                  DO 431 A=FA,LA
                     xa=map(a+no)
                     IA=ipq(max0(xi,xa))+min0(xi,xa)
                     beA=ipq(max0(xbe,xa))+min0(xbe,xa)
                     BEUIA=DOFF(MAX0(BEU,IA))+DADD(MIN0(BEU,IA))
                     BEaIu=cOFF(MAX0(BEa,Iu))+cADD(MIN0(BEa,Iu))
                     tmp=cin(beaiu)
                     if(uspn.eq.vspn)tmp=tmp-din(beuia)
                     XVAL(I)=XVAL(I)-tmp*T1O(V,A)
  431             CONTINUE
  434          CONTINUE
               DO 461 GA=FGA,LGA
                  VAL=0.0D0
                  DO 443 I=FI,LI
                     VAL=VAL+XVAL(I)*T1O(I,GA)
  443             CONTINUE
c                 UVBEGA=TOFF(v,u)+TADD(be,ga)
                  UVBEGA=TOFF(U,v)+TADD(be,ga)
                  T2(UVBEGA)=T2(UVBEGA)+VAL
  461          CONTINUE
  471       CONTINUE
  481    CONTINUE
  491 CONTINUE
c
      RETURN
      END
