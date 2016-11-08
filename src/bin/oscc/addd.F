
      SUBROUTINE ADDD  (T2,DIN,XVAL,
     .                  ORBSYM,FLOV,DOFF,DADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,
     .                  NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  T2(*),DIN(*),XVAL(*)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),
     .        DOFF(*),DADD(*),map(*),spin(*),ssov(2,4)
C
      DO 191 U=1,NO
         IF(FZO(U).EQ.1)GO TO 191
         USYM=ORBSYM(U)
         uspn=spin(u)
         xu=map(u)
         DO 181 V=1,u-1
            IF(FZO(V).EQ.1)GO TO 181
            VSYM=ORBSYM(V)
            VUSYM=IEOR(VSYM,USYM)
            vspn=spin(v)
            xv=map(v)
            fbe=ssov(uspn,3)-no
            lbe=ssov(uspn,4)-no
            DO 171 BE=fbe,lbe
               if(fzv(be).eq.1)go to 171
               xbe=map(be+no)
               BESYM=ORBSYM(BE+NO)
               GASYM=IEOR(BESYM,VUSYM)+1
               FGA=FLOV(GASYM,vspn,3)-NO
               LGA=FLOV(GASYM,vspn,4)-NO
               BEU=ipq(max0(xbe,xu))+min0(xbe,xu)
               beV=ipq(max0(xbe,xv))+min0(xbe,xv)
               IF(LGA.Ge.BE)LGA=BE-1
               DO 161 GA=FGA,LGA
                  xga=map(ga+no)
                  GAV=ipq(max0(xga,xv))+min0(xga,xv)
                  BEUGAV=DOFF(MAX0(BEU,GAV))+DADD(MIN0(BEU,GAV))
                  uvbega=toff(u,v)+tadd(be,ga)
                  t2(uvbega)=t2(uvbega)+DIN(BEUGAV)
                  if(uspn.eq.vspn)then
                  GAu=ipq(max0(xga,xu))+min0(xga,xu)
                  gaUbeV=DOFF(MAX0(gaU,beV))+DADD(MIN0(gaU,beV))
                  t2(uvbega)=t2(uvbega)-DIN(gaubev)
                  endif
  161          CONTINUE
  171       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
      RETURN
      END
