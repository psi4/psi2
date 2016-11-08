      SUBROUTINE FT2   (T1,T1O,T2O,FIN,AUX4,XVAL,
     .                  ORBSYM,FLOV,FOFF,FADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL1,VAL2,TMP,TMV
      REAL*8  T1(NO,NV),T1O(NO,NV),T2O(*),FIN(*),XVAL(*),AUX4(NO,NV)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        FOFF(*),FADD(*),map(*),spin(*),ssov(2,4)
C
      CALL ZERO(XVAL,NO)
C
      DO 191 I=1,NO
         IF(FZO(I).EQ.1)GO TO 191
         ISYM=ORBSYM(I)
         ispn=spin(i)
         xi=map(i)
         DO 181 A=1,nv 
            if(fzv(a).eq.1)go to 181
            xa=map(a+no)
            aspn=spin(a+no)
            ASYM=ORBSYM(A+NO)
            AISYM=IEOR(ASYM,ISYM)
            ia=ipq(max0(xi,xa))+min0(xi,xa)
            DO 171 B=1,a-1
               IF(FZV(B).EQ.1)GO TO 171
               bspn=spin(b+no)
               if(bspn.eq.aspn.and.bspn.ne.ispn)go to 171
               xb=map(b+no)
               ib=ipq(max0(xi,xb))+min0(xi,xb)
               BSYM=ORBSYM(B+NO)
               USYM=IEOR(BSYM,AISYM)+1
               if(ispn.eq.aspn)xspn=bspn
               if(ispn.eq.bspn)xspn=aspn
               FBE=FLOV(USYM,xspn,3)-NO
               LBE=FLOV(USYM,xspn,4)-NO
               FU=FLOV(USYM,xspn,1)
               LU=FLOV(USYM,xspn,2)
               lu2=lu
               if(lu.ge.i)lu2=i-1
               DO 111 U=FU,LU2
                  IUAB=TOFF(I,U)+TADD(A,B)
                  XVAL(U)=+T2O(IUAB)+T1O(I,A)*T1O(U,B)
                  if(aspn.eq.bspn)xval(u)=xval(u)-t1o(i,b)*t1o(u,a)
  111          CONTINUE
               fu2=fu
               if(fu.le.i)fu2=i+1
               DO 112 U=FU2,LU
                  uiAB=TOFF(u,i)+TADD(A,B)
                  XVAL(U)=+T2O(uiAB)+T1O(u,A)*T1O(i,B)
                  if(aspn.eq.bspn)xval(u)=xval(u)-t1o(u,b)*t1o(i,a)
  112          CONTINUE
               DO 161 BE=FBE,LBE
                  xbe=map(be+no)
                  BEB=ipq(MAX0(xBE,xB))+MIN0(xBE,xB)
                  IABEB=FOFF(max0(BEB,ia))+FADD(min0(beb,IA))
                  val1=fin(iabeb)
                  if(aspn.eq.bspn)then
                  BEa=ipq(MAX0(xBE,xa))+MIN0(xBE,xa)
                  IbBEa=FOFF(max0(BEa,ib))+FADD(min0(bea,Ib))
                  val1=val1-fin(ibbea)
                  endif
                  DO 151 U=FU,LU2
                     TMV=AUX4(U,BE)
                     TMV=TMV+val1*XVAL(U)
                     AUX4(U,BE)=TMV
  151             CONTINUE
  161          CONTINUE
               DO 162 BE=FBE,LBE
                  xbe=map(be+no)
                  BEa=ipq(MAX0(xBE,xa))+MIN0(xBE,xa)
                  IbBEa=FOFF(max0(BEa,ib))+FADD(min0(bea,Ib))
                  val1=fin(ibbea)
                  if(aspn.eq.bspn)then
                  BEb=ipq(MAX0(xBE,xb))+MIN0(xBE,xb)
                  IaBEb=FOFF(max0(BEb,ia))+FADD(min0(beb,Ia))
                  val1=val1-fin(iabeb)
                  endif
                  DO 152 U=FU2,LU
                     TMV=AUX4(U,BE)
                     TMV=TMV+val1*XVAL(U)
                     AUX4(U,BE)=TMV
  152             CONTINUE
  162          CONTINUE
  171       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
C     ADD INTO T1
C
      DO 291 U=1,NO
         IF(FZO(U).EQ.1)GO TO 291
         USYM=ORBSYM(U)
         uspn=spin(u)	
         FBE=FLOV(USYM+1,uspn,3)-NO
         LBE=FLOV(USYM+1,uspn,4)-NO
         DO 271 BE=FBE,LBE
            T1(U,BE)=T1(U,BE)+AUX4(U,BE)
c           write(6,675)uspn,map(u),map(be+no),t1(u,be)
  675  format(i3,2i6,f20.12)
  271    CONTINUE
  291 CONTINUE
C
      RETURN
      END
