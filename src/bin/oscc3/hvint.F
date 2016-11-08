      SUBROUTINE HVINT (T1O,T2O,DIN,FOCK,HV,XVAL,
     .                  ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .                  ZLX,NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,TMP
      REAL*8  T1O(NO,NV),T2O(*),DIN(*),FOCK(*),HV(NV,NV),XVAL(*)
      INTEGER ORBSYM(*),FLOV(8,2,4),DOFF(*),DADD(*),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        map(*),spin(*),ssov(2,4)
C
      call zero(hv,nv*nv)
c
      DO 132 BE=1,NV
         IF(FZV(BE).EQ.1)GO TO 132
         bespn=spin(be+no)
         BESYM=ORBSYM(BE+NO)
         FA=FLOV(BESYM+1,bespn,3)-NO
         LA=FLOV(BESYM+1,bespn,4)-NO
         DO 131 A=FA,LA
            VAL=0.0D0
            FABE=IPQ(MAX0(A,BE)+NO)+MIN0(A,BE)+NO
            IF(A.NE.BE)VAL=VAL+FOCK(FABE)
            HV(BE,A)=VAL
  131    CONTINUE
  132 CONTINUE
C
      DO 139 B=1,NV
         IF(FZV(B).EQ.1)GO TO 139
         bspn=spin(b+no)
         xb=map(b+no)
         BSYM=ORBSYM(B+NO)
         DO 138 i=1,NO
            IF(FZO(i).EQ.1)GO TO 138
            ispn=spin(i)
            xi=map(i)
            iSYM=ORBSYM(i)
            iBSYM=IEOR(iSYM,BSYM)
            iB=ipq(max0(xi,xb))+min0(xi,xb)
            DO 137 j=1,i-1
               IF(FZO(j).EQ.1)GO TO 137
               jspn=spin(j)
               if(ispn.eq.jspn.and.ispn.ne.bspn)go to 137
               xj=map(j)
               jSYM=ORBSYM(j)
               BESYM=IEOR(jSYM,iBSYM)+1
               if(bspn.eq.jspn)bespn=ispn
               if(bspn.eq.ispn)bespn=jspn
               FBE=FLOV(BESYM,bespn,3)-NO
               LBE=FLOV(BESYM,bespn,4)-NO
               jB=ipq(max0(xj,xb))+min0(xj,xb)
               fbe2=fbe
               if(fbe.le.b)fbe2=b+1
               DO 133 BE=FBE2,LBE
                  IJBEB=TOFF(I,J)+TADD(BE,B)
                  XVAL(BE)=T2O(IJBEB)+T1O(I,BE)*T1O(J,B)
                  if(ispn.eq.jspn)xval(be)=xval(be)-t1o(i,b)*t1o(j,be)
  133          CONTINUE
               DO 135 A=FBE,LBE
                  xa=map(a+no)
                  JA=ipq(max0(xj,xa))+min0(xj,xa)
                  IA=ipq(max0(xi,xa))+min0(xi,xa)
                  IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                  IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                  TMP=DIN(IAJB)
                  if(ispn.eq.jspn)tmp=tmp-din(ibja)
                  DO 134 BE=FBE2,LBE
                     HV(BE,A)=HV(BE,A)-TMP*XVAL(BE)
  134             CONTINUE
  135          CONTINUE
               lbe2=lbe
               if(lbe.ge.b)lbe2=b-1
               DO 143 BE=FBE,LBE2
                  IJbbe=TOFF(I,J)+TADD(B,Be)
                  XVAL(BE)=T2O(IJBBe)+T1O(I,B)*T1O(J,Be)
                  if(ispn.eq.jspn)xval(be)=xval(be)-t1o(i,be)*t1o(j,b)
  143          CONTINUE
               DO 145 A=FBE,LBE
                  xa=map(a+no)
                  JA=ipq(max0(xj,xa))+min0(xj,xa)
                  IA=ipq(max0(xi,xa))+min0(xi,xa)
                  IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                  IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                  TMP=DIN(IbJa)
                  if(ispn.eq.jspn)tmp=tmp-din(iajb)
                  DO 144 BE=FBE,LBE2
                     HV(BE,A)=HV(BE,A)-TMP*XVAL(BE)
  144             CONTINUE
  145          CONTINUE
  137       CONTINUE
  138    CONTINUE
  139 CONTINUE
C
      RETURN
      END
