      SUBROUTINE HOINT (T1O,T2O,DIN,FOCK,HO,XVAL,
     .                  ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .                  ZLX,NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,TMP
      REAL*8  T1O(NO,NV),T2O(*),DIN(*),FOCK(*),HO(NO,NO),XVAL(*)
      INTEGER ORBSYM(*),FLOV(8,2,4),DOFF(*),DADD(*),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        map(*),spin(*),ssov(2,4)
C
      call zero(ho,no*no)
c
      DO 132 U=1,NO
         IF(FZO(U).EQ.1)GO TO 132
         USYM=ORBSYM(U)
         uspn=spin(u)
         FI=FLOV(USYM+1,uspn,1)
         LI=FLOV(USYM+1,uspn,2)
         DO 131 I=FI,LI
            VAL=0.0D0
            FUI=IPQ(MAX0(U,I))+MIN0(U,I)
            IF(I.NE.U)VAL=VAL+FOCK(FUI)
            HO(U,I)=VAL
  131    CONTINUE
  132 CONTINUE
C
      DO 193 A=1,NV
         IF(FZV(A).EQ.1)GO TO 193
         ASYM=ORBSYM(A+NO)
         xa=map(a+no)
         aspn=spin(a+no)
         DO 183 B=1,a-1
            IF(FZV(B).EQ.1)GO TO 183
            xb=map(b+no)
            bspn=spin(b+no)
            BSYM=ORBSYM(B+NO)
            BASYM=IEOR(BSYM,ASYM)
            DO 173 J=1,NO
               IF(FZO(J).EQ.1)GO TO 173
               xj=map(j)
               jspn=spin(j)
               if(aspn.eq.bspn.and.jspn.ne.aspn)go to 173
               JSYM=ORBSYM(J)
               USYM=IEOR(JSYM,BASYM)+1
               if(bspn.eq.jspn)uspn=aspn
               if(aspn.eq.jspn)uspn=bspn
               FU=FLOV(USYM,uspn,1)
               LU=FLOV(USYM,uspn,2)
               JA=ipq(max0(xj,xa))+min0(xj,xa)
               JB=ipq(max0(xj,xb))+min0(xj,xb)
               fu2=fu
               if(fu.le.j)fu2=j+1
               DO 143 U=FU2,LU
                  UJAB=TOFF(U,J)+TADD(A,B)
                  XVAL(U)=T2O(UJAB)+T1O(U,A)*T1O(J,B)
                  if(aspn.eq.bspn)xval(u)=xval(u)-t1o(u,b)*t1o(j,a)
  143          CONTINUE
               DO 163 I=FU,LU
                  xi=map(i)
                  IA=ipq(max0(xi,xa))+min0(xi,xa)
                  IB=ipq(max0(xi,xb))+min0(xi,xb)
                  IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                  IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                  TMP=DIN(IAJB)
                  if(aspn.eq.bspn)tmp=tmp-din(ibja)
                  DO 153 U=FU2,LU
                     HO(U,I)=HO(U,I)+TMP*XVAL(U)
  153             CONTINUE
  163          CONTINUE
               lu2=lu
               if(lu.ge.j)lu2=j-1
               DO 144 U=FU,LU2
                  juAB=TOFF(j,u)+TADD(A,B)
                  XVAL(U)=T2O(juAB)+T1O(j,A)*T1O(u,B)
                  if(aspn.eq.bspn)xval(u)=xval(u)-t1o(j,b)*t1o(u,a)
  144          CONTINUE
               DO 164 I=FU,LU
                  xi=map(i)
                  IA=ipq(max0(xi,xa))+min0(xi,xa)
                  IB=ipq(max0(xi,xb))+min0(xi,xb)
                  IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
                  IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                  TMP=DIN(ibja)
                  if(aspn.eq.bspn)tmp=tmp-din(iajb)
                  DO 154 U=FU,LU2
                     HO(U,I)=HO(U,I)+TMP*XVAL(U)
  154             CONTINUE
  164          CONTINUE
c
  173       CONTINUE
  183    CONTINUE
  193 CONTINUE
C
      RETURN
      END
