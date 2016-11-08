      SUBROUTINE HIAINT(T1O,T2O,DIN,FOCK,HIA,XVAL,
     .                  ORBSYM,FLOV,DOFF,DADD,TOFF,TADD,IPQ,ITR,FZO,FZV,
     .                  ZLX,NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,TMP
      REAL*8  T1O(NO,NV),T2O(*),DIN(*),FOCK(*),HIA(NO,NV),XVAL(*)
      INTEGER ORBSYM(*),FLOV(8,2,4),DOFF(*),DADD(*),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        map(*),spin(*),ssov(2,4)
C
      call zero(hia,no*nv)
c
      DO 137 A=1,NV
         IF(FZV(A).EQ.1)GO TO 137
         ASYM=ORBSYM(A+NO)
         aspn=spin(a+no)
         xa=map(a+no)
         FI=FLOV(ASYM+1,aspn,1)
         LI=FLOV(ASYM+1,aspn,2)
         DO 136 I=FI,LI
            xi=map(i)
            VAL=0.0D0
            IA=ipq(max0(xa,xi))+min0(xa,xi)
            FIA=IPQ(A+NO)+I
            DO 134 J=1,NO
               IF(FZO(J).EQ.1)GO TO 134
               jspn=spin(j)
               xj=map(j)
               JSYM=ORBSYM(J)
               FB=FLOV(JSYM+1,jspn,3)-NO
               LB=FLOV(JSYM+1,jspn,4)-NO
               JA=ipq(max0(xj,xa))+min0(xj,xa)
               DO 130 B=FB,LB
                  xb=map(b+no)
                  Jb=ipq(max0(xj,xb))+min0(xj,xb)
                  ib=ipq(max0(xi,xb))+min0(xi,xb)
                  IAJB=DOFF(MAX0(JB,IA))+DADD(MIN0(JB,IA))
                  IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
                  tmp=din(iajb)
                  if(aspn.eq.jspn)tmp=tmp-din(ibja)
                  VAL=VAL+tmp*T1O(J,B)
  130          CONTINUE
  134       CONTINUE
            HIA(I,A)=VAL
     .                   +FOCK(FIA)
  136    CONTINUE
  137 CONTINUE
C
      RETURN
      END
