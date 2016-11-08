      SUBROUTINE INITT (T1,T2,xIN,FOCK,
     .                  ORBSYM,FLOV,iOFF,iADD,uOFF,vADD,IPQ,ITR,FZO,FZV,
     .                  NOs,NVs,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL
      REAL*8  T1(NOs,NVs),T2(*),xIN(*),FOCK(*)
      INTEGER ORBSYM(*),FLOV(8,2,4),iOFF(*),iADD(*),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),uOFF(NOs,NOs),vADD(NVs,NVs),spin(*),map(*),
     .        ssov(2,4)
C
      CALL ZERO(T1,NOs*NVs)
C
      DO 146 I=1,NOs
         IF(FZO(I).EQ.1)GO TO 146
         ISYM=ORBSYM(I)
         ispn=spin(i)
         xi=map(i)
         IY=IPQ(I)+I
         DO 136 J=1,i-1
            IF(FZO(J).EQ.1)GO TO 136
            JSYM=ORBSYM(J)
            jspn=spin(j)
            xj=map(j)
            JISYM=IEOR(JSYM,ISYM)
            JJ=IPQ(J)+J
            fa=ssov(ispn,3)-nos
            la=ssov(ispn,4)-nos
            DO 126 A=fa,la
               IF(FZV(A).EQ.1)GO TO 126
               xa=map(a+nos)
               ASYM=ORBSYM(A+NOs)
               BSYM=IEOR(ASYM,JISYM)+1
               FB=FLOV(BSYM,jspn,3)-NOs
               LB=FLOV(BSYM,jspn,4)-NOs
               AA=IPQ(A+NOs)+A+NOs
               xia=ipq(max0(xi,xa))+min0(xi,xa)
               IF(LB.Ge.A)LB=A-1
               DO 116 B=FB,LB
                  xb=map(b+nos)
                  BY=IPQ(B+NOs)+B+NOs
                  xJB=ipq(max0(xj,xb))+min0(xj,xb)
                  IAJB=iOFF(MAX0(xIA,xJB))+iADD(MIN0(xIA,xJB))
                  tIJAB=uOFF(I,J)+vADD(A,B)
                  VAL=FOCK(IY)+FOCK(JJ)-FOCK(AA)-FOCK(BY)
                  T2(tIJAB)=xIN(IAJB)/VAL
                  if(ispn.eq.jspn)then
                  xiB=ipq(max0(xi,xb))+min0(xi,xb)
                  xJa=ipq(max0(xj,xa))+min0(xj,xa)
                  IbJa=iOFF(MAX0(xIb,xJa))+iADD(MIN0(xIb,xJa))
                  T2(tIJAB)=t2(tijab)-xIN(IbJa)/VAL
                  endif
c         write(6,678)i,j,a,b,tijab,t2(tijab)
  678     format(4i4,i8,f20.12)
  116          CONTINUE
  126       CONTINUE
  136    CONTINUE
  146 CONTINUE
C
      RETURN
      END
