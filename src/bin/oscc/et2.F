      SUBROUTINE ET2   (T1O,T2O,EIN,AUX1,XVAL,
     .                  ORBSYM,FLOV,EOFF,EADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL1,VAL2,TMP,TMV
      REAL*8  T1O(NO,NV),T2O(*),EIN(*),XVAL(*),AUX1(NV,NO)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        EOFF(*),EADD(*),map(*),spin(*),ssov(2,4)
C
      CALL ZERO(AUX1,NO*NV)
      CALL ZERO(XVAL,NV)
C
      DO 191 I=1,NO
         IF(FZO(I).EQ.1)GO TO 191
         ISYM=ORBSYM(I)
         ispn=spin(i)
         xi=map(i)
         DO 181 J=1,i-1
            IF(FZO(J).EQ.1)GO TO 181
            jspn=spin(j)
            xj=map(j)
            JSYM=ORBSYM(J)
            JISYM=IEOR(JSYM,ISYM)
            DO 171 A=1,nv 
               if(fzv(a).eq.1)go to 171
               ASYM=ORBSYM(A+NO)
               xa=map(a+no)
               aspn=spin(a+no)
               if(ispn.eq.jspn.and.ispn.ne.aspn)go to 171
               ia=ipq(max0(xi,xa))+min0(xi,xa)
               ja=ipq(max0(xj,xa))+min0(xj,xa)
               BESYM=IEOR(ASYM,JISYM)+1
               if(aspn.eq.ispn)xspn=jspn
               if(aspn.eq.jspn)xspn=ispn
               FBE=FLOV(BESYM,xspn,3)-NO
               LBE=FLOV(BESYM,xspn,4)-NO
               FU=FLOV(BESYM,xspn,1)
               LU=FLOV(BESYM,xspn,2)
               lbe2=lbe
               if(lbe.ge.a)lbe2=a-1
               DO 101 BE=FBE,LBE2
                  IJABE=TOFF(I,J)+TADD(A,BE)
                  XVAL(BE)=T2O(IJABE)+T1O(I,A)*T1O(J,BE)
                  if(ispn.eq.jspn)xval(be)=xval(be)-T1O(J,A)*T1O(I,BE)
  101          CONTINUE
               fbe2=fbe
               if(fbe.le.a)fbe2=a+1
               DO 102 BE=FBE2,LBE
                  IJbea=TOFF(I,J)+TADD(be,a)
                  XVAL(BE)=T2O(IJbea)+T1O(I,be)*T1O(J,a)
                  if(ispn.eq.jspn)xval(be)=xval(be)-T1O(i,A)*T1O(j,BE)
  102          CONTINUE
               DO 161 U=FU,LU
                  xu=map(u)
                  JU=ipq(MAX0(xJ,xU))+MIN0(xJ,xU)
                  IAJU=eoff(max0(ia,ju))+eadd(min0(ia,ju))
                  TMP = -EIN(IAJU)
                  if(ispn.eq.jspn)then
                  iU=ipq(MAX0(xi,xU))+MIN0(xi,xU)
                  JAIU=eoff(max0(ja,iu))+eadd(min0(ja,iu))
                  TMP = tmp+EIN(JAIU)
                  endif
                  DO 151 BE=FBE,LBE2
                     TMV = AUX1(BE,U)
                     TMV = TMV + TMP *XVAL(BE)
                     AUX1(BE,U) = TMV
  151             CONTINUE
  161          CONTINUE
               DO 162 U=FU,LU
                  xu=map(u)
                  iU=ipq(MAX0(xi,xU))+MIN0(xi,xU)
                  JAIU=eoff(max0(ja,iu))+eadd(min0(ja,iu))
                  TMP = -EIN(jAiU)
                  if(ispn.eq.jspn)then
                  JU=ipq(MAX0(xJ,xU))+MIN0(xJ,xU)
                  IAJU=eoff(max0(ia,ju))+eadd(min0(ia,ju))
                  TMP = tmp+EIN(iAjU)
                  endif
                  DO 152 BE=FBE2,LBE
                     TMV = AUX1(BE,U)
                     TMV = TMV + TMP *XVAL(BE)
                     AUX1(BE,U) = TMV
  152             CONTINUE
  162          CONTINUE
  171       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
      RETURN
      END
