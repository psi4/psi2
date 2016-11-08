      SUBROUTINE HT1   (T1,T1O,T2O,CIN,DIN,FOCK,HO,HV,HIA,
     .                  ORBSYM,FLOV,COFF,CADD,DOFF,DADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,VAL1,VAL2,FVAL,tmp
      REAL*8  T1(NO,NV),T1O(NO,NV),T2O(*),CIN(*),DIN(*),FOCK(*),
     .        HO(NO,NO),HV(NV,NV),HIA(NO,NV)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        COFF(*),CADD(*),DOFF(*),DADD(*),map(*),spin(*),ssov(*)
C
      DO 191 U=1,NO
         IF(FZO(U).EQ.1)GO TO 191
         USYM=ORBSYM(U)
         uspn=spin(u)
         xu=map(u)
         FBE=FLOV(USYM+1,uspn,3)-NO
         LBE=FLOV(USYM+1,uspn,4)-NO
         FI =FLOV(USYM+1,uspn,1)
         LI =FLOV(USYM+1,uspn,2)
         DO 171 BE=FBE,LBE
            xbe=map(be+no)
            BEU=ipq(max0(xbe,xu))+min0(xbe,xu)
            FBEU=IPQ(BE+NO)+U
            VAL=0.0D0
            DO 135 I=FI,LI
               VAL=VAL-HO(U,I)*T1O(I,BE)
  135       CONTINUE
            DO 138 A=FBE,LBE
               VAL=VAL+HV(BE,A)*T1O(U,A)
  138       CONTINUE
            DO 148 I=1,NO
               IF(FZO(I).EQ.1)GO TO 148
               ispn=spin(i)
               xi=map(i)
               IU=ipq(MAX0(xI,xU))+MIN0(xI,xU)
               ISYM=ORBSYM(I)
               FA=FLOV(ISYM+1,ispn,3)-NO
               LA=FLOV(ISYM+1,ispn,4)-NO
               DO 145 A=FA,LA
                  xa=map(a+no)
                  Ia=ipq(MAX0(xI,xa))+MIN0(xI,xa)
                  BEUIA=DOFF(MAX0(BEU,IA))+DADD(MIN0(BEU,IA))
                  BEA=ipq(MAX0(xBE,xA))+MIN0(xBE,xA)
                  IUBEA=COFF(max0(iu,BEA))+CADD(min0(IU,bea))
                  tmp=din(beuia)
                  if(uspn.eq.ispn)tmp=tmp-cin(iubea)
                  VAL=VAL+tmp*T1O(I,A)
c
                  VAL1=+T1O(U,A)*T1O(I,BE)
                  VAL=VAL+HIA(I,A)*VAL1
c
                  FIA=IPQ(A+NO)+I
                  FVAL=FOCK(FIA)*T1O(U,A)*T1O(I,BE)
                  VAL=VAL-FVAL-FVAL
  145          CONTINUE
  148       CONTINUE
c
            DO 158 I=1,u-1
               IF(FZO(I).EQ.1)GO TO 158
               ispn=spin(i)
               ISYM=ORBSYM(I)
               FA=FLOV(ISYM+1,ispn,3)-NO
               LA=FLOV(ISYM+1,ispn,4)-NO
               la2=la
               if(la.ge.be)la2=be-1
               DO 155 A=FA,LA2
                  uibea=toff(u,i)+tadd(be,a)
                  VAL=VAL+HIA(I,A)*t2o(uibea)
  155          continue
               fa2=fa
               if(fa.le.be)fa2=be+1
               DO 156 A=FA2,LA
                  uibea=toff(u,i)+tadd(a,be)  
                  VAL=VAL-HIA(I,A)*t2o(uibea)
  156          continue
  158       continue
 
            DO 168 I=u+1,no
               IF(FZO(I).EQ.1)GO TO 168
               ispn=spin(i)
               ISYM=ORBSYM(I)
               FA=FLOV(ISYM+1,ispn,3)-NO
               LA=FLOV(ISYM+1,ispn,4)-NO
               la2=la
               if(la.ge.be)la2=be-1
               DO 165 A=FA,LA2
                  uibea=toff(i,u)+tadd(be,a)
                  VAL=VAL-HIA(I,A)*t2o(uibea)
  165          continue
               fa2=fa
               if(fa.le.be)fa2=be+1
               DO 166 A=FA2,LA
                  uibea=toff(i,u)+tadd(a,be)  
                  VAL=VAL+HIA(I,A)*t2o(uibea)
  166          continue
  168       continue
c
            VAL=VAL+FOCK(FBEU)
            T1(U,BE)=T1(U,BE)+VAL
  171    CONTINUE
  191 CONTINUE
c
      RETURN
      END
