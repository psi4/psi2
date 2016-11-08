      SUBROUTINE GOINT (T1O,T2O,EIN,FOCK,HO,GO,
     .                  ORBSYM,FLOV,EOFF,EADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,tmp
      REAL*8  T1O(NO,NV),T2O(*),EIN(*),FOCK(*),HO(NO,NO),GO(NO,NO)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        EOFF(*),EADD(*),map(*),spin(*),ssov(2,4)
C
      CALL ZERO(GO,NO*NO)
C
      DO 190 U=1,NO
         IF(FZO(U).EQ.1)GO TO 190
         USYM=ORBSYM(U)
         uspn=spin(u)
         xu=map(u)
         fi=ssov(uspn,1)
         li=ssov(uspn,2)
         DO 180 I=fi,li
            if(fzo(i).eq.1)go to 180
            xi=map(i)
            ISYM=ORBSYM(I)
            IF(USYM.NE.ISYM)GO TO 180
            VAL=0.0D0
            IU=ipq(MAX0(xI,xU))+MIN0(xI,xU)
            VAL=VAL+HO(U,I)
            FA=FLOV(ISYM+1,uspn,3)-NO
            LA=FLOV(ISYM+1,uspn,4)-NO
            DO 120 A=FA,LA
               FIA=IPQ(A+NO)+I
               VAL=VAL+FOCK(FIA)*T1O(U,A)
  120       CONTINUE
            DO 170 A=1,NV
               IF(FZV(A).EQ.1)GO TO 170
               aspn=spin(a+no)
               xa=map(a+no)
               ASYM=ORBSYM(A+NO)
               IA=ipq(max0(xi,xa))+min0(xi,xa)
               FJ=FLOV(ASYM+1,aspn,1)
               LJ=FLOV(ASYM+1,aspn,2)
               DO 160 J=FJ,LJ
                  xj=map(j)
                  jA=ipq(max0(xj,xa))+min0(xj,xa)
                  IUJA=EOFF(max0(iu,ja))+EADD(min0(iu,ja))
                  jU=ipq(MAX0(xj,xU))+MIN0(xj,xU)
                  jUiA=EOFF(max0(ju,ia))+EADD(min0(ju,ia))
                  tmp=ein(iuja)
                  if(uspn.eq.aspn)tmp=tmp-ein(juia)
                  VAL=VAL+tmp*T1O(J,A)
  160          CONTINUE
  170       CONTINUE
            GO(U,I)=VAL
  180    CONTINUE
  190 CONTINUE
C
      RETURN
      END
