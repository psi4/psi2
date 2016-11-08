      SUBROUTINE GVINT (T1O,T2O,FIN,FOCK,HV,GV,
     .                  ORBSYM,FLOV,FOFF,FADD,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  VAL,tmp
      REAL*8  T1O(NO,NV),T2O(*),FIN(*),FOCK(*),HV(NV,NV),GV(NV,NV)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        FOFF(*),FADD(*),map(*),spin(*),ssov(2,4)
C
      CALL ZERO(GV,NV*NV)
C
      DO 190 A=1,NV
         IF(FZV(A).EQ.1)GO TO 190
         ASYM=ORBSYM(A+NO)
         aspn=spin(a+no)
         xa=map(a+no)
         fbe=ssov(aspn,3)-no
         lbe=ssov(aspn,4)-no
         DO 180 BE=fbe,lbe
            if(fzv(be).eq.1)go to 180
            xbe=map(be+no)
            BESYM=ORBSYM(BE+NO)
            IF(BESYM.NE.ASYM)GO TO 180
            VAL=0.0D0
            BEA=ipq(MAX0(xA,xBE))+MIN0(xA,xBE)
            VAL=VAL+HV(BE,A)
            FI=FLOV(ASYM+1,aspn,1)
            LI=FLOV(ASYM+1,aspn,2)
            DO 110 I=FI,LI
               FIA=IPQ(A+NO)+I
               VAL=VAL-FOCK(FIA)*T1O(I,BE)
110         CONTINUE
            DO 170 I=1,NO
               IF(FZO(I).EQ.1)GO TO 170
               ISYM=ORBSYM(I)
               ispn=spin(i)
               xi=map(i)
               IA=ipq(max0(xi,xa))+min0(xi,xa)
               FB=FLOV(ISYM+1,ispn,3)-NO
               LB=FLOV(ISYM+1,ispn,4)-NO
               DO 160 B=FB,LB
                  xb=map(b+no)
                  IB=ipq(max0(xi,xb))+min0(xi,xb)
                  BEAIB=FOFF(max0(BEA,ib))+FADD(min0(bea,IB))
                  BEB=ipq(MAX0(xB,xBE))+MIN0(xB,xBE)
                  BEbIa=FOFF(max0(BEb,ia))+FADD(min0(beb,Ia))
                  tmp=fin(beaib)
                  if(aspn.eq.ispn)tmp=tmp-fin(bebia)
                  VAL=VAL+tmp*T1O(I,B)
  160          CONTINUE
  170       CONTINUE
            GV(BE,A)=VAL
  180    CONTINUE
  190 CONTINUE
C
      RETURN
      END
