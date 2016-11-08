      SUBROUTINE DIVIDE(T1,T2,FOCK,
     .                  ORBSYM,FLOV,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,spin)
      IMPLICIT INTEGER (A-Z)
      REAL*8  DELTA
      REAL*8  T1(NO,NV),T2(*),FOCK(*)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),spin(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV)
C
      DO 191 U=1,NO
         IF(FZO(U).EQ.1)GO TO 191
         USYM=ORBSYM(U)
         UU=IPQ(U)+U
         uspn=spin(u)
         DO 181 V=1,u-1
            IF(FZO(V).EQ.1)GO TO 181
            VSYM=ORBSYM(V)
            VUSYM=IEOR(VSYM,USYM)
            VV=IPQ(V)+V
            vspn=spin(v)
            DO 171 BE=1,NV
               IF(FZV(BE).EQ.1)GO TO 171
               bespn=spin(be+no)
               if(uspn.ne.bespn)go to 171
               BESYM=ORBSYM(BE+NO)
               GASYM=IEOR(BESYM,VUSYM)+1
               FGA=FLOV(GASYM,vspn,3)-NO
               LGA=FLOV(GASYM,vspn,4)-NO
               BEBE=IPQ(BE+NO)+BE+NO
               IF(LGA.Ge.BE)LGA=BE-1
               DO 161 GA=FGA,LGA
                  GAGA=IPQ(GA+NO)+GA+NO
                  DELTA=FOCK(UU)+FOCK(VV)-FOCK(BEBE)-FOCK(GAGA)
                  UVBEGA=TOFF(U,V)+TADD(BE,GA)
                  T2(UVBEGA)=T2(UVBEGA)/DELTA
  161          CONTINUE
  171       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
      DO 291 U=1,NO
         IF(FZO(U).EQ.1)GO TO 291
         USYM=ORBSYM(U)
         uspn=spin(u)
         FBE=FLOV(USYM+1,uspn,3)-NO
         LBE=FLOV(USYM+1,uspn,4)-NO
         UU=IPQ(U)+U
         DO 281 BE=FBE,LBE
            BEBE=IPQ(BE+NO)+BE+NO
            DELTA=FOCK(UU)-FOCK(BEBE)
            T1(U,BE)=T1(U,BE)/DELTA
  281    CONTINUE
  291 CONTINUE
C
      RETURN
      END
