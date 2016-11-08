      SUBROUTINE TDIFF(T2,T1,T2O,T1O,T1NORM,DELT,NO,NV,NDIMT2,spin,ssov,
     .                 ORBSYM,FLOV,NIRRED,TOFF,TADD,FZO,FZV,XE,XEO)
      IMPLICIT INTEGER(A-Z)
      REAL*8 XE,XEO,DELE,T1NORM,DELT
      REAL*8 T2(NDIMT2),T1(NO,NV),T2O(NDIMT2),T1O(NO,NV),t1sum(2)
      INTEGER ORBSYM(NO+NV),FLOV(8,2,4),TOFF(NO,NO),TADD(NV,NV),
     .        FZO(NO),FZV(NV),spin(*),ssov(2,4)
C
      T1NORM=0.0D+00
      t1sum(1)=0.0d0
      t1sum(2)=0.0d0
      DELT=0.0D+00
C
      iocc=0
      DO 10 U=1,NO
         IF(FZO(U).EQ.1)GO TO 10
         iocc=iocc+1
         USYM=ORBSYM(U)
         uspn=spin(u)
         fbe=ssov(uspn,3)-no
         lbe=ssov(uspn,4)-no
         DO 20 BE=fbe,lbe
            BESYM=ORBSYM(BE+NO)
            IF(USYM.EQ.BESYM)THEN
               DELT=DELT+(T1(U,BE)-T1O(U,BE))*(T1(U,BE)-T1O(U,BE))
               T1sum(uspn) =T1sum(uspn) +T1(U,BE)*T1(U,BE)
            ENDIF
            UBESYM=IEOR(USYM,BESYM)
            DO 30 V=1,U-1
               IF(FZO(V).EQ.1)GO TO 30
               vspn=spin(v)
               VSYM=ORBSYM(V)
               GASYM=IEOR(UBESYM,VSYM)+1
               FGA=FLOV(GASYM,vspn,3)-NO
               LGA=FLOV(GASYM,vspn,4)-NO
               if(lga.ge.be)lga=be-1
               DO 40 GA=FGA,LGA
                  UVBEGA=TOFF(U,V)+TADD(BE,GA)
                  DELT=DELT
     .     +      (T2O(UVBEGA)-T2(UVBEGA))*(T2O(UVBEGA)-T2(UVBEGA))
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      DELT=DSQRT(DELT)
      T1NORM=0.5d0*(DSQRT(T1sum(1))+dsqrt(t1sum(2)))/dsqrt(1.0d0*iocc)
C
      RETURN
      END
