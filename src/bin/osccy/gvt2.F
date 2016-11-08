      SUBROUTINE GVT2  (T2,T2O,GV,XVAL,AUX3,
     .                  ORBSYM,FLOV,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  TMP
      REAL*8  T2(*),T2O(*),GV(NV,NV),XVAL(*),AUX3(NV,NV)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        map(*),spin(*),ssov(2,4)
C
      DO 191 U=1,NO
         IF(FZO(U).EQ.1)GO TO 191
         USYM=ORBSYM(U)
         uspn=spin(u)
         DO 181 V=1,u-1
            IF(FZO(V).EQ.1)GO TO 181
            vspn=spin(v)
            VSYM=ORBSYM(V)
            VUSYM=IEOR(VSYM,USYM)
            fga=ssov(vspn,3)-no
            lga=ssov(vspn,4)-no
            DO 121 GA=fga,lga
               if(fzv(ga).eq.1)go to 121
               GASYM=ORBSYM(GA+NO)
               ASYM=IEOR(GASYM,VUSYM)+1
               FA=FLOV(ASYM,uspn,3)-NO
               LA=FLOV(ASYM,uspn,4)-NO
               fa2=fa
               if(fa.le.ga)fa2=ga+1
               DO 111 A=FA2,LA
                  UVAGA=TOFF(U,V)+TADD(A,GA)
                  AUX3(A,GA)=T2O(UVAGA)
  111          CONTINUE
               la2=la
               if(la.ge.ga)la2=ga-1
               DO 112 A=FA,LA2
                  UVAGA=TOFF(U,V)+TADD(GA,A)
                  AUX3(A,GA)=-T2O(UVAGA)
  112          CONTINUE
  121       CONTINUE
            fbe=ssov(uspn,3)-no
            lbe=ssov(uspn,4)-no
            DO 151 BE=fbe,lbe
               if(fzv(be).eq.1)go to 151
               BESYM=ORBSYM(BE+NO)
               FA=FLOV(BESYM+1,uspn,3)-NO
               LA=FLOV(BESYM+1,uspn,4)-NO
               DO 125 A=FA,LA
                  XVAL(A)=GV(BE,A)
  125          CONTINUE
               GASYM=IEOR(BESYM,VUSYM)+1
               FGA=FLOV(GASYM,vspn,3)-NO
               LGA=FLOV(GASYM,vspn,4)-NO
               if(lga.ge.be)lga=be-1
               DO 141 GA=FGA,LGA
                  TMP=0.0D0
                  DO 131 A=FA,LA
                     TMP=TMP+AUX3(A,GA)*XVAL(A)
  131             CONTINUE
                  UVBEGA=TOFF(U,V)+TADD(BE,GA)
                  T2(UVBEGA)=T2(UVBEGA)+TMP
  141          CONTINUE
  151       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
      DO 291 U=1,NO
         IF(FZO(U).EQ.1)GO TO 291
         USYM=ORBSYM(U)
         uspn=spin(u)
         DO 281 V=1,u-1
            IF(FZO(V).EQ.1)GO TO 281
            vspn=spin(v)
            VSYM=ORBSYM(V)
            VUSYM=IEOR(VSYM,USYM)
            fga=ssov(uspn,3)-no
            lga=ssov(uspn,4)-no
            DO 221 GA=fga,lga
               if(fzv(ga).eq.1)go to 221
               GASYM=ORBSYM(GA+NO)
               ASYM=IEOR(GASYM,VUSYM)+1
               FA=FLOV(ASYM,vspn,3)-NO
               LA=FLOV(ASYM,vspn,4)-NO
               fa2=fa
               if(fa.le.ga)fa2=ga+1
               DO 211 A=FA2,LA
                  UVAGA=TOFF(U,V)+TADD(A,GA)
                  AUX3(A,GA)=-T2O(UVAGA)
  211          CONTINUE
               la2=la
               if(la.ge.ga)la2=ga-1
               DO 212 A=FA,LA2
                  UVAGA=TOFF(U,V)+TADD(GA,A)
                  AUX3(A,GA)= T2O(UVAGA)
  212          CONTINUE
  221       CONTINUE
            fbe=ssov(vspn,3)-no
            lbe=ssov(vspn,4)-no
            DO 251 BE=fbe,lbe
               if(fzv(be).eq.1)go to 251
               BESYM=ORBSYM(BE+NO)
               FA=FLOV(BESYM+1,vspn,3)-NO
               LA=FLOV(BESYM+1,vspn,4)-NO
               DO 225 A=FA,LA
                  XVAL(A)=GV(BE,A)
  225          CONTINUE
               GASYM=IEOR(BESYM,VUSYM)+1
               FGA=FLOV(GASYM,uspn,3)-NO
               LGA=FLOV(GASYM,uspn,4)-NO
               if(fga.le.be)fga=be+1
               DO 241 GA=FGA,LGA
                  TMP=0.0D0
                  DO 231 A=FA,LA
                     TMP=TMP+AUX3(A,GA)*XVAL(A)
  231             CONTINUE
                  UVBEGA=TOFF(U,V)+TADD(ga,be)
                  T2(UVBEGA)=T2(UVBEGA)+TMP
  241          CONTINUE
  251       CONTINUE
  281    CONTINUE
  291 CONTINUE
      RETURN
      END
