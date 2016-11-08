      SUBROUTINE GOT2  (T2,T2O,GO,XVAL,AUX2,
     .                  ORBSYM,FLOV,
     .                  TOFF,TADD,IPQ,ITR,ITV,FZO,FZV,ZLX,
     .                  NO,NV,NIRRED,map,spin,ssov)
      IMPLICIT INTEGER (A-Z)
      REAL*8  TMP
      REAL*8  T2(*),T2O(*),GO(NO,NO),XVAL(*),AUX2(NO,NO)
      INTEGER ORBSYM(*),FLOV(8,2,4),FZO(*),FZV(*),
     .        IPQ(*),ITR(*),ITV(*),TOFF(NO,NO),TADD(NV,NV),ZLX(NV,NV),
     .        map(*),spin(*),ssov(2,4)
C
      call zero(aux2,no*no)
c
      DO 191 BE=1,NV
         IF(FZV(BE).EQ.1)GO TO 191
         BESYM=ORBSYM(BE+NO)
         bespn=spin(be+no) 
         DO 181 GA=1,be-1
            IF(FZV(GA).EQ.1)GO TO 181
            gaspn=spin(ga+no)
            GASYM=ORBSYM(GA+NO)
            BEGAS=IEOR(BESYM,GASYM)
            fv=ssov(gaspn,1)
            lv=ssov(gaspn,2)
            DO 121 V=fv,lv
               if(fzo(v).eq.1)go to 121
               VSYM=ORBSYM(V)
               ISYM=IEOR(VSYM,BEGAS)+1
               FI=FLOV(ISYM,bespn,1)
               LI=FLOV(ISYM,bespn,2)
               fi2=fi
               if(fi.le.v)fi2=v+1
               DO 111 I=FI2,LI
                  IVBEGA=TOFF(I,V)+TADD(BE,GA)
                  AUX2(I,V)=T2O(IVBEGA)
  111          CONTINUE
               li2=li
               if(li.ge.v)li2=v-1
               DO 112 I=FI,LI2
                  IVBEGA=TOFF(v,i)+TADD(BE,GA)
                  AUX2(I,V)=-T2O(IVBEGA)
  112          CONTINUE
  121       CONTINUE
            fu=ssov(bespn,1)
            lu=ssov(bespn,2)
            DO 151 U=fu,lu
               if(fzo(u).eq.1)go to 151
               USYM=ORBSYM(U)
               FI=FLOV(USYM+1,bespn,1)
               LI=FLOV(USYM+1,bespn,2)
               DO 125 I=FI,LI
                  XVAL(I)= GO(U,I)
  125          CONTINUE
               VSYM=IEOR(USYM,BEGAS)+1
               FV=FLOV(VSYM,gaspn,1)
               LV=FLOV(VSYM,gaspn,2)
               if(lv.ge.u)lv=u-1
               DO 141 V=FV,LV
                  TMP=0.0D0
                  DO 131 I=FI,LI
                     TMP=TMP+AUX2(I,V)*XVAL(I)
  131             CONTINUE
                    UVBEGA=TOFF(U,V)+TADD(BE,GA)
                    T2(UVBEGA)=T2(UVBEGA)-TMP
  141          CONTINUE
  151       CONTINUE
  181    CONTINUE
  191 CONTINUE
C
      call zero(aux2,no*no)
c
      DO 291 BE=1,NV
         IF(FZV(BE).EQ.1)GO TO 291
         BESYM=ORBSYM(BE+NO)
         bespn=spin(be+no) 
         DO 281 GA=1,be-1
            IF(FZV(GA).EQ.1)GO TO 281
            gaspn=spin(ga+no)
            GASYM=ORBSYM(GA+NO)
            BEGAS=IEOR(BESYM,GASYM)
            fv=ssov(bespn,1)
            lv=ssov(bespn,2)
            DO 221 V=fv,lv
               if(fzo(v).eq.1)go to 221
               VSYM=ORBSYM(V)
               ISYM=IEOR(VSYM,BEGAS)+1
               FI=FLOV(ISYM,gaspn,1)
               LI=FLOV(ISYM,gaspn,2)
               fi2=fi
               if(fi.le.v)fi2=v+1
               DO 211 I=FI2,LI
                  IVBEGA=TOFF(i,v)+TADD(BE,GA)
                  AUX2(I,V)=T2O(IVBEGA)
  211          CONTINUE
               li2=li
               if(li.ge.v)li2=v-1
               DO 212 I=FI,LI2
                  IVBEGA=TOFF(v,i)+TADD(BE,GA)
                  AUX2(I,V)=-T2O(IVBEGA)
  212          CONTINUE
  221       CONTINUE
            fu=ssov(gaspn,1)
            lu=ssov(gaspn,2)
            DO 251 U=fu,lu
               if(fzo(u).eq.1)go to 251
               USYM=ORBSYM(U)
               FI=FLOV(USYM+1,gaspn,1)
               LI=FLOV(USYM+1,gaspn,2)
               DO 225 I=FI,LI
                  XVAL(I)= GO(U,I)
  225          CONTINUE
               VSYM=IEOR(USYM,BEGAS)+1
               FV=FLOV(VSYM,bespn,1)
               LV=FLOV(VSYM,bespn,2)
               if(fv.le.u)fv=u+1
               DO 241 V=FV,LV
                  TMP=0.0D0
                  DO 231 I=FI,LI
                     TMP=TMP+AUX2(I,V)*XVAL(I)
  231             CONTINUE
                    UVBEGA=TOFF(v,u)+TADD(BE,GA)
                    T2(UVBEGA)=T2(UVBEGA)+TMP
  241          CONTINUE
  251       CONTINUE
  281    CONTINUE
  291 CONTINUE
c
      RETURN
      END
