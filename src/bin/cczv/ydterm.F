      subroutine ydterm(a,asym,i,isym,ai,aisym,ais,dim,tau,din,
     .                  fzo,fzv,orbsym,flov,nirred,no,nv,
     .                  itr,itv,ipq,doff,dadd,toff,tadd,zlx,uoff,
     .                  vadd,nt3,z1o)
      IMPLICIT INTEGER(A-Z)
      COMMON/NEWPAR/NAB(8),NAI(8),NIJ(8),MAXOO,MAXVV,MAXOV
      REAL*8 z1o(NO,NV),TAU(*),VAL,VAL1,VAL2,DIN(*)
      INTEGER ORBSYM(*),FLOV(NIRRED,4),TOFF(NO,NO,2,NIRRED),ZLX(NV,NV),
     .        TADD(NV,NV,NIRRED),FZO(NO),FZV(NV),UOFF(NO,NO,2),NT3(8),
     .        VADD(NV,NV),DOFF(*),DADD(*),ITR(*),ITV(*),IPQ(*)
c
      IASYM=aisym+1
      IA=ai
c
      DO 3180 J=1,NO
      if(fzo(j).eq.1)go to 3180
      JSYM=ORBSYM(J)
      JA=ITR(A)+J
      BSYM=IEOR(IASYM-1,JSYM)+1
      FB=FLOV(BSYM,3)-NO
      LB=FLOV(BSYM,4)-NO
      DO 3175 B=FB,LB
      IB=ITR(B)+I
      JB=ITR(B)+J
      IAJB=DOFF(MAX0(IA,JB))+DADD(MIN0(IA,JB))
      IBJA=DOFF(MAX0(IB,JA))+DADD(MIN0(IB,JA))
      DO 3170 K=1,NO
      if(fzo(k).eq.1)go to 3170
      KSYM=ORBSYM(K)
      FC=FLOV(KSYM+1,3)-NO
      LC=FLOV(KSYM+1,4)-NO
      LC2=LC
      IF(LC2.GT.B)LC2=B
      DO 3165 C=FC,LC2
      JKBC=TOFF(J,K,1,IASYM)+TADD(B,C,IASYM)
      TAU(JKBC)=TAU(JKBC)-(din(IAJB)+din(IAJB)-din(IBJA))*Z1O(K,C)
 3165 CONTINUE
      FC2=FC
      IF(FC2.LT.B)FC2=B
      DO 3166 C=FC2,LC
      KJCB=TOFF(K,J,1,IASYM)+TADD(C,B,IASYM)
      TAU(KJCB)=TAU(KJCB)-(din(IAJB)+din(IAJB)-din(IBJA))*Z1O(K,C)
 3166 CONTINUE
 3170 CONTINUE
 3175 CONTINUE
 
      CSYM=IEOR(IASYM-1,JSYM)+1
      FC=FLOV(CSYM,3)-NO
      LC=FLOV(CSYM,4)-NO
      DO 3775 C=FC,LC
      IC=ITR(C)+I
      JC=ITR(C)+J
      IAJC=DOFF(MAX0(IA,JC))+DADD(MIN0(IA,JC))
      ICJA=DOFF(MAX0(IC,JA))+DADD(MIN0(IC,JA))
      DO 3770 K=1,NO
      if(fzo(k).eq.1)go to 3770
      KSYM=ORBSYM(K)
      FB=FLOV(KSYM+1,3)-NO
      LB=FLOV(KSYM+1,4)-NO
      FB2=FB
      IF(FB.LT.C)FB2=C
      DO 3765 B=FB2,LB
      JKBC=TOFF(J,K,1,IASYM)+TADD(B,C,IASYM)
      TAU(JKBC)=TAU(JKBC)+(din(IAJC)+din(IAJC)-din(ICJA))*Z1O(K,B)
 3765 CONTINUE
      LB2=LB
      IF(LB.GT.C)LB2=C
      DO 3767 B=FB,LB2
      KJCB=TOFF(K,J,1,IASYM)+TADD(C,B,IASYM)
      TAU(KJCB)=TAU(KJCB)+(din(IAJC)+din(IAJC)-din(ICJA))*Z1O(K,B)
 3767 CONTINUE
 3770 CONTINUE
 3775 CONTINUE
 3180 CONTINUE
 
cS    IF(SORT.NE.'YES ')THEN
      DO 3380 J=1,NO
      if(fzo(j).eq.1)go to 3380
      JSYM=ORBSYM(J)
      JA=ITR(A)+J
      DO 3375 B=1,NV
      if(fzv(b).eq.1)go to 3375
      BSYM=ORBSYM(B+NO)
      IBSYM=IEOR(ISYM,BSYM)
      JBSYM=IEOR(JSYM,BSYM)
      IB=ITR(B)+I
      JB=ITR(B)+J
      IAJBS=IEOR(IASYM-1,JBSYM)
      DO 3370 K=1,NO
      if(fzo(k).eq.1)go to 3370
      KSYM=ORBSYM(K)
      KASYM=IEOR(KSYM,ASYM)
      KA=ITR(A)+K
      KB=ITR(B)+K
      CSYM=IEOR(IAJBS,KSYM)
      KCSYM=IEOR(KSYM,CSYM)
      FC=FLOV(CSYM+1,3)-NO
      LC=FLOV(CSYM+1,4)-NO
      LC2=LC
      IF(LC.GT.B)LC2=B
      DO 3365 C=FC,LC2
      IC=ITR(C)+I
      JC=ITR(C)+J
      KC=ITR(C)+K
      VAL=0.0D0
      JBKAS=IEOR(JBSYM,KASYM)
      IF(JBKAS.EQ.0)THEN
      JBKA=DOFF(MAX0(JB,KA))+DADD(MIN0(JB,KA))
      JAKB=DOFF(MAX0(JA,KB))+DADD(MIN0(JA,KB))
      VAL=VAL+(din(JBKA)+din(JBKA)-din(JAKB))*Z1O(I,C)
      ENDIF
      JBKCS=IEOR(JBSYM,KCSYM)
      IF(JBKCS.EQ.0)THEN
      JBKC=DOFF(MAX0(JB,KC))+DADD(MIN0(JB,KC))
      JCKB=DOFF(MAX0(JC,KB))+DADD(MIN0(JC,KB))
      VAL=VAL-(din(JBKC)+din(JBKC)-din(JCKB))*Z1O(I,A)
      ENDIF
      KCIBS=IEOR(KCSYM,IBSYM)
      IF(KCIBS.EQ.0)THEN
      KCIB=DOFF(MAX0(KC,IB))+DADD(MIN0(KC,IB))
      KBIC=DOFF(MAX0(KB,IC))+DADD(MIN0(KB,IC))
      VAL=VAL+(din(KCIB)+din(KCIB)-din(KBIC))*Z1O(J,A)
      ENDIF
      KCIAS=IEOR(KCSYM,IASYM-1)
      IF(KCIAS.EQ.0)THEN
      KCIA=DOFF(MAX0(KC,IA))+DADD(MIN0(KC,IA))
      KAIC=DOFF(MAX0(KA,IC))+DADD(MIN0(KA,IC))
      VAL=VAL-(din(KCIA)+din(KCIA)-din(KAIC))*Z1O(J,B)
      ENDIF
      JKBC=TOFF(J,K,1,IASYM)+TADD(B,C,IASYM)
      TAU(JKBC)=TAU(JKBC)+VAL
 3365 CONTINUE
      FC2=FC
      IF(FC.LT.B)FC2=B
      DO 3366 C=FC2,LC
      IC=ITR(C)+I
      JC=ITR(C)+J
      KC=ITR(C)+K
      VAL=0.0D0
      JBKAS=IEOR(JBSYM,KASYM)
      IF(JBKAS.EQ.0)THEN
      JBKA=DOFF(MAX0(JB,KA))+DADD(MIN0(JB,KA))
      JAKB=DOFF(MAX0(JA,KB))+DADD(MIN0(JA,KB))
      VAL=VAL+(din(JBKA)+din(JBKA)-din(JAKB))*Z1O(I,C)
      ENDIF
      JBKCS=IEOR(JBSYM,KCSYM)
      IF(JBKCS.EQ.0)THEN
      JBKC=DOFF(MAX0(JB,KC))+DADD(MIN0(JB,KC))
      JCKB=DOFF(MAX0(JC,KB))+DADD(MIN0(JC,KB))
      VAL=VAL-(din(JBKC)+din(JBKC)-din(JCKB))*Z1O(I,A)
      ENDIF
      KCIBS=IEOR(KCSYM,IBSYM)
      IF(KCIBS.EQ.0)THEN
      KCIB=DOFF(MAX0(KC,IB))+DADD(MIN0(KC,IB))
      KBIC=DOFF(MAX0(KB,IC))+DADD(MIN0(KB,IC))
      VAL=VAL+(din(KCIB)+din(KCIB)-din(KBIC))*Z1O(J,A)
      ENDIF
      KCIAS=IEOR(KCSYM,IASYM-1)
      IF(KCIAS.EQ.0)THEN
      KCIA=DOFF(MAX0(KC,IA))+DADD(MIN0(KC,IA))
      KAIC=DOFF(MAX0(KA,IC))+DADD(MIN0(KA,IC))
      VAL=VAL-(din(KCIA)+din(KCIA)-din(KAIC))*Z1O(J,B)
      ENDIF
      KJCB=TOFF(K,J,1,IASYM)+TADD(C,B,IASYM)
      TAU(KJCB)=TAU(KJCB)+VAL
 3366 CONTINUE
 3370 CONTINUE
 3375 CONTINUE
 3380 CONTINUE
cS    ENDIF
      return
      end
