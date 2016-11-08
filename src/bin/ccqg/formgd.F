      SUBROUTINE FORMGD(GDMAT,TAU,T2,T1,Z2,Z1,NO,NV,NO2,NV2,NTRO,NTRV,
     .                  IOFF,NTR,Q1,Q2,Q3,TZ2,VALO,VALV,JOUT,OPTION,
     .                  ORBSYM,FLOV,ZLX,ITR,ITV,NT3,NIRRED,NTAU,NOV,
     .                  T3OFF,T3ADD,W3,MOO,NVV,Y3,NDIMT2,TOFF,TADD,
     .                  DOFF,DADD,NSHOV,AUX,AUS,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      COMMON/TAPES/ITAP81,ITAP82,ITAP83,ITAP99,itap76,i76,itap74
      CHARACTER*4 OPTION
      REAL*8 VAL,VAL1,VAL2,VAL3,VAL4,VAL5,VAL6,VAL7,VAL8
      REAL*8 GDMAT(NSHOV),VALO(NO),
     .       T1(NO,NV),Z1(NO,NV),Q1(NO,NV),Q2(NO,NV),
     .       VALV(NV),Q3(NO,NV),TZ2(NO,NO),W3(NTAU),MOO(NO,NO),
     .       NVV(NV,NV),Y3(NTAU),AUX(NO,NO),AUS(NO,NO),
     .       T2(NDIMT2),Z2(NDIMT2),TAU(NDIMT2)
      INTEGER TOFF(NO,NO,2),TADD(NV,NV),
     .        IOFF(NTR),ITR(NOV),ITV(NV),
     .        ORBSYM(NO+NV),FLOV(NIRRED,4),ZLX(NV,NV),
     .        T3OFF(NO,NO,2,NIRRED),T3ADD(NV,NV,NIRRED),NT3(NIRRED),
     .        DOFF(NO*NV),DADD(NO*NV),FZO(NO),FZV(NV)
C
      CALL ZERO(GDMAT,NSHOV)
c
c     write(6,*)'in formgd: nvv'
c     call matout(nvv,nv,nv,nv,nv,6)
c     write(6,*)'in formgd: moo'
c     call matout(moo,no,no,no,no,6)
      if(option.eq.'CCT ')then 
c        write(6,*)'reading dden',nshov,'  long'
         call wreadw(itap76,gdmat,intowp(nshov),i76,i76)
c        write(6,'(8f10.6)')gdmat
      DO 991 A=1,NV
      IF(FZV(A).EQ.1)GO TO 991
      ASYM=ORBSYM(A+NO)
      DO 981 I=1,NO
      IF(FZO(I).EQ.1)GO TO 981
      ISYM=ORBSYM(I)
      IASYM=IEOR(ISYM,ASYM)
      AI=ITR(A)+I
      DO 971 B=1,nv
      IF(FZV(B).EQ.1)GO TO 971
      BSYM=ORBSYM(B+NO)
      JSYM=IEOR(IASYM,BSYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      DO 961 J=FJ,LJ
      BJ=ITR(B)+J
      AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
      IF(I.EQ.J)GDMAT(AIBJ)=GDMAT(AIBJ)-NVV(A,B)
  961 CONTINUE
  971 CONTINUE
  981 CONTINUE
  991 CONTINUE
      endif
C
C     FORM  MOO AND NVV
C
      IF(OPTION.EQ.'SDT1')THEN
      CALL ZERO(NVV,NV*NV)
      CALL ZERO(MOO,NO*NO)
      I81=1
      I82=1
      DO 4991 C=1,NV
      IF(FZV(C).EQ.1)GO TO 4991
      CSYM=ORBSYM(C+NO)
      DO 4981 K=1,NO
      IF(FZO(K).EQ.1)GO TO 4981
      KSYM=ORBSYM(K)
      KCSYM=IEOR(KSYM,CSYM)+1
      DIM=NT3(KCSYM)
      CALL WREADW(ITAP81,W3,INTOWP(DIM),I81,I81)
      CALL WREADW(ITAP82,Y3,INTOWP(DIM),I82,I82)
      DO 4971 B=1,NV
      IF(FZV(B).EQ.1)GO TO 4971
      BSYM=ORBSYM(B+NO)
      DO 4961 J=1,NO
      IF(FZO(J).EQ.1)GO TO 4961
      JSYM=ORBSYM(J)
      BJSYM=IEOR(BSYM,JSYM)
      KCBJS=IEOR(KCSYM-1,BJSYM)
      DO 4951 A=1,NV
      IF(FZV(A).EQ.1)GO TO 4951
      ASYM=ORBSYM(A+NO)
      ISYM=IEOR(KCBJS,ASYM)+1
      FI=FLOV(ISYM,1)
      LI=FLOV(ISYM,2)
      ZLAB=ZLX(A,B)

      FM=FLOV(ISYM,1)
      LM=FLOV(ISYM,2)
      DO 4941 I=FI,LI
      IJAB=T3OFF(I,J,ZLAB,KCSYM)+T3ADD(A,B,KCSYM)
      DO 4931 M=FM,LM
      MJAB=T3OFF(M,J,ZLAB,KCSYM)+T3ADD(A,B,KCSYM)
      MOO(I,M)=MOO(I,M)-Y3(IJAB)*W3(MJAB)
 4931 CONTINUE
 4941 CONTINUE
C
      FE=FLOV(ASYM+1,3)-NO
      LE=FLOV(ASYM+1,4)-NO
      DO 4942 I=FI,LI
      IJAB=T3OFF(I,J,ZLAB,KCSYM)+T3ADD(A,B,KCSYM)
      DO 4932 E=FE,LE
      ZLEB=ZLX(E,B)
      IJEB=T3OFF(I,J,ZLEB,KCSYM)+T3ADD(E,B,KCSYM)
      NVV(A,E)=NVV(A,E)-Y3(IJAB)*W3(IJEB)
 4932 CONTINUE
 4942 CONTINUE
c
 4951 CONTINUE
 4961 CONTINUE
 4971 CONTINUE
 4981 CONTINUE
 4991 CONTINUE
      DO 4993 I=1,NO
      DO 4994 J=1,NO
      MOO(I,J)=MOO(I,J)*0.5D0
 4994 CONTINUE
 4993 CONTINUE
      DO 4995 A=1,NV
      DO 4996 B=1,NV
      NVV(A,B)=NVV(A,B)*0.5D0
 4996 CONTINUE
 4995 CONTINUE
c     CALL MATOUT(MOO,NO,NO,NO,NO,6)
c     CALL MATOUT(NVV,NV,NV,NV,NV,6)
C
      I81=1
      DO 2991 A=1,NV
      IF(FZV(A).EQ.1)GO TO 2991
      ASYM=ORBSYM(A+NO)
      DO 2981 I=1,NO
      IF(FZO(I).EQ.1)GO TO 2981
      ISYM=ORBSYM(I)
      IASYM=IEOR(ISYM,ASYM)+1
      DIM=NT3(IASYM)
      CALL WREADW(ITAP81,W3,INTOWP(DIM),I81,I81)
      AI=ITR(A)+I
      DO 2971 B=1,nv
      IF(FZV(B).EQ.1)GO TO 2971
      BSYM=ORBSYM(B+NO)
      BI=ITR(B)+I
      JSYM=IEOR(IASYM-1,BSYM)+1
      FJ=FLOV(JSYM,1)
      LJ=FLOV(JSYM,2)
      DO 2961 J=FJ,LJ
      BJ=ITR(B)+J
      AJ=ITR(A)+J
      AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
      BIAJ=DOFF(MAX0(BI,AJ))+DADD(MIN0(BI,AJ))
      IF(I.EQ.J)GDMAT(AIBJ)=GDMAT(AIBJ)-NVV(A,B)
      DO 2951 C=1,NV
      IF(FZV(C).EQ.1)GO TO 2951
      ZLCB=ZLX(C,B)
      CSYM=ORBSYM(C+NO)
      FK=FLOV(CSYM+1,1)
      LK=FLOV(CSYM+1,2)
      DO 2941 K=FK,LK
      JKCB=T3OFF(J,K,ZLCB,IASYM)+T3ADD(C,B,IASYM)
      JKBC=T3OFF(K,J,ZLCB,IASYM)+T3ADD(B,C,IASYM)
      VAL= (W3(JKCB)-W3(JKBC))*Z1(K,C)
      GDMAT(AIBJ)=GDMAT(AIBJ)-VAL-VAL
      GDMAT(BIAJ)=GDMAT(BIAJ)+VAL
 2941 CONTINUE
 2951 CONTINUE
 2961 CONTINUE
 2971 CONTINUE
C
 2981 CONTINUE
 2991 CONTINUE
C
      ENDIF
C
C
      DO 10 A=1,NV
         IF(FZV(A).EQ.1)GO TO 10
         ASYM=ORBSYM(A+NO)
         DO 20 I=1,NO
            IF(FZO(I).EQ.1)GO TO 20
            ISYM=ORBSYM(I)
            AISYM=IEOR(ASYM,ISYM)
            AI=ITR(A)+I
            VAL5=0.0D+00
            VAL6=0.0D+00
            IF(ASYM.NE.ISYM)GO TO 199
            DO 19 K=1,NO
               IF(FZO(K).EQ.1)GO TO 19
               KSYM=ORBSYM(K)
               FC=FLOV(KSYM+1,3)-NO
               TC=FLOV(KSYM+1,4)-NO
               DO 18 C=FC,TC
                  ZLAC=ZLX(A,C)
                  IKAC=TOFF(I,K,ZLAC)+TADD(A,C)
                  KIAC=TOFF(K,I,ZLAC)+TADD(A,C)
                  VAL5=VAL5+T2(IKAC)*Z1(K,C)
                  VAL6=VAL6+TAU(KIAC)*Z1(K,C)
   18          CONTINUE
   19       CONTINUE
  199       CONTINUE
            DO 30 B=1,NV
               IF(FZV(B).EQ.1)GO TO 30
               BSYM=ORBSYM(B+NO)
               ZLAB=ZLX(A,B)
               BI=ITR(B)+I
               JSYM=IEOR(AISYM,BSYM)+1
               FJ=FLOV(JSYM,1)
               TJ=FLOV(JSYM,2)
               DO 40 J=FJ,TJ
                  IJAB=TOFF(I,J,ZLAB)+TADD(A,B)
                  JIAB=TOFF(J,I,ZLAB)+TADD(A,B)
                  BJ=ITR(B)+J
                  AJ=ITR(A)+J
                  AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  AJBI=DOFF(MAX0(AJ,BI))+DADD(MIN0(AJ,BI))
                  VAL1=Z2(IJAB)*0.5D+00
                  VAL2=TAU(IJAB)+TAU(IJAB)
                  VAL3=TAU(JIAB)
                  VAL4=T1(I,A)*Z1(J,B)
                  GDMAT(AIBJ)=GDMAT(AIBJ)+ VAL1+VAL2-VAL3+VAL4+VAL4
     .                        +T1(J,B)*(4.0D+00*VAL5-2.0D+00*VAL6)
                  GDMAT(AJBI)=GDMAT(AJBI)+ T1(J,B)*(VAL6-2.0D+00*VAL5)
   40          CONTINUE
   30       CONTINUE
C
C
      CALL ZERO(Q1,NO*NV)
      CALL ZERO(Q2,NO*NV)
      CALL ZERO(Q3,NO*NV)
            DO 330 L=1,NO
               IF(FZO(L).EQ.1)GO TO 330
               LSYM=ORBSYM(L)
               DSYM=IEOR(AISYM,LSYM)+1
               FD=FLOV(DSYM,3)-NO
               TD=FLOV(DSYM,4)-NO
               DO 340 D=FD,TD
                  DL=ITR(D)+L
                  Q1(L,D)=0.0D+00
                  Q2(L,D)=0.0D+00
                  Q3(L,D)=0.0D+00
                  DO 350 K=1,NO
                     IF(FZO(K).EQ.1)GO TO 350
                     KSYM=ORBSYM(K)
                     CSYM=IEOR(AISYM,KSYM)+1
                     FC=FLOV(CSYM,3)-NO
                     TC=FLOV(CSYM,4)-NO
                     DO 360 C=FC,TC
                        ZLAC=ZLX(A,C)
                        ZLDC=ZLX(D,C)
                        IKAC=TOFF(I,K,ZLAC)+TADD(A,C)
                        KIAC=TOFF(K,I,ZLAC)+TADD(A,C)
                        LKDC=TOFF(L,K,ZLDC)+TADD(D,C)
                        KLDC=TOFF(K,L,ZLDC)+TADD(D,C)
                        Q1(L,D)=Q1(L,D)+T2(IKAC)*Z2(LKDC)
                        Q2(L,D)=Q2(L,D)+T2(KIAC)*Z2(LKDC)
                        Q3(L,D)=Q3(L,D)+T2(KIAC)*Z2(KLDC)
  360                CONTINUE
  350             CONTINUE
                  AIDL=DOFF(MAX0(AI,DL))+DADD(MIN0(AI,DL))
                  GDMAT(AIDL)=GDMAT(AIDL)+Q1(L,D)+Q1(L,D)-Q2(L,D)
  340          CONTINUE
  330       CONTINUE
C
            CALL ZERO(TZ2,NO*NO)
            CALL ZERO(AUX,NO*NO)
            CALL ZERO(AUS,NO*NO)
            DO 359 L=1,NO
               IF(FZO(L).EQ.1)GO TO 359
               LSYM=ORBSYM(L)
               DSYM=IEOR(AISYM,LSYM)+1
               FD=FLOV(DSYM,3)-NO
               TD=FLOV(DSYM,4)-NO
               FJ=FLOV(DSYM,1)
               TJ=FLOV(DSYM,2)
               DO 357 D=FD,TD
                  ZLAD=ZLX(A,D)
                  ILAD=TOFF(I,L,ZLAD)+TADD(A,D)
                  DO 355 J=FJ,TJ
                     TZ2(L,J)=TZ2(L,J)+(Q1(L,D)+Q1(L,D)-Q2(L,D))*T1(J,D)
                     AUX(L,J)=AUX(L,J)+(Q1(L,D)+Q3(L,D))        *T1(J,D)
                     AUS(L,J)=AUS(L,J)+ Z2(ILAD)                *T1(J,D)
  355             CONTINUE
  357          CONTINUE
  359       CONTINUE
C
            DO 430 J=1,NO
               IF(FZO(J).EQ.1)GO TO 430
               JSYM=ORBSYM(J)
               AJ=ITR(A)+J
               BSYM=IEOR(AISYM,JSYM)+1
               FB=FLOV(BSYM,3)-NO
               TB=FLOV(BSYM,4)-NO
               FL=FLOV(BSYM,1)
               TL=FLOV(BSYM,2)
               DO 440 B=FB,TB
                  BI=ITR(B)+I
                  BJ=ITR(B)+J
                  AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  AJBI=DOFF(MAX0(AJ,BI))+DADD(MIN0(AJ,BI))
                  VAL1=0.0D+00
                  VAL2=0.0D+00
                  VAL3=0.0D+00
                  VAL4=0.0D+00
                  VAL5=0.0D0
                  VAL6=0.0D0
                  DO 365 L=FL,TL
                     VAL5=VAL5+(TZ2(L,J)+AUS(L,J))*T1(L,B)
                     VAL6=VAL6+AUX(L,J)*T1(L,B)
  365             CONTINUE
                  DO 450 L=1,NO
                     IF(FZO(L).EQ.1)GO TO 450
                     LSYM=ORBSYM(L)
                     DSYM=IEOR(AISYM,LSYM)+1
                     FD=FLOV(DSYM,3)-NO
                     TD=FLOV(DSYM,4)-NO
                     DO 460 D=FD,TD
                        ZLBD=ZLX(B,D)
                        JLBD=TOFF(J,L,ZLBD)+TADD(B,D)
                        LJBD=TOFF(L,J,ZLBD)+TADD(B,D)
                        VAL1=VAL1+(Q1(L,D)-Q2(L,D))*T2(JLBD)
                        VAL2=VAL2+(Q1(L,D)-Q2(L,D))*(T2(JLBD)-T2(LJBD))
                        VAL3=VAL3+Q1(L,D)*T2(JLBD)
                        VAL4=VAL4+Q3(L,D)*T2(LJBD)
  460                CONTINUE
  450             CONTINUE
                  GDMAT(AIBJ)=GDMAT(AIBJ)+VAL1+0.5D+00*(VAL2+VAL3)-VAL5
                  GDMAT(AJBI)=GDMAT(AJBI)-(VAL1-0.5D+00*VAL4-VAL6 )
  440          CONTINUE
  430       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      DO 810 J=1,NO
         IF(FZO(J).EQ.1)GO TO 810
         JSYM=ORBSYM(J)
         DO 820 K=1,NO
            IF(FZO(K).EQ.1)GO TO 820
            VALO(K)=0.0D+00
            DO 830 C=1,NV
               IF(FZV(C).EQ.1)GO TO 830
               VALO(K)=VALO(K)+T1(J,C)*Z1(K,C)
  830       CONTINUE
  820    CONTINUE
         DO 840 A=1,NV
            IF(FZV(A).EQ.1)GO TO 840
            ASYM=ORBSYM(A+NO)
            JASYM=IEOR(ASYM,JSYM)
            AJ=ITR(A)+J
            FK=FLOV(JSYM+1,1)
            TK=FLOV(JSYM+1,2)
            DO 850 I=1,NO
               IF(FZO(I).EQ.1)GO TO 850
               ISYM=ORBSYM(I)
               BSYM=IEOR(JASYM,ISYM)+1
               FB=FLOV(BSYM,3)-NO
               TB=FLOV(BSYM,4)-NO
               AI=ITR(A)+I
               DO 860 B=FB,TB
                  ZLAB=ZLX(A,B)
                  BI=ITR(B)+I
                  BJ=ITR(B)+J
                  AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  AJBI=DOFF(MAX0(AJ,BI))+DADD(MIN0(AJ,BI))
                  VAL1=0.0D+00
                  DO 870 K=FK,TK
                     IKAB=TOFF(I,K,ZLAB)+TADD(A,B)
                     VAL1=VAL1+VALO(K)*T2(IKAB)
  870             CONTINUE
                  GDMAT(AIBJ)=GDMAT(AIBJ)-(VAL1+VAL1)
                  GDMAT(AJBI)=GDMAT(AJBI)+VAL1
  860          CONTINUE
  850       CONTINUE
  840    CONTINUE
C
         CALL ZERO(VALO,NO)
         FL=FLOV(JSYM+1,1)
         TL=FLOV(JSYM+1,2)
         DO 809 L=FL,TL
            DO 819 K=1,NO
               IF(FZO(K).EQ.1)GO TO 819
               KSYM=ORBSYM(K)
               JKSYM=IEOR(JSYM,KSYM)
               DO 829 D=1,NV
                  IF(FZV(D).EQ.1)GO TO 829
                  DSYM=ORBSYM(D+NO)
                  CSYM=IEOR(JKSYM,DSYM)+1
                  FC=FLOV(CSYM,3)-NO
                  TC=FLOV(CSYM,4)-NO
                  DO 839 C=FC,TC
                     ZLDC=ZLX(D,C)
                     JKDC=TOFF(J,K,ZLDC)+TADD(D,C)
                     LKDC=TOFF(L,K,ZLDC)+TADD(D,C)
                     VALO(L)=VALO(L)+T2(JKDC)*Z2(LKDC)
  839             CONTINUE
  838             CONTINUE
  829          CONTINUE
  819       CONTINUE
  809    CONTINUE
         DO 817 B=1,NV
            IF(FZV(B).EQ.1)GO TO 817
            BSYM=ORBSYM(B+NO)
            JBSYM=IEOR(JSYM,BSYM)
            BJ=ITR(B)+J
            DO 827 I=1,NO
               IF(FZO(I).EQ.1)GO TO 827
               ISYM=ORBSYM(I)
               ASYM=IEOR(JBSYM,ISYM)+1
               FA=FLOV(ASYM,3)-NO
               TA=FLOV(ASYM,4)-NO
               DO 837 A=FA,TA
                  AI=ITR(A)+I
                  AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  ZLBA=ZLX(B,A)
                  VAL1=0.0D+00
                  DO 847 L=FL,TL
                     ILBA=TOFF(I,L,ZLBA)+TADD(B,A)
                     LIBA=TOFF(L,I,ZLBA)+TADD(B,A)
                     VAL1=VAL1+VALO(L)*(TAU(ILBA)-TAU(LIBA)-TAU(LIBA))
  847             CONTINUE
                  GDMAT(AIBJ)=GDMAT(AIBJ)+VAL1
  837          CONTINUE
  827       CONTINUE
  817    CONTINUE
  810 CONTINUE
C
      DO 910 B=1,NV
         IF(FZV(B).EQ.1)GO TO 910
         BSYM=ORBSYM(B+NO)
         DO 920 C=1,NV
            IF(FZV(C).EQ.1)GO TO 920
            VALV(C)=0.0D+00
            DO 930 K=1,NO
               IF(FZO(K).EQ.1)GO TO 930
               VALV(C)=VALV(C)+T1(K,B)*Z1(K,C)
  930       CONTINUE
  920    CONTINUE
         FC=FLOV(BSYM+1,3)-NO
         TC=FLOV(BSYM+1,4)-NO
         DO 940 A=1,NV
            IF(FZV(A).EQ.1)GO TO 940
            ASYM=ORBSYM(A+NO)
            BASYM=IEOR(BSYM,ASYM)
            DO 950 I=1,NO
               IF(FZO(I).EQ.1)GO TO 950
               ISYM=ORBSYM(I)
               JSYM=IEOR(BASYM,ISYM)+1
               FJ=FLOV(JSYM,1)
               TJ=FLOV(JSYM,2)
               AI=ITR(A)+I
               BI=ITR(B)+I
               DO 960 J=FJ,TJ
                  AJ=ITR(A)+J
                  BJ=ITR(B)+J
                  AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  AJBI=DOFF(MAX0(AJ,BI))+DADD(MIN0(AJ,BI))
                  VAL1=0.0D+00
                  DO 970 C=FC,TC
                     ZLAC=ZLX(A,C)
                     IJAC=TOFF(I,J,ZLAC)+TADD(A,C)
                     VAL1=VAL1+VALV(C)*T2(IJAC)
  970             CONTINUE
                  GDMAT(AIBJ)=GDMAT(AIBJ)-(VAL1+VAL1)
                  GDMAT(AJBI)=GDMAT(AJBI)+VAL1
  960          CONTINUE
  950       CONTINUE
  940    CONTINUE
C
         CALL ZERO(VALV,NV)
         FD=FLOV(BSYM+1,3)-NO
         TD=FLOV(BSYM+1,4)-NO
         DO 906 D=FD,TD
            DO 916 L=1,NO
               IF(FZO(L).EQ.1)GO TO 916
               LSYM=ORBSYM(L)
               BLSYM=IEOR(BSYM,LSYM)
               DO 926 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 926
                  KSYM=ORBSYM(K)
                  CSYM=IEOR(BLSYM,KSYM)+1
                  FC=FLOV(CSYM,3)-NO
                  TC=FLOV(CSYM,4)-NO
                  DO 936 C=FC,TC
                     ZLBC=ZLX(B,C)
                     ZLDC=ZLX(D,C)
                     LKBC=TOFF(L,K,ZLBC)+TADD(B,C)
                     LKDC=TOFF(L,K,ZLDC)+TADD(D,C)
                     VALV(D)=VALV(D)+T2(LKBC)*Z2(LKDC)
  936             CONTINUE
  926          CONTINUE
  916       CONTINUE
  906    CONTINUE
C
         DO 708 J=1,NO
            IF(FZO(J).EQ.1)GO TO 708
            JSYM=ORBSYM(J)
            BJSYM=IEOR(BSYM,JSYM)
            BJ=ITR(B)+J
            DO 718 A=1,NV
               IF(FZV(A).EQ.1)GO TO 718
               ASYM=ORBSYM(A+NO)
               ISYM=IEOR(BJSYM,ASYM)+1
               FI=FLOV(ISYM,1)
               TI=FLOV(ISYM,2)
               DO 728 I=FI,TI
                  AI=ITR(A)+I
                  AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  VAL1=0.0D+00
                  DO 738 D=FD,TD
                  ZLAD=ZLX(A,D)
                  IJAD=TOFF(I,J,ZLAD)+TADD(A,D)
                  JIAD=TOFF(J,I,ZLAD)+TADD(A,D)
                  VAL1=VAL1+VALV(D)*(TAU(JIAD)-TAU(IJAD)-TAU(IJAD))
  738             CONTINUE
                  GDMAT(AIBJ)=GDMAT(AIBJ)+VAL1
  728          CONTINUE
  718       CONTINUE
  708    CONTINUE
  910 CONTINUE
C
C
      DO 1001 I=1,NO
         IF(FZO(I).EQ.1)GO TO 1001
         ISYM=ORBSYM(I)
         DO 1011 J=1,NO
            IF(FZO(J).EQ.1)GO TO 1011
            JSYM=ORBSYM(J)
            IJSYM=IEOR(ISYM,JSYM)
            DO 1021 K=1,NO
               IF(FZO(K).EQ.1)GO TO 1021
               KSYM=ORBSYM(K)
               LSYM=IEOR(IJSYM,KSYM)+1
               FL=FLOV(LSYM,1)
               TL=FLOV(LSYM,2)
               DO 1031 L=FL,TL
                  TZ2(K,L)=0.0D+00
                  DO 1041 C=1,NV
                     IF(FZV(C).EQ.1)GO TO 1041
                     CSYM=ORBSYM(C+NO)
                     DSYM=IEOR(IJSYM,CSYM)+1
                     FD=FLOV(DSYM,3)-NO
                     TD=FLOV(DSYM,4)-NO
                     DO 1051 D=FD,TD
                        ZLCD=ZLX(C,D)
                        IJCD=TOFF(I,J,ZLCD)+TADD(C,D)
                        KLCD=TOFF(K,L,ZLCD)+TADD(C,D)
                        TZ2(K,L)=TZ2(K,L)+TAU(IJCD)*Z2(KLCD)
 1051                CONTINUE
 1041             CONTINUE
 1031          CONTINUE
 1021       CONTINUE
C
            DO 1022 A=1,NV
               IF(FZV(A).EQ.1)GO TO 1022
               ASYM=ORBSYM(A+NO)
               BSYM=IEOR(IJSYM,ASYM)+1
               FB=FLOV(BSYM,3)-NO
               TB=FLOV(BSYM,4)-NO
               AI=ITR(A)+I
               DO 1032 B=FB,TB
                  BJ=ITR(B)+J
                  AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  ZLAB=ZLX(A,B)
                  VAL1=0.0D+00
                  DO 1042 K=1,NO
                  IF(FZO(K).EQ.1)GO TO 1042
                  KSYM=ORBSYM(K)
                  LSYM=IEOR(IJSYM,KSYM)+1
                  FL=FLOV(LSYM,1)
                  TL=FLOV(LSYM,2)
                     DO 1054 L=FL,TL
                        KLAB=TOFF(K,L,ZLAB)+TADD(A,B)
                        VAL1=VAL1+TZ2(K,L)*TAU(KLAB)
 1054                CONTINUE
 1042             CONTINUE
                  GDMAT(AIBJ)=GDMAT(AIBJ)+0.5D+00*VAL1
 1032          CONTINUE
 1022       CONTINUE
 1011    CONTINUE
 1001 CONTINUE
C
C
      DO 2001 A=1,NV
         IF(FZV(A).EQ.1)GO TO 2001
         ASYM=ORBSYM(A+NO)
         FB=FLOV(ASYM+1,3)-NO
         TB=FLOV(ASYM+1,4)-NO
         DO 2011 B=FB,TB
            VAL1=0.0D+00
            VAL2=0.0D+00
            DO 2021 K=1,NO
               IF(FZO(K).EQ.1)GO TO 2021
               KSYM=ORBSYM(K)
               VAL2=VAL2+T1(K,A)*Z1(K,B)
               DO 2031 L=1,NO
                  IF(FZO(L).EQ.1)GO TO 2031
                  LSYM=ORBSYM(L)
                  KLSYM=IEOR(KSYM,LSYM)
                  CSYM=IEOR(ASYM,KLSYM)+1
                  FC=FLOV(CSYM,3)-NO
                  TC=FLOV(CSYM,4)-NO
                  DO 2041 C=FC,TC
                     ZLAC=ZLX(A,C)
                     ZLBC=ZLX(B,C)
                     KLAC=TOFF(K,L,ZLAC)+TADD(A,C)
                     KLBC=TOFF(K,L,ZLBC)+TADD(B,C)
                     VAL1=VAL1+T2(KLAC)*Z2(KLBC)
 2041             CONTINUE
 2031          CONTINUE
 2021       CONTINUE
            DO 2022 I=1,NO
               IF(FZO(I).EQ.1)GO TO 2022
               AI=ITR(A)+I
               DO 2032 J=I,I
                  IF(FZO(J).EQ.1)GO TO 2032
                  BJ=ITR(B)+J
                  AIBJ=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  GDMAT(AIBJ)=GDMAT(AIBJ)-(VAL1+VAL2)
 2032          CONTINUE
 2022       CONTINUE
 2011    CONTINUE
 2001 CONTINUE
C
      DO 2229 A=1,NV
         IF(FZV(A).EQ.1)GO TO 2229
         ASYM=ORBSYM(A+NO)
         DO 2228 I=1,NO
            IF(FZO(I).EQ.1)GO TO 2228
            ISYM=ORBSYM(I)
            AISYM=IEOR(ASYM,ISYM)
            AI=ITR(A)+I
            DO 2227 B=1,A
               IF(FZV(B).EQ.1)GO TO 2227
               BSYM=ORBSYM(B+NO)
               JSYM=IEOR(AISYM,BSYM)+1
               FJ=FLOV(JSYM,1)
               LJ=FLOV(JSYM,2)
               LJ2=LJ
               IF(B.EQ.A)LJ2=I
               DO 2226 J=FJ,LJ2
                  BJ=ITR(B)+J
                  AIBJ=DOFF(AI)+DADD(BJ)
                  IF(A.EQ.B.AND.I.EQ.J)GDMAT(AIBJ)=GDMAT(AIBJ)*2.0D0
 2226          CONTINUE
 2227       CONTINUE
 2228    CONTINUE
 2229 CONTINUE
C
    1 FORMAT (' GDMAT(',I3,',',I3,',',I3,',',I3,')=',F20.12)
      RETURN
      END
