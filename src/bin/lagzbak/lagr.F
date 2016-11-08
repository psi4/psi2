C
C----------------------------------------------------------------------
C
      SUBROUTINE LAGR (II,HI,GDM,QDM,IPQ,XLAG,NO,NV,NT,ORBSYM,FLOV,
     .                 NIRRED,IOFF,IADD,SCFTCC,CCTSCF,FOCK,
     .                 ENUC,ESCF,ETOT69,ITAP79,FZT,ZETO,ZETV,EX,ITR,
     .                 LFZO,LFZV,IFZO,IFZV,ITV,NRS,
     .                 SRTOUT,ITAP82,ITAP83,IMAT,GMAT,IAUX,GAUX,
     .                 I1OFF,I1ADD,X1INT,C2OFF,C2ADD,C2INT,
     .                 D2OFF,D2ADD,D2INT,AIN,CIN,DIN,EIN,
     .                 AOFF,AADD,COFF,CADD,DOFF,DADD,EOFF,EADD,
     .                 NSGOO,NSGOV,NSHOV,NSLOV,BUF,IBUF,LENBUF,
     .                 E3OFF,E3ADD,E3INT)
      IMPLICIT INTEGER (A-Z)
      REAL*8 II(*),HI(*),GDM(*),QDM(*),XLAG(NT,NT),FOCK(*),GAUX(*),
     .       ZETO(*),ZETV(*),EX(*),IMAT(NT,NT),GMAT(NT,NT),IAUX(*),
     .       X1INT(*),C2INT(*),D2INT(*),AIN(*),CIN(*),DIN(*),EIN(*),
     .       BUF(*),E3INT(*)
      REAL*8 VAL,ECC1,ECC2,ENUC,ESCF,ETOT69,XFAC,eref
      INTEGER IPQ(*),ITR(*),ORBSYM(NT),FLOV(NIRRED,2),IOFF(*),IADD(*),
     .        SCFTCC(*),CCTSCF(*),FZT(NT),IFZO(*),IFZV(*),ITV(*),NRS(*),
     .        I1OFF(*),I1ADD(*),C2OFF(*),C2ADD(*),EOFF(*),EADD(*),
     .        AOFF(*),AADD(*),COFF(*),CADD(*),DOFF(*),DADD(*),
     .        D2OFF(*),D2ADD(*),E3OFF(*),E3ADD(*),ibuf(*)
      CHARACTER*4 SRTOUT
C
      NT2=NT*(NT+1)/2
C
C     WRITE(6,*)' THE ONE-ELECTRON HAMILTONIAN IN SCF ORDERING'
C     CALL PRINT  (HI ,NT2,NT,6)
C
C >>> MODIFY THE ONE-E INTEGRALS  FOR ALL ORBITALS
C
      CALL RDINS (60,AIN,NSGOO,BUF,IBUF,LENBUF,NO,AOFF,AADD,ITR,ITR)
C
      DO 20 I=1,NO
         I1=CCTSCF(I)
         ISYM=ORBSYM(I1)
         DO 15 J=1,I
            J1=CCTSCF(J)
            JSYM=ORBSYM(J1)
            IF(ISYM.NE.JSYM)GO TO 15
            IJ=ITR(I)+J
            VAL=0.0D0
            DO 10 P=1,NO
               PP=CCTSCF(P)
               IF(FZT(PP).EQ.1)THEN
                  IP=ITR(MAX0(I,P))+MIN0(I,P)
                  JP=ITR(MAX0(J,P))+MIN0(J,P)
                  PP=ITR(P)+P
                  IJPP=AOFF(MAX0(IJ,PP))+AADD(MIN0(IJ,PP))
                  IPJP=AOFF(MAX0(IP,JP))+AADD(MIN0(IP,JP))
                  VAL=VAL+AIN(IJPP)+AIN(IJPP)-AIN(IPJP)
               ENDIF
   10       CONTINUE
            HIJ=IPQ(MAX0(I1,J1))+MIN0(I1,J1)
            HI(HIJ)=HI(HIJ)+VAL
   15    CONTINUE
   20 CONTINUE
C
      CALL RDINS (64,EIN,NSLOV,BUF,IBUF,LENBUF,NO,EOFF,EADD,ITR,ITR)
C
      DO 40 A=1,NV
         A1=CCTSCF(A+NO)
         ASYM=ORBSYM(A1)
         DO 35 J=1,NO
            J1=CCTSCF(J)
            JSYM=ORBSYM(J1)
            IF(ASYM.NE.JSYM)GO TO 35
            AJ=ITR(A)+J
            VAL=0.0D0
            DO 30 P=1,NO
               PP=CCTSCF(P)
               IF(FZT(PP).EQ.1)THEN
                  AP=ITR(A)+P
                  JP=ITR(MAX0(J,P))+MIN0(J,P)
                  PP=ITR(P)+P
                  AJPP=EOFF(AJ)+EADD(PP)
                  APJP=EOFF(AP)+EADD(JP)
                  VAL=VAL+EIN(AJPP)+EIN(AJPP)-EIN(APJP)
               ENDIF
   30       CONTINUE
            HAJ=IPQ(MAX0(A1,J1))+MIN0(A1,J1)
            HI(HAJ)=HI(HAJ)+VAL
   35    CONTINUE
   40 CONTINUE
C
      CALL RDINS (62,CIN,NSGOV,BUF,IBUF,LENBUF,NO,COFF,CADD,ITV,ITR)
      CALL RDINS (63,DIN,NSHOV,BUF,IBUF,LENBUF,NO,DOFF,DADD,ITR,ITR)
C
      DO 60 A=1,NV
         A1=CCTSCF(A+NO)
         ASYM=ORBSYM(A1)
         DO 55 B=1,A
            B1=CCTSCF(B+NO)
            BSYM=ORBSYM(B1)
            IF(ASYM.NE.BSYM)GO TO 55
            AB=ITV(A)+B
            VAL=0.0D0
            DO 50 P=1,NO
               PP=CCTSCF(P)
               IF(FZT(PP).EQ.1)THEN
                  AP=ITR(A)+P
                  BP=ITR(B)+P
                  PP=ITR(P)+P
                  ABPP=COFF(AB)+CADD(PP)
                  APBP=DOFF(AP)+DADD(BP)
                  VAL=VAL+CIN(ABPP)+CIN(ABPP)-DIN(APBP)
               ENDIF
   50       CONTINUE
            HAB=IPQ(MAX0(A1,B1))+MIN0(A1,B1)
            HI(HAB)=HI(HAB)+VAL
   55    CONTINUE
   60 CONTINUE
C
C     WRITE(6,*)'  '
C     WRITE(6,*)' THE MODIFIED ONE-HAMILTONIAN IN SCF ORDERING'
C     CALL PRINT  (HI ,NT2,NT,6)
C
      CALL ZERO(XLAG,NT*NT)
C
      ECC1=0.0D0
      ECC2=0.0D0
C
C >>> CONSTRUCT LAGRANGIAN FOR ACTIVE-ALL ORBITALS
C
C >>>  2-E  ACTIVE-ALL  BLOCK
C
      CALL RFILE(ITAP82)
      CALL RFILE(ITAP83)
      IPOINT=1
      GPOINT=1
      DO 190 TOTSYM=1,NIRRED
         DO 180 S=1,NT
            IF(FZT(S).EQ.1)GO TO 180
            SSYM=ORBSYM(S)
            DO 170 T=1,S
               IF(FZT(T).EQ.1)GO TO 170
               TSYM=ORBSYM(T)
               TSSYM=IEOR(TSYM,SSYM)
               IF(TSSYM.NE.TOTSYM-1)GO TO 170
               ST=IPQ(MAX0(S,T))+MIN0(S,T)
C
               IF(SRTOUT.EQ.'  NO')THEN
C
               DO 160 P=1,NT
                  IF(FZT(P).EQ.1)GO TO 160
                  PSYM=ORBSYM(P)
                  RSYM=IEOR(PSYM,TSSYM)+1
                  FR=FLOV(RSYM,1)
                  LR=FLOV(RSYM,2)
                  DO 150 Q=1,NT
                     QSYM=ORBSYM(Q)
                     IF(PSYM.NE.QSYM)GO TO 150
                     IF(FZT(Q).NE.1)THEN
                        VAL=0.0D0
                        DO 130 R=FR,LR
                           PR=IPQ(MAX0(P,R))+MIN0(P,R)
                           QR=IPQ(MAX0(Q,R))+MIN0(Q,R)
                           PRST=IOFF(MAX0(PR,ST))+IADD(MIN0(PR,ST))
                           QRST=IOFF(MAX0(QR,ST))+IADD(MIN0(QR,ST))
                           VAL=VAL+GDM(PRST)*II(QRST)
  130                   CONTINUE
                     ELSE
                        VAL=0.0D0
                        DO 140 R=FR,LR
                           PR=IPQ(MAX0(P,R))+MIN0(P,R)
                           QR=IPQ(MAX0(Q,R))+MIN0(Q,R)
                           PRST=IOFF(MAX0(PR,ST))+IADD(MIN0(PR,ST))
                           QRST=I1OFF(QR)+I1ADD(ST)
                           VAL=VAL+GDM(PRST)*X1INT(QRST)
  140                   CONTINUE
                     ENDIF
                     XFAC=2.0D0
                     IF(S.EQ.T)XFAC=1.0D0
                     XLAG(P,Q)=XLAG(P,Q)+(VAL+VAL)*XFAC
                     IF(P.EQ.Q)ECC2=ECC2+VAL*XFAC
  150             CONTINUE
  160          CONTINUE
C
               ELSE
C
               LRS=NRS(TOTSYM)
               CALL WREADW(ITAP82,IAUX,INTOWP(LRS),IPOINT,IPOINT)
               CALL WREADW(ITAP83,GAUX,INTOWP(LRS),GPOINT,GPOINT)
               ICNT=0
C
               DO 161 P=1,NT
                  IF(FZT(P).EQ.1)GO TO 161
                  PSYM=ORBSYM(P)
                  RSYM=IEOR(PSYM,TSSYM)+1
                  FR=FLOV(RSYM,1)
                  LR=FLOV(RSYM,2)
                  IF(LR.GT.P)LR=P
                  DO 131 R=FR,LR
                     ICNT=ICNT+1
                     IMAT(R,P)=IAUX(ICNT)
                     IMAT(P,R)=IAUX(ICNT)
                     GMAT(R,P)=GAUX(ICNT)
                     GMAT(P,R)=GAUX(ICNT)
  131             CONTINUE
  161          CONTINUE
               DO 163 P=1,NT
                  IF(FZT(P).EQ.1)GO TO 163
                  PSYM=ORBSYM(P)
                  RSYM=IEOR(PSYM,TSSYM)+1
                  FR=FLOV(RSYM,1)
                  LR=FLOV(RSYM,2)
                  DO 153 Q=1,NT
                     QSYM=ORBSYM(Q)
                     IF(PSYM.NE.QSYM)GO TO 153
                     IF(FZT(Q).NE.1)THEN
                        VAL=0.0D0
                        DO 133 R=FR,LR
                           VAL=VAL+GMAT(R,P)*IMAT(R,Q)
  133                   CONTINUE
                     ELSE
                        VAL=0.0D0
                        DO 135 R=FR,LR
                           QR=IPQ(MAX0(Q,R))+MIN0(Q,R)
                           QRST=I1OFF(QR)+I1ADD(ST)
                           VAL=VAL+GMAT(R,P)*X1INT(QRST)
  135                   CONTINUE
                     ENDIF
                     XFAC=2.0D0
                     IF(S.EQ.T)XFAC=1.0D0
                     XLAG(P,Q)=XLAG(P,Q)+(VAL+VAL)*XFAC
                     IF(P.EQ.Q)ECC2=ECC2+VAL*XFAC
  153             CONTINUE
  163          CONTINUE
C
               ENDIF
C
  170       CONTINUE
  180    CONTINUE
  190 CONTINUE
c
      call rclose(itap82,4)
c  itap83 is needed in back2
      call rclose(itap83,3)
C
C >>> 1-E  ACTIVE-ALL  BLOCK
C
      DO 260 P=1,NT
         IF(FZT(P).EQ.1)GO TO 260
         PSYM=ORBSYM(P)
         FR=FLOV(PSYM+1,1)
         LR=FLOV(PSYM+1,2)
         DO 250 Q=1,NT
COK         IF(FZT(Q).EQ.1)GO TO 250
            QSYM=ORBSYM(Q)
            IF(QSYM.NE.PSYM)GO TO 250
            VAL=0.0D0
            DO 230 R=FR,LR
               PR=IPQ(MAX0(P,R))+MIN0(P,R)
               QR=IPQ(MAX0(Q,R))+MIN0(Q,R)
               VAL=VAL+QDM(PR)*HI(QR)
  230       CONTINUE
            XLAG(P,Q)=XLAG(P,Q)+VAL
            IF(P.EQ.Q)ECC1=ECC1+VAL
  250    CONTINUE
  260 CONTINUE
C
C >>> ADD FZO-ALL CONTRIBUTION TO LAGRANGIAN
C
      DO 265 P=1,NT
         PP=SCFTCC(P)
         IF(FZT(P).NE.1.OR.PP.GT.NO)GO TO 265
         PSYM=ORBSYM(P)
         DO 255 Q=1,NT
COK         IF(FZT(Q).NE.1)GO TO 255
            QSYM=ORBSYM(Q)
            IF(PSYM.NE.QSYM)GO TO 255
            PQ=IPQ(MAX0(P,Q))+MIN0(P,Q)
            VAL=0.0D0
            DO 235 K=1,NT
               IF(FZT(K).EQ.1)GO TO 235
               KSYM=ORBSYM(K)
               FL=FLOV(KSYM+1,1)
               LL=FLOV(KSYM+1,2)
               DO 225 L=FL,LL
                  KL=IPQ(MAX0(K,L))+MIN0(K,L)
                  PK=IPQ(MAX0(P,K))+MIN0(P,K)
                  QL=IPQ(MAX0(Q,L))+MIN0(Q,L)
               IF(FZT(Q).EQ.0)THEN
                  PQKL=I1OFF(PQ)+I1ADD(KL)
                  PKQL=I1OFF(PK)+I1ADD(QL)
                  VAL=VAL+QDM(KL)*(X1INT(PQKL)+X1INT(PQKL)-X1INT(PKQL))
crsg           ELSE
crsg              PQKL=C2OFF(PQ)+C2ADD(KL)
crsg              PKQL=D2OFF(MAX0(PK,QL))+D2ADD(MIN0(PK,QL))
crsg              VAL=VAL+QDM(KL)*(C2INT(PQKL)+C2INT(PQKL)-D2INT(PKQL))
crsg           ENDIF
crsg non-canonical addressing scheme of fz-act-fz-act integrals
crsg in symcnt requires this additional check
c this line works for FZC and FZV in N2 STO (2fzc,1fzv)
               else
                  pqkl=c2off(pq)+c2add(kl)
                  if(p.lt.q.and.k.gt.l.and.k.gt.q.
     1                      or.q.lt.p.and.l.gt.p.and.l.gt.k) then
                      pkql=d2off(min0(pk,ql))+d2add(max0(pk,ql))
                  else
                      pkql=d2off(max0(pk,ql))+d2add(min0(pk,ql))
                  endif
                  val=val+qdm(kl)*(c2int(pqkl)+c2int(pqkl)-d2int(pkql))
               endif
  225          CONTINUE
  235       CONTINUE
            XLAG(P,Q)=XLAG(P,Q)+VAL
  255    CONTINUE
  265 CONTINUE
C
C     WRITE(6,*)'  '
C     WRITE(6,*)' (      CORRELATED ) LAGRANGIAN IN SCF ORDERING'
C     CALL MATOUT (XLAG,NT,NT,NT,NT,6)
C
C >>> READ one e ints from 79 to calculate reference energy
C
      NT2=NT*(NT+1)/2
      START=1
      CALL WREADW(ITAP79,FOCK,INTOWP(NT2),START,END)
      eref=0.0d0
      do 290 i=1,no
         fii=ipq(i)+i
         eref=eref+fock(fii)
  290 continue
C
C >>> READ FOCK MATRIX FROM FILE79
C
      START=1+INTOWP(NT2)
      CALL WREADW(ITAP79,FOCK,INTOWP(NT2),START,END)
      do 291 i=1,no
         fii=ipq(i)+i
         eref=eref+fock(fii)
  291 continue
      eref=eref+enuc 
C
C     WRITE(6,*)' '
C     WRITE(6,*)' FOCK MATRIX FROM TAPE79 IN CC ORDERING'
C     CALL PRINT(FOCK,NT2,NT,6)
C
C >>> ADD HF CONTRIBUTION TO LAGRANGIAN
C
      DO 410 Q=1,NT
         QSYM=ORBSYM(Q)
         QQ=SCFTCC(Q)
         DO 405 J=1,NO
             JJ=CCTSCF(J)
             JSYM=ORBSYM(JJ)
             IF(QSYM.NE.JSYM)GO TO 405
             FJQ=IPQ(MAX0(J,QQ))+MIN0(J,QQ)
             XLAG(JJ,Q)=XLAG(JJ,Q)+FOCK(FJQ)+FOCK(FJQ)
  405    CONTINUE
  410 CONTINUE
C
C >>> SCALE BY FACTOR 2
C
      DO 560 P=1,NT
         PSYM=ORBSYM(P)
         DO 550 Q=1,NT
            QSYM=ORBSYM(Q)
            IF(PSYM.NE.QSYM)GO TO 550
            XLAG(P,Q)=XLAG(P,Q)*2.0D0
  550    CONTINUE
  560 CONTINUE
C
C     WRITE(6,*)'  '
C     WRITE(6,*)' ( HF + CORRELATED ) LAGRANGIAN IN SCF ORDERING'
C     CALL MATOUT (XLAG,NT,NT,NT,NT,6)
C
C >>> CONSTRUCT Z VECTOR FOR FZO-DOC AND FZV-VIR ROTATIONS
C
      CALL ZERO(ZETO,LFZO*NO)
      CALL ZERO(ZETV,LFZV*NV)
      CALL ZERO(EX,NO*NV)
C
      DO 660 P=1,NO
         PP=CCTSCF(P)
         IF(FZT(PP).NE.1)GO TO 660
         FPP=IPQ(P)+P
         PSYM=ORBSYM(PP)
         DO 650 Q=1,NO
            QQ=CCTSCF(Q)
            IF(FZT(QQ).EQ.1)GO TO 650
            QSYM=ORBSYM(QQ)
            IF(PSYM.NE.QSYM)GO TO 650
Cbug        PQ=ITR(MAX0(PP,QQ))+MIN0(PP,QQ)
            PQ=ITR(MAX0(P ,Q ))+MIN0(P ,Q )
            PQZ=ITR(P)+Q
            FQQ=IPQ(Q)+Q
            ZPQ=IFZO(PQZ)
            ZETO(ZPQ)=(XLAG(PP,QQ)-XLAG(QQ,PP))/(FOCK(FPP)-FOCK(FQQ))
            DO 640 B=1,NV
               BB=CCTSCF(B+NO)
               BSYM=ORBSYM(BB)
               DO 630 J=1,NO
                  JJ=CCTSCF(J)
                  JSYM=ORBSYM(JJ)
                  IF(JSYM.NE.BSYM)GO TO 630
                  XBJ=ITR(B)+J
                  BJ=ITR(B)+J
                  BP=ITR(B)+P
                  BQ=ITR(B)+Q
                  JP=ITR(MAX0(J,P))+MIN0(J,P)
                  JQ=ITR(MAX0(J,Q))+MIN0(J,Q)
                  PQBJ=EOFF(BJ)+EADD(PQ)
                  JQBP=EOFF(BP)+EADD(JQ)
                  JPBQ=EOFF(BQ)+EADD(JP)
                  VAL=4.0D0*EIN(PQBJ)-EIN(JQBP)-EIN(JPBQ)
                  EX(XBJ)=EX(XBJ)+ZETO(ZPQ)*VAL
  630          CONTINUE
  640       CONTINUE
  650    CONTINUE
  660 CONTINUE
C
      DO 665 P=1,NV
         PP=CCTSCF(P+NO)
         IF(FZT(PP).NE.1)GO TO 665
         FPP=IPQ(P+NO)+P+NO
         PSYM=ORBSYM(PP)
         DO 655 Q=1,NV
            QQ=CCTSCF(Q+NO)
            IF(FZT(QQ).EQ.1)GO TO 655
            QSYM=ORBSYM(QQ)
            IF(PSYM.NE.QSYM)GO TO 655
            PQ=IPQ(MAX0(PP,QQ))+MIN0(PP,QQ)
            PQZ=ITV(P)+Q
            FQQ=IPQ(Q+NO)+Q+NO
            ZPQ=IFZV(PQZ)
            ZETV(ZPQ)=(XLAG(PP,QQ)-XLAG(QQ,PP))/(FOCK(FPP)-FOCK(FQQ))
            DO 645 B=1,NV
               BB=CCTSCF(B+NO)
               BSYM=ORBSYM(BB)
               DO 635 J=1,NO
                  JJ=CCTSCF(J)
                  JSYM=ORBSYM(JJ)
                  IF(JSYM.NE.BSYM)GO TO 635
                  XBJ=ITR(B)+J
                     BJ=IPQ(MAX0(BB,JJ))+MIN0(BB,JJ)
                     BP=IPQ(MAX0(BB,PP))+MIN0(BB,PP)
                     BQ=IPQ(MAX0(BB,QQ))+MIN0(BB,QQ)
                     JP=IPQ(MAX0(JJ,PP))+MIN0(JJ,PP)
                     JQ=IPQ(MAX0(JJ,QQ))+MIN0(JJ,QQ)
C
                  IF(FZT(BB).EQ.0.AND.FZT(JJ).EQ.0)THEN
                     PQBJ=I1OFF(PQ)+I1ADD(BJ)
                     JQBP=I1OFF(BP)+I1ADD(JQ)
                     JPBQ=I1OFF(JP)+I1ADD(BQ)
                     VAL=4.0D0*X1INT(PQBJ)-X1INT(JQBP)-X1INT(JPBQ)
                  ENDIF
C
crsg non-canonical addressing scheme of fz-act-fz-act integrals
crsg in symcnt requires additional checks below
c                 IF(FZT(BB).EQ.1.AND.FZT(JJ).EQ.0)THEN
c                    PQBJ=D2OFF(MAX0(BJ,PQ))+D2ADD(MIN0(BJ,PQ))
c                    JQBP=C2OFF(BP)+C2ADD(JQ)
c                    JPBQ=D2OFF(MAX0(JP,BQ))+D2ADD(MIN0(JP,BQ))
c                    VAL=4.0D0*D2INT(PQBJ)-C2INT(JQBP)-D2INT(JPBQ)
c                 ENDIF
C
c                 IF(FZT(BB).EQ.0.AND.FZT(JJ).EQ.1)THEN
c                    PQBJ=D2OFF(MAX0(BJ,PQ))+D2ADD(MIN0(BJ,PQ))
c                    JQBP=D2OFF(MAX0(JQ,BP))+D2ADD(MIN0(JQ,BP))
c                    JPBQ=C2OFF(JP)+C2ADD(BQ)
c                    VAL=4.0D0*D2INT(PQBJ)-D2INT(JQBP)-C2INT(JPBQ)
c                 ENDIF
C
                  if(fzt(bb).eq.1.and.fzt(jj).eq.0)then
                     if(pp.lt.bb.and.qq.gt.jj.and.qq.gt.bb.
     1                 or.bb.lt.pp.and.jj.gt.pp.and.jj.gt.qq) then
                        pqbj=d2off(min0(bj,pq))+d2add(max0(bj,pq))
                     else
                        pqbj=d2off(max0(bj,pq))+d2add(min0(bj,pq))
                     endif
                     jqbp=c2off(bp)+c2add(jq)
                     if(pp.lt.bb.and.jj.gt.bb.and.jj.gt.qq.
     1                 or.bb.lt.pp.and.qq.gt.jj.and.qq.gt.pp) then
                        jpbq=d2off(min0(jp,bq))+d2add(max0(jp,bq))
                     else
                        jpbq=d2off(max0(jp,bq))+d2add(min0(jp,bq))
                     endif
                     val=4.0d0*d2int(pqbj)-c2int(jqbp)-d2int(jpbq)
                  endif
C
                  if(fzt(bb).eq.0.and.fzt(jj).eq.1)then
                     if(pp.lt.jj.and.qq.gt.bb.and.qq.gt.jj.
     1                 or.jj.lt.pp.and.bb.gt.pp.and.bb.gt.qq) then
                        pqbj=d2off(min0(bj,pq))+d2add(max0(bj,pq))
                     else
                        pqbj=d2off(max0(bj,pq))+d2add(min0(bj,pq))
                     endif
                     if(pp.lt.jj.and.bb.gt.jj.and.bb.gt.qq.
     1                 or.jj.lt.pp.and.qq.gt.bb.and.qq.gt.pp) then
                        jqbp=d2off(min0(jq,bp))+d2add(max0(jq,bp))
                     else
                        jqbp=d2off(max0(jq,bp))+d2add(min0(jq,bp))
                     endif
                     jpbq=c2off(jp)+c2add(bq)
                     val=4.0d0*d2int(pqbj)-d2int(jqbp)-c2int(jpbq)
                  endif
crsg end changes
C
                  IF(FZT(BB).EQ.1.AND.FZT(JJ).EQ.1)THEN
                     PQBJ=E3OFF(PQ)+E3ADD(BJ)
                     JQBP=E3OFF(JQ)+E3ADD(BP)
                     JPBQ=E3OFF(BQ)+E3ADD(JP)
                     VAL=4.0D0*E3INT(PQBJ)-E3INT(JQBP)-E3INT(JPBQ)
                  ENDIF
C
                  EX(XBJ)=EX(XBJ)+ZETV(ZPQ)*VAL
  635          CONTINUE
  645       CONTINUE
  655    CONTINUE
  665 CONTINUE
C
C
      WRITE(6,*)'  '
      WRITE(6,*)' NUCLEAR REP. ENERGY = ', ENUC
      WRITE(6,*)' SCF          ENERGY = ', escf 
      WRITE(6,*)' REFERENCE    ENERGY = ', eref 
      WRITE(6,*)' CORRELATION  ENERGY = ', ECC1+ECC2
      WRITE(6,*)' TOTAL ENERGY LAGZBK = ', ECC1+ECC2+eref
      WRITE(6,*)' TOTAL ENERGY TAPE69 = ', ETOT69
C
C     WRITE(6,*)'  '
C     WRITE(6,*)' EX VECTOR'
C     CALL MATOUT (EX,NO,NV,NO,NV,6)
C
C     WRITE(6,*)'  '
C     WRITE(6,*)' ZETO'
C     CALL MATOUT (ZETO,NO,LFZO,NO,LFZO,6)
C
C     WRITE(6,*)'  '
C     WRITE(6,*)' ZETV'
C     CALL MATOUT (ZETV,NV,LFZV,NV,LFZV,6)
C
c     WRITE(6,*)'  '
c     WRITE(6,*)' THE ONE-PDM IN SCF ORDERING in ine lagr'
c     CALL PRINT  (QDM,NT2,NT,6)
C
      RETURN
      END
