C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZLAG (XLAG,FOCK,ZETA,AIN,CIN,DIN,EIN,
     .                 AOFF,AADD,COFF,CADD,DOFF,DADD,EOFF,EADD,
     .                 NSGOO,NSGOV,NSHOV,NSLOV,BUF,IBUF,LENBUF,
     .                 ORBSYM,NIRRED,FLOV,ITR,ITV,IPQ,NO,NV,NT,
     .                 CCTSCF,SCFTCC,
     .                 ZETO,ZETV,IFZO,IFZV,FZT,
     .                 AUX1,NSYMHF,ITYP,NFZO,NFZV)
      IMPLICIT INTEGER (A-Z)
      REAL*8 XLAG(NT,NT),FOCK(*),ZETA(*),AIN(*),CIN(*),DIN(*),
     .       EIN(*),BUF(*),ZETO(*),ZETV(*),AUX1(NT,NT)
      REAL*8 VAL,DELTA,XINT
      INTEGER ORBSYM(*),FLOV(NIRRED,2),ITR(*),ITV(*),IPQ(*),
     .        CCTSCF(*),SCFTCC(*),IBUF(*),FZT(*),EOFF(*),EADD(*),
     .        AOFF(*),AADD(*),COFF(*),CADD(*),DOFF(*),DADD(*),
     .        IFZO(*),IFZV(*),NFZO(*),NFZV(*)
      CHARACTER*4 ITYP(*)
C
C >>> ADD Z-CPHF TERMS TO LAGRANGIAN.
C
C     SOME NON-SCF TERMS ARE NOT INCLUDED !
C
C     ALL A (VIR+FZV), ALL I (FZO+DOC)
C
      DO 110 A=1,NV
         AA=CCTSCF(A+NO)
         ASYM=ORBSYM(AA)
         DO 105 I=1,NO
             II=CCTSCF(I)
             ISYM=ORBSYM(II)
             IF(ASYM.NE.ISYM)GO TO 105
             DELTA=XLAG(AA,II)-XLAG(II,AA)
             VAL=0.0D0
             DO 103 J=1,NO
                JJ=CCTSCF(J)
                JSYM=ORBSYM(JJ)
                IF(ASYM.NE.JSYM)GO TO 103
                AJ=ITR(A)+J
                FIJ=IPQ(MAX0(I,J))+MIN0(I,J)
                VAL=VAL+FOCK(FIJ)*ZETA(AJ)
  103       CONTINUE
            XLAG(II,AA)=XLAG(II,AA)+VAL +DELTA
            XLAG(AA,II)=XLAG(AA,II)+VAL
  105    CONTINUE
  110 CONTINUE
C
C     A: FZV  AND  B: VIR
C
      DO 210 A=1,NV
         AA=CCTSCF(A+NO)
         ASYM=ORBSYM(AA)
         DO 205 B=1,NV
             BB=CCTSCF(B+NO)
             BSYM=ORBSYM(BB)
             IF(ASYM.NE.BSYM)GO TO 205
             VAL=0.0D0
             DO 203 J=1,NO
                JJ=CCTSCF(J)
                JSYM=ORBSYM(JJ)
                IF(ASYM.NE.JSYM)GO TO 203
                AJ=ITR(A)+J
                FBJ=IPQ(B+NO)+J
                VAL=VAL+FOCK(FBJ)*ZETA(AJ)
  203       CONTINUE
            XLAG(AA,BB)=XLAG(AA,BB)+VAL
C
            IF(FZT(AA).EQ.1.AND.FZT(BB).NE.1)THEN
               ABZ=ITV(A)+B
               ZAB=IFZV(ABZ)
               FBB=IPQ(B+NO)+B+NO
               XLAG(AA,BB)=XLAG(AA,BB)+FOCK(FBB)*ZETV(ZAB)
               FAA=IPQ(A+NO)+A+NO
               XLAG(BB,AA)=XLAG(BB,AA)+FOCK(FAA)*ZETV(ZAB)
            ENDIF
C
  205    CONTINUE
  210 CONTINUE
C
C     ALL K (FZO+DOC), ALL I (FZO+DOC)
C
C
C     CALL RDINS (60,AIN,NSGOO,BUF,IBUF,LENBUF,NO,AOFF,AADD,ITR,ITR)
      CALL RDINS (63,DIN,NSHOV,BUF,IBUF,LENBUF,NO,DOFF,DADD,ITR,ITR)
C     CALL RDINS (64,EIN,NSLOV,BUF,IBUF,LENBUF,NO,EOFF,EADD,ITR,ITR)
C
      DO 390 I=1,NO
         II=CCTSCF(I)
         ISYM=ORBSYM(II)
         DO 385 K=1,NO
             KK=CCTSCF(K)
             KSYM=ORBSYM(KK)
             IF(ISYM.NE.KSYM)GO TO 385
C
             VAL=0.0D0
             DO 303 B=1,NV
                BB=CCTSCF(B+NO)
                BSYM=ORBSYM(BB)
                IF(BSYM.EQ.ISYM) THEN
                   BI=ITR(B)+I
                   FBK=IPQ(B+NO)+K
                   VAL=VAL+FOCK(FBK)*ZETA(BI)
                ENDIF
                DO 301 J=1,NO
                   JJ=CCTSCF(J)
                   JSYM=ORBSYM(JJ)
                   IF(JSYM.NE.BSYM)GO TO 301
                   BJ=ITR(B)+J
                   BI=ITR(B)+I
                   BK=ITR(B)+K
                   IJ=ITR(MAX0(I,J))+MIN0(I,J)
                   IK=ITR(MAX0(I,K))+MIN0(I,K)
                   JK=ITR(MAX0(J,K))+MIN0(J,K)
                   IKBJ=EOFF(BJ)+EADD(IK)
                   IJBK=EOFF(BK)+EADD(IJ)
                   JKBI=EOFF(BI)+EADD(JK)
                   XINT=4.0D0*EIN(IKBJ)-EIN(IJBK)-EIN(JKBI)
                   VAL=VAL+XINT*ZETA(BJ)
  301           CONTINUE
  303        CONTINUE
             XLAG(II,KK)=XLAG(II,KK)+VAL
C
             VAL=0.0D0
             DO 313 A=1,NV
                AA=CCTSCF(A+NO)
                IF(FZT(AA).NE.1)GO TO 313
                ASYM=ORBSYM(AA)
                DO 311 B=1,NV
                   BB=CCTSCF(B+NO)
                   IF(FZT(BB).EQ.1)GO TO 311
                   BSYM=ORBSYM(BB)
                   IF(ASYM.NE.BSYM)GO TO 311
                   BA=ITV(MAX0(B,A))+MIN0(B,A)
                   IK=ITR(MAX0(I,K))+MIN0(I,K)
                   IKBA=COFF(BA)+CADD(IK)
                   BI=ITR(B)+I
                   AK=ITR(A)+K
                   AKBI=DOFF(MAX0(AK,BI))+DADD(MIN0(AK,BI))
                   BK=ITR(B)+K
                   IA=ITR(A)+I
                   IABK=DOFF(MAX0(IA,BK))+DADD(MIN0(IA,BK))
                   XINT=4.0D0*CIN(IKBA)-DIN(IABK)-DIN(AKBI)
                   ABZ=ITV(A)+B
                   ZAB=IFZV(ABZ)
                   VAL=VAL+XINT*ZETV(ZAB)
  311           CONTINUE
  313        CONTINUE
             XLAG(II,KK)=XLAG(II,KK)+VAL
C
             VAL=0.0D0
             DO 323 J=1,NO
                JJ=CCTSCF(J)
                IF(FZT(JJ).NE.1)GO TO 323
                JSYM=ORBSYM(JJ)
                DO 321 L=1,NO
                   LL=CCTSCF(L)
                   IF(FZT(LL).EQ.1)GO TO 321
                   LSYM=ORBSYM(LL)
                   IF(JSYM.NE.LSYM)GO TO 321
                   KI=ITR(MAX0(K,I))+MIN0(K,I)
                   JL=ITR(MAX0(J,L))+MIN0(J,L)
                   KIJL=AOFF(MAX0(KI,JL))+AADD(MIN0(KI,JL))
                   KJ=ITR(MAX0(K,J))+MIN0(K,J)
                   IL=ITR(MAX0(I,L))+MIN0(I,L)
                   KJIL=AOFF(MAX0(KJ,IL))+AADD(MIN0(KJ,IL))
                   KL=ITR(MAX0(K,L))+MIN0(K,L)
                   JI=ITR(MAX0(J,I))+MIN0(J,I)
                   KLJI=AOFF(MAX0(KL,JI))+AADD(MIN0(KL,JI))
                   XINT=4.0D0*AIN(KIJL)-AIN(KJIL)-AIN(KLJI)
                   JLZ=ITR(J)+L
                   ZJL=IFZO(JLZ)
                   VAL=VAL+XINT*ZETO(ZJL)
  321           CONTINUE
  323        CONTINUE
             XLAG(II,KK)=XLAG(II,KK)+VAL
  325        CONTINUE
C
C      I: FZO  AND  K: DOC
C
             IF(FZT(II).EQ.1.AND.FZT(KK).NE.1)THEN
                IKZ=ITR(I)+K
                ZIK=IFZO(IKZ)
                FKK=IPQ(K)+K
                XLAG(II,KK)=XLAG(II,KK)+FOCK(FKK)*ZETO(ZIK)
                FII=IPQ(I)+I
                XLAG(KK,II)=XLAG(KK,II)+FOCK(FII)*ZETO(ZIK)
             ENDIF
 
  385    CONTINUE
  390 CONTINUE
C
c     WRITE(6,*)
c     WRITE(6,*)' TOTAL ( + Z-CPHF ) LAGRANGIAN '
c     WRITE(6,*)
C
      DO 510 ISYM=1,NSYMHF
         FI=FLOV(ISYM,1)-NFZO(ISYM)
         LI=FLOV(ISYM,2)+NFZV(ISYM)
         IS=0
         DO 505 I=FI,LI
            IS=IS+1
            JS=0
            DO 503 J=FI,LI
               JS=JS+1
               AUX1(IS,JS)=XLAG(I,J)
  503       CONTINUE
  505    CONTINUE
c          WRITE(6,*)
c          WRITE(6,*)' SYMMETRY   ', ITYP(ISYM),'*  ',ITYP(ISYM)
c          WRITE(6,*)
c          CALL MATOUT (AUX1,NT,NT,IS,IS,6)
  510 CONTINUE
C
C          CALL MATOUT (XLAG,NT,NT,NT,NT,6)
C
      RETURN
      END
