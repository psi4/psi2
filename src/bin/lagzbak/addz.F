C
C-----------------------------------------------------------------------
C
      SUBROUTINE ADDZ (QDMT,QDM,
     .                 ORBSYM,NIRRED,FLOV,ITR,ITV,IPQ,NO,NV,NT,
     .                 CCTSCF,SCFTCC,ZETA,ZETO,ZETV,FZT,IFZO,IFZV)
      IMPLICIT INTEGER (A-Z)
      REAL*8 QDMT(*),QDM(*),ZETA(*),ZETO(*),ZETV(*)
      INTEGER ORBSYM(*),FLOV(NIRRED,2),ITR(*),ITV(*),IPQ(*),
     .        CCTSCF(*),SCFTCC(*),FZT(*),IFZO(*),IFZV(*)
C
      NT2=NT*(NT+1)/2
      DO 210 I=1,NT2
         QDMT(I)=QDM(I)
  210 CONTINUE
C
C >>> ADD HF CONTRIBUTION TO ONEPDM IN MO BASIS
C
      DO 220 I=1,NO
         II=CCTSCF(I)
         IIII=IPQ(II)+II
         QDMT(IIII)=QDMT(IIII)+2.0D0
  220 CONTINUE
C
C     WRITE(6,*)
C     WRITE(6,*)'(HF + ACTIVE CORRELATED) ONEPDM IN MO BASIS SCF ORDER'
C     CALL PRINT(QDMT,NT2,NT,6)
C
C >>> ADD CPHF CONTRIBUTION TO ONEPDM IN MO BASIS
C
      DO 110 A=1,NV
         AA=CCTSCF(A+NO)
         ASYM=ORBSYM(AA)
         DO 105 I=1,NO
             II=CCTSCF(I)
             ISYM=ORBSYM(II)
             IF(ASYM.NE.ISYM)GO TO 105
             AI=ITR(A)+I
             AAII=IPQ(MAX0(AA,II))+MIN0(AA,II)
             QDMT(AAII)=QDMT(AAII)+0.5D0*ZETA(AI)
  105    CONTINUE
  110 CONTINUE
C
C >>> ADD Z VECTOR CONTRIBUTIONS FROM FZO-DOC AND FZV-VIR ROTATIONS
C
      DO 360 J=1,NO
         JJ=CCTSCF(J)
         IF(FZT(JJ).NE.1)GO TO 360
         JSYM=ORBSYM(JJ)
         DO 350 K=1,NO
            KK=CCTSCF(K)
            IF(FZT(KK).EQ.1)GO TO 350
            KSYM=ORBSYM(KK)
            IF(KSYM.NE.JSYM)GO TO 350
            JK=IPQ(MAX0(JJ,KK))+MIN0(JJ,KK)
            JKZ=ITR(J)+K
            ZJK=IFZO(JKZ)
            QDMT(JK)=QDMT(JK)+0.5D0*ZETO(ZJK)
  350    CONTINUE
  360 CONTINUE
C
      DO 460 A=1,NV
         AA=CCTSCF(A+NO)
         IF(FZT(AA).NE.1)GO TO 460
         ASYM=ORBSYM(AA)
         DO 450 B=1,NV
            BB=CCTSCF(B+NO)
            IF(FZT(BB).EQ.1)GO TO 450
            BSYM=ORBSYM(BB)
            IF(BSYM.NE.ASYM)GO TO 450
            AB=IPQ(MAX0(AA,BB))+MIN0(AA,BB)
            ABZ=ITV(A)+B
            ZAB=IFZV(ABZ)
            QDMT(AB)=QDMT(AB)+0.5D0*ZETV(ZAB)
  450    CONTINUE
  460 CONTINUE
C
C
C     WRITE(6,*)
C     WRITE(6,*)' TOTAL ONEPDM ( + CPHF TERMS) IN MO BASIS SCF ORDER'
C     CALL PRINT(QDMT,NT2,NT,6)
C
      RETURN
      END
