C
C----------------------------------------------------------------------
C
      SUBROUTINE CPHF (XLAG,CIN,DIN,ZETA,DELX,NO,NV,NT,
     .                 COFF,CADD,DOFF,DADD,ORBSYM,NIRRED,FLOV,
     .                 ITR,ITV,IPQ,BUF,IBUF,LENBUF,NSGOV,NSHOV,CCTSCF,
     .                 AUXOO,FOCK,EX)
      IMPLICIT INTEGER (A-Z)
      REAL*8 CIN(*),DIN(*),XLAG(NT,NT),DELX(*),ZETA(*),BUF(*),
     .       AUXOO(NO,NO),FOCK(*),EX(*)
      INTEGER ORBSYM(*),FLOV(NIRRED,2),ITR(*),ITV(*),IPQ(*),
     .        COFF(*),CADD(*),DOFF(*),DADD(*),IBUF(*),CCTSCF(*)
      REAL*8 EPS,VAL
      REAL*8 DET
C
C >>> CALCULATE PERTURBATION
C
      CALL ZERO(DELX,NO*NV)
      DO 120 A=1,NV
         AA=CCTSCF(A+NO)
         ASYM=ORBSYM(AA)
         DO 110 I=1,NO
            II=CCTSCF(I)
            ISYM=ORBSYM(II)
            IF(ASYM.NE.ISYM)GO TO 110
            AI=ITR(A)+I
            DELX(AI)=XLAG(AA,II)-XLAG(II,AA)-EX(AI)
  110    CONTINUE
  120 CONTINUE
C
C     WRITE(6,*)'  '
C     WRITE(6,*)' DELTA X '
C     CALL MATOUT (DELX,NO,NV,NO,NV,6)
C
C >>> READ C AND D INTS FROM FILES 62 AND 63
C
C     CALL RDINS (63,DIN,NSHOV,BUF,IBUF,LENBUF,NO,DOFF,DADD,ITR,ITR)
C     CALL RDINS (62,CIN,NSGOV,BUF,IBUF,LENBUF,NO,COFF,CADD,ITV,ITR)
C
C >>> CALCULATE HOW MANY A1 TRANSITIONS
C
      LIMA1=0
      DO 160 A=1,NV
         AA=CCTSCF(A+NO)
         ASYM=ORBSYM(AA)
         DO 150 I=1,NO
            II=CCTSCF(I)
            ISYM=ORBSYM(II)
            IF(ISYM.NE.ASYM)GO TO 150
            LIMA1=LIMA1+1
  150    CONTINUE
  160 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*)' NUMBER OF A1-A1 TRANSITIONS IS ',LIMA1
      WRITE(6,*)
C
C >>> FORM A1-A1 BLOCK OF A MATRIX
C
      DO 260 A=1,NV
         AA=CCTSCF(A+NO)
         ASYM=ORBSYM(AA)
         DO 250 B=1,A
            BB=CCTSCF(B+NO)
            BSYM=ORBSYM(BB)
            AB=ITV(A)+B
            CALL ZERO(AUXOO,NO*NO)
            DO 220 I=1,NO
               II=CCTSCF(I)
               ISYM=ORBSYM(II)
               IF(ISYM.NE.ASYM)GO TO 220
               AI=ITR(A)+I
               BI=ITR(B)+I
               LJ=NO
               IF(B.EQ.A)LJ=I
               DO 210 J=1,LJ
                  JJ=CCTSCF(J)
                  JSYM=ORBSYM(JJ)
                  IF(JSYM.NE.BSYM)GO TO 210
                  IJ=ITR(MAX0(I,J))+MIN0(I,J)
                  BJ=ITR(B)+J
                  AJ=ITR(A)+J
                  AIBJ=DOFF(AI)+DADD(BJ)
                  AJBI=DOFF(MAX0(AJ,BI))+DADD(MIN0(AJ,BI))
                  ABIJ=COFF(AB)+CADD(IJ)
C >>> TRIPLET PROPAGATOR
C                 VAL=-DIN(AJBI)-CIN(ABIJ)
C >>> SINGLET NR PROP
C                 VAL=+DIN(AJBI)-CIN(ABIJ)
C >>> SINGLET REAL PROPAGATOR
                  VAL=4.D0*DIN(AIBJ)-DIN(AJBI)-CIN(ABIJ)
                  AUXOO(I,J)=VAL
  210          CONTINUE
  220       CONTINUE
            DO 240 I=1,NO
               II=CCTSCF(I)
               ISYM=ORBSYM(II)
               IF(ISYM.NE.ASYM)GO TO 240
               AI=ITR(A)+I
               LJ=NO
               IF(B.EQ.A)LJ=I
               DO 230 J=1,LJ
                  IJ=ITR(I)+J
                  JJ=CCTSCF(J)
                  JSYM=ORBSYM(JJ)
                  IF(JSYM.NE.BSYM)GO TO 230
                  BJ=ITR(B)+J
                  AIBJ=DOFF(AI)+DADD(BJ)
                  DIN(AIBJ)=AUXOO(I,J)
  230          CONTINUE
  240       CONTINUE
C
  250    CONTINUE
  260 CONTINUE
C
      LIM=LIMA1*(LIMA1+1)/2
      WRITE(6,*)
      WRITE(6,*)' DIMENSION OF PACKED A MATRIX =',LIM
      WRITE(6,*)
C     WRITE(6,*)' A MATRIX INTS PART'
C     CALL PRINT(DIN,LIM,LIMA1,6)
C
      DO 460 A=1,NV
         AA=CCTSCF(A+NO)
         ASYM=ORBSYM(AA)
         DO 450 B=1,A
            BB=CCTSCF(B+NO)
            BSYM=ORBSYM(BB)
            FAB=IPQ(A+NO)+B+NO
            DO 420 I=1,NO
               II=CCTSCF(I)
               ISYM=ORBSYM(II)
               IF(ISYM.NE.ASYM)GO TO 420
               AI=ITR(A)+I
               LJ=NO
               IF(B.EQ.A)LJ=I
               DO 410 J=1,LJ
                  JJ=CCTSCF(J)
                  JSYM=ORBSYM(JJ)
                  IF(JSYM.NE.BSYM)GO TO 410
                  FIJ=IPQ(MAX0(I,J))+MIN0(I,J)
                  BJ=ITR(B)+J
                  AIBJ=DOFF(AI)+DADD(BJ)
                  IF(I.EQ.J)DIN(AIBJ)=DIN(AIBJ)+FOCK(FAB)
                  IF(A.EQ.B)DIN(AIBJ)=DIN(AIBJ)-FOCK(FIJ)
  410          CONTINUE
  420       CONTINUE
  450    CONTINUE
  460 CONTINUE
C
C     WRITE(6,*)
C     WRITE(6,*)' A MATRIX TOTAL'
C     CALL PRINT(DIN,LIM,LIMA1,6)
C
C >>> USE DIRECT METHOD TO INVERT  A MATRIX
C
      CALL DSINV(DIN,LIMA1,EPS,IER,DET)
C
C     WRITE(6,*)
C     WRITE(6,*)' A MATRIX INVERSE'
C     CALL PRINT(DIN,LIM,LIMA1,6)
C
C >>> CALCULATE  Z = DELTA X *  A ** -1
C
      CALL ZERO(ZETA,NO*NV)
      DO 360 A=1,NV
         AA=CCTSCF(A+NO)
         ASYM=ORBSYM(AA)
         DO 350 I=1,NO
            II=CCTSCF(I)
            ISYM=ORBSYM(II)
            IF(ISYM.NE.ASYM)GO TO 350
            AI=ITR(A)+I
            DO 320 B=1,NV
               BB=CCTSCF(B+NO)
               BSYM=ORBSYM(BB)
               DO 310 J=1,NO
                  JJ=CCTSCF(J)
                  JSYM=ORBSYM(JJ)
                  IF(JSYM.NE.BSYM)GO TO 310
                  BJ=ITR(B)+J
                  BJAI=DOFF(MAX0(AI,BJ))+DADD(MIN0(AI,BJ))
                  ZETA(AI)=ZETA(AI)+DELX(BJ)*DIN(BJAI)
  310          CONTINUE
  320       CONTINUE
  350    CONTINUE
  360 CONTINUE
C
C     WRITE(6,*)'  '
c     WRITE(6,*)' Z VECTOR'
c     CALL MATOUT (ZETA,NO,NV,NO,NV,6)
C
      RETURN
      END
