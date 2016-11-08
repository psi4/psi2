C
C---------------------------------------------------------------------
C
      SUBROUTINE FOCK2 (NT,HH,FOCK,EREF,ENUC,EIG,
     .                  IPQ,ITR,ITV,NO,NV,AIN,EIN,CIN,DIN,
     .                  AOFF,AADD,EOFF,EADD,COFF,CADD,DOFF,DADD,
     .                  ORBSYM,FLOV,NIRRED,I79,FZO,FZV)
      IMPLICIT INTEGER (A-Z)
      REAL*8 HH(*),FOCK(*),AIN(*),EIN(*),CIN(*),DIN(*),EIG(*)
      REAL*8 EREF,ENUC
      INTEGER IPQ(*),ITR(*),ITV(*),FZO(*),FZV(*),
     .        AOFF(*),AADD(*),EOFF(*),EADD(*),COFF(*),CADD(*),
     .        DOFF(*),DADD(*),ORBSYM(*),FLOV(NIRRED,4)
C
      integer prcntr
      NT2=NT*(NT+1)/2
C
C >>> READ ONE-E MO INTS FROM FILE79
C
      ITAP79=79
COK   CALL RFILE(ITAP79)
      CALL WREADW(ITAP79,HH,INTOWP(NT2),1,PON79)
C
C     WRITE(6,*)'ONE E HAMILTONIAN MO BASIS'
C     CALL PRINT(HH,NT2,NT,6)
C
      EREF=0.0D0
      DO 200 I=1,NO
         FII=IPQ(I)+I
         EREF=EREF+HH(FII)
  200 CONTINUE
C
C >>> ADD FROZEN ORBITALS CONTRIBUTION OBTAINED IN RDINTS
C
C
C     WRITE(6,*)'FROZEN MOS CONTRIBUTION '
C     WRITE(6,*)'  '
C     CALL PRINT(FOCK,NT2,NT,6)
C
      DO 220 IJ=1,NT2
         HH(IJ)=HH(IJ)+FOCK(IJ)
  220 CONTINUE
C
C >>> 2 ELECTRON CONTRIBUTION TO FOCK MATRIX
C
      CALL ZERO (FOCK,NT2)
      DO 248  I=1,NO
         IF(FZO(I).EQ.1)GO TO 248
         ISYM=ORBSYM(I)
         FM=FLOV(ISYM+1,1)
         LM=FLOV(ISYM+1,2)
         IF(LM.GT.I)LM=I
         DO 246  M=FM,LM
            IM=ITR(I)+M
            FIM=IPQ(I)+M
            DO 240 J=1,NO
               IF(FZO(J).EQ.1)GO TO 240
               IJ=ITR(MAX0(I,J))+MIN0(I,J)
               MJ=ITR(MAX0(M,J))+MIN0(M,J)
               JJ=ITR(J)+J
               IMJJ=AOFF(MAX0(IM,JJ))+AADD(MIN0(IM,JJ))
               IJMJ=AOFF(MAX0(IJ,MJ))+AADD(MIN0(IJ,MJ))
               FOCK(FIM)=FOCK(FIM)+2.D0*AIN(IMJJ)-AIN(IJMJ)
  240       CONTINUE
  246    CONTINUE
  248 CONTINUE
C
C     WRITE(6,*)'A INTS CONTRIBUTION '
C     WRITE(6,*)'  '
C     CALL PRINT(FOCK,NT2,NT,6)
C
      DO 221 IJ=1,NT2
         HH(IJ)=HH(IJ)+FOCK(IJ)
  221 CONTINUE
C
      CALL ZERO (FOCK,NT2)
      DO 448  I=1,NV
         IF(FZV(I).EQ.1)GO TO 448
         ISYM=ORBSYM(I+NO)
         FM=FLOV(ISYM+1,1)
         LM=FLOV(ISYM+1,2)
         DO 446  M=FM,LM
            IM=ITR(I)+M
            FIM=IPQ(I+NO)+M
            DO 440  J=1,NO
               IF(FZO(J).EQ.1)GO TO 440
               IJ=ITR(I)+J
               MJ=ITR(MAX0(M,J))+MIN0(M,J)
               JJ=ITR(J)+J
               IMJJ=EOFF(IM)+EADD(JJ)
               IJMJ=EOFF(IJ)+EADD(MJ)
               FOCK(FIM)=FOCK(FIM)+2.D0*EIN(IMJJ)-EIN(IJMJ)
  440       CONTINUE
  446    CONTINUE
  448 CONTINUE
C
      DO 230 IJ=1,NT2
         HH(IJ)=HH(IJ)+FOCK(IJ)
  230 CONTINUE
C
C     WRITE(6,*)'E INTS CONTRIBUTION '
C     WRITE(6,*)'  '
C     CALL PRINT(FOCK,NT2,NT,6)
C
      CALL ZERO (FOCK,NT2)
      DO 548  I=1,NV
         IF(FZV(I).EQ.1)GO TO 548
         ISYM=ORBSYM(I+NO)
         FM=FLOV(ISYM+1,3)-NO
         LM=FLOV(ISYM+1,4)-NO
         IF(LM.GT.I)LM=I
         DO 546  M=FM,LM
            FIM=IPQ(I+NO)+M+NO
            IM=ITV(I)+M
            DO 540  J=1,NO
               IF(FZO(J).EQ.1)GO TO 540
               JJ=ITR(J)+J
               IMJJ=COFF(IM)+CADD(JJ)
               IJ=ITR(I)+J
               MJ=ITR(M)+J
               IJMJ=DOFF(IJ)+DADD(MJ)
               FOCK(FIM)=FOCK(FIM)+2.D0*CIN(IMJJ)-DIN(IJMJ)
  540       CONTINUE
  546    CONTINUE
  548 CONTINUE
C
C     WRITE(6,*)'C-D INTS CONTRIBUTION '
C     WRITE(6,*)'  '
C     CALL PRINT(FOCK,NT2,NT,6)
C
      DO 242 IJ=1,NT2
         HH(IJ)=HH(IJ)+FOCK(IJ)
  242 CONTINUE
C
      DO 250 IJ=1,NT2
         FOCK(IJ)=HH(IJ)
  250 CONTINUE
C
C >>> PLUG MO ENERGIES INTO FOCK MATRIX
C
C     CALL ZERO(FOCK,NT2)
C     DO 620 I=1,NT
C        FII=IPQ(I)+I
C        FOCK(FII)=EIG(I)
C 620 CONTINUE
C
c     WRITE(6,*)'FOCK MATRIX IN MO BASIS'
c     WRITE(6,*)'  '
c     CALL PRINT(FOCK,NT2,NT,6)
C
      DO 700 I=1,NO
         FII=IPQ(I)+I
         EREF=EREF+FOCK(FII)
  700 CONTINUE
      EREF=EREF+ENUC
      if (prcntr('IS_ON BRIEF').eq.0) WRITE(6,638)EREF
  638 FORMAT(/,' REFERENCE ENERGY  =  ',F20.12)
C
C >>> WRITE FOCK MATRIX IN FILE79 AFTER THE ONE-E MO INTS
C
      PON79=1+INTOWP(NT2)
      CALL WWRITW(ITAP79,FOCK,INTOWP(NT2),PON79,PON79)
      I79=PON79
COK   CALL RCLOSE(ITAP79,3)
C
      RETURN
      END
