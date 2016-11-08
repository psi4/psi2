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
      NT2=NT*(NT+1)/2
C
C >>> READ ONE-E MO INTS FROM FILE79
C
c     ITAP79=79
COK   CALL RFILE(ITAP79)
c     CALL WREADW(ITAP79,HH,INTOWP(NT2),1,PON79)
C
C     WRITE(6,*)'ONE E HAMILTONIAN MO BASIS'
C     CALL PRINT(HH,NT2,NT,6)
C
c     EREF=0.0D0
c     DO 200 I=1,NO
c        FII=IPQ(I)+I
c        EREF=EREF+HH(FII)
c 200 CONTINUE
C
C >>> ADD FROZEN ORBITALS CONTRIBUTION OBTAINED IN RDINTS
C
C
C     WRITE(6,*)'FROZEN MOS CONTRIBUTION '
C     WRITE(6,*)'  '
C     CALL PRINT(FOCK,NT2,NT,6)
C
c     DO 220 IJ=1,NT2
c        HH(IJ)=HH(IJ)+FOCK(IJ)
c 220 CONTINUE
C
C >>> 2 ELECTRON CONTRIBUTION TO FOCK MATRIX
C
c     CALL ZERO (FOCK,NT2)
c     DO 248  I=1,NO
c        IF(FZO(I).EQ.1)GO TO 248
c        ISYM=ORBSYM(I)
c        FM=FLOV(ISYM+1,1)
c        LM=FLOV(ISYM+1,2)
c        IF(LM.GT.I)LM=I
c        DO 246  M=FM,LM
c           IM=ITR(I)+M
c           FIM=IPQ(I)+M
c           DO 240 J=1,NO
c              IF(FZO(J).EQ.1)GO TO 240
c              IJ=ITR(MAX0(I,J))+MIN0(I,J)
c              MJ=ITR(MAX0(M,J))+MIN0(M,J)
c              JJ=ITR(J)+J
c              IMJJ=AOFF(MAX0(IM,JJ))+AADD(MIN0(IM,JJ))
c              IJMJ=AOFF(MAX0(IJ,MJ))+AADD(MIN0(IJ,MJ))
c              FOCK(FIM)=FOCK(FIM)+2.D0*AIN(IMJJ)-AIN(IJMJ)
c 240       CONTINUE
c 246    CONTINUE
c 248 CONTINUE
c
c     WRITE(6,*)'A INTS CONTRIBUTION '
c     WRITE(6,*)'  '
c     CALL PRINT(FOCK,NT2,NT,6)
c
c     DO 221 IJ=1,NT2
c        HH(IJ)=HH(IJ)+FOCK(IJ)
c 221 CONTINUE
c
c     CALL ZERO (FOCK,NT2)
c     DO 448  I=1,NV
c        IF(FZV(I).EQ.1)GO TO 448
c        ISYM=ORBSYM(I+NO)
c        FM=FLOV(ISYM+1,1)
c        LM=FLOV(ISYM+1,2)
c        DO 446  M=FM,LM
c           IM=ITR(I)+M
c           FIM=IPQ(I+NO)+M
c           DO 440  J=1,NO
c              IF(FZO(J).EQ.1)GO TO 440
c              IJ=ITR(I)+J
c              MJ=ITR(MAX0(M,J))+MIN0(M,J)
c              JJ=ITR(J)+J
c              IMJJ=EOFF(IM)+EADD(JJ)
c              IJMJ=EOFF(IJ)+EADD(MJ)
c              FOCK(FIM)=FOCK(FIM)+2.D0*EIN(IMJJ)-EIN(IJMJ)
c 440       CONTINUE
c 446    CONTINUE
c 448 CONTINUE
c
c     DO 230 IJ=1,NT2
c        HH(IJ)=HH(IJ)+FOCK(IJ)
c 230 CONTINUE
c
c     WRITE(6,*)'E INTS CONTRIBUTION '
c     WRITE(6,*)'  '
c     CALL PRINT(FOCK,NT2,NT,6)
c
c     CALL ZERO (FOCK,NT2)
c     DO 548  I=1,NV
c        IF(FZV(I).EQ.1)GO TO 548
c        ISYM=ORBSYM(I+NO)
c        FM=FLOV(ISYM+1,3)-NO
c        LM=FLOV(ISYM+1,4)-NO
c        IF(LM.GT.I)LM=I
c        DO 546  M=FM,LM
c           FIM=IPQ(I+NO)+M+NO
c           IM=ITV(I)+M
c           DO 540  J=1,NO
c              IF(FZO(J).EQ.1)GO TO 540
c              JJ=ITR(J)+J
c              IMJJ=COFF(IM)+CADD(JJ)
c              IJ=ITR(I)+J
c              MJ=ITR(M)+J
c              IJMJ=DOFF(IJ)+DADD(MJ)
c              FOCK(FIM)=FOCK(FIM)+2.D0*CIN(IMJJ)-DIN(IJMJ)
c 540       CONTINUE
c 546    CONTINUE
c 548 CONTINUE
c
C     WRITE(6,*)'C-D INTS CONTRIBUTION '
C     WRITE(6,*)'  '
C     CALL PRINT(FOCK,NT2,NT,6)
C
c     DO 242 IJ=1,NT2
c        HH(IJ)=HH(IJ)+FOCK(IJ)
c 242 CONTINUE
C
c     DO 250 IJ=1,NT2
c        FOCK(IJ)=HH(IJ)
c 250 CONTINUE
C
C >>> PLUG MO ENERGIES INTO FOCK MATRIX
C
      CALL ZERO(FOCK,NT2)
      DO 620 I=1,NT
         FII=IPQ(I)+I
         FOCK(FII)=EIG(I)
  620 CONTINUE
C
c     WRITE(6,*)'FOCK MATRIX IN MO BASIS'
c     WRITE(6,*)'  '
c     CALL PRINT(FOCK,NT2,NT,6)
C
c     DO 700 I=1,NO
c        FII=IPQ(I)+I
c        EREF=EREF+FOCK(FII)
c 700 CONTINUE
c     EREF=EREF+ENUC
c     WRITE(6,*)' '
c     WRITE(6,638)EREF
c 638 FORMAT(/,' REFERENCE ENERGY  =  ',F20.12)
C
C >>> WRITE FOCK MATRIX IN FILE79 AFTER THE ONE-E MO INTS
C
c     PON79=1+INTOWP(NT2)
c     CALL WWRITW(ITAP79,FOCK,INTOWP(NT2),PON79,PON79)
c     I79=PON79
COK   CALL RCLOSE(ITAP79,3)
C
      RETURN
      END
