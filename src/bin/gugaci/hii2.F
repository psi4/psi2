      SUBROUTINE hii2(h2,isym,numi,INT,IJADD,KADD,LADD,H1)
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INTEGER SYMORB,ARR,ORBFRM,ASM,OS,BMAX
      REAL*8 INT
C
      COMMON /COEFS/  A,B,INTAL,INTAU,INTB,INTAD,INTBD
      COMMON /INTS/   NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
      COMMON /DIMS/   NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     #,               NROWOC,LEVFRM
     #,               NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /ALL/    VAL1,VAL2,VAL3,ARR,ITR1,ITR2,IA,JA,ITYPE,ISEGT
     #,               LVFRM1,NLWKI,NLWKJ,IMAX,IMIN
      COMMON /SYMM/   ASM,JS,IS,MX(8),MN(8),OS(8),NUMSYM(8)
C
      DIMENSION INT(NMAX),IJADD(NUMIJ),KADD(SYMORB),LADD(SYMORB)
C     DIMENSION HIIM (1),HIISM(1),HIJM (1),HIM(1),HSM(1)
      DIMENSION H1(1),H2(1)
C
      DATA ACRCY /1.0D-09/
C
cets110990      RETURN
C
C*********************************** HII2 ******************************
C
cets110990      ENTRY HII2(H2,ISYM,NUMI)
C
      LAD=ARR+KADD(OS(ASM+1)+MN(ISYM+1))
C
      IF (ABS(A).LT.ACRCY) GO TO 1110
      IF (ABS(B).LT.ACRCY) GO TO 1120
C
      JISV=1
      DO 1102 I=1,NUMI
         JI=JISV
         IJ=I
         DO 1101 J=1,I-1
            T=B*INT(LAD+INTB)
            H2(IJ)=T+A*INT(LAD+INTAL)
            H2(JI)=T+A*INT(LAD+INTAU)
            LAD=LAD+3
            IJ=IJ+NUMI
            JI=JI+1
 1101    CONTINUE
         H2(IJ)=B*INT(LAD+INTBD)+A*INT(LAD+INTAD)
         LAD=LAD+3
         JISV=JISV+NUMI
 1102 CONTINUE
      RETURN
C
 1110 CONTINUE
      JISV=1
      DO 1112 I=1,NUMI
         JI=JISV
         IJ=I
         DO 1111 J=1,I-1
            T=B*INT(LAD+INTB)
            H2(IJ)=T
            H2(JI)=T
            JI=JI+1
            IJ=IJ+NUMI
            LAD=LAD+3
 1111    CONTINUE
         H2(IJ)=B*INT(LAD+INTBD)
         JISV=JISV+NUMI
         LAD=LAD+3
 1112 CONTINUE
      RETURN
C
 1120 CONTINUE
      JISV=1
      DO 1122 I=1,NUMI
         JI=JISV
         IJ=I
         DO 1121 J=1,I-1
            H2(IJ)=A*INT(LAD+INTAL)
            H2(JI)=A*INT(LAD+INTAU)
            JI=JI+1
            IJ=IJ+NUMI
            LAD=LAD+3
 1121    CONTINUE
         H2(IJ)=A*INT(LAD+INTAD)
         JISV=JISV+NUMI
         LAD=LAD+3
 1122 CONTINUE
C
      RETURN
C
      END
