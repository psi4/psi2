      SUBROUTINE hiis(h1,isym,numi,INT,IJADD,KADD,LADD,H2)
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
C************************************ HIIS *****************************
C
cets110990      ENTRY HIIS(H1,ISYM,NUMI)
C
      LAD=ARR+KADD(OS(ASM+1)+MN(ISYM+1))
C
      IF (ABS(A).LT.ACRCY) GO TO 210
      IF (ABS(B).LT.ACRCY) GO TO 220
C
      JISV=1
      DO 202 I=1,NUMI
         JI=JISV
         IJ=I
         DO 201 J=1,I-1
            T=B*INT(LAD+INTB )+A*INT(LAD+INTAL)
            H1(IJ)=T
            H1(JI)=T
            JI=JI+1
            IJ=IJ+NUMI
            LAD=LAD+3
  201    CONTINUE
         H1(IJ)=0.0
         JISV=JISV+NUMI
         LAD=LAD+3
  202 CONTINUE
C
      RETURN
C
  210 CONTINUE
      JISV=1
      DO 212 I=1,NUMI
         JI=JISV
         IJ=I
         DO 211 J=1,I-1
            T=B*INT(LAD+INTB)
            H1(IJ)=T
            H1(JI)=T
            JI=JI+1
            IJ=IJ+NUMI
            LAD=LAD+3
  211    CONTINUE
         H1(IJ)=0.0
         JISV=JISV+NUMI
         LAD=LAD+3
  212 CONTINUE
C
      RETURN
C
  220 CONTINUE
      JISV=1
      DO 222 I=1,NUMI
         JI=JISV
         IJ=I
         DO 221 J=1,I-1
            T=A*INT(LAD+INTAL)
            H1(IJ)=T
            H1(JI)=T
            JI=JI+1
            IJ=IJ+NUMI
            LAD=LAD+3
  221    CONTINUE
         H1(IJ)=0.0
         JISV=JISV+NUMI
         LAD=LAD+3
  222 CONTINUE
C
      RETURN
C
      END
