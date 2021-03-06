cets The following is called by a C program which handles dynamic
cets allocation of memory.
      subroutine fentry(cc,ia,maxcor)
cets  PROGRAM RAMANC
C   ANALYTICAL RAMAN INTENSITIES
C   FOR CLOSED SCF WAVE FUNCTIONS
C*****************************************************
C***LAST UPDATED ON MAY 28, 1985 BY YUKIO YAMAGUCHI***
C*****************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
cets  DIMENSION CC(400000),IA(1)
      dimension cc(maxcor),ia(1)
      DIMENSION AFG(6,150)
      COMMON/BASIS/NBASIS,NBFAO,NTRI,NTRI2,NBASQ
      COMMON/CALIF/LPARA(1024),APARA(1024)
      COMMON/DIPOL/DIPDA(135),DIPDN(135)
      COMMON/EAFGS/AFG1(3,3,150),AFG2(3,3,150),AFGT(3,3,150)
      COMMON/FUNCS/NATOMS,N3N,N3NP3,N3NX3
      COMMON/MFSEC/MASTER,NSECT
      COMMON/SECTF/ISF(153),IHF(153),IFF(153),IBF(153),IUF(153)
      COMMON/SIGNS/IOFF(466),IPRNT
      COMMON/TAPES/ITAP43,ITAP44
      COMMON/CI101/IHOMO,ILUMO
cets  EQUIVALENCE (CC,IA)
      DATA A00 / 0.52917706D+00 /
    1 FORMAT(//,2X,' RAMAN INTENSITY PROGRAM FOR CLOSED SHELL SCF IN AO
     1BASIS'/)
    2 FORMAT(5I5)
    3 FORMAT(2X,8I6)
    4 FORMAT(//,2X,' PARAMETERS'/
     * 2X,' NBASIS = ',I8/
     * 2X,' NBFAO  = ',I8/
     * 2X,' NTRI   = ',I8/
     * 2X,' NATOMS = ',I8/
     * 2X,' N3N    = ',I8/
     * 2X,' N3NP3  = ',I8/
     * 2X,' N3NX3  = ',I8/
     * 2X,' IHOMO  = ',I8/
     * 2X,' ILUMO  = ',I8/
     * 2X,' IPRNT  = ',I8/)
    5 FORMAT(//,2X,' SCF EIGENVECTOR (AO BASIS) MATRIX'/)
    6 FORMAT(2I5)
    7 FORMAT(3F20.10)
    8 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
    9 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
   10 FORMAT(//,2X,' AFG MATRIX IN (BOHR+2)'/)
   11 FORMAT(//,2X,' AFG MATRIX IN (A+2)'/)
   12 FORMAT(3F20.10)
C
      call drum
      CALL TSTART(6)
cets  CALL NOUNFL
      CALL INITMF(1)
C
      ITAPE3=3
      ITAPE6=6
      INPUT=5
      ITAP18=18
      ITAP43=43
      ITAP44=44
cets  MAXCOR=400000
cets  open (unit=18,file='file18',status='unknown')
      call ffile(itap18,' ',0)
C
      CALL LOCATE(INPUT,'# RAMAN ##',IERR)
C
C   LOCATION OF MATRICES IN JTAPE1
C   (1) SA (N3NP3*NTRIL)
C   (2) HA (N3NP3*NTRIL)
C   (3) FA (N3NP3*NTRIL)
C   (4) BA (N3NP3*NTRIL)
C   (5) UA (N3NP3*NBASQL)
C
      IOFF(1)=0
      DO 101 I=1,465
  101 IOFF(I+1)=IOFF(I)+I
C
      WRITE(6,1)
      WRITE(3,1)
C
C   READIN PARAMETERS
      READ(5,2) IPRNT
C
C   READIN SCF INFORMATION
      NBASIS=LPARA(3)
      NTRI=LPARA(4)
      NTRI2=NTRI*2
      NTOT=IOFF(NTRI+1)
      IHOMO=LPARA(7)
      NATOMS=LPARA(10)
      NBFAO=LPARA(11)
      ILUMO=IHOMO+1
      N3N=LPARA(17)
      N3NP3=N3N+3
      N3NX3=N3N*3
      NTRIL=((NTRI2-1)/NSECT+1)
      NBASQ=NBASIS*NBASIS*2
      NBASQL=((NBASQ-1)/NSECT+1)
      N3TRL=N3NP3*NTRIL
      DO 102 I=1,N3NP3
      ISF(I)=(I-1)*NTRIL+1
      IHF(I)=(I-1)*NTRIL+N3TRL+1
      IFF(I)=(I-1)*NTRIL+N3TRL+N3TRL+1
      IBF(I)=(I-1)*NTRIL+N3TRL+N3TRL+N3TRL+1
      IUF(I)=(I-1)*NBASQL+N3TRL+N3TRL+N3TRL+N3TRL+1
      WRITE(3,3) I,ISF(I),IHF(I),IFF(I),IBF(I),IUF(I)
  102 CONTINUE
C
      WRITE(6,4) NBASIS,NBFAO,NTRI,NATOMS,N3N,N3NP3,N3NX3,IHOMO,ILUMO,
     1           IPRNT
C
C   READ IN SO-MO EIGENVECTORS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IC1=1
      IC2=IC1+NBASIS
      IC3=IC2+NBASIS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      CALL MREAD(CC(IC1),17)
      CALL MREAD(CC(IC2),18)
      CALL MREAD(CC(IC3),20)
      IF(IPRNT.LE.2) GO TO 301
      WRITE(6,5)
      CALL EIGOUT(CC(IC3),CC(IC1),CC(IC2),NBFAO,NBASIS,NBFAO,NBASIS,6)
C
  301 CONTINUE
      WRITE(3,21)
   21 FORMAT(//,2X,' NOW YOU ARE IN SSMAT'/)
      IC4=IC3+NBFAO*NBASIS
      IC5=IC4+NBASIS*NBASIS*N3N
      IC6=IC5+NBASIS*NBASIS*3
      IC7=IC6+NTRI*N3NX3
      ICMAX=IC7+NTRI*9
      WRITE(3,8) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................UA      UF      SSAF    SSFG....
      CALL SSMAT(CC(IC4),CC(IC5),CC(IC6),CC(IC7))
C
      WRITE(3,22)
   22 FORMAT(//,2X,' NOW YOU ARE IN POLDER'/)
      IC4=IC3+NBFAO*NBASIS
      IC5=IC4+NBASIS*NBASIS*N3N
      IC6=IC5+NBASIS*NBASIS*3
      IC7=IC6+NTRI*N3NX3
      IC8=IC7+NTRI*9
      IC9=IC8+NTRI*N3N
      IC10=IC9+NTRI*3
      ICMAX=IC10+NTRI*N3NX3
      WRITE(3,8) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................EIG     UA      UF      SSAF    SSFG    BA......
      CALL POLDER(CC(IC1),CC(IC4),CC(IC5),CC(IC6),CC(IC7),CC(IC8),
C.................BF      EAF......
     1            CC(IC9),CC(IC10))
C
C   STORE POLARIZABILITY DERIVATIVES (PART I) ON THE TAPE18
      DO 105 IABC=1,N3N
      IXYZ=0
      DO 104 I=1,3
      DO 104 J=1,I
      IXYZ=IXYZ+1
      AFG(IXYZ,IABC)=AFGT(I,J,IABC)
  104 CONTINUE
  105 CONTINUE
      WRITE(6,10)
      CALL MATOUT(AFG,6,150,6,N3N,6)
C
CCC   A22=A00*A00
CCC   DO 106 IABC=1,N3N
CCC   DO 106 IXYZ=1,6
CCC   AFG(IXYZ,IABC)=AFG(IXYZ,IABC)*A22
CC106 CONTINUE
CCC   WRITE(6,11)
CCC   CALL MATOUT(AFG,6,150,6,N3N,6)
C
      REWIND ITAP18
      WRITE(ITAP18,2) NATOMS,N3N
      WRITE(ITAP18,12) ((AFG(IXYZ,IABC),IABC=1,N3N),IXYZ=1,6)
      REWIND ITAP18
      GO TO 400
C
  399 WRITE(6,9) ICMAX,MAXCOR
  400 CONTINUE
C
      CALL RCLOSE(ITAP43,3)
      CALL RCLOSE(ITAP44,3)
      CALL RCLOSE(MASTER,3)
      CALL TSTOP(6)
C
      STOP
      END
