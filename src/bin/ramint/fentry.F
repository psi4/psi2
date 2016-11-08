cets The following is called by a C program which handles dynamic
cets allocation of memory.
      subroutine fentry(cc,ia,maxcor)
cets  PROGRAM RAMINT
C    THE RAMAN INTENSITY PROGRAM
c
c    edward seidl (cets030991)
c    march 9, 1991
c    use new input format
c
C*****************************************************
C***LAST UPDATED ON MAY 24, 1985 BY YUKIO YAMAGUCHI***
C*****************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
cets030991
#include <error.h>
      integer errcod, frdc, frdi
      integer ip, prcntr, prdup, idegn
      character*8 optyp
cets  DIMENSION CC(400000),IA(1)
      dimension cc(maxcor),ia(1)
      CHARACTER*8  SCFTYP,CITYP,NOSYM
      COMMON/BASIS/NBASIS,NBFAO,NBATRI,NBATR2
      COMMON/CALCS/ISTYPE,ICTYPE
      COMMON/COORD/ZAN(50),X(3,50)
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOM,N3N,NATRI,NSYM
      COMMON/GAUSA/EX(200),CS(200)
      COMMON/GAUSB/KATOM(200),KTYPE(200),KSTART(200),KLOC(200),
     1             KMIN(200),KMAX(200),KNG(200)
      COMMON/GAUSC/T(36),PTR(3,3,8),DTR(6,6,8)
      COMMON/GAUSD/INVT(8),ISO(200,8),LBLAT(2,50),ICT(50,8)
      COMMON/GAUSE/MXSPDF(50),MINSH(3,50),MAXSH(3,50)
      COMMON/MAXDS/NMAXD
      COMMON/MFSEC/MASTER,NSECT
      COMMON/NDERS/NDERIV
      COMMON/PARA1/NSHELL,NT,MXT,NISO
      COMMON/PARA2/NPRIM
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/TAPES/ITAP30,ITAP44
cets  EQUIVALENCE (CC,IA)
    1 FORMAT(//,2X,' THE RAMAN INTENSITY PROGRAM'/)
    2 FORMAT(3(A8,2X))
    3 FORMAT(//,3X,2(A8,2X))
    4 FORMAT(5I5)
    5 FORMAT(A8,2X)
    6 FORMAT(//,2X,' WRONG INPUT INFORMATION'/
     1          2X,' SCFTYP = ',A8,5X,' CITYP = ',A8,5X/)
    7 FORMAT(//,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
    8 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
C
cets030691      call drum
      call psinit('APPEND')
      CALL TSTART(6)
cets  CALL NOUNFL
C
      ITAPE3=3
      INPUT=5
      ITAP30=30
      ITAP44=44
      NMAXD=10
cets  MAXCOR=400000
C
      CALL LOCATE(INPUT,'# RAMINT #',IERR)
C
      IOFF(1)=0
      DO 101 I=1,255
  101 IOFF(I+1)=IOFF(I)+I
C
      WRITE(6,1)
C
C   CLSCF--+
C   HSOP---+
C   OPSNG--+--SCF
C   GRSCF--+
C   TCSCF--+
      if(ierr.eq.0) then
        READ(5,2) SCFTYP,CITYP
        WRITE(6,3) SCFTYP,CITYP
        READ(5,4) IPRNT
        READ(5,5) NOSYM
      else
        nosym = 'NOSYM   '

        cityp='SCF     '
        errcod = frdc('WFN',cityp)
        if(cityp.ne.'SCF     ') then
          write(6,*) 'this code only does scf wavefunctions'
          call qabort
        endif

        scftyp='CLSCF   '
        optyp = 'NONE    '
        errcod= frdc('OPENTYPE',optyp)
        if(optyp.ne.'NONE    ') scftyp='GRSCF   '

        errcod = frdboo('DEGEN',idegn)
        if(idegn.eq.1) scftyp = 'GRSCF   '

        write(6,3) scftyp,cityp
      endif

C
      ISTYPE=0
      IF(SCFTYP.EQ.'CLSCF   ') ISTYPE=1
      IF(SCFTYP.EQ.'HSOP    ') ISTYPE=2
      IF(SCFTYP.EQ.'OPSNG   ') ISTYPE=2
      IF(SCFTYP.EQ.'GRSCF   ') ISTYPE=2
      IF(SCFTYP.EQ.'TCSCF   ') ISTYPE=3
      IF(SCFTYP.EQ.'MCSCF   ') ISTYPE=4
      ICTYPE=0
      IF(CITYP.EQ.'SCF     ') ICTYPE=1
      IF(CITYP.EQ.'MCSCF   ') ICTYPE=2
      IF(CITYP.EQ.'CI      ') ICTYPE=3
      NDERIV=1
C
      IF(ISTYPE*ICTYPE.EQ.0) GO TO 399
C
C================================================
C===READ IN BASIS SET AND GEOMETRY FROM TAPE30===
C================================================
C
      CALL GETBAS
C
C===========================================
C===FORCE CALCULATION TO NOSYM IF DESIRED===
C===========================================
C
      IF(NOSYM.NE.'NOSYM   ') GO TO 201
      CALL SETNOS
C
C========================
C===FORM LIMITS ARRAYS===
C========================
C
  201 CONTINUE
      CALL LIMITS
C
C==================================================
C===FORM DENSITY AND LAGRANGIAN MATRICES FOR SCF===
C==================================================
C
      CALL READMF
C
      IC1=1
      IC2=IC1+NBATRI*NTYPEP
      IC3=IC2+NBATRI
      IC4=IC3+NBASIS
      IC5=IC4+NBASIS
      ICMAX=IC5+NBFAO*NBASIS
      WRITE(3,7) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 499
C.................DS      DST     EIG     OCC     EAO.....
      CALL DENS40(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5))
C
C   CALCULATE ONE-ELECTRON INTEGRALS AS A TEST
      IC1=1
      IC2=IC1+NBATRI*NTYPEP
      IC3=IC2+NBATRI
      IC4=IC3+NBASIS
      IC5=IC4+NBASIS
      IC6=IC5+NBFAO*NBASIS
      IC7=IC6+NBATRI
      IC8=IC7+NBATRI
      IC9=IC8+NBATRI
      ICMAX=IC9+NBATRI
      WRITE(3,7) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 499
C.................DST     SS      TT      VV      HH......
      CALL ONEHGT(CC(IC2),CC(IC6),CC(IC7),CC(IC8),CC(IC9))
C
C   CALCULATE FIELD WEIGHTED DENSITY MATRICES
      IC1=1
      IC2=IC1+NBATRI*NTYPEP
      IC3=IC2+NBATRI
      IC4=IC3+NBASIS
      IC5=IC4+NBASIS
      IC6=IC5+NBFAO*NBASIS
      IC7=IC6+NBATRI*3*NTYPEP
      IC8=IC7+NBATRI*3
      ICMAX=IC8+NBASIS*NBASIS*3
      WRITE(3,7) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 499
C.................EAO     DF      DFT     UF......
      CALL DENSEF(CC(IC5),CC(IC6),CC(IC7),CC(IC8))
C
C   CALCULATE TWO-ELECTRON FIRST DERIVATIVES
      IC1=1
      IC2=IC1+NBATRI*NTYPEP
      IC3=IC2+NBATRI
      IC4=IC3+NBASIS
      IC5=IC4+NBASIS
      IC6=IC5+NBFAO*NBASIS
      IC7=IC6+NBATRI*3*NTYPEP
      IC8=IC7+NBATRI*3
      ICMAX=IC8+N3N*NMAXD
      WRITE(3,7) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 499
C.................DS      DST     DF      DFT     D1T.....
      CALL JK1DER(CC(IC1),CC(IC2),CC(IC6),CC(IC7),CC(IC8))
C
C   SYMMETRIZE THE FIRST DERIVATIVE MATRIX
C...............D1T.....
      CALL SYMD(CC(IC8))
C
C   SAVE POLARIZABILITY DERIVATIVES ON THE TAPE18
C.................D1T.....
      CALL FORM18(CC(IC8))
C
      GO TO 500
C
  399 CONTINUE
      WRITE(6,6) SCFTYP,CITYP
      GO TO 500
  499 CONTINUE
      WRITE(6,8) ICMAX,MAXCOR
  500 CONTINUE
C
      CALL RCLOSE(ITAP30,3)
      CALL RCLOSE(ITAP44,3)
      CALL RCLOSE(MASTER,3)
      CALL TSTOP(6)
      STOP
      END
