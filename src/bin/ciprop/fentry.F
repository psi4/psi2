cets The following is called by a C program which handles dynamic
cets allocation of memory.
      subroutine fentry(cc,ia,maxcor)
cets  PROGRAM PROPER
C   THE MULLIKEN POPULATON AND DIPOLE MOMENT PROGRAM
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
c   by: Roger Grev, December 20, 1991             search:  crsg
c       changes to basset.f and popul.f so it works for 5d7f9g.
c       basically changing NBASIS to nbfao.
c
c       by: Edward Seidl (cets072891)
c     date: July 28, 1991
c   reason: change dimension of ex,ey and ez in common exyzs and
c           contr in basset to allow for f&g functions
c
c       by: Edward Seidl (cets030291)
c     date: March 2, 1991
c   reason: replace drum with psinit and use new input format
c
c       by: Edward Seidl (cets)
c     date: August 11, 1988
c   reason: change to fentry to allow dynamic core allocation
c**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c3-24-88  *
C*  DATE:  MARCH  24,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C*  LAST UPDATED ON MARCH 20, 1985 BY TJL                             *
C**********************************************************************
C*  LAST UPDATED ON MAY 27, 1985  BY RICHARD REMINGTON                *
C*   CALL TO RFILE(ITAP30) WAS REPLACED BY SREW AND THE               *
C*   RFILE WAS MOVED TO AFTER ITAP30=30, ALSO DEFAULTS                *
C*   WERE ADDED IN CASE NO INPUT TO PROPER IS FOUND.                  *
C**********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
cets030291
#include <error.h>
      integer errcod, frdi
c
c3-24-88  DIMENSION CC(360000),IA(1)
cets  DIMENSION CC(360000),IA(1)
      dimension cc(maxcor),ia(1)
      DIMENSION I30(200),A30(100)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI
      COMMON/COORD/CHGS(112),X(3,112),W(112)
      COMMON/FLAGS/MPOINT,MCONST
      COMMON/FUNCS/NATOM,N3N,NSYM
      COMMON/GAUSS/ZS(200),CS(200)
      COMMON/NORBS/NAST(112),NAED(112)
      COMMON/NSHEL/KNUC(200),KTYPE(200),KPRIM(200),KLOC(200),
     1 KMIN(200),KMAX(200),KNG(200),NSHELL,NPRIM
      COMMON/PARAS/MXCOEF,LOCVEC
      COMMON/QVALS/QORB(256),QATM(112),QNET(112)
      COMMON/SIGNS/IOFF(500),IPRNT
      COMMON/CI101/IOPEN,IOCC,JOCC,KOCC
cets  EQUIVALENCE (CC,IA)
      EQUIVALENCE (I30,A30)
    1 FORMAT(//,2X,' THE MULLIKEN POPULATION AND DIPOLE MOMENT PROGRAM
     1       '/)
    2 FORMAT(3I5)
    3 FORMAT(/,2X,' PARAMETERS'/
     * 2X,' ISCFCI = ',I8/
     * 2X,' ICENT  = ',I8/
     * 2X,' IPRNT  = ',I8/
     * 2X,' NBASIS = ',I8/
     * 2X,' NBFAO  = ',I8/
     * 2X,' NBFSO  = ',I8/
     * 2X,' NBATRI = ',I8/
     * 2X,' NATOM  = ',I8/
     * 2X,' N3N    = ',I8/
     * 2X,' NSYM   = ',I8/
     * 2X,' IOPEN  = ',I8/
     * 2X,' NSHELL = ',I8/
     * 2X,' NPRIM  = ',I8/)
    4 FORMAT(/,3X,20A4)
    5 FORMAT(//,2X,' ENUC  = ',F20.10/
     1          2X,' ESCF  = ',F20.10/)
    6 FORMAT(//,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
    7 FORMAT(//,2X,' WARNING! NO INPUT TO PROPER HAS BEEN FOUND. '/
     1          2X,'          DEFAULT PARAMETERS WILL BE USED. ',/)
C
cets  CALL NOUNFL
cets030291      CALL DRUM
      call psinit('APPEND')
      CALL TSTART(6)
C
      ITAPE3=3
      INPUT=5
      ITAPE6=6
      ITAP30=30
c3-24-88  MAXCOR=360000
cets  MAXCOR=360000
      CALL RFILE(ITAP30)
C
      WRITE(6,1)
      write(6,*) '       F & G function version '
      WRITE(3,1)
C
C   ISCFCI IS A PARAMETER FOR TYPE OF WAVE FUNCTION
C     ISCFCI = 0 OR 1  FOR SCF
C     ISCFCI = 2       FOR CI
C     ISCFCI = 3       FOR MCSCF
C   ICENT IS A PARAMETER FOR REFERENCE COORDINATE OF DIPOLE MOMENT
C     ICENT  = 0 OR 1  FOR CENTER OF MASS
C     ICENT  = 2       FOR ORIGIN OF SPACE FIXED COORDINATE
C     ICENT  = 3       FOR CENTER OF ELECTRONIC CHARGE
C                          BASED ON MULLIKEN POPULATION
C     ICENT  = 4       FOR CENTER OF NUCLEAR CHARGE
C     ICENT  = 5       FOR CENTER OF NET CHARGE
C   IPRINT IS A PARAMETER FOR SIZE OF OUTPUT
C
C   SET UP DEFAULT VALUES IF INPUT TO PROPER IS NOT FOUND
C
      ISCFCI = 2
      ICENT  = 1
      iprnt = 0
C
      CALL LOCATE(INPUT,'# PROPER #',IERR)
C
cets030291      IF(IERR .NE. 0)  WRITE(6,7)
cets030291      IF(IERR .NE. 0)  GOTO 100
      if(ierr.ne.0) then
        errcod = frdi('ICENT',icent)
        errcod = frdi('IPRINT',iprnt)
      else
        READ(5,2) ISCFCI,ICENT,IPRNT
      endif
      IF(ISCFCI.EQ.0) ISCFCI=2
      IF(ICENT.EQ.0) ICENT=1
  100 CONTINUE
      IOFF(1)=0
      DO 101 I=1,499
  101 IOFF(I+1)=IOFF(I)+I
C
      CALL SREW(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)
      MPOINT=I30(2)
      MCONST=I30(3)
      NCALCS=I30(5)
      NBFSO=I30(18)
      NATOM=I30(19)
      NBFAO=I30(22)
      NSHELL=I30(27)
      NPRIM=I30(32)
      NSYM=I30(41)
      MXCOEF=I30(42)
      IOPEN=I30(43)
      NBASIS=NBFSO
      N3N=NATOM*3
      NBATRI=IOFF(NBFAO+1)
C
      WRITE(6,3) ISCFCI,ICENT,IPRNT,NBASIS,NBFAO,NBFSO,NBATRI,NATOM,
     1           N3N,NSYM,IOPEN,NSHELL,NPRIM
C
C   READ IN LABELS
      JUNK=101+MCONST+MPOINT+NCALCS-1
      CALL WREADW(ITAP30,LOCCAL,1,JUNK,JUNK)
      CALL WREADW(ITAP30,I30,60,LOCCAL,LOCCAL)
      WRITE(6,4) (I30(I),I=1,20)
      CALL WREADW(ITAP30,I30,20,LOCCAL,JUNK)
      LOCVEC=I30(1)
C
C   READ IN GEOMETRY
      LOCCAL=LOCCAL+20
      CALL WREADW(ITAP30,X,NATOM*6,LOCCAL,LOCCAL)
C
C   READ IN ENUC AND ESCF
      CALL WREADW(ITAP30,A30,20,LOCCAL,LOCCAL)
      ENUC=A30(1)
      ESCF=A30(2)
      WRITE(6,5) ENUC,ESCF
C
C   READ IN BASIS SET INFORMATION
CC    WRITE(6,21)
      WRITE(3,21)
   21 FORMAT(//,2X,' NOW YOU ARE IN BASSET'/)
      IC1=1
C.................UAS.....
      CALL BASSET(CC(IC1))
C
C   READ IN EIGENVECTORS
CC    WRITE(6,22)
      WRITE(3,22)
   22 FORMAT(//,2X,'NOW YOU ARE IN EIGVEC'/)
      IC1=1
      IC2=IC1+NBFAO*NBFSO
      IC3=IC2+NBASIS
      IC4=IC3+NBASIS
      IC5=IC4+NBFSO*NBFSO
      IC6=IC5+NBFAO*NBFSO
      ICMAX=IC6+NBFAO*NBFSO
      IF(ICMAX.GT.MAXCOR) GO TO 399
C................UAS     EIG     OCC     ESO     EAO     EE......
      CALL EIVEC(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6))
C
C   FORM A DENSITY MATRIX
CC    WRITE(6,23)
      WRITE(3,23)
   23 FORMAT(//,2X,' NOW YOU ARE IN DENSTY'/)
      IC1=1
      IC2=IC1+NBFAO*NBFSO
      IC3=IC2+NBASIS
      IC4=IC3+NBASIS
      IC5=IC4+NBFSO*NBFSO
      IC6=IC5+NBFAO*NBFSO
      ICMAX=IC6+NBATRI
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................OCC     EAO     P.......
      CALL DENSTY(CC(IC3),CC(IC5),CC(IC6),ISCFCI)
C
C   CALCULATE ONE ELECTRON INTEGRALS FOR A TEST
CC    WRITE(6,24)
      WRITE(3,24)
   24 FORMAT(//,2X,' NOW YOU ARE IN ONEHGT'/)
      IC1=1
      IC2=IC1+NBFAO*NBFSO
      IC3=IC2+NBASIS
      IC4=IC3+NBASIS
      IC5=IC4+NBFSO*NBFSO
      IC6=IC5+NBFAO*NBFSO
      IC7=IC6+NBATRI
      IC8=IC7+NBATRI
      IC9=IC8+NBATRI
      IC10=IC9+NBATRI
      ICMAX=IC10+NBATRI
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................P       SS      TT      VV      HH.......
      CALL ONEHGT(CC(IC6),CC(IC7),CC(IC8),CC(IC9),CC(IC10))
C
C   CALCULATE MULLIKEN POPULATION
CC    WRITE(6,25)
      WRITE(3,25)
   25 FORMAT(//,2X,' NOW YOU ARE IN POPUL'/)
      IC1=1
      IC2=IC1+NBFAO*NBFSO
      IC3=IC2+NBASIS
      IC4=IC3+NBASIS
      IC5=IC4+NBFSO*NBFSO
      IC6=IC5+NBFAO*NBFSO
      IC7=IC6+NBATRI
      IC8=IC7+NBATRI
      IC9=IC8+NBATRI
      ICMAX=IC9+NATOM*NATOM
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................P       S       AOPOP   BATM....
      CALL  POPUL(CC(IC6),CC(IC7),CC(IC8),CC(IC9))
C
C   FIND A REFERENCE CENTER IF NECESSARY
CC    WRITE(6,26)
      WRITE(3,26)
   26 FORMAT(//,2X,' NOW YOU ARE IN CENTER'/)
      CALL CENTER(ICENT)
C
C   CALCULATE DIPOLE MOMENT
CC    WRITE(6,27)
      WRITE(3,27)
   27 FORMAT(//,2X,' NOW YOU ARE IN DIPOLE'/)
      IC1=1
      IC2=IC1+NBFAO*NBFSO
      IC3=IC2+NBASIS
      IC4=IC3+NBASIS
      IC5=IC4+NBFSO*NBFSO
      IC6=IC5+NBFAO*NBFSO
      IC7=IC6+NBATRI
      IC8=IC7+NBATRI
      IC9=IC8+NBATRI
      IC10=IC9+NBATRI
      ICMAX=IC10+NBATRI
      IF(ICMAX.GT.MAXCOR) GO TO 399
C.................P       DIPX    DIPY    DIPZ...........
      CALL DIPOLE(CC(IC6),CC(IC8),CC(IC9),CC(IC10),ICENT,ISCFCI)
      GO TO 400
  399 WRITE(6,6) ICMAX,MAXCOR
  400 CONTINUE
      CALL TSTOP(6)
C
      STOP
      END
