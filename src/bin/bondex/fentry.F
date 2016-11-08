cets The following is called by a C program which handles dynamic
cets allocation of memory.
      subroutine fentry(cc,ia,maxcor)
cets  PROGRAM BONDEX
C     THE BOND-INDEX AND VALENCE  PROGRAM
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
c   by:  Edward Seidl (cets072891)
c   date: July 28, 1991
c   reason: change dimension of ex,ey,and ez for f&g functions
c           make ityfc=5 the default
c
c   by:  Edward Seidl (cets030991)
c   date: March 9, 1991
c   reason: use new input format
C**********************************************************************
c   by:  Edward Seidl (cets)
c   date: August 11, 1988
c   reason: change to fentry for dynamic core allocation
c**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c3-18-88  *
C*  DATE:  MARCH  18,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C*          and remove variable "GY" in block data section            *
C**********************************************************************
C*  BY:  GUSTAVO E. SCUSERIA                                          *
C*  DATE:  APRIL 30, 1987                                             *
C*  REASON:  BUG FOR TCSCF CASE FIXED                                 *
C**********************************************************************
C*  BY:  GUSTAVO E. SCUSERIA                                          *
C*  DATE:  MAY 18, 1987                                               *
C*  REASON:  PRINT OUT GEOMETRY                                       *
C**********************************************************************
C
C
C     *****************************************************
C     *****************************************************
C     **  A PROGRAM TO CALCULATE SOPHISTICATED NUMBERS ? **
C     **                                                 **
C     **            BROUGHT TO YOU BY                    **
C     **                                                 **
C     **                 GES                             **
C     *****************************************************
C     *****************************************************
C
C     PLEASE NOTIFY ANY WEIRD RESULTS OR BUGS.
C
C     INPUT: SET UP # BONDEX ## IN YOUR INPUT FILE
C            ISCFCI, IPRNT, ITYFC (3I5)
C            ISCFCI= 0,1 SCF CALCULATION
C                    2   CI CALCULATION
C            IPRNT = 0,1 MINIMUM PRINT OUT
C                    2   EXTENDED PRINT OUT INCLUDING DENSITY MATRIX
C            ITYFC = 3   SPD FUNCTIONS (INPUT$ CODES)
C                  = 5   IF F OR G ARE PRESENT (PINPUT$ CODES)
C
C     NOTE!! IF USING  ISCFCI=2 , 3RD FLAG IN # ONEPDM ## MUST BE .GE. 1
C            AND FIRST FLAG .EQ. 1 (YES, THE PRINT OPTION!)
C            YOU'LL ALSO NEED TO RUN MASTER$ BEFORE ONEPDM$
C            (AND INTS$ BEFORE MASTER$, SORRY)
C
C     IF USING MORE THAN 200 BFS OR 112 ATOMS TAKE A LOOK AT THE FORTRAN
C     LAST UPDATE : APRIL 7, 1987
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      integer wpadti
cets030991
#include <error.h>
      integer errcod, frdc, frdi, frdboo
      character*8 wfn
      integer ip, prcntr, prdup
      integer angmom
c3-18-88  DIMENSION CC(360000),IA(1)
cets  dimension cc(maxcor),ia(intowp(maxcor))
      dimension cc(maxcor),ia(1)
cets  DIMENSION CC(360000),IA(1)
      DIMENSION I30(200),A30(100)
      COMMON/BASIS/NBASIS,NBFAO,NBFSO,NBATRI
      COMMON/COORD/CHGS(112),X(3,112),WW(112)
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
    1 FORMAT(//,2X,' THE BOND-INDEX AND VALENCE PROGRAM'/
     .         ,2X,'       GES VERSION 05.18.87        '/)
    2 FORMAT(16I5)
    3 FORMAT(/,2X,' PARAMETERS'/
     * 2X,' ISCFCI = ',I8/
     * 2X,' IPRNT  = ',I8/
     * 2X,' ITYFC  = ',I8/
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
    7 FORMAT(//,2X,' WARNING! NO INPUT TO BONDEX HAS BEEN FOUND. '/
     1          2X,'          DEFAULT PARAMETERS WILL BE USED. ',/)
C
c     CALL TSTART(6)
cets  CALL NOUNFL
cets030991      call drum
      call psinit('APPEND')
C
      ITAPE3=3
      INPUT=5
      ITAPE6=6
      ITAP30=30
cets  MAXCOR=360000
c3-18-88  MAXCOR=360000
      call tstart(itape6)
      CALL RFILE(ITAP30)
C
      WRITE(6,1)
      WRITE(3,1)
C
C   ISCFCI IS A PARAMETER FOR TYPE OF WAVE FUNCTION
C     ISCFCI = 0 OR 1  FOR SCF
C     ISCFCI = 2       FOR CI
C   IPRINT IS A PARAMETER FOR SIZE OF OUTPUT
C   ITYFC  IS A PARAMETER FOR DEALING WITH SPD (3) OR FG (5) BFS.
C
C   SET UP DEFAULT VALUES IF INPUT TO PROPER IS NOT FOUND
C
      ISCFCI = 1
      IPRINT = 0
      ITYFC  = 5
C
      CALL LOCATE(INPUT,'# BONDEX #',IERR)
C
cets030991      IF(IERR .NE. 0)  WRITE(6,7)
cets030991      IF(IERR .NE. 0)  GOTO 100
      if(ierr.eq.0) then
        READ(5,2) ISCFCI,IPRNT, ITYFC
        IF(ISCFCI.EQ.0)ISCFCI=1
        IF(ITYFC .EQ.0)ITYFC=5
      else
        iprnt=0
        errcod = frdi('IPRINT',iprnt)
    
        errcod = frdc('WFN',wfn)
        if(wfn.ne.'SCF     ') iscfci=2

        angmom=1
        errcod = frdboo('F_OR_G',angmom)
        if(angmom.eq.1) ityfc=5
      endif

cets030991 print control, not used for now, someone else can fix this up
c          later

      errcod = prcntr('ON BRIEF')
      errcod = prcntr('OFF VERBOSE DEBUG')
      errcod = prdup()
      
      if(prcntr('IS_ON VERBOSE').eq.1.or.prcntr('IS_ON DEBUG').eq.1)
     +  errcod = prcntr('OFF BRIEF')
c
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
      WRITE(6,3) ISCFCI,IPRNT,ITYFC,NBASIS,NBFAO,NBFSO,NBATRI,NATOM,
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
      WRITE(6,994)
  994 FORMAT(//,'GEOMETRY FROM FILE30',/)
      DO 66 I=1,NATOM
      WRITE(6,993)I,X(1,I),X(2,I),X(3,I)
  993 FORMAT(1X,I4,5X,3(F12.8,5X))
   66 CONTINUE
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
      CALL BASSET(CC(IC1),ITYFC)
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
C   CALCULATE THE DAM INDEX
CC    WRITE(6,27)
      WRITE(3,27)
   27 FORMAT(//,2X,' NOW YOU ARE IN VINDEX'/)
      IC1=1
      IC2=IC1+NBFAO*NBFSO
      IC3=IC2+NBASIS
      IC4=IC3+NBASIS
      IC5=IC4+NBFSO*NBFSO
      IC6=IC5+NBFAO*NBFSO
      IOV=IC6+NBATRI
      IBON=IOV+NBATRI
      IPRO=IBON+NBASIS*NBASIS
      IVIN=IPRO+NBASIS*NBASIS
      IV1 =IVIN+NATOM*NATOM
      IE1 =IV1 +NBASIS*NBASIS
      IGAMMA=IE1+NBASIS
      IPOSV=WPADTI(IGAMMA)
CCC     IPOSV=IGAMMA*2-1
      IBETA=IGAMMA+NBASIS
      IVPOS=WPADTI(IBETA)
      IBETAQ=IBETA+NBASIS
      IORD=WPADTI(IBETAQ)
      IW=IBETAQ+NBASIS
      IMINX=IW+NBASIS
      ICMAX=IMINX+NBASIS
      IF(ICMAX.GT.MAXCOR) GO TO 399
      CALL DIPORB(CC(IOV),CC(IC6),CC(IBON),CC(IPRO),CC(IVIN),ISCFCI,
     .            CC(IV1),CC(IE1),CC(IMINX),
     .            CC(IGAMMA),CC(IBETA),CC(IBETAQ),CC(IW),
     .            IA(IPOSV),IA(IVPOS),IA(IORD))
      GO TO 400
  399 WRITE(6,6) ICMAX,MAXCOR
  400 CONTINUE
      CALL TSTOP(itape6)
C
      STOP
      END

