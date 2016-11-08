cges  PROGRAM DERIVY
      subroutine fentry(cc,ia,maxcor)
C   THE DERIVATIVE PROGRAM
C   THIS PROGRAM IS ABLE TO CALCULATE FIRST AND SECOND DERIVATIVES
C   (FOR ALL KINDS OF WAVE FUNCTIONS ??)
c    to run g-functions derivtaives some trivial changes for size of
c    common blocks are needed. Look in NOTE.
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C   by David Sherrill
C   July 8, 1995
C   add tcgrfl flag so that wavefunction 'TCSCF-GR' can be specified
C   using new input style.
C**********************************************************************
C   BY:  RICHARD REMINGTON                  search: c010692            
C   DATE:  January  06, 1992                                           
C   REASON: Modify new input format to handle open shells correctly    
C           (maybe).
C**********************************************************************
c   by Curtis Janssen
c   march 26, 1991
c   removed some blank lines and various bits of garbage from the
c   output to save trees
C**********************************************************************
c   by Edward Seidl (cets053191)
c   may 31, 1991
c   change way mopn is determined for high spin cases
C**********************************************************************
c   by Edward Seidl (cets030691)
c   march 6, 1991
c   change wrder to write wfn, dertyp, and number of fzv&fzc to file11
c   also add DELETE55 option
c   also change derfk1 to write file42 differently for ci gradients
c   add printing control also
c   change dens40 to read file49 for openshell ci gradients
c
c   by Edward Seidl (cets030491)
c   march 4, 1991
c   replace drum with psinit, and add new input
C**********************************************************************
c   by GES
c   September 27, 1990
c   Bug fix reported by Tracy in initialization of B matrix in jk2der
c   implemented.
C**********************************************************************
c   by GES
c   may 22, 1990
c   get rid of file55 only in ci or cc cases.
c   Also: file42 is created! for scf2nd. Need by cpgraoth.
c         this modification done previously.
C**********************************************************************
C*  BY:  GUSTAVO E. SCUSERIA                                          *
C*  DATE:  november 12, 1989                                          *
c   reason: avoid creating file42 (used by NGRCPHF)
C**********************************************************************
C*  BY:  GUSTAVO E. SCUSERIA                                          *
C*  DATE:  MAY 08, 1989                                               *
C*  REASON: FOR COUPLED CLUSTER OPTION  5D, 7F 9G FIX NBATRI/NTRI     *
C**********************************************************************
C*  BY:  TRACY HAMILTON                                               *
C*  DATE:  APRIL 21, 1989                                             *
C*  REASON: READ IN AND PASS NPERSH ARRAY TO DERFK1 FOR PROPER        *
C           CONSTRUCTION OF THE IBFTAT ARRAY WRITTEN IN FILE42.       *
C**********************************************************************
C*  BY:  GUSTAVO E. SCUSERIA                                          *
C*  DATE:  FEBRUARY 22, 1989                                          *
C*  REASON: INCLUDE COUPLED CLUSTER OPTION  'CC' IN # DERIV ###       *
C           REQUIRES FILE79. COMPATIBLE ONLY WITH LAGZBAK CODE.       *
C*     FOR CHANGES SEARCH CGES                                        *
C**********************************************************************
C*  BY:  TRACY HAMILTON                                               *
C*  DATE:  DECEMBER 1, 1988                                           *
C*  REASON: PUT IN CORRELATED DERIVATIVES                             *
C*     FOR CHANGES SEARCH C12-1-88.                                   *
C*     NTRI,NTRI2 IN COMMON/BASIS/                                    *
C*  NOTE THAT B MATRIX REQUIRES NTYPEP MATRICES (NTYPES+1)            *
C**********************************************************************
C*  BY:  TRACY HAMILTON                                               *
C*  DATE:  NOVEMBER 11, 1988                                          *
C*  REASON: REDUCE CORE, AND CHANGE GETBAS TO ONLY ALLOW UP TO F      *
C*      FUNCTIONS.                                                    *
C**********************************************************************
C*  BY:  TRACY HAMILTON                         SEARCH CTPH           *
C*  DATE:  NOVEMBER 11, 1988                                          *
C*  REASON: PUT IN DERIVATIVES OF D5,D6,F7,F10,G9,G15 FUNCTIONS.      *
C**********************************************************************
C*  BY:  RICHARD REMINGTON                     search: MachineSize    *
C*  DATE:  MARCH  12,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370  AND  ADD A WARNING   *
C*  WITH STOP AND EXIT IF LABEL "# DERIV ##" IS NOT FOUND BY LOCATE.  *
C**********************************************************************
C*  BY:  RICHARD REMINGTON                                            *
C*  DATE:  NOVEMBER 11, 1987                                          *
C*  REASON: IMPROVE PRECISION BY CHANGING DLIMIT = 1.0E-10 TO 1.0E-16 *
C*          CHANGE DLIMIT TO DLIM AND PUT IT INTO COMMON TOLER        *
C**********************************************************************
C* LAST UPDATED ON JUNE 20, 1987 BY RBR---WRITE ESCF TO 11 FOR SCF2ND *
C**********************************************************************
C*   BY: WESLEY ALLEN                                                 *
C*   DATE: NOVEMBER 5, 1986                                           *
C*   REASON:  FIX BUG IN THE FORMATION OF THE OCCUPATION VECTOR FOR   *
C*            TCSCF CASES WHERE THE TWO ORBITALS ARE OF THE SAME      *
C*            SYMMETRY (THIS INCLUDES ALL SCFX CASES).                *
C**********************************************************************
C*  LAST UPDATED ON MAY 15, 1985 BY YUKIO YAMAGUCHI                   *
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C3-12-88  DIMENSION CC(1200000),IA(1),I30(200)
CMachineSize (For 15Mb)  DIMENSION CC(850000),IA(1),I30(200)
CMACHINESIZE (FOR  7MB)  DIMENSION CC(350000),IA(1),I30(200)
      DIMENSION CC(maxcor),IA(*),I30(200),NPERSH(5)
cets030491
#include <error.h>
      integer errcod, frdc, frdi, frdboo, frdia1, fcnt
      integer prdup, prcntr
      character*8 opentp
      logical oldinp
      integer tcgrfl
c
      CHARACTER*8 SCFTYP,CITYP,DERTYP,RAMAN
      CHARACTER*8 NOSYM,AFOCK,DFOCK,DERINT
      CHARACTER*8 OPTYPE
      integer opshls(20),nopshl
      COMMON/BASIS/NBASIS,NBFAO,NBATRI,NBATR2,NTRI,NTRI2
      COMMON/CALCS/ISTYPE,ICTYPE,IRAMAN,IAFOCK,IDFOCK,IDINT
      COMMON/COORD/ZAN(50),X(3,50)
      COMMON/ENRGY/ENUC,ESCF
      COMMON/FUNCS/NTYPES,NTYPEP,NTREAD,NATOM,N3N,NATRI,NSYM
      COMMON/GAUSA/EX(200),CS(200)
      COMMON/GAUSB/KATOM(200),KTYPE(200),KSTART(200),KLOC(200),
     1             KMIN(200),KMAX(200),KNG(200)
CTPH  COMMON/GAUSC/T(36),PTR(3,3,8),DTR(6,6,8)
      COMMON/GAUSC/T(225),PTR(3,3,8),DTR(6,6,8),FTR(10,10,8),
     1             GTR(15,15,8)
      COMMON/GAUSD/INVT(8),ISO(200,8),LBLAT(2,50),ICT(50,8)
CTPH  COMMON/GAUSE/MXSPDF(50),MINSH(3,50),MAXSH(3,50)
      COMMON/GAUSE/MXSPDF(50),MINSH(5,50),MAXSH(5,50)
      COMMON/NDERS/NDERIV
      COMMON/OPEN1/NUNIQ,LOPN(5),MOPN(5,5)
      COMMON/OPEN2/GOCC(10)
      COMMON/PARA1/NSHELL,NT,MXT,NISO
      COMMON/PARA2/NPRIM
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/TAPES/ITAP42
      COMMON/TOLER/DLIM
cges  EQUIVALENCE (CC,IA)
    1 FORMAT(/,2X,' THE DERIVATIVE PROGRAM')
    2 FORMAT(4(A8,2X))
    3 FORMAT(/,3X,4(A8,2X))
    4 FORMAT(5I5)
    5 FORMAT(A8,2X,I5)
    6 FORMAT(F10.5,10I5)
    7 FORMAT(4(A8,2X),I5)
    8 FORMAT(/,2X,' WRONG INPUT '/
     1          2X,' SCFTYP = ',A8,5X,' CITYP = ',A8,5X,' DERTYP = ',A8,
     2        /)
    9 FORMAT(/,2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
   10 FORMAT(/,2X,' REQUIRED MEMORY EXCEEDS MAXCOR'/
     1          2X,' ICMAX = ',I10,5X,' MAXCOR = ',I10/)
   11 FORMAT(/,2X,' WRONG INPUT '/
     1          2X,' THE LABEL "# DERIV ##" WAS NOT FOUND IN THE INPUT,'
     2     ,' DECK.',/)
   12 FORMAT(/,' USING NEW FILE55 IO (INTEGER WORD ADDRESSING)',/)
C
cets030491      call drum
      call psinit('APPEND')
      CALL TSTART(6)
clj   write(6,*)'entering derivth'
cges  CALL NOUNFL
C12-1-88
clj   WRITE(6,12)
C
      DLIM=1.0D-16
      ITAPE3 =  3
      ICHEK  =  3
      INPUT  =  5
      IOUT   =  6
      ITAP30 = 30
      ITAP42 = 42
C3-12-88   MAXCOR=1200000
CMachineSize (For 15Mb)  MAXCOR=850000
CMACHINESIZE (FOR  7MB)  MAXCOR=350000
cges  MAXCOR=350000
C
C  READ IN PARAMETERS FROM TAPE30
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)
      IOPEN=I30(43)
      CALL SREW(ITAP30)
C
      CALL LOCATE(INPUT,'# DERIV ##',IERR)
C
      oldinp=.true.
      IF(IERR .EQ. 0) GOTO 100
      oldinp = .false.
cets030491      WRITE(IOUT,11)
cets030491      WRITE(ICHEK,11)
cets030491      STOP 99100
C
  100 CONTINUE
      IOFF(1)=0
      DO 101 I=1,255
  101 IOFF(I+1)=IOFF(I)+I
C
      WRITE(6,1)
C
C   CLSCF--+
C   HSOP---+
C   OPSNG--+--SCF----SECOND   FOR SCF 2ND DERIVATIVE
C   GRSCF--+
C   TCSCF--+
C
C   CLSCF--+
C   HSOP---+
C   OPSNG--+--CI-----SECOND   FOR SCF-CI 2ND DERIVATIVE
C   GRSCF--+
C   TCSCF--+
C
C   MCSCF--+--MCSCF--SECOND   FOR MCSCF 2ND DERIVATIVE
C
C   MCSCF--+--CI-----SECOND   FOR MCSCF-CI 2ND DERIVATIVE
C
cets030491 add new input here
      if(oldinp) then
        READ(5,2) SCFTYP,CITYP,DERTYP,RAMAN
        WRITE(6,3) SCFTYP,CITYP,DERTYP,RAMAN
        READ(5,4) IPRNT
        NUNIQ=0
        IF(IOPEN.NE.0.AND.SCFTYP.EQ.'GRSCF   ') THEN
          READ(5,5) OPTYPE,NUNIQ
          DO 102 I=1,NUNIQ
            READ(5,6) GOCC(I),LL,(MOPN(I,J),J=1,LL)
            LOPN(I)=LL
  102     CONTINUE
        END IF
        READ(5,7) NOSYM,AFOCK,DFOCK,DERINT
      else
        errcod = frdc('WFN',cityp)
        if(errcod.ne.EOK) then
          write(6,*) ' no WFN input found in deriv'
          call qabort
        endif

        dertyp = 'FIRST   '
        errcod = frdc('DERTYPE',dertyp)

        iraman=0
        errcod = frdboo('RAMAN_INT',iraman)

        if(iraman.eq.1) raman = 'RAMAN   '

        opentp = 'NONE    '
        errcod = frdc('OPENTYPE',opentp)

        iprnt=0
        errcod = frdi('IPRINT',iprnt)

cets the following are added for compatability with old input
c    it is advised that they never be used
        nosym = ' '
        errcod = frdc('NOSYM',nosym)

        afock = ' '
        errcod = frdc('AFOCK',afock)

        dfock = ' '
        errcod = frdc('DFOCK',dfock)

        derint = ' '
        errcod = frdc('DERINT',derint)

cets
        scftyp = 'GRSCF   '
        if(cityp.eq.'SCF     ') then
          if(opentp.eq.'SINGLET ') then
            if(iopen.lt.0) then
              write(6,*) ' file30 says this is TWOCON, but'
              write(6,*) ' OPENTYPE = SINGLET '
              write(6,*) ' '
              write(6,*) ' most peculiar'
              call pabort
              call qabort
            endif
            nuniq=2
            gocc(1)=1.0d0
            gocc(2)=1.0d0
            lopn(1)=1
            lopn(2)=1
            mopn(1,1) = 1
            mopn(2,1) = 2
          else if(opentp.eq.'HIGHSPIN') then
            nuniq=1
            gocc(1)=1.0d0
            if(iopen.lt.0) then
              write(6,*) ' file30 says this is TWOCON, but'
              write(6,*) ' OPENTYPE = HIGHSPIN '
              write(6,*) ' '
              write(6,*) ' most peculiar'
              call pabort
              call qabort
            endif
            errcod = fcnt('SOCC',nopshl)
            errcod = frdia1('SOCC',opshls,nopshl)
            lopn(1)=0
c010692  5001        lopn(1)=lopn(1)+opshls(1) change 1 to i
            do 5001 i=1,nopshl
 5001        lopn(1)=lopn(1)+opshls(i)
            do 5000 i=1,lopn(1)
 5000         mopn(1,i) = i
          else if(opentp.eq.'SPECIAL ') then
            write(6,*) 'cannot do special cases'
            call qabort
          else
            scftyp = 'CLSCF   '
          endif
        endif

        if(opentp.eq.'TWOCON  ') then
           scftyp = 'TCSCF   '
           tcgrfl = 0
           errcod = frdboo('TCSCF_GR',tcgrfl)
           if (tcgrfl.eq.1) then
              scftyp = 'TCSCF-GR'
           endif
        endif

        if(cityp(1:2).eq.'CC') cityp = 'CC      '

        if(cityp.eq.'CC      '.and.opentp.eq.'NONE    ') then
          scftyp = 'CLSCF   '
        endif
c
        write(6,3) scftyp,cityp,dertyp,raman
      endif
c
cets030691 add print control options
c
c have two modes for now, DEBUG, and BRIEF
c       BRIEF  prints only the final gradient (no 2derivs printed)
c       DEBUG  will continue to change...prints everything I can find
c
c       if !BRIEF print all contributions to gradient and 2nd derivs
c
      errcod = prcntr('ON BRIEF')
      errcod = prcntr('OFF VERBOSE DEBUG')
      errcod = prdup()
      
      if(prcntr('IS_ON VERBOSE').eq.1.or.prcntr('IS_ON DEBUG').eq.1)
     +  errcod = prcntr('OFF BRIEF')

cets0030691
      iopen = abs(iopen)
C
      ISTYPE=0
      IF(SCFTYP.EQ.'CLSCF   ') ISTYPE=1
      IF(SCFTYP.EQ.'HSOP    ') ISTYPE=2
      IF(SCFTYP.EQ.'OPSNG   ') ISTYPE=3
      IF(SCFTYP.EQ.'GRSCF   ') ISTYPE=4
      IF(SCFTYP.EQ.'TCSCF-GR') ISTYPE=5
      IF(SCFTYP.EQ.'TCSCF   ') ISTYPE=6
      IF(SCFTYP.EQ.'MCSCF   ') ISTYPE=7
      ICTYPE=0
      IF(CITYP.EQ.'SCF     ') ICTYPE=1
      IF(CITYP.EQ.'CI      ') ICTYPE=2
      IF(CITYP.EQ.'MCSCF   ') ICTYPE=3
      IF(CITYP.EQ.'GVBCI   ') ICTYPE=4
CGES
      IF(CITYP.EQ.'CC      ') ICTYPE=9
C
      NDERIV=0
      IF(DERTYP.EQ.'FIRST   ') NDERIV=1
      IF(DERTYP.EQ.'SECOND  ') NDERIV=2
      IF(DERTYP.EQ.'THIRD   ') NDERIV=3
      IF(DERTYP.EQ.'FOURTH  ') NDERIV=4
      NTYPES=0
      IF(SCFTYP.EQ.'CLSCF   ') NTYPES=1
      IF(SCFTYP.EQ.'HSOP    ') NTYPES=2
      IF(SCFTYP.EQ.'OPSNG   ') NTYPES=3
      IF(SCFTYP.EQ.'GRSCF   ') NTYPES=NUNIQ+1
      IF(SCFTYP.EQ.'TCSCF-GR') NTYPES=3
      IF(SCFTYP.EQ.'TCSCF   ') NTYPES=3
      IF(SCFTYP.EQ.'MCSCF   ') NTYPES=1
      NTYPEP=NTYPES+1
      IF(NDERIV.GE.2.AND.ISTYPE.EQ.2) ISTYPE=4
      IF(NDERIV.GE.2.AND.ISTYPE.EQ.3) ISTYPE=4
      IF(NDERIV.GE.2.AND.ISTYPE.EQ.6) ISTYPE=5
C
      IF(ISTYPE*ICTYPE*NDERIV.EQ.0) GO TO 399
C
      if(oldinp) IRAMAN=0
      IF(RAMAN.EQ.'RAMAN   ') IRAMAN=1
      IDFOCK=0
      IF(NDERIV.GE.2) IDFOCK=1
      IF(ICTYPE.GE.2) IDFOCK=1
      IF(DFOCK.EQ.'FIRST   ') IDFOCK=1
      IF(DFOCK.EQ.'SECOND  ') IDFOCK=2
      IF(DFOCK.EQ.'BOTH    ') IDFOCK=3
      IAFOCK=0
      IF(AFOCK.EQ.'AFOCK   ') IAFOCK=1
      IF(ICTYPE.EQ.2) IAFOCK=1
      IDINT=0
      IF(DERINT.EQ.'FIRST   ') IDINT=1
      IF(DERINT.EQ.'SECOND  ') IDINT=2
      IF(DERINT.EQ.'BOTH    ') IDINT=3
C
C================================================
C===READ IN BASIS SET AND GEOMETRY FROM TAPE30===
C================================================
C
      WRITE(3,21)
clj   WRITE(6,21)
   21 FORMAT(//,2X,' NOW YOU ARE IN GETBAS'/)
      CALL GETBAS(NPERSH)
C
C===========================================
C===FORCE CALCULATION TO NOSYM IF DESIRED===
C===========================================
C
      IF(NOSYM.NE.'NOSYM   ') GO TO 201
      WRITE(3,22)
   22 FORMAT(//,2X,' NOW YOU ARE IN SETNOS'/)
      CALL SETNOS
C
C========================
C===FORM LIMITS ARRAYS===
C========================
C
  201 CONTINUE
      WRITE(3,23)
   23 FORMAT(//,2X,' NOW YOU ARE IN LIMITS'/)
      CALL LIMITS
C
C:::::::::::::::::::::::::::
C:::::NUCLEAR REPULSION:::::
C:::::::::::::::::::::::::::
C
C   CALCULATE DERIVATIVES OF NUCLEAR REPULSION
      WRITE(3,24)
   24 FORMAT(//,2X,' NOW YOU ARE IN NUCTWO'/)
      IC1=1
      IC2=IC1+N3N
      IC3=IC2+N3N
      IC4=IC3+N3N
      IC5=IC4+N3N
      IC6=IC5+N3N
      IC7=IC6+N3N
      IC8=IC7+N3N*N3N
      IC9=IC8+N3N*N3N
      IC10=IC9+N3N*N3N
      IC11=IC10+N3N*N3N
      IC12=IC11+N3N*N3N
      ICMAX=IC12+N3N*N3N
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 499
C.................D1N     D2N.....
      CALL NUCTWO(CC(IC1),CC(IC7))
C
C=========================================================
C===FORM DENSITY AND LAGRANGIAN MATRICES FOR SCF AND CI===
C=========================================================
C
c12-1-88  note change in go to
      IF(NDERIV.GE.2) GO TO 202
      IF(CITYP.EQ.'CI      ') GO TO 202
CGES
      IF(CITYP.EQ.'CC      ') GO TO 2020
C
        WRITE(3,25)
   25   FORMAT(//,2X,' NOW YOU ARE IN SCFMAT'/)
        IC13=IC12+N3N*N3N
        IC14=IC13+NBATRI*NTYPEP
        IC15=IC14+NBATRI
        IC16=IC15+NBATRI
        IC17=IC16+NBASIS
        IC18=IC17+NBASIS
        IC19=IC18+NBFAO*NBASIS
        IC20=IC19+NBFAO*NBASIS
        IC21=IC20+NBASIS*NBASIS
        IC22=IC21+NBFAO*NBFAO
        IC23=IC22+NBFAO*NBFAO
        ICMAX=IC23+NBFAO*NBFAO
        WRITE(3,9) ICMAX,MAXCOR
        IF(ICMAX.GT.MAXCOR) GO TO 499
C.................WAO      EIG      OCC      EAO      UAS      ESO......
      CALL SCFMAT(CC(IC15),CC(IC16),CC(IC17),CC(IC18),CC(IC19),CC(IC20),
C.................EE       U        T........
     1            CC(IC21),CC(IC22),CC(IC23))
C
C   FORM DENSITY MATRICES
        WRITE(3,26)
   26   FORMAT(//,2X,' NOW YOU ARE IN DENS30'/)
        IC13=IC12+N3N*N3N
        IC14=IC13+NBATRI*NTYPEP
        IC15=IC14+NBATRI
        IC16=IC15+NBATRI
        IC17=IC16+NBASIS
        IC18=IC17+NBASIS
        ICMAX=IC18+NBFAO*NBASIS
        WRITE(3,9) ICMAX,MAXCOR
        IF(ICMAX.GT.MAXCOR) GO TO 499
C...................DENS     DENT     EIG      OCC      EAO......
        CALL DENS30(CC(IC13),CC(IC14),CC(IC16),CC(IC17),CC(IC18))
        GO TO 203
C
C   READ IN PARAMETERS FROM THE MASTER FILE
  202 CONTINUE
      WRITE(3,27)
   27 FORMAT(//,2X,' NOW YOU ARE IN READMF'/)
      CALL READMF
C
C   FORM DENSITY MATRICES
      WRITE(3,28)
   28 FORMAT(//,2X,' NOW YOU ARE IN DENS40'/)
      IC13=IC12+N3N*N3N
      IC14=IC13+NBATRI*NTYPEP
      IC15=IC14+NBATRI
      IC16=IC15+NBATRI
      IC17=IC16+NBATRI
      IC18=IC17+NBASIS
      IC19=IC18+NBASIS
      IC20=IC19+NBFAO*NBASIS
      IC21=IC20+NBFAO*NBASIS
      ICMAX=IC21+NBASIS*NBASIS
C12-1 ICMAX=IC19+NBFAO*NBASIS
      WRITE(3,9) ICMAX,MAXCOR
C12-1  ADD DEN1 AFER WAO
      IF(ICMAX.GT.MAXCOR) GO TO 499
C.................DENS     DENT     WAO      DEN1     EIG      OCC......
      CALL DENS40(CC(IC13),CC(IC14),CC(IC15),CC(IC16),CC(IC17),CC(IC18),
C.................EAO      U        T........
     1            CC(IC19),CC(IC20),CC(IC21),iopen)
C
C:::::::::::::::::::::::::::::::::
C:::::SCF OR CI WAVE FUNCTION:::::
C:::::::::::::::::::::::::::::::::
C
C   CALCULATE SCF DERIVATIVES
  203 CONTINUE
 2020 CONTINUE
      WRITE(3,29)
   29 FORMAT(//,2X,' NOW YOU ARE IN SCFDER'/)
C#################################################
C###NOTICE THE FOLLOWING RESTRICTIONS#############
      NTREAD=NTYPEP
      IF(ISTYPE.EQ.5.OR.ISTYPE.EQ.6) NTREAD=7
      IF(IAFOCK.EQ.0) NTREAD=NTREAD-1
      IF(IDFOCK.LE.1) NATRI=1
C###END OF NOTICE#################################
C#################################################
C12-1-88  PUT IN SPACE FOR CI ONEPDM
      IC13=IC12+N3N*N3N
      IC14=IC13+NBATRI*NTYPEP
      IC15=IC14+NBATRI
      IC16=IC15+NBATRI
      IC17=IC16+NBATRI
      IA17=IC17+IC17-1
      MAXLFT=MAXCOR-IC17
      ICMAX=IC17+NBATRI*N3N+NBATRI*3+NBATRI*NATRI+NBATRI*6+NBATRI*2
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 499
CGES
      IF(CITYP.EQ.'CC      ') THEN
      ITAP79=79
      CALL RFILE(ITAP79)
      I79=1+INTOWP(NTRI+NTRI)
      CALL WREADW(ITAP79,ECC,INTOWP(1),I79,I79)
      ESCF=ECC
      WRITE(6,663)ESCF
  663 FORMAT(1X,//,' TOTAL ENERGY FROM FILE79  =',F20.12)
      I79=I79+INTOWP(NTRI+NBASIS*NBASIS)
      CALL WREADW(ITAP79,CC(IC16),INTOWP(NBATRI),I79,I79)
C
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        WRITE(6,*)' '
        WRITE(6,*)'AO ONE-PDM in fentry from 79 '
        CALL PRINT(CC(IC16),NBATRI,NBASIS,6)
      endif
c
      CALL WREADW(ITAP79,CC(IC15),INTOWP(NBATRI),I79,I79)
C
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        WRITE(6,*)' '
        WRITE(6,*)'AO SYMMETRIZED LAGRANGIAN in fentry from 79'
        CALL PRINT(CC(IC15),NBATRI,NBASIS,6)
      endif
C
      CALL RCLOSE(ITAP79,3)
      ENDIF
C
C.................D1N     D1O     D1F     D1W     D1T     D1TOT...
      CALL SCFDER(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),
C.................D2N     D2O     D2F     D2W      D2T      D2TOT....
     1            CC(IC7),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),
C.................DENS     DENT     WAO      DEN1.....
     2            CC(IC13),CC(IC14),CC(IC15),CC(IC16),
C.................IA       CC..............
     3            IA(IC17),CC(IC17),MAXLFT,NPERSH)
cges this is a bug !!!!!!!!!!!!  ia(ia17) is equivalent to cc(ic17)
      GO TO 210
C
C:::::::::::::::::::::::::::::
C:::::MCSCF WAVE FUNCTION:::::
C:::::::::::::::::::::::::::::
C
C   CALCULATE MCSCF DERIVATIVES
  206 CONTINUE
      IF(CITYP.NE.'MCSCF   ') GO TO 207
      WRITE(3,31)
   31 FORMAT(//,2X,' NOW YOU ARE IN MCSCF'/)
      IC13=IC12+N3N*N3N
      IC14=IC13+NBATRI
      IC15=IC14+NBATRI
      IA15=IC15+IC15-1
      MAXLFT=MAXCOR-IC15
      ICMAX=IC16+NBATRI*N3N+NBATRI*3+NBATRI*NATRI+NBATRI*6+NBATRI
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 499
C................D1N     D1O     D1F     D1W     D1T     D1TOT...
      CALL MCSCF(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),
C................D2N     D2O     D2F     D2W      D2T      D2TOT....
     1           CC(IC7),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),
C................DENT     WAO      IA       CC..............
     2           CC(IC13),CC(IC14),IA(IC15),CC(IC15),MAXLFT)
      GO TO 210
C
C::::::::::::::::::::::::::::::::
C:::::MCSCF-CI WAVE FUNCTION:::::
C::::::::::::::::::::::::::::::::
C
C   CALCULATE MCSCF-CI DERIVATIVES
  207 CONTINUE
      IF(CITYP.NE.'CI      ') GO TO 210
      WRITE(3,32)
   32 FORMAT(//,2X,' NOW YOU ARE IN MCSCI'/)
      IC13=IC12+N3N*N3N
      IC14=IC13+NBATRI*NTYPEP
      IC15=IC14+NBATRI
      IA15=IC15+IC15-1
      MAXLFT=MAXCOR-IC15
      ICMAX=IC16+NBATRI*N3N+NBATRI*3+NBATRI*NATRI+NBATRI*6+NBATRI
      WRITE(3,9) ICMAX,MAXCOR
      IF(ICMAX.GT.MAXCOR) GO TO 499
C................D1N     D1O     D1F     D1W     D1T     D1TOT...
      CALL MCSCI(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),
C................D2N     D2O     D2F     D2S      D2T      D2TOT....
     1           CC(IC7),CC(IC8),CC(IC9),CC(IC10),CC(IC11),CC(IC12),
C................DENT     WAO      IA       CC..............
     2           CC(IC13),CC(IC14),IA(IC15),CC(IC15),MAXLFT)
C
  210 CONTINUE
      GO TO 500
  399 CONTINUE
      WRITE(3,8) SCFTYP,CITYP,DERTYP
      WRITE(6,8) SCFTYP,CITYP,DERTYP
      GO TO 500
  499 CONTINUE
      WRITE(3,10) ICMAX,MAXCOR
      WRITE(6,10) ICMAX,MAXCOR
  500 CONTINUE
      CALL TSTOP(6)
      STOP
      END
