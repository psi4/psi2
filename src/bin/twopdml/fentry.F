cets The following is called by a C program which handles dynamic
cets allocation of memory.
      subroutine fentry(z,id,maxcor)
cets  PROGRAM TWOPDML
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
c   Last updated on July 09, 1998                                     *
c            by Dr. Yukio Yamaguchi                                   *
c   Modification for an integer array (IA)                            *
c      `Maxi' is defined as a parameter (currently 200000)            *
c       This `Maxi' value can be increased if required.               *
c**********************************************************************
c       by: Edward Seidl (cets030291)                                 *
c     date: March 3, 1991                                             *
c   reason: add new input, and replace drum with psinit               *
c                                                                     *
C*  BY:  RICHARD REMINGTON                         SEARCH:  C4-29-88  *
C*  DATE:  APRIL  29,  1988                                           *
C*  REASON: DECREASE CORE: ARRAY Z(720000)  TO  Z(680000)             *
C**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c3-24-88  *
C*  DATE:  MARCH  24,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C*          AND COMMENT OUT UNREACHABLE CODE                          *
C**********************************************************************
C      BY:  TRASH                                                     *
C    DATE:  JUNE 29, 1987                                             *
C  REASON:  FIX SMALL BUG----MAKE MAXI=IA(90000)  COMPATIBLE WITH DRT *
C**********************************************************************
C      BY:  WESLEY ALLEN                                              *
C    DATE:  NOVEMBER 5, 1986                                          *
C  REASON:  ACTIVATE THE ROOT SELECTION OPTION.                       *
C           (CHANGES IN TWOPDMLI)                                     *
C**********************************************************************
C      BY:  TJL                                                       *
C    DATE:  AUGUST 4, 1985                                            *
C  REASON:  CORRECT CALLS TO NCDLBL AND DCDLBL                        *
C**********************************************************************
C      BY:  RICHARD REMINGTON                                         *
C    DATE:  FEBRUARY 14, 1985                                         *
C  REASON:  CLEANUP OUTPUT AND STANDARDIZE UNIT 6 TO OUTPUT AND       *
C     UNIT 3 TO CHECK                                                 *
C  REASON:  CHANGE UNIT NUMBERS !! 2=>52, 82=>56                      *
C  REASON:  CHANGE UNIT NUMBERS !! 8=>58, 20=>53                      *
C  REASON:  CHANGE UNIT NUMBERS !!12=>54,                             *
C**********************************************************************
C    CONVERTED TO THE IBM BY DJF AND TJL  OCTOBER 1984                *
C**********************************************************************
C  *******  LAST UPDATED ON DECEMBER 13 1982  *******                 *
C  * MODIFIED TO SPLIT I-BLOCKS ON JAN. 20 1983 DJF *                 *
C**********************************************************************
C       PROGRAM OPTIONS
C
C       1.  IGUESS=ROOT OF CI TO BE USED IN CREATION OF DENSITY MATRIX
C                  DEFAULT=1.
C       2.  IPRINT=PRINT OPTION (OPERATIONAL 5-1-86 RBR)
C                  DEFAULT=0.
C
C***********************************************************************
C     FILE ASSIGNMENTS
C
C         5=INPUT FILE, SINGLE CARD I5 FORMAT
C         6=OUTPUT FILE
C        53=TWO-PARTICLE DENSITY MATRIX STORED LIKE THE INTEGRALS
C        54=FINAL VECTOR(S),112 WORD LABEL + NUMBER OF WALKS LONG
C        58=DISTINCT ROW TABLE, A FEW 10,000 LONG AT MOST.
C
C***********************************************************************
      IMPLICIT INTEGER (A-Z)
cets030291
#include <error.h>
      integer frdi, errcod
      INTEGER NUMINT
      LOGICAL PAGEIN
c3-24-88  COMMON Z(1200000)
C4-29-88  COMMON Z(720000)
cets  COMMON Z(680000)
cets the following did not work
cets  common z(maxcor)
      real*8 z(maxcor)
      integer id(1)
c=============================================
cyy 07-08-1998
c   Modification for a larger integer array
c=============================================
      parameter (maxi=200000)
      dimension ia(maxi)
c=============================================
cets  REAL*8 Z,X,REP,FZCORE,EGUESS,ECI,CNVERG,ACRCY,SQCDIF
      REAL*8   X,REP,FZCORE,EGUESS,ECI,CNVERG,ACRCY,SQCDIF
      REAL*8 CZERO,DRTVER,CIVER,INPUT,IVER
      REAL*8 TIME,TOTAL
      COMMON /IENT/ IDEN,IEXT,LOPIN,IALL,LOPEX,INTOD
C
cyy 07-08-1998  COMMON /DRT/  IA(90000)
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP03,ITAP04
     *,             ITAPE3,ITAP05,ITAP06,ITAP56,IT53
      COMMON /DIMS/ NBF,NSYM,NORBS,NROWS,NROWS4,NWKS,NWKS2,NLEVS
     *,             NROWOC,NROW4O,NWKSOC,NLEVOC,NORBOC,LEVFRM
     *,             NWKSMX,NLWKMX,NUWKMX,BMAX,NROOTS,ORBFRM
      COMMON /INTS/ NMAX,NMAX2,NGROUP,NBLKOC,NUMIJ,SYMORB,INTSRT
      COMMON /OPS/  IGUESS,IRSTRT,IROOTI,IROOTF,IPRINT
      COMMON /ALL/  X(3),ARR,IX(6),LVFRM1,NLWKI,NLWKJ,IMAX,IMIN
      COMMON /SYMM/ JSM,JFSYM,IFSYM,MAXSYM(8),MINSYM(8),ISMOFF(8)
     #,             NUMSYM(8)
CTJL  COMMON /DIAG/ REP,FZCORE,EGUESS,ECI,SQCDIF,CNVERG
CTJL *,             CZERO,MXITER,ITER,ICNVG,NROOT,REFWLK
      COMMON /DIAG/ REP,FZCORE,EGUESS,ECI,REFWLK,MXITER,CNVERG,ICNVG
     *,             ITER,SQCDIF,CZERO,NROOT
      COMMON /PAGE/ IWORD3,IWORD4,PAGEIN,NOFFI,NOFFJ,NPASS,NWORDL,IWORD6
      COMMON /LBLS/ CIVER,DRTVER,LBLINT(26),LBLDRT(26)
C*
      COMMON /XX4X/   NIJVIR,NREFS
C
      DIMENSION ILBL(26),ILBM(26),ILBN(26),III(3),IDENT(26)
      DIMENSION IFCB8(16),IFCB12(16)
C
cets030291      call drum
      call psinit('APPEND')
c     CALL TSTART(ITAPE3)
      CALL TSTART(ITAPE6)
C
C
C     MAXI IS NOT AN INTEGER EQUIVALENT OF MAXR
cets080593      MAXI= 90000
c===============================================
cyy-07-08-1998  following statement is an error
cyy     maxi = intowp(maxcor)
c===============================================
C
C
c3-24-88  MAXR=1200000
C4-29-88  MAXR= 720000
cets  MAXR= 680000
      maxr = maxcor
C
C
C     WRITE(ITAPE3,'(" PLEASE ENTER MAXR")')
C     READ (3,*) MAXR
C     IF (MAXR.EQ.0) MAXR=160000
      PDMVER=4.830208
      WRITE(ITAPE6,1)
    1 FORMAT(' ')
    4 FORMAT(///,13X,'     >>>>> TWOPDM CALCULATION FINISHED <<<<< ',/)
    5 FORMAT(/)
      WRITE(ITAPE6,10)
   10 FORMAT(' ',12X,'SHAPE DRIVEN GRAPHICAL UNITARY GROUP TPDM PROGRAM'
     1)
   11 FORMAT(' ',12X,'************ MATRIX REFORMULATION *************')
   12 FORMAT(' ',25X,'WRITTEN BY PAUL SAXE')
      WRITE(ITAPE6,13)
   13 FORMAT(' ',15X,'BASED ON ORIGINAL SDGUGA PROGRAM WRITTEN BY')
      WRITE(ITAPE6,14)
   14 FORMAT(' ',16X,'PAUL SAXE, DOUGLAS J. FOX AND N. C. HANDY',//
     , 11X,'MODIFIED TO CONSTRUCT TWO PARTICLE DENSITY MATRIX BY'
     , /22X,'GEORGE FITZGERALD AND DOUG FOX',/
     , /22X,'Modified for a larger basis set by'
     , /31X,'Yukio Yamaguchi'
     , /25x,' (Version of July 09, 1998)',///)
      WRITE(ITAPE6,15) MAXI,MAXR
   15 FORMAT(' ','MAXIMUM DIMENSION OF INTEGER COMMON = ',I10,/
     *       ' ','MAXIMUM DIMENSION OF REAL COMMON    = ',I10,/)
C   READ IN OPTIONS FOR THIS RUN
C
      IGUESS = 1
      IPRINT = 0
      CALL LOCATE(ITAPE5,'# TWOPDM #',IERROR)
cets030291
      if(ierror.ne.0) then
c           WRITE(ITAPE6,16)
c           GO TO 17
c   16 FORMAT(1X,'WARNING! NO TWOPDM INPUT FOUND ',
c     +' **** DEFAULT OPTIONS BEING USED **** ')
c
        errcod = frdi('ROOT1',iguess)
        errcod = frdi('IPRINT',iprint)
      else
        READ(ITAPE5,18) IGUESS,IPRINT
      END IF
      IF(IGUESS.EQ.0) IGUESS = 1
   18 FORMAT(10I5)
   17 CONTINUE
C
      WRITE(ITAPE6,19) IGUESS,IPRINT
   19 FORMAT(' ','OPTIONS',11X,'IGUESS   IPRINT',/,15X,6I8)
      CALL RFILE(ITAP54)
      CALL SREW(ITAP54)
      CALL DCDLBL(ITAP54,LBLINT,LBLDRT,REFWLK,IX,IX,EGUESS,ECI,CZERO)
CTJL  CALL DCDLBL(ITAP56,LBLINT,LBLDRT,IX,IX,IX,EGUESS,ECI,CZERO)
      CALL SREW(ITAP54)
CTJL  CALL SREW(ITAP56)
      IRE = REFWLK
      CALL NCDLBL(ITAP54,LBLINT,LBLDRT,IRE,IX,IGUESS,EGUESS,ECI,CZERO)
CTJL  CALL NCDLBL(ITAP56,LBLINT,LBLDRT,IX,IX,IGUESS,EGUESS,ECI,CZERO)
      WRITE(ITAPE6,4321) EGUESS,ECI
 4321 FORMAT(/,' EGUESS = ',F20.12,/,' ECI    = ',F20.12,/)
       WRITE(ITAPE3,20)
   20 FORMAT(' ','READ IN DISTINCT ROW TABLE')
      CALL RFILE(ITAP58)
CIBM  CALL SREW(ITAP58)
CH    READ (ITAP58) DRTVER
C*
CH    IF (DRTVER.LT.3.0) STOP
CH    READ (ITAP58) LBLDRT
CH    READ (ITAP58) NBF,NSYM,NORBS,NROWS,NROWS4
CH   #,              NLEVS,NREFS,NROWOC,NROW4O,NWKS,NWKSOC,NLEVOC
CH   #,              ORBFRM,SYMORB,NUMIJ,NGROUP,INTNUM,NMAX,BMAX
CH   #,              NIJVIR
      CALL WREADW(ITAP58,IVER,INTOWP(1),1,END)
      CALL WREADW(ITAP58,LBLDRT,26,END,END)
      CALL WREADW(ITAP58,NBF,1,END,END)
      CALL WREADW(ITAP58,NSYM,1,END,END)
      CALL WREADW(ITAP58,NORBS,1,END,END)
      CALL WREADW(ITAP58,NROWS,1,END,END)
      CALL WREADW(ITAP58,NROWS4,1,END,END)
      CALL WREADW(ITAP58,NLEVS,1,END,END)
      CALL WREADW(ITAP58,NREFS,1,END,END)
      CALL WREADW(ITAP58,NROWOC,1,END,END)
      CALL WREADW(ITAP58,NROW4O,1,END,END)
      CALL WREADW(ITAP58,NWKS,1,END,END)
      CALL WREADW(ITAP58,NWKSOC,1,END,END)
      CALL WREADW(ITAP58,NLEVOC,1,END,END)
      CALL WREADW(ITAP58,ORBFRM,1,END,END)
      CALL WREADW(ITAP58,SYMORB,1,END,END)
      CALL WREADW(ITAP58,NUMIJ,1,END,END)
      CALL WREADW(ITAP58,NGROUP,1,END,END)
      CALL WREADW(ITAP58,NUMINT,1,END,END)
      CALL WREADW(ITAP58,NMAX,1,END,END)
      CALL WREADW(ITAP58,BMAX,1,END,END)
      CALL WREADW(ITAP58,NIJVIR,1,END,END)
      WRITE(ITAPE3,*) ' NIJVIR = ',NIJVIR
      END=END+1
      NMAX2=NMAX+NMAX
      NWKS2=NWKS+NWKS
      LVFRM1=ORBFRM
      LEVFRM=LVFRM1+1
      NORBOC=NLEVOC-1
 
      WRITE(ITAPE6,31)LBLDRT
   31 FORMAT(' LABEL FROM DRT :   ',26A3,/)
      WRITE(ITAPE6,32)NWKS
   32 FORMAT(' NUMBER OF CONFIGURATIONS    = ',I10)
      WRITE(ITAPE6,33)NBF
   33 FORMAT(' NUMBER OF BASIS FUNCTIONS   = ',I10)
      WRITE(ITAPE6,34) NORBS
   34 FORMAT(' NUMBER OF ALLOWED ORBITALS  = ',I10)
      WRITE(ITAPE6,35)LVFRM1
   35 FORMAT(' FERMI LEVEL IS ',I5,/)
      WRITE(ITAPE3,*) ' NMAX = ',NMAX
C
C     CALCULATE POINTERS INTO BLANK COMMON
      ARC=1
      IJXX=ARC+NROWS4
      KLXX=IJXX+NUMIJ
      NKLXX=KLXX+NIJVIR
      IJWW=NKLXX+NSYM*ORBFRM
      KLWW=IJWW+NUMIJ
      NKLWW=KLWW+NIJVIR
      IST4X=NKLWW+NSYM*ORBFRM
      LEN4X=IST4X+NORBS
      WEIGHT=LEN4X+NORBS
      WTW=WEIGHT+NROWS4
      WTX=WTW+NSYM*LVFRM1
      WTY=WTX+NSYM*LVFRM1
      WAB=WTY+LVFRM1
      NLWKS=WAB+LVFRM1
      AVAL=NLWKS+NROWS
      BVAL=AVAL+NROWS
      ROWSYM=BVAL+NROWS
      IJADD=ROWSYM+NROWS
      IJGRP=IJADD+NUMIJ
      KADD=IJGRP+NUMIJ
      LADD=KADD+SYMORB
      JMXNXT=LADD+SYMORB
      JMNNXT=JMXNXT+NORBS
      ININT=JMNNXT+NORBS
      INEXT=ININT+NORBS
      NINGRP=INEXT+NORBS
      IMNGRP=NINGRP+NGROUP
      IMXGRP=IMNGRP+NGROUP
      JMNGRP=IMXGRP+NGROUP
      JMXGRP=JMNGRP+NGROUP
      ORBSYM=JMXGRP+NGROUP
      ISEGM=ORBSYM+NORBS
      JSEGM=ISEGM+NLEVS
      IMAIN=JSEGM+NLEVS
      ISUB=IMAIN+NLEVS
      IUWKMN=ISUB+NLEVS
      IUWKSB=IUWKMN+NLEVS
      ITRACK=IUWKSB+NLEVS
      LMIN=ITRACK+NLEVS
      HDWGT=LMIN+NLEVS
      HDSEG=HDWGT+NLEVS
      HDPNT=HDSEG+NLEVS
      ACOEF=HDPNT+NLEVS
cets added on advice of the kind people at ardent
      if (mod(acoef,2).eq.0) then
         write(0,*) 'TWOPDM: acoef = ',acoef,'must be odd'
         acoef = acoef + 1
      endif
      BCOEF=ACOEF+NLEVS*2
      LENI=BCOEF+NLEVS*2
      WRITE(ITAPE6,41) LENI
   41 FORMAT(1X,'LENGTH OF INTEGER COMMON    = ',I10,/)
      IF(LENI.GT.MAXI) GO TO 930
      WRITE(ITAPE6,42)
   42 FORMAT(' READ DRT FILE')
      CALL GETDRT(IA(ISEGM),IA(AVAL),IA(BVAL)
     #,        IA(ROWSYM),IA(ARC),IA(NLWKS)
     #,        IA(WEIGHT),IA(WAB),IA(WTW),IA(WTX),IA(WTY)
     #,                   IA(KADD),IA(LADD),IA(IJADD),IA(IJGRP)
     #,        IA(ININT),IA(INEXT),IA(JMNNXT),IA(JMXNXT),IA(NINGRP)
     #,        IA(ORBSYM),IA(IMNGRP),IA(IMXGRP),IA(JMNGRP),IA(JMXGRP)
     #,        IA(IJXX),IA(KLXX),IA(NKLXX),IA(IJWW),IA(KLWW),IA(NKLWW)
     #,        END)
      NWKSMX=NWKS
      LN=0
      DO 43 I=1,NSYM
      IF (MINSYM(I).GE.1000) NUMSYM(I)=0
      IF (MINSYM(I).LT.1000) NUMSYM(I)=MAXSYM(I)-MINSYM(I)+1
      IF (NUMSYM(I).GT.LN) LN=NUMSYM(I)
   43 CONTINUE
      LNSQ=LN**2
      INTS=1
      CI=INTS+NMAX
      INTSRT=1
      WRITE(ITAPE3,*) ' INTSRT = ',INTSRT
      IF(INTSRT/2.GT.NMAX) CI=INTS+INTSRT/2
      SI=CI+LNSQ
      CJ=SI+LNSQ
      SJ=CJ+LNSQ
      H1=SJ+LNSQ
      H2=H1+LNSQ
      VECTOR=H2+LNSQ
      LENR=VECTOR+NWKS
      WRITE(ITAPE3,*) ' LENR,VECTOR = ',LENR,VECTOR
C ALLOCATE CORE
      IF (LENR.GT.MAXR) THEN
      NASK=190
      ELSE
      NASK=INT(LENR/1000)+1
      END IF
      IF(LENR.LT.MAXR) GO TO 50
      WRITE(ITAPE6,45)
   45 FORMAT(1H ,'INSUFFICIENT CORE TO HOLD ALL VECTORS. PAGING WILL'
     *,' BE TRIED,')
      NWKSMX=(MAXR-VECTOR)
      IF(2*NLWKMX.GT.NWKSMX) GO TO 940
      WRITE(ITAPE6,46) NWKSMX
  46  FORMAT(1X,I8,' ELEMENTS OF THE VECTOR CAN BE HELD AT ONCE.')
      LENR=VECTOR+NWKSMX
      WRITE (ITAPE6,49) LENR
  49  FORMAT (1X,' AMOUNT OF REAL CORE USED:',I8)
  50  CONTINUE
      CALL INITPG(Z(VECTOR))
      CALL INITLP(Z(INTS),Z(VECTOR),Z(LENR),IA(IJADD),IA(KADD)
     *,           IA(LADD),IA(WTW),IA(WTX),IA(WTY),IA(WAB),IA(ORBSYM))
      CALL INITDG(Z(INTS),Z(VECTOR),IA(IJADD),IA(KADD),IA(LADD)
     *,           IA(WTW),IA(WTX),IA(WTY),IA(WAB),IA(ORBSYM))
      CALL INITEX(Z(INTS),Z(VECTOR),IA(IJADD),IA(KADD)
     *,           IA(LADD),IA(WTW),IA(WTX),IA(WTY),IA(WAB),IA(ORBSYM))
      CALL INITXX(Z(INTS),Z(VECTOR),IA(IJADD),IA(KADD)
     *,           IA(LADD),IA(WTW),IA(WTX),IA(WTY),IA(WAB),IA(ORBSYM))
      CALL INITX2(Z(INTS),Z(VECTOR),Z(LENR),IA(IJADD),IA(KADD)
     *,           IA(LADD),IA(WTW),IA(WTX),IA(WTY),IA(WAB),IA(ORBSYM))
C     CALL INITST(Z(INTS),Z(VECTOR),IA(IJADD),IA(KADD)
C    *,           IA(LADD),IA(WTW),IA(WTX),IA(WTY),IA(WAB),IA(ORBSYM)
C    *,           Z(CI),Z(SI),Z(CJ),Z(SJ),Z(H1),Z(H2)
C    #,           IA(IJXX),IA(NKLXX),IA(IJWW),IA(NKLWW))
C     CALL INIMAT(Z(INTS),IA(IJADD),IA(KADD),IA(LADD))
      CALL INIT4X(Z(INTS),Z(VECTOR),Z(LENR),IA(IJADD),IA(KADD)
     *,           IA(LADD),IA(WTW),IA(WTX),IA(WTY),IA(WAB),IA(ORBSYM))
      CALL LOOPYD(IA(BVAL),IA(NLWKS),IA(WEIGHT),ia(aval),ia(ijadd)
     *,           IA(IJGRP),IA(KADD),IA(LADD),ia(inext),ia(arc)
     *,           IA(JMNNXT),IA(JMXNXT),IA(ROWSYM),ia(orbsym),ia(isegm)
     *,           IA(JSEGM),IA(IMAIN),IA(ISUB),IA(IUWKMN),IA(IUWKSB)
     *,           IA(ITRACK),IA(LMIN),IA(ACOEF),IA(BCOEF)
     *,           IA(HDWGT),IA(HDSEG),IA(HDPNT),IA(IMXGRP),IA(IMNGRP)
     #,           IA(JMXGRP),IA(JMNGRP),Z(INTS),z(vector),z(lenr)
     #,           ia(wtw),ia(wtx),ia(wty),ia(wab))
      WRITE(ITAPE6,1001) IDEN
 1001 FORMAT(/' ENTRIES INTO DIAGONAL   = ',I12)
      WRITE(ITAPE6,1002) INTOD
 1002 FORMAT(/' INTERNAL OFF-DIAGONAL   = ',I12)
      WRITE(ITAPE6,1003) IEXT
 1003 FORMAT(/' ENTRIES INTO EXTERNAL   = ',I12)
      WRITE(ITAPE6,1004) LOPEX
 1004 FORMAT(/' ENTRIES INTO LOOPEX     = ',I12)
      WRITE(ITAPE6,1005) LOPIN
 1005 FORMAT(/' ENTRIES INTO LOOPIN     = ',I12)
      WRITE(ITAPE6,1006) IALL
 1006 FORMAT(/' ENTRIES INTO ALLEXT     = ',I12)
C
      WRITE(ITAPE6,4)
C
      JCODE = 3
      CALL RCLOSE(IT53,JCODE)
      CALL RCLOSE(ITAP54,JCODE)
      CALL RCLOSE(ITAP58,JCODE)
C
c     CALL TSTOP(ITAPE3)
      CALL TSTOP(ITAPE6)
C
      stop
  900 WRITE(ITAPE6,905)NWKS,NWKSB
  905 FORMAT(' NUMBER OF CONFIGURATIONS DOES NOT MATCH ',2I10)
      CALL DUMP
  910 WRITE(ITAPE6,915)NSYM,NSYM2
  915 FORMAT(' NUMBER OF SYMMETRY TYPES DOES NOT MATCH ',2I10)
      CALL DUMP
  930 WRITE(ITAPE6,935)LENI
  935 FORMAT(' INSUFFICIENT SPACE IN INTEGER BLANK COMMON NEED',I10,
     *' WORDS.')
      CALL DUMP
  940 WRITE(ITAPE6,945)4*NLWKMX+VECTOR
  945 FORMAT(' INSUFFICIENT SPACE IN BLANK COMMON NEED',I10,' WORDS.')
      CALL DUMP
      stop
      END
