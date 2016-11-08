cges  PROGRAM CCSR
      subroutine fentry(cc,ic,maxcor)
C
C  >>>  LAST UPDATE :  april 29, 1990  <<<
C
C***********************************************************************
C  NOTICE OF PROGRAM MODIFICATIONS
c     by:  Edward Seidl                  search: cets070991
c   date:  July 9,1991
c reason:  hardwire to use <= 8MB of core.  option LARGE_CORE added
c          will use amount of core specified by MEMORY if LARGE_CORE is
c          true
c
C      BY:  GUSTAVO E. SCUSERIA (unless otherwise noted)
c
c    date:  Thu Feb 28 11:56:43 EST 1991
c      by:  Curtis L Janssen
c  reason:  added new style input
c
c    date:  april 29, 1990
c  reason:  add # CCSR #### with code71 (I5) .ne.0 to avoid
c           creating file71 in soppa calculations.
C
C    DATE:  MAY   10, 1989
C  REASON:  CHANGE NAME TO CCSR.
C
C    DATE:  MARCH 27, 1989
C  REASON:  CHANGE READ30 ROUTINE TO NEW VERSION IN NCC9.
C
C    DATE:  NOVEMBER 16, 1988
C  REASON:  WRITE OUT ONE-E MO INTS (CC LABELS) TO FILE 79.
C
C    DATE:  NOVEMBER 1, 1988
C  REASON:  IF(INTBUF.LE.0) STOP!
C           READS BUFFER SIZE FROM 78, COMPATIBLE W/ CCTRANSL ONLY.
C
C    DATE:  SEPTEMBER 16, 1988
C  REASON:  FORCE NTAO=NBF FOR NEW GCINTS CODE.
C***********************************************************************
      IMPLICIT INTEGER (A-Z)
      REAL*8 CC(maxcor)
      INTEGER IC(*)
cges  EQUIVALENCE (CC,IC)
      REAL*8 ENUC,ESCF,ECORR,ETOT
      CHARACTER* 4 LABEL(20),OPTION,PRTS,RSTR,CHAR,SORT
      INTEGER NLAMDA(8),NCL(8),ITYP(8),FLOV(8,4),
     .        NOSYM(8),NVSYM(8),NFZO(8),NFZV(8),NDOC(8),NUOC(8),
     .        NT3(8),NOP(8)
c
#include <error.h>
      integer errcod, frdboo, icore
      integer prcntr, prdup
C
c     CALL NOUNFL
C
c     MAXCOR= 500000
C
      call psinit('APPEND')
      CALL TSTART(6)
clj Set up print control
      errcod = prcntr('ON BRIEF')
      errcod = prdup()
      if (prcntr('IS_ON DEBUG').ne.0) then
        errcod = prcntr('OFF BRIEF')
        endif
C
      WRITE(6,6000)
clj save paper
c6000 FORMAT(//,47('-'),/,47('-'),/,
c    .' CCSR  : RESORTED B INTS FOR NEW CC CODES        ',/,
c    .'     WRITTEN BY GUSTAVO SCUSERIA  ',/,
c    .'     LAST UPDATE : MAY 10, 1989   ',/,47('-'),/,47('-'),//)
 6000 FORMAT(/,
     .' CCSR  : RESORTED B INTS FOR NEW CC CODES        ',/,
     .'     WRITTEN BY GUSTAVO SCUSERIA  ',/,
     .'     LAST UPDATE : MAY 10, 1989   ',/)
C
      INPUT=5
      JOUT=6
      ITAP30=30
      ITAP60=60
      ITAP61=61
      ITAP62=62
      ITAP63=63
      ITAP64=64
      ITAP65=65
      ITAP66=66
      ITAP71=71
      ITAP78=78
C
      CALL RFILE(ITAP60)
      CALL RFILE(ITAP61)
      CALL RFILE(ITAP62)
      CALL RFILE(ITAP63)
      CALL RFILE(ITAP64)
      CALL RFILE(ITAP65)
      CALL RFILE(ITAP66)
      CALL RFILE(ITAP78)
C
C >>> GET BUFFER SIZE FROM FILE78
C
C     INTBUF=SEC2I(100)/INTOWP(1)
      CALL SREAD(ITAP78,INTBUF,1)
      WRITE (6,*) ' INTBUF=',INTBUF
      IF(INTBUF.LE.0)THEN
         WRITE(6,*)' PROBLEMS WITH INTBUF. PROGRAM QUITS ...'
         STOP
      ENDIF
c
      call locate(5,'# CCSR ####',ierr)
      if(ierr.ne.0) then
        errcod = frdboo('NO71',code71)
        if (errcod.eq.EType) then
          write(6,*) 'ERROR: NO71 must be boolean'
     &               //' (TRUE, FALSE, YES, NO, 1, or 0)'
          call pabort
          call qabort
          endif
        if (errcod.ne.EOK) code71=0
      else
        write(6,*) 'WARNING: using old style ''# CCSR ####'' input'
        read(5,621)code71
  621   format(16i5)
      endif
      if(code71.ne.0)then
      write(6,*)
      write(6,*)'code71 =',code71
      write(6,*)'code71.ne.0 ===>  file71 is NOT created'
      write(6,*)
      endif
c
cets070991
c    test to see if using more than 8mb of core
c    if MEMORY is > 8mb, then reduce it unless LARGE_CORE is true
c
#ifdef AIX
      icore=0
      errcod = frdboo('LARGE_CORE',icore)
      if(maxcor.gt.1000000.and.icore.eq.0) then
        write(6,*) ' '
        write(6,*) 'You have requested more than 8 MB of memory'
        write(6,*) 'which will cause poor performance on the'
        write(6,*) 'rs6000s, thus I will only use 8 MB of core.'
        write(6,*) 'If this is not enough, may I suggest using'
        write(6,*) 'the LARGE_CORE option (see psi manual)'
        write(6,*) ' '
        maxcor=1000000
      endif
#endif
cets070991
C
C >>> ZERO OUT ALL THE CORE
C
      CALL ZERO(CC,MAXCOR)
C
C >>> READ FILE30
C
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,IC,200,101,LWORD)
      MPOINT = IC(2)
      MCONST = IC(3)
      MCALCS = IC(4)
      NCALCS = IC(5)
      NT     = IC(18)
      NAT    = IC(19)
      NTAO   = IC(22)
      NSHELL = IC(27)
      NIRRED = IC(28)
      NPRIM  = IC(32)
      NSYMHF = IC(41)
      MXCOEF = IC(42)
C
      NBF=NT
      NT2=NT*(NT+1)/2
C
C >>> START CORE ALLOCATION
C
      WTEMP  = 1
      COORD  = WTEMP + NT
      ZAN    = COORD + 3*NAT
      DVEC   = ZAN + NAT
      WVEC   = DVEC + NT*NT
      XVEC   = WVEC + NTAO*NT
      XARR   = XVEC + NTAO*NT
      EEE    = XARR + NT
      ONEI   = EEE + NT
      SCFTCC = WPADTI(ONEI+NT2)
      ORBSYM = SCFTCC + NT
      ILAST  = IADTWP(ORBSYM+NT)
C
      CALL READ30(ITAP30,MPOINT,MCONST,MCALCS,NCALCS,NT,NTAO,
     .            NSYMHF,MXCOEF,CC(EEE),NLAMDA,6,NCL,NO,
     .            CC(WTEMP),IC(SCFTCC),ITYP,IC(ORBSYM),FLOV,NIRRED,
     .            ENUC,ESCF,NAT,NSHELL,NPRIM,
     .            CHAR,NORD,CC(ZAN),CC(COORD),
     .            IOPEN,NOP)
C
      NV=NT-NO
C
      LTYP=7
      NIOBF=100
      LENGTH=NIOBF*(SEC2I(1)/INTOWP(1))
      IBKTSP=LENGTH*LTYP
C
C >>> CORE ALLOCATION FOR SORT
C
      JBKT   = ILAST
      JIBKT  = WPADTI(JBKT)
      JBUF   = JBKT + IBKTSP
      JIBUF  = WPADTI(JBUF)
      JTOP   = JBUF + INTBUF
C
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)' REQUIRED  MEMORY FOR CCSORT IS',JTOP
      LEFT=MAXCOR-JTOP
      IF(LEFT.LT.0)THEN
         WRITE(6,*)'NOT ENOUGH CC, REQUIRES',JTOP, '  REAL*8 WORDS'
         call qabort
      ENDIF
C
      CALL SORTI(CC(JBKT),IC(JIBKT),CC(JBUF),IC(JIBUF),LENGTH,
     .           LTYP,IBKTSP,INTBUF,NBF,NO,ITAP60,ITAP61,ITAP62,ITAP63,
     .           ITAP64,ITAP65,ITAP66,ITAP78,JOUT,CC(ONEI))
C
c >>> resort b ints to file71
c
      if(code71.eq.0)then
      BUF   = ILAST
      IBUF  = WPADTI(BUF)
      IPQ   = WPADTI(BUF + LENGTH)
      ITV   = IPQ + NBF
      BOFF  = ITV + NT
      BADD  = BOFF + NV*NV
      GBUCK = IADTWP(BADD + NV*NV)
      JTOP  = GBUCK + NV*NV
C
      WRITE(6,*)' REQUIRED  MEMORY FOR RE-SORTING B INTS IS ',JTOP
      LEFT=MAXCOR-JTOP
      WRITE(6,*)' REMAINING MEMORY IS ',LEFT
      IF(LEFT.LT.0)THEN
         WRITE(6,*)'NOT ENOUGH MEMORY, REQUIRES',JTOP, '  REAL*8 WORDS'
         call qabort
      ENDIF
C
      CALL SORTB(CC(BUF),IC(IBUF),LENGTH,
     .           NBF,NO,ITAP61,ITAP71,JOUT,
     .           NIRRED,FLOV,IC(ORBSYM),IC(IPQ),IC(ITV),NV,IC(BOFF),
     .           IC(BADD),LEFT,CC(JTOP),CC(GBUCK))
      endif
C
      CALL RCLOSE(ITAP78,3)
      CALL RCLOSE(ITAP60,3)
      CALL RCLOSE(ITAP61,3)
      CALL RCLOSE(ITAP62,3)
      CALL RCLOSE(ITAP63,3)
      CALL RCLOSE(ITAP64,3)
      CALL RCLOSE(ITAP65,3)
      CALL RCLOSE(ITAP66,3)
C
      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE(6,*)
        WRITE(6,*)'  >>>>>   SORTING OF INTEGRALS  COMPLETED   <<<<<'
        WRITE(6,*)
        endif
      CALL TSTOP(6)
      call psdone()
      stop
C
      END
