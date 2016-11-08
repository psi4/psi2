cets The following is called by a C program which handles dynamic
cets allocation of memory.
      subroutine fentry(rz,iz,maxcr)
cets  PROGRAM TRANS
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C   Last updated on Febuary 03, 1995 by Evan Bolton                   *
C   Reason:  Addition of CASSCF WFN type required that the WFN is to  *
C              be read in and when the WFN = CASSCF, then file34      *
C              must NOT be deleted.                                   *
C              See comments labeled CeebCAS.                          *
C**********************************************************************
c     by:  Edward Seidl                  search: cets070891
c   date:  July 8,1991                                      
c reason:  hardwire to use <= 8MB of core.  option LARGE_CORE added
c          will use amount of core specified by MEMORY if LARGE_CORE is
c          true
c
c     by:  Edward Seidl                  search: cets030291
c   date:  March 2,1991                                      
c reason:  add new input format, and replace drum with psinit
c
c     by:  Edward Seidl                              search: cets
c   date:  October 14,1988                                      
c reason: remove it34 after getint to save disk space if       
c           itest is set = 1 in trans input                   
c           also remove file95                               
c**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c4-19-88  *
C*  DATE:  April  19,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C      BY:  RICHARD REMINGTON                                         *
C    DATE:  FEBRUARY  4, 1987                                         *
C  REASON:  MODIFY TO HANDLE LARGER CASES (MAYBE---3M INTEGER WORDS)  *
C**********************************************************************
C      BY:  TJL                                                       *
C    DATE:  AUGUST 3,1985                                             *
C  REASON:  MODIFY TO HANDLE LARGER CASES                             *
C**********************************************************************
C      BY:  RICHARD REMINGTON                                         *
C    DATE:  FEBRUARY 12, 1985                                         *
C  REASON:  CLEANUP OUTPUT AND STANDARDIZE UNIT 6 TO OUTPUT AND       *
C     UNIT 3 TO CHECK                                                 *
C**********************************************************************
C      BY:  RICHARD REMINGTON                                         *
C    DATE:  FEBRUARY 13, 1985                                         *
C  REASON:  CHANGE UNIT NUMBERS !! 2=>52, 8=>58, 44=>34               *
C           CHANGE UNIT NUMBERS !! 10=>30, 105=85=>95, 101=81=>91     *
C           CHANGE UNIT NUMBERS !! 93=>97, 103=83=>93,                *
C**********************************************************************
C
C
      COMMON /TAPES/  IT34,IT58,IT47,IT52,IT91,IT93,IT95
      COMMON /OUT/    IOUT6,ICHECK,IDBG,ICNT
      COMMON /TEST/   ITEST
C
CUCB  REAL*8 RZ(1500000)
c4-19-88  REAL*8 RZ(1200000)
cets  REAL*8 RZ(725000)
cets  INTEGER IZ(1)
      real*8 rz(maxcr)
      integer iz(1)
cets  EQUIVALENCE (RZ,IZ)
cets030291
#include <error.h>
      integer errcod, frdboo, frdc
CeebCAS
      character*8 cityp
C
C     WRITE(ICHECK,101)
C     READ(1,*) NASK
C     IF (NASK.LE.0  .OR. NASK.GT.190) NASK=190
cets030291      call drum
      call psinit('APPEND')
      ICHECK = 3
      IOUT6  = 6
c     CALL TSTART(ICHECK)
      CALL TSTART(IOUT6)
      itest=0
CUCB  NASK=1500
c4-19-88  NASK=1200
      NASK=725     
cets  ITEST = 0
C     WRITE(ICHECK,102) NASK
 101  FORMAT(/'  ENTER NASK'/)
 102  FORMAT(/'  NASK=',I5)
CTJL  CALL X:GDSPCE(NASK,NGET,,)
C
cets  MAXCOR=NASK*2000
      input = 5
      icore=0
      call locate (input,'# TRANS ##',IERR)
      if (ierr.eq.0) then
        read(5,199) itest
  199   format(i5)
        write(6,*) 'itest =',itest
cets030291
      else
        itest=1
        errcod = frdboo('DELETE34',itest)

CeebCAS  Begin
        errcod = frdc('WFN',cityp)
        if (cityp.eq.'CASSCF  ') itest=0

        write(6,*) '  '
        write(6,*) 'WFN = ',cityp
        write(6,*) '  '
CeebCAS  End

        errcod = frdboo('LARGE_CORE',icore)
      endif
c
cets070891
#ifdef AIX
      if(maxcr.gt.1000000.and.icore.eq.0) then
        write(6,*) ' '
        write(6,*) 'You have requested more than 8 MB of memory'
        write(6,*) 'which will cause poor performance on the'
        write(6,*) 'rs6000s, thus I will only use 8 MB of core.'
        write(6,*) 'If this is not enough, may I suggest using'
        write(6,*) 'the LARGE_CORE option (see psi manual)'
        write(6,*) ' '
        maxcr=1000000
      endif
#endif
      maxcor=intowp(maxcr)
      CALL NEOTRN(RZ,IZ,MAXCOR)
C
      JCODE = 3
      IT30 = 30
      IT97 = 97
      CALL RCLOSE(IT30,JCODE)
cets got rid of file34 after getint in neotrn
cets  CALL RCLOSE(IT34,JCODE)
      CALL RCLOSE(IT52,JCODE)
      CALL RCLOSE(IT58,JCODE)
      JCODE = 4
CTJL  CALL RCLOSE(IT91,JCODE)
      CALL RCLOSE(IT93,JCODE)
cets file95 closed in neotrn
cets  if(itest.eq.0) CALL RCLOSE(IT95,JCODE)
C
      WRITE(IOUT6,105)
  105 FORMAT(/,1H ,12X,'>>>>>  TRANSFORMATION ',
     +'CALCULATION COMPLETED  <<<<< ')
C
c     CALL TSTOP(ICHECK)
      CALL TSTOP(IOUT6)
C
      STOP
      END
