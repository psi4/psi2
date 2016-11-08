cets The following is called by a C program which handles dynamic
cets allocation of memory.
      subroutine fentry(rz,iz,maxcor)
cets  PROGRAM CISORT
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
c
c     by:  Russ Thomas & Justin Fermann
c   date:  March 15, 1994
c reason:  Put in explicit abort when there's not enough memory
c          to sort the integrals.  search cjrt in main.F
c**********************************************************************
c     by:  Edward Seidl                  search: cets070991
c   date:  July 9,1991
c reason:  hardwire to use <= 8MB of core.  option LARGE_CORE added
c          will use amount of core specified by MEMORY if LARGE_CORE is
c          true
c
c      by:  Edward Seidl (cets030291)
c    date:  March 3, 1991
c  reason:  replace drum with psinit
C**********************************************************************
c      by:  Edward Seidl (cets)
c    date:  August 11, 1988
c  reason:  Add fentry to allow dynamic core allocation
c**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c3-24-88  *
C*  DATE:  MARCH  24,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C      BY:  RICHARD REMINGTON                                         *
C    DATE:  FEBRUARY 12, 1985                                         *
C  REASON:  CLEANUP OUTPUT AND STANDARDIZE UNIT 6 TO OUTPUT AND       *
C     UNIT 3 TO CHECK                                                 *
C**********************************************************************
C      BY:  RICHARD REMINGTON                                         *
C    DATE:  FEBRUARY 13, 1985                                         *
C  REASON:  CHANGE UNIT NUMBERS !! 2=>52, 8=>58, 03=99=>99            *
C**********************************************************************
C
C
      IMPLICIT INTEGER (A-Z)
#include <error.h>
      integer errcod, frdboo, icore
C
CTJL  EXTENDED BLOCK /CORE/ IZ(1)
CTJL  EXTENDED BASE  /CORE/ 256
      COMMON /TAPES/  ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP12,ITAP99,ITAP04
     *,               ITAPE3,ITAP05,ITAP06
C
C  AT THE PRESENT TIME, THE FOLLOWING TAPES ARE "NOT" USED:
C         ITAPE5, ITAP12, ITAP04, ITAP05, ITAP06
C
      COMMON /TEST/ ITEST
C
c3-24-88  REAL*8 RZ(600000)
cets  REAL*8 RZ(600000)
cets  INTEGER IZ(1)
cets  EQUIVALENCE (IZ,RZ)
      real*8 rz(maxcor)
      integer iz(1)
C
CTJL  NASK=200
CTJL  CALL X:GDSPCE(NASK,NGET,,)
C
cets030291      CALL DRUM
      call psinit('APPEND')
      ITAPE3 = 3
      ITAPE6 = 6
c     CALL TSTART(ITAPE3)
      CALL TSTART(ITAPE6)
cets070991
      icore=1
      errcod = frdboo('LARGE_CORE',icore)
C
c3-24-88  MAXR = 600000
cets  MAXR = 600000
cets070991
#ifdef AIX
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
      maxr = maxcor
      MAXI = INTOWP(MAXR)
      ITEST = 0
C
      CALL MAINST(IZ,RZ,MAXI,MAXR)
C
      JCODE = 3
      CALL RCLOSE(ITAP52,JCODE)
      CALL RCLOSE(ITAP58,JCODE)
      JCODE = 4
      IF(ITEST.EQ.1) CALL RCLOSE(ITAP99,JCODE)
C
c     CALL TSTOP(ITAPE3)
      CALL TSTOP(ITAPE6)
C
      stop
C
      END
