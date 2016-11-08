cets110990      PROGRAM GUGACI
      subroutine fentry(z,iz,maxcor)
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
c       by: Edward Seidl (cets030291)
c     date: March 3, 1991
c   reason: add new input format, and replace drum with psinit
c
c       by: Edward Seidl (cets111090)
c     date: Nov. 10, 1990
c   reason: split up loopyd and loopy, add common block 'sqrts'
c           and common block 'ets'
c
c       by: Edward Seidl (cets110990)
c     date: Nov. 9, 1990
c   reason: convert Oct version of gugaci to UNIX
c           also delete file94 when finished
c           no longer allow restarts
c           convert as many entry points to subroutine calls as possible
c**********************************************************************
C   BY:  RICHARD REMINGTON                         SEARCH: C10-13-90  *
C   DATE:  OCTOBER  13, 1990                                          *
C   REASON: FIX ERRORS FOR LARGE CASES.  ARRAYS IN STATPG INCREASED;  *
C           SUBROUTINE WFILE MODIFIED AND RENAMED WFILE0.             *
C**********************************************************************
C   BY:  Curtis Janssen                            key: clj120788     *
C   DATE:  Dec. 7, 1988                                               *
C   REASON: Program searches for '# GUGAREF ' in file 5.  If found,   *
C           one free format integer is read which is the drt walk     *
C           be used as the initial guess vector. If 0 is read, the    *
C           lowest energy configuration will be used as the reference.*
C**********************************************************************
C   BY:  RICHARD REMINGTON                         search:  c3-24-88  *
C   DATE:  MARCH  24,  1988                                           *
C   REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C   NOTICE OF PROGRAM MODIFICATION                                    *
C      BY:  RICHARD REMINGTON                                         *
C    DATE:  MAY   2,  1986                                            *
C  REASON:  ADD PRINT OPTION TO REDUCE OUTPUT DATA (THE 9TH OPTION)   *
C               A DEFAULT OF IPRINT = 0 GIVES MINIMUM OUTPUT          *
C**********************************************************************
C      BY:  RICHARD REMINGTON                                         *
C    DATE:  FEBRUARY 13, 1985                                         *
C  REASON:  CHANGE UNIT NUMBERS !! 2=>52, 03=99=>99, 06=96=>96/94     *
C  REASON:  CHANGE UNIT NUMBERS !! 8=>58, 04=94=>94                   *
C  REASON:  CHANGE UNIT NUMBERS !!12=>54, 05=95=>95                   *
C**********************************************************************
C CHANGES MADE IN 2 SEPARATE PLACES BY WESLEY ALLEN 2-27-85 TO        *
C REACTIVATE THE CALCULATION OF MULTIPLE ROOTS IN SUBROUTINE DAVIDS.  *
C**********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /TAPES/ITAP52,ITAPE5,ITAPE6,ITAP58,ITAP54,ITAP99,ITAP94
     *,             ITAPE3,ITAP95,ITAP96
C
c3-24-88 (For  7Mb)  REAL*8 Z(750000) ,CUTOF1,CUTOF2,CUTOF3
c3-24-88 (For 15Mb)  REAL*8 Z(1600000) ,CUTOF1,CUTOF2,CUTOF3
cets110990      REAL*8 Z(1400000) ,CUTOF1,CUTOF2,CUTOF3
cets110990      INTEGER IZ(1)
cets110990      EQUIVALENCE (Z,IZ)
      real*8 z(maxcor),cutof1,cutof2
      integer iz(1)
C
      COMMON /MATRIX/ CUTOF1,CUTOF2,NMULTP,NMULTD
C
cets030291      call drum
      call psinit('APPEND')

      ITAPE3 = 3
      ITAPE6 = 6
C      CALL TSTART(ITAPE3)
      CALL TSTART(ITAPE6)
cets110990      CALL NOUNFL
C
      CALL RFILE(54)
C
      CUTOF1=10.0D+00**(-10)
      CUTOF2=10.0D+00**(-11)
      NMULTP=0
      NMULTD=0
c3-24-88 (For  7Mb)  MAXI=INTOWP(750000)
c3-24-88 (For 15Mb)  MAXI=INTOWP(1600000)
cets110990      MAXI=INTOWP(1400000)
      maxi=intowp(maxcor)
C
      CALL CIMAIN(IZ,Z,MAXI)
C
      JCODE = 3
      CALL RCLOSE(ITAP52,JCODE)
      CALL RCLOSE(ITAP58,JCODE)
      CALL RCLOSE(ITAP54,JCODE)
      JCODE = 4
      call rclose(itap94,jcode)
      CALL RCLOSE(ITAP99,JCODE)
C
C      CALL TSTOP(ITAPE3)
      CALL TSTOP(ITAPE6)
C
cets110990      CALL EXIT
      stop
C
      END
