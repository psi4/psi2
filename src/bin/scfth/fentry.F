      subroutine fentry(z,iz,maxcr)
C  THIS PROGRAM IS REALLY SCFTFK FORT-OLD, SCFTFK AND NO OPEN SHELL DIIS
C***********************************************************************
C SCF PROGRAM WITH TFOCK SUBROUTINE PRESENTLY NEEDED FOR CIGRADS 3-7-85
C***********************************************************************
C   NOTICE OF PROGRAM MODIFICATION                                     *
C***********************************************************************
c      by:  Edward Seidl (cets070891)
c    date:  July 8, 1991
c  reason:  hardwire for only 1 megabyte of core on rs6000. more core
c           than this results in cpu/wall time ratios of about 9%
c
c      by:  Edward Seidl (cets050391)
c    date:  May 3, 1991
c  reason:  double "psize" in pitzer.f to allow for large ano basis sets
c
c      by:  C Janssen
c    date:  March 9, 1991
c  reason:  modified virtul.F for variable linewidth and removed
c           all calls to PRINTS and replaced these with calls to PRINT
c           in iomr.
c
c      by:  Edward Seidl (cets030291)
c    date:  March 3, 1991
c  reason:  1) merge scf and scftfk. now does orbital rotation if
c              corrtype is not equal to 'SCF'
c           2) do not call tfock if not a ci gradient
c           3) write a 1 to i10(45) so scf knows an old vector is in
c              file30.  this way can keep RESTART set to 'yes' all
c              the time
c
c      by:  Edward Seidl (cets022791)
c    date:  Feb. 27, 1991
c  reason:  add new-style input to input.f, scf1.f, and pitzer.f
c
c      by:  Edward Seidl (cets022291)
c    date:  Feb. 22, 1991
c  reason:  replace call to drum with call to psinit
c
C      BY:  RICHARD REMINGTON                                          *
C    DATE:  April    24, 1989                                          *
C  REASON:  Add  "common/tapes/..." to subroutines DMAT, ERDUW, CONTR2,*
C           SQTOTQ, and SQTOTR  so that variables ICHECK and IOUT will *
C           be defined in write statements. Add write(iout, ) in some. *
C***********************************************************************
C      BY:  RICHARD REMINGTON                     SEARCH: MACHINESIZE  *
C    DATE:  APRIL    21, 1989                                          *
C  REASON:  DOUBLE Z() AND MAXCOR FROM 200000 REAL*8 TO 400000 REAL*8. *
C***********************************************************************
C      BY:  RICHARD REMINGTON                     Search: C01-16-89    *
C    DATE:  January  16, 1989                                          *
C  REASON:  OPEN AND CLOSE FILE30 BETWEEN ITERATIONS to write to file6.*
C***********************************************************************
C      BY:  RICHARD REMINGTON                                          *
C    DATE:  APRIL   3,  1988                                           *
C  REASON:  DUE TO PROBLEMS WITH OPEN SHELLS AND DIIS, THE CURRENT SCF *
C           PROGRAM HAS BEEN RENAMED  "SCFTFK FORTDIIS"  AND THE OLD   *
C           VERSION NAMED  "SCFTFK FORT-OLD"  HAS BEEN COPIED OVER AS  *
C           "SCFTFK FORTRAN".    THIS PROGRAM IS THE SAME AS           *
C                                   "SCFTFK FORT-OLD"                  *
C***********************************************************************
C      BY:  RICHARD REMINGTON                                          *
C    DATE:  OCTOBER   4, 1987                                          *
C  REASON:  INCREASE VARIABLE "PSIZE" TO HANDLE MORE D-SHELLS  16=>32  *
C           AND CHANGE PRINT FORMAT STATEMENTS WITH 1P2E... TO 1P,2E...*
C***********************************************************************
C      BY:  WESLEY ALLEN                                               *
C    DATE:  NOVEMBER 5, 1986                                           *
C  REASON:  FIX BUG IN THE REORDERING OF VIRTUAL ORBITALS              *
C***********************************************************************
C      BY:  TJL AND ACS                                                *
C    DATE:  FEBRUARY 24, 1986                                          *
C  REASON:  MODIFY FOR TCSCF-CI GRADIENTS AND CI SECOND DERIVATIVES    *
C***********************************************************************
C      BY:  TJL                                                        *
C    DATE:  JULY 13, 1985                                              *
C  REASON:  CORRECT THE DROL ARRAY USED IN GRCPHF                      *
C***********************************************************************
C      BY:  ROBERT MACLAGAN                                            *
C    DATE:  FEBRUARY 27, 1985                                          *
C  REASON:  PROBLEMS IN "TFOCK" FOR "TCSCF"                            *
C***********************************************************************
C      BY:  RICHARD REMINGTON                                          *
C    DATE:  FEBRUARY  1, 1985                                          *
C  REASON:  CLEANUP OUTPUT AND STANDARDIZE UNIT 6 TO OUTPUT AND        *
C     UNIT 3 TO CHECK; REMOVE DUPLICATE MATOUT, RENAME SUBROUTINE      *
C     "PRINT" AS "PRINTS" AND "FLIN" AS "FLIN2"                        *
C     ALSO VARIABLE "NXTFLG" IS PRINTED OUT AS THE OPTION "DIIS"       *
C***********************************************************************
C    DATE:  FEBRUARY  4, 1985                                          *
C  REASON:  DYNAMICALLY DIMENSION WMO AND WMO2 IN TFOCK AND ALPHA AND  *
C           BETA IN INPUT AND LBLI IN TWOI                             *
C***********************************************************************
C      BY:  RICHARD REMINGTON                                          *
C    DATE:  FEBRUARY  6, 1985                                          *
C  REASON:  CHANGE UNIT NUMBERS !!10=>30,44=>34,57=>47 ! ! ! ! !       *
C           AND STANDARDIZE "SCF" INPUT TO "I5" FORMAT                 *
C           "IDFIX", "ITPREP", AND "ITCON" ARE NOT READ FROM SCF INPUT *
C           BUT ATTAIN THE DEFAULT VALUES GIVEN BY THE PROGRAM         *
C***********************************************************************
C      BY:  RICHARD REMINGTON                                          *
C    DATE:  FEBRUARY 12, 1985                                          *
C  REASON:  ADD LINES TO SUBROUTINE TFOCK FOR OPEN/CLOSED SCF2ND CALC  *
C***********************************************************************
C      BY:  RICHARD REMINGTON                                          *
C    DATE:  FEBRUARY 18, 1985                                          *
C  REASON:  MOVE RFILE(49) TO TOP OF PROGRAM FROM SCF1 AND LEAVE A SREW*
C***********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
CMACHINESIZE  (For  7MB)  REAL*8 Z(200000)
CMACHINESIZE  (For 15MB)  REAL*8 Z(400000)
      REAL*8 Z(maxcr)
      INTEGER IZ(1)
c     EQUIVALENCE (Z,IZ)
C
CMACHINESIZE  (For  7MB)  MAXCOR=400000
CMACHINESIZE  (For 15MB)  MAXCOR=800000
C
CTJL  NASK=200
CTJL  CALL X:GDSPCE(NASK,NGET,,)
C
cets022291      call drum
      call psinit('APPEND')
C      CALL TSTART(3)
      CALL TSTART(6)
C
      I30 = 30
      I34 = 34
      I49 = 49
      I92 = 92
      CALL RFILE(I30)
      CALL RFILE(I49)
      CALL RFILE(I92)
cets070891
#ifdef AIX
      if(maxcr.gt.1000000) then
        write(6,*) ' '
        write(6,*) 'You have requested more than 8 MB of memory'
        write(6,*) 'which will cause poor performance on the'
        write(6,*) 'rs6000s, thus I will only use 8 MB of core'
        write(6,*) 'if this is not enough, may I suggest using'
        write(6,*) 'the program CSCF.'
        write(6,*) ' '
        maxcr=1000000
      endif
#endif
      MAXCOR=intowp(maxcr)
C
      CALL PITZER(Z,IZ,MAXCOR)
C
      JCODE = 3
      CALL RCLOSE(I30,JCODE)
      CALL RCLOSE(I34,JCODE)
      JCODE = 4
      CALL RCLOSE(I92,JCODE)
C
C      CALL TSTOP(3)
      CALL TSTOP(6)
C
      STOP
      END
