      subroutine fentry(z,iz,maxcr)
C  THIS PROGRAM IS REALLY SCFTFK FORT-OLD, SCFTFK AND NO OPEN SHELL DIIS
C***********************************************************************
C SCF PROGRAM WITH TFOCK SUBROUTINE PRESENTLY NEEDED FOR CIGRADS 3-7-85
C***********************************************************************
C   NOTICE OF PROGRAM MODIFICATION                                     *
C***********************************************************************
c      by:  Gustavo E. Scuseria
c    date:  03 jul 1991
c  reason:  modify code to write out fock alpha and beta matrices
c           into file79 (leaving space for one-el ham)
c           intwo was modified for case=4 and now both val1 and val2
c                 are written to disk when the supermatrix is formed
c           contr2 is modified to accomodate new super-matrix format
c                  for case=4
c           contr2g same as new contr2 but q-terms are added for fock
c                   matrices
c           (Search for cges)
C***********************************************************************
c      by:  Gustavo E. Scuseria
c    date:  12 dec 1990
c  reason:  Rotate mo's for open-shell CC. Add variable # PSIZE #
c           (Search for cges)
C***********************************************************************
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
      MAXCOR=intowp(maxcr)
C
CTJL  NASK=200
CTJL  CALL X:GDSPCE(NASK,NGET,,)
C
      call drum
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
