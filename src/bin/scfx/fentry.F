cets083089      PROGRAM SCFEX
      subroutine fentry(z,iz,mxcor)
C***********************************************************************
C   NOTICE OF PROGRAM MODIFICATION                                     *
C***********************************************************************
c Moved to PSI distribution disk on 020289 - clj.
C***********************************************************************
C      BY:  RICHARD REMINGTON                                          *
C    DATE:  OCTOBER   4, 1987                                          *
C  REASON:  INCREASE VARIABLE "PSIZE" TO HANDLE MORE D-SHELLS  16=>32  *
C           AND CHANGE PRINT FORMAT STATEMENTS WITH 1P2E... TO 1P,2E...*
C***********************************************************************
C     BY: WESLEY ALLEN                                                 *
C     DATE: JANUARY 5, 1987                                            *
C     REASON: PUT IN PLACE A MORE SOPHISTICATED DAMPING (AVERAGING)    *
C             PROCEDURE.                                               *
C***********************************************************************
C     BY: WESLEY ALLEN                                                 *
C     DATE: DECEMBER 4, 1986                                           *
C     REASON: A BUG WHICH AFFECTED ONLY CERTAIN CASES WAS FIXED.       *
C             THE ERROR CONCERNED THE CONSTRUCTION OF THE CLOSED       *
C             SHELL FOCK OPERATOR USED IN THE OPTIMIZATION OF THE      *
C             DOUBLY OCCUPIED ORBITALS IN THE SAME SYMMETRY BLOCK      *
C             AS THE TWO OPEN SHELL ORBITALS.                          *
C***********************************************************************
C     BY: WESLEY ALLEN                                                 *
C     DATE: NOVEMBER 5, 1986                                           *
C     A NEW VERSION OF SCFX WHICH INCLUDES SEVERAL CHANGES             *
C     FOR CI GRADIENTS AND OTHER THINGS HAS BEEN ADDED.                *
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
cets083089      REAL*8 Z(200000)
      real*8 z(mxcor)
      INTEGER IZ(1)
cets083089      EQUIVALENCE (Z,IZ)
cets083089      MAXCOR=400000
      maxcor=mxcor*2
cets053193      call drum
      call psinit('APPEND')
C      CALL TSTART(3)
      CALL TSTART(6)
cets083089      CALL NOUNFL
      I30 = 30
      I34 = 34
      I48 = 48
      I49 = 49
      I92 = 92
      CALL RFILE(I30)
      CALL RFILE(I48)
      CALL RFILE(I49)
      CALL RFILE(I92)
C
      CALL PITZER(Z,IZ,MAXCOR)
C
      JCODE = 3
      CALL RCLOSE(I30,JCODE)
      CALL RCLOSE(I34,JCODE)
      CALL RCLOSE(I48,JCODE)
      JCODE = 4
      CALL RCLOSE(I92,JCODE)
C
C      CALL TSTOP(3)
      CALL TSTOP(6)
C
      STOP
      END
