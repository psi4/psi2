cets The following is called by a C program which handles dynamic
cets allocation of memory.
      subroutine fentry(a,id,maxcor)
cets  PROGRAM ONEPDM
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
c      by: Roger Grev (crsg021594)
c    date: February 15, 1994
c  reason: to correct the dimension of the array BFCODE, so that
c          multireference CI onepdm output would have the correct
c          symmetry and number codes.
c
c      by: Edward Seidl (cets071791)
c    date: July 17, 1991
c  reason: pass integer array id to mainop, get rid of common block /drt/
c          large numbers of f functions cause leni to exceed maxi which
c          was hardwired.  Now all is dynamic (hopefully no overwrites)
c
c      by: Edward Seidl (cets030291)
c    date: March 2, 1991
c  reason: add new input, and use psinit instead of drum
c
C*  BY:  RICHARD REMINGTON                         search:  c3-24-88  *
C*  DATE:  MARCH  24,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C*  NOTICE OF PROGRAM MODIFICATION                                    *
C*       BY: WESLEY ALLEN                                             *
C*     DATE: NOVEMBER 5, 1986                                         *
C*   REASON: ADD A WARNING TO THE USER WHO WRITES NATURAL ORBITALS    *
C*           OVER THE MOLECULAR ORBITALS IN FILE30.                   *
C**********************************************************************
C*  NOTICE OF PROGRAM MODIFICATION                                    *
C*       BY: WESLEY ALLEN                                             *
C*     DATE: MAY 27, 1985                                             *
C*   REASON: SET UP DEFAULTS IF "ONEPDM" INPUT NOT IN INPUT DECK      *
C**********************************************************************
C*  NOTICE OF PROGRAM MODIFICATION                                    *
C*       BY: WESLEY ALLEN                                             *
C*     DATE: AUGUST 16, 1985                                          *
C*   REASON: MODIFICATIONS TO LABEL THE IRREDUCIBLE REPRESENTATIONS   *
C*           OF DIFFERENT POINT GROUPS CORRECTLY. ALSO, ITAP'S WERE   *
C*           CHANGED TO IFIL'S WITH THE CORRECT NUMBERS.              *
C**********************************************************************
C*  NOTICE OF PROGRAM MODIFICATION                                    *
C*       BY: WESLEY ALLEN                                             *
C*     DATE: SEPTEMBER 17, 1985                                       *
C*   REASON: OPTIONS ADDED TO WRITE THE AO 1-PDM TO THE MASTER FILE   *
C*           AND THE NATURAL ORBITALS OVER THE SCF VECTOR IN FILE30.  *
C**********************************************************************
C
c3-24-88  REAL*8 A(100000)
cets  REAL*8 A(100000)
      real*8 a(maxcor)
cets added dummy
cets  integer id(intowp(maxcor))
      integer id(1)
C
    1 FORMAT(/,1H ,16X,'>>>>>  ONEPDM CALCULATION COMPLETED  <<<<< ')
C
c3-24-88  MAXCOR=100000
cets  MAXCOR=100000
cets030291      call drum
      call psinit('APPEND')
c     CALL TSTART(3)
      CALL TSTART(6)
      CALL MAINOP(A,id,MAXCOR)
      WRITE(6,1)
c     CALL TSTOP(3)
      CALL TSTOP(6)
C
      STOP
      END
