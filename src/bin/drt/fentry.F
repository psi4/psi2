cets103089      PROGRAM DRT
      subroutine fentry(rz,z,maxcor)
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C      BY:  Evan Bolton  (See CeebCAS)                                *
C    DATE:  Febuary 3, 1995                                           *
C  REASON:  Made changes necessary for WFN = CASSCF.  Namely if       *
C           cityp.eq.'CASSCF  ' then res=.true.                       *
C           I also updated the manual page to reflect that DRT reads  *
C           in the derivative type and WFN type.                      *
C**********************************************************************
C      BY:  Edward Seidl (cets050191)
C    DATE:  May 1, 1991
C  REASON:  add routine scancj using CLJ's new input routines
c           hopefully open-shell singlets are done sensibly now.
C**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c3-24-88  *
C*  DATE:  MARCH  24,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C      BY:  Curtis Janssen                                            *
C    DATE:  Aug. 25, 1987 at UGA                                      *
C  REASON:  Remove a branch to the inside of a do loop in WEIGHT.     *
C**********************************************************************
C      BY:  RICHARD REMINGTON                                         *
C    DATE:  FEBRUARY 13, 1985                                         *
C  REASON:  CLEANUP OUTPUT AND STANDARDIZE UNIT 6 TO OUTPUT AND       *
C     UNIT 3 TO CHECK                                                 *
C**********************************************************************
C      BY:  RICHARD REMINGTON                                         *
C    DATE:  FEBRUARY 12, 1985                                         *
C  REASON:  CHANGE UNIT NUMBERS !! 10=>30, 8=>58                      *
C**********************************************************************
C
C
      IMPLICIT INTEGER (A-Z)
C
c3-24-88  COMMON  Z(90000)
cets103089      COMMON  Z(90000)
      dimension z(maxcor)
C
cets050191      call drum 
      call psinit('APPEND')
cets  CALL TSTART (3)
      CALL TSTART (6)
C
c3-24-88  CALL MNDRT(Z,90000)
cets103089      CALL MNDRT(Z,90000)
      call mndrt(z,maxcor)
C
cets  CALL TSTOP (3)
      CALL TSTOP (6)
C
C     CALL EXIT
C
      STOP
      END














