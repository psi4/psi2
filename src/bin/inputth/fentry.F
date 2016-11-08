cges  PROGRAM INPUT
      subroutine fentry(z,ia,mxcor)
C    (Program INPUTTH)
C**********************************************************************
C* NOTICE OF PROGRAM MODIFICATION                                     *
C**********************************************************************
c
c	by: croag
c	date: Mon Mar  7 12:37:26 EST 1994
c	where: below and in the function intrb1
c	reason: allow for generic local and global basis set
c		definitions.  path should now be:
c	./input.dat -- basis:
c	./basis.dat -- basis:
c	User_set_path/user_set_name -- basis:
c	-or- User_set_path/basis.dat -- basis: 
c	PSI_LIB/pbasis.dat -- basis:
c
c666 666 666 666 666 666 666 666 666 666 666 666 666 666 666 666 666 666 
C     BY:CChaos
C   DATE: December 3, 1992
C REASON: Write Nuc. Rep. Energy 
C         call to enuc function in subroutine mole.
C**********************************************************************
C     BY: ETS                                Search: cets052791       *
C   DATE: May 27, 1991                                                *
C REASON: diddle with the way the basis set is read in                *
c         also add blkatm.f, common block with atom names and symbols
c         used to figure out charge of atom if CHARGES vector not
c         given in input.  Now can also use atomic symbols to specify
c         basis set.  I only added atoms up through lanthanum, someone
c         else can go up to number 103 if they so desire.
C**********************************************************************
C     BY: CLJ                                Search: clj021091        *
C   DATE: March 10, 1991                                              *
C REASON: New input format and print options                          *
C**********************************************************************
C     BY: GES                                Search: cges             *
C   DATE: may 24, 1990                                                *
C REASON: Augment maximum number of primitives to 35                  *
C**********************************************************************
C     BY: GES                                Search: cges             *
C   DATE: mar 27, 1990                                                *
C REASON: If keyword 'ANO' is used with 'GET' change gralized contr.  *
c         format to MOLECULE (Nasa ANOs may be used straightforwardly)*
C**********************************************************************
C     BY: GES                                Search: cges             *
C   DATE: dec 12, 1989                                                *
C REASON: Augment maximum number of primitives to 30                  *
c         and maximum number of generalized contractions to 9         *
C**********************************************************************
C     BY: Brenda Colegrove                   Search: c06-09-89        *
C   DATE: June 9, 1989                                                *
C REASON: Get rid of references to TIMLIM and NE                      *
C**********************************************************************
C     BY: Richard Remington                  Search: c03-14-89        *
C   DATE: March     14, 1989                                          *
C REASON: CHANGE CORE TO RUN IN 7MB ON 9370 OR 15MB ON 3090           *
***********************************************************************
C     BY: Richard Remington                  Search: c01-24-89        *
C   DATE: January 24, 1989                                            *
C REASON: Remove printing of angles for atoms and diatomics.          *
C**********************************************************************
C   CHANGE THE PROGRAM SO THAT YOU CAN HAVE SOME SETS OF BASIS        *
C   FUNCTIONS AS D6 AND SOME AS D5, ETC.  FIRST AND SECOND SCF        *
C   DERIVATIVES NOW WORK FOR SEGMENTED CONTRACTIONS D5,D6,F7,F10,G9,  *
C   AND G15.  Change the defaults back to d6, f10, g15.               *
***********************************************************************
C   THE PROGRAM ONLY READS AND PRINTS 6 GENERAL CONTRACTIONS NOW.     *
C   MRCRU = 6.  THE READING AND PRINTING IS ALSO SPECIFIC TO A        *
C   MAXIMUM OF 6.  THERE IS A MAXIMUM OF 20 PRIMITIVES PER SHELL.     *
***********************************************************************
C   SO FAR THIS SHOULD WORK FOR SCF, CI, AND CC ENERGY POINTS         *
C   WORK ON DERIVATIVES IS IN PROGRESS.  This also should work with   *
C   the old gradient codes for d6 functions.                          *
C**********************************************************************
C     BY:  RUSS PITZER AND TRACY HAMILTON                             *
C   DATE:  AUGUST 13, 1988                                            *
C REASON:  MODIFICATION FOR 5 COMPONENT D AND 7 COMPONENT F           *
C          FUNCTIONS, GENERAL CONTRACTIONS AND CLEAN UP PROGRAM       *
C**********************************************************************
C         MODIFICATION FOR F AND G FUNCTIONS                          *
C     BY:  PAUL SAXE AND JULIA RICE                                   *
C   DATE:  NOVEMBER 12,1986                                           *
C**********************************************************************
C NOTICE OF PROGRAM CREATION                                          *
C     BY:  ROBERT MACLAGAN                                            *
C          FROM FEBRUARY 5, 1985 VERSION OF 'INPUT'                   *
C   DATE:  FEBRUARY 10, 1984                                          *
C ENTERED INTO SCHAEFER AS INPUT ON 3-7-85                            *
C**********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
c03-14-89  (For either 7Mb or 15Mb)  DIMENSION Z(400000),IA(1)
cges  DIMENSION Z(400000),IA(1)
      DIMENSION Z(mxcor),IA(*)
#include <error.h>
#include <misc.h>
      integer inold,ierr,prcntr,prdup,isf7,isf8,isf9
      integer errcod,fioopn,fiocls
clj021091 This common will be used to pass the input format style around
      common /instyl/ inold, isf9, isf8
clj021091 Use psinit and fioopn instead of ffile
      errcod = fioopn(6, 'output.dat', 'USER', 'WRITE')
      errcod = fiocls(6, 'DELETE')
      call psinit('WRITE')
      CALL TSTART(6)

clj021091 Regardless of the input style, use the new style print control
c     The BRIEF prcntr flag is ON by default.
      ierr = prcntr('ON BRIEF ANGSTROM DEGREE SHELLS')
      ierr = prdup()
      if (prcntr('IS_ON DEBUG').eq.1) then
        ierr = prcntr('OFF BRIEF')
        endif
      if (prcntr('IS_ON BRIEF').eq.0) then
        ierr = prcntr('ON BOHR')
        ierr = prcntr('ON RADIAN')
        endif
c call prdup a second time to allow users to override the defaults set up above
      ierr = prdup()

clj021091 Determine if we are using the old style input
      call locate(5,'# INPUT ##',ierr)
      if (ierr.ne.0) inold = 0
      if (ierr.eq.0) inold = 1

      if (inold.eq.1) then
c       open basis.dat in the current directory
        write(6,*) 'WARNING: using old style ''# INPUT ##'''
        if (prcntr('IS_ON BRIEF').eq.0)
     &    write(6,*) 'basis set data is in the user file basis.dat'
croag        call fioopn(7,'basis.dat','USER','READ')
        errcod =  fioopn(7,'basis.dat','USER','READ')
      else
c       open pbasis.dat in the psi library directory

croag   open a hell of alot more than just that.  use the function
croag	return wisely.  use isf7 & 8 are check flags used in readbasis
croag	file7 in the LIB area should always exist, but put in the flag anyway.
croag	these are is-f(ile)-7, etc.  1 = exists, 0 = non.
croag	unit9 is ./basis.dat
croag	unit8 is User_Defined basis file
croag	unit7 is pbasis.dat in the psi library directory

        if (prcntr('IS_ON BRIEF').eq.0)
     &    write(6,*) 'basis set data is in the library file pbasis.dat'
croag        call fioopn(7,'pbasis.dat','LIB','READ')
	isf7 = 1
	isf8 = 1
	isf9 = 1
        errcod = fioopn(9,'basis.dat','USER','READ')
        if (errcod .ne. EOK) isf9 = 0
c	write(6,*) 'isf9 is ',isf9
        errcod = fioopn(8,'basis.dat','SUSER','READ')
        if (errcod .ne. EOK) isf8 = 0
c	write(6,*) 'isf8 is ',isf8
        errcod = fioopn(7,'pbasis.dat','LIB','READ')
        if (errcod .ne. EOK) isf7 = 0
c	write(6,*) 'ierr7 is ',isf7
        endif

C
C***********************************************************************
C           CALL THE MAIN SUBROUTINE, PASSING IN CORE                  *
C***********************************************************************
      maxcor=intowp(mxcor)
c530 unbelievable! you can't call main a subroutine
      CALL xMAIN(Z,IA,MAXCOR)
      ITAP30 = 30
      JCODE = 3
      CALL RCLOSE(ITAP30,JCODE)
      CALL TSTOP(6)
      call psdone()
cges  CALL EXIT
      stop
C
      END
