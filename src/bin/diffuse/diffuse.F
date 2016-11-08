      subroutine fentry(cc,ia,maxcor)
c**********************************************************************
c     PROGRAM DIFFUSE                                                 *
c**********************************************************************
c   Last updated on May 14, 2008 by Dr. Yukio Yamaguchi               *
c**********************************************************************
c  This program computes exponents of the diffuse functions           *
c   for a given basis set                                             *
c**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension cc(maxcor),ia(maxcor*2)
      parameter (max_cal=100)
      character*5 ATOM
      character*5 ATMSHELL(max_cal)
      character*15 BASIS
      dimension EXPO(max_cal)
      dimension NUM_DIFF(max_cal)
      data half / 0.5d+00 /
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(/,2X,' THE DIFFUSE PROGRAM'/
     1         2x,' to dertermine exponents of the diffuse functions'/
     2         2x,' for a given basis set'/
     3         2x,'      (Version of 05/14/2008)'/ 
     4         2x,'      (16:30 PM)'/)
    2 FORMAT(/,2X,' PARAMETERS'/
     1 2X,' ATOM       = ',2x,A5/
     2 2X,' BASIS SET  = ',2x,A15/
     3 2X,' NUM_SHELLS = ',I3/
     4 2X,' PRINT      = ',I3/)
    3 FORMAT(A3,2x,I5)
    4 FORMAT(2x,' Atomic Shell      = ',4x,A3)
    5 FORMAT(2x,' Number of diffuse = ',I5)
    6 FORMAT(/,2X,' NO.',5x,' EXPONENT')
    7 FORMAT(3x,'----------------------')
    8 FORMAT(f10.6)
    9 FORMAT(1x,i5,1x,f15.10)
C
      call psinit('APPEND')
      CALL TSTART(6)     
      input=5
      itap67=68
      call ffile(itap68,'file68',0)
c
      write(6,1)
      write(6,*) '  '
      write(6,*) '  maxcor     = ',maxcor
C
      errcod = fwkclr()
      errcod = fwkadd(5,'DIFFUSE')
      errcod = fwkadd(5,'DEFAULT')
      iprnt  = 0
      errcod = frdi('PRINT',iprnt)
      errcod = frdi('NUM_SHELLS',num_shells)
      errcod = frdc('ATOM',atom)
      errcod = frdc('BASIS',basis)
C
      WRITE(6,2) ATOM,BASIS,NUM_SHELLS,IPRNT
c     WRITE(6,2) ATOM,NUM_SHELLS,IPRNT
C
C====================================================
C  Read in two frequencies for each level of theory
C  and calculate the Renner parameters
C====================================================
      write(6,*) ' '
      DO 103 ISHL=1,NUM_SHELLS
      read(itap68,3) ATMSHELL(ISHL),NUM_DIFF(ISHL)
      NUM_DFS=NUM_DIFF(ISHL)
      WRITE(6,4) ATMSHELL(ISHL)
      WRITE(6,5) NUM_DIFF(ISHL)
      WRITE(6,6)     
      WRITE(6,7)     
      ii = 3
      DO 101 I=1,3
      READ(itap68,8) EXPO(I)
      WRITE(6,9) I,EXPO(I)
  101 CONTINUE          
      DO 102 IDIF=1,NUM_DFS
      ii = ii + 1
      iexp    = ii - 3
      expo1   = expo(iexp)
      expo2   = expo(iexp+1)
      expo3   = expo(iexp+2)
      expo32  = (expo3/expo2)
      expo21  = (expo2/expo1)
      expo4   = ( expo32 + expo21 ) * (expo3) * half
      expo(iexp+3) = expo4
      WRITE(6,9) ii,expo4
  102 continue
      write(6,*) ' '
      write(6,*) ' '
  103 CONTINUE
C
      call tstop(6)
      stop
      end
