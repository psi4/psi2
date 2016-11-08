      subroutine fentry(cc,ia,maxcor)
c**********************************************************************
c     PROGRAM DIFFUSE                                                 *
c**********************************************************************
c   Last updated on May 22, 2008 by Dr. Yukio Yamaguchi               *
c**********************************************************************
c  This program computes Rydberg series                               *
c   for a given atom                                                  *
c**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension cc(maxcor),ia(maxcor*2)
      parameter (max_cal=100)
      character*5 ATOM
      character*1 signal
      data rconst / 109737.0d+00 /
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(/,2X,' THE RYDBERG PROGRAM'/
     1         2x,' to dertermine excitation energies'/
     2         2x,' of Rydberg states for a given atom'/
     3         2x,'      (Version of 05/23/2008)'/ 
     4         2x,'      (12:15 PM)'/)
    2 FORMAT(/,2X,' PARAMETERS'/
     1 2X,' ATOM       = ',2x,A5/
     2 2X,' PRINT      = ',I3/ 
     3 2X,' NQ_START   = ',I3/ 
     4 2X,' NQ_END     = ',I3/ 
     5 2X,' NUM_DEL    = ',I3/ 
     6 2X,' I.P.       = ',f12.3/ 
     7 2X,' RCONSTANT  = ',f12.3/ 
     8 2X,' DEL_INC    = ',f12.3/)
    3 FORMAT(A3,2x,I5)
    4 FORMAT(2x,' Principal Quantum Number = ',4x,I3)
    5 FORMAT(/,2X,' No.',5x,' n',5X,' delta',6X,' n-delta',4X,
     1            ' (n-delta)^2',4X,' term[R/(n-delta)^2]',3X,
     2            ' energy[I.P.-term]')
    6 FORMAT(3x,'-------------------------------------------------------
     1------------------------------------------')
    7 FORMAT(2x,I3,5X,I3,1X,f12.5,f12.5,f12.5,6X,f12.3,11X,f12.3,1X,a1)
C
      call psinit('APPEND')
      CALL TSTART(6)     
      input=5
c
      write(6,1)
      write(6,*) '  '
      write(6,*) '  maxcor     = ',maxcor
      write(6,*) '  '
C
      errcod = fwkclr()
      errcod = fwkadd(5,'RYDBERG')
      errcod = fwkadd(5,'DEFAULT')
      iprnt  = 0
      errcod = frdi('PRINT',iprnt)
      errcod = frdi('NQ_START',nq_start)
      errcod = frdi('NQ_END',nq_end)
      errcod = frdi('NUM_DEL',num_del)
      errcod = frdc('ATOM',atom)
      errcod = frdf('IP',potent)
      delinc = 0.01d+00
      errcod = frdf('DEL_INC',delinc)
C
      WRITE(6,2) ATOM,IPRNT,nq_start,nq_end,num_del,potent,rconst,delinc
C
C====================================================
C  Read in two frequencies for each level of theory
C  and calculate the Renner parameters
C====================================================
      write(6,*) ' '
      DO 102 INQ=NQ_START,NQ_END
      WRITE(6,4) INQ
      WRITE(6,5)
      WRITE(6,6)
      rnq = inq
      DO 101 IDEL=1,NUM_DEL+1
      rdel = idel
      delta = (rdel - 1) * delinc
      diff = rnq - delta
      if(diff .le. 0.0d+00) go to 101
      diffsq = diff * diff
      term = (rconst) / (diffsq)
      if(term .gt. rconst) go to 101
      signal = ' '
      if(term .gt. potent) then
        signal = '*'
      end if
      ex_ener = potent - term
      WRITE(6,7) IDEL,INQ,delta,diff,diffsq,term,ex_ener,signal
  101 continue
      write(6,*) ' '
      write(6,*) ' '
  102 CONTINUE
C
      call tstop(6)
      stop
      end
