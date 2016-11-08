      subroutine fentry(cc,ia,maxcor)
c**********************************************************************
c     PROGRAM RENNER                                                  *
c**********************************************************************
c   Last updated on March 19, 2004 by Dr. Yukio Yamaguchi             *
c**********************************************************************
c  This program determines Renner parameters                          *
c   for linear nondegenerate bending mode                             *
c**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension cc(maxcor),ia(maxcor*2)
      parameter (max_cal=100)
      character*20 LEVEL(max_cal), BASIS(max_cal)
      dimension FREQ1(max_cal),FREQ2(max_cal)
      dimension NPTS(max_cal)
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(/,2X,' THE RENNER PROGRAM'/
     1         2x,' to dertermine Renner parameters'/
     2         2x,' for a linear triatomic molecule'/
     3         2x,'      (Version of 03/22/2004)'/ 
     4         2x,'      (10:40 AM)'/)
    2 FORMAT(/,2X,' PARAMETERS'/
     1 2X,' NUM_LEVELS = ',I5/
     2 2X,' IPRNT      = ',I5/)
    3 FORMAT(A20,I5)
    4 FORMAT(2x,' Level of sophistication = ',3x,A20)
    5 FORMAT(2x,' Number of points        = ',3x,I5)
    6 FORMAT(/,2X,' NO.',5x,' FREQ1',5x,' FREQ2',5x,' BASIS SET'/)
    7 FORMAT(f10.1,f10.1,A20)
    8 FORMAT(1x,i5,f10.1,5x,f10.1,5x,A20)
    9 FORMAT(/,' Renner Parameters'/)
   10 FORMAT(/,2X,' NO.',5x,' EPSILON',3x,' OMEGA',5x,' BASIS SET'/)
   11 FORMAT(1x,i5,f10.4,5x,f10.1,5x,A20)
C                                                                       MK100690
      call psinit('APPEND')
      CALL TSTART(6)                                                    NOR01250
      input=5
      itap67=67
      call ffile(itap67,'file67',0)
c
      write(6,1)
      write(6,*) '  '
      write(6,*) '  maxcor     = ',maxcor
C
      errcod = fwkclr()
      errcod = fwkadd(5,'RENNER')
      errcod = fwkadd(5,'DEFAULT')
      iprnt  = 0
      errcod = frdi('PRINT',iprnt)
      errcod = frdi('NUM_LEVELS',num_lvls)
C
      WRITE(6,2) NUM_LVLS,IPRNT
C
C====================================================
C  Read in two frequencies for each level of theory
C====================================================
      write(6,*) ' '
      II=0
      DO 102 ILEV=1,NUM_LVLS                                       
      read(itap67,3) LEVEL(ILEV),NPTS(ILEV)
      NUM_PTS=NPTS(ILEV)
      WRITE(6,4) LEVEL(ILEV)
      WRITE(6,5) NPTS(ILEV)
      WRITE(6,6)                                                        MK100910
      DO 101 I=1,NUM_PTS                                       
      II=II+1
      READ(itap67,7) FREQ1(II),FREQ2(II),BASIS(II)
      WRITE(6,8) I,FREQ1(II),FREQ2(II),BASIS(II)
  101 CONTINUE                                                          MK100950
      write(6,*) ' '
      write(6,*) ' '
  102 CONTINUE
C
C====================================================
C  Calculate the Renner parameters
C====================================================
      WRITE(6,9)                                                        MK100910
      II=0
      DO 104 ILEV=1,NUM_LVLS
      NUM_PTS=NPTS(ILEV)
      WRITE(6,4) LEVEL(ILEV)
      WRITE(6,5) NPTS(ILEV)
      WRITE(6,10) 
      NUM_PTS=NPTS(ILEV)
      DO 103 I=1,NUM_PTS
      II=II+1
      fr1sq = freq1(ii) * freq1(ii)
      fr2sq = freq2(ii) * freq2(ii)
      f1pf2 = fr1sq + fr2sq
      f1mf2 = fr1sq - fr2sq
      epsil = f1mf2/f1pf2
      omega = dsqrt(f1pf2*0.50d+00)
      WRITE(6,11) I,epsil,omega,BASIS(II)
  103 CONTINUE
      write(6,*) '  '
      write(6,*) '  '
  104 CONTINUE
      write(6,*) '  '
C
      call tstop(6)
      stop
      end
