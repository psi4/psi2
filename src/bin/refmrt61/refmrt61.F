      subroutine fentry(cc,ia,maxcor)
c**********************************************************************
C     PROGRAM REFMRT61                                                *
c     This program will construct an analytic potential               *
c       for a linear triatomic molecule                               *
c     Adjusted to the MOLPRO grogram                                  *
c**********************************************************************
C   Last updated on 04-23-2007 by Dr. Yukio Yamaguchi                 *
c**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (max_cal=10000)
      parameter (max_ijk=1000)
      dimension cc(maxcor),ia(maxcor*2)
      dimension RI0(max_cal),RJ0(max_cal),PK0(max_cal)
      dimension ENER0(max_cal)
      character*80 atitle
      character*80 title1,title2
      data a00,one,two / 0.0d+00 , 1.0d+00 , 2.0d+00 /
      data bohr / 0.529177249d+00 /
      data pi / 3.1415926536d+00 /
      data wave / 219474.63d+00 /
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(//,2X,' ================================================ '/
     1          2X,'     THE ANALYTIC POTENTIAL PROGRAM               '/
     2          2X,'     FOR A LINEAR TRIATOMIC MOLECULE              '/
     3          2X,' ------------------------------------------------ '/
     4          2X,'     to determine analytic potential              '/
     5          2X,'     for a linear triatomic molecule              '/
     6          2X,'     adjusted to the MOLPRO program package       '/
     7          2X,' ------------------------------------------------ '/
     8          2x,'  <also with SPF, Legendre, and Sin coordinates>  '/
     9          2x,'             (Version of 04/23/2007)              '/
     a          2X,' ================================================ '/
     b      )
    2 FORMAT(A80)                                                       FDG00200
    3 FORMAT(3x,a80/3x,'(title in file61)')
    4 FORMAT(//,2X,' PARAMETERS'/
     1 2X,' NUM_PTS   = ',I6/
     2 2X,' IPRNT     = ',I6/)
    5 FORMAT(I5)        
    6 format(//,2x,' PERTURBED GEOMETRIES AND ENERGIES'/
     1 4X,4H NO.,6X,4H R_i,14X,4H R_j,14X,4H P_k,19X,7H Energy/
     3 14X,7H (Bohr),11X,7H (Bohr),11X,9H (radian),14X,10H (hartree),/)
c   7 format(f7.1,5x,3f15.10,5x,f20.12)
    7 format(f7.1,3f15.10,f20.12)
    8 format(2x,i5,3x,f15.10,3x,f15.10,3x,f15.10,5x,f20.12)
    9 format(f7.1,5x,3f15.10,5x,f20.12)
C                                                                       FDG00570
      INPUT=5                                                           FDG00580
      ITAPE6=6                                                          FDG00590
      ITAP61=61
      ITAP62=62
C                                                                       FDG00620
      call psinit('APPEND')
      CALL TSTART(6)                                                    FDG00550
      call ffile(itap61,'file61',0)
      call ffile(itap62,'refmrt61',0)
      WRITE(6,1)                                                        FDG00630
      write(6,*) '  maxcor    = ',maxcor
      write(6,*) '  max_cal   = ',max_cal
      write(6,*) '  max_ijk   = ',max_ijk
C
      errcod    = fwkclr()
      errcod    = fwkadd(5,'REFMRT61')
      errcod    = fwkadd(5,'DEFAULT')
      errcod    = frdi('PRINT',iprnt)
C
C  Read in information from file61
      REWIND ITAP61                                                     FDG00800
      read(itap61,2) atitle
      write(6,*) ' '
      write(6,3) atitle
      read(itap61,2) atitle
      read(itap61,2) atitle
      npoints = 0
  201 continue
      read(itap61,2,end=202) atitle
      npoints = npoints + 1
      go to 201
  202 continue
C
      num_pts = npoints
      WRITE(6,4) IPRNT,NUM_PTS
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C   Read in perturbed geometrical parameters and energies 
C      from "file61.dat" file
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      REWIND ITAP61                                                     FDG00800
      READ(ITAP61,2) title1
      READ(ITAP61,2) title2
      READ(ITAP61,5) num_temp
C
      if(num_temp.ne.num_pts) then
        write(6,*) '  '
        write(6,*) '  *********************************'
        write(6,*) '  NUM_TEMP not equal to NUM_PTS !!!'
        write(6,*) '  NUM_TEPM = ',num_temp
        write(6,*) '  NUM_PTS  = ',num_pts
        write(6,*) '  ---> STOP !!!'
        write(6,*) '  *********************************'
        write(6,*) '  '
        stop
      end if
C
      write(6,6)
      do 101 i=1,num_pts
      read(itap61,7) xnn,rri,rrj,ppk,evalu
      nn = int(xnn)
      write(6,8) nn,rri,rrj,ppk,evalu
      ri0(i)=rri
      rj0(i)=rrj
      pk0(i)=ppk
      ener0(i)=evalu
  101 continue
      write(6,*) '  '
      REWIND ITAP61          
C
      REWIND ITAP62
      WRITE(ITAP62,2) title1
      WRITE(ITAP62,2) title2
      WRITE(ITAP62,5) num_temp
      do 102 i=1,num_pts
      xnn = i
      rri   = ri0(i)
      rrj   = rj0(i)
      ppk   = pk0(i) * 0.50d+00
      evalu = ener0(i)
      write(itap62,9) xnn,rri,rrj,ppk,evalu
  102 continue
      REWIND ITAP62
C
      CALL TSTOP(6)                                                     FDG01730
      STOP                                                              FDG01750
      END                                                               FDG01760
