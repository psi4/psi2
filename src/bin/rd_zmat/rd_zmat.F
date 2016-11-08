      subroutine fentry(cc,ia,maxcor)
c**********************************************************************
C     PROGRAM RD_ZMAT                                                 *
c     This program will update the perturbed geometry                 *
c          in terms of internal coordinates for XYZ type molecules    *
c**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
c**********************************************************************
C   Last updated on May 19, 2003 by Dr. Yukio Yamaguchi               *
c**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension cc(maxcor),ia(maxcor*2)
      character*80 title
      character*80 record(200)
      character*80 gparas(100)
      character*8 wfn,moltyp
      dimension istp(100)
      dimension valu(100)
c==================================================================
c  For psi 2.0 input
#include <error.h>
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(//,2X,' A PROGRAM TO READ INFORMATION FROM ZMAT INPUT'/
     2          2x,'        for a XYZ type molecule  (Version of 05/19/2  
     3003)'/)
    2 FORMAT(//,2X,' PARAMETERS'/
     1 2X,' IPRNT    = ',I5/
     2 2X,' WFN      = ',A8/
     3 2X,' MOLTYP   = ',A8/
     4 2X,' NATOMS   = ',I5/
     5 2X,' N3N      = ',I5/)
    3 FORMAT(A80)
    4 FORMAT(2x,I5,2x,A75)
    5 FORMAT(f20.16)
    6 FORMAT(A7,1x,f20.16)
C                                                                       FDG00570
      INPUT=5                                                           FDG00580
      ITAPE6=6                                                          FDG00590
      ITAP18=18                                                         FDG00600
      ITAP20=20                                                         FDG00600
      ITAP21=21                                                         FDG00600
      ITAP30=30                                                         FDG00610
C                                                                       FDG00620
      call psinit('APPEND')
      CALL TSTART(6)                                                    FDG00550
      call ffile(itap18,'inpnew',0)
      call ffile(itap20,'zmat',0)
      call ffile(itap21,'newzmat',1)
      WRITE(6,1)                                                        FDG00630
C
C======================================================
C   New input
C======================================================
      errcod = fwkclr()
      errcod = fwkadd(5,'RD_ZMAT')
      errcod = fwkadd(5,'DEFAULT')
      errcod = frdi('PRINT',iprnt)
      errcod = frdc('WFN',wfn)
      errcod = frdc('MOLTYP',moltyp)
      write(6,2) iprnt,wfn,moltyp,natoms,n3n
C
      rewind itap20
C  First block
      nlineT = 0
      nlineA = 0
  201 continue
      read(itap20,3) title
      nlineT = nlineT + 1
      if(iprnt.ge.5) then
        write(6,3) title
      end if
      record(nlineT) = title
C::::::::::::::::::::::::123456789012345:::::::
      if(title(1:15).eq.'               ') then
        go to 301
      end if
        nlineA = nlineA + 1
      go to 201
  301 continue
      if(iprnt.ge.5) then
        write(6,*) '  nlineA = ',nlineA
        write(6,*) ' '
      end if
C
C  Second block
      nlineB = 0
  202 continue
      read(itap20,3) title
      nlineT = nlineT + 1
      if(iprnt.ge.5) then
        write(6,3) title
      end if
      record(nlineT) = title
      if(title(1:15).eq.'               ') then
         go to 302
      end if
      nlineB = nlineB + 1
      gparas(nlineB) = title
      istp(nlineB)=index(title(1:),'=')
      go to 202
  302 continue
      if(iprnt.ge.5) then
        write(6,*) '  nlineB = ',nlineB
        write(6,*) ' '
      end if
C
C  Third block
      nlineC = 0
  203 continue
      read(itap20,3) title
      nlineT = nlineT + 1
      if(iprnt.ge.5) then
        write(6,3) title
      end if
      record(nlineT) = title
      if(title(1:15).eq.'               ') then
         go to 303
      end if
      nlineC = nlineC + 1
      go to 203
  303 continue
      if(iprnt.ge.5) then
        write(6,*) '  nlineC = ',nlineC
        write(6,*) ' '
      end if
C
C  Fourth block
      nlineD = 0
  204 continue
      read(itap20,3,end=400) title
      nlineT = nlineT + 1
      if(iprnt.ge.5) then
        write(6,3) title
      end if
      record(nlineT) = title
      if(title(1:15).eq.'               ') then
         go to 304
      end if
      nlineD = nlineD + 1
      go to 204
  304 continue
      if(iprnt.ge.5) then
        write(6,*) '  nlineD = ',nlineD
        write(6,*) ' '
      end if
C
  400 continue
      write(6,*) '  nlineT = ',nlineT
      write(6,*) '  End of ZMAT input'
C
      if(iprnt.ge.3) then
        write(6,*) ' '
        write(6,*) ' ===================='
        write(6,*) '  Rewrite all records'
        write(6,*) ' ===================='
        write(6,*) ' '
        do 101 i=1,nlineT
        write(6,3) record(i)
  101   continue
        write(6,*) ' '
      endif
C
      if(iprnt.ge.5) then
        write(6,*) ' '
        write(6,*) ' :::::::::::::::::::::::::::::::'
        write(6,*) '  Rewrite geometrical parameters'
        write(6,*) ' :::::::::::::::::::::::::::::::'
        write(6,*) ' '
        do 102 i=1,nlineB
        write(6,4) istp(i),gparas(i)
  102   continue
        write(6,*) ' '
      end if
C
      rewind itap18
      npara=0
  205 continue
      read(itap18,5,end=500) temp
      npara=npara+1
      valu(npara)=temp
      go to 205
  500 continue
      rewind itap18
C
      rewind itap21
      do 103 i=1,nlineA+1
      write(itap21,3) record(i)
  103 continue
C
      if(npara.ne.nlineB) then
        write(6,*) ' '
        write(6,*) '   WARNING   WARNING  '
        write(6,*) '   npara .ne. nlineB  '
        write(6,*) '   npara  = ',npara
        write(6,*) '   nlineB = ',nlineB
        write(6,*) '   First nlineB will be used '
        write(6,*) ' '
      end if
      do 104 i=1,nlineB
      title = gparas(i)
      nstp  = istp(i)
      write(6,6) title(1:7),valu(i)
      write(itap21,6) title(1:7),valu(i)
cc    write(6,6) title(1:nstp),valu(i)
cc    write(itap21,6) title(1:nstp),valu(i)
  104 continue
C
      nstart = nlineA + nlineB + 2
      do 105 i=nstart,nlineT
      write(itap21,3) record(i)
  105 continue
      rewind itap21
C
      CALL TSTOP(6)                                                     FDG01730
      STOP                                                              FDG01750
      END                                                               FDG01760
