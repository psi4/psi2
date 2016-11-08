C**********************************************************************
C     PROGRAM GEOMUP                                                  *
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C   Last updated on April 25, 2003 by Dr. Yukio Yamaguchi             *
c   Updated for the LINUX system                                      *
C**********************************************************************
C   Last updated on June 12, 1996 by Dr. Yukio Yamaguchi              *
c   Modification for power2-aix.41                                    *
c   Need "ropen" and "rclose" statements after call "locate"          *
C**********************************************************************
C   Last updated on June 28, 1995 by Dr. Yukio Yamaguchi              *
c   Modification to read in "chg" from file30                         *
C**********************************************************************
C   Last updated on April 04, 1995 by Dr. Yukio Yamaguchi             *
c   Lines for "subroutine fentry" were added                          *
c   Compiled using "gmake" file                                       *
C**********************************************************************
C   THE SOURCE CODE OF GEOMETRY UPDATING PROGRAM                        FDG00070
C   THIS PROGRAM IS ABLE TO TREAT UP TO `50' ATOMS                      FDG00090
      subroutine fentry(cc,ia,maxcor)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension cc(maxcor),ia(maxcor*2)
      dimension i30(200)
      dimension x(50),y(50),z(50)
      dimension chg(50),ichg(100)
      dimension coord(3,50),icoord(300)
      equivalence (chg,ichg)
      equivalence (coord,icoord)
    1 FORMAT(//,2X,' THE GEOMETRY UPDATING PROGRAM (Version of 4/25/2003
     1)'/)
    2 FORMAT(//,2X,' PARAMETERS'/                                       FDG00270
     1 2X,' MPOINT = ',I5/                                              FDG00280
     2 2X,' MCONST = ',I5/                                              FDG00290
     3 2X,' NCALCS = ',I5/                                              FDG00300
     4 2X,' NATOM  = ',I5/                                              FDG00310
     5 2X,' N3N    = ',I5/                                              FDG00320
     6 2X,' IGEOP  = ',I5/)
    3 FORMAT(//,2X,' OLD COORDINATES IN TAPE30'/                        FDG00440
     1 39X,' OLD X IN A.U.'/                                            FDG00450
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               FDG00460
    4 FORMAT(3X,I5,1X,4F16.10)                                          FDG00470
    5 format(3f20.10)
    6 FORMAT(//,2X,' NEW COORDINATES IN TAPE30'/                        FDG00480
     1 39X,' NEW X IN A.U.'/                                            FDG00490
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)               FDG00500
    7 FORMAT(//,2X,' GEOMETRY IN TAPE30 IS UPDATED'/)                   FDG00510
C                                                                       FDG00570
      INPUT=5                                                           FDG00580
      ITAPE6=6                                                          FDG00590
      ITAP30=30                                                         FDG00610
C                                                                       FDG00620
      call psinit('APPEND')
      CALL TSTART(6)                                                    FDG00550
      WRITE(6,1)                                                        FDG00630
C                                                                       FDG00640
C************************************                                   FDG00650
C***UPDATE TAPE30 FOR NEW GEOMETRY***                                   FDG00660
C************************************                                   FDG00670
C                                                                       FDG00680
C   READ IN GEOMETRY FROM TAPE30                                        FDG00690
      CALL RFILE(ITAP30)                                                FDG00700
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              FDG00710
      MPOINT=I30(2)                                                     FDG00720
      MCONST=I30(3)                                                     FDG00730
      NCALCS=I30(5)                                                     FDG00740
      NATOM=I30(19)                                                     FDG00750
      N3N=NATOM*3                                                       FDG00760
      IGEOP=100+MCONST+MPOINT+NCALCS                                    FDG00770
C                                                                       FDG00780
      WRITE(6,2) MPOINT,MCONST,NCALCS,NATOM,N3N,IGEOP
C                                                                       FDG01130
C   READ IN NUCLEAR CHARGE                                              MAS11050
      JUNK=101+MCONST                                                   MAS11020
      CALL WREADW(ITAP30,I30,MPOINT,JUNK,JUNK)                          MAS11030
Cyy   CALL WREADW(ITAP30,CHG,NATOM*2,I30(1),JUNK)
      CALL WREADW(ITAP30,ICHG,NATOM*2,I30(1),JUNK)
C                                                                       FDG01130
C   READ IN THE PRESENT GEOMETRY                                        FDG01260
      CALL WREADW(ITAP30,LOCCAL,1,IGEOP,JUNK)                           FDG01270
      IGEOP=LOCCAL+80                                                   FDG01280
Cyy   CALL WREADW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        FDG01290
      CALL WREADW(ITAP30,ICOORD,N3N*2,IGEOP,JUNK)                        FDG01290
      WRITE(6,3)
      do 101 i=1,natom
      WRITE(6,4) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)
  101 continue
      CALL RCLOSE(ITAP30,3)
c
c   Read in a set of new cartesian coordinates
      CALL LOCATE(INPUT,'# GEOMUP #',IERR)
      do 102 i=1,natom
      read(5,5) x(i),y(i),z(i)
  102 continue
c
      WRITE(6,6)
      do 103 i=1,natom
      COORD(1,I)=X(I)
      COORD(2,I)=Y(I)
      COORD(3,I)=Z(I)
      WRITE(6,4) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)
  103 continue
      call rfile(itap30)
      CALL WWRITW(ITAP30,COORD,N3N*2,IGEOP,JUNK)
      call rclose(itap30,3)
      WRITE(6,7)
c
      CALL TSTOP(6)
      stop
      end
