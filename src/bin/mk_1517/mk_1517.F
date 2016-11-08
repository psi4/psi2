      subroutine fentry(cc,ia,maxcor)
ccc   PROGRAM MK_1517                                               
C**********************************************************************
C*  Last updated on June 02, 2003 by Dr. Yukio Yamaguchi              *
C*      LINUX system                                                  *
C**********************************************************************
C*  Last updated on October 03, 2000 by Dr. Yukio Yamaguchi           *
C*    Modification for all type of wave functions                     *
C*      and for a simplified input                                    *
C=====================================================================*
C  (1)  print = a number                                              *
C     A larger number porvides more output                            *
C  (2)  dipole = YES or NO                                            *
C     YES means dipole moment derivatives will be calculated          *
C     NO  means dipole moment derivatives will not be calculated      *
C     Default is YES                                                  *
C  (3)  wfn = wave function type                                      *
C     SCF and CI wave functions require file11                        *
C     CCSD and CCSD(T) wave functions require file11 and file28       *
C     Default is WFN given in default section                         *
C  (4)  geom_up = YES or NO                                           *
C     Update file30 with the optimized geometry given in file11       *
C     Default is YES                                                  *
C=====================================================================*
C**********************************************************************
C*  Last updated on January 31, 2000 by Dr. Yukio Yamaguchi           *
C*  Modification for the CCSD and CCSD(T)                             *
C*    Closed-Shell Wave Functions                                     *
C*    to use file28 for dipole moment derivatives                     *
C**********************************************************************
C*  Last updated on October 16, 1991 by Dr. Yukio Yamaguchi           *
C*  Modification for the Unix Version                                 *
C********************************************************************** MK100020
C*  LAST UPDATED ON DECEMBER 23, 1990 BY YUKIO YAMAGUCHI              * MK100030
C********************************************************************** MK100040
C   FORM THE HESSIAN (FILE15) AND DIPOLE MOMENT DERIVATIVE (FILE17)   * MK100050
C   FILES FROM THE CUMMULATED FILE11                                  * MK100060
C   CARTESIAN COORDINATE SYSTEM                                       * MK100070
C********************************************************************** MK100080
      IMPLICIT REAL*8 (A-H,O-Z)                                         MK100090
      dimension cc(maxcor),ia(maxcor*2)
      DIMENSION FF(150,150),AA(150,150)                                 MK100100
      DIMENSION XX(150),YY(150),ZZ(150)                                 MK100110
      CHARACTER*80 ATITLE                                               MK100120
      COMMON/COM101/NATOM,N3N                                           MK100130
      COMMON/COM102/CHG(150),X(150),Y(150),Z(150)                       MK100140
      COMMON/COM103/R(1275)                                             MK100150
      COMMON/COM104/NT,NCASE,ISPOT(150)                                 MK100160
      COMMON/COM105/IOFF(256),IPRNT                                     MK100170
      COMMON/COM106/PTR(3,3,48)                                         MK100180
      COMMON/COM107/ICT(150,48)                                         MK100190
      COMMON/COM201/ENGX(150),ENGY(150),ENGZ(150)                       MK100200
      COMMON/COM202/GRADX(150),GRADY(150),GRADZ(150)                    MK100210
      COMMON/COM203/FP(150,150),FM(150,150),FX(150,150)                 MK100220
      COMMON/COM204/DP(3,150),DM(3,150),DIPDT(3,150)                    MK100230
      DATA A00,ONE,TWO / 0.0D+00 , 1.0D+00 , 2.0D+00 /                  MK100240
      DATA ITAP11,ITAP15,ITAP17,ITAP28 / 11 , 15 , 17 , 28 /   
      DATA BOHR / 0.52917706D+00 /                                      MK100260
      DATA XLIMIT / 1.0D-6 /                                            MK100270
      DATA DIFLIM / 1.0D-5 /                                            MK100270
c
c==================================================================
c  For psi 2.0 input
#include <error.h>
      character*6 wfn, dipole, geomup
      integer frdc, frdi, frdf, fwkclr, fwkadd, frdboo, frdpos
      integer fcnt, frdia1
      integer errcod
c==================================================================
    1 FORMAT(//,2X,' THE PROGRAM TO FORM FILE15 AND FILE17'/            MK100280
     1          2X,' IN THE CARTESIAN COORDINATE SYSTEM'/     
     2          2X,' (VERSION OF 06/14/2001)'/)                         MK100300
c    2          2X,' (VERSION OF 06/02/2003)'/)                         MK100300
    2 FORMAT(5I5)                                                       MK100310
    3 FORMAT(A80)                                                       MK100320
    4 FORMAT(3X,A80)                                                    MK100330
    5 FORMAT(I5,F20.10,F15.10,2F20.10)                                  MK100340
    6 FORMAT(/,2X,' PARAMETERS FOR THE CALCULATION'//                   MK100350
     * 2X,' IPRNT  = ',I5/                                              MK100360
     * 2X,' WFN    = ',A6/ 
     * 2X,' DIPOLE = ',A6/ 
     * 2X,' GEOMUP = ',A6/ 
     * 2X,' NATOM  = ',I5/    
     * 2X,' N3N    = ',I5/   
     * 2X,' ENERGY = ',F20.10/)
    7 FORMAT(//,2X,' REFERENCE GEOMETRY'/                               MK100400
     1 5X,4H NO.,11X,7H CHARGE,13X,2H X,18X,2H Y,18X,2H Z/)             MK100410
    8 FORMAT(4F20.10)                                                   MK100420
    9 FORMAT(2X,I7,5X,4F20.10)                                          MK100430
   10 FORMAT(/,2X,' DIPOLE MOMENTS',17X,3F20.10/)                       MK100440
   11 FORMAT(//,2X,' CARTESIAN ENERGY GRADIENTS'/                       MK100450
     1 5X,4H NO.,11X,7H CHARGE,13X,2H X,18X,2H Y,18X,2H Z/)             MK100460
   12 FORMAT(20X,3F20.10)                                               MK100470
   13 FORMAT(/,2X,' NCASE = ',I5/)                                      MK100480
   14 FORMAT(//,2X,' PERTURBED GEOMETRY'/                               MK100490
     1 5X,4H NO.,11X,7H CHARGE,13X,2H X,18X,2H Y,18X,2H Z/)             MK100500
   15 FORMAT(//,2X,' NUMBER OF GRADIENT AND DIPOLE MOMENT SETS',I5/)    MK100510
   16 FORMAT(//,2X,' DP MATRIX BEFORE SYMMETRIZATION'/)                 MK100520
   17 FORMAT(//,2X,' DM MATRIX BEFORE SYMMETRIZATION'/)                 MK100530
   18 FORMAT(//,2X,' FP MATRIX BEFORE SYMMETRIZATION'/)                 MK100540
   19 FORMAT(//,2X,' FM MATRIX BEFORE SYMMETRIZATION'/)                 MK100550
   20 FORMAT(//,2X,' DELX = ',F20.10/)                                  MK100560
   21 FORMAT(//,3X,52('*')                                              MK100570
     *        /,2X,' DIFFERENT DELX IS USED !   NO. OF OCCURENCES = ',I5MK100580
     *        /,3X,52('*'),/)                                           MK100590
   22 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS (IN DEBYE/BOHR)MK100600
     *'/)                                                               MK100610
   23 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS (IN DEBYE/A)'/)MK100620
   24 FORMAT(//,2X,' FX MATRIX (IN A.U.)'/)                             MK100630
   25 FORMAT(//,2X,' AVERAGED FX MATRIX (UPPER TRI. = FX, LOWER TRI. = DMK100640
     1IFFERENCE)'/)                                                     MK100650
   26 FORMAT(//,2X,' THE LARGEST DIFFERENCE IN DIAGONAL ELEMENTS IS = ',
     1 2x,F15.10/)
   27 FORMAT(//,2X,' FX MATRIX (IN HARTREE/BOHR+2)'/)                   MK100660
   28 FORMAT(2I5)                                                       MK100670
   29 FORMAT(3F20.10)                                                   MK100680
C                                                                       MK100690
      call psinit('APPEND')
      CALL TSTART(6)                                                    MK100700
      INPUT=5                                                           MK100720
      ITAPE6=6                                                          MK100730
      NDIM=150                                                          MK100740
      NDIM3=NDIM*3                                                      MK100750
      NDIMSQ=NDIM*NDIM                                                  MK100760
      RBOHR=ONE/BOHR                                                    MK100770
C                                                                       MK100780
      call ffile(itap11,'file11',0)
      call ffile(itap15,'file15',0)
C                                                                       MK100800
      WRITE(6,1)                                                        MK100810
      write(6,*) '  '
      write(6,*) '  maxcor = ',maxcor
C======================================================
C   Old input
C======================================================
c     READ(5,2) IPRNT      
C======================================================
C   New input
C======================================================
      errcod = fwkclr()
      errcod = fwkadd(5,'MK_1517')
      errcod = fwkadd(5,'DEFAULT')
      errcod = frdi('PRINT',iprnt)
      errcod = frdc('WFN',wfn)
      dipole = 'YES   '
      errcod = frdc('DIPOLE',dipole)
      if(dipole.eq.'YES   ') then
        call ffile(itap17,'file17',0)
        if(wfn.eq.'CCSD  '.or.wfn.eq.'CCSD_T') then
          call ffile(itap28,'file28',0)
        end if
      end if
      geomup = 'YES   '
      errcod = frdc('GEOM_UP',geomup)
C======================================================
C                                                                       MK100830
C   READ IN THE REFERENCE GEOMETRY                                      MK100840
      REWIND ITAP11                                                     MK100850
      if(dipole.eq.'YES   ') then
        if(wfn.eq.'CCSD  '.or.wfn.eq.'CCSD_T') then
          rewind itap28
        end if
      end if
c
      READ(ITAP11,3) ATITLE                                             MK100860
      WRITE(6,4) ATITLE                                                 MK100870
      READ(ITAP11,5) NATOM,ENERGY,dx,dy,dz
      N3N=NATOM*3                                                       MK100890
      WRITE(6,6) iprnt,wfn,dipole,geomup,
     *           NATOM,N3N,ENERGY
      WRITE(6,7)                                                        MK100910
      DO 101 I=1,NATOM                                                  MK100920
      READ(ITAP11,8) CHG(I),X(I),Y(I),Z(I)                              MK100930
      WRITE(6,9) I,CHG(I),X(I),Y(I),Z(I)                                MK100940
  101 CONTINUE                                                          MK100950
      if(dipole.eq.'YES   ') then
        if(wfn.eq.'CCSD  '.or.wfn.eq.'CCSD_T') then
          read(itap28,29) dx,dy,dz
        end if
      end if
      WRITE(6,10) DX,DY,DZ
      WRITE(6,11)                                                       MK100970
      DO 102 I=1,NATOM                                                  MK100980
      READ(ITAP11,12) ENGX(I),ENGY(I),ENGZ(I)                           MK100990
      WRITE(6,9) I,CHG(I),ENGX(I),ENGY(I),ENGZ(I)                       MK101000
  102 CONTINUE                                                          MK101010
      CALL DIST                                                         MK101020
C                                                                       MK101030
C   READ IN THE GRADIENT SET WITH PERTURBATIONS                         MK101040
      CALL ZERO(DP,NDIM3)                                               MK101050
      CALL ZERO(DM,NDIM3)                                               MK101060
      CALL ZERO(FP,NDIMSQ)                                              MK101070
      CALL ZERO(FM,NDIMSQ)                                              MK101080
C                                                                       MK101090
      IDELTA=0                                                          MK101100
      NCASE=0                                                           MK101110
  200 NCASE=NCASE+1                                                     MK101120
      IF(IPRNT.GT.3) THEN                                               MK101130
        WRITE(6,13) NCASE                                               MK101140
      END IF                                                            MK101150
      READ(ITAP11,3,END=205) ATITLE                                     MK101160
      IF(IPRNT.GT.3) THEN                                               MK101170
        WRITE(6,4) ATITLE                                               MK101180
      END IF                                                            MK101190
      READ(ITAP11,5) NATOM,ENERGY,dx,dy,dz
      if(dipole.eq.'YES   ') then
        if(wfn.eq.'CCSD  '.or.wfn.eq.'CCSD_T') then
          read(itap28,29) dx,dy,dz
        end if
      end if
      IF(IPRNT.GT.3) THEN                                               MK101210
        WRITE(6,14)                                                     MK101220
      END IF                                                            MK101230
      DO 103 I=1,NATOM                                                  MK101240
      READ(ITAP11,8) CHG(I),XX(I),YY(I),ZZ(I)                           MK101250
      IF(IPRNT.LE.2) GO TO 201                                          MK101260
      WRITE(6,9) I,CHG(I),XX(I),YY(I),ZZ(I)                             MK101270
  201 CONTINUE                                                          MK101280
      IADD=(I-1)*3                                                      MK101290
      DIFX=XX(I)-X(I)                                                   MK101300
      DIFY=YY(I)-Y(I)                                                   MK101310
      DIFZ=ZZ(I)-Z(I)                                                   MK101320
      IF(DABS(DIFX).GT.XLIMIT) GO TO 301                                MK101330
      IF(DABS(DIFY).GT.XLIMIT) GO TO 302                                MK101340
      IF(DABS(DIFZ).GT.XLIMIT) GO TO 303                                MK101350
      GO TO 103                                                         MK101360
C                                                                       MK101370
C  X PERTURBATION                                                       MK101380
  301 CONTINUE                                                          MK101390
      IPOS=IADD+1                                                       MK101400
      IF(DIFX.LT.A00) IPOS=IPOS+N3N                                     MK101410
      IF(NCASE.EQ.1) DELX=DABS(DIFX)                                    MK101420
      XDIF=DELX-DABS(DIFX)                                              MK101430
      IF(DABS(XDIF).GT.XLIMIT) IDELTA=IDELTA+1                          MK101440
      GO TO 103                                                         MK101450
C                                                                       MK101460
C  Y PERTURBATION                                                       MK101470
  302 CONTINUE                                                          MK101480
      IPOS=IADD+2                                                       MK101490
      IF(DIFY.LE.A00) IPOS=IPOS+N3N                                     MK101500
      IF(NCASE.EQ.1) DELX=DABS(DIFY)                                    MK101510
      XDIF=DELX-DABS(DIFY)                                              MK101520
      IF(DABS(XDIF).GT.XLIMIT) IDELTA=IDELTA+1                          MK101530
      GO TO 103                                                         MK101540
C                                                                       MK101550
C  Z PERTURBATION                                                       MK101560
  303 CONTINUE                                                          MK101570
      IPOS=IADD+3                                                       MK101580
      IF(DIFZ.LE.A00) IPOS=IPOS+N3N                                     MK101590
      IF(NCASE.EQ.1) DELX=DABS(DIFZ)                                    MK101600
      XDIF=DELX-DABS(DIFZ)                                              MK101610
      IF(DABS(XDIF).GT.XLIMIT) IDELTA=IDELTA+1                          MK101620
  103 CONTINUE                                                          MK101630
      IF(IPRNT.GT.3) THEN                                               MK101640
        WRITE(6,10) DX,DY,DZ                                            MK101650
        WRITE(6,11)                                                     MK101660
      END IF                                                            MK101670
      DO 104 I=1,NATOM                                                  MK101680
      READ(ITAP11,12) GRADX(I),GRADY(I),GRADZ(I)                        MK101690
      IF(IPRNT.LE.2) GO TO 104                                          MK101700
      WRITE(6,9) I,CHG(I),GRADX(I),GRADY(I),GRADZ(I)                    MK101710
  104 CONTINUE                                                          MK101720
      ISPOT(NCASE)=IPOS                                                 MK101730
      IF(IPOS.GT.N3N) GO TO 202                                         MK101740
C                                                                       MK101750
C   POSITIVE PERTURBATION                                               MK101760
      DP(1,IPOS)=DX                                                     MK101770
      DP(2,IPOS)=DY                                                     MK101780
      DP(3,IPOS)=DZ                                                     MK101790
      DO 105 I=1,NATOM                                                  MK101800
      IX=(I-1)*3+1                                                      MK101810
      IY=IX+1                                                           MK101820
      IZ=IY+1                                                           MK101830
      FP(IPOS,IX)=GRADX(I)                                              MK101840
      FP(IPOS,IY)=GRADY(I)                                              MK101850
      FP(IPOS,IZ)=GRADZ(I)                                              MK101860
  105 CONTINUE                                                          MK101870
      GO TO 203                                                         MK101880
C                                                                       MK101890
C   NEGATIVE PERTURBATION                                               MK101900
  202 CONTINUE                                                          MK101910
      DM(1,IPOS-N3N)=DX                                                 MK101920
      DM(2,IPOS-N3N)=DY                                                 MK101930
      DM(3,IPOS-N3N)=DZ                                                 MK101940
      DO 106 I=1,NATOM                                                  MK101950
      IX=(I-1)*3+1                                                      MK101960
      IY=IX+1                                                           MK101970
      IZ=IY+1                                                           MK101980
      FM(IPOS-N3N,IX)=GRADX(I)                                          MK101990
      FM(IPOS-N3N,IY)=GRADY(I)                                          MK102000
      FM(IPOS-N3N,IZ)=GRADZ(I)                                          MK102010
  106 CONTINUE                                                          MK102020
  203 CONTINUE                                                          MK102030
      GO TO 200                                                         MK102040
C                                                                       MK102050
C   END OF GRADIENT SETS                                                MK102060
  205 CONTINUE                                                          MK102070
      rewind itap11
      if(dipole.eq.'YES   ') then
        if(wfn.eq.'CCSD  '.or.wfn.eq.'CCSD_T') then
          rewind itap28
        end if
      end if
      NCASE=NCASE-1                                                     MK102080
      WRITE(6,15) NCASE                                                 MK102090
      if(dipole.eq.'YES   ') then
        WRITE(6,16) 
        CALL MATOUT(DP,3,NDIM,3,N3N,6)    
        WRITE(6,17)                      
        CALL MATOUT(DM,3,NDIM,3,N3N,6)  
      end if
      WRITE(6,18)                                                       MK102140
      CALL MATOUT(FP,NDIM,NDIM,N3N,N3N,6)                               MK102150
      WRITE(6,19)                                                       MK102160
      CALL MATOUT(FM,NDIM,NDIM,N3N,N3N,6)                               MK102170
      WRITE(6,20) DELX                                                  MK102180
      IF(IDELTA.GE.1) THEN                                              MK102190
        WRITE(6,21) IDELTA                                              MK102200
      END IF                                                            MK102210
C
C   Update geometry in file30 using the optimized geometry in file11
C     if required
      if(geomup.eq.'YES   ') then
          CALL UPDATE
      end if
C                                                                       MK102220
C   GET SYMMETRY INFORMATION FROM TAPE30                                MK102230
      CALL GETBAS                                                       MK102240
C                                                                       MK102250
C   SYMMETRIZE THE GRADIENT AND DIPOLE MOMENT MATRICES                  MK102260
      IF(NT.EQ.1) GO TO 206                                             MK102270
      if(dipole.eq.'YES   ') then
        CALL SYMDP(DP,DM,NDIM)
      end if
      CALL SYMFX(FP,FM,NDIM)                                            MK102290
C                                                                       MK102300
C   CALCULATE THE DIPOLE DERIVATIVE MATRIX                              MK102310
  206 CONTINUE                                                          MK102320
      FACT=DELX*TWO                                                     MK102330
      if(dipole.eq.'NO    ') go to 207
      DO 107 I=1,3                                                      MK102340
      DO 107 J=1,N3N                                                    MK102350
      DIPDT(I,J)=(DP(I,J)-DM(I,J))/FACT                                 MK102360
  107 CONTINUE                                                          MK102370
      WRITE(6,22)                                                       MK102380
      CALL MATOUT(DIPDT,3,NDIM,3,N3N,6)                                 MK102390
      DO 108 I=1,3                                                      MK102400
      DO 108 J=1,N3N                                                    MK102410
      DIPDT(I,J)=DIPDT(I,J)*RBOHR                                       MK102420
  108 CONTINUE                                                          MK102430
      WRITE(6,23)                                                       MK102440
      CALL MATOUT(DIPDT,3,NDIM,3,N3N,6)                                 MK102450
C                                                                       MK102460
C   CALCULATE THE FORCE CONSTANT MATRIX                                 MK102470
  207 continue
      DO 109 I=1,N3N                                                    MK102480
      DO 109 J=1,N3N                                                    MK102490
      AA(I,J)=(FP(I,J)-FM(I,J))/FACT                                    MK102500
  109 CONTINUE                                                          MK102510
      WRITE(6,24)                                                       MK102520
      CALL MATOUT(AA,NDIM,NDIM,N3N,N3N,6)                               MK102530
      DO 110 I=1,N3N                                                    MK102540
      DO 110 J=I,N3N                                                    MK102550
      FF(I,J)=(AA(I,J)+AA(J,I))/TWO                                     MK102560
      IF(I.EQ.J) GO TO 110                                              MK102570
      FF(J,I)=(AA(I,J)-AA(J,I))/TWO                                     MK102580
  110 CONTINUE                                                          MK102590
      WRITE(6,25)                                                       MK102600
      CALL MATOUT(FF,NDIM,NDIM,N3N,N3N,6)                               MK102610
C
C  Find the largest difference in off-diagonal elements
      difmax = 0.0
      DO 111 I=1,N3N                                                    MK102620
      DO 111 J=I,N3N                                                    MK102630
      IF(I.EQ.J) GO TO 111                                              MK102650
      temp = FF(J,I)
      if(dabs(temp).gt.difmax) then
        difmax = dabs(temp)
      end if
  111 CONTINUE                                        
      write(6,*) '  '
      write(6,26) difmax
      if(difmax.gt.diflim) then
        write(6,*) '  '
        write(6,*) '  ******************************************'
        write(6,*) '  WARNING!         WARNING!         WARNING!'
        write(6,*) '  Large difference in off-diagonal elements '
        write(6,*) '  ***************************************** '
      end if
C
C  Symmetrized the FX matrix
      DO 112 I=1,N3N                                                    MK102620
      DO 112 J=I,N3N                                                    MK102630
      FX(I,J)=FF(I,J)                                                   MK102640
      IF(I.EQ.J) GO TO 112                                              MK102650
      FX(J,I)=FF(I,J)                                                   MK102660
  112 CONTINUE                                                          MK102670
      WRITE(6,27)                                                       MK102680
      CALL MATOUT(FX,NDIM,NDIM,N3N,N3N,6)                               MK102690
      WRITE(6,*) '  '
C                                                                       MK102700
C   FORM TAPE15                                                         MK102710
      REWIND ITAP15                                                     MK102720
      WRITE(ITAP15,28) NATOM,N3N*2                                      MK102730
      WRITE(ITAP15,29) ((FX(I,J),J=1,N3N),I=1,N3N)                      MK102740
      REWIND ITAP15                                                     MK102750
C                                                                       MK102760
C   FORM TAPE17                                                         MK102770
      if(dipole.eq.'YES   ') then
        REWIND ITAP17  
        WRITE(ITAP17,28) NATOM,N3N 
        WRITE(ITAP17,29) ((DIPDT(I,J),J=1,N3N),I=1,3) 
        REWIND ITAP17                                                     MK102810
      end if
C                                                                       MK102820
      CALL TSTOP(6)                                                     MK102830
      STOP                                                              MK102840
      END                                                               MK102850
      SUBROUTINE DIST                                                   MK102860
      IMPLICIT REAL*8 (A-H,O-Z)                                         MK102870
      COMMON/COM101/NATOM,N3N                                           MK102880
      COMMON/COM102/CHG(150),X(150),Y(150),Z(150)                       MK102890
      COMMON/COM103/R(1275)                                             MK102900
    1 FORMAT(//,2X,' INTERATOMIC DISTANCE MATRIX'/)                     MK102910
C                                                                       MK102920
      CALL ZERO(R,1275)                                                 MK102930
      IJ=0                                                              MK102940
      DO 101 I=1,NATOM                                                  MK102950
      DO 101 J=1,I                                                      MK102960
      IJ=IJ+1                                                           MK102970
      IF(I.EQ.J) GO TO 101                                              MK102980
      XD=X(I)-X(J)                                                      MK102990
      YD=Y(I)-Y(J)                                                      MK103000
      ZD=Z(I)-Z(J)                                                      MK103010
      R(IJ)=DSQRT(XD*XD+YD*YD+ZD*ZD)                                    MK103020
  101 CONTINUE                                                          MK103030
      WRITE(6,1)                                                        MK103040
      CALL PRINT(R,1275,NATOM,6)                                        MK103050
      RETURN                                                            MK103060
      END                                                               MK103070
      SUBROUTINE UPDATE
      IMPLICIT REAL*8 (A-H,O-Z)   
      COMMON/COM101/NATOM,N3N
      COMMON/COM102/CHG(150),X(150),Y(150),Z(150)
      DIMENSION COORD(3,50)
      DIMENSION I30(200)
    1 FORMAT(//,2X,' THE OPTIMIZED GEOMETRY IN TAPE11'/      
     1 39X,' OPT. X IN A.U.'/                              
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)  
    2 FORMAT(3X,I5,1X,4F16.10)                             
    3 FORMAT(//,2X,' OLD COORDINATES IN TAPE30'/          
     1 39X,' OLD X IN A.U.'/                                
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)    
    4 FORMAT(//,2X,' NEW COORDINATES IN TAPE30'/            
     1 39X,' NEW X IN A.U.'/                                 
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z/)  
    5 FORMAT(//,2X,' GEOMETRY IN TAPE30 IS UPDATED'/) 
C                                                         
C   PRINT OUT THE OPTIMIZED GEOMETRY IN FILE11 
      WRITE(6,1) 
      DO 101 I=1,NATOM                        
      WRITE(6,2) I,CHG(I),X(I),Y(I),Z(I)   
  101 CONTINUE                               
C
C   READ IN THE PRESENT GEOMETRY FROM TAPE30
      itap30 = 30
      CALL RFILE(ITAP30)                                                FDG00700
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              FDG00710
      MPOINT=I30(2)                                                     FDG00720
      MCONST=I30(3)                                                     FDG00730
      NCALCS=I30(5)                                                     FDG00740
      NATOM=I30(19)                                                     FDG00750
      N3N=NATOM*3                                                       FDG00760
      IGEOP=100+MCONST+MPOINT+NCALCS                                    FDG00770
      CALL WREADW(ITAP30,LOCCAL,1,IGEOP,JUNK)                           FDG01270
      IGEOP=LOCCAL+80                                                   FDG01280
      CALL WREADW(ITAP30,COORD,N3N*2,IGEOP,JUNK)                        FDG01290
      WRITE(6,3) 
      DO 102 I=1,NATOM                                                  FDG01310
      WRITE(6,2) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)
  102 CONTINUE                                                          FDG01330
C                                                                       FDG01340
C   UPDATE THE GEOMETRY IN TAPE30                                       FDG01350
      WRITE(6,4) 
      DO 103 I=1,NATOM                                                  FDG01510
      IX=I*3-2                                                          FDG01520
      IY=I*3-1                                                          FDG01530
      IZ=I*3                                                            FDG01540
      COORD(1,I)=X(I)
      COORD(2,I)=Y(I)
      COORD(3,I)=Z(I)
      WRITE(6,2) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)
  103 CONTINUE   
      CALL WWRITW(ITAP30,COORD,N3N*2,IGEOP,JUNK)    
      WRITE(6,5)                                                 
      CALL RCLOSE(ITAP30,3)
C
      RETURN
      END
      SUBROUTINE GETBAS                                                 MK103080
      IMPLICIT REAL*8 (A-H,O-Z)                                         MK103090
      DIMENSION I30(200),LICT(400)                                      MK103100
      COMMON/COM101/NATOM,N3N                                           MK103110
      COMMON/COM104/NT,NCASE,ISPOT(150)                                 MK103120
      COMMON/COM105/IOFF(256),IPRNT                                     MK103130
      COMMON/COM106/PTR(3,3,48)                                         MK103140
      COMMON/COM107/ICT(150,48)                                         MK103150
C                                                                       MK103160
    1 FORMAT(/,2X,' NT               = ',I5/)   
    2 FORMAT(/,2X,' ICT MATRIX, ISYM = ',I5/)  
    3 FORMAT(2X,10I5)                         
    4 FORMAT(/,2X,' PTR MATRIX, ISYM = ',I5/)
C                                                                       MK103200
      IOFF(1)=0                                                         MK103210
      DO 101 I=1,255                                                    MK103220
      IOFF(I+1)=IOFF(I)+I                                               MK103230
  101 CONTINUE                                                          MK103240
C                                                                       MK103250
C   GET CONSTANTS FROM TAPE30                                           MK103260
      ITAP30=30                                                         MK103270
      CALL RFILE(ITAP30)                                                MK103280
      CALL SREW(ITAP30)                                                 MK103290
      CALL WREADW(ITAP30,I30,200,101,JUNK)                              MK103300
C                                                                       MK103310
      MPOINT  = I30(2)                                                  MK103320
      MCONST  = I30(3)                                                  MK103330
      NATOM   = I30(19)                                                 MK103340
      NT      = I30(29)                                                 MK103350
      N3N=NATOM*3                                                       MK103360
      WRITE(6,1) NT
C                                                                       MK103380
C   READ POINTERS FROM TAPE30                                           MK103390
      IPOS=101+MCONST                                                   MK103400
      CALL WREADW(ITAP30,I30,MPOINT,IPOS,JUNK)                          MK103410
C                                                                       MK103420
C   READ SYMMETRY INFORMATION                                           MK103430
      CALL WREADW(ITAP30,LICT,NATOM*NT,I30(2),JUNK)                     MK103440
      CALL WREADW(ITAP30,PTR,3*3*NT*2,I30(31),JUNK)                     MK103450
C                                                                       MK103460
      DO 103 I=1,NT                                                     MK103470
      II=(I-1)*NATOM                                                    MK103480
      DO 102 IATOM=1,NATOM                                              MK103490
      II=II+1                                                           MK103500
      ICT(IATOM,I)=LICT(II)                                             MK103510
  102 CONTINUE                                                          MK103520
  103 CONTINUE                                                          MK103530
C                                                                       MK103540
      DO 104 ISYM=1,NT                                                  MK103550
      WRITE(6,2) ISYM                                                   MK103560
      WRITE(6,3) (ICT(I,ISYM),I=1,NATOM)                                MK103570
      WRITE(6,4) ISYM                                                   MK103580
      CALL MATOUT(PTR(1,1,ISYM),3,3,3,3,6)                              MK103590
  104 CONTINUE                                                          MK103600
C                                                                       MK103610
      CALL RCLOSE(ITAP30,3)                                             MK103620
      RETURN                                                            MK103630
      END                                                               MK103640
      SUBROUTINE SYMDP(DP,DM,NDIM)                                      MK103650
      IMPLICIT REAL*8 (A-H,O-Z)                                         MK103660
      DIMENSION DP(3,NDIM),DM(3,NDIM)                                   MK103670
      DIMENSION D1T(3),DNU(3)                                           MK103680
      COMMON/COM101/NATOM,N3N                                           MK103690
      COMMON/COM104/NT,NCASE,ISPOT(150)                                 MK103700
      COMMON/COM105/IOFF(256),IPRNT                                     MK103710
      COMMON/COM106/PTR(3,3,48)                                         MK103720
      COMMON/COM107/ICT(150,48)                                         MK103730
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                MK103740
      DATA XLIM / 1.0D-4 /                                              MK103750
    1 FORMAT(//,2X,' DP MATRIX BEFORE SYMMETRIZATION'/)                 MK103760
    2 FORMAT(//,2X,' DM MATRIX BEFORE SYMMETRIZATION'/)                 MK103770
    3 FORMAT(//,2X,' DP BEFORE SYMMETRIZATION, IABC = ',I5/)            MK103780
    4 FORMAT(//,2X,' DM BEFORE SYMMETRIZATION, IABC = ',I5/)            MK103790
    5 FORMAT(//,2X,' GRADIENT AFTER SYMMETRIZATION, IT = ',I5/)         MK103800
    6 FORMAT(//,2X,' DP MATRIX AFTER SYMMETRIZATION'/)                  MK103810
    7 FORMAT(//,2X,' DM MATRIX AFTER SYMMETRIZATION'/)                  MK103820
C                                                                       MK103830
C*    IF(IPRNT.LE.2) GO TO 201                                          MK103840
C*    WRITE(6,1)                                                        MK103850
C*    CALL MATOUT(DP,3,NDIM,3,N3N,6)                                    MK103860
C*    WRITE(6,2)                                                        MK103870
C*    CALL MATOUT(DM,3,NDIM,3,N3N,6)                                    MK103880
C                                                                       MK103890
  201 CONTINUE                                                          MK103900
      DO 120 IABC=1,NCASE                                               MK103910
      IPOS=ISPOT(IABC)                                                  MK103920
      IATOM=(IPOS+2)/3                                                  MK103930
      IXYZ=IPOS-(IATOM-1)*3                                             MK103940
      IF(IPOS.GT.N3N) THEN                                              MK103950
        IATOM=(IPOS-N3N+2)/3                                            MK103960
        IXYZ=IPOS-N3N-(IATOM-1)*3                                       MK103970
      END IF                                                            MK103980
C                                                                       MK103990
C   POSITIVE PERTURBATIONS                                              MK104000
      IF(IPOS.GT.N3N) GO TO 301                                         MK104010
      D1T(1)=DP(1,IPOS)                                                 MK104020
      D1T(2)=DP(2,IPOS)                                                 MK104030
      D1T(3)=DP(3,IPOS)                                                 MK104040
      IF(IPRNT.LE.2) GO TO 202                                          MK104050
      WRITE(6,3) IABC                                                   MK104060
      CALL MATOUT(D1T,3,1,3,1,6)                                        MK104070
  202 CONTINUE                                                          MK104080
      GO TO 302                                                         MK104090
C                                                                       MK104100
C   NEGATIVE PERTURBATIONS                                              MK104110
  301 CONTINUE                                                          MK104120
      D1T(1)=DM(1,IPOS-N3N)                                             MK104130
      D1T(2)=DM(2,IPOS-N3N)                                             MK104140
      D1T(3)=DM(3,IPOS-N3N)                                             MK104150
      IF(IPRNT.LE.2) GO TO 302                                          MK104160
      WRITE(6,4) IABC                                                   MK104170
      CALL MATOUT(D1T,3,1,3,1,6)                                        MK104180
C                                                                       MK104190
  302 CONTINUE                                                          MK104200
      DO 110 IIT=1,NT                                                   MK104210
      JATOM=ICT(IATOM,IIT)                                              MK104220
      XNU=PTR(IXYZ,IXYZ,IIT)                                            MK104230
      IF(XNU.NE.ONE) GO TO 303                                          MK104240
      IF(JATOM.NE.IATOM) GO TO 303                                      MK104250
      GO TO 110                                                         MK104260
C                                                                       MK104270
C   SYMMETRIZE GRADIENT VECTOR                                          MK104280
  303 CONTINUE                                                          MK104290
      IT=IIT                                                            MK104300
      DEDX=A00                                                          MK104310
      DEDY=A00                                                          MK104320
      DEDZ=A00                                                          MK104330
      DEDXP=D1T(1)                                                      MK104340
      DEDYP=D1T(2)                                                      MK104350
      DEDZP=D1T(3)                                                      MK104360
      DEDX=DEDX+DEDXP*PTR(1,1,IT)+DEDYP*PTR(2,1,IT)+DEDZP*PTR(3,1,IT)   MK104370
      DEDY=DEDY+DEDXP*PTR(1,2,IT)+DEDYP*PTR(2,2,IT)+DEDZP*PTR(3,2,IT)   MK104380
      DEDZ=DEDZ+DEDXP*PTR(1,3,IT)+DEDYP*PTR(2,3,IT)+DEDZP*PTR(3,3,IT)   MK104390
      DNU(1)=DEDX                                                       MK104400
      DNU(2)=DEDY                                                       MK104410
      DNU(3)=DEDZ                                                       MK104420
      IF(IPRNT.LE.2) GO TO 203                                          MK104430
      WRITE(6,5) IT                                                     MK104440
      CALL MATOUT(DNU,3,1,3,1,6)                                        MK104450
C                                                                       MK104460
  203 CONTINUE                                                          MK104470
      DO 105 I=1,3                                                      MK104480
      XNU=PTR(IXYZ,I,IT)                                                MK104490
      IF(DABS(XNU).GT.XLIM) GO TO 304                                   MK104500
  105 CONTINUE                                                          MK104510
  304 CONTINUE                                                          MK104520
      JXYZ=I                                                            MK104530
      JPOS=(JATOM-1)*3+JXYZ                                             MK104540
      IF(IPOS.LE.N3N.AND.XNU.EQ.-ONE) JPOS=JPOS+N3N                     MK104550
      IF(IPOS.GT.N3N.AND.XNU.EQ.ONE) JPOS=JPOS+N3N                      MK104560
      IF(JPOS.GT.N3N) GO TO 305                                         MK104570
C                                                                       MK104580
C   POSITIVE PERTURBATION                                               MK104590
      DP(1,JPOS)=DNU(1)                                                 MK104600
      DP(2,JPOS)=DNU(2)                                                 MK104610
      DP(3,JPOS)=DNU(3)                                                 MK104620
  107 CONTINUE                                                          MK104630
      GO TO 110                                                         MK104640
C                                                                       MK104650
C   NEGATIVE PERTURBATION                                               MK104660
  305 CONTINUE                                                          MK104670
      DM(1,JPOS-N3N)=DNU(1)                                             MK104680
      DM(2,JPOS-N3N)=DNU(2)                                             MK104690
      DM(3,JPOS-N3N)=DNU(3)                                             MK104700
C                                                                       MK104710
  110 CONTINUE                                                          MK104720
  120 CONTINUE                                                          MK104730
C                                                                       MK104740
      WRITE(6,6)                                                        MK104750
      CALL MATOUT(DP,3,NDIM,3,N3N,6)                                    MK104760
      WRITE(6,7)                                                        MK104770
      CALL MATOUT(DM,3,NDIM,3,N3N,6)                                    MK104780
C                                                                       MK104790
      RETURN                                                            MK104800
      END                                                               MK104810
      SUBROUTINE SYMFX(FP,FM,NDIM)                                      MK104820
      IMPLICIT REAL*8 (A-H,O-Z)                                         MK104830
      DIMENSION FP(NDIM,NDIM),FM(NDIM,NDIM)                             MK104840
      DIMENSION D1T(3,150),DNU(3,150)                                   MK104850
      COMMON/COM101/NATOM,N3N                                           MK104860
      COMMON/COM104/NT,NCASE,ISPOT(150)                                 MK104870
      COMMON/COM105/IOFF(256),IPRNT                                     MK104880
      COMMON/COM106/PTR(3,3,48)                                         MK104890
      COMMON/COM107/ICT(150,48)                                         MK104900
      DATA A00,ONE / 0.0D+00 , 1.0D+00 /                                MK104910
      DATA XLIM / 1.0D-4 /                                              MK104920
    1 FORMAT(//,2X,' FP MATRIX BEFORE SYMMETRIZATION'/)                 MK104930
    2 FORMAT(//,2X,' FM MATRIX BEFORE SYMMETRIZATION'/)                 MK104940
    3 FORMAT(//,2X,' FP BEFORE SYMMETRIZATION, IABC = ',I5/)            MK104950
    4 FORMAT(//,2X,' FM BEFORE SYMMETRIZATION, IABC = ',I5/)            MK104960
    5 FORMAT(//,2X,' GRADIENT AFTER SYMMETRIZATION, IT = ',I5/)         MK104970
    6 FORMAT(//,2X,' FP MATRIX AFTER SYMMETRIZATION'/)                  MK104980
    7 FORMAT(//,2X,' FM MATRIX AFTER SYMMETRIZATION'/)                  MK104990
C                                                                       MK105000
C*    IF(IPRNT.LE.2) GO TO 201                                          MK105010
C*    WRITE(6,1)                                                        MK105020
C*    CALL MATOUT(FP,NDIM,NDIM,N3N,N3N,6)                               MK105030
C*    WRITE(6,2)                                                        MK105040
C*    CALL MATOUT(FM,NDIM,NDIM,N3N,N3N,6)                               MK105050
C                                                                       MK105060
  201 CONTINUE                                                          MK105070
      DO 120 IABC=1,NCASE                                               MK105080
      IPOS=ISPOT(IABC)                                                  MK105090
      IATOM=(IPOS+2)/3                                                  MK105100
      IXYZ=IPOS-(IATOM-1)*3                                             MK105110
      IF(IPOS.GT.N3N) THEN                                              MK105120
        IATOM=(IPOS-N3N+2)/3                                            MK105130
        IXYZ=IPOS-N3N-(IATOM-1)*3                                       MK105140
      END IF                                                            MK105150
C                                                                       MK105160
C   POSITIVE PERTURBATIONS                                              MK105170
      IF(IPOS.GT.N3N) GO TO 301                                         MK105180
      DO 102 I=1,NATOM                                                  MK105190
      IADD=(I-1)*3                                                      MK105200
      D1T(1,I)=FP(IPOS,IADD+1)                                          MK105210
      D1T(2,I)=FP(IPOS,IADD+2)                                          MK105220
      D1T(3,I)=FP(IPOS,IADD+3)                                          MK105230
  102 CONTINUE                                                          MK105240
      IF(IPRNT.LE.2) GO TO 202                                          MK105250
      WRITE(6,3) IABC                                                   MK105260
      CALL MATOUT(D1T,3,NDIM,3,NATOM,6)                                 MK105270
  202 CONTINUE                                                          MK105280
      GO TO 302                                                         MK105290
C                                                                       MK105300
C   NEGATIVE PERTURBATIONS                                              MK105310
  301 CONTINUE                                                          MK105320
      DO 103 I=1,NATOM                                                  MK105330
      IADD=(I-1)*3                                                      MK105340
      D1T(1,I)=FM(IPOS-N3N,IADD+1)                                      MK105350
      D1T(2,I)=FM(IPOS-N3N,IADD+2)                                      MK105360
      D1T(3,I)=FM(IPOS-N3N,IADD+3)                                      MK105370
  103 CONTINUE                                                          MK105380
      IF(IPRNT.LE.2) GO TO 302                                          MK105390
      WRITE(6,4) IABC                                                   MK105400
      CALL MATOUT(D1T,3,NDIM,3,NATOM,6)                                 MK105410
C                                                                       MK105420
  302 CONTINUE                                                          MK105430
      DO 110 IIT=1,NT                                                   MK105440
      JATOM=ICT(IATOM,IIT)                                              MK105450
      XNU=PTR(IXYZ,IXYZ,IIT)                                            MK105460
      IF(XNU.NE.ONE) GO TO 303                                          MK105470
      IF(JATOM.NE.IATOM) GO TO 303                                      MK105480
      GO TO 110                                                         MK105490
C                                                                       MK105500
C   SYMMETRIZE GRADIENT VECTOR                                          MK105510
  303 CONTINUE                                                          MK105520
      IT=IIT                                                            MK105530
      DO 104 IC=1,NATOM                                                 MK105540
      DEDX=A00                                                          MK105550
      DEDY=A00                                                          MK105560
      DEDZ=A00                                                          MK105570
      ICNU=ICT(IC,IT)                                                   MK105580
      DEDXP=D1T(1,IC)                                                   MK105590
      DEDYP=D1T(2,IC)                                                   MK105600
      DEDZP=D1T(3,IC)                                                   MK105610
      DEDX=DEDX+DEDXP*PTR(1,1,IT)+DEDYP*PTR(2,1,IT)+DEDZP*PTR(3,1,IT)   MK105620
      DEDY=DEDY+DEDXP*PTR(1,2,IT)+DEDYP*PTR(2,2,IT)+DEDZP*PTR(3,2,IT)   MK105630
      DEDZ=DEDZ+DEDXP*PTR(1,3,IT)+DEDYP*PTR(2,3,IT)+DEDZP*PTR(3,3,IT)   MK105640
      DNU(1,ICNU)=DEDX                                                  MK105650
      DNU(2,ICNU)=DEDY                                                  MK105660
      DNU(3,ICNU)=DEDZ                                                  MK105670
  104 CONTINUE                                                          MK105680
      IF(IPRNT.LE.2) GO TO 203                                          MK105690
      WRITE(6,5) IT                                                     MK105700
      CALL MATOUT(DNU,3,NDIM,3,NATOM,6)                                 MK105710
C                                                                       MK105720
  203 CONTINUE                                                          MK105730
      DO 105 I=1,3                                                      MK105740
      XNU=PTR(IXYZ,I,IT)                                                MK105750
      IF(DABS(XNU).GT.XLIM) GO TO 304                                   MK105760
  105 CONTINUE                                                          MK105770
  304 CONTINUE                                                          MK105780
      JXYZ=I                                                            MK105790
      JPOS=(JATOM-1)*3+JXYZ                                             MK105800
      IF(IPOS.LE.N3N.AND.XNU.EQ.-ONE) JPOS=JPOS+N3N                     MK105810
      IF(IPOS.GT.N3N.AND.XNU.EQ.ONE) JPOS=JPOS+N3N                      MK105820
      IF(JPOS.GT.N3N) GO TO 305                                         MK105830
C                                                                       MK105840
C   POSITIVE PERTURBATION                                               MK105850
      DO 107 I=1,NATOM                                                  MK105860
      IADD=(I-1)*3                                                      MK105870
      FP(JPOS,IADD+1)=DNU(1,I)                                          MK105880
      FP(JPOS,IADD+2)=DNU(2,I)                                          MK105890
      FP(JPOS,IADD+3)=DNU(3,I)                                          MK105900
  107 CONTINUE                                                          MK105910
      GO TO 110                                                         MK105920
C                                                                       MK105930
C   NEGATIVE PERTURBATION                                               MK105940
  305 CONTINUE                                                          MK105950
      DO 108 I=1,NATOM                                                  MK105960
      IADD=(I-1)*3                                                      MK105970
      FM(JPOS-N3N,IADD+1)=DNU(1,I)                                      MK105980
      FM(JPOS-N3N,IADD+2)=DNU(2,I)                                      MK105990
      FM(JPOS-N3N,IADD+3)=DNU(3,I)                                      MK106000
  108 CONTINUE                                                          MK106010
C                                                                       MK106020
  110 CONTINUE                                                          MK106030
  120 CONTINUE                                                          MK106040
C                                                                       MK106050
      WRITE(6,6)                                                        MK106060
      CALL MATOUT(FP,NDIM,NDIM,N3N,N3N,6)                               MK106070
      WRITE(6,7)                                                        MK106080
      CALL MATOUT(FM,NDIM,NDIM,N3N,N3N,6)                               MK106090
C                                                                       MK106100
      RETURN                                                            MK106110
      END                                                               MK106120
