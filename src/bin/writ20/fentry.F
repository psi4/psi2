      SUBROUTINE FENTRY(RARRAY,IARRAY,MAXCOR)
C**********************************************************************
C* THIS PROGRAM TAKES AN INPUT FILE OF CARTESIAN SECOND DERIVATIVES   *
C* (FILETYPE = "TOTAL15" AND WRITES A "FILE20" OF THIRD DERIVATIVES   *
C**********************************************************************
C*  MODIFICATION FOR IMS VERSION                                      *
C*  BY: YUKIO YAMAGUCHI                                               *
C*  DATE: FEBRUARY 21, 1989                                           *
C**********************************************************************
c Moved to PSI distribution disk on 020389 - clj.                     *
C**********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RARRAY(MAXCOR),IARRAY(1)
C
      DIMENSION X2NDP(22500),X2NDM(22500),TMP3RD(22500),FD3RD(573800)
C
      DATA ANGBHR,ZERO / 0.52917706D+00, 0.0D+00 /
C
    1 FORMAT(10I5)
    2 FORMAT(//,2X,' THE CARTESIAN FINITE DIFFERENCE PROGRAM FOR '
     1   ,'THIRD DERIVATIVES: WRIT20 ',//,
     2 2X,' ISYM30 = ',I5/
     3 2X,' I2PNTD = ',I5/
     4 2X,' IFLAG3 = ',I5/
     5 2X,' IPRINT = ',I5/
     6 2X,' FINDIF = ',F10.6/
     7 2X,' NATOMS = ',I5/
     8 2X,' NN3N   = ',I5/
     9 2X,' N6     = ',I5/)
    3 FORMAT(/,3I5)
    4 FORMAT(3F20.10)
    5 FORMAT(/,' NATOMS,ATOM,AXIS ',3I5)
    6 FORMAT(/,10X,' >>>>> END <<<<< ',/)
C
cets  CALL NOUNFL
      call drum
      CALL TSTART(6)
C
      ICHEK  =  3
      INPUT  =  5
      IOUT   =  6
      ITAP15 = 15
      ITAP20 = 20
      ITAP21 = 21
      call ffile(itap15,'total15',0)
      call ffile(itap20,' ',0)
      call ffile(itap21,' ',0)
C
      ISYM30  =  0
      I2PNTD  =  2
      IFLAG3  =  0
      IPRINT  =  0
      FINDIF  =  0.001D+00
C
      CALL LOCATE(INPUT,'# FDGEOM #',IERR)
      IF(IERR .NE. 0) GOTO 100
          READ(INPUT,1) ISYM30,I2PNTD,IFLAG3,IPRINT
          READ(INPUT,*) FINDIF
              IF(ISYM30 .EQ. 0) ISYM30 = 0
              IF(I2PNTD .EQ. 0) I2PNTD = 2
              IF(IFLAG3 .EQ. 0) IFLAG3 = 0
              IF(IPRINT .EQ. 0) IPRINT = 0
              IF(FINDIF .EQ. 0.0) FINDIF = 0.001D+00
  100 CONTINUE
C
      REWIND(ITAP15)
C
      READ(ITAP15,3) NATOM,NAXIS,NDISPS
      READ(ITAP15,1) NATOMS,N6
C
      N3N   = NATOMS * 3
      NN3N  = NATOMS * 3 * NATOMS
      N3TOT = N3N*(N3N+1)*(N3N+2)/6
      DIS2  = 2 * FINDIF
C
      REWIND(ITAP15)
      WRITE(IOUT,2) ISYM30,I2PNTD,IFLAG3,IPRINT,FINDIF,NATOMS,NN3N,N6
      WRITE(ITAP20,1) NATOMS,N3TOT
C
C
      DO 120  II=1,N3TOT
          FD3RD(II) = ZERO
  120 CONTINUE
C
C
      DO 700 NFD=1,NDISPS/I2PNTD
C
C
      DO 210 I=1,NN3N*3
          X2NDP(I)  = ZERO
          X2NDM(I)  = ZERO
          TMP3RD(I) = ZERO
  210 CONTINUE
C
      READ(ITAP15,3) IATOM,IAXIS
      READ(ITAP15,1) IATOMS,I6
C
      DO 220  I=1,NN3N
          II = (I-1)*3
          READ(ITAP15,4) X2NDP(II+1),X2NDP(II+2),X2NDP(II+3)
  220 CONTINUE
C
      READ(ITAP15,3) JATOM,JAXIS
      READ(ITAP15,1) JATOMS,J6
C
      IF(IATOM.NE.JATOM .OR. IAXIS.NE.-JAXIS) STOP 99999
C
      DO 230  I=1,NN3N
          II = (I-1)*3
          READ(ITAP15,4) X2NDM(II+1),X2NDM(II+2),X2NDM(II+3)
  230 CONTINUE
C
      DO 240  I=1,NN3N*3
          TMP3RD(I) = (X2NDP(I) - X2NDM(I))/DIS2
  240 CONTINUE
C
CCC   WRITE(ITAP21,5) IATOMS,IATOM,IAXIS
CCC
CCC   DO 250  I=1,NN3N
CCC       II = (I-1)*3
CCC       WRITE(ITAP21,4) TMP3RD(II+1),TMP3RD(II+2),TMP3RD(II+3)
CC250 CONTINUE
C
C
      INDX   = (IATOM-1)*3 + IAXIS
      NSTART = (INDX-1)*(INDX)*(INDX+1)/6
      NSTOP  = (INDX)*(INDX+1)*(INDX+2)/6
      NLEN   = NSTOP - NSTART
      JJ     = NSTART
C
      DO 300  IATM=1,INDX
          IND  = (IATM-1)*N3N
          DO 310  II=1,IATM
              JJ  =  JJ + 1
              IND  = IND + 1
              FD3RD(JJ) = TMP3RD(IND)
  310     CONTINUE
  300 CONTINUE
C
  700 CONTINUE
C
C
      DO 600  I=1,(N3TOT+2)/3
          II = (I-1)*3
          WRITE(ITAP20,4) FD3RD(II+1),FD3RD(II+2),FD3RD(II+3)
  600 CONTINUE
C
      WRITE(IOUT,6)
      CALL TSTOP(6)
C
C2-21-89 STOP
      RETURN
      END
