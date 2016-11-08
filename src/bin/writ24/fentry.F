      SUBROUTINE FENTRY(RARRAY,IARRAY,MAXCOR)
C***********************************************************************
C*  THIS PROGRAM TAKES AN INPUT FILE OF CARTESIAN THIRD DERIVATIVES    *
C* (FILETYPE = "TOTAL20" AND WRITES A "FILE24" OF FOURTH DERIVATIVES.  *
C   PROGRAM IS SETUP TO HANDLE 20 ATOMS.                               *
C***********************************************************************
c Moved to PSI distribution disk on 020389 - clj.
C***********************************************************************
C*   MODIFICATION FOR IMS VERSIO                                       *
C*   BY: YUKIO YAMAGUCHI                                               *
C*   DATE: FEBRUARY 21, 1989                                           *
C***********************************************************************
C*  PROGRAM LAST MODIFIED ON: OCT.  01, 1985 BY:  RBR                  *
C*      REASON:  FIRST IMPLEMENTATION                                  *
C***********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RARRAY(MAXCOR),IARRAY(1)
C
      DIMENSION X3RDP(37820),X3RDM(37820),TMP4TH(37820),FD4TH(595665)
      DIMENSION I3LEN(20*3),I4LEN(20*3)
C
      DATA ANGBHR,ZERO / 0.52917706D+00, 0.0D+00 /
C
    1 FORMAT(10I5)
    2 FORMAT(//,2X,' THE CARTESIAN FINITE DIFFERENCE PROGRAM FOR '
     1   ,'FOURTH DERIVATIVES: WRIT24 ',//,
     2 2X,' ISYM30 = ',I5/
     3 2X,' I2PNTD = ',I5/
     4 2X,' IFLAG3 = ',I5/
     5 2X,' IPRINT = ',I5/
     6 2X,' FINDIF = ',F10.6/
     7 2X,' NATOMS = ',I5/
     8 2X,' NN3N   = ',I5/
     9 2X,' N3TOT  = ',I5/
     A 2X,' N4TOT  = ',I5/)
    3 FORMAT(/,3I5)
    4 FORMAT(3F20.10)
    5 FORMAT(/,' NATOMS,ATOM,AXIS ',3I5)
    6 FORMAT(/,15X,' >>>>> PROGRAM WRIT24 COMPLETED <<<<< ',/)
C
      call drum
cets  CALL NOUNFL
      CALL TSTART(6)
C
      ICHEK  =  3
      INPUT  =  5
      IOUT   =  6
      ITAP20 = 20
      ITAP24 = 24
      call ffile(itap20,'total20',0)
      call ffile(itap24,' ',0)
C
C     TAPE 20 IS "TOTAL20" THE CUMULATIVE TOTAL OF FILE20'S BY FIN-DIFF
C
      ISYM30  =  0
      I2PNTD  =  2
      IFLAG3  =  0
      IPRINT  =  0
      FINDIF  =  0.01D+00
C
      CALL LOCATE(INPUT,'# FDGEOM #',IERR)
      IF(IERR .NE. 0) GOTO 100
          READ(INPUT,1) ISYM30,I2PNTD,IFLAG3,IPRINT
          READ(INPUT,*) FINDIF
              IF(ISYM30 .EQ. 0) ISYM30 = 0
              IF(I2PNTD .EQ. 0) I2PNTD = 2
              IF(IFLAG3 .EQ. 0) IFLAG3 = 0
              IF(IPRINT .EQ. 0) IPRINT = 0
              IF(FINDIF .EQ. 0.0) FINDIF = 0.01D+00
  100 CONTINUE
C
      REWIND(ITAP20)
C
      READ(ITAP20,3) NATOM,NAXIS,NDISPS
      READ(ITAP20,1) NATOMS,N3CHK
C
      N3N   = NATOMS * 3
      NN3N  = NATOMS * 3 * NATOMS
      N3TOT = N3N*(N3N+1)*(N3N+2)/6
      N4TOT = N3N*(N3N+1)*(N3N+2)*(N3N+3)/(2*3*4)
      DIS2  = 2 * FINDIF
C
      DO 110  II=1,3*NATOMS
          I3LEN(II) = II*(II+1)*(II+2)/6
          I4LEN(II) = I3LEN(II)*(II+3)/4
  110 CONTINUE
C
      REWIND(ITAP20)
      WRITE(IOUT,2) ISYM30,I2PNTD,IFLAG3,IPRINT,FINDIF,NATOMS,
     1              NN3N,N3TOT,N4TOT
      WRITE(ITAP24,1) NATOMS,N4TOT
C
      IF(N3CHK.NE.N3TOT) STOP 99999
C
      DO 120  II=1,N4TOT
          FD4TH(II) = ZERO
  120 CONTINUE
C
C     READ IN THE FILE 20'S ---- BOTH PLUS AND MINUS DISPLACEMENTS
C                                FOR A SINGLE ATOM AND AXIS
C     I2PNTD = 1 FOR SINGLE POINT, AND =2 FOR DOUBLE POINT DISPLACEMENTS
C
      DO 700 NFD=1,NDISPS/I2PNTD
C
C
      DO 210 I=1,N3TOT
          X3RDP(I)  = ZERO
          X3RDM(I)  = ZERO
          TMP4TH(I) = ZERO
  210 CONTINUE
C
      READ(ITAP20,3) IATOM,IAXIS
      READ(ITAP20,1) IATOMS,I3CHK
      READ(ITAP20,4) (X3RDP(JJ),JJ=1,N3TOT)
C
C
      READ(ITAP20,3) JATOM,JAXIS
      READ(ITAP20,1) JATOMS,J6
C
CCC   IF(I2PTND.EQ.1) DO SOMETHING (BACKUP 20 AND GET DER4 FROM REFF)
      IF(IATOM.NE.JATOM .OR. IAXIS.NE.-JAXIS) STOP 99999
C
      READ(ITAP20,4) (X3RDM(JJ),JJ=1,N3TOT)
C
      DO 240  I=1,N3TOT
          TMP4TH(I) = (X3RDP(I) - X3RDM(I))/DIS2
  240 CONTINUE
C
C
      INDX   = (IATOM-1)*3 + IAXIS
      IOFF4  = I4LEN(INDX)-I3LEN(INDX)
C
      DO 300  II=1,I3LEN(INDX)
              FD4TH(IOFF4+II) = TMP4TH(II)
  300 CONTINUE
C
  700 CONTINUE
C
C
      WRITE(ITAP24,4) (FD4TH(II),II=1,N4TOT)
C
      WRITE(IOUT,6)
      CALL TSTOP(6)
C
C2-21-89 STOP
      RETURN
      END
