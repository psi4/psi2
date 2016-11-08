      SUBROUTINE FENTRY(CC,IA,MAXCOR)
C*********************************************************************
C  PROGRAM HVIB FOR 20 OR LESS ATOMS---GIVES FREQUENCIES FROM        *
C  GRADIENTS OF CARTESIAN DISPLACED SYMMETRY DISTINCT ATOMS          *
C*********************************************************************
c Moved to PSI distribution disk on 020389 - clj.
C*********************************************************************
C  NOTICE OF PROGRAM MODIFICATION                                    *
C        BY:  YUKIO YAMAGUCHI                                        *
C      DATE:  FEBRUARY 23, 1989                                      *
C    REASON:  MODIFICATION FOR IMS VERSION                           *
C*********************************************************************
C  NOTICE OF PROGRAM MODIFICATION              IBM-4381 VERSION      *
C        BY:  RICHARD REMINGTON                                      *
C      DATE:  NOVEMBER 18, 1984                                      *
C    REASON:  MODIFY  TO HANDLE MORE ATOMS                           *
C*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(MAXCOR),IA(1)
      COMMON A(30000)
      REAL*8 MDYNE
      INTEGER SYM,FLAG,GEOM,YWGHT,YFLAG
      COMMON /PARAM/VECT(60),ANGS,FRCS,MDYNE
     1,SYM(60),ICEN,N,ISYM,ITYPE(20),IATNO(20),ITAPE5
      COMMON /MASSES/ AMASS(30)
C
      call drum
      CALL TSTART(6)
cets  CALL NOUNFL
C
      ITAPE5=15
      call ffile(itape5,' ',0)
      READ(ITAPE5,9011) ICEN
      DO 20 I=1,30000
      A(I)=0.0D+00
  20  CONTINUE
 1    FORMAT(4X,I1)
 9011 FORMAT(4X,I2)
      N = 3 * ICEN
      IH1 = 1
      IH2 = IH1 + N * N * 2
      IH3 = IH2 + N * N * 2
      IH4 = IH3 + N * N * 2
      IH5 = IH4 + N * N * 2
      IH6 = IH5 + N * N * 2
      YWGHT = IH6 + ICEN * 2
      GEOM = YWGHT + N * 2
CRBR  PRINT  *, 'ICEN BEFORE CALL INPUT ',ICEN
      CALL INPUT(A(IH1),A(IH2),A(IH4),A(IH6),A(GEOM)
     *,          FLAG,YFLAG)
      CALL SYMM(A(IH1),A(IH2),A(IH6),FLAG)
      IF(FLAG.GT.0) GO TO 10
       CALL CRUNCH(A(IH1),A(IH2),A(IH3),A(IH4),A(IH5),A(IH6)
     *,            A(YWGHT),A(GEOM),YFLAG)
   10  CONTINUE
C
      CALL TSTOP(6)
C
C2-23-89 STOP
       RETURN
       END
