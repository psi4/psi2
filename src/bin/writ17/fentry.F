      SUBROUTINE FENTRY(RARRAY,IARRAY,MAXCOR)
C***********************************************************************
C* THIS PROGRAM TAKES AN INPUT FILE OF CARTESIAN DIPOLE MOMENT VECTORS *
C* AND WRITES A "FILE17" WHICH IS SUITABLE FOR USE WITH NORMCO$        *
C* It is just a finite difference and scaling calculation.   RBR       *
C***********************************************************************
C*   MODIFICATION FOR IMS VERSION                                      *
C*   BY: YUKIO YAMAGUCHI                                               *
C*   DATE: FEBRUARY 21, 1989                                           *
C***********************************************************************
c Moved to PSI distribution disk on 020389 - clj.
C***********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RARRAY(MAXCOR),IARRAY(1)
C
      DIMENSION X1(100),X2(100),Y1(100),Y2(100),Z1(100),Z2(100)
      DIMENSION XX(100),YY(100),ZZ(100)
C
      DATA ANGBHR,ZERO / 0.52917706D+00, 0.0D+00 /
C
    1 FORMAT(10I5)
    2 FORMAT(3F20.10)
    3 FORMAT(//,3F20.10)
C
      call drum
      ICHEK  =  3
      INPUT  =  5
      IOUT   =  6
      ITAP17 = 17
      ITAP18 = 18
      call ffile(itap17,' ',0)
      call ffile(itap18,' ',0)
C
      REWIND(ITAP18)
C
      DO 100 I=1,100
          X1(I) = ZERO
          Y1(I) = ZERO
          Z1(I) = ZERO
          X2(I) = ZERO
          Y2(I) = ZERO
          Z2(I) = ZERO
          XX(I) = ZERO
          YY(I) = ZERO
          ZZ(I) = ZERO
  100 CONTINUE
C
      READ(ITAP18,1) NATOM,IFLAG
      N3N  = 3 * NATOM
      READ(ITAP18,2) DIS
      DIS2 = 2.0D+00 * DIS * ANGBHR
      WRITE(ITAP17,1) NATOM,N3N
C
      DO 110  I=1,N3N
          READ(ITAP18,3) X1(I),Y1(I),Z1(I)
          READ(ITAP18,3) X2(I),Y2(I),Z2(I)
          XX(I) = (X1(I) - X2(I))/DIS2
          YY(I) = (Y1(I) - Y2(I))/DIS2
          ZZ(I) = (Z1(I) - Z2(I))/DIS2
  110 CONTINUE
C
      DO 120  I=1,NATOM
          II  = I + NATOM
          III = I + NATOM*2
          WRITE(ITAP17,2) XX(I),XX(II),XX(III)
  120 CONTINUE
C
      DO 130  I=1,NATOM
          II  = I + NATOM
          III = I + NATOM*2
          WRITE(ITAP17,2) YY(I),YY(II),YY(III)
  130 CONTINUE
C
      DO 140  I=1,NATOM
          II  = I + NATOM
          III = I + NATOM*2
          WRITE(ITAP17,2) ZZ(I),ZZ(II),ZZ(III)
  140 CONTINUE
C
C2-21-89 STOP
      RETURN
      END
