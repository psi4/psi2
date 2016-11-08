      SUBROUTINE SETNOS
      IMPLICIT REAL*8 (A-H,O-Z)
cets030691
#include <error.h>
      integer ip, prcntr
CTPH  COMMON/GAUSC/T(36),PTR(3,3,8),DTR(6,6,8)
      COMMON/GAUSC/T(225),PTR(3,3,8),DTR(6,6,8),FTR(10,10,8),
     1             GTR(15,15,8)
      COMMON/GAUSD/INVT(8),ISO(200,8),LBLAT(2,50),ICT(50,8)
      COMMON/ISPAC/INDIN(120),INDOUT(30)
      COMMON/PARA1/NSHELL,NT,MXT,NISO
      COMMON/SIGNS/IOFF(256),IPRNT
      COMMON/TOLER/DLIM
      DATA ZERO,ONE / 0.0D+00 , 1.0D+00 /
C
C   FORCE THE NO-SYMMETRY CALCULATION FOR CI SECOND DERIVATIVE
    1 FORMAT(/'   $$$$$ SET NO-SYMMETRY FLAG $$$$$'/
     1,       '   VALUES BEFORE SETNOS    NTRANS=',I3)
    2 FORMAT( '  TRANS=',9F10.6)
    3 FORMAT( '   ISO =',10I10)
    4 FORMAT(/'   %%%%%%%%%% AFTER SETNO  NTRANS=',I3)
C
cets030691      IF(IPRNT.LE.3) GO TO 201
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        WRITE(6,1) NT
        WRITE(6,2) (T(I),I=1,MXT)
        WRITE(6,3) ((ISO(I,J),I=1,NSHELL),J=1,NISO)
      endif
C
  201 CONTINUE
      NT=1
      MXT=9
      DO 101 I=1,9
  101 T(I)=ZERO
      T(1)=ONE
      T(5)=ONE
      T(9)=ONE
      DO 102 I=1,NSHELL
      INDIN(1)=I
cges  CALL ISOIN
      call ysoin
      ISO(I,1)=INDOUT(1)
  102 CONTINUE
cets030691      IF(IPRNT.LE.3) GO TO 202
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        WRITE(6,4) NT
        WRITE(6,2) (T(I),I=1,MXT)
        WRITE(6,3) (ISO(I,1),I=1,NSHELL)
      endif
C
  202 CONTINUE
      RETURN
      END
