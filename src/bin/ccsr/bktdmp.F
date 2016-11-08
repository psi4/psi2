C
C----------------------------------------------------------------------
C
      SUBROUTINE BKTDMP(BKT,LENGTH,IDISK)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BKT(LENGTH)
      CALL SWRIT(IDISK,BKT,INTOWP(LENGTH))
      RETURN
      END
