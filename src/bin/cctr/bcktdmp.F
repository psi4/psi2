C
C                           ******************
C---------------------------***** BKTDMP *****-------------------------
C                           ******************
C
      SUBROUTINE BKTDMP(BKT,LENGTH,IDISK,ICHAN)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BKT(LENGTH)
      CALL RGETSA(IDISK,ICHAN)
      CALL SWRIT(IDISK,BKT,INTOWP(LENGTH))
      RETURN
      END
