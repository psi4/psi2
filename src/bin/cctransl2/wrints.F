C
C                            ********************
C----------------------------****** WRINTS ******--------------------
C                            ********************
C
      SUBROUTINE WRINTS (IDISK,MOINTS,LENINT,IAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MOINTS(LENINT)
      CALL RGETSA(IDISK,IAD)
      CALL SWRIT (IDISK,MOINTS,INTOWP(LENINT))
      RETURN
      END
