C
C                            ********************
C----------------------------****** REINTS ******--------------------
C                            ********************
C
      SUBROUTINE REINTS (IDISK,MOINTS,LENINT,IDIM,JOUT,IAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MOINTS(LENINT)
      CALL RSETSA(IDISK,IAD)
      CALL SREAD (IDISK,MOINTS,INTOWP(LENINT))
      CALL MATOUT(MOINTS,IDIM,IDIM,IDIM,IDIM,JOUT)
      RETURN
      END
