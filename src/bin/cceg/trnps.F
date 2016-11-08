C
C ----------------------------------------------------------------------
C
      SUBROUTINE TRNPS (TT1,T1,NO,NV)
      IMPLICIT INTEGER (A-Z)
      REAL*8 TT1(NV,NO),T1(NO,NV)
      DO 102 I = 1, NO
          DO 202 J = 1, NV
             TT1(J,I) = T1(I,J)
202          CONTINUE
102       CONTINUE
      RETURN
      END
