C     //////////////////////////////////////////////////////
      SUBROUTINE EXPAND(BIG,SMAL,NB,NS)
 
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DOUBLE PRECISION BIG(NB,NB),SMAL(NS),VAL,TOL,CUT,A0
C
      DATA CUT /1.0D-09/
      DATA A0 /0.0D+00/
C
      IJ=0
      DO 10 I=1,NB
        DO 20 J=1,I
          IJ=IJ+1
          VAL=SMAL(IJ)
CWA       IF (ABS(VAL) .LT. CUT) VAL=A0
          BIG(I,J)=VAL
          BIG(J,I)=VAL
 20     CONTINUE
 10   CONTINUE
      RETURN
      END
