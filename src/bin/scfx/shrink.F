C     ////////////////////////////////////
      SUBROUTINE SHRINK(BIG,SMAL,NB,NS)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DOUBLE PRECISION BIG(NB,NB),SMAL(NS),VAL,TOL,CUT,A0
      COMMON /TAPES/NIN,IR,IOUT,ICHECK,IFLAG3,ITAPE,ITERML
CWA   DATA CUT /1.0D+09/
      DATA CUT /1.0D-10/
      DATA A0 /0.0D+00/
      DATA TOL /1.0D-07/
C
      IJ=0
      DO 30 I=1,NB
        DO 40 J=1,I
          IJ=IJ+1
          VAL=BIG(I,J)
          IF (ABS(VAL) .LT. CUT) VAL=A0
          SMAL(IJ)=VAL
          IF(ABS(BIG(J,I)) .LT. CUT) GOTO 50
          IF (VAL .EQ. 0.0) THEN
            IF (ABS(VAL-BIG(J,I)) .GT. TOL) WRITE(ICHECK,101) I,J,
     *      VAL,BIG(J,I)
          ELSE
            IF (ABS((VAL-BIG(J,I))/VAL) .GT. TOL)WRITE(ICHECK,101) I
     *      ,J,VAL,BIG(J,I)
          END IF
 50       CONTINUE
 40     CONTINUE
 30   CONTINUE
 101  FORMAT('  ASYMMETRY IN SHRINK: I,J,VAL1,VAL2=',2I5,2G18.8)
      RETURN
      END
