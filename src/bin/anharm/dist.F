      SUBROUTINE DIST(IPRNT)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/VIB101/NATOM,N3N,NATRI,ILIN,NVIB
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50),W(50)
      COMMON/VIB202/R(1275)
      DATA A00 / 0.0D+00 /
    1 FORMAT(//,2X,' INTERATOMIC DISTANCE MATRIX'/)
C
      IJ=0
      DO 101 I=1,NATOM
      DO 101 J=1,I
      IJ=IJ+1
      R(IJ)=A00
      IF(I.EQ.J) GO TO 101
      XD=X(I)-X(J)
      YD=Y(I)-Y(J)
      ZD=Z(I)-Z(J)
      R(IJ)=DSQRT(XD*XD+YD*YD+ZD*ZD)
  101 CONTINUE
      IF(IPRNT.EQ.0) GO TO 201
      WRITE(6,1)
      CALL PRINT(R,1275,NATOM,6)
  201 RETURN
      END
