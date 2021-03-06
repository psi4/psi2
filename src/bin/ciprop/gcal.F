      SUBROUTINE GCAL(PCX,PCY,PCZ,VX,VY,VZ,A22,NI,NJ,MAXT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION GX(10),GY(10),GZ(10)
cets072891      COMMON/EXYZS/EX(7,7,13),EY(7,7,13),EZ(7,7,13)
      common/exyzs/ex(7,7,14),ey(7,7,14),ez(7,7,14)
      DATA ZERO,ONE / 0.0D+00 , 1.0D+00 /
C
      GX(1)=ONE
      GY(1)=ONE
      GZ(1)=ONE
      GX(2)=PCX*A22
      GY(2)=PCY*A22
      GZ(2)=PCZ*A22
C
      IF(MAXT.LE.2) GO TO 201
      DO 101 I=2,MAXT-1
      XI1=FLOAT(I-1)
      GX(I+1)=(PCX*GX(I)+XI1*GX(I-1))*A22
      GY(I+1)=(PCY*GY(I)+XI1*GY(I-1))*A22
      GZ(I+1)=(PCZ*GZ(I)+XI1*GZ(I-1))*A22
  101 CONTINUE
C
  201 CONTINUE
      VX=ZERO
      VY=ZERO
      VZ=ZERO
      DO 102 I=1,MAXT
      VX=VX+EX(NI,NJ,I)*GX(I)
      VY=VY+EY(NI,NJ,I)*GY(I)
      VZ=VZ+EZ(NI,NJ,I)*GZ(I)
  102 CONTINUE
C
      RETURN
      END
