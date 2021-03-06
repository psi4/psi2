      SUBROUTINE FRMIJ
      IMPLICIT REAL*8 (A-H,O-Z)
CTPH  COMMON/ABDEN/ABDENS(1296)
c11-14-88  (For  7Mb)  COMMON/ABDEN/ABDENS(10000)
c11-14-88  (For 15Mb)  COMMON/ABDEN/ABDENS(50625)
      COMMON/ABDEN/ABDENS(10000)
CTPH  COMMON/DINDX/IJX(100),IJY(100),IJZ(100),KLX(100),KLY(100),KLZ(100)
      COMMON/DINDX/IJX(225),IJY(225),IJZ(225),KLX(225),KLY(225),KLZ(225)
      COMMON/D1INT/XIDER,YIDER,ZIDER,XJDER,YJDER,ZJDER,
     1             XKDER,YKDER,ZKDER
      COMMON/PARA4/MIJ,MKL,ISAME(100)
      COMMON/ROOTS/NONZER
      COMMON/TOLER/DLIM
CTPH  COMMON/TWOD0/XIN(2304),YIN(2304),ZIN(2304)
CTPH  COMMON/TWODI/DXI(2304),DYI(2304),DZI(2304)
CTPH  COMMON/TWODJ/DXJ(2304),DYJ(2304),DZJ(2304)
      COMMON/TWOD0/XIN(6250),YIN(6250),ZIN(6250)
      COMMON/TWODI/DXI(6250),DYI(6250),DZI(6250)
      COMMON/TWODJ/DXJ(6250),DYJ(6250),DZJ(6250)
C
      N=0
      GO TO (100,200,300,400,500,600,700,800,900),NONZER
C
  100 CONTINUE
      DO 101 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 101 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DENSTY=ABDENS(N)
      IF(DABS(DENSTY).LT.DLIM) GO TO 101
      XY=XIN(MX) * YIN(MY) * DENSTY
      XZ=XIN(MX) * ZIN(MZ) * DENSTY
      YZ=YIN(MY) * ZIN(MZ) * DENSTY
C
      XIDER=XIDER + DXI(MX) * YZ
      YIDER=YIDER + DYI(MY) * XZ
      ZIDER=ZIDER + DZI(MZ) * XY
C
      XJDER=XJDER + DXJ(MX) * YZ
      YJDER=YJDER + DYJ(MY) * XZ
      ZJDER=ZJDER + DZJ(MZ) * XY
C
  101 CONTINUE
      RETURN
C
  200 CONTINUE
      DO 102 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 102 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DENSTY=ABDENS(N)
      IF(DABS(DENSTY).LT.DLIM) GO TO 102
      XY1 = XIN(MX) * YIN(MY) * DENSTY
      XZ1 = XIN(MX) * ZIN(MZ) * DENSTY
      YZ1 = YIN(MY) * ZIN(MZ) * DENSTY
C
      XY2 = XIN(MX+625) * YIN(MY+625) * DENSTY
      XZ2 = XIN(MX+625) * ZIN(MZ+625) * DENSTY
      YZ2 = YIN(MY+625) * ZIN(MZ+625) * DENSTY
C
      XIDER = XIDER + DXI(MX) * YZ1 + DXI(MX+625) * YZ2
      YIDER = YIDER + DYI(MY) * XZ1 + DYI(MY+625) * XZ2
      ZIDER = ZIDER + DZI(MZ) * XY1 + DZI(MZ+625) * XY2
C
      XJDER = XJDER + DXJ(MX) * YZ1 + DXJ(MX+625) * YZ2
      YJDER = YJDER + DYJ(MY) * XZ1 + DYJ(MY+625) * XZ2
      ZJDER = ZJDER + DZJ(MZ) * XY1 + DZJ(MZ+625) * XY2
C
  102 CONTINUE
      RETURN
C
  300 CONTINUE
      DO 103 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 103 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DENSTY=ABDENS(N)
      IF(DABS(DENSTY).LT.DLIM) GO TO 103
      XY1 = XIN(MX     ) * YIN(MY     ) * DENSTY
      XZ1 = XIN(MX     ) * ZIN(MZ     ) * DENSTY
      YZ1 = YIN(MY     ) * ZIN(MZ     ) * DENSTY
C
      XY2 = XIN(MX+ 625) * YIN(MY+ 625) * DENSTY
      XZ2 = XIN(MX+ 625) * ZIN(MZ+ 625) * DENSTY
      YZ2 = YIN(MY+ 625) * ZIN(MZ+ 625) * DENSTY
C
      XY3 = XIN(MX+  1250) * YIN(MY+  1250) * DENSTY
      XZ3 = XIN(MX+  1250) * ZIN(MZ+  1250) * DENSTY
      YZ3 = YIN(MY+  1250) * ZIN(MZ+  1250) * DENSTY
C
      XIDER = XIDER + DXI(MX     ) * YZ1 + DXI(MX+ 625) * YZ2
     1              + DXI(MX+  1250) * YZ3
      YIDER = YIDER + DYI(MY     ) * XZ1 + DYI(MY+ 625) * XZ2
     1              + DYI(MY+  1250) * XZ3
      ZIDER = ZIDER + DZI(MZ     ) * XY1 + DZI(MZ+ 625) * XY2
     1              + DZI(MZ+  1250) * XY3
C
      XJDER = XJDER + DXJ(MX     ) * YZ1 + DXJ(MX+ 625) * YZ2
     1              + DXJ(MX+  1250) * YZ3
      YJDER = YJDER + DYJ(MY     ) * XZ1 + DYJ(MY+ 625) * XZ2
     1              + DYJ(MY+  1250) * XZ3
      ZJDER = ZJDER + DZJ(MZ     ) * XY1 + DZJ(MZ+ 625) * XY2
     1              + DZJ(MZ+  1250) * XY3
C
  103 CONTINUE
      RETURN
C
  400 CONTINUE
      DO 104 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 104 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DENSTY=ABDENS(N)
      IF(DABS(DENSTY).LT.DLIM) GO TO 104
      XY1 = XIN(MX     ) * YIN(MY     ) * DENSTY
      XZ1 = XIN(MX     ) * ZIN(MZ     ) * DENSTY
      YZ1 = YIN(MY     ) * ZIN(MZ     ) * DENSTY
C
      XY2 = XIN(MX+ 625) * YIN(MY+ 625) * DENSTY
      XZ2 = XIN(MX+ 625) * ZIN(MZ+ 625) * DENSTY
      YZ2 = YIN(MY+ 625) * ZIN(MZ+ 625) * DENSTY
C
      XY3 = XIN(MX+  1250) * YIN(MY+  1250) * DENSTY
      XZ3 = XIN(MX+  1250) * ZIN(MZ+  1250) * DENSTY
      YZ3 = YIN(MY+  1250) * ZIN(MZ+  1250) * DENSTY
C
      XY4 = XIN(MX+  1875) * YIN(MY+  1875) * DENSTY
      XZ4 = XIN(MX+  1875) * ZIN(MZ+  1875) * DENSTY
      YZ4 = YIN(MY+  1875) * ZIN(MZ+  1875) * DENSTY
C
      XIDER = XIDER + DXI(MX     ) * YZ1 + DXI(MX+ 625) * YZ2
     1              + DXI(MX+  1250) * YZ3 + DXI(MX+  1875) * YZ4
      YIDER = YIDER + DYI(MY     ) * XZ1 + DYI(MY+ 625) * XZ2
     1              + DYI(MY+  1250) * XZ3 + DYI(MY+  1875) * XZ4
      ZIDER = ZIDER + DZI(MZ     ) * XY1 + DZI(MZ+ 625) * XY2
     1              + DZI(MZ+  1250) * XY3 + DZI(MZ+  1875) * XY4
C
      XJDER = XJDER + DXJ(MX     ) * YZ1 + DXJ(MX+ 625) * YZ2
     1              + DXJ(MX+  1250) * YZ3 + DXJ(MX+  1875) * YZ4
      YJDER = YJDER + DYJ(MY     ) * XZ1 + DYJ(MY+ 625) * XZ2
     1              + DYJ(MY+  1250) * XZ3 + DYJ(MY+  1875) * XZ4
      ZJDER = ZJDER + DZJ(MZ     ) * XY1 + DZJ(MZ+ 625) * XY2
     1              + DZJ(MZ+  1250) * XY3 + DZJ(MZ+  1875) * XY4
C
  104 CONTINUE
      RETURN
C
  500 CONTINUE
      DO 105 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 105 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DENSTY=ABDENS(N)
      IF(DABS(DENSTY).LT.DLIM) GO TO 105
      XY1 = XIN(MX     ) * YIN(MY     ) * DENSTY
      XZ1 = XIN(MX     ) * ZIN(MZ     ) * DENSTY
      YZ1 = YIN(MY     ) * ZIN(MZ     ) * DENSTY
C
      XY2 = XIN(MX+ 625) * YIN(MY+ 625) * DENSTY
      XZ2 = XIN(MX+ 625) * ZIN(MZ+ 625) * DENSTY
      YZ2 = YIN(MY+ 625) * ZIN(MZ+ 625) * DENSTY
C
      XY3 = XIN(MX+  1250) * YIN(MY+  1250) * DENSTY
      XZ3 = XIN(MX+  1250) * ZIN(MZ+  1250) * DENSTY
      YZ3 = YIN(MY+  1250) * ZIN(MZ+  1250) * DENSTY
C
      XY4 = XIN(MX+  1875) * YIN(MY+  1875) * DENSTY
      XZ4 = XIN(MX+  1875) * ZIN(MZ+  1875) * DENSTY
      YZ4 = YIN(MY+  1875) * ZIN(MZ+  1875) * DENSTY
C
      XY5 = XIN(MX+2500) * YIN(MY+2500) * DENSTY
      XZ5 = XIN(MX+2500) * ZIN(MZ+2500) * DENSTY
      YZ5 = YIN(MY+2500) * ZIN(MZ+2500) * DENSTY
C
      XIDER = XIDER + DXI(MX     ) * YZ1 + DXI(MX+ 625) * YZ2
     1              + DXI(MX+  1250) * YZ3 + DXI(MX+  1875) * YZ4
     2              + DXI(MX+2500) * YZ5
      YIDER = YIDER + DYI(MY     ) * XZ1 + DYI(MY+ 625) * XZ2
     1              + DYI(MY+  1250) * XZ3 + DYI(MY+  1875) * XZ4
     2              + DYI(MY+2500) * XZ5
      ZIDER = ZIDER + DZI(MZ     ) * XY1 + DZI(MZ+ 625) * XY2
     1              + DZI(MZ+  1250) * XY3 + DZI(MZ+  1875) * XY4
     2              + DZI(MZ+2500) * XY5
C
      XJDER = XJDER + DXJ(MX     ) * YZ1 + DXJ(MX+ 625) * YZ2
     1              + DXJ(MX+  1250) * YZ3 + DXJ(MX+  1875) * YZ4
     2              + DXJ(MX+2500) * YZ5
      YJDER = YJDER + DYJ(MY     ) * XZ1 + DYJ(MY+ 625) * XZ2
     1              + DYJ(MY+  1250) * XZ3 + DYJ(MY+  1875) * XZ4
     2              + DYJ(MY+2500) * XZ5
      ZJDER = ZJDER + DZJ(MZ     ) * XY1 + DZJ(MZ+ 625) * XY2
     1              + DZJ(MZ+  1250) * XY3 + DZJ(MZ+  1875) * XY4
     2              + DZJ(MZ+2500) * XY5
C
  105 CONTINUE
      RETURN
C
  600 CONTINUE
      DO 106 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 106 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DENSTY=ABDENS(N)
      IF(DABS(DENSTY).LT.DLIM) GO TO 106
      XY1 = XIN(MX     ) * YIN(MY     ) * DENSTY
      XZ1 = XIN(MX     ) * ZIN(MZ     ) * DENSTY
      YZ1 = YIN(MY     ) * ZIN(MZ     ) * DENSTY
C
      XY2 = XIN(MX+ 625) * YIN(MY+ 625) * DENSTY
      XZ2 = XIN(MX+ 625) * ZIN(MZ+ 625) * DENSTY
      YZ2 = YIN(MY+ 625) * ZIN(MZ+ 625) * DENSTY
C
      XY3 = XIN(MX+  1250) * YIN(MY+  1250) * DENSTY
      XZ3 = XIN(MX+  1250) * ZIN(MZ+  1250) * DENSTY
      YZ3 = YIN(MY+  1250) * ZIN(MZ+  1250) * DENSTY
C
      XY4 = XIN(MX+  1875) * YIN(MY+  1875) * DENSTY
      XZ4 = XIN(MX+  1875) * ZIN(MZ+  1875) * DENSTY
      YZ4 = YIN(MY+  1875) * ZIN(MZ+  1875) * DENSTY
C
      XY5 = XIN(MX+2500) * YIN(MY+2500) * DENSTY
      XZ5 = XIN(MX+2500) * ZIN(MZ+2500) * DENSTY
      YZ5 = YIN(MY+2500) * ZIN(MZ+2500) * DENSTY
C
      XY6 = XIN(MX+3125) * YIN(MY+3125) * DENSTY
      XZ6 = XIN(MX+3125) * ZIN(MZ+3125) * DENSTY
      YZ6 = YIN(MY+3125) * ZIN(MZ+3125) * DENSTY
C
      XIDER = XIDER + DXI(MX     ) * YZ1 + DXI(MX+ 625) * YZ2
     1              + DXI(MX+  1250) * YZ3 + DXI(MX+  1875) * YZ4
     2              + DXI(MX+2500) * YZ5  + DXI(MX+ 3125) *YZ6
      YIDER = YIDER + DYI(MY     ) * XZ1 + DYI(MY+ 625) * XZ2
     1              + DYI(MY+  1250) * XZ3 + DYI(MY+  1875) * XZ4
     2              + DYI(MY+2500) * XZ5  + DYI(MY+ 3125) * XZ6
      ZIDER = ZIDER + DZI(MZ     ) * XY1 + DZI(MZ+ 625) * XY2
     1              + DZI(MZ+  1250) * XY3 + DZI(MZ+  1875) * XY4
     2              + DZI(MZ+2500) * XY5  + DZI(MZ+ 3125) * XY6
C
      XJDER = XJDER + DXJ(MX     ) * YZ1 + DXJ(MX+ 625) * YZ2
     1              + DXJ(MX+  1250) * YZ3 + DXJ(MX+  1875) * YZ4
     2              + DXJ(MX+2500) * YZ5  + DXJ(MX+ 3125) * YZ6
      YJDER = YJDER + DYJ(MY     ) * XZ1 + DYJ(MY+ 625) * XZ2
     1              + DYJ(MY+  1250) * XZ3 + DYJ(MY+  1875) * XZ4
     2              + DYJ(MY+2500) * XZ5  + DYJ(MY+ 3125) * XZ6
      ZJDER = ZJDER + DZJ(MZ     ) * XY1 + DZJ(MZ+ 625) * XY2
     1              + DZJ(MZ+  1250) * XY3 + DZJ(MZ+  1875) * XY4
     2              + DZJ(MZ+2500) * XY5  + DZJ(MZ+ 3125) * XY6
C
  106 CONTINUE
      RETURN
C
  700 CONTINUE
      DO 107 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 107 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DENSTY=ABDENS(N)
      IF(DABS(DENSTY).LT.DLIM) GO TO 107
      XY1 = XIN(MX     ) * YIN(MY     ) * DENSTY
      XZ1 = XIN(MX     ) * ZIN(MZ     ) * DENSTY
      YZ1 = YIN(MY     ) * ZIN(MZ     ) * DENSTY
C
      XY2 = XIN(MX+ 625) * YIN(MY+ 625) * DENSTY
      XZ2 = XIN(MX+ 625) * ZIN(MZ+ 625) * DENSTY
      YZ2 = YIN(MY+ 625) * ZIN(MZ+ 625) * DENSTY
C
      XY3 = XIN(MX+  1250) * YIN(MY+  1250) * DENSTY
      XZ3 = XIN(MX+  1250) * ZIN(MZ+  1250) * DENSTY
      YZ3 = YIN(MY+  1250) * ZIN(MZ+  1250) * DENSTY
C
      XY4 = XIN(MX+  1875) * YIN(MY+  1875) * DENSTY
      XZ4 = XIN(MX+  1875) * ZIN(MZ+  1875) * DENSTY
      YZ4 = YIN(MY+  1875) * ZIN(MZ+  1875) * DENSTY
C
      XY5 = XIN(MX+2500) * YIN(MY+2500) * DENSTY
      XZ5 = XIN(MX+2500) * ZIN(MZ+2500) * DENSTY
      YZ5 = YIN(MY+2500) * ZIN(MZ+2500) * DENSTY
C
      XY6 = XIN(MX+3125) * YIN(MY+3125) * DENSTY
      XZ6 = XIN(MX+3125) * ZIN(MZ+3125) * DENSTY
      YZ6 = YIN(MY+3125) * ZIN(MZ+3125) * DENSTY
C
      XY7 = XIN(MX+3750) * YIN(MY+3750) * DENSTY
      XZ7 = XIN(MX+3750) * ZIN(MZ+3750) * DENSTY
      YZ7 = YIN(MY+3750) * ZIN(MZ+3750) * DENSTY
C
      XIDER = XIDER + DXI(MX     ) * YZ1 + DXI(MX+ 625) * YZ2
     1              + DXI(MX+  1250) * YZ3 + DXI(MX+  1875) * YZ4
     2              + DXI(MX+2500) * YZ5  + DXI(MX+ 3125) *YZ6
     3              + DXI(MX+3750) * YZ7
      YIDER = YIDER + DYI(MY     ) * XZ1 + DYI(MY+ 625) * XZ2
     1              + DYI(MY+  1250) * XZ3 + DYI(MY+  1875) * XZ4
     2              + DYI(MY+2500) * XZ5  + DYI(MY+ 3125) * XZ6
     3              + DYI(MY+3750) * XZ7
      ZIDER = ZIDER + DZI(MZ     ) * XY1 + DZI(MZ+ 625) * XY2
     1              + DZI(MZ+  1250) * XY3 + DZI(MZ+  1875) * XY4
     2              + DZI(MZ+2500) * XY5  + DZI(MZ+ 3125) * XY6
     3              + DZI(MZ+3750) * XY7
C
      XJDER = XJDER + DXJ(MX     ) * YZ1 + DXJ(MX+ 625) * YZ2
     1              + DXJ(MX+  1250) * YZ3 + DXJ(MX+  1875) * YZ4
     2              + DXJ(MX+2500) * YZ5  + DXJ(MX+ 3125) * YZ6
     3              + DXJ(MX+3750) * YZ7
      YJDER = YJDER + DYJ(MY     ) * XZ1 + DYJ(MY+ 625) * XZ2
     1              + DYJ(MY+  1250) * XZ3 + DYJ(MY+  1875) * XZ4
     2              + DYJ(MY+2500) * XZ5  + DYJ(MY+ 3125) * XZ6
     3              + DYJ(MY+3750) * XZ7
      ZJDER = ZJDER + DZJ(MZ     ) * XY1 + DZJ(MZ+ 625) * XY2
     1              + DZJ(MZ+  1250) * XY3 + DZJ(MZ+  1875) * XY4
     2              + DZJ(MZ+2500) * XY5  + DZJ(MZ+ 3125) * XY6
     3              + DZJ(MZ+3750) * XY7
C
  107 CONTINUE
      RETURN
C
  800 CONTINUE
      DO 108 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 108 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DENSTY=ABDENS(N)
      IF(DABS(DENSTY).LT.DLIM) GO TO 108
      XY1 = XIN(MX     ) * YIN(MY     ) * DENSTY
      XZ1 = XIN(MX     ) * ZIN(MZ     ) * DENSTY
      YZ1 = YIN(MY     ) * ZIN(MZ     ) * DENSTY
C
      XY2 = XIN(MX+ 625) * YIN(MY+ 625) * DENSTY
      XZ2 = XIN(MX+ 625) * ZIN(MZ+ 625) * DENSTY
      YZ2 = YIN(MY+ 625) * ZIN(MZ+ 625) * DENSTY
C
      XY3 = XIN(MX+ 1250) * YIN(MY+ 1250) * DENSTY
      XZ3 = XIN(MX+ 1250) * ZIN(MZ+ 1250) * DENSTY
      YZ3 = YIN(MY+ 1250) * ZIN(MZ+ 1250) * DENSTY
C
      XY4 = XIN(MX+ 1875) * YIN(MY+ 1875) * DENSTY
      XZ4 = XIN(MX+ 1875) * ZIN(MZ+ 1875) * DENSTY
      YZ4 = YIN(MY+ 1875) * ZIN(MZ+ 1875) * DENSTY
C
      XY5 = XIN(MX+2500) * YIN(MY+2500) * DENSTY
      XZ5 = XIN(MX+2500) * ZIN(MZ+2500) * DENSTY
      YZ5 = YIN(MY+2500) * ZIN(MZ+2500) * DENSTY
C
      XY6 = XIN(MX+3125) * YIN(MY+3125) * DENSTY
      XZ6 = XIN(MX+3125) * ZIN(MZ+3125) * DENSTY
      YZ6 = YIN(MY+3125) * ZIN(MZ+3125) * DENSTY
C
      XY7 = XIN(MX+3750) * YIN(MY+3750) * DENSTY
      XZ7 = XIN(MX+3750) * ZIN(MZ+3750) * DENSTY
      YZ7 = YIN(MY+3750) * ZIN(MZ+3750) * DENSTY
C
      XY8 = XIN(MX+4375) * YIN(MY+4375) * DENSTY
      XZ8 = XIN(MX+4375) * ZIN(MZ+4375) * DENSTY
      YZ8 = YIN(MY+4375) * ZIN(MZ+4375) * DENSTY
C
      XIDER = XIDER + DXI(MX     ) * YZ1 + DXI(MX+ 625) * YZ2
     1              + DXI(MX+ 1250) * YZ3 + DXI(MX+ 1875) * YZ4
     2              + DXI(MX+2500) * YZ5  + DXI(MX+ 3125) * YZ6
     3              + DXI(MX+3750) * YZ7  + DXI(MX+ 4375) * YZ8
      YIDER = YIDER + DYI(MY     ) * XZ1 + DYI(MY+ 625) * XZ2
     1              + DYI(MY+ 1250) * XZ3 + DYI(MY+ 1875) * XZ4
     2              + DYI(MY+2500) * XZ5  + DYI(MY+ 3125) * XZ6
     3              + DYI(MY+3750) * XZ7  + DYI(MY+ 4375) * XZ8
      ZIDER = ZIDER + DZI(MZ     ) * XY1 + DZI(MZ+ 625) * XY2
     1              + DZI(MZ+ 1250) * XY3 + DZI(MZ+ 1875) * XY4
     2              + DZI(MZ+2500) * XY5  + DZI(MZ+ 3125) * XY6
     3              + DZI(MZ+3750) * XY7  + DZI(MZ+ 4375) * XY8
C
      XJDER = XJDER + DXJ(MX     ) * YZ1 + DXJ(MX+ 625) * YZ2
     1              + DXJ(MX+ 1250) * YZ3 + DXJ(MX+ 1875) * YZ4
     2              + DXJ(MX+2500) * YZ5  + DXJ(MX+ 3125) * YZ6
     3              + DXJ(MX+3750) * YZ7  + DXJ(MX+ 4375) * YZ8
      YJDER = YJDER + DYJ(MY     ) * XZ1 + DYJ(MY+ 625) * XZ2
     1              + DYJ(MY+ 1250) * XZ3 + DYJ(MY+ 1875) * XZ4
     2              + DYJ(MY+2500) * XZ5  + DYJ(MY+ 3125) * XZ6
     3              + DYJ(MY+3750) * XZ7  + DYJ(MY+ 4375) * XZ8
      ZJDER = ZJDER + DZJ(MZ     ) * XY1 + DZJ(MZ+ 625) * XY2
     1              + DZJ(MZ+ 1250) * XY3 + DZJ(MZ+ 1875) * XY4
     2              + DZJ(MZ+2500) * XY5  + DZJ(MZ+ 3125) * XY6
     3              + DZJ(MZ+3750) * XY7  + DZJ(MZ+ 4375) * XY8
C
  108 CONTINUE
      RETURN
C
  900 CONTINUE
      DO 109 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 109 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DENSTY=ABDENS(N)
      IF(DABS(DENSTY).LT.DLIM) GO TO 109
      XY1 = XIN(MX     ) * YIN(MY     ) * DENSTY
      XZ1 = XIN(MX     ) * ZIN(MZ     ) * DENSTY
      YZ1 = YIN(MY     ) * ZIN(MZ     ) * DENSTY
C
      XY2 = XIN(MX+ 625) * YIN(MY+ 625) * DENSTY
      XZ2 = XIN(MX+ 625) * ZIN(MZ+ 625) * DENSTY
      YZ2 = YIN(MY+ 625) * ZIN(MZ+ 625) * DENSTY
C
      XY3 = XIN(MX+ 1250) * YIN(MY+ 1250) * DENSTY
      XZ3 = XIN(MX+ 1250) * ZIN(MZ+ 1250) * DENSTY
      YZ3 = YIN(MY+ 1250) * ZIN(MZ+ 1250) * DENSTY
C
      XY4 = XIN(MX+ 1875) * YIN(MY+ 1875) * DENSTY
      XZ4 = XIN(MX+ 1875) * ZIN(MZ+ 1875) * DENSTY
      YZ4 = YIN(MY+ 1875) * ZIN(MZ+ 1875) * DENSTY
C
      XY5 = XIN(MX+2500) * YIN(MY+2500) * DENSTY
      XZ5 = XIN(MX+2500) * ZIN(MZ+2500) * DENSTY
      YZ5 = YIN(MY+2500) * ZIN(MZ+2500) * DENSTY
C
      XY6 = XIN(MX+3125) * YIN(MY+3125) * DENSTY
      XZ6 = XIN(MX+3125) * ZIN(MZ+3125) * DENSTY
      YZ6 = YIN(MY+3125) * ZIN(MZ+3125) * DENSTY
C
      XY7 = XIN(MX+3750) * YIN(MY+3750) * DENSTY
      XZ7 = XIN(MX+3750) * ZIN(MZ+3750) * DENSTY
      YZ7 = YIN(MY+3750) * ZIN(MZ+3750) * DENSTY
C
      XY8 = XIN(MX+4375) * YIN(MY+4375) * DENSTY
      XZ8 = XIN(MX+4375) * ZIN(MZ+4375) * DENSTY
      YZ8 = YIN(MY+4375) * ZIN(MZ+4375) * DENSTY
C
      XY9 = XIN(MX+5000) * YIN(MY+5000) * DENSTY
      XZ9 = XIN(MX+5000) * ZIN(MZ+5000) * DENSTY
      YZ9 = YIN(MY+5000) * ZIN(MZ+5000) * DENSTY
C
      XIDER = XIDER + DXI(MX     ) * YZ1 + DXI(MX+ 625) * YZ2
     1              + DXI(MX+ 1250) * YZ3 + DXI(MX+ 1875) * YZ4
     2              + DXI(MX+2500) * YZ5  + DXI(MX+ 3125) * YZ6
     3              + DXI(MX+3750) * YZ7  + DXI(MX+ 4375) * YZ8
     4              + DXI(MX+5000) * YZ9
      YIDER = YIDER + DYI(MY     ) * XZ1 + DYI(MY+ 625) * XZ2
     1              + DYI(MY+ 1250) * XZ3 + DYI(MY+ 1875) * XZ4
     2              + DYI(MY+2500) * XZ5  + DYI(MY+ 3125) * XZ6
     3              + DYI(MY+3750) * XZ7  + DYI(MY+ 4375) * XZ8
     4              + DYI(MY+5000) * XZ9
      ZIDER = ZIDER + DZI(MZ     ) * XY1 + DZI(MZ+ 625) * XY2
     1              + DZI(MZ+ 1250) * XY3 + DZI(MZ+ 1875) * XY4
     2              + DZI(MZ+2500) * XY5  + DZI(MZ+ 3125) * XY6
     3              + DZI(MZ+3750) * XY7  + DZI(MZ+ 4375) * XY8
     4              + DZI(MZ+5000) * XY9
C
      XJDER = XJDER + DXJ(MX     ) * YZ1 + DXJ(MX+ 625) * YZ2
     1              + DXJ(MX+ 1250) * YZ3 + DXJ(MX+ 1875) * YZ4
     2              + DXJ(MX+2500) * YZ5  + DXJ(MX+ 3125) * YZ6
     3              + DXJ(MX+3750) * YZ7  + DXJ(MX+ 4375) * YZ8
     4              + DXJ(MX+5000) * YZ9
      YJDER = YJDER + DYJ(MY     ) * XZ1 + DYJ(MY+ 625) * XZ2
     1              + DYJ(MY+ 1250) * XZ3 + DYJ(MY+ 1875) * XZ4
     2              + DYJ(MY+2500) * XZ5  + DYJ(MY+ 3125) * XZ6
     3              + DYJ(MY+3750) * XZ7  + DYJ(MY+ 4375) * XZ8
     4              + DYJ(MY+5000) * XZ9
      ZJDER = ZJDER + DZJ(MZ     ) * XY1 + DZJ(MZ+ 625) * XY2
     1              + DZJ(MZ+ 1250) * XY3 + DZJ(MZ+ 1875) * XY4
     2              + DZJ(MZ+2500) * XY5  + DZJ(MZ+ 3125) * XY6
     3              + DZJ(MZ+3750) * XY7  + DZJ(MZ+ 4375) * XY8
     4              + DZJ(MZ+5000) * XY9
C
  109 CONTINUE
      RETURN
C
      END
