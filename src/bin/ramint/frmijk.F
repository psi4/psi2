      SUBROUTINE FRMIJK
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ABDEN/ABDENS(1296,10)
      COMMON/DINDX/IJX(100),IJY(100),IJZ(100),KLX(100),KLY(100),KLZ(100)
      COMMON/D1INT/XIDER(10),YIDER(10),ZIDER(10),
     1             XJDER(10),YJDER(10),ZJDER(10),
     2             XKDER(10),YKDER(10),ZKDER(10)
      COMMON/MAXDS/NMAXD
      COMMON/PARA4/MIJ,MKL,ISAME(100)
      COMMON/ROOTS/NONZER
      COMMON/TWOD0/XIN(2304),YIN(2304),ZIN(2304)
      COMMON/TWODI/DXI(2304),DYI(2304),DZI(2304)
      COMMON/TWODJ/DXJ(2304),DYJ(2304),DZJ(2304)
      COMMON/TWODK/DXK(2304),DYK(2304),DZK(2304)
      DATA DLIM / 1.0D-10 /
C
      N=0
      GO TO (100,200,300,400,500),NONZER
C
  100 CONTINUE
      DO 111 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 111 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DO 101 IXYZ=1,NMAXD
      DENSTY=ABDENS(N,IXYZ)
      IF(DABS(DENSTY).LT.DLIM) GO TO 101
      XY=XIN(MX) * YIN(MY) * DENSTY
      XZ=XIN(MX) * ZIN(MZ) * DENSTY
      YZ=YIN(MY) * ZIN(MZ) * DENSTY
C
      XIDER(IXYZ)=XIDER(IXYZ)+DXI(MX)*YZ
      YIDER(IXYZ)=YIDER(IXYZ)+DYI(MY)*XZ
      ZIDER(IXYZ)=ZIDER(IXYZ)+DZI(MZ)*XY
C
      XJDER(IXYZ)=XJDER(IXYZ)+DXJ(MX)*YZ
      YJDER(IXYZ)=YJDER(IXYZ)+DYJ(MY)*XZ
      ZJDER(IXYZ)=ZJDER(IXYZ)+DZJ(MZ)*XY
C
      XKDER(IXYZ)=XKDER(IXYZ)+DXK(MX)*YZ
      YKDER(IXYZ)=YKDER(IXYZ)+DYK(MY)*XZ
      ZKDER(IXYZ)=ZKDER(IXYZ)+DZK(MZ)*XY
C
  101 CONTINUE
  111 CONTINUE
      RETURN
C
  200 CONTINUE
      DO 112 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 112 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DO 102 IXYZ=1,NMAXD
      DENSTY=ABDENS(N,IXYZ)
      IF(DABS(DENSTY).LT.DLIM) GO TO 102
      XY1 = XIN(MX) * YIN(MY) * DENSTY
      XZ1 = XIN(MX) * ZIN(MZ) * DENSTY
      YZ1 = YIN(MY) * ZIN(MZ) * DENSTY
      XY2 = XIN(MX+256) * YIN(MY+256) * DENSTY
      XZ2 = XIN(MX+256) * ZIN(MZ+256) * DENSTY
      YZ2 = YIN(MY+256) * ZIN(MZ+256) * DENSTY
C
      XIDER(IXYZ)=XIDER(IXYZ)+DXI(MX)*YZ1+DXI(MX+256)*YZ2
      YIDER(IXYZ)=YIDER(IXYZ)+DYI(MY)*XZ1+DYI(MY+256)*XZ2
      ZIDER(IXYZ)=ZIDER(IXYZ)+DZI(MZ)*XY1+DZI(MZ+256)*XY2
C
      XJDER(IXYZ)=XJDER(IXYZ)+DXJ(MX)*YZ1+DXJ(MX+256)*YZ2
      YJDER(IXYZ)=YJDER(IXYZ)+DYJ(MY)*XZ1+DYJ(MY+256)*XZ2
      ZJDER(IXYZ)=ZJDER(IXYZ)+DZJ(MZ)*XY1+DZJ(MZ+256)*XY2
C
      XKDER(IXYZ)=XKDER(IXYZ)+DXK(MX)*YZ1+DXK(MX+256)*YZ2
      YKDER(IXYZ)=YKDER(IXYZ)+DYK(MY)*XZ1+DYK(MY+256)*XZ2
      ZKDER(IXYZ)=ZKDER(IXYZ)+DZK(MZ)*XY1+DZK(MZ+256)*XY2
C
  102 CONTINUE
  112 CONTINUE
      RETURN
C
  300 CONTINUE
      DO 113 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 113 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DO 103 IXYZ=1,NMAXD
      DENSTY=ABDENS(N,IXYZ)
      IF(DABS(DENSTY).LT.DLIM) GO TO 103
      XY1 = XIN(MX     ) * YIN(MY     ) * DENSTY
      XZ1 = XIN(MX     ) * ZIN(MZ     ) * DENSTY
      YZ1 = YIN(MY     ) * ZIN(MZ     ) * DENSTY
      XY2 = XIN(MX+ 256) * YIN(MY+ 256) * DENSTY
      XZ2 = XIN(MX+ 256) * ZIN(MZ+ 256) * DENSTY
      YZ2 = YIN(MY+ 256) * ZIN(MZ+ 256) * DENSTY
      XY3 = XIN(MX+ 512) * YIN(MY+ 512) * DENSTY
      XZ3 = XIN(MX+ 512) * ZIN(MZ+ 512) * DENSTY
      YZ3 = YIN(MY+ 512) * ZIN(MZ+ 512) * DENSTY
C
      XIDER(IXYZ)=XIDER(IXYZ)+DXI(MX     )*YZ1+DXI(MX+ 256)*YZ2
     1                       +DXI(MX+ 512)*YZ3
      YIDER(IXYZ)=YIDER(IXYZ)+DYI(MY     )*XZ1+DYI(MY+ 256)*XZ2
     1                       +DYI(MY+ 512)*XZ3
      ZIDER(IXYZ)=ZIDER(IXYZ)+DZI(MZ     )*XY1+DZI(MZ+ 256)*XY2
     1                       +DZI(MZ+ 512)*XY3
C
      XJDER(IXYZ)=XJDER(IXYZ)+DXJ(MX     )*YZ1+DXJ(MX+ 256)*YZ2
     1                       +DXJ(MX+ 512)*YZ3
      YJDER(IXYZ)=YJDER(IXYZ)+DYJ(MY     )*XZ1+DYJ(MY+ 256)*XZ2
     1                       +DYJ(MY+ 512)*XZ3
      ZJDER(IXYZ)=ZJDER(IXYZ)+DZJ(MZ     )*XY1+DZJ(MZ+ 256)*XY2
     1                       +DZJ(MZ+ 512)*XY3
C
      XKDER(IXYZ)=XKDER(IXYZ)+DXK(MX     )*YZ1+DXK(MX+ 256)*YZ2
     1                       +DXK(MX+ 512)*YZ3
      YKDER(IXYZ)=YKDER(IXYZ)+DYK(MY     )*XZ1+DYK(MY+ 256)*XZ2
     1                       +DYK(MY+ 512)*XZ3
      ZKDER(IXYZ)=ZKDER(IXYZ)+DZK(MZ     )*XY1+DZK(MZ+ 256)*XY2
     1                       +DZK(MZ+ 512)*XY3
C
  103 CONTINUE
  113 CONTINUE
      RETURN
C
  400 CONTINUE
      DO 114 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 114 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DO 104 IXYZ=1,NMAXD
      DENSTY=ABDENS(N,IXYZ)
      IF(DABS(DENSTY).LT.DLIM) GO TO 104
      XY1 = XIN(MX     ) * YIN(MY     ) * DENSTY
      XZ1 = XIN(MX     ) * ZIN(MZ     ) * DENSTY
      YZ1 = YIN(MY     ) * ZIN(MZ     ) * DENSTY
C
      XY2 = XIN(MX+ 256) * YIN(MY+ 256) * DENSTY
      XZ2 = XIN(MX+ 256) * ZIN(MZ+ 256) * DENSTY
      YZ2 = YIN(MY+ 256) * ZIN(MZ+ 256) * DENSTY
C
      XY3 = XIN(MX+ 512) * YIN(MY+ 512) * DENSTY
      XZ3 = XIN(MX+ 512) * ZIN(MZ+ 512) * DENSTY
      YZ3 = YIN(MY+ 512) * ZIN(MZ+ 512) * DENSTY
C
      XY4 = XIN(MX+ 768) * YIN(MY+ 768) * DENSTY
      XZ4 = XIN(MX+ 768) * ZIN(MZ+ 768) * DENSTY
      YZ4 = YIN(MY+ 768) * ZIN(MZ+ 768) * DENSTY
C
      XIDER(IXYZ)=XIDER(IXYZ)+DXI(MX     )*YZ1+DXI(MX+ 256)*YZ2
     1                       +DXI(MX+ 512)*YZ3+DXI(MX+ 768)*YZ4
      YIDER(IXYZ)=YIDER(IXYZ)+DYI(MY     )*XZ1+DYI(MY+ 256)*XZ2
     1                       +DYI(MY+ 512)*XZ3+DYI(MY+ 768)*XZ4
      ZIDER(IXYZ)=ZIDER(IXYZ)+DZI(MZ     )*XY1+DZI(MZ+ 256)*XY2
     1                       +DZI(MZ+ 512)*XY3+DZI(MZ+ 768)*XY4
C
      XJDER(IXYZ)=XJDER(IXYZ)+DXJ(MX     )*YZ1+DXJ(MX+ 256)*YZ2
     1                       +DXJ(MX+ 512)*YZ3+DXJ(MX+ 768)*YZ4
      YJDER(IXYZ)=YJDER(IXYZ)+DYJ(MY     )*XZ1+DYJ(MY+ 256)*XZ2
     1                       +DYJ(MY+ 512)*XZ3+DYJ(MY+ 768)*XZ4
      ZJDER(IXYZ)=ZJDER(IXYZ)+DZJ(MZ     )*XY1+DZJ(MZ+ 256)*XY2
     1                       +DZJ(MZ+ 512)*XY3+DZJ(MZ+ 768)*XY4
C
      XKDER(IXYZ)=XKDER(IXYZ)+DXK(MX     )*YZ1+DXK(MX+ 256)*YZ2
     1                       +DXK(MX+ 512)*YZ3+DXK(MX+ 768)*YZ4
      YKDER(IXYZ)=YKDER(IXYZ)+DYK(MY     )*XZ1+DYK(MY+ 256)*XZ2
     1                       +DYK(MY+ 512)*XZ3+DYK(MY+ 768)*XZ4
      ZKDER(IXYZ)=ZKDER(IXYZ)+DZK(MZ     )*XY1+DZK(MZ+ 256)*XY2
     1                       +DZK(MZ+ 512)*XY3+DZK(MZ+ 768)*XY4
C
  104 CONTINUE
  114 CONTINUE
      RETURN
C
  500 CONTINUE
      DO 115 I=1,MIJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DO 115 K=1,MKL
      MX=NX+KLX(K)
      MY=NY+KLY(K)
      MZ=NZ+KLZ(K)
      N=N+1
      DO 105 IXYZ=1,NMAXD
      DENSTY=ABDENS(N,IXYZ)
      IF(DABS(DENSTY).LT.DLIM) GO TO 105
      XY1 = XIN(MX     ) * YIN(MY     ) * DENSTY
      XZ1 = XIN(MX     ) * ZIN(MZ     ) * DENSTY
      YZ1 = YIN(MY     ) * ZIN(MZ     ) * DENSTY
C
      XY2 = XIN(MX+ 256) * YIN(MY+ 256) * DENSTY
      XZ2 = XIN(MX+ 256) * ZIN(MZ+ 256) * DENSTY
      YZ2 = YIN(MY+ 256) * ZIN(MZ+ 256) * DENSTY
C
      XY3 = XIN(MX+ 512) * YIN(MY+ 512) * DENSTY
      XZ3 = XIN(MX+ 512) * ZIN(MZ+ 512) * DENSTY
      YZ3 = YIN(MY+ 512) * ZIN(MZ+ 512) * DENSTY
C
      XY4 = XIN(MX+ 768) * YIN(MY+ 768) * DENSTY
      XZ4 = XIN(MX+ 768) * ZIN(MZ+ 768) * DENSTY
      YZ4 = YIN(MY+ 768) * ZIN(MZ+ 768) * DENSTY
C
      XY5 = XIN(MX+1024) * YIN(MY+1024) * DENSTY
      XZ5 = XIN(MX+1024) * ZIN(MZ+1024) * DENSTY
      YZ5 = YIN(MY+1024) * ZIN(MZ+1024) * DENSTY
C
      XIDER(IXYZ)=XIDER(IXYZ)+DXI(MX     )*YZ1+DXI(MX+ 256)*YZ2
     1                       +DXI(MX+ 512)*YZ3+DXI(MX+ 768)*YZ4
     2                       +DXI(MX+1024)*YZ5
      YIDER(IXYZ)=YIDER(IXYZ)+DYI(MY     )*XZ1+DYI(MY+ 256)*XZ2
     1                       +DYI(MY+ 512)*XZ3+DYI(MY+ 768)*XZ4
     2                       +DYI(MY+1024)*XZ5
      ZIDER(IXYZ)=ZIDER(IXYZ)+DZI(MZ     )*XY1+DZI(MZ+ 256)*XY2
     1                       +DZI(MZ+ 512)*XY3+DZI(MZ+ 768)*XY4
     2                       +DZI(MZ+1024)*XY5
C
      XJDER(IXYZ)=XJDER(IXYZ)+DXJ(MX     )*YZ1+DXJ(MX+ 256)*YZ2
     1                       +DXJ(MX+ 512)*YZ3+DXJ(MX+ 768)*YZ4
     2                       +DXJ(MX+1024)*YZ5
      YJDER(IXYZ)=YJDER(IXYZ)+DYJ(MY     )*XZ1+DYJ(MY+ 256)*XZ2
     1                       +DYJ(MY+ 512)*XZ3+DYJ(MY+ 768)*XZ4
     2                       +DYJ(MY+1024)*XZ5
      ZJDER(IXYZ)=ZJDER(IXYZ)+DZJ(MZ     )*XY1+DZJ(MZ+ 256)*XY2
     1                       +DZJ(MZ+ 512)*XY3+DZJ(MZ+ 768)*XY4
     2                       +DZJ(MZ+1024)*XY5
C
      XKDER(IXYZ)=XKDER(IXYZ)+DXK(MX     )*YZ1+DXK(MX+ 256)*YZ2
     1                       +DXK(MX+ 512)*YZ3+DXK(MX+ 768)*YZ4
     2                       +DXK(MX+1024)*YZ5
      YKDER(IXYZ)=YKDER(IXYZ)+DYK(MY     )*XZ1+DYK(MY+ 256)*XZ2
     1                       +DYK(MY+ 512)*XZ3+DYK(MY+ 768)*XZ4
     2                       +DYK(MY+1024)*XZ5
      ZKDER(IXYZ)=ZKDER(IXYZ)+DZK(MZ     )*XY1+DZK(MZ+ 256)*XY2
     1                       +DZK(MZ+ 512)*XY3+DZK(MZ+ 768)*XY4
     2                       +DZK(MZ+1024)*XY5
C
  105 CONTINUE
  115 CONTINUE
      RETURN
C
      END
