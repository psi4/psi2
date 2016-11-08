      subroutine stvint
c     ----- gauss-hermite quadrature using minimum-point formula -----
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      common/stv/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,ni,nj
      parameter ( zero=0.0d0, one=1.0d0 )
      dimension h(21), w(21), ilo(6), ihi(6)
      data h /0.0d0, -0.70710678118654752d0, 0.70710678118654752d0,
     1        -1.2247448713915890d0, 0.0d0, 1.2247448713915890d0,
     2        -1.6506801238857846d0, -0.52464762327529032d0,
     3        0.52464762327529032d0, 1.6506801238857846d0,
     4        -2.0201828704560856d0, -0.95857246461381851d0,
     5        0.0d0, 0.95857246461381851d0, 2.0201828704560856d0,
     6        -2.3506049736744922d0, -1.3358490740136969d0,
     7        -0.43607741192761651d0, 0.43607741192761651d0,
     8        1.3358490740136969d0, 2.3506049736744922d0/
      data w /1.0d0, 0.5d0, 0.5d0, 0.16666666666666667d0,
     1        0.66666666666666667d0, 0.16666666666666667d0,
     2        0.045875854768068492d0, 0.45412414523193151d0,
     3        0.45412414523193151d0, 0.045875854768068492d0,
     4        0.011257411327720689d0, 0.22207592200561264d0,
     5        0.53333333333333333d0, 0.22207592200561264d0,
     6        0.011257411327720689d0, 0.0025557844020562464d0,
     7        0.088615746041914527d0, 0.40882846955602923d0,
     8        0.40882846955602923d0, 0.088615746041914527d0,
     9        0.0025557844020562464d0/
      data ilo /1,2,4,7,11,16/
      data ihi /1,3,6,10,15,21/
c
      xint=zero
      yint=zero
      zint=zero
      npts=(ni+nj)/2
      imin=ilo(npts)
      imax=ihi(npts)
      do 13 i=imin,imax
      px=one
      py=one
      pz=one
      dum=h(i)/t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (6,5,4,3,2,1),ni
    1 px=   ax
      py=   ay
      pz=   az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 px=px*ax
      py=py*ay
      pz=pz*az
    6 go to (12,11,10,9,8,7),nj
    7 px=px*bx
      py=py*by
      pz=pz*bz
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 dum=w(i)
      xint=xint+dum*px
      yint=yint+dum*py
      zint=zint+dum*pz
   13 continue
      return
      end
