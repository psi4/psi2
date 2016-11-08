c Arguments:
c   r = real*8 geometry (3 by mxat)
c   mxat = maximum number of atoms
c   nat = the actual number of atoms
c   scale = a scaling factor for the distances
c           scale = 1.d0 will output the distances in radian
      subroutine bonda(r,mxat,nat,scale)
#if defined(UNDEFINED)
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      integer mxat,nat
      real*8 r(3,mxat),scale

c     if a bond length is > long the angle is not printed.
      real*8 long
      parameter(long=8.d0)

      real*8 vij1,vij2,vij3,vik1,vik2,vik3
      real*8 dot,angle,rij,rik
      integer fmt
      integer i,j,k,ij

  11  format(1X,I1,'-',I1,'-',I1,4X,F18.8)
  12  format(1X,I2,'-',I2,'-',I2,4X,F18.8)
  15  format(1X,I5,'-',I5,'-',I5,4X,F18.8)

      if (nat.lt.10) then
        assign 11 to fmt
      else if (nat.lt.100) then
        assign 12 to fmt
      else
        assign 15 to fmt
        endif

c     CENTER atom first
      do 101 i = 1, nat
        do 201 j = 1, nat
          if (j.eq.i) goto 201
          vij1 = r(1,i)-r(1,j)
          vij2 = r(2,i)-r(2,j)
          vij3 = r(3,i)-r(3,j)
          rij = dsqrt(vij1**2 + vij2**2 + vij3**2)
          if (rij.lt.0.001) goto 901
          if (rij.gt.4.d0) goto 201
          do 301 k = 1, j - 1
            if ((k.eq.i).or.(k.eq.j)) goto 301
            vik1 = r(1,i)-r(1,k)
            vik2 = r(2,i)-r(2,k)
            vik3 = r(3,i)-r(3,k)
            rik = dsqrt(vik1**2 + vik2**2 + vik3**2)
c           If the atoms are distant, then don't need the angle
            if (rik.gt.4.d0) goto 301
            if (rik.lt.0.001) goto 901

            dot = (vij1*vik1 + vij2*vik2 + vij3*vik3)/(rij*rik)
            if (dot.gt.1.0) dot = 1.0
            if (dot.lt.-1.0) dot = -1.0
            angle = dacos(dot)
            write(6,fmt) j,i,k,angle*scale
  301       continue
  201     continue
  101   continue

        return

  901 continue
      write(6,*) 'ERROR: bond length was too small in bonda'
      call qabort
      return
      end
