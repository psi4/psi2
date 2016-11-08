      subroutine forms(gout)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
c     ----- form integrals over functions -----
ctph  common/root/xx,u(9),w(9),nroots
      common/root/xx,u(10),w(10),nroots
      common/xyz/xin(5625),yin(5625),zin(5625)
      common /index/ ijgt(225), ijx(225), ijy(225), ijz(225), ik(225),
     1               klgt(225), klx(225), kly(225), klz(225)
      common /shlnos/lit,ljt,lkt,llt,ij,ijkl
      dimension gout(*)
      go to (500,1500,2500,3500,4500,5500,6500,7500,8500),nroots
  500 do 1000 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      n1=ijgt(i)
      kmax=ik(i)
      do 1000 k=1,kmax
      mx=nx+klx(k)
      my=ny+kly(k)
      mz=nz+klz(k)
      n=n1+klgt(k)
 1000 gout(n)=(       xin(mx     )*yin(my     )*zin(mz     ))+gout(n)
      return
 1500 do 2000 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      n1=ijgt(i)
      kmax=ik(i)
      do 2000 k=1,kmax
      mx=nx+klx(k)
      my=ny+kly(k)
      mz=nz+klz(k)
      n=n1+klgt(k)
 2000 gout(n)=(       xin(mx     )*yin(my     )*zin(mz     )
     1               +xin(mx+ 625)*yin(my+ 625)*zin(mz+ 625))+gout(n)
      return
 2500 do 3000 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      n1=ijgt(i)
      kmax=ik(i)
      do 3000 k=1,kmax
      mx=nx+klx(k)
      my=ny+kly(k)
      mz=nz+klz(k)
      n=n1+klgt(k)
 3000 gout(n)=(       xin(mx     )*yin(my     )*zin(mz     )
     1               +xin(mx+ 625)*yin(my+ 625)*zin(mz+ 625)
     1               +xin(mx+1250)*yin(my+1250)*zin(mz+1250))+gout(n)
      return
 3500 do 4000 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      n1=ijgt(i)
      kmax=ik(i)
      do 4000 k=1,kmax
      mx=nx+klx(k)
      my=ny+kly(k)
      mz=nz+klz(k)
      n=n1+klgt(k)
 4000 gout(n)=(       xin(mx     )*yin(my     )*zin(mz     )
     1               +xin(mx+ 625)*yin(my+ 625)*zin(mz+ 625)
     1               +xin(mx+1250)*yin(my+1250)*zin(mz+1250)
     1               +xin(mx+1875)*yin(my+1875)*zin(mz+1875))+gout(n)
      return
 4500 do 5000 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      n1=ijgt(i)
      kmax=ik(i)
      do 5000 k=1,kmax
      mx=nx+klx(k)
      my=ny+kly(k)
      mz=nz+klz(k)
      n=n1+klgt(k)
 5000 gout(n)=(       xin(mx     )*yin(my     )*zin(mz     )
     1               +xin(mx+ 625)*yin(my+ 625)*zin(mz+ 625)
     1               +xin(mx+1250)*yin(my+1250)*zin(mz+1250)
     1               +xin(mx+1875)*yin(my+1875)*zin(mz+1875)
     1               +xin(mx+2500)*yin(my+2500)*zin(mz+2500))+gout(n)
      return
 5500 do 6000 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      n1=ijgt(i)
      kmax=ik(i)
      do 6000 k=1,kmax
      mx=nx+klx(k)
      my=ny+kly(k)
      mz=nz+klz(k)
      n=n1+klgt(k)
 6000 gout(n)=(       xin(mx     )*yin(my     )*zin(mz     )
     1               +xin(mx+ 625)*yin(my+ 625)*zin(mz+ 625)
     1               +xin(mx+1250)*yin(my+1250)*zin(mz+1250)
     1               +xin(mx+1875)*yin(my+1875)*zin(mz+1875)
     1               +xin(mx+2500)*yin(my+2500)*zin(mz+2500)
     1               +xin(mx+3125)*yin(my+3125)*zin(mz+3125))+gout(n)
      return
 6500 do 7000 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      n1=ijgt(i)
      kmax=ik(i)
      do 7000 k=1,kmax
      mx=nx+klx(k)
      my=ny+kly(k)
      mz=nz+klz(k)
      n=n1+klgt(k)
 7000 gout(n)=(       xin(mx     )*yin(my     )*zin(mz     )
     1               +xin(mx+ 625)*yin(my+ 625)*zin(mz+ 625)
     1               +xin(mx+1250)*yin(my+1250)*zin(mz+1250)
     1               +xin(mx+1875)*yin(my+1875)*zin(mz+1875)
     1               +xin(mx+2500)*yin(my+2500)*zin(mz+2500)
     1               +xin(mx+3125)*yin(my+3125)*zin(mz+3125)
     1               +xin(mx+3750)*yin(my+3750)*zin(mz+3750))+gout(n)
      return
 7500 do 8000 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      n1=ijgt(i)
      kmax=ik(i)
      do 8000 k=1,kmax
      mx=nx+klx(k)
      my=ny+kly(k)
      mz=nz+klz(k)
      n=n1+klgt(k)
 8000 gout(n)=(       xin(mx     )*yin(my     )*zin(mz     )
     1               +xin(mx+ 625)*yin(my+ 625)*zin(mz+ 625)
     1               +xin(mx+1250)*yin(my+1250)*zin(mz+1250)
     1               +xin(mx+1875)*yin(my+1875)*zin(mz+1875)
     1               +xin(mx+2500)*yin(my+2500)*zin(mz+2500)
     1               +xin(mx+3125)*yin(my+3125)*zin(mz+3125)
     1               +xin(mx+3750)*yin(my+3750)*zin(mz+3750)
     1               +xin(mx+4375)*yin(my+4375)*zin(mz+4375))+gout(n)
      return
 8500 do 9000 i=1,ij
      nx=ijx(i)
      ny=ijy(i)
      nz=ijz(i)
      n1=ijgt(i)
      kmax=ik(i)
      do 9000 k=1,kmax
      mx=nx+klx(k)
      my=ny+kly(k)
      mz=nz+klz(k)
      n=n1+klgt(k)
 9000 gout(n)=(       xin(mx     )*yin(my     )*zin(mz     )
     1               +xin(mx+ 625)*yin(my+ 625)*zin(mz+ 625)
     1               +xin(mx+1250)*yin(my+1250)*zin(mz+1250)
     1               +xin(mx+1875)*yin(my+1875)*zin(mz+1875)
     1               +xin(mx+2500)*yin(my+2500)*zin(mz+2500)
     1               +xin(mx+3125)*yin(my+3125)*zin(mz+3125)
     1               +xin(mx+3750)*yin(my+3750)*zin(mz+3750)
     1               +xin(mx+4375)*yin(my+4375)*zin(mz+4375)
     1               +xin(mx+5000)*yin(my+5000)*zin(mz+5000))+gout(n)
      return
      end
