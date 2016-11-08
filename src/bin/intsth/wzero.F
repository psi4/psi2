      subroutine wzero(n,a,inca)
*mdc*on vax ibm cray
      implicit real*8 (a-h,o-z)
*mdc*off
      parameter (zero=0.0d0)
      dimension a(inca,*)
c
      do 10 i=1,n
10    a(1,i)=zero
c
      return
      end
