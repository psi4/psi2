      subroutine aclear(n,a)
      real*8 a(n)
c
c  A = 0
c
      do 10 i=1,n
 10   a(i) = 0.0d0
      return
c
      end
