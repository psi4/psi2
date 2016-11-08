      subroutine convgd(A,B,ireslt)
      real*8 a,b
      character*3 ireslt
c
c  A < B ?
c
      if(a.lt.b) then
          ireslt = 'YES'
          else
          ireslt = ' NO'
          endif
      return
c
      end
