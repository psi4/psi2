      subroutine bottom(iunit)
c  Skip to end of file 'iunit'
      character*1 dummy
1000  format(a1)

10    read(iunit,1000,end=20) dummy
      go to 10
20    continue
      return
      end
