      program oscc
      parameter maxcor=10500000
      real*8 cc(maxcor)
      integer ic(1)
      equivalence(cc,ic)
c
      call fentry(cc,ic,maxcor)
c    
      end
  
