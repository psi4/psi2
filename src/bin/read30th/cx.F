      subroutine cx(cspdfg,mrcru,mxpr,i,iend)
      real*8 cspdfg
      dimension cspdfg(mrcru,mxpr)
      itap30=30
c        0         0         0         0         0         0         0
        call wreadw(itap30,cspdfg,intowp(mrcru*mxpr),iend,iend)
      do 876 k=1,mrcru
            WRITE (6,*) '          '
            WRITE (6,*) '    k =',k
            WRITE (6,*) '          '
        write(6,31) (cspdfg(k,j),j=1,mxpr)
  876 continue
   31 format(3f20.10)
      return
      end
