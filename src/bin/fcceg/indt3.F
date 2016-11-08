
       subroutine indt3(k1,k2,k3,ind,i1off,i2off,i3off,trsym)
       implicit integer (a-z)
       common/sheet/nm,nirred
       integer i1off(*),i2off(nm,nirred),i3off(*),trsym(*)
c      
       i1=max0(k1,k2,k3)
       i3=min0(k1,k2,k3)
       i2=k1+k2+k3-i1-i3
c    
       i3sym=trsym(i3)
       ind=i1off(i1)+i2off(i2,i3sym)+i3off(i3)
c
       return
       end
