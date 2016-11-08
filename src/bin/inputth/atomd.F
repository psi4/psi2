c Arguments:
c   r = real*8 geometry (3 by mxat)
c   mxat = maximum number of atoms
c   nat = the actual number of atoms
c   dist = a real*8 scratch array (size nat*(nat+1)/2)
c   scale = a scaling factor for the distances
c           scale = 1.d0 will output the distances in bohr
      subroutine atomd(r,mxat,nat,dist,scale)
#if defined(UNDEFINED)
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      integer mxat,nat
      real*8 r(3,mxat),scale
      real*8 dist(*)

      integer i,j,ij

      ij = 0
      do 101 i = 1, nat
        do 201 j = 1, i
          ij = ij + 1
          dist(ij) = dsqrt(  (r(1,i)-r(1,j))**2
     &                     + (r(2,i)-r(2,j))**2
     &                     + (r(3,i)-r(3,j))**2
     &                     )*scale
  201     continue
  101   continue

      call print(dist,ij,nat,6)
      return
      end
