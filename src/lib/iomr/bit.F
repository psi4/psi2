c_*BEGIN_FILE bit.F

      integer function iand(i,j)
      implicit integer*4 (a-z)
      iand=and(i,j)
      return
      end

      integer function ishft(i,j)
      implicit integer*4 (a-z)
      if (j .le. 0) then
          ishft=rshift(i,-j)
      else
          ishft=lshift(i, j)
          endif
      return
      end

      integer function ior(i,j)
      implicit integer*4 (a-z)
      ior = or(i,j)
      return
      end

      integer function ieor(i,j)
      implicit integer*4 (a-z)
      ieor = xor(i,j)
      return
      end

c_##ifdef APOLLO
c_#      integer function rshift(i,j)
c_#      implicit integer*4 (a-z)
c_#      rshift = rshft(i,j)
c_#      return
c_#      end
c_#
c_#      integer function lshift(i,j)
c_#      implicit integer*4 (a-z)
c_#      lshift = lshft(i,j)
c_#      return
c_#      end
c_##endif
