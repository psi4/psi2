c Reads the geometry from the input
      subroutine intrg(c,maxcor)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      integer maxcor
      integer c(maxcor)

#include <basis.h>

      integer coallc,isizof

      pgeom = coallc(c,3*nucent*isizof('REAL'))

      call intrg1(c(pgeom),nucent)

      return
      end
c-----------------------------------------------------------------------
      subroutine intrg1(geom,nucent)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include <error.h>
      integer nucent
      real*8 geom(3,nucent)

      integer errcod, frdfa2

      errcod = frdfa2('GEOMETRY',geom,3,nucent)
      if (errcod.ne.EOK) then
        write(6,*) 'FATAL ERROR: couldn''t read ''GEOMETRY:0'''
        call pabort
        call qabort
        endif

      return
      end
c-----------------------------------------------------------------------
c This routine prints out the geometry of the unique centers.
c style tells how to print the info out, unit gives the units
      subroutine intpg(c,maxcor,style,unit)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      integer maxcor
      integer c(maxcor)
      character*(*) style,unit

#include <basis.h>

      call intpg1(c(pgeom),nucent,style,unit)

      return
      end
c-----------------------------------------------------------------------
      subroutine intpg1(geom,nucent,style,unit)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include <error.h>
      integer nucent
      real*8 geom(3,nucent)
      character*(*) style,unit

      integer i,j
      real*8 conv

      if (style.ne.'CARTESIAN') then
        write(6,*) 'ERROR: intpg1 can only print CARTESIAN geometries'
        call qabort
        endif

      if (unit.eq.'BOHR') then
        conv = 1.0d0
      else if (unit.eq.'ANGSTROM') then
        conv = 0.529177249d0
      else
        write(6,*) 'ERROR: intpg1 got an unrecognized unit'
        call qabort
        endif

      write(6,1) 'Center','x','y','z'
    1 format(a10,2x,a10,10x,2x,a10,10x,2x,a10)
      write(6,3) '    ------','--------------------',
     &           '--------------------',
     &           '--------------------'
    3 format(a10,2x,a20,2x,a20,2x,a20)

      do 101 i = 1, nucent
        write(6,2) i,(geom(j,i)*conv,j=1,3)
    2   format(i10,2x,f20.10,2x,f20.10,2x,f20.10)
 101    continue

      return
      end

