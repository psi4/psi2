c This routine reads in the charges from the input.
c The nucent value in integral.h must have already been initialized
c before this routine is called.
      subroutine intrc(c,maxcor)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      integer maxcor
      integer c(maxcor)

#include <basis.h>

      integer coallc,isizof

      pcharg = coallc(c,nucent*isizof('REAL'))

      call intrc1(c(pcharg),nucent)

      return
      end
c-----------------------------------------------------------------------
      subroutine intrc1(charg,nucent)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include <error.h>
      integer nucent
      real*8 charg(nucent)

      integer errcod, frdfa1

      errcod = frdfa1('CHARGES',charg,nucent)
      if(errcod.ne.EOK) call defchg(charg,nucent)
cets052791      if (errcod.ne.EOK) then
c        write(6,*) 'FATAL ERROR: couldn''t read ''CHARGES'''
c        call pabort
c        call qabort
c        endif

      return
      end
c-----------------------------------------------------------------------
c This routine prints out the charges on the unique centers.
      subroutine intpc(c,maxcor)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      integer maxcor
      integer c(maxcor)

#include <basis.h>

      call intpc1(c(pcharg),nucent)

      return
      end
c-----------------------------------------------------------------------
      subroutine intpc1(charg,nucent)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include <error.h>
      integer nucent
      real*8 charg(nucent)

      integer i

      write(6,1) 'Center','Charge'
    1 format(a10,2x,a10)
      write(6,1) '------','------'

      do 101 i = 1, nucent
        write(6,2) i,charg(i)
    2   format(i10,2x,f10.4)
 101    continue

      return
      end
c------------------------------------------------------------------
cets052791
c this subroutine puts the default charges in the charg array
      subroutine defchg(charg,nucent)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include <error.h>
      integer errcod, nucent, i, j
      integer frdc1
      real*8 charg(nucent)
      character*2 symb
      character*40 atom,atoms
      common/atoml/symb(103),atoms(103)

      do 100 i=1,nucent
        errcod = frdc1('ATOMS',atom,i)
        if (errcod.ne.0) then
          write(6,*) 'in defchg'
          write(6,*) 'ERROR: unexpected problem getting one of'
     &             // ' the ''ATOMS'''
          call pabort
          call qabort
        endif
        do 101 j=1,58
 101      if(atom.eq.atoms(j).or.atom.eq.symb(j)) goto 199
 199    continue
        if(j.eq.58) then
          write(6,*)'in defchg, unknown atom type',atom
          call pabort
          call qabort
        endif
        charg(i)= dfloat(j)
 100  continue
      return
      end
