      subroutine rdsimp(stype,satom,nsim,nstre,nbnd,ntor,nlin1,nlin2,
     1 ninv,nout)
      implicit integer(a-z)
#include <error.h>
      character*4 stype(nsim)
      character*80 keywrd
      dimension satom(nsim,4)

      if(nstre.ne.0) then
        keywrd = ':INTCO:STRE'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 10 i=1,nstre
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            stype(ij)='STRE'
            errcod=frdint(satom(ij,1))
            errcod=frdint(satom(ij,2))
            errcod=frdint(dummy)
            if(errcod.ne.Eevect) errcod = frdrpn()
  10      continue
        endif
      endif

      if(ninv.ne.0) then
        keywrd = ':INTCO:INVR'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 15 i=1,ninv
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            stype(ij)='INVR'
            errcod=frdint(satom(ij,1))
            errcod=frdint(satom(ij,2))
            errcod=frdint(dummy)
            if(errcod.ne.Eevect) errcod = frdrpn()
  15      continue
        endif
      endif

      if(nbnd.ne.0) then
        keywrd = ':INTCO:BEND'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 11 i=1,nbnd
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            stype(ij)='BEND'
            errcod=frdint(satom(ij,1))
            errcod=frdint(satom(ij,2))
            errcod=frdint(satom(ij,3))
            errcod=frdint(dummy)
            if(errcod.ne.Eevect) errcod = frdrpn()
  11      continue
        endif
      endif

      if(ntor.ne.0) then
        keywrd = ':INTCO:TORS'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 12 i=1,ntor
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            stype(ij)='TORS'
            errcod=frdint(satom(ij,1))
            errcod=frdint(satom(ij,2))
            errcod=frdint(satom(ij,3))
            errcod=frdint(satom(ij,4))
            errcod=frdint(dummy)
            if(errcod.ne.Eevect) errcod = frdrpn()
  12      continue
        endif
      endif

      if(nlin1.ne.0) then
        keywrd = ':INTCO:LIN1'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 13 i=1,nlin1
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            stype(ij)='LIN1'
            errcod=frdint(satom(ij,1))
            errcod=frdint(satom(ij,2))
            errcod=frdint(satom(ij,3))
            errcod=frdint(satom(ij,4))
            errcod=frdint(dummy)
            if(errcod.ne.Eevect) errcod = frdrpn()
  13      continue
        endif
      endif

      if(nlin2.ne.0) then
        keywrd = ':INTCO:LIN2'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 16 i=1,nlin2
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            stype(ij)='LIN2'
            errcod=frdint(satom(ij,1))
            errcod=frdint(satom(ij,2))
            errcod=frdint(satom(ij,3))
            errcod=frdint(satom(ij,4))
            errcod=frdint(dummy)
            if(errcod.ne.Eevect) errcod = frdrpn()
  16      continue
        endif
      endif

      keywrd = ':INTCO:OUT'
      errcod = frdpos(keywrd)
      if(errcod.eq.EOK) then
        errcod = frdlpn()
        if(errcod.ne.EOK) goto 901
        do 14 i=1,nout
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          errcod = frdint(ij)
          stype(ij)=' OUT'
          errcod=frdint(satom(ij,1))
          errcod=frdint(satom(ij,2))
          errcod=frdint(satom(ij,3))
          errcod=frdint(satom(ij,4))
          errcod=frdint(dummy)
          if(errcod.ne.Eevect) errcod = frdrpn()
  14    continue
      endif

      return

  901 continue
      write(6,*) 'trouble parsing simple internal coordinates'
      call pabort
      call qabort
      end
