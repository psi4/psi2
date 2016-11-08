      subroutine rdsymm(symms,symmc,nsym,nasym,lbls)
      implicit integer(a-z)
#include <error.h>
      character*80 lbls(nsym)
      real*8 symmc(nsym,100),cff
      dimension symms(nsym,0:100)

      errcod = frdpos(':INTCO:SYMM')
      if(errcod.ne.EOK) goto 901

      errcod = frdlpn()
      do 10 symco=1,nsym-nasym
        errcod = frdlpn()
        errcod = frdchr(lbls(symco))
        if(errcod.ne.EOK) goto 901
        errcod = frdlpn()
        if(errcod.ne.EOK) goto 901
        nsimp=0
 201    errcod = frdint(simno)
        if(errcod.ne.Eevect) then
          nsimp=nsimp+1
          if(simno.lt.0) then
            symms(symco,nsimp) = -simno
            symmc(symco,nsimp) = -1.0d0
          else
            symms(symco,nsimp) = simno
            symmc(symco,nsimp) = 1.0d0
          endif
          goto 201
        endif
        errcod = frdlpn()
        if(errcod.eq.EOK) then
          do 100 sim=1,nsimp
            errcod = frdflt(cff)
            symmc(symco,sim)=symmc(symco,sim)*cff
 100      continue
          errcod=frdrpn()
          errcod=frdrpn()
        endif
        symms(symco,0)=nsimp
  10  continue

      if (nasym.eq.0) return

      errcod = frdpos(':INTCO:ASYMM')

      errcod = frdlpn()
      do 20 symco=nsym-nasym+1,nsym
        errcod = frdlpn()
        errcod = frdchr(lbls(symco))
        if(errcod.ne.EOK) goto 902
        errcod = frdlpn()
        if(errcod.ne.EOK) goto 902
        nsimp=0
 221    errcod = frdint(simno)
        if(errcod.ne.Eevect) then
          nsimp=nsimp+1
          if(simno.lt.0) then
            symms(symco,nsimp) = -simno
            symmc(symco,nsimp) = -1.0d0
          else
            symms(symco,nsimp) = simno
            symmc(symco,nsimp) = 1.0d0
          endif
          goto 221
        endif
        errcod = frdlpn()
        if(errcod.eq.EOK) then
          do 200 sim=1,nsimp
            errcod = frdflt(cff)
            symmc(symco,sim)=symmc(symco,sim)*cff
 200      continue
          errcod=frdrpn()
          errcod=frdrpn()
        endif
        symms(symco,0)=nsimp
  20  continue
      return 

 901  write(6,*) 'error parsing SYMM vector'
      call pabort
      call qabort
 902  write(6,*) 'error parsing ASYMM vector'
      call pabort
      call qabort
      end
