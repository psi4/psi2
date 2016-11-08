      subroutine rddisp(indi,ccin,ndisp)
      implicit integer(a-z)
#include <error.h>
      character*80 lbls
      real*8 ccin(4,ndisp),cff
      dimension indi(4,ndisp)

      errcod = frdpos('DISPLACEMENTS')
      if(errcod.ne.EOK) goto 901

      errcod = frdlpn()
      do 10 disno=1,ndisp
        errcod = frdlpn()
        if(errcod.ne.EOK) goto 901
        ndp=0
 201    errcod = frdint(dipno)
        if(errcod.ne.Eevect) then
          ndp=ndp+1
          if(ndp.gt.4) goto 902
          indi(ndp,disno)=dipno
          errcod = frdflt(ccin(ndp,disno))
          if(errcod.ne.EOK) goto 901
          goto 201
        endif
  10  continue
      return 

 901  write(6,*) 'error parsing DISPLACEMENTS vector'
      call pabort
      call qabort
 902  write(6,*) 'can only handle 4 coordinates per displacement'
      call pabort
      call qabort
      end
