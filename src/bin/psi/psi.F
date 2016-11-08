      program psi
#if defined(UNDEFINED)
      implicit undefined(a-z)
#else
      implicit integer (a-z)
#endif
#include <error.h>
      integer mxexec
      parameter (mxexec = 100)
      character*80 exec(mxexec)
      character*80 wfn,derty,openty,proced
      integer i,j
      integer check
      integer nexec
      integer rdepth
      integer stability, disp, opt, large

      integer errcod, frdc, chrlen, chrend, frdboo, fwkadd, frdpos

      call psinit('TERMINAL!')

      write(6,1) 'The Psi Execution Driver'
  1   format(20X,A)

      call fioopn(7,'psi.dat','LIB','READ')
      call fwkdef()
      errcod = fwkadd(7,'PSI')
      if (errcod.ne.0) then
        write(6,*) 'INTERNAL ERROR: fwkadd failed'
        endif

c--------------------------------------------- Read in input parameters
      errcod = frdc('WFN',wfn)
      if (errcod.ne.EOK) then
        write(6,*) 'ERROR: a problem arose while reading the required'
     &           // ' string valued keyword ''WFN'''
        call qabort
        endif

      errcod = frdc('DERTYPE',derty)
      if (errcod.eq.EEOD) then
        derty = 'NONE'
      else if (errcod.ne.EOK) then
        write(6,*) 'ERROR: a problem arose while reading the optional'
     &           // ' string valued keyword ''DERTYPE'''
        call qabort
        endif

      errcod = frdc('OPENTYPE',openty)
      if (errcod.eq.EEOD) then
        openty = 'NONE'
      else if (errcod.ne.EOK) then
        write(6,*) 'ERROR: a problem arose while reading the optional'
     &           // ' string valued keyword ''OPENTYPE'''
        call qabort
        endif

      errcod = frdboo('CHECK',check)
      if (errcod.eq.EEOD) then
        check = 0
      else if (errcod.ne.EOK) then
        write(6,*) 'ERROR: a problem arose while reading the optional'
     &           // ' boolean valued keyword ''CHECK'''
        call pabort
        call qabort
        endif

      errcod = frdboo('DISP',disp)
      if (errcod.eq.EEOD) then
        disp = 0
      else if (errcod.ne.EOK) then
        write(6,*) 'ERROR: a problem arose while reading the optional'
     &           // ' boolean valued keyword ''DISP'''
        call pabort
        call qabort
        endif

      errcod = frdboo('OPT',opt)
      if (errcod.eq.EEOD) then
        opt = 0
      else if (errcod.ne.EOK) then
        write(6,*) 'ERROR: a problem arose while reading the optional'
     &           // ' boolean valued keyword ''OPT'''
        call pabort
        call qabort
        endif

      errcod = frdboo('LARGE',large)
      if (errcod.eq.EEOD) then
        large = 0
      else if (errcod.ne.EOK) then
        write(6,*) 'ERROR: a problem arose while reading the optional'
     &           // ' boolean valued keyword ''LARGE'''
        call pabort
        call qabort
        endif

      errcod = frdboo('STABILITY',stability)
      if (errcod.eq.EEOD) then
        stability = 0
      else if (errcod.ne.EOK) then
        write(6,*) 'ERROR: a problem arose while reading the optional'
     &           // ' boolean valued keyword ''STABILITY'''
        call pabort
        call qabort
        endif

      write(6,*) 'WFN       = ',wfn(1:chrlen(wfn))
      write(6,*) 'DERTYPE   = ',derty(1:chrlen(derty))
      write(6,*) 'OPENTYPE  = ',openty(1:chrlen(openty))
      if (check.ne.0) write(6,*) 'CHECK     = YES'
      if (disp.ne.0) write(6,*)  'DISP      = YES'
      if (opt.ne.0) write(6,*)   'OPT       = YES'
      if (large.ne.0) write(6,*) 'LARGE     = YES'
      if (stability.ne.0) write(6,*) 'STABILITY = YES'

c ---------------------------------
c The basic restrictions are checked for here, more may be caught below.

      if (     (openty.ne.'NONE')
     &    .and.(openty.ne.'HIGHSPIN')
     &    .and.(openty.ne.'TWOCON')
     &    .and.(openty.ne.'SPECIAL')
     &    .and.(openty.ne.'SINGLET')) then
        write(6,*) 'ERROR: bad ''OPENTYPE'': must be one of'
        write(6,*) ' NONE HIGHSPIN TWOCON SINGLET SPECIAL'
        write(6,*) 'Further restrictions may apply, depending on'
     &              // ' other option settings.'
        call qabort
        endif

      if (     (derty.ne.'NONE')
     &    .and.(derty.ne.'FIRST')
     &    .and.(derty.ne.'SECOND')
     &    .and.(derty.ne.'THIRD')
     &    .and.(derty.ne.'FOURTH')) then
        write(6,*) 'ERROR: bad ''DERTYPE'': must be one of'
        write(6,*) ' NONE FIRST SECOND THIRD FOURTH'
        write(6,*) 'Further restrictions may apply, depending on'
     &              // ' other option settings.'
        call qabort
        endif

      if ((disp.ne.0).and.(opt.ne.0)) then
        write(6,*) 'ERROR: both ''DISP'' and ''OPT'' are YES,'
        write(6,*) '  only one of these options must be active'
        call qabort
        endif

      if (     (stability.ne.0).and.(derty.ne.'SECOND')
     &    .and.(wfn.ne.'SCF')) then
        write(6,*) 'ERROR: if ''STABILITY'' is YES,'
        write(6,*) '  then ''DERTYPE'' must be SECOND'
        write(6,*) '  and ''WFN'' must be SCF'
        call qabort
        endif

      if (     (stability.ne.0).and.(openty.ne.'NONE')
     &    .and.(openty.ne.'HIGHSPIN')
     &    .and.(openty.ne.'SINGLET')) then
        write(6,*) 'ERROR: if ''STABILITY'' is YES,'
        write(6,*) '  then ''OPENTYPE'' must be one of'
        write(6,*) ' NONE HIGHSPIN SINGLET'
        call qabort
        endif


c ---------------------------------
c  See if a user defined procedure is to be executed
      errcod = frdpos('EXEC')
      if (errcod.eq.EEOD) then
c     Construct the name of the procedure from the derty, openty, and wfn
        proced = ' '
        if (large.ne.0) then
          call scat(proced,'LARGE')
          endif
        call scat(proced,wfn)
        if (openty.eq.'NONE') then
          call scat(proced,'CLOSEDSHELL')
        else
          call scat(proced,openty)
          endif
        if (derty.eq.'NONE') then
          call scat(proced,'ENERGY')
        else
          call scat(proced,derty)
          endif
        if (disp.ne.0) then
          call scat(proced,'DISP')
        else if (opt.ne.0) then
          call scat(proced,'OPT')
          endif
        if (stability.ne.0) then
          call scat(proced,'STABILITY')
          endif
        call lexec(exec,nexec,mxexec,proced)
      else
c     Use the user provided execution list
        write(6,*) 'Using the user provided execution list from'
     &           // ' ''EXEC'''
        call lexec(exec,nexec,mxexec,'EXEC')
        endif

      write(6,*)
      if (check.ne.0) then
        write(6,*) '''CHECK'' is YES, so nothing will be executed.'
        write(6,*) 'The following programs would otherwise be executed:'
      else
        write(6,*) 'The following programs will be executed:'
        endif

      rdepth = 0
      do 101 i = 1, nexec
        if (exec(i) .eq. 'END') rdepth = rdepth - 1
        write(6,*) ('  ',j=1,rdepth),exec(i)(1:chrend(exec(i)))
        if (exec(i)(1:6).eq.'REPEAT') rdepth = rdepth + 1
 101    continue

      write(6,*)
      if (check.eq.0) call execut(exec,nexec)

      call psdone()
      stop
      end
