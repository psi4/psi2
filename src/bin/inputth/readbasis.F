c This reads in the basis set using the libpsi input parsing utilities.
c These utilities must be properly initialized before this routine is
c called.
      subroutine intrb(c,maxcor)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include <error.h>
      common /instyl/ inold, isf9, isf8
      integer maxcor
      integer c(maxcor)

#include <basis.h>

      integer errcod, dummy
      integer mxprim

      integer intrb1, coallc, frdc, isizof, fcnt

c The number of unique centers
      errcod = fcnt('ATOMS',nucent)
      if (errcod .ne. EOK) then
        write(6,*) 'FATAL ERROR: looking for ''ATOMS'''
        call pabort
        call qabort
        endif

c     This is the maximum number of primitives that can be held on the
c     first attempt to read the basis set.
c     Use a huge number to begin with, the excess will be freed later on.
      mxprim = 1000
  801 continue
        dummy = coallc(c,mxprim*isizof('REAL'))
        pexpon = coallc(c,mxprim*isizof('REAL'))
        pconco = coallc(c,mxprim*isizof('REAL'))
        ppshc  = coallc(c,nucent*isizof('INT'))
        pnshc  = coallc(c,nucent*isizof('INT'))
        pnprim = coallc(c,mxprim*isizof('INT'))
        ppprim = coallc(c,mxprim*isizof('INT'))
        ptypes = coallc(c,mxprim*isizof('INT'))
        errcod = intrb1(c(pexpon),c(pconco),c(ppshc),c(pnprim),
     &                  c(ppprim),c(ptypes),c(pnshc),
     &                  mxprim,nshell,nprimt,nucent)
        if (errcod.ne.0) then
          call cofree(c,dummy)
          call cofree(c,pexpon)
          call cofree(c,pconco)
          call cofree(c,ppshc)
          call cofree(c,pnshc)
          call cofree(c,pnprim)
          call cofree(c,ppprim)
          call cofree(c,ptypes)
          mxprim = mxprim * 2
          goto 801
          endif

c Repack the arrays so the minimum amount of space is used.
      call cofree(c,dummy)
      call corllc(c,pexpon,nprimt,'REAL')
      call corllc(c,pconco,nprimt,'REAL')
      call corllc(c,ppshc,nucent,'INT')
      call corllc(c,pnshc,nucent,'INT')
      call corllc(c,pnprim,nshell,'INT')
      call corllc(c,ppprim,nshell,'INT')
      call corllc(c,ptypes,nshell,'INT')

c Set up ntri and nprim2
      nprim2 = nprimt**2

      return
      end
c---------------------------------------------------------------------
      integer function intrb1(expon,conco,pshc,nprim,pprim,types,nshc,
     &                        mxprim,nshell,nprimt,nucent)
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include <error.h>
      integer mxprim,nshell,nprimt,nucent
      real*8 expon(mxprim),conco(mxprim)
      integer pshc(nucent),nshc(nucent),nprim(mxprim),pprim(mxprim)
      integer types(mxprim)

      integer errcod,i,depth,just1
      character*10 junk, type
      character*80 atom,basis,keywrd,sbasis
      character*40 atoms
      character*2 symb

      common/atoml/symb(103),atoms(103)
      common /instyl/ inold, isf9, isf8

      integer inold,isf9,isf8
      integer chrlen, frdpos, frdc1, frdrpn, frdlpn, frdc
      integer frdchr, frdflt, fstck

c Look in the input and compute the number of shells and primitives so
c that the arrays can be allocated.

      nprimt = 0
      nshell = 0

cets052791 if BASIS is a scalar, then use the same basis set for each atom

      errcod = fcnt('BASIS',just1)
      if(just1.lt.0) then
        errcod = frdc('BASIS',sbasis)
        if (errcod .ne. EOK) then
          write(6,*) 'FATAL ERROR: could not read ''BASIS'' vector'
          call pabort
          call qabort
        endif
      endif

c Loop thru the unique atoms
      do 101 i = 1, nucent
        pshc(i) = nshell + 1
        nshc(i) = 0
        errcod = frdc1('ATOMS',atom,i)
        if (errcod .ne. EOK) then
          write(6,*) 'FATAL ERROR: could not read ''ATOMS'' vector'
          call pabort
          call qabort
          endif
        if(just1.gt.0) then
          errcod = frdc1('BASIS',basis,i)
          if (errcod .ne. EOK) then
            write(6,*) 'FATAL ERROR: could not read ''BASIS'' vector'
            call pabort
            call qabort
            endif
        else
          basis=sbasis
          endif

cets052791 addition for atomic symbols

        do 201 j=1,58
  201     if(atom.eq.atoms(j).or.atom.eq.symb(j)) goto 202
  202   atom=atoms(j)

c Set up the keyword path to the value needed to get the basis sets.
croag	added two more files to the working path.
croag	note that only basis sections are read, not default.
        call fwkclr()
        errcod = fwkadd(5,'BASIS')
	if (isf9.eq.1) errcod = fwkadd(9,'BASIS')
	if (isf8.eq.1) errcod = fwkadd(8,'BASIS')
        errcod = fwkadd(7,'BASIS')

        keywrd = ' '
        call scat(keywrd,atom)
        call scat(keywrd,':')
        call scat(keywrd,basis)
        errcod = frdpos(keywrd)
        if (errcod.ne.EOK) then
          write(6,*) 'FATAL ERROR: could not find this basis set:'
          write(6,*) keywrd(1:chrlen(keywrd))
          write(6,*) 'The working keyword path:'
          call fwkswk()
          rewind 7
          keywrd = ' '
          read(7,'(A)') keywrd
          write(6,'(A)') keywrd
          read(7,'(A)') keywrd
          write(6,'(A)') keywrd
          read(7,'(A)') keywrd
          write(6,'(A)') keywrd
          read(7,'(A)') keywrd
          write(6,'(A)') keywrd
          read(7,'(A)') keywrd
          write(6,'(A)') keywrd
          read(7,'(A)') keywrd
          write(6,'(A)') keywrd
          read(7,'(A)') keywrd
          write(6,'(A)') keywrd
          call pabort
          call qabort
          endif

c Begin reading the basis set.  May need to recursively read the
c data since a basis set may have a 'get' in place of a primitive.
        depth = 0

        errcod = frdlpn()
        if (errcod .ne. EOK) goto 901

c         Loop over the shells (and get's)
  801   continue

          errcod = frdchr(junk)
          if (errcod .eq. Eevect) then
c           We are done with this set of shells
            if (depth.eq.0) goto 102
            depth = depth - 1
            errcod = fstck('POP')
            if (errcod.ne.EOK) goto 901
            goto 801
          else if (errcod .eq. Ebvect) then
c           Start a new shell
            nshell = nshell + 1
            nshc(i) = nshc(i) + 1
            if (nshell.gt.mxprim) then
              intrb1 = 1
              return
              endif

c           Update the information for this shell
            nprim(nshell) = 0
            pprim(nshell) = nprimt + 1

 803        continue

            errcod = frdchr(type)
            if (errcod .ne. EOK) goto 901

            if (type.eq.'GET') then
              errcod = frdchr(basis)
              if (errcod.ne.EOK) goto 901
              errcod = frdrpn()
              if (errcod.ne.EOK) goto 901
              depth = depth + 1
              errcod = fstck('PUSH')
              if (errcod.ne.EOK) goto 901
              keywrd = ':BASIS:'
              call scat(keywrd,atom)
              call scat(keywrd,':')
              call scat(keywrd,basis)
              errcod = frdpos(keywrd)
              if (errcod.ne.EOK) goto 901
              errcod = frdlpn()
              if (errcod.ne.EOK) goto 901
              errcod = frdlpn()
              if (errcod.ne.EOK) goto 901
              goto 803
              endif

            call inttyp('C_TO_I',type,types(nshell))

 802        continue
c           Now get the primitives
            errcod = frdchr(junk)
c           If we have Eevect then this shell is done; go to next shell
            if (errcod.eq.Eevect) goto 801
            if (errcod.ne.Ebvect) goto 901
            nprimt = nprimt + 1
            nprim(nshell) = nprim(nshell) + 1
            if (nprimt .gt. mxprim) then
              intrb1 = 1
              return
              endif
            errcod = frdflt(expon(nprimt))
            if (errcod.ne.EOK) goto 901
            errcod = frdflt(conco(nprimt))
            if (errcod.ne.EOK) goto 901
            errcod = frdchr(junk)
            if (errcod.ne.Eevect) goto 901
            goto 802
          else
            goto 901
            endif

  102   continue
c Reset the keyword path to the default path
        call fwkdef()

  101   continue

c   NORMAL EXIT
      intrb1 = 0
      return

c   SYNTAX ERROR IN BASIS SPECIFICATION
  901 continue
      write(6,*) 'FATAL ERROR: bad basis set specification(',errcod,')'
      call pabort
      call qabort
      return
      end   
