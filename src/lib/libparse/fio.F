c_*BEGIN_FILE fio.F
c The formatted io file management routines.
c-----------------------------------------------------------------------------
c  Formatted I/O file Open.
c name is a symbolic name or a filename.
c type is one of
c   USER: a user file is opened with name name
c   LIB: a library file is opened with name name
c   SUSER: a user file is opened with symbolic name name
c   (SLIB: is not permitted)
c access if one of
c   READ
c   WRITE (whether or not the file exists)
c   OLD   (only if the file exists)
c   NEW
c   APPEND
c   (If some of these access modes are not supported then other modes
c    may be quietly substituted)
c Symbolic file names are resolved by looking in the FILES input.
c Error recovery mode: returns a value from "error.h"
c
c	modifications made by croag on 5 march 1994
c	the SUSER option had an error, and I corrected both
c	that, and the general search path.
c	comments also added for readability.
      integer function fioopn(unit,name,type,access)
c_##ifdef sun
c_#      implicit undefined (a-z)
c_##else
      implicit integer (a-z)
c_##endif
c_##include "error.h"
c_*BEGIN_INCLUDE ./error.h
      integer EOK, Etype, Ekey, Elngth, Esyntx, Einflp, Efind, EEOF
      integer Etrunc, Eseg, Eparen, Eevect, Erdrct, EEOD, Earg, Ebvect
      integer Eopen, Eclose, Elimit, Emrpar
      parameter (
c                    OK, no error.
     & EOK = 0,
c                    Bad type in a formatted read.
     & Etype = 1,
c                    Could not find the keyword.
     & Ekey  = 2,
c                    Bad array length.
     & Elngth  = 3,
c                    Syntax error.
     & Esyntx  = 4,
c                    Infinite loop
     & Einflp  = 5,
c                    Efind (not used, 6 is available for replacement)
     & Efind = 6,
c                    End of File
     & EEOF = 7,
c                    Character value has been truncated
     & Etrunc = 8,
c                    Too many segments in a keyword
     & Eseg = 9,
c                    Too many levels of parenthesis
     & Eparen = 10,
c                    End of vector (not really an error)
     & Eevect = 11,
c                    Redirection depth is too deep.
     & Erdrct = 12,
c                    End of Data is like EOF but expected.
     & EEOD = 13,
c                    Bad value in argument list
     & Earg = 14,
c                    Could not open a file.
     & Eopen = 15,
c                    Beginning of vector (not really an error)
     & Ebvect = 16,
c                    Some error occurred while closing a file.
     & Eclose = 17,
c                    An internal limit has been exceeded.
     & Elimit = 18,
c                    A right parenthesis is missing
     & Emrpar = 19 )
c_*END_INCLUDE ./error.h
c_##include "params.h"
c_*BEGIN_INCLUDE ./params.h
c Length of the error message buffer
c (this is currently hardwired to 80 in some places)
      integer lemess
      parameter (lemess = 80)
c Maximum length of a token
      integer lentok
      parameter (lentok = 80)
c Maximum number of segments in a keyword
      integer maxseg
      parameter (maxseg = 20)
c Maximum number of current working keywords
      integer maxcwk
      parameter (maxcwk = 10)
c Maximum length of a filename
      integer maxnam
      parameter (maxnam = 100)
c Maximum number of units (if this is less than 100 then problems may occur,
c if this is more than one hundred then some portability problems may
c arise)
      integer maxunt
      parameter (maxunt = 100)
c_*END_INCLUDE ./params.h
      integer unit
      character*(*) name,type,access

      character*(maxnam) fname
      character*10 stat
      integer ios
      character*10 frmttd, int
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess
      character*1 junk
      integer gfilnm,trunit
      data frmttd/'FORMATTED'/

c Get the precise filename
      fioopn = gfilnm(name,type,fname)
      if (fioopn.ne.EOK) return

      if (access .eq. 'READ') then
        stat = 'OLD'
      else if (access .eq. 'OLD') then
        stat = 'OLD'
      else if (access .eq. 'NEW') then
        stat = 'NEW'
      else if (access .eq. 'APPEND') then
c_##ifdef AIX
c_#        stat = 'UNKNOWN'
c_##else
        stat = 'APPEND'
c_##endif
      else if (access .eq. 'WRITE') then
        stat = 'UNKNOWN'
        endif

c Assign unit number and record information about the file.
      fioopn = trunit(unit,name,frmttd,'INIT')
      if (fioopn.ne.EOK) return

c Open the file 
      open(unit,file=fname,form=frmttd,status=stat,iostat=ios)
      if (ios .ne. 0) then
        errcod = Eopen
        fioopn = errcod
        emess = 'Formatted open failed for '''
        call scat(emess,name)
        call scat(emess,''' iostat=')
        call intchr(ios,int)
        call scat(emess,int)
        call scatb(emess,'(fioopn)')
        return
        endif

      if (access .eq. 'APPEND') then
 801    read(unit,'(A1)',end=901) junk
          goto 801
 901      continue
      else if (access .eq. 'WRITE') then
        rewind unit
        endif

      return
      end
c----------------------------------------------------------------------------
c  The formatted I/O close routine
      integer function fiocls(unit,status)
c_##ifdef sun
c_#      implicit undefined (a-z)
c_##else
      implicit integer (a-z)
c_##endif
c_##include "error.h"
c_*BEGIN_INCLUDE ./error.h
      integer EOK, Etype, Ekey, Elngth, Esyntx, Einflp, Efind, EEOF
      integer Etrunc, Eseg, Eparen, Eevect, Erdrct, EEOD, Earg, Ebvect
      integer Eopen, Eclose, Elimit, Emrpar
      parameter (
c                    OK, no error.
     & EOK = 0,
c                    Bad type in a formatted read.
     & Etype = 1,
c                    Could not find the keyword.
     & Ekey  = 2,
c                    Bad array length.
     & Elngth  = 3,
c                    Syntax error.
     & Esyntx  = 4,
c                    Infinite loop
     & Einflp  = 5,
c                    Efind (not used, 6 is available for replacement)
     & Efind = 6,
c                    End of File
     & EEOF = 7,
c                    Character value has been truncated
     & Etrunc = 8,
c                    Too many segments in a keyword
     & Eseg = 9,
c                    Too many levels of parenthesis
     & Eparen = 10,
c                    End of vector (not really an error)
     & Eevect = 11,
c                    Redirection depth is too deep.
     & Erdrct = 12,
c                    End of Data is like EOF but expected.
     & EEOD = 13,
c                    Bad value in argument list
     & Earg = 14,
c                    Could not open a file.
     & Eopen = 15,
c                    Beginning of vector (not really an error)
     & Ebvect = 16,
c                    Some error occurred while closing a file.
     & Eclose = 17,
c                    An internal limit has been exceeded.
     & Elimit = 18,
c                    A right parenthesis is missing
     & Emrpar = 19 )
c_*END_INCLUDE ./error.h
c_##include "params.h"
c_*BEGIN_INCLUDE ./params.h
c Length of the error message buffer
c (this is currently hardwired to 80 in some places)
      integer lemess
      parameter (lemess = 80)
c Maximum length of a token
      integer lentok
      parameter (lentok = 80)
c Maximum number of segments in a keyword
      integer maxseg
      parameter (maxseg = 20)
c Maximum number of current working keywords
      integer maxcwk
      parameter (maxcwk = 10)
c Maximum length of a filename
      integer maxnam
      parameter (maxnam = 100)
c Maximum number of units (if this is less than 100 then problems may occur,
c if this is more than one hundred then some portability problems may
c arise)
      integer maxunt
      parameter (maxunt = 100)
c_*END_INCLUDE ./params.h
      integer unit
      character*(*) status

      integer errcod,ios
      character*30 int
      character*80 emess
      common /lparer/ errcod, emess

      integer trunit

      close(unit,status=status,err=901,iostat=ios)

c Tell the unit tracker that this unit number is once again available
      fiocls = trunit(unit,' ',' ','DELETE')
      return

 901  continue
      errcod = Eclose
      fiocls = errcod
      emess = 'Close on unit'
      call intchr(unit,int)
      call scatb(emess,int)
      call scatb(emess,'gave iostat')
      call intchr(ios,int)
      call scatb(emess,int)
      call scatb(emess,'(fiocls)')
      return
      end
c----------------------------------------------------------------------------
c  Get File Name
c This determines the name for a file based on name and type.
c name is the name given by the application
c type is one of LIB USER SUSER
c fname is set to the computed filename.
croag  - I modified the erroneous SUSER branch.  should be spiffy now.
      integer function gfilnm(name,type,fname)
c_##ifdef sun
c_#      implicit undefined (a-z)
c_##else
      implicit integer (a-z)
c_##endif
c_##include "error.h"
c_*BEGIN_INCLUDE ./error.h
      integer EOK, Etype, Ekey, Elngth, Esyntx, Einflp, Efind, EEOF
      integer Etrunc, Eseg, Eparen, Eevect, Erdrct, EEOD, Earg, Ebvect
      integer Eopen, Eclose, Elimit, Emrpar
      parameter (
c                    OK, no error.
     & EOK = 0,
c                    Bad type in a formatted read.
     & Etype = 1,
c                    Could not find the keyword.
     & Ekey  = 2,
c                    Bad array length.
     & Elngth  = 3,
c                    Syntax error.
     & Esyntx  = 4,
c                    Infinite loop
     & Einflp  = 5,
c                    Efind (not used, 6 is available for replacement)
     & Efind = 6,
c                    End of File
     & EEOF = 7,
c                    Character value has been truncated
     & Etrunc = 8,
c                    Too many segments in a keyword
     & Eseg = 9,
c                    Too many levels of parenthesis
     & Eparen = 10,
c                    End of vector (not really an error)
     & Eevect = 11,
c                    Redirection depth is too deep.
     & Erdrct = 12,
c                    End of Data is like EOF but expected.
     & EEOD = 13,
c                    Bad value in argument list
     & Earg = 14,
c                    Could not open a file.
     & Eopen = 15,
c                    Beginning of vector (not really an error)
     & Ebvect = 16,
c                    Some error occurred while closing a file.
     & Eclose = 17,
c                    An internal limit has been exceeded.
     & Elimit = 18,
c                    A right parenthesis is missing
     & Emrpar = 19 )
c_*END_INCLUDE ./error.h
c_##include "params.h"
c_*BEGIN_INCLUDE ./params.h
c Length of the error message buffer
c (this is currently hardwired to 80 in some places)
      integer lemess
      parameter (lemess = 80)
c Maximum length of a token
      integer lentok
      parameter (lentok = 80)
c Maximum number of segments in a keyword
      integer maxseg
      parameter (maxseg = 20)
c Maximum number of current working keywords
      integer maxcwk
      parameter (maxcwk = 10)
c Maximum length of a filename
      integer maxnam
      parameter (maxnam = 100)
c Maximum number of units (if this is less than 100 then problems may occur,
c if this is more than one hundred then some portability problems may
c arise)
      integer maxunt
      parameter (maxunt = 100)
c_*END_INCLUDE ./params.h
      character*(*) name,type,fname

      logical dot
      character*(lentok) keywrd
      character*(lentok) value
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess
      character*7 files
      integer frdc

croag the trailing colon is an error.  the lead colon messes up the
croag	path too.  trailing colon meant search was :files::path
croag	which always failed.  better is :files:path which parses
croag 	just that.  best is files:path, which parses
croag	(input.dat- input:files:path and then default:files:path)
croag      data files/':FILES:'/
      data files/'FILES'/

      errcod = EOK
      gfilnm = errcod
      emess = ' '

      if (type .eq. 'USER') then
croag	This file has been exactly specified by the programmer in the function
croag	call.  remains in in the job directory.
        fname = name
        if (name.ne.fname) goto 901
      else if (type .eq. 'LIB') then
c This file resides in the library directory
croag	and uses the name specified in the function call.
c        fname = '/home/xerxes/psi/lib/'
        fname = "/usr/local/psi/share/"
        call scat(fname,name)
        if (fname(len(fname):len(fname)).ne.' ') goto 901
      else if (type .eq. 'SUSER') then
c This file has been symbolically specified, read the input to find its name
c       Check for a path input
croag	if there is a path, it is used, if not then defaults to job directory.
        keywrd = files
        call scat(keywrd,':PATH')
        gfilnm = frdc(keywrd,fname)
        if (gfilnm.ne.EOK) then
          if (gfilnm.eq.EEOD) then
croag 	got to end of data, and no path was specified, so use job directory
            fname = ' '
          else
            return
            endif
          endif
c       See if there is a name input for this file
croag	rewrite starting here.  if no name, defaults to name parameter
croag	in the function call.
croag	if name is given by the user, then it looks for an extension
        keywrd = files
        call scat(keywrd,':NAME')
        gfilnm = frdc(keywrd,value)
        if (gfilnm.ne.EOK) then
          if (gfilnm.eq.EEOD) then
croag	no name was specified by user
croag  	change within this conditional
croag            value = ' '
            value = name
            call tolowr(value)
            gfilnm = EOK
            emess = ' '
croag	endchange
            dot = .false.
          else
            return
            endif
        else
croag	i think this is bogus, but the parser probably can't handle
croag	a period, if it has special meaning.  
          dot = .true.
          endif
        call scat(fname,value)
croag	many changes follow, within the "dot" conditional
c       If a name was found, then a '.' is needed between the name and the ext.
croag	it only needs the dot if an extension exists.
croag        if (dot) call scat(fname,'.')
c       Get the extension
	if (dot) then
          keywrd = files
croag 	learn to spell curt.  :-)
croag          call scat(keywrd,':EXTENTION')
          call scat(keywrd,':EXTENSION')
          gfilnm = frdc(keywrd,value)
          if (gfilnm.ne.EOK) return
croag	there was a name, but there is no extension.  fine, then we've got 
croag	the whole name already.
croag          if (gfilnm.ne.EOK) then
croag            if (gfilnm.eq.EEOD) then
croag              value = name
croag              call tolowr(value)
croag              errcod = EOK
croag              gfilnm = errcod
croag              emess = ' '
croag	      return
croag            else
croag              return
croag              endif
croag            endif
          call scat(fname,'.')
          call scat(fname,value)
	endif
c       We should have the complete filename now
        if (fname(len(fname):len(fname)) .ne. ' ') goto 901
      else
croag  it's not USER, LIB or SUSER
        errcod = Earg
        gfilnm = errcod
        emess = 'Bad file type argument:'
        call scatb(emess,type)
        call scatb(emess,'(gfilnm)')
        endif

      return

c Here I handle truncation errors
 901  continue
      errcod = Etrunc
      gfilnm = errcod
      emess = 'A filename was truncated (gfilnm)'
      return
      end
c--------------------------------------------------------------------------
c  Formatted I/O delete
c  This takes a filename and type and deletes it.  The file should NOT
c already be open.  If it is already open then use its unit number.
c The formatted file open routine is used, but I don't think that will
c cause any problem, since the file is immediately closed after being
c deleted.
c Error recovery: ignores all errors.
      subroutine fiodel(name,type)
c_##ifdef sun
c_#      implicit undefined (a-z)
c_##else
      implicit integer (a-z)
c_##endif
c_##include "error.h"
c_*BEGIN_INCLUDE ./error.h
      integer EOK, Etype, Ekey, Elngth, Esyntx, Einflp, Efind, EEOF
      integer Etrunc, Eseg, Eparen, Eevect, Erdrct, EEOD, Earg, Ebvect
      integer Eopen, Eclose, Elimit, Emrpar
      parameter (
c                    OK, no error.
     & EOK = 0,
c                    Bad type in a formatted read.
     & Etype = 1,
c                    Could not find the keyword.
     & Ekey  = 2,
c                    Bad array length.
     & Elngth  = 3,
c                    Syntax error.
     & Esyntx  = 4,
c                    Infinite loop
     & Einflp  = 5,
c                    Efind (not used, 6 is available for replacement)
     & Efind = 6,
c                    End of File
     & EEOF = 7,
c                    Character value has been truncated
     & Etrunc = 8,
c                    Too many segments in a keyword
     & Eseg = 9,
c                    Too many levels of parenthesis
     & Eparen = 10,
c                    End of vector (not really an error)
     & Eevect = 11,
c                    Redirection depth is too deep.
     & Erdrct = 12,
c                    End of Data is like EOF but expected.
     & EEOD = 13,
c                    Bad value in argument list
     & Earg = 14,
c                    Could not open a file.
     & Eopen = 15,
c                    Beginning of vector (not really an error)
     & Ebvect = 16,
c                    Some error occurred while closing a file.
     & Eclose = 17,
c                    An internal limit has been exceeded.
     & Elimit = 18,
c                    A right parenthesis is missing
     & Emrpar = 19 )
c_*END_INCLUDE ./error.h
c_##include "params.h"
c_*BEGIN_INCLUDE ./params.h
c Length of the error message buffer
c (this is currently hardwired to 80 in some places)
      integer lemess
      parameter (lemess = 80)
c Maximum length of a token
      integer lentok
      parameter (lentok = 80)
c Maximum number of segments in a keyword
      integer maxseg
      parameter (maxseg = 20)
c Maximum number of current working keywords
      integer maxcwk
      parameter (maxcwk = 10)
c Maximum length of a filename
      integer maxnam
      parameter (maxnam = 100)
c Maximum number of units (if this is less than 100 then problems may occur,
c if this is more than one hundred then some portability problems may
c arise)
      integer maxunt
      parameter (maxunt = 100)
c_*END_INCLUDE ./params.h
      character*(*) name, type

      integer errcod,unit

      integer fioopn,fiocls

      if (type.eq.'LIB') return

      unit = 0
      errcod = fioopn(unit,name,type,'WRITE')
      if (errcod.ne.EOK) return

      errcod = fiocls(unit,'DELETE')
      return
      end

