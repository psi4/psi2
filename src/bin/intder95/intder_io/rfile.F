CC//////////////////////////////////////////////////////////////////////

      SUBROUTINE RFILE(ITAPE)
C
C USE THIS TO OPEN ALL DIRECT ACCESS FILES
C
      IMPLICIT INTEGER (A-Z)
	
C
      parameter (reclen = 1024, numbuf = 10)
      common /pointr/ wptr(128),tptr(128)
      COMMON /SECT/ SECTOR
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &   untbuf(128), bufrec(numbuf), buffer(reclen,numbuf)
      character*1 digit
      character*71 file
      character*80 prefix
      character*71 scrtch
      data first/1/
C
      ITAP3 = 3
      ITAP6 = 6

c Initialize on the first call to rfile.
      if (first .eq. 1) then
        first = 0
        do 101 i = 1, numbuf
          bufunt(i) = 0
          bufrec(i) = 0
          wrtbuf(i) = 0
  101     continue
        do 102 I = 1, 128
          untbuf(i) = 0
  102     continue
        nxtbuf = 1
        endif
c
c  Some programs require SECTOR = 1024
      SECTOR = 1024

      call locate(5,'# FILES ##',ierr)
      if (ierr .ne. 0) then
        iprfx = 0
      else
 804    continue
        prefix = ' '
        read(5,'(A80)',err=904) prefix

 904    continue
        iprfx = 0
 805    continue
        iprfx = iprfx + 1
        if (iprfx .eq. 81) goto 905
        if (prefix(iprfx:iprfx).ne.' ') goto 805
 905    continue
        iprfx = iprfx - 1
        endif
      file='file'//digit(itape/10)//digit(itape-10*(itape/10))
      if (iprfx .ne. 0) then
        scrtch = file
        file=prefix(1:iprfx)//'.'//scrtch
        endif



      open (unit=itape,access='DIRECT',recl=reclen*4,form='UNFORMATTED',
     &      iostat=ierr,status='UNKNOWN',file=file)
C
      IF (IERR .ne. 0) then
        WRITE(*,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        WRITE(6,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        write(*,*) ' rfile: filename=',file
        write(6,*) ' rfile: filename=',file
        call ioerr(ierr)
        call mabort
        endif
C
C PLACE THE POINTER OF ITAPE AT THE BEGINNING
C
      CALL SREW(ITAPE)
C
      RETURN
      END
