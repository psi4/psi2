      program main
      integer status
      integer errcod, frdc, chrlen
      character*80 name
c The output must be to the terminal for this program to work with
c the shell script that actually deletes the file.
      call psinit('TERMINAL!')
      call gfname(status)
      if (status.eq.1) return
      errcod = frdc(':DEFAULT:FILES:DEFAULT:NAME',name)
      if (errcod.eq.0) write(6,'(A)') name(1:chrlen(name))
      call psdone()
      end

      subroutine gprgid(progid)
      character*(*) progid
      progid = 'FILESNAME'
      return
      end
