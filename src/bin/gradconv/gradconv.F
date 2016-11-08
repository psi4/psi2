      program gradcv
cets110289      program gradcv
clj      subroutine fentry(zdum,izdum,maxcor)
      implicit real*8 (a-h, o-z)
      parameter (maxatm=300)
      character*10 aster
      character*1 junk
      dimension x1(maxatm,3), x2(maxatm,3)
      logical onegrd
 
      call psinit('APPEND')
      call ffile(11,' ',0)

      call tstart(6)
      write(6,*)
      write(6,fmt='(t20, a)') 'Geometry Convergence Test Program'
      write(6,*)
 
      call locate(5,'# GRADCV #',ierror)
      if (ierror .ne. 0) then
c Here is where the defaults are set:
          grmsi = -1.0d-8
          gmaxi =  0.0d0
          grmsc =  0.2d-5
          gmaxc =  0.4d-5
          goto 701
          endif
    3 format(f20.10)
      read(5,*) grmsi
      read(5,*) gmaxi
      read(5,*) grmsc
      read(5,*) gmaxc
  701 write(6,*) 'Gradconv will use the following input:'
      write(6,*) '     grmsi = ', grmsi
      write(6,*) '     gmaxi = ', gmaxi
      write(6,*) '     grmsc = ', grmsc
      write(6,*) '     gmaxc = ', gmaxc
      write(6,*)
 
      length = 0
  801 read (11,'(1A1)',end=901) junk
          length = length + 1
          goto 801
  901     continue
 
  802 continue
          backspace 11
          backspace 11
          read(11,'(I5)') natoms
          if (natoms .eq. 0) goto 802
      if (natoms .gt. maxatm) then
          write(6,*) 'Too many atoms for gradconv. Aborting.'
#if defined(NO_EXIT_ROUTINE)
          stop 111
#else
          call exit(111)
#endif
          endif
c See if only one gradient run is in file11.
      if (length .eq. 2*natoms + 2) then
          onegrd = .true.
      else
          onegrd = .false.
          endif
 
      write(6,*) 'Read natoms = ',natoms
      if (.not. onegrd) then
          do 101 i = 1, natoms + 2
  101         backspace 11
 
          do 102 i = 1, natoms
  102         read (11,1) (x1(i,j), j=1, 3)
    1         format(20x, 3f20.10)
 
          do 103 i = 1, 2
  103         read (11,*)
          endif
      do 109 i = 1, natoms
          read (11,*)
  109     continue
 
      do 104 i = 1, natoms
          read (11,1) (x2(i,j), j=1, 3)
  104     continue
 
      if (.not. onegrd) then
          grms1 = 0.d0
          gmax1 = 0.d0
          do 107 i = 1, natoms
              do 207 j = 1, 3
                  if ( abs(x1(i,j)) .gt. gmax1 ) gmax1 = abs(x1(i,j))
                  grms1 = x1(i,j)**2 + grms1
  207             continue
  107         continue
          grms1 = dsqrt(grms1/(natoms*3.d0))
      else
c Note: This sets the limit on how lousy the first run can be:
          grms1 = 0.75d0
          gmax1 = 1.00d0
          endif
 
      grms2 = 0.d0
      gmax2 = 0.d0
      do 108 i = 1, natoms
          do 208 j = 1, 3
              if ( abs(x2(i,j)) .gt. gmax2 ) gmax2 = abs(x2(i,j))
              grms2 = grms2 + x2(i,j)**2
  208         continue
  108     continue
      grms2 = dsqrt(grms2/(natoms*3.d0))
 
      if (.not. onegrd) then
      write(6,*)
      write(6,*) 'For the next to last gradient run:'
      write(6,*) '   grms = ',grms1
      write(6,*) '   gmax = ',gmax1
      endif
      write(6,*)
      write(6,*) 'For the last gradient run:'
      write(6,*) '   grms = ',grms2
      write(6,*) '   gmax = ',gmax2
      write(6,*)
      write(6,*) ' where grms = the root-mean-square of the cartesian'
      write(6,*) '              gradients'
      write(6,*) '   and gmax = the maximum absolute value of the'
      write(6,*) '              cartesian gradients.'
      write(6,*) 'The units are as in file11.'
      write(6,*)
 
      if ( (gmax2 .lt. gmaxc) .and. (grms2 .lt. grmsc) ) then
          write(6,*) 'Converged.'
          write(6,*)
          call tstop(6)
          call psdone()
#if defined(NO_EXIT_ROUTINE)
          stop 11
#else
          call exit(11)
#endif
      else  if (       (gmax2 .gt. gmax1+gmaxi )
     &           .or.  (grms2 .gt. grms1+grmsi ) ) then
          write(6,*) 'Calculation aborted.'
          write(6,*)
          call tstop(6)
          call psdone()
#if defined(NO_EXIT_ROUTINE)
          stop 111
#else
          call exit(111)
#endif
      else
          write(6,*) 'Calculation continuing.'
          write(6,*)
          call tstop(6)
          call psdone()
#if defined(NO_EXIT_ROUTINE)
          stop 0
#else
          call exit(0)
#endif
          endif

      end
